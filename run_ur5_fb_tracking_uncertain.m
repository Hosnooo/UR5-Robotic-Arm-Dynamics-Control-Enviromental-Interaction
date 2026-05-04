clc; clear; close all;
addpath('maple_gen')
addpath('utils')

%% ========================= Setup =========================
assert(exist('UR5_params','file') == 2, 'UR5_params.m not found.');
assert(exist('UR5_M','file')      == 2, 'UR5_M.m not found.');
assert(exist('UR5_C','file')      == 2, 'UR5_C.m not found.');
assert(exist('UR5_h','file')      == 2, 'UR5_h.m not found.');
assert(exist('UR5_G','file')      == 2, 'UR5_G.m not found.');
assert(exist('UR5_fdyn','file')   == 2, 'UR5_fdyn.m not found.');

USE_FIG_HELPERS = (exist('figureoptscall','file') == 2) && ...
                  (exist('saveFigureAsPDF','file') == 2);

if USE_FIG_HELPERS
    figureoptscall;
end

outdir = 'ur5_fb_uncertainty_results';
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

%% ========================= Robot parameters =========================
[~, Pi_true] = UR5_params();
n = 6;

%% ========================= Simulation settings =========================
Tf = 15.0;
t_plot = linspace(0, Tf, 3001);

ode_opts = odeset( ...
    'RelTol', 1e-6, ...
    'AbsTol', 1e-8, ...
    'MaxStep', 1e-2);

%% ========================= Trajectory selection =========================
% Choose:
%   'quintic'  -> one smooth point-to-point move, then hold
%   'sine'     -> continuous sinusoidal motion
TRAJ_MODE = 'sine';

traj = struct();
traj.mode = lower(TRAJ_MODE);

switch traj.mode
    case 'quintic'
        traj.Tmove   = 10.0;
        traj.q_start = [ 0.3; -0.8;  0.6; -0.4;  0.2; -0.1];
        traj.q_final = [ 0.0; -0.5;  0.8;  0.0;  0.3;  0.0];

    case 'sine'
        traj.q_bias  = [ 0.20; -0.70;  0.70; -0.20;  0.15;  0.00];
        traj.q_amp   = [ 0.35;  0.25;  0.20;  0.30;  0.15;  0.20];
        traj.w_ref   = 2*pi*[0.12; 0.10; 0.14; 0.09; 0.11; 0.13];
        traj.phi_ref = [ 0.0;  pi/4;  pi/2;  pi/6;  pi/3;  pi/5];

    otherwise
        error('Unknown TRAJ_MODE. Use ''quintic'' or ''sine''.');
end

%% ========================= Initial state =========================
% Start from the desired trajectory at t = 0, then add a chosen bias.
[q_ref0, dq_ref0, ~] = desired_trajectory_general(0, traj);

q0_bias  = [ 0.15; -0.10;  0.08;  0.12; -0.06;  0.10];
dq0_bias = zeros(n,1);

q0  = q_ref0  + q0_bias;
dq0 = dq_ref0 + dq0_bias;

%% ========================= Controller gains =========================
% Baseline FFW+PD gains
Kp_vec = [10; 10;  8; 5; 3; 2];
Kd_vec = [ 4;  4;  3; 2; 1; 1];

% Updated-law gains
Lambda_vec = [3.0; 3.0; 2.5; 2.0; 1.5; 1.2];
Kr_vec     = [8.0; 8.0; 6.0; 4.0; 2.5; 2.0];

settling_tol_abs = 1e-3;

%% ========================= Uncertainty settings =========================
% Use the same perturbation pattern across all uncertainty levels.
rng(2);
delta_pattern = 2*rand(size(Pi_true)) - 1;
delta_pattern(end) = 0;   % keep gravity constant exact

% Keep 5% in the list. Add another level to examine the effect away from 5%.
UNC_LEVELS = [0.05, 0.10];

%% ========================= Case list =========================
cases = struct([]);

cases(1).name      = 'Exact_FFW_PD';
cases(1).label     = '$\mathrm{Exact\ FFW+PD}$';
cases(1).mode      = 'EXACT_FFW_PD';
cases(1).unc_level = 0.0;

idx = 2;
for k = 1:numel(UNC_LEVELS)
    u = UNC_LEVELS(k);

    cases(idx).name      = sprintf('Uncertain_FFW_PD_%02d', round(100*u));
    cases(idx).label     = sprintf('$\\mathrm{Uncertain\\ FFW+PD}\\ (%.0f\\%%)$', 100*u);
    cases(idx).mode      = 'UNCERTAIN_FFW_PD';
    cases(idx).unc_level = u;
    idx = idx + 1;

    cases(idx).name      = sprintf('UpdatedLaw_%02d', round(100*u));
    cases(idx).label     = sprintf('$\\mathrm{Updated\\ law}\\ (%.0f\\%%)$', 100*u);
    cases(idx).mode      = 'UPDATED_LAW';
    cases(idx).unc_level = u;
    idx = idx + 1;
end

nCases = numel(cases);

%% ========================= Run simulations =========================
results = struct([]);
x0 = [q0; dq0];

for c = 1:nCases
    fprintf('Running %s ...\n', cases(c).name);

    Pi_hat = build_uncertain_params(Pi_true, cases(c).unc_level, delta_pattern);

    rhs = @(t, x) closed_loop_rhs_uncertain( ...
        t, x, cases(c), Pi_true, Pi_hat, traj, ...
        Kp_vec, Kd_vec, Lambda_vec, Kr_vec);

    sol = ode45(rhs, [0 Tf], x0, ode_opts);

    sim = postprocess_solution_uncertain( ...
        sol, t_plot, cases(c), Pi_true, Pi_hat, traj, ...
        Kp_vec, Kd_vec, Lambda_vec, Kr_vec);

    metrics = compute_metrics_uncertain(sim, settling_tol_abs, traj);

    results(c).name      = cases(c).name;
    results(c).label     = cases(c).label;
    results(c).mode      = cases(c).mode;
    results(c).unc_level = cases(c).unc_level;
    results(c).sim       = sim;
    results(c).metrics   = metrics;

    write_timeseries_csv_uncertain(results(c), outdir);
end

%% ========================= Summary table =========================
summaryTbl = build_summary_table_uncertain(results);

summaryTbl = sortrows(summaryTbl, ...
    {'RMS_PosErr_Norm', 'RMS_VelErr_Norm', 'Control_Energy'}, ...
    {'ascend',          'ascend',          'ascend'});

disp(' ');
disp(summaryTbl);

writetable(summaryTbl, fullfile(outdir, 'summary_fb_uncertainty.csv'));
save(fullfile(outdir, 'summary_fb_uncertainty.mat'), ...
    'results', 'summaryTbl', ...
    't_plot', 'Tf', 'traj', ...
    'q0', 'dq0', ...
    'Kp_vec', 'Kd_vec', 'Lambda_vec', 'Kr_vec', ...
    'UNC_LEVELS', 'delta_pattern');

%% ========================= Plot styles =========================
styles = {'-','--','-.',':','-','--','-.'};
legend_labels = {results.label};

%% ========================= Figure 1: Overview =========================
fig1 = figure('Color','w', 'Name', 'FB Uncertainty Overview');
tl1 = tiledlayout(fig1, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

hLeg1 = gobjects(1, nCases);

ax1 = nexttile(tl1); hold(ax1, 'on');
for c = 1:nCases
    h = plot(ax1, results(c).sim.t, results(c).sim.pos_err_norm, ...
        'LineStyle', styles{c}, 'LineWidth', 1.5);
    hLeg1(c) = h;
end
grid(ax1, 'on');
xlabel(ax1, 'Time [s]');
ylabel(ax1, '$\|e_q\|_2$ [rad]');

ax2 = nexttile(tl1); hold(ax2, 'on');
for c = 1:nCases
    plot(ax2, results(c).sim.t, results(c).sim.vel_err_norm, ...
        'LineStyle', styles{c}, 'LineWidth', 1.5);
end
grid(ax2, 'on');
xlabel(ax2, 'Time [s]');
ylabel(ax2, '$\|e_{\dot q}\|_2$ [rad/s]');

ax3 = nexttile(tl1); hold(ax3, 'on');
for c = 1:nCases
    plot(ax3, results(c).sim.t, results(c).sim.tau_norm, ...
        'LineStyle', styles{c}, 'LineWidth', 1.5);
end
grid(ax3, 'on');
xlabel(ax3, 'Time [s]');
ylabel(ax3, '$\|\tau\|_2$ [Nm]');

lgd1 = legend(hLeg1, legend_labels, 'Orientation', 'horizontal');
lgd1.Layout.Tile = 'north';
lgd1.NumColumns = 3;

if USE_FIG_HELPERS
    saveFigureAsPDF(fig1, fullfile(outdir, 'fig_overview.pdf'));
end

%% ========================= Figure 2: Joint position errors =========================
fig2 = figure('Color','w', 'Name', 'Joint Position Errors');
tl2 = tiledlayout(fig2, 3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

hLeg2 = gobjects(1, nCases);

for j = 1:n
    ax = nexttile(tl2); hold(ax, 'on');
    for c = 1:nCases
        h = plot(ax, results(c).sim.t, results(c).sim.e_log(j,:), ...
            'LineStyle', styles{c}, 'LineWidth', 1.35);
        if j == 1
            hLeg2(c) = h;
        end
    end
    yline(ax, 0, 'k:', 'LineWidth', 1.2);
    grid(ax, 'on');
    xlabel(ax, 'Time [s]');
    ylabel(ax, sprintf('$e_{q,%d}$ [rad]', j));
end

lgd2 = legend(hLeg2, legend_labels, 'Orientation', 'horizontal');
lgd2.Layout.Tile = 'north';
lgd2.NumColumns = 3;

if USE_FIG_HELPERS
    saveFigureAsPDF(fig2, fullfile(outdir, 'fig_joint_position_errors.pdf'));
end

%% ========================= Figure 3: Joint velocity errors =========================
fig3 = figure('Color','w', 'Name', 'Joint Velocity Errors');
tl3 = tiledlayout(fig3, 3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

hLeg3 = gobjects(1, nCases);

for j = 1:n
    ax = nexttile(tl3); hold(ax, 'on');
    for c = 1:nCases
        h = plot(ax, results(c).sim.t, results(c).sim.de_log(j,:), ...
            'LineStyle', styles{c}, 'LineWidth', 1.35);
        if j == 1
            hLeg3(c) = h;
        end
    end
    yline(ax, 0, 'k:', 'LineWidth', 1.2);
    grid(ax, 'on');
    xlabel(ax, 'Time [s]');
    ylabel(ax, sprintf('$e_{\\dot q,%d}$ [rad/s]', j));
end

lgd3 = legend(hLeg3, legend_labels, 'Orientation', 'horizontal');
lgd3.Layout.Tile = 'north';
lgd3.NumColumns = 3;

if USE_FIG_HELPERS
    saveFigureAsPDF(fig3, fullfile(outdir, 'fig_joint_velocity_errors.pdf'));
end

%% ========================= Figure 4: Joint positions =========================
fig4 = figure('Color','w', 'Name', 'Joint Positions');
tl4 = tiledlayout(fig4, 3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

hLeg4 = gobjects(1, nCases);
hqd4 = gobjects(1,1);

for j = 1:n
    ax = nexttile(tl4); hold(ax, 'on');
    for c = 1:nCases
        h = plot(ax, results(c).sim.t, results(c).sim.q_log(j,:), ...
            'LineStyle', styles{c}, 'LineWidth', 1.35);
        if j == 1
            hLeg4(c) = h;
        end
    end
    hqd = plot(ax, results(1).sim.t, results(1).sim.qd_log(j,:), ...
        'k:', 'LineWidth', 1.2);
    if j == 1
        hqd4 = hqd;
    end
    grid(ax, 'on');
    xlabel(ax, 'Time [s]');
    ylabel(ax, sprintf('$q_%d$ [rad]', j));
end

legend_entries4 = [legend_labels, {'$q_d$'}];
lgd4 = legend([hLeg4, hqd4], legend_entries4, 'Orientation', 'horizontal');
lgd4.Layout.Tile = 'north';
lgd4.NumColumns = 3;

if USE_FIG_HELPERS
    saveFigureAsPDF(fig4, fullfile(outdir, 'fig_joint_positions.pdf'));
end

%% ========================= Figure 5: Joint torques =========================
fig5 = figure('Color','w', 'Name', 'Joint Torques');
tl5 = tiledlayout(fig5, 3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

hLeg5 = gobjects(1, nCases);

for j = 1:n
    ax = nexttile(tl5); hold(ax, 'on');
    for c = 1:nCases
        h = plot(ax, results(c).sim.t, results(c).sim.tau_log(j,:), ...
            'LineStyle', styles{c}, 'LineWidth', 1.35);
        if j == 1
            hLeg5(c) = h;
        end
    end
    grid(ax, 'on');
    xlabel(ax, 'Time [s]');
    ylabel(ax, sprintf('$\\tau_%d$ [Nm]', j));
end

lgd5 = legend(hLeg5, legend_labels, 'Orientation', 'horizontal');
lgd5.Layout.Tile = 'north';
lgd5.NumColumns = 3;

if USE_FIG_HELPERS
    saveFigureAsPDF(fig5, fullfile(outdir, 'fig_joint_torques.pdf'));
end

%% ========================= Figure 6: Sliding variable norm =========================
fig6 = figure('Color','w', 'Name', 'Sliding Variable Norm');
tl6 = tiledlayout(fig6, 1, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

hLeg6 = gobjects(1, nCases);

ax6 = nexttile(tl6); hold(ax6, 'on');
for c = 1:nCases
    h = plot(ax6, results(c).sim.t, results(c).sim.s_norm, ...
        'LineStyle', styles{c}, 'LineWidth', 1.5);
    hLeg6(c) = h;
end
grid(ax6, 'on');
xlabel(ax6, 'Time [s]');
ylabel(ax6, '$\|s\|_2$');

lgd6 = legend(hLeg6, legend_labels, 'Orientation', 'horizontal');
lgd6.Layout.Tile = 'north';
lgd6.NumColumns = 3;

if USE_FIG_HELPERS
    saveFigureAsPDF(fig6, fullfile(outdir, 'fig_s_norm.pdf'));
end

fprintf('\nSaved results to:\n%s\n', fullfile(pwd, outdir));

%% ========================= Local functions =========================

function Pi_hat = build_uncertain_params(Pi_true, unc_level, delta_pattern)
Pi_hat = Pi_true .* (1 + unc_level * delta_pattern);
Pi_hat(end) = Pi_true(end);   % keep gravity exact
end

function dx = closed_loop_rhs_uncertain(t, x, caseDef, Pi_true, Pi_hat, traj, ...
                                        Kp_vec, Kd_vec, Lambda_vec, Kr_vec)

n = 6;

q  = x(1:n);
dq = x(n+1:2*n);

[qd, dqd, ddqd] = desired_trajectory_general(t, traj);

[~, ~, tau] = uncertain_control_law( ...
    q, dq, qd, dqd, ddqd, caseDef, Pi_hat, ...
    Kp_vec, Kd_vec, Lambda_vec, Kr_vec);

ddq = UR5_fdyn(q, dq, tau, Pi_true);

dx = [dq; ddq];
end

function sim = postprocess_solution_uncertain(sol, t_plot, caseDef, Pi_true, Pi_hat, traj, ...
                                              Kp_vec, Kd_vec, Lambda_vec, Kr_vec)

X = deval(sol, t_plot);
n = 6;
N = numel(t_plot);

q_log    = X(1:n, :);
dq_log   = X(n+1:2*n, :);

qd_log   = zeros(n, N);
dqd_log  = zeros(n, N);
ddqd_log = zeros(n, N);
ddq_log  = zeros(n, N);
tau_log  = zeros(n, N);
e_log    = zeros(n, N);
de_log   = zeros(n, N);
s_log    = zeros(n, N);

for k = 1:N
    t  = t_plot(k);
    q  = q_log(:,k);
    dq = dq_log(:,k);

    [qd, dqd, ddqd] = desired_trajectory_general(t, traj);

    [e, de, tau, s] = uncertain_control_law( ...
        q, dq, qd, dqd, ddqd, caseDef, Pi_hat, ...
        Kp_vec, Kd_vec, Lambda_vec, Kr_vec);

    ddq = UR5_fdyn(q, dq, tau, Pi_true);

    qd_log(:,k)   = qd;
    dqd_log(:,k)  = dqd;
    ddqd_log(:,k) = ddqd;
    ddq_log(:,k)  = ddq;
    tau_log(:,k)  = tau;
    e_log(:,k)    = e;
    de_log(:,k)   = de;
    s_log(:,k)    = s;
end

sim.t            = t_plot(:).';
sim.q_log        = q_log;
sim.dq_log       = dq_log;
sim.ddq_log      = ddq_log;
sim.qd_log       = qd_log;
sim.dqd_log      = dqd_log;
sim.ddqd_log     = ddqd_log;
sim.tau_log      = tau_log;
sim.e_log        = e_log;
sim.de_log       = de_log;
sim.s_log        = s_log;
sim.pos_err_norm = vecnorm(e_log, 2, 1);
sim.vel_err_norm = vecnorm(de_log, 2, 1);
sim.tau_norm     = vecnorm(tau_log, 2, 1);
sim.s_norm       = vecnorm(s_log, 2, 1);
end

function [e, de, tau, s] = uncertain_control_law( ...
    q, dq, qd, dqd, ddqd, caseDef, Pi_hat, ...
    Kp_vec, Kd_vec, Lambda_vec, Kr_vec)

e  = qd  - q;
de = dqd - dq;
s  = zeros(6,1);

switch caseDef.mode
    case 'EXACT_FFW_PD'
        tau_ff = UR5_M(qd, Pi_hat) * ddqd + UR5_h(qd, dqd, Pi_hat) + UR5_G(qd, Pi_hat);
        tau_fb = Kp_vec .* e + Kd_vec .* de;
        tau = tau_ff + tau_fb;

    case 'UNCERTAIN_FFW_PD'
        tau_ff = UR5_M(qd, Pi_hat) * ddqd + UR5_h(qd, dqd, Pi_hat) + UR5_G(qd, Pi_hat);
        tau_fb = Kp_vec .* e + Kd_vec .* de;
        tau = tau_ff + tau_fb;

    case 'UPDATED_LAW'
        qrd   = dqd  + Lambda_vec .* e;
        q_rdd = ddqd + Lambda_vec .* de;
        s     = dq - qrd;

        tau = UR5_M(q, Pi_hat) * q_rdd ...
            + UR5_C(q, dq, Pi_hat) * qrd ...
            + UR5_G(q, Pi_hat) ...
            - Kr_vec .* s;

    otherwise
        error('Unknown case mode: %s', caseDef.mode);
end
end

function [qd, dqd, ddqd] = desired_trajectory_general(t, traj)

switch lower(traj.mode)
    case 'quintic'
        q0 = traj.q_start;
        qf = traj.q_final;
        Tmove = traj.Tmove;

        n = numel(q0);
        qd   = zeros(n,1);
        dqd  = zeros(n,1);
        ddqd = zeros(n,1);

        if t <= 0
            qd = q0;
            return;
        elseif t >= Tmove
            qd = qf;
            return;
        end

        s = t / Tmove;

        sigma   = 10*s^3 - 15*s^4 + 6*s^5;
        dsigma  = (30*s^2 - 60*s^3 + 30*s^4) / Tmove;
        ddsigma = (60*s - 180*s^2 + 120*s^3) / (Tmove^2);

        dq_total = qf - q0;

        qd   = q0 + dq_total .* sigma;
        dqd  = dq_total .* dsigma;
        ddqd = dq_total .* ddsigma;

    case 'sine'
        qd   = traj.q_bias + traj.q_amp .* sin(traj.w_ref .* t + traj.phi_ref);
        dqd  = traj.q_amp  .* traj.w_ref .* cos(traj.w_ref .* t + traj.phi_ref);
        ddqd = -traj.q_amp .* (traj.w_ref.^2) .* sin(traj.w_ref .* t + traj.phi_ref);

    otherwise
        error('Unknown trajectory mode.');
end
end

function metrics = compute_metrics_uncertain(sim, settling_tol_abs, traj)

t            = sim.t;
pos_err_norm = sim.pos_err_norm;
vel_err_norm = sim.vel_err_norm;
tau_norm     = sim.tau_norm;
s_norm       = sim.s_norm;

metrics.final_pos_err_norm = pos_err_norm(end);
metrics.rms_pos_err_norm   = sqrt(trapz(t, pos_err_norm.^2) / (t(end) - t(1)));
metrics.peak_pos_err_norm  = max(pos_err_norm);

metrics.final_vel_err_norm = vel_err_norm(end);
metrics.rms_vel_err_norm   = sqrt(trapz(t, vel_err_norm.^2) / (t(end) - t(1)));
metrics.peak_vel_err_norm  = max(vel_err_norm);

metrics.final_tau_norm     = tau_norm(end);
metrics.rms_tau_norm       = sqrt(trapz(t, tau_norm.^2) / (t(end) - t(1)));
metrics.peak_tau_norm      = max(tau_norm);

metrics.final_s_norm       = s_norm(end);
metrics.rms_s_norm         = sqrt(trapz(t, s_norm.^2) / (t(end) - t(1)));

metrics.control_energy     = trapz(t, tau_norm.^2);

if strcmpi(traj.mode, 'quintic')
    settle_thresh = max(0.02 * pos_err_norm(1), settling_tol_abs);

    settle_idx = NaN;
    for k = 1:numel(t)
        if all(pos_err_norm(k:end) <= settle_thresh)
            settle_idx = k;
            break;
        end
    end

    if isnan(settle_idx)
        metrics.settling_time = NaN;
    else
        metrics.settling_time = t(settle_idx);
    end
else
    metrics.settling_time = NaN;
end
end

function tbl = build_summary_table_uncertain(results)

nCases = numel(results);

Controller        = cell(nCases,1);
Label             = cell(nCases,1);
Uncertainty_pct   = zeros(nCases,1);
Final_PosErr_Norm = zeros(nCases,1);
RMS_PosErr_Norm   = zeros(nCases,1);
Peak_PosErr_Norm  = zeros(nCases,1);
Final_VelErr_Norm = zeros(nCases,1);
RMS_VelErr_Norm   = zeros(nCases,1);
Peak_VelErr_Norm  = zeros(nCases,1);
Settling_Time_s   = NaN(nCases,1);
Final_Tau_Norm    = zeros(nCases,1);
RMS_Tau_Norm      = zeros(nCases,1);
Peak_Tau_Norm     = zeros(nCases,1);
Final_s_Norm      = zeros(nCases,1);
RMS_s_Norm        = zeros(nCases,1);
Control_Energy    = zeros(nCases,1);

for c = 1:nCases
    Controller{c}         = results(c).name;
    Label{c}              = results(c).label;
    Uncertainty_pct(c)    = 100 * results(c).unc_level;
    Final_PosErr_Norm(c)  = results(c).metrics.final_pos_err_norm;
    RMS_PosErr_Norm(c)    = results(c).metrics.rms_pos_err_norm;
    Peak_PosErr_Norm(c)   = results(c).metrics.peak_pos_err_norm;
    Final_VelErr_Norm(c)  = results(c).metrics.final_vel_err_norm;
    RMS_VelErr_Norm(c)    = results(c).metrics.rms_vel_err_norm;
    Peak_VelErr_Norm(c)   = results(c).metrics.peak_vel_err_norm;
    Settling_Time_s(c)    = results(c).metrics.settling_time;
    Final_Tau_Norm(c)     = results(c).metrics.final_tau_norm;
    RMS_Tau_Norm(c)       = results(c).metrics.rms_tau_norm;
    Peak_Tau_Norm(c)      = results(c).metrics.peak_tau_norm;
    Final_s_Norm(c)       = results(c).metrics.final_s_norm;
    RMS_s_Norm(c)         = results(c).metrics.rms_s_norm;
    Control_Energy(c)     = results(c).metrics.control_energy;
end

tbl = table(Controller, Label, Uncertainty_pct, ...
    Final_PosErr_Norm, RMS_PosErr_Norm, Peak_PosErr_Norm, ...
    Final_VelErr_Norm, RMS_VelErr_Norm, Peak_VelErr_Norm, ...
    Settling_Time_s, ...
    Final_Tau_Norm, RMS_Tau_Norm, Peak_Tau_Norm, ...
    Final_s_Norm, RMS_s_Norm, ...
    Control_Energy);
end

function write_timeseries_csv_uncertain(result, outdir)

t            = result.sim.t(:);
q_log        = result.sim.q_log.';
dq_log       = result.sim.dq_log.';
ddq_log      = result.sim.ddq_log.';
qd_log       = result.sim.qd_log.';
dqd_log      = result.sim.dqd_log.';
ddqd_log     = result.sim.ddqd_log.';
tau_log      = result.sim.tau_log.';
e_log        = result.sim.e_log.';
de_log       = result.sim.de_log.';
s_log        = result.sim.s_log.';
pos_err_norm = result.sim.pos_err_norm(:);
vel_err_norm = result.sim.vel_err_norm(:);
tau_norm     = result.sim.tau_norm(:);
s_norm       = result.sim.s_norm(:);

varNames = [{'time'}, ...
    arrayfun(@(i) sprintf('q%d', i),     1:6, 'UniformOutput', false), ...
    arrayfun(@(i) sprintf('dq%d', i),    1:6, 'UniformOutput', false), ...
    arrayfun(@(i) sprintf('ddq%d', i),   1:6, 'UniformOutput', false), ...
    arrayfun(@(i) sprintf('qd%d', i),    1:6, 'UniformOutput', false), ...
    arrayfun(@(i) sprintf('dqd%d', i),   1:6, 'UniformOutput', false), ...
    arrayfun(@(i) sprintf('ddqd%d', i),  1:6, 'UniformOutput', false), ...
    arrayfun(@(i) sprintf('tau%d', i),   1:6, 'UniformOutput', false), ...
    arrayfun(@(i) sprintf('eq%d', i),    1:6, 'UniformOutput', false), ...
    arrayfun(@(i) sprintf('edq%d', i),   1:6, 'UniformOutput', false), ...
    arrayfun(@(i) sprintf('s%d', i),     1:6, 'UniformOutput', false), ...
    {'pos_err_norm', 'vel_err_norm', 'tau_norm', 's_norm'}];

data = [t, q_log, dq_log, ddq_log, qd_log, dqd_log, ddqd_log, ...
        tau_log, e_log, de_log, s_log, ...
        pos_err_norm, vel_err_norm, tau_norm, s_norm];

T = array2table(data, 'VariableNames', varNames);
writetable(T, fullfile(outdir, ['timeseries_' result.name '.csv']));
end