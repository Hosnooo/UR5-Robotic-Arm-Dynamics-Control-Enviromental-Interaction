clc; clear; close all;
addpath('maple_gen')
addpath('utils')

%% ========================= Setup =========================
assert(exist('UR5_params','file') == 2, 'UR5_params.m not found.');
assert(exist('UR5_M','file')      == 2, 'UR5_M.m not found.');
assert(exist('UR5_h','file')      == 2, 'UR5_h.m not found.');
assert(exist('UR5_G','file')      == 2, 'UR5_G.m not found.');
assert(exist('UR5_fdyn','file')   == 2, 'UR5_fdyn.m not found.');
assert(exist('UR5_fkine','file')   == 2, 'UR5_fkine.m not found.');

USE_FIG_HELPERS = (exist('figureoptscall','file') == 2) && ...
                  (exist('saveFigureAsPDF','file') == 2);

if USE_FIG_HELPERS
    figureoptscall;
end

outdir = 'ur5_fb_tracking_results';
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

%% ========================= Robot parameters =========================
[~, Pi_true] = UR5_params();
Pi_hat = Pi_true;   % Full-knowledge case
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
%   'sine'     -> continuous sinusoidal motion for the whole simulation
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
% Start from the desired trajectory at t = 0, then add a bias.
[q_ref0, dq_ref0, ~] = desired_trajectory_general(0, traj);

q0_bias  = [ 0.15; -0.10;  0.08;  0.12; -0.06;  0.10];
dq0_bias = zeros(n,1);   % or set a velocity mismatch too

q0  = q_ref0  + q0_bias;
dq0 = dq_ref0 + dq0_bias;

%% ========================= Controller gains =========================
Kp_vec = [10; 10;  8; 5; 3; 2];
Kd_vec = [ 4;  4;  3; 2; 1; 1];
Ki_vec = [ 0.6; 1.8; 1.6; 0.7; 0.15; 0.10];

% Variable-PD scheduling coefficients
alphaKp = [8; 8; 6; 4; 2; 1];
alphaKd = [2; 2; 1.5; 1; 0.5; 0.3];

settling_tol_abs = 1e-3;

%% ========================= Controller list =========================
controllers = struct([]);

controllers(1).name  = 'FFW_PD';
controllers(1).label = '$\mathrm{FFW+PD}$';
controllers(1).mode  = 'FFW_PD';

controllers(2).name  = 'FFW_varPD';
controllers(2).label = '$\mathrm{FFW+variable\ PD}$';
controllers(2).mode  = 'FFW_varPD';

controllers(3).name  = 'FBL_PDFFW';
controllers(3).label = '$\mathrm{FBL+[PD+FFW]}$';
controllers(3).mode  = 'FBL_PDFFW';

controllers(4).name  = 'FBL_PIDFFW';
controllers(4).label = '$\mathrm{FBL+[PID+FFW]}$';
controllers(4).mode  = 'FBL_PIDFFW';

nCtrl = numel(controllers);

%% ========================= Run simulations =========================
results = struct([]);

x0 = [q0; dq0; zeros(n,1)];   % x = [q; dq; eint]

for c = 1:nCtrl
    fprintf('Running %s ...\n', controllers(c).name);

    rhs = @(t, x) closed_loop_rhs_tracking( ...
        t, x, controllers(c), Pi_true, Pi_hat, ...
        traj, Kp_vec, Kd_vec, Ki_vec, alphaKp, alphaKd);

    sol = ode45(rhs, [0 Tf], x0, ode_opts);

    sim = postprocess_solution_tracking( ...
        sol, t_plot, controllers(c), Pi_true, Pi_hat, ...
        traj, Kp_vec, Kd_vec, Ki_vec, alphaKp, alphaKd);

    metrics = compute_metrics_tracking(sim, settling_tol_abs, traj);

    results(c).name    = controllers(c).name;
    results(c).label   = controllers(c).label;
    results(c).mode    = controllers(c).mode;
    results(c).sim     = sim;
    results(c).metrics = metrics;

    write_timeseries_csv_tracking(results(c), outdir);
end

%% ========================= Summary table =========================
summaryTbl = build_summary_table_tracking(results);

summaryTbl = sortrows(summaryTbl, ...
    {'RMS_PosErr_Norm', 'RMS_VelErr_Norm', 'Control_Energy'}, ...
    {'ascend',          'ascend',          'ascend'});

disp(' ');
disp(summaryTbl);

writetable(summaryTbl, fullfile(outdir, 'summary_fb_tracking.csv'));
save(fullfile(outdir, 'summary_fb_tracking.mat'), ...
    'results', 'summaryTbl', ...
    't_plot', 'Tf', 'traj', ...
    'q0', 'dq0', ...
    'Kp_vec', 'Kd_vec', 'Ki_vec', ...
    'alphaKp', 'alphaKd');

%% ========================= Plot styles =========================
styles = {'-','--','-.',':'};
legend_labels = {results.label};

%% ========================= Figure 1: Overview =========================
fig1 = figure('Color','w', 'Name', 'FB Tracking Overview');
tl1 = tiledlayout(fig1, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

hLeg1 = gobjects(1, nCtrl);

ax1 = nexttile(tl1); hold(ax1, 'on');
for c = 1:nCtrl
    h = plot(ax1, results(c).sim.t, results(c).sim.pos_err_norm, ...
        'LineStyle', styles{c}, 'LineWidth', 1.5);
    hLeg1(c) = h;
end
grid(ax1, 'on');
xlabel(ax1, 'Time [s]');
ylabel(ax1, '$\|e_q\|_2$ [rad]');

ax2 = nexttile(tl1); hold(ax2, 'on');
for c = 1:nCtrl
    plot(ax2, results(c).sim.t, results(c).sim.vel_err_norm, ...
        'LineStyle', styles{c}, 'LineWidth', 1.5);
end
grid(ax2, 'on');
xlabel(ax2, 'Time [s]');
ylabel(ax2, '$\|e_{\dot q}\|_2$ [rad/s]');

ax3 = nexttile(tl1); hold(ax3, 'on');
for c = 1:nCtrl
    plot(ax3, results(c).sim.t, results(c).sim.tau_norm, ...
        'LineStyle', styles{c}, 'LineWidth', 1.5);
end
grid(ax3, 'on');
xlabel(ax3, 'Time [s]');
ylabel(ax3, '$\|\tau\|_2$ [Nm]');

lgd1 = legend(hLeg1, legend_labels, 'Orientation', 'horizontal');
lgd1.Layout.Tile = 'north';
lgd1.NumColumns = 2;

if USE_FIG_HELPERS
    saveFigureAsPDF(fig1, fullfile(outdir, 'fig_overview.pdf'));
end

%% ========================= Figure 2: Joint position errors =========================
fig2 = figure('Color','w', 'Name', 'Joint Position Errors');
tl2 = tiledlayout(fig2, 3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

hLeg2 = gobjects(1, nCtrl);

for j = 1:n
    ax = nexttile(tl2); hold(ax, 'on');
    for c = 1:nCtrl
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
lgd2.NumColumns = 2;

if USE_FIG_HELPERS
    saveFigureAsPDF(fig2, fullfile(outdir, 'fig_joint_position_errors.pdf'));
end

%% ========================= Figure 3: Joint velocity errors =========================
fig3 = figure('Color','w', 'Name', 'Joint Velocity Errors');
tl3 = tiledlayout(fig3, 3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

hLeg3 = gobjects(1, nCtrl);

for j = 1:n
    ax = nexttile(tl3); hold(ax, 'on');
    for c = 1:nCtrl
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
lgd3.NumColumns = 2;

if USE_FIG_HELPERS
    saveFigureAsPDF(fig3, fullfile(outdir, 'fig_joint_velocity_errors.pdf'));
end

%% ========================= Figure 4: Joint positions =========================
fig4 = figure('Color','w', 'Name', 'Joint Positions');
tl4 = tiledlayout(fig4, 3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

hLeg4 = gobjects(1, nCtrl);
hqd4 = gobjects(1,1);

for j = 1:n
    ax = nexttile(tl4); hold(ax, 'on');
    for c = 1:nCtrl
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

hLeg5 = gobjects(1, nCtrl);

for j = 1:n
    ax = nexttile(tl5); hold(ax, 'on');
    for c = 1:nCtrl
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
lgd5.NumColumns = 2;

if USE_FIG_HELPERS
    saveFigureAsPDF(fig5, fullfile(outdir, 'fig_joint_torques.pdf'));
end

%% ========================= Figure 6: Integral torque for FBL+[PID+FFW] =========================
pid_idx = find(strcmp({results.name}, 'FBL_PIDFFW'), 1);

if ~isempty(pid_idx)
    fig6 = figure('Color','w', 'Name', 'FBL PID Integral Torque');
    tl6 = tiledlayout(fig6, 3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    for j = 1:n
        ax = nexttile(tl6); hold(ax, 'on');
        plot(ax, results(pid_idx).sim.t, results(pid_idx).sim.tauI_log(j,:), ...
            'LineWidth', 1.35);
        grid(ax, 'on');
        xlabel(ax, 'Time [s]');
        ylabel(ax, sprintf('$\\tau_{I,%d}$ [Nm]', j));
    end

    if USE_FIG_HELPERS
        saveFigureAsPDF(fig6, fullfile(outdir, 'fig_fbl_pid_integral_torque.pdf'));
    end
end

fprintf('\nSaved results to:\n%s\n', fullfile(pwd, outdir));

%% ========================= Local functions =========================

function dx = closed_loop_rhs_tracking(t, x, controller, Pi_true, Pi_hat, ...
                                       traj, Kp_vec, Kd_vec, Ki_vec, alphaKp, alphaKd)

n = 6;

q    = x(1:n);
dq   = x(n+1:2*n);
eint = x(2*n+1:3*n);

[qd, dqd, ddqd] = desired_trajectory_general(t, traj);

[e, de, tau, deint] = tracking_control_law( ...
    q, dq, eint, qd, dqd, ddqd, controller, Pi_hat, ...
    Kp_vec, Kd_vec, Ki_vec, alphaKp, alphaKd);

ddq = UR5_fdyn(q, dq, tau, Pi_true);

dx = [dq; ddq; deint];
end

function sim = postprocess_solution_tracking(sol, t_plot, controller, Pi_true, Pi_hat, ...
                                             traj, Kp_vec, Kd_vec, Ki_vec, alphaKp, alphaKd)

X = deval(sol, t_plot);
n = 6;
N = numel(t_plot);

q_log    = X(1:n, :);
dq_log   = X(n+1:2*n, :);
eint_log = X(2*n+1:3*n, :);

qd_log   = zeros(n, N);
dqd_log  = zeros(n, N);
ddqd_log = zeros(n, N);
ddq_log  = zeros(n, N);
tau_log  = zeros(n, N);
tauI_log = zeros(n, N);
e_log    = zeros(n, N);
de_log   = zeros(n, N);

for k = 1:N
    t    = t_plot(k);
    q    = q_log(:,k);
    dq   = dq_log(:,k);
    eint = eint_log(:,k);

    [qd, dqd, ddqd] = desired_trajectory_general(t, traj);

    [e, de, tau, ~, tauI] = tracking_control_law( ...
        q, dq, eint, qd, dqd, ddqd, controller, Pi_hat, ...
        Kp_vec, Kd_vec, Ki_vec, alphaKp, alphaKd);

    ddq = UR5_fdyn(q, dq, tau, Pi_true);

    qd_log(:,k)   = qd;
    dqd_log(:,k)  = dqd;
    ddqd_log(:,k) = ddqd;
    ddq_log(:,k)  = ddq;
    tau_log(:,k)  = tau;
    tauI_log(:,k) = tauI;
    e_log(:,k)    = e;
    de_log(:,k)   = de;
end

sim.t            = t_plot(:).';
sim.q_log        = q_log;
sim.dq_log       = dq_log;
sim.ddq_log      = ddq_log;
sim.qd_log       = qd_log;
sim.dqd_log      = dqd_log;
sim.ddqd_log     = ddqd_log;
sim.tau_log      = tau_log;
sim.tauI_log     = tauI_log;
sim.e_log        = e_log;
sim.de_log       = de_log;
sim.pos_err_norm = vecnorm(e_log, 2, 1);
sim.vel_err_norm = vecnorm(de_log, 2, 1);
sim.tau_norm     = vecnorm(tau_log, 2, 1);
sim.tauI_norm    = vecnorm(tauI_log, 2, 1);
end

function [e, de, tau, deint, tauI] = tracking_control_law( ...
    q, dq, eint, qd, dqd, ddqd, controller, Pi_hat, ...
    Kp_vec, Kd_vec, Ki_vec, alphaKp, alphaKd)

e  = qd  - q;
de = dqd - dq;

tauI  = zeros(6,1);
deint = zeros(6,1);

switch controller.mode
    case 'FFW_PD'
        tau_ff = UR5_M(qd, Pi_hat) * ddqd + UR5_h(qd, dqd, Pi_hat) + UR5_G(qd, Pi_hat);
        tau_fb = Kp_vec .* e + Kd_vec .* de;
        tau = tau_ff + tau_fb;

    case 'FFW_varPD'
        Kp_var = Kp_vec + alphaKp .* abs(e);
        Kd_var = Kd_vec + alphaKd .* abs(de);

        tau_ff = UR5_M(qd, Pi_hat) * ddqd + UR5_h(qd, dqd, Pi_hat) + UR5_G(qd, Pi_hat);
        tau_fb = Kp_var .* e + Kd_var .* de;
        tau = tau_ff + tau_fb;

    case 'FBL_PDFFW'
        v = ddqd + Kd_vec .* de + Kp_vec .* e;
        tau = UR5_M(q, Pi_hat) * v + UR5_h(q, dq, Pi_hat) + UR5_G(q, Pi_hat);

    case 'FBL_PIDFFW'
        tauI = Ki_vec .* eint;
        v = ddqd + Kd_vec .* de + Kp_vec .* e + Ki_vec .* eint;
        tau = UR5_M(q, Pi_hat) * v + UR5_h(q, dq, Pi_hat) + UR5_G(q, Pi_hat);
        deint = e;

    otherwise
        error('Unknown controller mode: %s', controller.mode);
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

function metrics = compute_metrics_tracking(sim, settling_tol_abs, traj)

t            = sim.t;
pos_err_norm = sim.pos_err_norm;
vel_err_norm = sim.vel_err_norm;
tau_norm     = sim.tau_norm;
tauI_norm    = sim.tauI_norm;

metrics.initial_pos_err_norm = pos_err_norm(1);
metrics.final_pos_err_norm   = pos_err_norm(end);
metrics.rms_pos_err_norm     = sqrt(trapz(t, pos_err_norm.^2) / (t(end) - t(1)));
metrics.peak_pos_err_norm    = max(pos_err_norm);

metrics.initial_vel_err_norm = vel_err_norm(1);
metrics.final_vel_err_norm   = vel_err_norm(end);
metrics.rms_vel_err_norm     = sqrt(trapz(t, vel_err_norm.^2) / (t(end) - t(1)));
metrics.peak_vel_err_norm    = max(vel_err_norm);

metrics.final_tau_norm       = tau_norm(end);
metrics.rms_tau_norm         = sqrt(trapz(t, tau_norm.^2) / (t(end) - t(1)));
metrics.peak_tau_norm        = max(tau_norm);

metrics.final_tauI_norm      = tauI_norm(end);
metrics.peak_tauI_norm       = max(tauI_norm);

metrics.control_energy       = trapz(t, tau_norm.^2);

if strcmpi(traj.mode, 'quintic')
    settle_thresh = max(0.02 * pos_err_norm(1), settling_tol_abs);
    metrics.settling_threshold = settle_thresh;

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
    metrics.settling_threshold = NaN;
    metrics.settling_time = NaN;
end
end

function tbl = build_summary_table_tracking(results)

nCtrl = numel(results);

Controller        = cell(nCtrl,1);
Label             = cell(nCtrl,1);
Final_PosErr_Norm = zeros(nCtrl,1);
RMS_PosErr_Norm   = zeros(nCtrl,1);
Peak_PosErr_Norm  = zeros(nCtrl,1);
Final_VelErr_Norm = zeros(nCtrl,1);
RMS_VelErr_Norm   = zeros(nCtrl,1);
Peak_VelErr_Norm  = zeros(nCtrl,1);
Settling_Time_s   = NaN(nCtrl,1);
Final_Tau_Norm    = zeros(nCtrl,1);
RMS_Tau_Norm      = zeros(nCtrl,1);
Peak_Tau_Norm     = zeros(nCtrl,1);
Final_TauI_Norm   = zeros(nCtrl,1);
Peak_TauI_Norm    = zeros(nCtrl,1);
Control_Energy    = zeros(nCtrl,1);

for c = 1:nCtrl
    Controller{c}         = results(c).name;
    Label{c}              = results(c).label;
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
    Final_TauI_Norm(c)    = results(c).metrics.final_tauI_norm;
    Peak_TauI_Norm(c)     = results(c).metrics.peak_tauI_norm;
    Control_Energy(c)     = results(c).metrics.control_energy;
end

tbl = table(Controller, Label, ...
    Final_PosErr_Norm, RMS_PosErr_Norm, Peak_PosErr_Norm, ...
    Final_VelErr_Norm, RMS_VelErr_Norm, Peak_VelErr_Norm, ...
    Settling_Time_s, ...
    Final_Tau_Norm, RMS_Tau_Norm, Peak_Tau_Norm, ...
    Final_TauI_Norm, Peak_TauI_Norm, ...
    Control_Energy);
end

function write_timeseries_csv_tracking(result, outdir)

t            = result.sim.t(:);
q_log        = result.sim.q_log.';
dq_log       = result.sim.dq_log.';
ddq_log      = result.sim.ddq_log.';
qd_log       = result.sim.qd_log.';
dqd_log      = result.sim.dqd_log.';
ddqd_log     = result.sim.ddqd_log.';
tau_log      = result.sim.tau_log.';
tauI_log     = result.sim.tauI_log.';
e_log        = result.sim.e_log.';
de_log       = result.sim.de_log.';
pos_err_norm = result.sim.pos_err_norm(:);
vel_err_norm = result.sim.vel_err_norm(:);
tau_norm     = result.sim.tau_norm(:);
tauI_norm    = result.sim.tauI_norm(:);

varNames = [{'time'}, ...
    arrayfun(@(i) sprintf('q%d', i),     1:6, 'UniformOutput', false), ...
    arrayfun(@(i) sprintf('dq%d', i),    1:6, 'UniformOutput', false), ...
    arrayfun(@(i) sprintf('ddq%d', i),   1:6, 'UniformOutput', false), ...
    arrayfun(@(i) sprintf('qd%d', i),    1:6, 'UniformOutput', false), ...
    arrayfun(@(i) sprintf('dqd%d', i),   1:6, 'UniformOutput', false), ...
    arrayfun(@(i) sprintf('ddqd%d', i),  1:6, 'UniformOutput', false), ...
    arrayfun(@(i) sprintf('tau%d', i),   1:6, 'UniformOutput', false), ...
    arrayfun(@(i) sprintf('tauI%d', i),  1:6, 'UniformOutput', false), ...
    arrayfun(@(i) sprintf('eq%d', i),    1:6, 'UniformOutput', false), ...
    arrayfun(@(i) sprintf('edq%d', i),   1:6, 'UniformOutput', false), ...
    {'pos_err_norm', 'vel_err_norm', 'tau_norm', 'tauI_norm'}];

data = [t, q_log, dq_log, ddq_log, qd_log, dqd_log, ddqd_log, ...
        tau_log, tauI_log, e_log, de_log, ...
        pos_err_norm, vel_err_norm, tau_norm, tauI_norm];

T = array2table(data, 'VariableNames', varNames);
writetable(T, fullfile(outdir, ['timeseries_' result.name '.csv']));
end