clc; clear; close all;
addpath('maple_gen')
addpath('utils')

%% ========================= Setup =========================
assert(exist('UR5_params','file') == 2, 'UR5_params.m not found.');
assert(exist('UR5_fdyn','file')   == 2, 'UR5_fdyn.m not found.');
assert(exist('UR5_fkine','file')  == 2, 'UR5_fkine.m not found.');

USE_FIG_HELPERS = (exist('figureoptscall','file') == 2) && ...
                  (exist('saveFigureAsPDF','file') == 2);

if USE_FIG_HELPERS
    figureoptscall;
end

outdir = 'ur5_pid_sweep_results';
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

%% ========================= Robot parameters =========================
[~, Pi] = UR5_params();
n = 6;

%% ========================= Simulation settings =========================
Tf = 30.0;
t_plot = linspace(0, Tf, 3001);

ode_opts = odeset( ...
    'RelTol', 1e-6, ...
    'AbsTol', 1e-8, ...
    'MaxStep', 1e-2);

%% ========================= Initial conditions =========================
q0  = [ 0.3; -0.8;  0.6; -0.4;  0.2; -0.1];
dq0 = zeros(n,1);

%% ========================= Desired regulation target =========================
qd_des  = [0; -0.5; 0.8; 0; 0.3; 0];
dqd_des = zeros(n,1);

%% ========================= Fixed PD gains =========================
Kp_vec = [10; 10;  8; 5; 3; 2];
Kd_vec = [ 4;  4;  3; 2; 1; 1];

%% ========================= Base integral-gain pattern =========================
Ki_base_vec = [ 3;  3;  2.5; 2.2; 1; 0.8];

%% ========================= Integral-gain sweep =========================
Ki_scales = [0.5, 1.0, 1.5, 2.0, 2.2];
nCases = numel(Ki_scales);

%% ========================= Metrics settings =========================
settling_tol_abs = 1e-3;

%% ========================= Run sweep =========================
results = struct([]);
x0 = [q0; dq0; zeros(n,1)];   % x = [q; dq; eint]

for c = 1:nCases
    scale = Ki_scales(c);
    Ki_vec = scale * Ki_base_vec;

    fprintf('Running PID with Ki scale = %.2f ...\n', scale);

    rhs = @(t, x) closed_loop_rhs_pid( ...
        t, x, Pi, qd_des, dqd_des, Kp_vec, Kd_vec, Ki_vec);

    sol = ode45(rhs, [0 Tf], x0, ode_opts);

    sim = postprocess_solution_pid( ...
        sol, t_plot, Pi, qd_des, dqd_des, Kp_vec, Kd_vec, Ki_vec);

    metrics = compute_metrics_pid(sim, settling_tol_abs);

    results(c).name    = sprintf('PID_KiScale_%0.2f', scale);
    results(c).label   = sprintf('$K_i = %.1fK_{i,\\mathrm{base}}$', scale);
    results(c).scale   = scale;
    results(c).Ki_vec  = Ki_vec;
    results(c).sim     = sim;
    results(c).metrics = metrics;

    write_timeseries_csv_pid(results(c), outdir);
end

%% ========================= Summary table =========================
summaryTbl = build_summary_table_pid(results);

summaryTbl = sortrows(summaryTbl, ...
    {'RMS_Error_Norm', 'Settling_Time_s', 'Control_Energy'}, ...
    {'ascend',         'ascend',          'ascend'});

disp(' ');
disp(summaryTbl);

writetable(summaryTbl, fullfile(outdir, 'summary_pid_igain_sweep.csv'));
save(fullfile(outdir, 'summary_pid_igain_sweep.mat'), ...
    'results', 'summaryTbl', ...
    't_plot', 'Tf', ...
    'q0', 'dq0', 'qd_des', 'dqd_des', ...
    'Kp_vec', 'Kd_vec', 'Ki_base_vec', 'Ki_scales');

%% ========================= Plot styles =========================
styles = {'-','--','-.',':','-'};
legend_labels = {results.label};

%% ========================= Figure 1: Overview =========================
fig1 = figure('Color','w', 'Name', 'PID I-Gain Sweep Overview');
tl1 = tiledlayout(fig1, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

hLeg1 = gobjects(1, nCases);

ax1 = nexttile(tl1); hold(ax1, 'on');
for c = 1:nCases
    h = plot(ax1, results(c).sim.t, results(c).sim.err_norm, ...
        'LineStyle', styles{c}, 'LineWidth', 1.5);
    hLeg1(c) = h;
end
grid(ax1, 'on');
xlabel(ax1, 'Time [s]');
ylabel(ax1, '$\|e\|_2$ [rad]');

ax2 = nexttile(tl1); hold(ax2, 'on');
for c = 1:nCases
    plot(ax2, results(c).sim.t, results(c).sim.tau_norm, ...
        'LineStyle', styles{c}, 'LineWidth', 1.5);
end
grid(ax2, 'on');
xlabel(ax2, 'Time [s]');
ylabel(ax2, '$\|\tau\|_2$ [Nm]');

lgd1 = legend(hLeg1, legend_labels, 'Orientation', 'horizontal');
lgd1.Layout.Tile = 'north';
lgd1.NumColumns = 3;

if USE_FIG_HELPERS
    saveFigureAsPDF(fig1, fullfile(outdir, 'fig_pid_igain_overview.pdf'));
end

%% ========================= Figure 2: Joint errors =========================
fig2 = figure('Color','w', 'Name', 'PID I-Gain Sweep Joint Errors');
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
    ylabel(ax, sprintf('$e_%d$ [rad]', j));
end

lgd2 = legend(hLeg2, legend_labels, 'Orientation', 'horizontal');
lgd2.Layout.Tile = 'north';
lgd2.NumColumns = 3;

if USE_FIG_HELPERS
    saveFigureAsPDF(fig2, fullfile(outdir, 'fig_pid_igain_joint_errors.pdf'));
end

%% ========================= Figure 3: Joint positions =========================
fig3 = figure('Color','w', 'Name', 'PID I-Gain Sweep Joint Positions');
tl3 = tiledlayout(fig3, 3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

hLeg3 = gobjects(1, nCases);
hqd3 = gobjects(1,1);

for j = 1:n
    ax = nexttile(tl3); hold(ax, 'on');
    for c = 1:nCases
        h = plot(ax, results(c).sim.t, results(c).sim.q_log(j,:), ...
            'LineStyle', styles{c}, 'LineWidth', 1.35);
        if j == 1
            hLeg3(c) = h;
        end
    end
    hqd = yline(ax, qd_des(j), 'k:', 'LineWidth', 1.2);
    if j == 1
        hqd3 = hqd;
    end
    grid(ax, 'on');
    xlabel(ax, 'Time [s]');
    ylabel(ax, sprintf('$q_%d$ [rad]', j));
end

legend_entries3 = [legend_labels, {'$q_d$'}];
lgd3 = legend([hLeg3, hqd3], legend_entries3, 'Orientation', 'horizontal');
lgd3.Layout.Tile = 'north';
lgd3.NumColumns = 3;

if USE_FIG_HELPERS
    saveFigureAsPDF(fig3, fullfile(outdir, 'fig_pid_igain_joint_positions.pdf'));
end

%% ========================= Figure 4: Joint torques =========================
fig4 = figure('Color','w', 'Name', 'PID I-Gain Sweep Joint Torques');
tl4 = tiledlayout(fig4, 3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

hLeg4 = gobjects(1, nCases);

for j = 1:n
    ax = nexttile(tl4); hold(ax, 'on');
    for c = 1:nCases
        h = plot(ax, results(c).sim.t, results(c).sim.tau_log(j,:), ...
            'LineStyle', styles{c}, 'LineWidth', 1.35);
        if j == 1
            hLeg4(c) = h;
        end
    end
    grid(ax, 'on');
    xlabel(ax, 'Time [s]');
    ylabel(ax, sprintf('$\\tau_%d$ [Nm]', j));
end

lgd4 = legend(hLeg4, legend_labels, 'Orientation', 'horizontal');
lgd4.Layout.Tile = 'north';
lgd4.NumColumns = 3;

if USE_FIG_HELPERS
    saveFigureAsPDF(fig4, fullfile(outdir, 'fig_pid_igain_joint_torques.pdf'));
end

%% ========================= Figure 5: Integral torque contribution =========================
fig5 = figure('Color','w', 'Name', 'PID I-Gain Sweep Integral Torque');
tl5 = tiledlayout(fig5, 3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

hLeg5 = gobjects(1, nCases);

for j = 1:n
    ax = nexttile(tl5); hold(ax, 'on');
    for c = 1:nCases
        h = plot(ax, results(c).sim.t, results(c).sim.tauI_log(j,:), ...
            'LineStyle', styles{c}, 'LineWidth', 1.35);
        if j == 1
            hLeg5(c) = h;
        end
    end
    grid(ax, 'on');
    xlabel(ax, 'Time [s]');
    ylabel(ax, sprintf('$\\tau_{I,%d}$ [Nm]', j));
end

lgd5 = legend(hLeg5, legend_labels, 'Orientation', 'horizontal');
lgd5.Layout.Tile = 'north';
lgd5.NumColumns = 3;

if USE_FIG_HELPERS
    saveFigureAsPDF(fig5, fullfile(outdir, 'fig_pid_igain_integral_torque.pdf'));
end

fprintf('\nSaved results to:\n%s\n', fullfile(pwd, outdir));

%% ========================= Local functions =========================

function dx = closed_loop_rhs_pid(~, x, Pi, qd_des, dqd_des, Kp_vec, Kd_vec, Ki_vec)

n = 6;

q    = x(1:n);
dq   = x(n+1:2*n);
eint = x(2*n+1:3*n);

[e, de, tau, deint] = pid_control_law(q, dq, eint, qd_des, dqd_des, Kp_vec, Kd_vec, Ki_vec);
ddq = UR5_fdyn(q, dq, tau, Pi);

dx = [dq; ddq; deint];
end

function sim = postprocess_solution_pid(sol, t_plot, Pi, qd_des, dqd_des, Kp_vec, Kd_vec, Ki_vec)

X = deval(sol, t_plot);
n = 6;
N = numel(t_plot);

q_log    = X(1:n, :);
dq_log   = X(n+1:2*n, :);
eint_log = X(2*n+1:3*n, :);

ddq_log  = zeros(n, N);
tau_log  = zeros(n, N);
tauI_log = zeros(n, N);
e_log    = zeros(n, N);
de_log   = zeros(n, N);

for k = 1:N
    q    = q_log(:,k);
    dq   = dq_log(:,k);
    eint = eint_log(:,k);

    [e, de, tau, ~, tauI] = pid_control_law(q, dq, eint, qd_des, dqd_des, Kp_vec, Kd_vec, Ki_vec);
    ddq = UR5_fdyn(q, dq, tau, Pi);

    e_log(:,k)    = e;
    de_log(:,k)   = de;
    tau_log(:,k)  = tau;
    tauI_log(:,k) = tauI;
    ddq_log(:,k)  = ddq;
end

sim.t         = t_plot(:).';
sim.q_log     = q_log;
sim.dq_log    = dq_log;
sim.ddq_log   = ddq_log;
sim.tau_log   = tau_log;
sim.tauI_log  = tauI_log;
sim.e_log     = e_log;
sim.de_log    = de_log;
sim.err_norm  = vecnorm(e_log, 2, 1);
sim.tau_norm  = vecnorm(tau_log, 2, 1);
sim.tauI_norm = vecnorm(tauI_log, 2, 1);
end

function [e, de, tau, deint, tauI] = pid_control_law(q, dq, eint, qd_des, dqd_des, Kp_vec, Kd_vec, Ki_vec)

e  = qd_des  - q;
de = dqd_des - dq;

tauPD = Kp_vec .* e + Kd_vec .* de;
tauI  = Ki_vec .* eint;
tau   = tauPD + tauI;

deint = e;
end

function metrics = compute_metrics_pid(sim, settling_tol_abs)

t         = sim.t;
err_norm  = sim.err_norm;
tau_norm  = sim.tau_norm;
tauI_norm = sim.tauI_norm;

metrics.initial_error_norm = err_norm(1);
metrics.final_error_norm   = err_norm(end);
metrics.rms_error_norm     = sqrt(trapz(t, err_norm.^2) / (t(end) - t(1)));
metrics.peak_error_norm    = max(err_norm);

metrics.final_tau_norm     = tau_norm(end);
metrics.rms_tau_norm       = sqrt(trapz(t, tau_norm.^2) / (t(end) - t(1)));
metrics.peak_tau_norm      = max(tau_norm);

metrics.final_tauI_norm    = tauI_norm(end);
metrics.peak_tauI_norm     = max(tauI_norm);

metrics.control_energy     = trapz(t, tau_norm.^2);

settle_thresh = max(0.02 * err_norm(1), settling_tol_abs);
metrics.settling_threshold = settle_thresh;

settle_idx = NaN;
for k = 1:numel(t)
    if all(err_norm(k:end) <= settle_thresh)
        settle_idx = k;
        break;
    end
end

if isnan(settle_idx)
    metrics.settling_time = NaN;
else
    metrics.settling_time = t(settle_idx);
end
end

function tbl = build_summary_table_pid(results)

nCases = numel(results);

CaseName          = cell(nCases,1);
Label             = cell(nCases,1);
KiScale           = zeros(nCases,1);
Final_Error_Norm  = zeros(nCases,1);
RMS_Error_Norm    = zeros(nCases,1);
Peak_Error_Norm   = zeros(nCases,1);
Settling_Time_s   = NaN(nCases,1);
Final_Tau_Norm    = zeros(nCases,1);
RMS_Tau_Norm      = zeros(nCases,1);
Peak_Tau_Norm     = zeros(nCases,1);
Final_TauI_Norm   = zeros(nCases,1);
Peak_TauI_Norm    = zeros(nCases,1);
Control_Energy    = zeros(nCases,1);

for c = 1:nCases
    CaseName{c}         = results(c).name;
    Label{c}            = results(c).label;
    KiScale(c)          = results(c).scale;
    Final_Error_Norm(c) = results(c).metrics.final_error_norm;
    RMS_Error_Norm(c)   = results(c).metrics.rms_error_norm;
    Peak_Error_Norm(c)  = results(c).metrics.peak_error_norm;
    Settling_Time_s(c)  = results(c).metrics.settling_time;
    Final_Tau_Norm(c)   = results(c).metrics.final_tau_norm;
    RMS_Tau_Norm(c)     = results(c).metrics.rms_tau_norm;
    Peak_Tau_Norm(c)    = results(c).metrics.peak_tau_norm;
    Final_TauI_Norm(c)  = results(c).metrics.final_tauI_norm;
    Peak_TauI_Norm(c)   = results(c).metrics.peak_tauI_norm;
    Control_Energy(c)   = results(c).metrics.control_energy;
end

tbl = table(CaseName, Label, KiScale, ...
    Final_Error_Norm, RMS_Error_Norm, Peak_Error_Norm, Settling_Time_s, ...
    Final_Tau_Norm, RMS_Tau_Norm, Peak_Tau_Norm, ...
    Final_TauI_Norm, Peak_TauI_Norm, ...
    Control_Energy);
end

function write_timeseries_csv_pid(result, outdir)

t         = result.sim.t(:);
q_log     = result.sim.q_log.';
dq_log    = result.sim.dq_log.';
ddq_log   = result.sim.ddq_log.';
tau_log   = result.sim.tau_log.';
tauI_log  = result.sim.tauI_log.';
e_log     = result.sim.e_log.';
de_log    = result.sim.de_log.';
err_norm  = result.sim.err_norm(:);
tau_norm  = result.sim.tau_norm(:);
tauI_norm = result.sim.tauI_norm(:);

varNames = [{'time'}, ...
    arrayfun(@(i) sprintf('q%d', i),    1:6, 'UniformOutput', false), ...
    arrayfun(@(i) sprintf('dq%d', i),   1:6, 'UniformOutput', false), ...
    arrayfun(@(i) sprintf('ddq%d', i),  1:6, 'UniformOutput', false), ...
    arrayfun(@(i) sprintf('tau%d', i),  1:6, 'UniformOutput', false), ...
    arrayfun(@(i) sprintf('tauI%d', i), 1:6, 'UniformOutput', false), ...
    arrayfun(@(i) sprintf('e%d', i),    1:6, 'UniformOutput', false), ...
    arrayfun(@(i) sprintf('de%d', i),   1:6, 'UniformOutput', false), ...
    {'err_norm', 'tau_norm', 'tauI_norm'}];

data = [t, q_log, dq_log, ddq_log, tau_log, tauI_log, e_log, de_log, ...
        err_norm, tau_norm, tauI_norm];

T = array2table(data, 'VariableNames', varNames);
writetable(T, fullfile(outdir, ['timeseries_' result.name '.csv']));
end