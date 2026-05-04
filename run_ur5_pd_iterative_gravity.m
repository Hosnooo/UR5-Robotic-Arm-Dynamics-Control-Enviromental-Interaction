clc; clear; close all;
addpath('maple_gen')
addpath('utils')

%% ========================= Setup =========================
assert(exist('UR5_params','file') == 2, 'UR5_params.m not found.');
assert(exist('UR5_G','file')      == 2, 'UR5_G.m not found.');
assert(exist('UR5_fdyn','file')   == 2, 'UR5_fdyn.m not found.');

USE_FIG_HELPERS = (exist('figureoptscall','file') == 2) && ...
                  (exist('saveFigureAsPDF','file') == 2);

if USE_FIG_HELPERS
    figureoptscall;
end

outdir = 'ur5_pd_iterative_gravity_results';
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

%% ========================= Initial and desired states =========================
q0  = [ 0.3; -0.8;  0.6; -0.4;  0.2; -0.1];
dq0 = zeros(n,1);

qd_des  = [0; -0.5; 0.8; 0; 0.3; 0];
dqd_des = zeros(n,1);

%% ========================= Two-case ILC settings =========================
N_trials = 10;

Kp_base = [10; 10; 8.0; 5.0; 3.0; 2.0];
Kd_vec  = [ 4;  4; 3.0; 2.0; 1.0; 1.0];

% Case 1: intended to converge
Gamma_good = diag([0.10 0.10 0.08 0.06 0.04 0.03]);
beta_good = 0.98;
delta_ghat_max_good = [4.0; 4.0; 3.0; 2.0; 1.2; 0.8];

% Case 2: intended to improve initially but lose convergence
Gamma_bad = diag([0.22 0.22 0.18 0.14 0.10 0.08]);
beta_bad = 0.995;
delta_ghat_max_bad = [8.0; 8.0; 6.0; 4.0; 2.5; 1.5];

ghat_limit = [80; 120; 80; 40; 20; 10];
learn_start_time = 10.0;

cases = struct([]);

cases(1).name            = 'ILC_convergent';
cases(1).label           = '$\mathrm{ILC}$ convergent';
cases(1).Kp              = Kp_base;
cases(1).Gamma           = Gamma_good;
cases(1).beta            = beta_good;
cases(1).delta_ghat_max  = delta_ghat_max_good;

cases(2).name            = 'ILC_aggressive';
cases(2).label           = '$\mathrm{ILC}$ aggressive';
cases(2).Kp              = Kp_base;
cases(2).Gamma           = Gamma_bad;
cases(2).beta            = beta_bad;
cases(2).delta_ghat_max  = delta_ghat_max_bad;

nCases = numel(cases);

%% ========================= Run trial batches =========================
all_results = struct([]);

for c = 1:nCases
    fprintf('Running case: %s\n', cases(c).name);

    ghat_k = zeros(n,1);
    case_results = struct([]);

    for k = 1:N_trials
        fprintf('  Trial %d / %d\n', k, N_trials);

        x0 = [q0; dq0];

        rhs = @(t, x) closed_loop_rhs( ...
            t, x, Pi, qd_des, dqd_des, ...
            cases(c).Kp, Kd_vec, ghat_k);

        sol = ode45(rhs, [0 Tf], x0, ode_opts);

        sim = postprocess_solution( ...
            sol, t_plot, Pi, qd_des, dqd_des, ...
            cases(c).Kp, Kd_vec, ghat_k, learn_start_time);

        metrics = compute_metrics(sim, settling_tol_abs(err_norm0_from_state(q0, qd_des)));

        case_results(k).trial   = k;
        case_results(k).ghat    = ghat_k;
        case_results(k).sim     = sim;
        case_results(k).metrics = metrics;

        ghat_k = update_gravity_estimate( ...
            ghat_k, sim.e_bias_integral, ...
            cases(c).Gamma, cases(c).beta, ...
            cases(c).delta_ghat_max, ghat_limit);
    end

    all_results(c).name    = cases(c).name;
    all_results(c).label   = cases(c).label;
    all_results(c).Kp      = cases(c).Kp;
    all_results(c).trials  = case_results;
end

%% ========================= Summary table =========================
trialTbl = build_trial_table(all_results);
disp(' ');
disp(trialTbl);

writetable(trialTbl, fullfile(outdir, 'trial_metrics.csv'));
save(fullfile(outdir, 'trial_metrics.mat'), ...
    'all_results', 'trialTbl', ...
    'Tf', 't_plot', 'q0', 'dq0', 'qd_des', 'dqd_des', ...
    'Kp_base', 'Kd_vec', ...
    'Gamma_good', 'Gamma_bad', ...
    'beta_good', 'beta_bad', ...
    'delta_ghat_max_good', 'delta_ghat_max_bad', ...
    'ghat_limit', 'learn_start_time', 'N_trials');

%% ========================= Plot styles =========================
styles = {'-','--'};

%% ========================= Figure 1: Final error vs trial =========================
fig1 = figure('Color','w', 'Name', 'Final Error vs Trial');
ax1 = axes(fig1); hold(ax1, 'on');

for c = 1:nCases
    vals = arrayfun(@(s) s.metrics.final_error_norm, all_results(c).trials);
    plot(ax1, 1:N_trials, vals, 'LineStyle', styles{c}, 'LineWidth', 1.8);
end

grid(ax1, 'on');
xlabel(ax1, 'Trial');
ylabel(ax1, 'Final $\|e\|_2$ [rad]');
legend({all_results.label}, 'Location', 'northeast');

if USE_FIG_HELPERS
    saveFigureAsPDF(fig1, fullfile(outdir, 'fig_final_error_vs_trial.pdf'));
end

%% ========================= Figure 2: RMS error vs trial =========================
fig2 = figure('Color','w', 'Name', 'RMS Error vs Trial');
ax2 = axes(fig2); hold(ax2, 'on');

for c = 1:nCases
    vals = arrayfun(@(s) s.metrics.rms_error_norm, all_results(c).trials);
    plot(ax2, 1:N_trials, vals, 'LineStyle', styles{c}, 'LineWidth', 1.8);
end

grid(ax2, 'on');
xlabel(ax2, 'Trial');
ylabel(ax2, 'RMS $\|e\|_2$ [rad]');
legend({all_results.label}, 'Location', 'northeast');

if USE_FIG_HELPERS
    saveFigureAsPDF(fig2, fullfile(outdir, 'fig_rms_error_vs_trial.pdf'));
end

%% ========================= Figure 3: Learned gravity estimate =========================
fig3 = figure('Color','w', 'Name', 'Learned Gravity Estimate');
tl3 = tiledlayout(fig3, 3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for j = 1:n
    ax = nexttile(tl3); hold(ax, 'on');
    for c = 1:nCases
        Ghist = zeros(N_trials,1);
        for k = 1:N_trials
            Ghist(k) = all_results(c).trials(k).ghat(j);
        end
        plot(ax, 1:N_trials, Ghist, 'LineStyle', styles{c}, 'LineWidth', 1.5);
    end
    grid(ax, 'on');
    xlabel(ax, 'Trial');
    ylabel(ax, sprintf('$\\hat g_%d$ [Nm]', j));
end

lgd3 = legend({all_results.label}, 'Orientation', 'horizontal');
lgd3.Layout.Tile = 'north';
lgd3.NumColumns = 2;

if USE_FIG_HELPERS
    saveFigureAsPDF(fig3, fullfile(outdir, 'fig_learned_gravity_vs_trial.pdf'));
end

%% ========================= Figure 4: Last-trial joint errors =========================
fig4 = figure('Color','w', 'Name', 'Last Trial Error Histories');
tl4 = tiledlayout(fig4, 3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

hLeg4 = gobjects(1, nCases);

for j = 1:n
    ax = nexttile(tl4); hold(ax, 'on');
    for c = 1:nCases
        sim_last = all_results(c).trials(end).sim;
        h = plot(ax, sim_last.t, sim_last.e_log(j,:), ...
            'LineStyle', styles{c}, 'LineWidth', 1.4);
        if j == 1
            hLeg4(c) = h;
        end
    end
    yline(ax, 0, 'k:', 'LineWidth', 1.0);
    grid(ax, 'on');
    xlabel(ax, 'Time [s]');
    ylabel(ax, sprintf('$e_%d$ [rad]', j));
end

lgd4 = legend(hLeg4, {all_results.label}, 'Orientation', 'horizontal');
lgd4.Layout.Tile = 'north';
lgd4.NumColumns = 2;

if USE_FIG_HELPERS
    saveFigureAsPDF(fig4, fullfile(outdir, 'fig_last_trial_joint_errors.pdf'));
end

fprintf('\nSaved results to:\n%s\n', fullfile(pwd, outdir));

%% ========================= Local functions =========================
function dx = closed_loop_rhs(~, x, Pi, qd_des, dqd_des, Kp_vec, Kd_vec, ghat)
n = 6;

q  = x(1:n);
dq = x(n+1:2*n);

[e, de, tau] = control_law(q, dq, qd_des, dqd_des, Kp_vec, Kd_vec, ghat);
ddq = UR5_fdyn(q, dq, tau, Pi);

dx = [dq; ddq];
end

function sim = postprocess_solution(sol, t_plot, Pi, qd_des, dqd_des, Kp_vec, Kd_vec, ghat, learn_start_time)
X = deval(sol, t_plot);
n = 6;
N = numel(t_plot);

q_log   = X(1:n, :);
dq_log  = X(n+1:2*n, :);

ddq_log = zeros(n, N);
tau_log = zeros(n, N);
e_log   = zeros(n, N);
de_log  = zeros(n, N);

for k = 1:N
    q  = q_log(:,k);
    dq = dq_log(:,k);

    [e, de, tau] = control_law(q, dq, qd_des, dqd_des, Kp_vec, Kd_vec, ghat);
    ddq = UR5_fdyn(q, dq, tau, Pi);

    e_log(:,k)   = e;
    de_log(:,k)  = de;
    tau_log(:,k) = tau;
    ddq_log(:,k) = ddq;
end

sim.t       = t_plot(:).';
sim.q_log   = q_log;
sim.dq_log  = dq_log;
sim.ddq_log = ddq_log;
sim.tau_log = tau_log;
sim.e_log   = e_log;
sim.de_log  = de_log;

sim.err_norm = vecnorm(e_log, 2, 1);
sim.tau_norm = vecnorm(tau_log, 2, 1);

idx_learn = sim.t >= learn_start_time;
if ~any(idx_learn)
    idx_learn = true(size(sim.t));
end
sim.e_bias_integral = trapz(sim.t(idx_learn), e_log(:, idx_learn), 2);
end

function [e, de, tau] = control_law(q, dq, qd_des, dqd_des, Kp_vec, Kd_vec, ghat)
e  = qd_des  - q;
de = dqd_des - dq;

tau = Kp_vec .* e + Kd_vec .* de + ghat;
end

function metrics = compute_metrics(sim, settling_tol_abs)
t        = sim.t;
err_norm = sim.err_norm;
tau_norm = sim.tau_norm;

metrics.initial_error_norm = err_norm(1);
metrics.final_error_norm   = err_norm(end);
metrics.rms_error_norm     = sqrt(trapz(t, err_norm.^2) / (t(end) - t(1)));
metrics.peak_error_norm    = max(err_norm);

metrics.final_tau_norm     = tau_norm(end);
metrics.rms_tau_norm       = sqrt(trapz(t, tau_norm.^2) / (t(end) - t(1)));
metrics.peak_tau_norm      = max(tau_norm);

settle_thresh = max(0.02 * err_norm(1), settling_tol_abs);

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

function tbl = build_trial_table(all_results)
Case            = {};
Trial           = [];
Final_Error     = [];
RMS_Error       = [];
Peak_Error      = [];
Settling_Time_s = [];
Peak_Tau        = [];

for c = 1:numel(all_results)
    for k = 1:numel(all_results(c).trials)
        Case{end+1,1}            = all_results(c).name;
        Trial(end+1,1)           = all_results(c).trials(k).trial;
        Final_Error(end+1,1)     = all_results(c).trials(k).metrics.final_error_norm;
        RMS_Error(end+1,1)       = all_results(c).trials(k).metrics.rms_error_norm;
        Peak_Error(end+1,1)      = all_results(c).trials(k).metrics.peak_error_norm;
        Settling_Time_s(end+1,1) = all_results(c).trials(k).metrics.settling_time;
        Peak_Tau(end+1,1)        = all_results(c).trials(k).metrics.peak_tau_norm;
    end
end

tbl = table(Case, Trial, Final_Error, RMS_Error, Peak_Error, Settling_Time_s, Peak_Tau);
end

function ghat_next = update_gravity_estimate(ghat_k, e_bias_integral, Gamma, beta, delta_lim, ghat_lim)
delta = Gamma * e_bias_integral;
delta = min(max(delta, -delta_lim), delta_lim);

ghat_next = beta * ghat_k + delta;
ghat_next = min(max(ghat_next, -ghat_lim), ghat_lim);
end

function y = err_norm0_from_state(q0, qd_des)
y = norm(qd_des - q0, 2);
end

function y = settling_tol_abs(err0)
y = max(1e-3, 0.02 * err0);
end
