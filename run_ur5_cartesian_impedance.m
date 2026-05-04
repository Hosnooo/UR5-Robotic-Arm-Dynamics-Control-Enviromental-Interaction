clc; clear; close all;
addpath('maple_gen')
addpath('utils')

%% ========================= Setup =========================
assert(exist('UR5_params','file')              == 2, 'UR5_params.m not found.');
assert(exist('UR5_M','file')                   == 2, 'UR5_M.m not found.');
assert(exist('UR5_h','file')                   == 2, 'UR5_h.m not found.');
assert(exist('UR5_G','file')                   == 2, 'UR5_G.m not found.');
assert(exist('UR5_fkine','file')               == 2, 'UR5_fkine.m not found.');
assert(exist('UR5_jacobian_geometric','file')  == 2, 'UR5_jacobian_geometric.m not found.');

USE_FIG_HELPERS = (exist('figureoptscall','file') == 2) && ...
                  (exist('saveFigureAsPDF','file') == 2);

if USE_FIG_HELPERS
    figureoptscall;
end

outdir = 'ur5_cartesian_impedance_results';
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

%% ========================= Robot parameters =========================
[~, Pi] = UR5_params();
n = 6;

%% ========================= Simulation settings =========================
Tf     = 18.0;
t_plot = linspace(0, Tf, 3601);

ode_opts = odeset( ...
    'RelTol', 1e-6, ...
    'AbsTol', 1e-8, ...
    'MaxStep', 1e-2);

%% ========================= Initial state =========================
q0  = [ 0.35; -1.05; 1.20; -1.30; 1.00; 0.20 ];
dq0 = zeros(n,1);

[p0, ~] = ee_pos_jac(q0, Pi);

%% ========================= Cartesian trajectory =========================
traj = struct();
traj.mode  = 'quintic3d';
traj.Tmove = 12.0;
traj.p0    = p0;
traj.pf    = p0 + [0.06; 0.05; 0.03];

traj.pd_ref = zeros(3, numel(t_plot));
for k = 1:numel(t_plot)
    [traj.pd_ref(:,k), ~, ~] = desired_cartesian_path(t_plot(k), traj);
end

%% ========================= External wrench =========================
extw = struct();
extw.mode   = 'half_sine';
extw.t_on   = 6.0;
extw.t_off  = 10.0;
extw.f_axis = [0; 1; 0];
extw.f_peak = 18.0;

%% ========================= Controller settings =========================
ctrl = struct();
ctrl.Kq_damp = diag([3.0, 3.0, 2.0, 1.5, 1.0, 0.8]);
ctrl.tau_max = [180; 180; 140; 90; 60; 40];

%% ========================= Numerical settings =========================
con = struct();
con.fd_dt    = 1e-6;
con.jac_damp = 1e-4;

%% ========================= Impedance cases =========================
cases = struct([]);

cases(1).name        = 'Nominal_NoForce';
cases(1).label       = 'Nominal, no force';
cases(1).apply_force = false;
cases(1).Md          = diag([2.5, 2.5, 2.5]);
cases(1).Dd          = diag([35, 35, 35]);
cases(1).Kd          = diag([120, 120, 120]);

cases(2).name        = 'Nominal_WithForce';
cases(2).label       = 'Nominal, with force';
cases(2).apply_force = true;
cases(2).Md          = diag([2.5, 2.5, 2.5]);
cases(2).Dd          = diag([35, 35, 35]);
cases(2).Kd          = diag([120, 120, 120]);

cases(3).name        = 'Soft_WithForce';
cases(3).label       = 'Soft, with force';
cases(3).apply_force = true;
cases(3).Md          = diag([3.0, 3.0, 3.0]);
cases(3).Dd          = diag([24, 24, 24]);
cases(3).Kd          = diag([70, 70, 70]);

cases(4).name        = 'Stiff_WithForce';
cases(4).label       = 'Stiff, with force';
cases(4).apply_force = true;
cases(4).Md          = diag([2.0, 2.0, 2.0]);
cases(4).Dd          = diag([50, 50, 50]);
cases(4).Kd          = diag([220, 220, 220]);

nCases = numel(cases);

%% ========================= Run simulations =========================
results = struct([]);
x0 = [q0; dq0];

for c = 1:nCases
    fprintf('Running %s ...\n', cases(c).name);

    rhs = @(t, x) closed_loop_rhs_impedance(t, x, Pi, traj, extw, ctrl, cases(c), con);
    sol = ode45(rhs, [0 Tf], x0, ode_opts);

    sim = postprocess_solution_impedance(sol, t_plot, Pi, traj, extw, ctrl, cases(c), con);

    results(c).name  = cases(c).name;
    results(c).label = cases(c).label;
    results(c).Md    = cases(c).Md;
    results(c).Dd    = cases(c).Dd;
    results(c).Kd    = cases(c).Kd;
    results(c).sim   = sim;

    write_timeseries_csv_impedance(results(c), outdir);
end

results   = add_impedance_metrics(results);
summaryTbl = build_summary_table_impedance(results);

disp(' ');
disp(summaryTbl);

writetable(summaryTbl, fullfile(outdir, 'summary_impedance.csv'));
save(fullfile(outdir, 'summary_impedance.mat'), ...
    'results', 'summaryTbl', 'traj', 'extw', 'ctrl', 'con', ...
    'q0', 'dq0', 'Tf', 't_plot');

%% ========================= Plot styles =========================
styles = {'-','--','-.',':'};
legend_labels = {results.label};

%% ========================= Figure 1: Overview =========================
fig1 = figure('Color','w', 'Name', 'Cartesian Impedance Overview');
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
ylabel(ax1, '$\|e_p\|_2$ [m]');

ax2 = nexttile(tl1); hold(ax2, 'on');
for c = 1:nCases
    plot(ax2, results(c).sim.t, results(c).sim.fext_norm, ...
        'LineStyle', styles{c}, 'LineWidth', 1.5);
end
grid(ax2, 'on');
xlabel(ax2, 'Time [s]');
ylabel(ax2, '$\|f_{ext}\|_2$ [N]');

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
lgd1.NumColumns = 2;

if USE_FIG_HELPERS
    saveFigureAsPDF(fig1, fullfile(outdir, 'fig_impedance_overview.pdf'));
end

%% ========================= Figure 2: 3D paths =========================
fig2 = figure('Color','w', 'Name', 'Cartesian Impedance 3D Paths');
hold on; grid on; axis equal; view(135, 25);

plot3(traj.pd_ref(1,:), traj.pd_ref(2,:), traj.pd_ref(3,:), ...
    'k:', 'LineWidth', 1.4, 'DisplayName', 'Desired path');

for c = 1:nCases
    plot3(results(c).sim.p_log(1,:), results(c).sim.p_log(2,:), results(c).sim.p_log(3,:), ...
        'LineStyle', styles{c}, 'LineWidth', 1.6, 'DisplayName', results(c).label);
end

plot3(traj.pd_ref(1,1), traj.pd_ref(2,1), traj.pd_ref(3,1), ...
    'go', 'MarkerFaceColor', 'g', 'DisplayName', 'Start');

plot3(traj.pd_ref(1,end), traj.pd_ref(2,end), traj.pd_ref(3,end), ...
    'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'End');

xlabel('x [m]');
ylabel('y [m]');
zlabel('z [m]');
legend('Location', 'best');

if USE_FIG_HELPERS
    saveFigureAsPDF(fig2, fullfile(outdir, 'fig_impedance_3d_paths.pdf'));
end

%% ========================= Figure 3: Cartesian position components =========================
fig3 = figure('Color','w', 'Name', 'Cartesian Position Components');
tl3 = tiledlayout(fig3, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

coord_labels = {'$p_x$ [m]', '$p_y$ [m]', '$p_z$ [m]'};
hLeg3 = gobjects(1, nCases);
hPd3  = gobjects(1,1);

for i = 1:3
    ax = nexttile(tl3); hold(ax, 'on');
    for c = 1:nCases
        h = plot(ax, results(c).sim.t, results(c).sim.p_log(i,:), ...
            'LineStyle', styles{c}, 'LineWidth', 1.35);
        if i == 1
            hLeg3(c) = h;
        end
    end
    hpd = plot(ax, t_plot, traj.pd_ref(i,:), 'k:', 'LineWidth', 1.25);
    if i == 1
        hPd3 = hpd;
    end
    grid(ax, 'on');
    xlabel(ax, 'Time [s]');
    ylabel(ax, coord_labels{i});
end

lgd3 = legend([hLeg3, hPd3], [legend_labels, {'Desired'}], 'Orientation', 'horizontal');
lgd3.Layout.Tile = 'north';
lgd3.NumColumns = 3;

if USE_FIG_HELPERS
    saveFigureAsPDF(fig3, fullfile(outdir, 'fig_impedance_position_components.pdf'));
end

%% ========================= Figure 4: Force-induced deviation =========================
fig4 = figure('Color','w', 'Name', 'Force-Induced Deviation');
tl4 = tiledlayout(fig4, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

force_case_idx = find(contains({results.name}, 'WithForce'));
hLeg4 = gobjects(1, numel(force_case_idx));

ax1 = nexttile(tl4); hold(ax1, 'on');
for k = 1:numel(force_case_idx)
    c = force_case_idx(k);
    h = plot(ax1, results(c).sim.t, results(c).sim.path_dev_vs_nominal, ...
        'LineStyle', styles{c}, 'LineWidth', 1.5);
    hLeg4(k) = h;
end
grid(ax1, 'on');
xlabel(ax1, 'Time [s]');
ylabel(ax1, '$\|p - p_{nom}\|_2$ [m]');

ax2 = nexttile(tl4); hold(ax2, 'on');
for k = 1:numel(force_case_idx)
    c = force_case_idx(k);
    plot(ax2, results(c).sim.t, results(c).sim.pos_err_norm, ...
        'LineStyle', styles{c}, 'LineWidth', 1.5);
end
grid(ax2, 'on');
xlabel(ax2, 'Time [s]');
ylabel(ax2, '$\|e_p\|_2$ [m]');

lgd4 = legend(hLeg4, {results(force_case_idx).label}, 'Orientation', 'horizontal');
lgd4.Layout.Tile = 'north';
lgd4.NumColumns = 3;

if USE_FIG_HELPERS
    saveFigureAsPDF(fig4, fullfile(outdir, 'fig_impedance_force_deviation.pdf'));
end

%% ========================= Figure 5: Joint torques (nominal with force) =========================
idx_nom_force = find(strcmp({results.name}, 'Nominal_WithForce'), 1);

if ~isempty(idx_nom_force)
    fig5 = figure('Color','w', 'Name', 'Joint Torques - Nominal With Force');
    tl5 = tiledlayout(fig5, 3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    for j = 1:n
        ax = nexttile(tl5); hold(ax, 'on');
        plot(ax, results(idx_nom_force).sim.t, results(idx_nom_force).sim.tau_log(j,:), ...
            'LineWidth', 1.35);
        grid(ax, 'on');
        xlabel(ax, 'Time [s]');
        ylabel(ax, sprintf('$\\tau_%d$ [Nm]', j));
    end

    if USE_FIG_HELPERS
        saveFigureAsPDF(fig5, fullfile(outdir, 'fig_impedance_joint_torques.pdf'));
    end
end

fprintf('\nSaved results to:\n%s\n', fullfile(pwd, outdir));

%% ========================= Local functions =========================

function dx = closed_loop_rhs_impedance(t, x, Pi, traj, extw, ctrl, imp, con)

n = 6;
q  = x(1:n);
dq = x(n+1:2*n);

[tau, ~] = cartesian_impedance_control(t, q, dq, Pi, traj, extw, ctrl, imp, con);
qdd = forward_dynamics_with_wrench(t, q, dq, tau, Pi, extw, imp);

dx = [dq; qdd];
end

function sim = postprocess_solution_impedance(sol, t_plot, Pi, traj, extw, ctrl, imp, con)

X = deval(sol, t_plot);

n = 6;
N = numel(t_plot);

q_log    = X(1:n, :);
dq_log   = X(n+1:2*n, :);
qdd_log  = zeros(n, N);
tau_log  = zeros(n, N);
p_log    = zeros(3, N);
pd_log   = zeros(3, N);
v_log    = zeros(3, N);
vd_log   = zeros(3, N);
err_log  = zeros(3, N);
fext_log = zeros(3, N);

for k = 1:N
    t  = t_plot(k);
    q  = q_log(:,k);
    dq = dq_log(:,k);

    [tau, ~] = cartesian_impedance_control(t, q, dq, Pi, traj, extw, ctrl, imp, con);
    qdd = forward_dynamics_with_wrench(t, q, dq, tau, Pi, extw, imp);

    [p, Jv] = ee_pos_jac(q, Pi);
    v = Jv * dq;
    [pd, vd, ~] = desired_cartesian_path(t, traj);
    fext = external_force_profile(t, extw, imp.apply_force);

    qdd_log(:,k) = qdd;
    tau_log(:,k) = tau;
    p_log(:,k)   = p;
    pd_log(:,k)  = pd;
    v_log(:,k)   = v;
    vd_log(:,k)  = vd;
    err_log(:,k) = pd - p;
    fext_log(:,k)= fext;
end

sim.t        = t_plot(:).';
sim.q_log    = q_log;
sim.dq_log   = dq_log;
sim.qdd_log  = qdd_log;
sim.tau_log  = tau_log;
sim.p_log    = p_log;
sim.pd_log   = pd_log;
sim.v_log    = v_log;
sim.vd_log   = vd_log;
sim.err_log  = err_log;
sim.fext_log = fext_log;

sim.pos_err_norm = vecnorm(err_log, 2, 1);
sim.fext_norm    = vecnorm(fext_log, 2, 1);
sim.tau_norm     = vecnorm(tau_log, 2, 1);
end

function [tau, data] = cartesian_impedance_control(t, q, dq, Pi, traj, extw, ctrl, imp, con)

M_dyn = UR5_M(q, Pi);
h_dyn = UR5_h(q, dq, Pi);
g_dyn = UR5_G(q, Pi);

[p, Jv, Jdot_dq] = ee_pos_jac_jdotdq(q, dq, Pi, con.fd_dt);
v = Jv * dq;

[pd, vd, ddpd] = desired_cartesian_path(t, traj);
fext = external_force_profile(t, extw, imp.apply_force);

e  = pd - p;
de = vd - v;

% M_d (xdd_d - xdd) + D_d (xd_d - xd) + K_d (x_d - x) = f_ext
acc_task = ddpd + imp.Md \ (imp.Dd * de + imp.Kd * e - fext);

Jv_pinv = damped_pinv(Jv, con.jac_damp);
qdd_cmd = Jv_pinv * (acc_task - Jdot_dq);

tau = M_dyn * qdd_cmd + h_dyn + g_dyn - ctrl.Kq_damp * dq;
tau = clamp_vec(tau, ctrl.tau_max);

data = struct();
data.p        = p;
data.v        = v;
data.pd       = pd;
data.vd       = vd;
data.ddpd     = ddpd;
data.fext     = fext;
data.acc_task = acc_task;
end

function qdd = forward_dynamics_with_wrench(t, q, dq, tau, Pi, extw, imp)

M_dyn = UR5_M(q, Pi);
h_dyn = UR5_h(q, dq, Pi);
g_dyn = UR5_G(q, Pi);

[Jv, ~] = jacobian_pos_only(q, Pi);
fext = external_force_profile(t, extw, imp.apply_force);

qdd = M_dyn \ (tau + Jv.' * fext - h_dyn - g_dyn);
end

function [pd, vd, ddpd] = desired_cartesian_path(t, traj)

switch lower(traj.mode)
    case 'quintic3d'
        p0    = traj.p0;
        pf    = traj.pf;
        Tmove = traj.Tmove;

        pd   = zeros(3,1);
        vd   = zeros(3,1);
        ddpd = zeros(3,1);

        if t <= 0
            pd = p0;
        elseif t >= Tmove
            pd = pf;
        else
            tau = t / Tmove;

            sigma   = 10*tau^3 - 15*tau^4 + 6*tau^5;
            dsigma  = (30*tau^2 - 60*tau^3 + 30*tau^4) / Tmove;
            ddsigma = (60*tau - 180*tau^2 + 120*tau^3) / (Tmove^2);

            dp = pf - p0;

            pd   = p0 + dp * sigma;
            vd   = dp * dsigma;
            ddpd = dp * ddsigma;
        end

    otherwise
        error('Unknown Cartesian trajectory mode.');
end
end

function fext = external_force_profile(t, extw, apply_force)

if ~apply_force
    fext = zeros(3,1);
    return;
end

switch lower(extw.mode)
    case 'half_sine'
        if (t < extw.t_on) || (t > extw.t_off)
            amp = 0.0;
        else
            tau = (t - extw.t_on) / (extw.t_off - extw.t_on);
            amp = extw.f_peak * sin(pi * tau);
        end
        fext = extw.f_axis * amp;

    otherwise
        error('Unknown external force profile.');
end
end

function results = add_impedance_metrics(results)

idx_nom = find(strcmp({results.name}, 'Nominal_NoForce'), 1);
assert(~isempty(idx_nom), 'Nominal_NoForce case is required as the reference case.');

p_nom = results(idx_nom).sim.p_log;
t     = results(idx_nom).sim.t;

for c = 1:numel(results)
    sim = results(c).sim;

    sim.path_dev_vs_nominal = vecnorm(sim.p_log - p_nom, 2, 1);
    results(c).sim = sim;

    metrics.final_pos_err_norm    = sim.pos_err_norm(end);
    metrics.rms_pos_err_norm      = sqrt(trapz(t, sim.pos_err_norm.^2) / (t(end) - t(1)));
    metrics.peak_pos_err_norm     = max(sim.pos_err_norm);
    metrics.peak_path_dev_nominal = max(sim.path_dev_vs_nominal);
    metrics.peak_fext_norm        = max(sim.fext_norm);
    metrics.final_tau_norm        = sim.tau_norm(end);
    metrics.rms_tau_norm          = sqrt(trapz(t, sim.tau_norm.^2) / (t(end) - t(1)));
    metrics.peak_tau_norm         = max(sim.tau_norm);
    metrics.control_energy        = trapz(t, sim.tau_norm.^2);

    results(c).metrics = metrics;
end
end

function summaryTbl = build_summary_table_impedance(results)

nCases = numel(results);

CaseName  = cell(nCases,1);
Label     = cell(nCases,1);
Final_PosErr_m = zeros(nCases,1);
RMS_PosErr_m   = zeros(nCases,1);
Peak_PosErr_m  = zeros(nCases,1);
Peak_PathDev_vs_NoForce_m = zeros(nCases,1);
Peak_ExtForce_N = zeros(nCases,1);
Peak_TauNorm_Nm = zeros(nCases,1);
Control_Energy  = zeros(nCases,1);

for c = 1:nCases
    CaseName{c} = results(c).name;
    Label{c}    = results(c).label;

    Final_PosErr_m(c) = results(c).metrics.final_pos_err_norm;
    RMS_PosErr_m(c)   = results(c).metrics.rms_pos_err_norm;
    Peak_PosErr_m(c)  = results(c).metrics.peak_pos_err_norm;
    Peak_PathDev_vs_NoForce_m(c) = results(c).metrics.peak_path_dev_nominal;
    Peak_ExtForce_N(c) = results(c).metrics.peak_fext_norm;
    Peak_TauNorm_Nm(c) = results(c).metrics.peak_tau_norm;
    Control_Energy(c)  = results(c).metrics.control_energy;
end

summaryTbl = table(CaseName, Label, ...
    Final_PosErr_m, RMS_PosErr_m, Peak_PosErr_m, ...
    Peak_PathDev_vs_NoForce_m, Peak_ExtForce_N, Peak_TauNorm_Nm, ...
    Control_Energy);
end

function write_timeseries_csv_impedance(result, outdir)

sim = result.sim;

t        = sim.t(:);
q_log    = sim.q_log.';
dq_log   = sim.dq_log.';
qdd_log  = sim.qdd_log.';
tau_log  = sim.tau_log.';
p_log    = sim.p_log.';
pd_log   = sim.pd_log.';
v_log    = sim.v_log.';
vd_log   = sim.vd_log.';
err_log  = sim.err_log.';
fext_log = sim.fext_log.';

pos_err_norm = sim.pos_err_norm(:);
fext_norm    = sim.fext_norm(:);
tau_norm     = sim.tau_norm(:);

varNames = [{'time'}, ...
    arrayfun(@(i) sprintf('q%d', i),   1:6, 'UniformOutput', false), ...
    arrayfun(@(i) sprintf('dq%d', i),  1:6, 'UniformOutput', false), ...
    arrayfun(@(i) sprintf('qdd%d', i), 1:6, 'UniformOutput', false), ...
    arrayfun(@(i) sprintf('tau%d', i), 1:6, 'UniformOutput', false), ...
    {'px','py','pz'}, {'pdx','pdy','pdz'}, {'vx','vy','vz'}, {'vdx','vdy','vdz'}, ...
    {'ex','ey','ez'}, {'fextx','fexty','fextz'}, {'pos_err_norm','fext_norm','tau_norm'}];

data = [t, q_log, dq_log, qdd_log, tau_log, ...
        p_log, pd_log, v_log, vd_log, err_log, fext_log, ...
        pos_err_norm, fext_norm, tau_norm];

T = array2table(data, 'VariableNames', varNames);
writetable(T, fullfile(outdir, sprintf('timeseries_%s.csv', result.name)));
end

function [p, Jv] = ee_pos_jac(q, Pi)
T0e = UR5_fkine(q, Pi);
p   = T0e(1:3,4);
[Jv, ~, ~] = UR5_jacobian_geometric(q, Pi);
end

function [p, Jv, Jdot_dq] = ee_pos_jac_jdotdq(q, dq, Pi, fd_dt)
[p, Jv] = ee_pos_jac(q, Pi);

[Jv_plus, ~]  = jacobian_pos_only(q + fd_dt * dq, Pi);
[Jv_minus, ~] = jacobian_pos_only(q - fd_dt * dq, Pi);

Jdot    = (Jv_plus - Jv_minus) / (2 * fd_dt);
Jdot_dq = Jdot * dq;
end

function [Jv, p] = jacobian_pos_only(q, Pi)
[p, Jv] = ee_pos_jac(q, Pi);
end

function Ainv = damped_pinv(A, damp)
Ainv = A' / (A*A' + (damp^2) * eye(size(A,1)));
end

function y = clamp_vec(x, lim)
y = min(max(x, -lim), lim);
end