clc; clear; close all;
addpath('maple_gen')
addpath('utils')

%% ========================= Setup =========================
assert(exist('UR5_params','file')              == 2, 'UR5_params.m not found.');
assert(exist('UR5_M','file')                   == 2, 'UR5_M.m not found.');
assert(exist('UR5_h','file')                   == 2, 'UR5_h.m not found.');
assert(exist('UR5_G','file')                   == 2, 'UR5_G.m not found.');
assert(exist('UR5_fkine','file')               == 2, 'UR5_fkine.m not found.');
assert(exist('UR5_fkall','file')               == 2, 'UR5_fkall.m not found.');
assert(exist('UR5_jacobian_geometric','file')  == 2, 'UR5_jacobian_geometric.m not found.');

USE_FIG_HELPERS = (exist('figureoptscall','file') == 2) && ...
                  (exist('saveFigureAsPDF','file') == 2);

if USE_FIG_HELPERS
    figureoptscall;
end

outdir = 'ur5_hybrid_plane_results';
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

ANIMATE_ROBOT = false;
ANIM_STRIDE   = 20;

%% ========================= Initial state =========================
q0  = [ 0.35; -1.05; 1.20; -1.30; 1.00; 0.20 ];
dq0 = zeros(n,1);

T0 = UR5_fkine(q0, Pi);
p0 = T0(1:3,4);

%% ========================= Plane environment / task frame =========================
env = struct();

% Rigid horizontal plane passing through the initial EE point
env.n_hat = [0; 0; 1];
env.n_hat = env.n_hat / norm(env.n_hat);

env.t1 = [1; 0; 0];
env.t2 = [0; 1; 0];

env.p_plane = p0;                 % point on the plane
env.B       = [env.t1, env.t2, env.n_hat];

% For visualization only
env.plane_half_width = 0.12;

%% ========================= Hybrid task definition =========================
% Motion-controlled tangent coordinates s = [s1; s2]
% Force-controlled normal scalar lambda

traj = struct();
traj.mode  = 'quintic2d';
traj.Tmove = 12.0;
traj.s0    = [0.00; 0.00];
traj.sf    = [0.08; 0.05];

%% ========================= Controller settings =========================
ctrl = struct();

% Motion loop in tangent coordinates: sdd = as
ctrl.Kp_s = diag([60, 60]);
ctrl.Kd_s = diag([18, 18]);

% Force loop in normal direction: lambda_cmd = lambda_d + Ki * integral error
ctrl.lambda_d   = 25.0;      % desired normal contact force [N]
ctrl.Ki_lambda  = 18.0;      % integral gain for force loop
ctrl.lambda_max = 120.0;     % clamp for commanded normal load [N]

% Mild joint damping
ctrl.Kq_damp = diag([3.0, 3.0, 2.0, 1.5, 1.0, 0.8]);

% Torque saturation
ctrl.tau_max = [180; 180; 140; 90; 60; 40];

%% ========================= Numerical settings =========================
con = struct();

% Finite differences
con.fd_dt    = 1e-6;
con.jac_damp = 1e-4;

% Constraint stabilization
con.alpha_baum = 18.0;
con.beta_baum  = 18.0;

% Force-measurement low-pass filter bandwidth
con.lambda_bw = 60.0;

%% ========================= Run simulation =========================
% State:
% x = [ q(6); dq(6); eInt_lambda(1); lambda_hat(1) ]
x0 = [q0; dq0; 0; 0];

rhs = @(t, x) closed_loop_rhs_hybrid_plane(t, x, Pi, env, traj, ctrl, con);

sol = ode45(rhs, [0 Tf], x0, ode_opts);

sim = postprocess_solution_hybrid_plane(sol, t_plot, Pi, env, traj, ctrl, con);
metrics = compute_metrics_hybrid_plane(sim, ctrl);

disp(' ');
disp(struct2table(metrics));

writetable(struct2table(metrics), fullfile(outdir, 'summary_hybrid_plane.csv'));
write_timeseries_csv_hybrid_plane(sim, outdir);

save(fullfile(outdir, 'summary_hybrid_plane.mat'), ...
    'sim', 'metrics', 'env', 'traj', 'ctrl', 'con', 'q0', 'dq0', 'Tf', 't_plot');

%% ========================= Figure 1: Hybrid overview =========================
fig1 = figure('Color','w', 'Name', 'Hybrid Plane Overview');
tl1 = tiledlayout(fig1, 4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile(tl1); hold(ax1, 'on');
plot(ax1, sim.t, sim.s_err_norm, 'LineWidth', 1.5);
grid(ax1, 'on');
xlabel(ax1, 'Time [s]');
ylabel(ax1, '$\|e_s\|_2$ [m]');
%title(ax1, 'Tangential motion error norm');

ax2 = nexttile(tl1); hold(ax2, 'on');
plot(ax2, sim.t, sim.lambda_log,     'LineWidth', 1.5);
plot(ax2, sim.t, sim.lambda_hat_log, '--', 'LineWidth', 1.2);
yline(ax2, ctrl.lambda_d, 'k:', 'LineWidth', 1.2);
grid(ax2, 'on');
xlabel(ax2, 'Time [s]');
ylabel(ax2, 'Force [N]');
%title(ax2, 'Normal contact force');
legend(ax2, {'$\lambda$ actual', '$\hat{\lambda}$ filtered', '$\lambda_d$'}, 'Location', 'northeast');

ax3 = nexttile(tl1); hold(ax3, 'on');
plot(ax3, sim.t, sim.h_log, 'LineWidth', 1.5);
yline(ax3, 0, 'k:', 'LineWidth', 1.2);
grid(ax3, 'on');
xlabel(ax3, 'Time [s]');
ylabel(ax3, '$h(q)$ [m]');
%title(ax3, 'Plane constraint function');

ax4 = nexttile(tl1); hold(ax4, 'on');
plot(ax4, sim.t, sim.tau_norm, 'LineWidth', 1.5);
grid(ax4, 'on');
xlabel(ax4, 'Time [s]');
ylabel(ax4, '$\|\tau\|_2$ [Nm]');
%title(ax4, 'Torque norm');

if USE_FIG_HELPERS
    saveFigureAsPDF(fig1, fullfile(outdir, 'fig_hybrid_overview.pdf'));
end

%% ========================= Figure 2: Tangential coordinates =========================
fig2 = figure('Color','w', 'Name', 'Tangential Coordinates');
tl2 = tiledlayout(fig2, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile(tl2); hold(ax1, 'on');
plot(ax1, sim.t, sim.s_log(1,:),  'LineWidth', 1.5);
plot(ax1, sim.t, sim.sd_log(1,:), 'k--', 'LineWidth', 1.2);
grid(ax1, 'on');
xlabel(ax1, 'Time [s]');
ylabel(ax1, '$s_1$ [m]');
%title(ax1, 'First tangential coordinate');
legend(ax1, {'$s_1$', '$s_{1,d}$'}, 'Location', 'best');

ax2 = nexttile(tl2); hold(ax2, 'on');
plot(ax2, sim.t, sim.s_log(2,:),  'LineWidth', 1.5);
plot(ax2, sim.t, sim.sd_log(2,:), 'k--', 'LineWidth', 1.2);
grid(ax2, 'on');
xlabel(ax2, 'Time [s]');
ylabel(ax2, '$s_2$ [m]');
%title(ax2, 'Second tangential coordinate');
legend(ax2, {'$s_2$', '$s_{2,d}$'}, 'Location', 'best');

if USE_FIG_HELPERS
    saveFigureAsPDF(fig2, fullfile(outdir, 'fig_tangent_coordinates.pdf'));
end

%% ========================= Figure 3: Joint torques =========================
fig3 = figure('Color','w', 'Name', 'Joint Torques');
tl3 = tiledlayout(fig3, 3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for j = 1:n
    ax = nexttile(tl3); hold(ax, 'on');
    plot(ax, sim.t, sim.tau_log(j,:), 'LineWidth', 1.35);
    grid(ax, 'on');
    xlabel(ax, 'Time [s]');
    ylabel(ax, sprintf('$\\tau_%d$ [Nm]', j));
end

if USE_FIG_HELPERS
    saveFigureAsPDF(fig3, fullfile(outdir, 'fig_joint_torques.pdf'));
end

%% ========================= Figure 4: 3D trajectory on plane =========================
fig4 = plot_hybrid_plane_3d(sim, env);

if USE_FIG_HELPERS
    saveFigureAsPDF(fig4, fullfile(outdir, 'fig_3d_hybrid_plane.pdf'));
end

%% ========================= Optional simple animation =========================
if ANIMATE_ROBOT
    animate_ur5_hybrid_plane(sim, Pi, env, ANIM_STRIDE);
end

fprintf('\nSaved results to:\n%s\n', fullfile(pwd, outdir));

%% ========================= Local functions =========================

function dx = closed_loop_rhs_hybrid_plane(t, x, Pi, env, traj, ctrl, con)

n = 6;

q          = x(1:n);
dq         = x(n+1:2*n);
eIntLambda = x(2*n+1);
lambda_hat = x(2*n+2);

[tau, ~] = hybrid_control_law_plane(t, q, dq, eIntLambda, lambda_hat, Pi, env, traj, ctrl, con);
[qdd, lambda_act, ~] = constrained_fdyn_plane(q, dq, tau, Pi, env, con);

eLambda = ctrl.lambda_d - lambda_hat;

dx = [dq;
      qdd;
      eLambda;
      con.lambda_bw * (lambda_act - lambda_hat)];
end

function sim = postprocess_solution_hybrid_plane(sol, t_plot, Pi, env, traj, ctrl, con)

X = deval(sol, t_plot);

n = 6;
N = numel(t_plot);

q_log          = X(1:n, :);
dq_log         = X(n+1:2*n, :);
eIntLambda_log = X(2*n+1, :);
lambda_hat_log = X(2*n+2, :);

qdd_log        = zeros(n, N);
tau_log        = zeros(n, N);
tau_motion_log = zeros(n, N);
p_log          = zeros(3, N);
pd_log         = zeros(3, N);
v_log          = zeros(3, N);
vd_log         = zeros(3, N);
s_log          = zeros(2, N);
sd_log         = zeros(2, N);
ds_log         = zeros(2, N);
dsd_log        = zeros(2, N);
lambda_log     = zeros(1, N);
lambda_cmd_log = zeros(1, N);
h_log          = zeros(1, N);
hdot_log       = zeros(1, N);

for k = 1:N
    t          = t_plot(k);
    q          = q_log(:,k);
    dq         = dq_log(:,k);
    eIntLambda = eIntLambda_log(k);
    lambda_hat = lambda_hat_log(k);

    [tau, ctrl_data] = hybrid_control_law_plane(t, q, dq, eIntLambda, lambda_hat, Pi, env, traj, ctrl, con);
    [qdd, lambda_act, con_data] = constrained_fdyn_plane(q, dq, tau, Pi, env, con);

    [p, Jv] = ee_pos_jac(q, Pi);
    v = Jv * dq;

    [sd, dsd, ~, pd, vd] = desired_tangent_path(t, env, traj);
    [s, ds] = task_coords_plane(p, v, env);

    qdd_log(:,k)        = qdd;
    tau_log(:,k)        = tau;
    tau_motion_log(:,k) = ctrl_data.tau_motion;
    p_log(:,k)          = p;
    pd_log(:,k)         = pd;
    v_log(:,k)          = v;
    vd_log(:,k)         = vd;
    s_log(:,k)          = s;
    sd_log(:,k)         = sd;
    ds_log(:,k)         = ds;
    dsd_log(:,k)        = dsd;
    lambda_log(k)       = lambda_act;
    lambda_cmd_log(k)   = ctrl_data.lambda_cmd;
    h_log(k)            = con_data.h;
    hdot_log(k)         = con_data.hdot;
end

sim.t               = t_plot(:).';
sim.q_log           = q_log;
sim.dq_log          = dq_log;
sim.qdd_log         = qdd_log;
sim.tau_log         = tau_log;
sim.tau_motion_log  = tau_motion_log;
sim.p_log           = p_log;
sim.pd_log          = pd_log;
sim.v_log           = v_log;
sim.vd_log          = vd_log;
sim.s_log           = s_log;
sim.sd_log          = sd_log;
sim.ds_log          = ds_log;
sim.dsd_log         = dsd_log;
sim.lambda_log      = lambda_log;
sim.lambda_cmd_log  = lambda_cmd_log;
sim.lambda_hat_log  = lambda_hat_log;
sim.eIntLambda_log  = eIntLambda_log;
sim.h_log           = h_log;
sim.hdot_log        = hdot_log;

sim.s_err_log       = sd_log - s_log;
sim.ds_err_log      = dsd_log - ds_log;
sim.lambda_err_log  = ctrl.lambda_d - lambda_log;

sim.s_err_norm      = vecnorm(sim.s_err_log,  2, 1);
sim.ds_err_norm     = vecnorm(sim.ds_err_log, 2, 1);
sim.tau_norm        = vecnorm(tau_log, 2, 1);
end

function [tau, data] = hybrid_control_law_plane(t, q, dq, eIntLambda, lambda_hat, Pi, env, traj, ctrl, con)

M_dyn = UR5_M(q, Pi);
h_dyn = UR5_h(q, dq, Pi);
g_dyn = UR5_G(q, Pi);

[p, Jv, Jdot_dq] = ee_pos_jac_jdotdq(q, dq, Pi, con.fd_dt);
v = Jv * dq;

[s, ds]                = task_coords_plane(p, v, env);
[sd, dsd, ddsd, ~, ~] = desired_tangent_path(t, env, traj);

e_s  = sd  - s;
de_s = dsd - ds;

a_s = ddsd + ctrl.Kd_s * de_s + ctrl.Kp_s * e_s;

% Convert desired tangent acceleration back to 3D Cartesian space
a_xyz = env.t1 * a_s(1) + env.t2 * a_s(2);

% Joint-space inverse-dynamics-like motion term
Jv_pinv = damped_pinv(Jv, con.jac_damp);
qdd_cmd = Jv_pinv * (a_xyz - Jdot_dq);

tau_motion = M_dyn * qdd_cmd + h_dyn + g_dyn - ctrl.Kq_damp * dq;

% Force loop in normal direction
[~, A_con, ~] = plane_constraint_geom(q, Pi, env);

lambda_cmd = ctrl.lambda_d + ctrl.Ki_lambda * eIntLambda;
lambda_cmd = clamp_scalar(lambda_cmd, 0, ctrl.lambda_max);

% To increase actual contact reaction, push against the plane normal
tau = tau_motion - A_con.' * lambda_cmd;
tau = clamp_vec(tau, ctrl.tau_max);

data = struct();
data.tau_motion = tau_motion;
data.lambda_cmd = lambda_cmd;
end

function [qdd, lambda, con_data] = constrained_fdyn_plane(q, dq, tau, Pi, env, con)

M_dyn = UR5_M(q, Pi);
h_dyn = UR5_h(q, dq, Pi);
g_dyn = UR5_G(q, Pi);

[h_con, A_con, ~] = plane_constraint_geom(q, Pi, env);
hdot = A_con * dq;

A_plus   = plane_constraint_jacobian(q + con.fd_dt * dq, Pi, env);
A_minus  = plane_constraint_jacobian(q - con.fd_dt * dq, Pi, env);
A_dot    = (A_plus - A_minus) / (2 * con.fd_dt);
A_dot_dq = A_dot * dq;

gamma = A_dot_dq + 2 * con.alpha_baum * hdot + (con.beta_baum^2) * h_con;

KKT = [M_dyn,   -A_con.';
       A_con,    0      ];

rhs = [tau - h_dyn - g_dyn;
       -gamma];

sol = KKT \ rhs;

qdd    = sol(1:6);
lambda = sol(7);

con_data = struct();
con_data.h    = h_con;
con_data.hdot = hdot;
end

function [h_con, A_con, p] = plane_constraint_geom(q, Pi, env)

[p, Jv] = ee_pos_jac(q, Pi);

% Plane constraint h(q) = n^T (p - p_plane) = 0
h_con = env.n_hat.' * (p - env.p_plane);

% Constraint Jacobian A(q) = dh/dq = n^T Jv
A_con = env.n_hat.' * Jv;   % 1x6
end

function A_con = plane_constraint_jacobian(q, Pi, env)
[~, A_con, ~] = plane_constraint_geom(q, Pi, env);
end

function [s, ds] = task_coords_plane(p, v, env)

dp = p - env.p_plane;

s  = [env.t1.' * dp;
      env.t2.' * dp];

ds = [env.t1.' * v;
      env.t2.' * v];
end

function [sd, dsd, ddsd, pd, vd] = desired_tangent_path(t, env, traj)

switch lower(traj.mode)
    case 'quintic2d'
        s0    = traj.s0;
        sf    = traj.sf;
        Tmove = traj.Tmove;

        sd   = zeros(2,1);
        dsd  = zeros(2,1);
        ddsd = zeros(2,1);

        if t <= 0
            sd = s0;
        elseif t >= Tmove
            sd = sf;
        else
            tau = t / Tmove;

            sigma   = 10*tau^3 - 15*tau^4 + 6*tau^5;
            dsigma  = (30*tau^2 - 60*tau^3 + 30*tau^4) / Tmove;
            ddsigma = (60*tau - 180*tau^2 + 120*tau^3) / (Tmove^2);

            ds_total = sf - s0;

            sd   = s0 + ds_total * sigma;
            dsd  = ds_total * dsigma;
            ddsd = ds_total * ddsigma;
        end

    otherwise
        error('Unknown tangent path mode.');
end

pd = env.p_plane + env.t1 * sd(1) + env.t2 * sd(2);
vd = env.t1 * dsd(1) + env.t2 * dsd(2);
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

function metrics = compute_metrics_hybrid_plane(sim, ctrl)

t = sim.t;

metrics.final_s_err_norm      = sim.s_err_norm(end);
metrics.rms_s_err_norm        = sqrt(trapz(t, sim.s_err_norm.^2) / (t(end) - t(1)));
metrics.peak_s_err_norm       = max(sim.s_err_norm);

metrics.final_lambda          = sim.lambda_log(end);
metrics.final_lambda_err      = ctrl.lambda_d - sim.lambda_log(end);
metrics.rms_lambda_err        = sqrt(trapz(t, sim.lambda_err_log.^2) / (t(end) - t(1)));
metrics.peak_abs_lambda_err   = max(abs(sim.lambda_err_log));

metrics.max_abs_h             = max(abs(sim.h_log));
metrics.max_abs_hdot          = max(abs(sim.hdot_log));

metrics.final_tau_norm        = sim.tau_norm(end);
metrics.rms_tau_norm          = sqrt(trapz(t, sim.tau_norm.^2) / (t(end) - t(1)));
metrics.peak_tau_norm         = max(sim.tau_norm);

metrics.control_energy        = trapz(t, sim.tau_norm.^2);
end

function write_timeseries_csv_hybrid_plane(sim, outdir)

t               = sim.t(:);
q_log           = sim.q_log.';
dq_log          = sim.dq_log.';
qdd_log         = sim.qdd_log.';
tau_log         = sim.tau_log.';
p_log           = sim.p_log.';
pd_log          = sim.pd_log.';
v_log           = sim.v_log.';
vd_log          = sim.vd_log.';
s_log           = sim.s_log.';
sd_log          = sim.sd_log.';
ds_log          = sim.ds_log.';
dsd_log         = sim.dsd_log.';
lambda_log      = sim.lambda_log(:);
lambda_cmd_log  = sim.lambda_cmd_log(:);
lambda_hat_log  = sim.lambda_hat_log(:);
eIntLambda_log  = sim.eIntLambda_log(:);
h_log           = sim.h_log(:);
hdot_log        = sim.hdot_log(:);
s_err_norm      = sim.s_err_norm(:);
ds_err_norm     = sim.ds_err_norm(:);
tau_norm        = sim.tau_norm(:);

varNames = [{'time'}, ...
    arrayfun(@(i) sprintf('q%d', i),      1:6, 'UniformOutput', false), ...
    arrayfun(@(i) sprintf('dq%d', i),     1:6, 'UniformOutput', false), ...
    arrayfun(@(i) sprintf('qdd%d', i),    1:6, 'UniformOutput', false), ...
    arrayfun(@(i) sprintf('tau%d', i),    1:6, 'UniformOutput', false), ...
    {'px','py','pz'}, ...
    {'pdx','pdy','pdz'}, ...
    {'vx','vy','vz'}, ...
    {'vdx','vdy','vdz'}, ...
    {'s1','s2'}, ...
    {'sd1','sd2'}, ...
    {'ds1','ds2'}, ...
    {'dsd1','dsd2'}, ...
    {'lambda', 'lambda_cmd', 'lambda_hat', 'eIntLambda', ...
     'h', 'hdot', 's_err_norm', 'ds_err_norm', 'tau_norm'}];

data = [t, ...
    q_log, dq_log, qdd_log, tau_log, ...
    p_log, pd_log, v_log, vd_log, ...
    s_log, sd_log, ds_log, dsd_log, ...
    lambda_log, lambda_cmd_log, lambda_hat_log, eIntLambda_log, ...
    h_log, hdot_log, s_err_norm, ds_err_norm, tau_norm];

T = array2table(data, 'VariableNames', varNames);
writetable(T, fullfile(outdir, 'timeseries_hybrid_plane.csv'));
end

function fig = plot_hybrid_plane_3d(sim, env)

fig = figure('Color','w', 'Name', '3D Hybrid Plane Task');
hold on; grid on; axis equal; view(135, 25);

all_pts = [sim.p_log, sim.pd_log];
mins = min(all_pts, [], 2);
maxs = max(all_pts, [], 2);

xc = 0.5 * (mins(1) + maxs(1));
yc = 0.5 * (mins(2) + maxs(2));
zc = env.p_plane(3);

hw = max(env.plane_half_width, 0.5 * max(maxs(1:2) - mins(1:2)) + 0.05);

Xp = [xc-hw, xc+hw; xc-hw, xc+hw];
Yp = [yc-hw, yc-hw; yc+hw, yc+hw];
Zp = zc * ones(2,2);

surf(Xp, Yp, Zp, ...
    'FaceAlpha', 0.10, ...
    'EdgeAlpha', 0.15, ...
    'HandleVisibility', 'off');

plot3(sim.pd_log(1,:), sim.pd_log(2,:), sim.pd_log(3,:), ...
    'k--', 'LineWidth', 1.25, 'DisplayName', 'Desired path');

plot3(sim.p_log(1,:), sim.p_log(2,:), sim.p_log(3,:), ...
    'b-', 'LineWidth', 1.8, 'DisplayName', 'Actual path');

plot3(sim.p_log(1,1), sim.p_log(2,1), sim.p_log(3,1), ...
    'go', 'MarkerFaceColor', 'g', 'DisplayName', 'Start');

plot3(sim.p_log(1,end), sim.p_log(2,end), sim.p_log(3,end), ...
    'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'End');

xlabel('x [m]');
ylabel('y [m]');
zlabel('z [m]');
%title('Hybrid force/motion task on a rigid plane');
legend('Location', 'best');
end

function animate_ur5_hybrid_plane(sim, Pi, env, frame_stride)

fig = figure('Color','w', 'Name', 'UR5 Hybrid Plane Animation');

all_pts = [sim.p_log, sim.pd_log];
mins = min(all_pts, [], 2) - 0.35;
maxs = max(all_pts, [], 2) + 0.35;

for k = 1:frame_stride:numel(sim.t)
    clf(fig);
    hold on; grid on; axis equal; view(135, 25);

    xc = 0.5 * (mins(1) + maxs(1));
    yc = 0.5 * (mins(2) + maxs(2));
    zc = env.p_plane(3);
    hw = env.plane_half_width + 0.20;

    Xp = [xc-hw, xc+hw; xc-hw, xc+hw];
    Yp = [yc-hw, yc-hw; yc+hw, yc+hw];
    Zp = zc * ones(2,2);

    surf(Xp, Yp, Zp, 'FaceAlpha', 0.08, 'EdgeAlpha', 0.08, 'HandleVisibility', 'off');

    plot3(sim.pd_log(1,:), sim.pd_log(2,:), sim.pd_log(3,:), 'k--', 'LineWidth', 1.0);
    plot3(sim.p_log(1,1:k), sim.p_log(2,1:k), sim.p_log(3,1:k), 'b-', 'LineWidth', 1.6);

    [Tlist, ~] = UR5_fkall(sim.q_log(:,k), Pi);
    P = zeros(3, numel(Tlist));
    for i = 1:numel(Tlist)
        P(:,i) = Tlist{i}(1:3,4);
    end

    plot3(P(1,:), P(2,:), P(3,:), 'o-', ...
        'LineWidth', 1.6, ...
        'MarkerSize', 5, ...
        'DisplayName', 'UR5');

    plot3(sim.p_log(1,k), sim.p_log(2,k), sim.p_log(3,k), ...
        'ro', 'MarkerFaceColor', 'r', 'HandleVisibility', 'off');

    xlim([mins(1), maxs(1)]);
    ylim([mins(2), maxs(2)]);
    zlim([mins(3), maxs(3)]);

    xlabel('x [m]');
    ylabel('y [m]');
    zlabel('z [m]');
    %title(sprintf('t = %.2f s', sim.t(k)));

    drawnow;
end
end

function Ainv = damped_pinv(A, damp)
Ainv = A' / (A*A' + (damp^2) * eye(size(A,1)));
end

function y = clamp_vec(x, lim)
y = min(max(x, -lim), lim);
end

function y = clamp_scalar(x, xmin, xmax)
y = min(max(x, xmin), xmax);
end