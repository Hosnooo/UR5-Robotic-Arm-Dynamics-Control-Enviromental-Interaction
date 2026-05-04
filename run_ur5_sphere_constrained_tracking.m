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

outdir = 'ur5_sphere_constrained_results';
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

%% ========================= Robot parameters =========================
[~, Pi] = UR5_params();
n = 6;

%% ========================= Simulation settings =========================
Tf     = 20.0;
t_plot = linspace(0, Tf, 4001);

ode_opts = odeset( ...
    'RelTol', 1e-6, ...
    'AbsTol', 1e-8, ...
    'MaxStep', 1e-2);

ANIMATE_ROBOT = false;     % set true if you want simple robot animation
ANIM_STRIDE   = 20;

%% ========================= Initial state =========================
q0  = [ 0.35; -1.10;  1.25; -1.35;  1.05;  0.25 ];
dq0 = zeros(n,1);

T0 = UR5_fkine(q0, Pi);
r0 = T0(1:3, 4);

%% ========================= Sphere environment =========================
env = struct();

env.R = 0.18;                         % sphere radius [m]

% Choose the initial radial direction so that q0 starts exactly on the sphere.
radial0 = [0; 0; 1];
radial0 = radial0 / norm(radial0);

env.r_c = r0 - env.R * radial0;       % sphere center so that r0 lies on the sphere

% Great-circle path basis on the sphere:
env.a = (r0 - env.r_c) / norm(r0 - env.r_c);   % ensures r_d(0) = r0
env.b = build_perp_unit(env.a);                 % tangent direction at t = 0

env.omega = 0.30;                    % angular speed along the sphere [rad/s]

%% ========================= Controller settings =========================
ctrl = struct();

% Tangential task-space PD action
ctrl.Kp_t = diag([220, 220, 220]);   % [N/m]
ctrl.Kd_t = diag([ 40,  40,  40]);   % [N s/m]

% Mild joint damping for numerical smoothness
ctrl.Kq_damp = diag([3.0, 3.0, 2.0, 1.5, 1.0, 0.8]);

% Torque saturation
ctrl.tau_max = [150; 150; 120; 80; 50; 35];

%% ========================= Constraint solve settings =========================
con = struct();

% Numerical differentiation step for A_dot * dq
con.fd_dt = 1e-6;

% Baumgarte stabilization for numerical drift reduction.
% Set both to zero if you want the pure slide equations only.
con.alpha_baum = 12.0;
con.beta_baum  = 12.0;

%% ========================= Run simulation =========================
x0 = [q0; dq0];

rhs = @(t, x) closed_loop_rhs_sphere(t, x, Pi, env, ctrl, con);

sol = ode45(rhs, [0 Tf], x0, ode_opts);

sim = postprocess_solution_sphere(sol, t_plot, Pi, env, ctrl, con);
metrics = compute_metrics_sphere(sim);

disp(' ');
disp(struct2table(metrics));

writetable(struct2table(metrics), fullfile(outdir, 'summary_sphere_constrained.csv'));
write_timeseries_csv_sphere(sim, outdir);

save(fullfile(outdir, 'summary_sphere_constrained.mat'), ...
    'sim', 'metrics', 'env', 'ctrl', 'con', 'q0', 'dq0', 'Tf', 't_plot');

%% ========================= Figure 1: constraint and control overview =========================
fig1 = figure('Color','w', 'Name', 'Sphere-Constrained Overview');
tl1 = tiledlayout(fig1, 4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile(tl1); hold(ax1, 'on');
plot(ax1, sim.t, sim.h_log, 'LineWidth', 1.5);
grid(ax1, 'on');
xlabel(ax1, 'Time [s]');
ylabel(ax1, '$h(q)$');
%title(ax1, 'Constraint function');
ylim([-1 1])
ax2 = nexttile(tl1); hold(ax2, 'on');
plot(ax2, sim.t, sim.radius_err_log, 'LineWidth', 1.5);
yline(ax2, 0, 'k:', 'LineWidth', 1.2);
grid(ax2, 'on');
xlabel(ax2, 'Time [s]');
ylabel(ax2, '$\|r-r_c\|-R$ [m]');
%title(ax2, 'Radial constraint error');
ylim([-1 1])

ax3 = nexttile(tl1); hold(ax3, 'on');
plot(ax3, sim.t, sim.lambda_log, 'LineWidth', 1.5);
grid(ax3, 'on');
xlabel(ax3, 'Time [s]');
ylabel(ax3, '$\lambda$');
%title(ax3, 'Constraint multiplier');

ax4 = nexttile(tl1); hold(ax4, 'on');
plot(ax4, sim.t, sim.tau_norm, 'LineWidth', 1.5);
grid(ax4, 'on');
xlabel(ax4, 'Time [s]');
ylabel(ax4, '$\|\tau\|_2$ [Nm]');
%title(ax4, 'Torque norm');

if USE_FIG_HELPERS
    saveFigureAsPDF(fig1, fullfile(outdir, 'fig_constraint_overview.pdf'));
end

%% ========================= Figure 2: tangential tracking =========================
fig2 = figure('Color','w', 'Name', 'Tangential Tracking');
tl2 = tiledlayout(fig2, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile(tl2); hold(ax1, 'on');
plot(ax1, sim.t, sim.e_t_norm, 'LineWidth', 1.5);
grid(ax1, 'on');
xlabel(ax1, 'Time [s]');
ylabel(ax1, '$\|e_t\|_2$ [m]');
%title(ax1, 'Tangential position tracking error');

ax2 = nexttile(tl2); hold(ax2, 'on');
plot(ax2, sim.t, sim.de_t_norm, 'LineWidth', 1.5);
grid(ax2, 'on');
xlabel(ax2, 'Time [s]');
ylabel(ax2, '$\|\dot e_t\|_2$ [m/s]');
%title(ax2, 'Tangential velocity tracking error');

if USE_FIG_HELPERS
    saveFigureAsPDF(fig2, fullfile(outdir, 'fig_tangential_tracking.pdf'));
end

%% ========================= Figure 3: joint torques =========================
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

%% ========================= Figure 4: 3D constrained trajectory =========================
fig4 = plot_sphere_trajectory_3d(sim, env);

if USE_FIG_HELPERS
    saveFigureAsPDF(fig4, fullfile(outdir, 'fig_3d_sphere_trajectory.pdf'));
end

%% ========================= Optional simple animation =========================
if ANIMATE_ROBOT
    animate_ur5_on_sphere(sim, Pi, env, ANIM_STRIDE);
end

fprintf('\nSaved results to:\n%s\n', fullfile(pwd, outdir));

%% ========================= Local functions =========================

function dx = closed_loop_rhs_sphere(t, x, Pi, env, ctrl, con)

n  = 6;
q  = x(1:n);
dq = x(n+1:2*n);

tau = sphere_control_law(t, q, dq, Pi, env, ctrl);
[qdd, ~, ~] = constrained_fdyn_sphere(q, dq, tau, Pi, env, con);

dx = [dq; qdd];
end

function sim = postprocess_solution_sphere(sol, t_plot, Pi, env, ctrl, con)

X = deval(sol, t_plot);

n = 6;
N = numel(t_plot);

q_log          = X(1:n, :);
dq_log         = X(n+1:2*n, :);

qdd_log        = zeros(n, N);
tau_log        = zeros(n, N);
r_log          = zeros(3, N);
dr_log         = zeros(3, N);
rd_log         = zeros(3, N);
drd_log        = zeros(3, N);
h_log          = zeros(1, N);
hdot_log       = zeros(1, N);
lambda_log     = zeros(1, N);
radius_err_log = zeros(1, N);
e_t_log        = zeros(3, N);
de_t_log       = zeros(3, N);

for k = 1:N
    q  = q_log(:,k);
    dq = dq_log(:,k);

    tau = sphere_control_law(t_plot(k), q, dq, Pi, env, ctrl);
    [qdd, lambda, con_data] = constrained_fdyn_sphere(q, dq, tau, Pi, env, con);
    [rd, drd, ~] = desired_sphere_path(t_plot(k), env);

    [T0e, Jv] = ee_pos_jac(q, Pi);
    r  = T0e(1:3, 4);
    dr = Jv * dq;

    n_hat = (r - env.r_c) / norm(r - env.r_c);
    P_t   = eye(3) - n_hat * n_hat.';

    e_t  = P_t * (rd  - r);
    de_t = P_t * (drd - dr);

    qdd_log(:,k)        = qdd;
    tau_log(:,k)        = tau;
    r_log(:,k)          = r;
    dr_log(:,k)         = dr;
    rd_log(:,k)         = rd;
    drd_log(:,k)        = drd;
    h_log(k)            = con_data.h;
    hdot_log(k)         = con_data.hdot;
    lambda_log(k)       = lambda;
    radius_err_log(k)   = con_data.radius_err;
    e_t_log(:,k)        = e_t;
    de_t_log(:,k)       = de_t;
end

sim.t              = t_plot(:).';
sim.q_log          = q_log;
sim.dq_log         = dq_log;
sim.qdd_log        = qdd_log;
sim.tau_log        = tau_log;
sim.r_log          = r_log;
sim.dr_log         = dr_log;
sim.rd_log         = rd_log;
sim.drd_log        = drd_log;
sim.h_log          = h_log;
sim.hdot_log       = hdot_log;
sim.lambda_log     = lambda_log;
sim.radius_err_log = radius_err_log;
sim.e_t_log        = e_t_log;
sim.de_t_log       = de_t_log;

sim.tau_norm  = vecnorm(tau_log, 2, 1);
sim.e_t_norm  = vecnorm(e_t_log,  2, 1);
sim.de_t_norm = vecnorm(de_t_log, 2, 1);
end

function tau = sphere_control_law(t, q, dq, Pi, env, ctrl)

[T0e, Jv] = ee_pos_jac(q, Pi);

r  = T0e(1:3, 4);
dr = Jv * dq;

[rd, drd, ~] = desired_sphere_path(t, env);

n_hat = (r - env.r_c) / norm(r - env.r_c);
P_t   = eye(3) - n_hat * n_hat.';

e_t  = P_t * (rd  - r);
de_t = P_t * (drd - dr);

F_t = ctrl.Kp_t * e_t + ctrl.Kd_t * de_t;

% Joint torque input u = tau
tau = Jv.' * F_t + UR5_G(q, Pi) - ctrl.Kq_damp * dq;
tau = clamp_vec(tau, ctrl.tau_max);
end

function [qdd, lambda, con_data] = constrained_fdyn_sphere(q, dq, tau, Pi, env, con)

M_dyn = UR5_M(q, Pi);
c_dyn = UR5_h(q, dq, Pi);
g_dyn = UR5_G(q, Pi);

[h_con, A_con, r] = sphere_constraint_geom(q, Pi, env);
hdot = A_con * dq;

% Numerical approximation of A_dot * dq
A_plus  = sphere_constraint_jacobian(q + con.fd_dt * dq, Pi, env);
A_minus = sphere_constraint_jacobian(q - con.fd_dt * dq, Pi, env);
A_dot   = (A_plus - A_minus) / (2 * con.fd_dt);
A_dot_dq = A_dot * dq;

% Baumgarte stabilization term for numerical robustness
gamma = A_dot_dq + 2 * con.alpha_baum * hdot + (con.beta_baum^2) * h_con;

% KKT solve:
% [ M  -A' ] [qdd   ] = [tau - c - g]
% [ A   0  ] [lambda]   [ -gamma   ]
KKT = [M_dyn,   -A_con.';
       A_con,    0     ];

rhs = [tau - c_dyn - g_dyn;
       -gamma];

sol = KKT \ rhs;

qdd    = sol(1:6);
lambda = sol(7);

radius_err = norm(r - env.r_c) - env.R;

con_data = struct();
con_data.h          = h_con;
con_data.hdot       = hdot;
con_data.radius_err = radius_err;
end

function [h_con, A_con, r] = sphere_constraint_geom(q, Pi, env)

[T0e, Jv] = ee_pos_jac(q, Pi);
r = T0e(1:3, 4);

delta = r - env.r_c;

% h(q) = ||r-r_c||^2 - R^2 = 0
h_con = delta.' * delta - env.R^2;

% A(q) = dh/dq = 2 (r-r_c)^T J_r(q)
A_con = 2 * delta.' * Jv;      % 1x6 row vector
end

function A_con = sphere_constraint_jacobian(q, Pi, env)
[~, A_con, ~] = sphere_constraint_geom(q, Pi, env);
end

function [rd, drd, ddrd] = desired_sphere_path(t, env)

theta = env.omega * t;

u    =  cos(theta) * env.a + sin(theta) * env.b;
du   = env.omega * (-sin(theta) * env.a + cos(theta) * env.b);
ddu  = env.omega^2 * (-cos(theta) * env.a - sin(theta) * env.b);

rd   = env.r_c + env.R * u;
drd  = env.R * du;
ddrd = env.R * ddu;
end

function [T0e, Jv] = ee_pos_jac(q, Pi)
T0e = UR5_fkine(q, Pi);
[Jv, ~, ~] = UR5_jacobian_geometric(q, Pi);
end

function metrics = compute_metrics_sphere(sim)

t = sim.t;

metrics.initial_tangent_err_norm = sim.e_t_norm(1);
metrics.final_tangent_err_norm   = sim.e_t_norm(end);
metrics.rms_tangent_err_norm     = sqrt(trapz(t, sim.e_t_norm.^2)  / (t(end) - t(1)));
metrics.rms_tangent_vel_err_norm = sqrt(trapz(t, sim.de_t_norm.^2) / (t(end) - t(1)));

metrics.max_abs_h                = max(abs(sim.h_log));
metrics.max_abs_radius_error     = max(abs(sim.radius_err_log));

metrics.final_lambda             = sim.lambda_log(end);
metrics.peak_abs_lambda          = max(abs(sim.lambda_log));

metrics.final_tau_norm           = sim.tau_norm(end);
metrics.rms_tau_norm             = sqrt(trapz(t, sim.tau_norm.^2) / (t(end) - t(1)));
metrics.peak_tau_norm            = max(sim.tau_norm);

metrics.control_energy           = trapz(t, sim.tau_norm.^2);
end

function write_timeseries_csv_sphere(sim, outdir)

t              = sim.t(:);
q_log          = sim.q_log.';
dq_log         = sim.dq_log.';
qdd_log        = sim.qdd_log.';
tau_log        = sim.tau_log.';
r_log          = sim.r_log.';
dr_log         = sim.dr_log.';
rd_log         = sim.rd_log.';
drd_log        = sim.drd_log.';
e_t_log        = sim.e_t_log.';
de_t_log       = sim.de_t_log.';
h_log          = sim.h_log(:);
hdot_log       = sim.hdot_log(:);
lambda_log     = sim.lambda_log(:);
radius_err_log = sim.radius_err_log(:);
tau_norm       = sim.tau_norm(:);
e_t_norm       = sim.e_t_norm(:);
de_t_norm      = sim.de_t_norm(:);

varNames = [{'time'}, ...
    arrayfun(@(i) sprintf('q%d', i),    1:6, 'UniformOutput', false), ...
    arrayfun(@(i) sprintf('dq%d', i),   1:6, 'UniformOutput', false), ...
    arrayfun(@(i) sprintf('qdd%d', i),  1:6, 'UniformOutput', false), ...
    arrayfun(@(i) sprintf('tau%d', i),  1:6, 'UniformOutput', false), ...
    {'rx','ry','rz'}, ...
    {'drx','dry','drz'}, ...
    {'rdx','rdy','rdz'}, ...
    {'drdx','drdy','drdz'}, ...
    {'etx','ety','etz'}, ...
    {'detx','dety','detz'}, ...
    {'h', 'hdot', 'lambda', 'radius_err', 'tau_norm', 'e_t_norm', 'de_t_norm'}];

data = [ ...
    t, ...
    q_log, dq_log, qdd_log, tau_log, ...
    r_log, dr_log, rd_log, drd_log, ...
    e_t_log, de_t_log, ...
    h_log, hdot_log, lambda_log, radius_err_log, tau_norm, e_t_norm, de_t_norm];

T = array2table(data, 'VariableNames', varNames);
writetable(T, fullfile(outdir, 'timeseries_sphere_constrained.csv'));
end

function fig = plot_sphere_trajectory_3d(sim, env)

fig = figure('Color','w', 'Name', '3D Sphere-Constrained Trajectory');
hold on; grid on; axis equal; view(135, 25);

[Xs, Ys, Zs] = sphere(60);
surf(env.r_c(1) + env.R * Xs, ...
     env.r_c(2) + env.R * Ys, ...
     env.r_c(3) + env.R * Zs, ...
     'FaceAlpha', 0.12, ...
     'EdgeAlpha', 0.08, ...
     'HandleVisibility', 'off');

plot3(sim.rd_log(1,:), sim.rd_log(2,:), sim.rd_log(3,:), ...
    'k--', 'LineWidth', 1.25, 'DisplayName', 'Desired path');

plot3(sim.r_log(1,:), sim.r_log(2,:), sim.r_log(3,:), ...
    'b-', 'LineWidth', 1.8, 'DisplayName', 'Actual path');

plot3(sim.r_log(1,1), sim.r_log(2,1), sim.r_log(3,1), ...
    'go', 'MarkerFaceColor', 'g', 'DisplayName', 'Start');

plot3(sim.r_log(1,end), sim.r_log(2,end), sim.r_log(3,end), ...
    'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'End');

xlabel('x [m]');
ylabel('y [m]');
zlabel('z [m]');
%title('UR5 end-effector trajectory on a spherical constraint');
legend('Location', 'northeast');
end

function animate_ur5_on_sphere(sim, Pi, env, frame_stride)

fig = figure('Color','w', 'Name', 'UR5 Sphere-Constrained Animation');

all_pts = [sim.r_log, sim.rd_log];
mins = min(all_pts, [], 2) - 0.35;
maxs = max(all_pts, [], 2) + 0.35;

for k = 1:frame_stride:numel(sim.t)
    clf(fig);
    hold on; grid on; axis equal; view(135, 25);

    [Xs, Ys, Zs] = sphere(40);
    surf(env.r_c(1) + env.R * Xs, ...
         env.r_c(2) + env.R * Ys, ...
         env.r_c(3) + env.R * Zs, ...
         'FaceAlpha', 0.10, ...
         'EdgeAlpha', 0.06, ...
         'HandleVisibility', 'off');

    plot3(sim.rd_log(1,:), sim.rd_log(2,:), sim.rd_log(3,:), ...
        'k--', 'LineWidth', 1.0);
    plot3(sim.r_log(1,1:k), sim.r_log(2,1:k), sim.r_log(3,1:k), ...
        'b-', 'LineWidth', 1.6);

    [Tlist, ~] = UR5_fkall(sim.q_log(:,k), Pi);
    P = zeros(3, numel(Tlist));
    for i = 1:numel(Tlist)
        P(:,i) = Tlist{i}(1:3,4);
    end

    plot3(P(1,:), P(2,:), P(3,:), 'o-', ...
        'LineWidth', 1.6, ...
        'MarkerSize', 5, ...
        'DisplayName', 'UR5');

    plot3(sim.r_log(1,k), sim.r_log(2,k), sim.r_log(3,k), ...
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

function b = build_perp_unit(a)

a = a / norm(a);

candidates = eye(3);
proj = abs(candidates.' * a);
[~, idx] = min(proj);

v = candidates(:, idx);
b = v - a * (a.' * v);
b = b / norm(b);
end

function y = clamp_vec(x, lim)
y = min(max(x, -lim), lim);
end