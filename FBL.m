clc; clear; close all;

% Apply global LaTeX + font formatting to all figures
figureoptscall;

outdir = 'q2_figures';
if ~exist(outdir, 'dir'), mkdir(outdir); end

%% ---- Parameters --------------------------------------------
[~, Pi_true] = UR5_params();
n = 6;

%% ---- Simulation settings -----------------------------------
Tf     = 15.0;
t_plot = linspace(0, Tf, 3001);
ode_opts = odeset('RelTol',1e-7,'AbsTol',1e-9,'MaxStep',5e-3);

%% ---- Trajectory --------------------------------------------
traj.q_bias  = [ 0.20; -0.70;  0.70; -0.20;  0.15;  0.00];
traj.q_amp   = [ 0.30;  0.20;  0.18;  0.25;  0.12;  0.15];
traj.w_ref   = 2*pi*[0.10; 0.08; 0.12; 0.07; 0.09; 0.11];
traj.phi_ref = [ 0.0;  pi/4;  pi/2;  pi/6;  pi/3;  pi/5];

%% ---- Initial state -----------------------------------------
[q_ref0, dq_ref0, ~] = joint_trajectory(0, traj);
q0  = q_ref0  + [0.12; -0.08; 0.06; 0.10; -0.05; 0.08];
dq0 = dq_ref0;

%% ---- Gains -------------------------------------------------
Kp_fbl     = diag([200; 200; 200;  80;  80;  80]);
Kd_fbl     = diag([ 28;  28;  28;  18;  18;  18]);
Kp_pd      = diag([ 12;  12;  10;   6;   4;   3]);
Kd_pd      = diag([  5;   5;   4;   2;   1;   1]);
Kp_fbl_pid = Kp_fbl;
Kd_fbl_pid = Kd_fbl;
Ki_fbl_pid = diag([8; 8; 8; 3; 3; 3]);

%% ---- Controller list ---------------------------------------
CTRL_LIST   = {'FFW_PD', 'FBL_PD', 'FBL_PID'};
CTRL_LABELS = {'A: $\mathrm{FFW}+PD_q$', 'B: $\mathrm{FBL}+PD_x$', 'C: $\mathrm{FBL}+PID_x$'};
nCtrl = numel(CTRL_LIST);

%% ---- Simulate ----------------------------------------------
results = struct([]);

for c = 1:nCtrl
    fprintf('Running %s ...\n', CTRL_LIST{c});
    if strcmp(CTRL_LIST{c},'FBL_PID')
        x0 = [q0; dq0; zeros(6,1)];
    else
        x0 = [q0; dq0];
    end
    rhs = @(t,x) ode_rhs(t, x, Pi_true, CTRL_LIST{c}, traj, ...
                          Kp_fbl, Kd_fbl, Kp_pd, Kd_pd, ...
                          Kp_fbl_pid, Kd_fbl_pid, Ki_fbl_pid);
    sol = ode45(rhs, [0 Tf], x0, ode_opts);
    sim = postprocess(sol, t_plot, Pi_true, CTRL_LIST{c}, traj, ...
                      Kp_fbl, Kd_fbl, Kp_pd, Kd_pd, ...
                      Kp_fbl_pid, Kd_fbl_pid, Ki_fbl_pid);
    results(c).name  = CTRL_LIST{c};
    results(c).label = CTRL_LABELS{c};
    results(c).sim   = sim;
end

%% ---- Plot styles -------------------------------------------
cols = {[0.20 0.20 0.80], [0.85 0.33 0.10], [0.47 0.67 0.19]};
styles = {':', '-', '--'};
lw = 1.6;
fbl_idx = find(strcmp(CTRL_LIST,'FBL_PD'));

%% 
%  FIGURE 1 — Overview

fig1 = figure('Color','w');

tl1 = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

ax1 = nexttile; hold on; grid on; box on;
for c = 1:nCtrl
    plot(ax1, results(c).sim.t, results(c).sim.pos_err*1e3, ...
        styles{c}, 'LineWidth', lw, 'Color', cols{c});
end
ylabel(ax1, '$\|e_p\|_2$ [mm]');
xlabel(ax1, 'Time [s]');
legend(ax1, CTRL_LABELS, 'Location','northeast');

ax2 = nexttile; hold on; grid on; box on;
for c = 1:nCtrl
    plot(ax2, results(c).sim.t, rad2deg(results(c).sim.ori_err), ...
        styles{c}, 'LineWidth', lw, 'Color', cols{c});
end
ylabel(ax2, '$\|e_o\|_2$ [deg]');
xlabel(ax2, 'Time [s]');
% legend(ax2, CTRL_LABELS, 'Location','northeast');

ax3 = nexttile; hold on; grid on; box on;
for c = 1:nCtrl
    plot(ax3, results(c).sim.t, results(c).sim.tau_norm, ...
        styles{c}, 'LineWidth', lw, 'Color', cols{c});
end
ylabel(ax3, '$\|\tau\|_2$ [Nm]');
xlabel(ax3, 'Time [s]');
% legend(ax3, CTRL_LABELS, 'Location','northeast');

saveFigureAsPDF(fig1, fullfile(outdir,'q2_fig1_overview.pdf'));

%%
%  FIGURE 2 — Task-space error components

fig2 = figure('Color','w');

tl2 = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

pos_lbl = {'$e_{p,x}$','$e_{p,y}$','$e_{p,z}$'};
ori_lbl = {'$e_{o,1}$','$e_{o,2}$','$e_{o,3}$'};

for j = 1:3
    ax = nexttile(tl2); hold on; grid on; box on;
    for c = 1:nCtrl
        plot(ax, results(c).sim.t, results(c).sim.ex_log(j,:)*1e3, ...
            styles{c}, 'LineWidth', lw, 'Color', cols{c});
    end
    yline(ax, 0, 'k:', 'LineWidth', 1.0);
    ylabel(ax, [pos_lbl{j} ' [mm]']);
    xlabel(ax, 'Time [s]');
    if j == 1, legend(ax, CTRL_LABELS, 'Location','best'); end
end

for j = 1:3
    ax = nexttile(tl2); hold on; grid on; box on;
    for c = 1:nCtrl
        plot(ax, results(c).sim.t, rad2deg(results(c).sim.ex_log(3+j,:)), ...
            styles{c}, 'LineWidth', lw, 'Color', cols{c});
    end
    yline(ax, 0, 'k:', 'LineWidth', 1.0);
    ylabel(ax, [ori_lbl{j} ' [deg]']);
    xlabel(ax, 'Time [s]');
end

saveFigureAsPDF(fig2, fullfile(outdir,'q2_fig2_error_components.pdf'));

%%
%  FIGURE 3 — 3D path + position components

fig3 = figure('Color','w');

tl3 = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

% 3D path
ax3a = nexttile(tl3);
hold(ax3a,'on'); grid(ax3a,'on'); box(ax3a,'on');
axis(ax3a,'equal'); view(ax3a,3);

plot3(ax3a, ...
    results(fbl_idx).sim.pos_ref(1,:)*1e3, ...
    results(fbl_idx).sim.pos_ref(2,:)*1e3, ...
    results(fbl_idx).sim.pos_ref(3,:)*1e3, ...
    'k--', 'LineWidth', 1.4);
plot3(ax3a, ...
    results(fbl_idx).sim.pos_log(1,:)*1e3, ...
    results(fbl_idx).sim.pos_log(2,:)*1e3, ...
    results(fbl_idx).sim.pos_log(3,:)*1e3, ...
    '-', 'LineWidth', lw, 'Color', cols{fbl_idx});

xlabel(ax3a,'$x$ [mm]'); ylabel(ax3a,'$y$ [mm]'); zlabel(ax3a,'$z$ [mm]');
legend(ax3a, {'$p_d$','$p$'}, 'Location','best');

% x,y,z vs time
ax3b = nexttile(tl3);
hold(ax3b,'on'); grid(ax3b,'on'); box(ax3b,'on');

clr_xyz = {[0.8 0.1 0.1],[0.1 0.6 0.1],[0.1 0.1 0.8]};
axis_lbl = {'$x$','$y$','$z$'};
t = results(fbl_idx).sim.t;
leg_entries = {};

for j = 1:3
    plot(ax3b, t, results(fbl_idx).sim.pos_log(j,:)*1e3, '-', ...
        'LineWidth', lw, 'Color', clr_xyz{j});
    plot(ax3b, t, results(fbl_idx).sim.pos_ref(j,:)*1e3, 'k--', ...
        'LineWidth', 1.2);
    leg_entries{end+1} = axis_lbl{j};
    leg_entries{end+1} = ['$' char('x'+(j-1)) '_d$'];
end

xlabel(ax3b, 'Time [s]');
ylabel(ax3b, '$p$ [mm]');
legend(ax3b, leg_entries, 'Location','best', 'NumColumns', 3);


saveFigureAsPDF(fig3, fullfile(outdir,'q2_fig3_cartesian_tracking.pdf'));

%% 
%  FIGURE 4 — Decoupled error dynamics

fig4 = figure('Color','w');

tl4 = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

fbl_mask   = find(~strcmp(CTRL_LIST,'FFW_PD'));
fbl_labels = CTRL_LABELS(fbl_mask);
fbl_styles = styles(fbl_mask);
fbl_cols   = cols(fbl_mask);

for j = 1:3
    ax = nexttile(tl4); hold on; grid on; box on;
    for ci = 1:numel(fbl_mask)
        c = fbl_mask(ci);
        plot(ax, results(c).sim.t, results(c).sim.ex_log(j,:)*1e3, ...
            fbl_styles{ci}, 'LineWidth', lw, 'Color', fbl_cols{ci});
    end
    yline(ax, 0, 'k:', 'LineWidth', 1.0);
    ylabel(ax, [pos_lbl{j} ' [mm]']);
    xlabel(ax, 'Time [s]');
    if j == 1, legend(ax, fbl_labels, 'Location','best'); end
end

for j = 1:3
    ax = nexttile(tl4); hold on; grid on; box on;
    for ci = 1:numel(fbl_mask)
        c = fbl_mask(ci);
        plot(ax, results(c).sim.t, rad2deg(results(c).sim.ex_log(3+j,:)), ...
            fbl_styles{ci}, 'LineWidth', lw, 'Color', fbl_cols{ci});
    end
    yline(ax, 0, 'k:', 'LineWidth', 1.0);
    ylabel(ax, [ori_lbl{j} ' [deg]']);
    xlabel(ax, 'Time [s]');
end

saveFigureAsPDF(fig4, fullfile(outdir,'q2_fig4_decoupled_errors.pdf'));

%% 
%  FIGURE 5 — Joint torques

fig5 = figure('Color','w');

tl5 = tiledlayout(3,2,'TileSpacing','compact','Padding','compact');

for j = 1:6
    ax = nexttile(tl5); hold on; grid on; box on;
    for c = 1:nCtrl
        plot(ax, results(c).sim.t, results(c).sim.tau_log(j,:), ...
            styles{c}, 'LineWidth', lw, 'Color', cols{c});
    end
    ylabel(ax, sprintf('$\\tau_%d$ [Nm]', j));
    xlabel(ax, 'Time [s]');
    if j == 1, legend(ax, CTRL_LABELS, 'Location','best'); end
end

saveFigureAsPDF(fig5, fullfile(outdir,'q2_fig5_joint_torques.pdf'));

fprintf('\nAll Q2 figures saved to: %s\n', fullfile(pwd, outdir));

%% 
%  LOCAL FUNCTIONS

function dx = ode_rhs(t, x, Pi, mode, traj, ...
                       Kp_fbl, Kd_fbl, Kp_pd, Kd_pd, ...
                       Kp_fbl_pid, Kd_fbl_pid, Ki_fbl_pid)
    n  = 6;
    q  = x(1:n);   dq = x(n+1:2*n);
    [qd, dqd, ddqd]          = joint_trajectory(t, traj);
    [y_ref, yd_ref, ydd_ref] = cartesian_reference(qd, dqd, ddqd, Pi);
    switch mode
        case 'FFW_PD'
            e  = qd - q;   de = dqd - dq;
            tau_ff = UR5_M(qd,Pi)*ddqd + UR5_h(qd,dqd,Pi) + UR5_G(qd,Pi);
            tau    = tau_ff + Kp_pd*e + Kd_pd*de;
            dx = [dq; UR5_fdyn(q,dq,tau,Pi)];
        case 'FBL_PD'
            tau = fbl_torque(q,dq,Pi,y_ref,yd_ref,ydd_ref, ...
                             Kp_fbl,Kd_fbl,zeros(6),zeros(6,1));
            dx  = [dq; UR5_fdyn(q,dq,tau,Pi)];
        case 'FBL_PID'
            eint = x(2*n+1:3*n);
            tau  = fbl_torque(q,dq,Pi,y_ref,yd_ref,ydd_ref, ...
                              Kp_fbl_pid,Kd_fbl_pid,Ki_fbl_pid,eint);
            ddq  = UR5_fdyn(q,dq,tau,Pi);
            T_cur = UR5_fkine(q,Pi);  T_ref = UR5_fkine(qd,Pi);
            ex    = task_error(T_cur(1:3,4),T_cur(1:3,1:3), ...
                               T_ref(1:3,4),T_ref(1:3,1:3));
            dx = [dq; ddq; ex];
    end
end

function tau = fbl_torque(q, dq, Pi, y_ref, yd_ref, ydd_ref, Kp, Kd, Ki, eint)
    T_cur = UR5_fkine(q, Pi);
    p_cur = T_cur(1:3,4);   R_cur = T_cur(1:3,1:3);
    [~,~,Jg] = UR5_jacobian_geometric(q, Pi);
    M = UR5_M(q,Pi);   h = UR5_h(q,dq,Pi);   G = UR5_G(q,Pi);
    y_cur  = [p_cur; rot2euler(R_cur)];
    dy_cur = Jg * dq;
    ex     = y_ref - y_cur;
    dt     = 1e-6;
    [~,~,Jg_p]  = UR5_jacobian_geometric(q + dq*dt, Pi);
    Jdot_dq = (Jg_p - Jg)/dt * dq;
    v   = ydd_ref + Kd*(yd_ref - dy_cur) + Kp*ex + Ki*eint;
    tau = M * (Jg \ (v - Jdot_dq)) + h + G;
end

function [y, dy, ddy] = cartesian_reference(qd, dqd, ddqd, Pi)
    T   = UR5_fkine(qd, Pi);
    p   = T(1:3,4);   R = T(1:3,1:3);
    [~,~,Jg] = UR5_jacobian_geometric(qd, Pi);
    dy  = Jg * dqd;
    dt  = 1e-6;
    [~,~,Jg_p] = UR5_jacobian_geometric(qd + dqd*dt, Pi);
    ddy = Jg*ddqd + (Jg_p - Jg)/dt * dqd;
    y   = [p; rot2euler(R)];
end

function [qd, dqd, ddqd] = joint_trajectory(t, traj)
    qd   = traj.q_bias + traj.q_amp .* sin(traj.w_ref.*t + traj.phi_ref);
    dqd  = traj.q_amp  .* traj.w_ref .* cos(traj.w_ref.*t + traj.phi_ref);
    ddqd = -traj.q_amp .* (traj.w_ref.^2) .* sin(traj.w_ref.*t + traj.phi_ref);
end

function ex = task_error(p_cur, R_cur, p_des, R_des)
    ep    = p_des - p_cur;
    R_err = R_des * R_cur.';
    sk    = 0.5*(R_err - R_err.');
    eo    = [sk(3,2); sk(1,3); sk(2,1)];
    ex    = [ep; eo];
end

function euler = rot2euler(R)
    psi   =  atan2(R(2,1), R(1,1));
    theta =  atan2(-R(3,1), sqrt(R(1,1)^2 + R(2,1)^2));
    phi   =  atan2(R(3,2), R(3,3));
    euler = [phi; theta; psi];
end

function sim = postprocess(sol, t_plot, Pi, mode, traj, ...
                            Kp_fbl, Kd_fbl, Kp_pd, Kd_pd, ...
                            Kp_fbl_pid, Kd_fbl_pid, Ki_fbl_pid)
    n = 6;
    X = deval(sol, t_plot);
    N = numel(t_plot);

    q_log = X(1:n,:);
    dq_log = X(n+1:2*n,:);
    tau_log = zeros(n,N);
    ex_log = zeros(6,N);
    pos_log = zeros(3,N);
    pos_ref = zeros(3,N);

    for k = 1:N
        q = q_log(:,k);
        dq = dq_log(:,k);
        t = t_plot(k);

        [qd, dqd, ddqd] = joint_trajectory(t, traj);
        [y_ref, yd_ref, ydd_ref] = cartesian_reference(qd, dqd, ddqd, Pi);
        pos_ref(:,k) = y_ref(1:3);

        T_cur = UR5_fkine(q, Pi);
        pos_log(:,k) = T_cur(1:3,4);
        ex_log(:,k) = task_error(T_cur(1:3,4), T_cur(1:3,1:3), y_ref(1:3), euler2rot(y_ref(4:6)));

        if strcmp(mode,'FBL_PID') && size(X,1) > 2*n
            eint = X(2*n+1:3*n,k);
        else
            eint = zeros(6,1);
        end

        switch mode
            case 'FFW_PD'
                e  = qd  - q;
                de = dqd - dq;
                tau_ff = UR5_M(qd,Pi)*ddqd + UR5_h(qd,dqd,Pi) + UR5_G(qd,Pi);
                tau_log(:,k) = tau_ff + Kp_pd*e + Kd_pd*de;
            case 'FBL_PD'
                tau_log(:,k) = fbl_torque(q,dq,Pi,y_ref,yd_ref,ydd_ref,Kp_fbl,Kd_fbl,zeros(6),eint);
            case 'FBL_PID'
                tau_log(:,k) = fbl_torque(q,dq,Pi,y_ref,yd_ref,ydd_ref,Kp_fbl_pid,Kd_fbl_pid,Ki_fbl_pid,eint);
        end
    end

    sim.t = t_plot;
    sim.q_log = q_log;
    sim.dq_log = dq_log;
    sim.tau_log = tau_log;
    sim.ex_log = ex_log;
    sim.pos_log = pos_log;
    sim.pos_ref = pos_ref;
    sim.pos_err = vecnorm(ex_log(1:3,:),2,1);
    sim.ori_err = vecnorm(ex_log(4:6,:),2,1);
    sim.tau_norm = vecnorm(tau_log,2,1);
end

function m = compute_metrics(sim)
    t = sim.t;
    dt = t(end)-t(1);
    m.rms_pos = sqrt(trapz(t, sim.pos_err.^2)/dt);
    m.rms_ori = sqrt(trapz(t, sim.ori_err.^2)/dt);
    m.pk_pos = max(sim.pos_err);
    m.pk_ori = max(sim.ori_err);
    m.energy = trapz(t, sim.tau_norm.^2);
end

function R = euler2rot(euler)
    phi = euler(1);
    theta = euler(2);
    psi = euler(3);
    Rz = [cos(psi) -sin(psi) 0; sin(psi) cos(psi) 0; 0 0 1];
    Ry = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
    Rx = [1 0 0; 0 cos(phi) -sin(phi); 0 sin(phi) cos(phi)];
    R = Rz*Ry*Rx;
end