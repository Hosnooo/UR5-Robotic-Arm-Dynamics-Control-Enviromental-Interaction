clc; clear; close all;
 
% Apply global LaTeX + font formatting to all figures
figureoptscall;
 
outdir = 'q1_figures';
if ~exist(outdir, 'dir'), mkdir(outdir); end
 
%% Parameters
[~, Pi] = UR5_params();
n = 6;
 
%% Simulation settings 
Tf     = 20.0;
t_plot = linspace(0, Tf, 4001);
ode_opts = odeset('RelTol',1e-7,'AbsTol',1e-9,'MaxStep',5e-3);
 
%% Initial condition 
q0  = [0.4; -1.0; 0.9; -0.5; 0.4; -0.2];
dq0 = zeros(n,1);
 
%% Desired target 
q_des = [0.0; -0.7; 0.6; 0.0; 0.3; 0.0];
T_des = UR5_fkine(q_des, Pi);
p_des = T_des(1:3,4);
R_des = T_des(1:3,1:3);
 
%% Gains 
Kp   = diag([120; 120; 120;  50;  50;  50]);
Kd_A = diag([ 18;  18;  15;  10;   5;   4]);
Kd_B = diag([ 35;  35;  35;  18;  18;  18]);
 
%% Simulate
fprintf('Running Controller A ...\n');
rhs_A = @(t,x) rhs_ctrl(t,x,Pi,'A',p_des,R_des,Kp,Kd_A,Kd_B);
sol_A = ode45(rhs_A, [0 Tf], [q0;dq0], ode_opts);
 
fprintf('Running Controller B ...\n');
rhs_B = @(t,x) rhs_ctrl(t,x,Pi,'B',p_des,R_des,Kp,Kd_A,Kd_B);
sol_B = ode45(rhs_B, [0 Tf], [q0;dq0], ode_opts);
 
simA = postprocess(sol_A, t_plot, Pi, 'A', p_des, R_des, Kp, Kd_A, Kd_B);
simB = postprocess(sol_B, t_plot, Pi, 'B', p_des, R_des, Kp, Kd_A, Kd_B);
 
metricsA = compute_metrics(simA);
metricsB = compute_metrics(simB);
 
%% Plot styles
colA   = [0.12 0.47 0.71];   % blue
colB   = [0.89 0.10 0.11];   % red
lw     = 1.6;
styles = {'-','--'};
sims   = {simA, simB};
labels = {'Ctrl A: $P$-task, $D$-joint', 'Ctrl B: $PD$-task'};
 
%% 
%  FIGURE 1 — Overview

fig1 = figure('Color','w','Units','inches','Position',[1 1 8 7]);
 
tl = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
 
allcolor=[colA; colB];

% --- Position error ---
ax1 = nexttile; hold on; grid on; box on;
for k = 1:2
    plot(ax1, sims{k}.t, sims{k}.pos_err*1e3, ...
        styles{k}, 'LineWidth', lw, 'Color', allcolor(k,:));
end
yline(ax1, 0, 'k:', 'LineWidth', 1.0);
ylabel(ax1, 'Position error [mm]');
xlabel(ax1, 'Time [s]');
legend(ax1, labels, 'Location','northeast');
%title(ax1, 'End-Effector Position Error');
 
% --- Orientation error ---
ax2 = nexttile; hold on; grid on; box on;
for k = 1:2
    plot(ax2, sims{k}.t, rad2deg(sims{k}.ori_err), ...
        styles{k}, 'LineWidth', lw, 'Color', allcolor(k,:));
end
yline(ax2, 0, 'k:', 'LineWidth', 1.0);
ylabel(ax2, 'Orientation error [deg]');
xlabel(ax2, 'Time [s]');
legend(ax2, labels, 'Location','northeast');
%title(ax2, 'End-Effector Orientation Error');
 
% --- Torque norm ---
ax3 = nexttile; hold on; grid on; box on;
for k = 1:2
    plot(ax3, sims{k}.t, sims{k}.tau_norm, ...
        styles{k}, 'LineWidth', lw, 'Color', allcolor(k,:));
end
ylabel(ax3, '$\|\tau\|_2$ [Nm]');
xlabel(ax3, 'Time [s]');
legend(ax3, labels, 'Location','northeast');
%title(ax3, 'Joint Torque Norm');
 
%%title(tl, 'Question 1 --- Cartesian Regulation Controllers: Overview', ...
%%'FontSize', 13, 'FontWeight', 'bold');
 
saveFigureAsPDF(fig1, fullfile(outdir,'q1_fig1_overview.pdf'), 8, 7, 600);
 
%% 
%  FIGURE 2 — Task-space error components
fig2 = figure('Color','w','Units','inches','Position',[1 1 11 7]);
 
tl2 = tiledlayout(3,2,'TileSpacing','compact','Padding','compact');
 
pos_lbl = {'$e_x$','$e_y$','$e_z$'};
ori_lbl = {'$e_{\phi}$','$e_{\theta}$','$e_{\psi}$'};
 
for j = 1:3
    ax = nexttile(tl2); hold on; grid on; box on;
    for k = 1:2
        plot(ax, sims{k}.t, sims{k}.ex_log(j,:)*1e3, ...
            styles{k}, 'LineWidth', lw, 'Color', allcolor(k,:));
    end
    yline(ax, 0, 'k:', 'LineWidth', 1.0);
    ylabel(ax, [pos_lbl{j} ' [mm]']);
    xlabel(ax, 'Time [s]');
    if j == 1
        legend(ax, labels, 'Location','best');
    end
    grid(ax,'on');
end
 
for j = 1:3
    ax = nexttile(tl2); hold on; grid on; box on;
    for k = 1:2
        plot(ax, sims{k}.t, rad2deg(sims{k}.ex_log(3+j,:)), ...
            styles{k}, 'LineWidth', lw, 'Color', allcolor(k,:));
    end
    yline(ax, 0, 'k:', 'LineWidth', 1.0);
    ylabel(ax, [ori_lbl{j} ' [deg]']);
    xlabel(ax, 'Time [s]');
    grid(ax,'on');
end
 
%%title(tl2, 'Question 1 --- Task-Space Error Components', ...
%%'FontSize', 13, 'FontWeight', 'bold');
 
saveFigureAsPDF(fig2, fullfile(outdir,'q1_fig2_task_errors.pdf'), 11, 7, 600);
 
%% 
%  FIGURE 3 — Joint torques

fig3 = figure('Color','w','Units','inches','Position',[1 1 11 7]);
 
tl3 = tiledlayout(3,2,'TileSpacing','compact','Padding','compact');
 
for j = 1:6
    ax = nexttile(tl3); hold on; grid on; box on;
    for k = 1:2
        plot(ax, sims{k}.t, sims{k}.tau_log(j,:), ...
            styles{k}, 'LineWidth', lw, 'Color', allcolor(k,:));
    end
    ylabel(ax, sprintf('$\\tau_%d$ [Nm]', j));
    xlabel(ax, 'Time [s]');
    if j == 1
        legend(ax, labels, 'Location','best');
    end
end
 
%%title(tl3, 'Question 1 --- Joint Torques', ...
%%'FontSize', 13, 'FontWeight', 'bold');
 
saveFigureAsPDF(fig3, fullfile(outdir,'q1_fig3_joint_torques.pdf'), 11, 7, 600);
 
%% 
%  FIGURE 4 — 3D end-effector path

fig4 = figure('Color','w','Units','inches','Position',[1 1 7 6]);
 
ax4 = axes('Parent', fig4);
hold(ax4,'on'); grid(ax4,'on'); box(ax4,'on'); axis(ax4,'equal'); view(ax4,3);
 
for k = 1:2
    plot3(ax4, ...
        sims{k}.pos_log(1,:)*1e3, ...
        sims{k}.pos_log(2,:)*1e3, ...
        sims{k}.pos_log(3,:)*1e3, ...
        styles{k}, 'LineWidth', lw, 'Color', allcolor(k,:));
end
 
p0_val = sims{1}.pos_log(:,1)*1e3;
plot3(ax4, p0_val(1), p0_val(2), p0_val(3), ...
    'o', 'MarkerSize', 9, 'LineWidth', 2, ...
    'MarkerEdgeColor', [0.2 0.6 0.2], 'MarkerFaceColor', [0.2 0.8 0.2]);
plot3(ax4, p_des(1)*1e3, p_des(2)*1e3, p_des(3)*1e3, ...
    '*', 'MarkerSize', 11, 'LineWidth', 2, 'Color', [0.8 0.1 0.1]);
 
xlabel(ax4, '$x$ [mm]'); ylabel(ax4, '$y$ [mm]'); zlabel(ax4, '$z$ [mm]');
legend(ax4, [labels, {'Start','Target'}], 'Location','best');
%%title(ax4, 'Question 1 --- End-Effector Trajectory (3D)', ...
%%'FontSize', 13, 'FontWeight', 'bold');
 
saveFigureAsPDF(fig4, fullfile(outdir,'q1_fig4_ee_path.pdf'), 7, 6, 600);

 
%% 
%  LOCAL FUNCTIONS (same as run_cartesian_regulation.m)

function dx = rhs_ctrl(~, x, Pi, mode, p_des, R_des, Kp, Kd_A, Kd_B)
    n  = 6;
    q  = x(1:n);   dq = x(n+1:2*n);
    tau = compute_tau(q, dq, Pi, mode, p_des, R_des, Kp, Kd_A, Kd_B);
    ddq = UR5_fdyn(q, dq, tau, Pi);
    dx  = [dq; ddq];
end
 
function tau = compute_tau(q, dq, Pi, mode, p_des, R_des, Kp, Kd_A, Kd_B)
    T_cur = UR5_fkine(q, Pi);
    p_cur = T_cur(1:3,4);
    R_cur = T_cur(1:3,1:3);
    [~,~,Jg] = UR5_jacobian_geometric(q, Pi);
    G  = UR5_G(q, Pi);
    ex = task_error(p_cur, R_cur, p_des, R_des);
    switch mode
        case 'A'
            tau = Jg.' * (Kp * ex) - Kd_A * dq + G;
        case 'B'
            dx_task = Jg * dq;
            tau = Jg.' * (Kp * ex - Kd_B * dx_task) + G;
    end
end
 
function ex = task_error(p_cur, R_cur, p_des, R_des)
    ep    = p_des - p_cur;
    R_err = R_des * R_cur.';
    sk    = 0.5 * (R_err - R_err.');
    eo    = [sk(3,2); sk(1,3); sk(2,1)];
    ex    = [ep; eo];
end
 
function sim = postprocess(sol, t_plot, Pi, mode, p_des, R_des, Kp, Kd_A, Kd_B)
    n = 6;  X = deval(sol, t_plot);  N = numel(t_plot);
    q_log   = X(1:n,:);   dq_log = X(n+1:2*n,:);
    tau_log = zeros(n,N); ex_log = zeros(6,N); pos_log = zeros(3,N);
    for k = 1:N
        q = q_log(:,k); dq = dq_log(:,k);
        T_cur = UR5_fkine(q,Pi);
        pos_log(:,k)  = T_cur(1:3,4);
        ex_log(:,k)   = task_error(T_cur(1:3,4),T_cur(1:3,1:3),p_des,R_des);
        tau_log(:,k)  = compute_tau(q,dq,Pi,mode,p_des,R_des,Kp,Kd_A,Kd_B);
    end
    sim.t        = t_plot;
    sim.q_log    = q_log;   sim.dq_log  = dq_log;
    sim.tau_log  = tau_log; sim.ex_log  = ex_log;
    sim.pos_log  = pos_log;
    sim.pos_err  = vecnorm(ex_log(1:3,:),2,1);
    sim.ori_err  = vecnorm(ex_log(4:6,:),2,1);
    sim.tau_norm = vecnorm(tau_log,2,1);
end
 
function m = compute_metrics(sim)
    t = sim.t;  dt = t(end)-t(1);
    m.rms_pos   = sqrt(trapz(t,sim.pos_err.^2)/dt);
    m.rms_ori   = sqrt(trapz(t,sim.ori_err.^2)/dt);
    m.final_pos = sim.pos_err(end);
    m.final_ori = sim.ori_err(end);
    m.energy    = trapz(t, sim.tau_norm.^2);
    thresh = max(0.02*sim.pos_err(1), 1e-4);
    m.t_settle = NaN;
    for k = 1:numel(t)
        if all(sim.pos_err(k:end) <= thresh)
            m.t_settle = t(k); break;
        end
    end
end