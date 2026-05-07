clc; clear; close all;

% Apply global LaTeX + font formatting to all figures
figureoptscall;

%% ---- Symbolic model ----------------------------------------
syms t1 t2 t3 t4 t5 t6 d1 d4 d5 d6 a2 a3 real
t=[t1;t2;t3;t4;t5;t6];
%{
DH parameters
i alpha_i-1 a_i-1 d_i theta_i
1 |  0  | 0  | d1 | t1
2 |  90 | 0  | 0  | t2
3 |  0  | a2 | 0  | t3
4 |  0  | a3 | d4 | t4
5 |  90 | 0  | d5 | t5
6 | -90 | 0  | d6 | t6
%}
%{
 T = [ cos(t),          -sin(t),         0,          a
       sin(t)*cos(a),   cos(t)*cos(a), -sin(a),  -sin(a)*d
       sin(t)*sin(a),   cos(t)*sin(a),  cos(a),   cos(a)*d
             0,               0,         0,          1    ]
%}
T01=[cos(t1) -sin(t1) 0 0; sin(t1) cos(t1) 0 0; 0 0 1 d1; 0 0 0 1];
T12=[cos(t2) -sin(t2) 0 0; 0 0 -1 0; sin(t2) cos(t2) 0 0; 0 0 0 1];
T23=[cos(t3) -sin(t3) 0 a2; sin(t3) cos(t3) 0 0; 0 0 1 0; 0 0 0 1];
T34=[cos(t4) -sin(t4) 0 a3; sin(t4) cos(t4) 0 0; 0 0 1 d4; 0 0 0 1];
T45=[cos(t5) -sin(t5) 0 0; 0 0 -1 -d5; sin(t5) cos(t5) 0 0; 0 0 0 1];
T56=[cos(t6) -sin(t6) 0 0; 0 0 1 d6; -sin(t6) -cos(t6) 0 0; 0 0 0 1];

T02=simplify(T01*T12);
T03=simplify(T02*T23);
T04=simplify(T03*T34);
T05=simplify(T04*T45);
T06=simplify(T05*T56);

Px=T06(1,4);
Py=T06(2,4);
Pz=T06(3,4);
P=[Px;Py;Pz];

J=sym(zeros(3,6));
for j=1:6
    J(1,j)=diff(P(1,1),t(j,1));
    J(2,j)=diff(P(2,1),t(j,1));
    J(3,j)=diff(P(3,1),t(j,1));
end

%% ---- Numeric model -----------------------------------------
DH_sym=[d1 d4 d5 d6 a2 a3];
DH_num=[0.089159 0.10915 0.09465 0.0823 0.425 0.39225];
J_num=subs(J,DH_sym,DH_num);
P_num=subs(P,DH_sym,DH_num);
T06_num=subs(T06,DH_sym,DH_num);

J_func=matlabFunction(J_num,"Vars",{t1,t2,t3,t4,t5,t6});
P_func=matlabFunction(P_num,"Vars",{t1,t2,t3,t4,t5,t6});
T06_func=matlabFunction(T06_num,"Vars",{t1,t2,t3,t4,t5,t6}); %#ok<NASGU>

%% ---- Labels ------------------------------------------------
lbl_std='Standard pinv';
lbl_dls='DLS pinv';
lbl_ns='Null-space manip.';
lbl_nc='Task aug. (no conflict)';
lbl_c='Task aug. (conflict)';

%% ---- Initial pose ------------------------------------------
q0=[0;0;pi/2;0;1.4;0]; % t5=1.4, pi/2 to show null space

P0=P_func(q0(1),q0(2),q0(3),q0(4),q0(5),q0(6));
fprintf("\nInitial EE position:\n");
fprintf("  x0 = %.4f m\n", P0(1));
fprintf("  y0 = %.4f m\n", P0(2));
fprintf("  z0 = %.4f m\n", P0(3));

JP0=J_func(q0(1),q0(2),q0(3),q0(4),q0(5),q0(6));
JP0=JP0(1:3,:);
w0=sqrt(abs(det(JP0 * JP0')));
fprintf("Initial manipulability: w0 = %.6f\n", w0);
fprintf("(w0 ~ 0 confirms wrist singularity at t5=1.4, DLS will activate)\n");

%% ---- Simulation parameters ---------------------------------
dt=0.01;
omega=0.5;
T_end=4*pi/omega;
t_vec=0:dt:T_end;
N=length(t_vec);

%% ---- Desired trajectory ------------------------------------
Ax=0.15;
Ay=0.10;
Az=0.10;
phi=pi/4;

x0=P0(1);
y0=P0(2);
z0=P0(3);

x_d=@(tk) x0+Ax*sin(omega*tk);
y_d=@(tk) y0+Ay*sin(2*omega*tk);
z_d=@(tk) z0+Az*sin(omega*tk+phi);
dx_d=@(tk) Ax*omega*cos(omega*tk);
dy_d=@(tk) 2*Ay*omega*cos(2*omega*tk);
dz_d=@(tk) Az*omega*cos(omega*tk+phi);

Kp=10*eye(3);

%% ---- Method 1: standard pinv and DLS pinv -----------------
lambda_max_m1=0.5;
epsilon=0.05;

q_std=q0;
Q_std=zeros(6,N);
E_std=zeros(3,N);
W_std=zeros(1,N);
P_actual_std=zeros(3,N);

fprintf("standard pseudo inverse start\n")

for k=1:N
    tk=t_vec(k);

    P_act_std=P_func(q_std(1),q_std(2),q_std(3),q_std(4),q_std(5),q_std(6));
    J_eval_std=J_func(q_std(1),q_std(2),q_std(3),q_std(4),q_std(5),q_std(6));
    JP_std=J_eval_std(1:3,:);

    xd_std=[x_d(tk);y_d(tk);z_d(tk)];
    dxd_std=[dx_d(tk);dy_d(tk);dz_d(tk)];

    e_std=xd_std-P_act_std;
    w_std=sqrt(abs(det(JP_std*JP_std')));

    JP_inv_std=(JP_std')/(JP_std*(JP_std'));
    dq_std=JP_inv_std*(dxd_std+Kp*e_std);

    q_std=q_std+dq_std*dt;

    Q_std(:,k)=q_std;
    E_std(:,k)=e_std;
    W_std(k)=w_std;
    P_actual_std(:,k)=P_act_std;
end

fprintf("standard pseudo inverse done, starting DLS pseudo inverse\n")

q_dls=q0;
Q_dls=zeros(6, N);
E_dls=zeros(3, N);
W_dls=zeros(1, N);
P_actual_dls=zeros(3, N);
Lam_hist_dls=zeros(1, N);

for k=1:N
    tk=t_vec(k);

    P_act_dls=P_func(q_dls(1),q_dls(2),q_dls(3),q_dls(4),q_dls(5),q_dls(6));
    J_eval_dls=J_func(q_dls(1),q_dls(2),q_dls(3),q_dls(4),q_dls(5),q_dls(6));
    JP_dls=J_eval_dls(1:3,:);

    xd_dls=[x_d(tk);y_d(tk);z_d(tk)];
    dxd_dls=[dx_d(tk);dy_d(tk);dz_d(tk)];

    e_dls=xd_dls-P_act_dls;
    w_dls=sqrt(abs(det(JP_dls*JP_dls')));

    if w_dls>=epsilon
        lam2_dls=0;
    else
        lam2_dls=lambda_max_m1^2*(1-w_dls/epsilon)^2;
    end

    JP_inv_dls=(JP_dls')/(JP_dls*(JP_dls')+lam2_dls*eye(3));
    dq_dls=JP_inv_dls*(dxd_dls+Kp*e_dls);

    q_dls=q_dls+dq_dls*dt;

    Q_dls(:,k)=q_dls;
    E_dls(:,k)=e_dls;
    W_dls(k)=w_dls;
    P_actual_dls(:,k)=P_act_dls;
    Lam_hist_dls(k)=sqrt(lam2_dls);
end

fprintf("DLS pseudo inverse done, starting repeatability test (3 cycles)\n")

t_rep=0:dt:3*T_end;
N_rep=length(t_rep);
q_rep=q0;
Q_rep=zeros(6, N_rep);

for k=1:N_rep
    tk=t_rep(k);
    P_act_rep=P_func(q_rep(1),q_rep(2),q_rep(3),q_rep(4),q_rep(5),q_rep(6));
    J_eval_rep=J_func(q_rep(1),q_rep(2),q_rep(3),q_rep(4),q_rep(5),q_rep(6));
    JP_rep=J_eval_rep(1:3,:);
    xd_rep=[x_d(tk);y_d(tk);z_d(tk)];
    dxd_rep=[dx_d(tk);dy_d(tk);dz_d(tk)];
    e_rep=xd_rep-P_act_rep;
    w_rep=sqrt(abs(det(JP_rep*JP_rep')));
    if w_rep>=epsilon
        lam2_rep=0;
    else
        lam2_rep=lambda_max_m1^2*(1-w_rep/epsilon)^2;
    end
    JP_inv_rep=JP_rep'/(JP_rep*JP_rep'+lam2_rep*eye(3));
    dq_rep=JP_inv_rep*(dxd_rep+Kp*e_rep);
    q_rep=q_rep+dq_rep*dt;
    Q_rep(:,k)=q_rep;
end

fprintf("Repeatability test done, starting null space method\n");

%% ---- Method 2: null-space manipulability ------------------
q_min=[-2*pi;-2*pi;-pi;-2*pi;-2*pi;-2*pi];
q_max=[2*pi;2*pi;pi;2*pi;2*pi;2*pi];
q_mid=(q_min+q_max)/2; %#ok<NASGU>

K0=0.5;
fd_delta=1e-4;
lambda_max_m2=5;

q_ns=q0;
Q_ns=zeros(6,N);
E_ns=zeros(3,N);
W_ns=zeros(1,N);
P_actual_ns=zeros(3,N);
NS_mag=zeros(1,N);

for k=1:N
    tk=t_vec(k);

    P_act_ns=P_func(q_ns(1),q_ns(2),q_ns(3),q_ns(4),q_ns(5),q_ns(6));
    J_eval_ns=J_func(q_ns(1),q_ns(2),q_ns(3),q_ns(4),q_ns(5),q_ns(6));
    JP_ns=J_eval_ns(1:3,:);

    xd_ns=[x_d(tk);y_d(tk);z_d(tk)];
    dxd_ns=[dx_d(tk);dy_d(tk);dz_d(tk)];

    e_ns=xd_ns-P_act_ns;
    w_ns=sqrt(abs(det(JP_ns*JP_ns')));

    if w_ns>=epsilon
        lam2_ns=0;
    else
        lam2_ns=lambda_max_m2^2*(1-w_ns/epsilon)^2;
    end

    JP_inv_ns=(JP_ns')/(JP_ns*(JP_ns')+lam2_ns*eye(3));
    P_null=eye(6)-JP_inv_ns*JP_ns;

    grad_H=grad_mani(q_ns,J_func,fd_delta);
    dq_ns_term=P_null*(K0*grad_H);

    dq_ns=JP_inv_ns*(dxd_ns+Kp*e_ns)+dq_ns_term;

    q_ns=q_ns+dq_ns*dt;

    Q_ns(:,k)=q_ns;
    E_ns(:,k)=e_ns;
    W_ns(k)=w_ns;
    P_actual_ns(:,k)=P_act_ns;
    NS_mag(k)=norm(dq_ns_term);
end

fprintf("null space method end, starting repeatability\n")

t_rep2=0:dt:3*T_end;
N_rep2=length(t_rep2);
q_rep2=q0;
Q_rep2=zeros(6, N_rep2);

for k=1:N_rep2
    tk=t_rep2(k);
    P_act_rep2=P_func(q_rep2(1),q_rep2(2),q_rep2(3),q_rep2(4),q_rep2(5),q_rep2(6));
    J_eval_rep2=J_func(q_rep2(1),q_rep2(2),q_rep2(3),q_rep2(4),q_rep2(5),q_rep2(6));
    JP_rep2=J_eval_rep2(1:3,:);
    xd_rep2=[x_d(tk);y_d(tk);z_d(tk)];
    dxd_rep2=[dx_d(tk);dy_d(tk);dz_d(tk)];
    e_rep2=xd_rep2-P_act_rep2;
    w_rep2=sqrt(abs(det(JP_rep2*JP_rep2')));
    if w_rep2>=epsilon
        lam2_rep2=0;
    else
        lam2_rep2=lambda_max_m2^2*(1-w_rep2/epsilon)^2;
    end
    JP_inv_rep2=JP_rep2'/(JP_rep2*JP_rep2'+lam2_rep2*eye(3));
    P_null=eye(6)-JP_inv_rep2*JP_rep2;
    grad_H_rep2=grad_mani(q_rep2,J_func,fd_delta);
    dq_rep2=JP_inv_rep2*(dxd_rep2+Kp*e_rep2)+P_null*(K0*grad_H_rep2);
    q_rep2=q_rep2+dq_rep2*dt;
    Q_rep2(:,k)=q_rep2;
end

fprintf("Repeatability test done\n")

%% ---- Method 3: task augmentation ---------------------------
z_const=P0(3);
Kp2=10;
Kp1=Kp;
lambda_max_m3=0.001;

x_d_nc=@(tk) P0(1) + Ax*sin(omega*tk);
y_d_nc=@(tk) P0(2) + Ay*sin(2*omega*tk);
z_d_nc=@(tk) z_const;
dx_d_nc=@(tk) Ax*omega*cos(omega*tk);
dy_d_nc=@(tk) 2*Ay*omega*cos(2*omega*tk);
dz_d_nc=@(tk) 0;

x_d_c=x_d;
y_d_c=y_d;
z_d_c=z_d;
dx_d_c=dx_d;
dy_d_c=dy_d;
dz_d_c=dz_d;

q_aug_nc=q0;
Q_aug_nc=zeros(6, N);
E1_aug_nc=zeros(3, N);
E2_aug_nc=zeros(1, N);
W_aug_nc=zeros(1, N);
P_actual_nc=zeros(3, N);
COND_nc=zeros(1, N);

for k=1:N
    tk=t_vec(k);

    P_act_nc=P_func(q_aug_nc(1),q_aug_nc(2),q_aug_nc(3),q_aug_nc(4),q_aug_nc(5),q_aug_nc(6));
    J_eval_nc=J_func(q_aug_nc(1),q_aug_nc(2),q_aug_nc(3),q_aug_nc(4),q_aug_nc(5),q_aug_nc(6));
    JP_nc=J_eval_nc(1:3,:);

    xd1=[x_d_nc(tk); y_d_nc(tk); z_d_nc(tk)];
    dxd1=[dx_d_nc(tk); dy_d_nc(tk); dz_d_nc(tk)];
    e1=xd1-P_act_nc;

    x2_actual=P_act_nc(3);
    x2_des=z_const;
    e2=x2_des-x2_actual;
    dx2_d=0;
    J2=JP_nc(3,:);

    Ja=[JP_nc;J2];
    COND_nc(k)=cond(Ja);

    Ja_inv=aug_dls_pinv(Ja,lambda_max_m3,epsilon);
    dxa=[dxd1 + Kp1*e1; dx2_d + Kp2*e2];

    dq=Ja_inv*dxa;
    w_nc=sqrt(abs(det(JP_nc*JP_nc')));

    q_aug_nc=q_aug_nc+dq*dt;

    Q_aug_nc(:,k)=q_aug_nc;
    E1_aug_nc(:,k)=e1;
    E2_aug_nc(k)=e2;
    W_aug_nc(k)=w_nc;
    P_actual_nc(:,k)=P_act_nc;
end

fprintf("no conflict done, starting with conflict\n")

q_aug_c=q0;
Q_aug_c=zeros(6, N);
E1_aug_c=zeros(3, N);
E2_aug_c=zeros(1, N);
W_aug_c=zeros(1, N);
P_actual_c=zeros(3, N);
COND_c=zeros(1, N);

for k=1:N
    tk=t_vec(k);

    P_act_c=P_func(q_aug_c(1),q_aug_c(2),q_aug_c(3),q_aug_c(4),q_aug_c(5),q_aug_c(6));
    J_eval_c=J_func(q_aug_c(1),q_aug_c(2),q_aug_c(3),q_aug_c(4),q_aug_c(5),q_aug_c(6));
    JP_c=J_eval_c(1:3,:);

    xd1_c=[x_d_c(tk);y_d_c(tk);z_d_c(tk)];
    dxd1_c=[dx_d_c(tk);dy_d_c(tk);dz_d_c(tk)];
    e1_c=xd1_c-P_act_c;

    x2_actual_c=P_act_c(3);
    e2_c=z_const-x2_actual_c;
    dx2_d_c=0;
    J2_c=JP_c(3,:);

    Ja_c=[JP_c;J2_c];
    COND_c(k)=cond(Ja_c);

    Ja_inv_c=aug_dls_pinv(Ja_c,lambda_max_m3,epsilon);
    dxa_c=[dxd1_c+Kp1*e1_c;dx2_d_c+Kp2*e2_c];

    dq_c=Ja_inv_c*dxa_c;
    w_c=sqrt(abs(det(JP_c*JP_c')));

    q_aug_c=q_aug_c+dq_c*dt;

    Q_aug_c(:,k)=q_aug_c;
    E1_aug_c(:,k)=e1_c;
    E2_aug_c(k)=e2_c;
    W_aug_c(k)=w_c;
    P_actual_c(:,k)=P_act_c;
end

fprintf("conflict done\n")

%% ---- Plot helper data --------------------------------------
xd_plot=arrayfun(@(tk) x_d(tk), t_vec);
yd_plot=arrayfun(@(tk) y_d(tk), t_vec);
zd_plot=arrayfun(@(tk) z_d(tk), t_vec);

%% FIGURE 1: 3D Cartesian Trajectory
figure('Name','Fig 1 - 3D Trajectory','NumberTitle','off');
plot3(xd_plot, yd_plot, zd_plot, 'k--','LineWidth',1.5); hold on;
plot3(P_actual_std(1,:), P_actual_std(2,:), P_actual_std(3,:), ...
      'b','LineWidth',1.5);
plot3(P_actual_dls(1,:), P_actual_dls(2,:), P_actual_dls(3,:), ...
      'g','LineWidth',1.5);
xlabel('$x$ [m]'); ylabel('$y$ [m]'); zlabel('$z$ [m]');
legend('$p_d$', lbl_std, lbl_dls, 'Location','best');
grid on; axis equal;

%% FIGURE 2: Position Error Norm
figure('Name','Fig 2 - Position Error','NumberTitle','off');
subplot(2,1,1);
plot(t_vec, vecnorm(E_std,2,1),'b','LineWidth',1.5);
xlabel('Time [s]'); ylabel('$\|e_p\|_2$ [m]');
grid on;

subplot(2,1,2);
plot(t_vec, vecnorm(E_dls,2,1),'g','LineWidth',1.5);
xlabel('Time [s]'); ylabel('$\|e_p\|_2$ [m]');
grid on;

%% FIGURE 3: Joint Trajectories
figure('Name','Fig 3 - Joint Trajectories','NumberTitle','off');
subplot(2,1,1);
plot(t_vec, Q_std','LineWidth',1.2);
xlabel('Time [s]'); ylabel('$q$ [rad]');
legend('$q_1$','$q_2$','$q_3$','$q_4$','$q_5$','$q_6$','Location','best');
grid on;

subplot(2,1,2);
plot(t_vec, Q_dls','LineWidth',1.2);
xlabel('Time [s]'); ylabel('$q$ [rad]');
legend('$q_1$','$q_2$','$q_3$','$q_4$','$q_5$','$q_6$','Location','best');
grid on;

%% FIGURE 4: Manipulability Index
figure('Name','Fig 4 - Manipulability','NumberTitle','off');
plot(t_vec, W_std,'b','LineWidth',1.5); hold on;
plot(t_vec, W_dls,'g','LineWidth',1.5);
yline(epsilon,'r--','LineWidth',1.5,'Label','\epsilon threshold');
xlabel('Time [s]'); ylabel('$w(q)$');
legend(lbl_std, lbl_dls, 'Location','best');
grid on;

%% FIGURE 5: DLS Damping Factor
figure('Name','Fig 5 - DLS Damping Factor','NumberTitle','off');
plot(t_vec, Lam_hist_dls,'k','LineWidth',1.5);
xlabel('Time [s]'); ylabel('$\lambda$');
grid on;

%% FIGURE 6: Repeatability Test
figure('Name','Fig 6 - Repeatability','NumberTitle','off');
plot(t_rep, Q_rep','LineWidth',1.2); hold on;
xline(T_end,   '--r','LineWidth',1.5,'Label','End Cycle 1');
xline(2*T_end, '--r','LineWidth',1.5,'Label','End Cycle 2');
xlabel('Time [s]'); ylabel('$q$ [rad]');
legend('$q_1$','$q_2$','$q_3$','$q_4$','$q_5$','$q_6$','Location','best');
grid on;

%% FIGURE 7: Cartesian Error Components
figure('Name','Fig 7 - Error Components','NumberTitle','off');
labels = {'$e_x$ [m]','$e_y$ [m]','$e_z$ [m]'};
for i = 1:3
    subplot(3,1,i);
    plot(t_vec, E_std(i,:),'b','LineWidth',1.2); hold on;
    plot(t_vec, E_dls(i,:),'g','LineWidth',1.2);
    ylabel(labels{i}); grid on;
    if i == 1
        legend(lbl_std, lbl_dls, 'Location','best');
    end
end
xlabel('Time [s]');

%% FIGURE 8: 3D Trajectory
figure('Name','Fig 8 - Method 2 3D Trajectory','NumberTitle','off');
plot3(xd_plot, yd_plot, zd_plot, 'k--', 'LineWidth', 1.5); hold on;
plot3(P_actual_dls(1,:), P_actual_dls(2,:), P_actual_dls(3,:), ...
      'b', 'LineWidth', 1.5);
plot3(P_actual_ns(1,:),  P_actual_ns(2,:),  P_actual_ns(3,:), ...
      'g', 'LineWidth', 1.5);
xlabel('$x$ [m]'); ylabel('$y$ [m]'); zlabel('$z$ [m]');
legend('$p_d$', lbl_dls, lbl_ns, 'Location','best');
grid on; axis equal;

%% FIGURE 9: Manipulability Comparison
figure('Name','Fig 9 - Manipulability','NumberTitle','off');
plot(t_vec, W_dls, 'b', 'LineWidth', 1.5); hold on;
plot(t_vec, W_ns,  'g', 'LineWidth', 1.5);
yline(epsilon, 'r--', 'LineWidth', 1.5, 'Label', '\epsilon threshold');
xlabel('Time [s]'); ylabel('$w(q)$');
legend(lbl_dls, lbl_ns, 'Location','best');
grid on;

%% FIGURE 10: Position Error
figure('Name','Fig 10 - Position Error','NumberTitle','off');
plot(t_vec, vecnorm(E_dls,2,1), 'b', 'LineWidth', 1.5); hold on;
plot(t_vec, vecnorm(E_ns, 2,1), 'g', 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('$\|e_p\|_2$ [m]');
legend(lbl_dls, lbl_ns, 'Location','best');
grid on;

%% FIGURE 11: Joint Trajectories
figure('Name','Fig 11 - Joint Trajectories','NumberTitle','off');
subplot(2,1,1);
plot(t_vec, Q_dls', 'LineWidth', 1.2);
xlabel('Time [s]'); ylabel('$q$ [rad]');
legend('$q_1$','$q_2$','$q_3$','$q_4$','$q_5$','$q_6$','Location','best');
grid on;

subplot(2,1,2);
plot(t_vec, Q_ns', 'LineWidth', 1.2);
xlabel('Time [s]'); ylabel('$q$ [rad]');
legend('$q_1$','$q_2$','$q_3$','$q_4$','$q_5$','$q_6$','Location','best');
grid on;

%% FIGURE 12: Null-Space Term Magnitude
figure('Name','Fig 12 - Null-Space Magnitude','NumberTitle','off');
plot(t_vec, NS_mag, 'k', 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('$\|\dot q_{ns}\|_2$ [rad/s]');
grid on;

%% FIGURE 13: Repeatability Comparison
figure('Name','Fig 13 - Repeatability Comparison','NumberTitle','off');
subplot(2,1,1);
plot(t_rep, Q_rep',  'LineWidth', 1.2); hold on;
xline(T_end,   '--r', 'LineWidth', 1.5);
xline(2*T_end, '--r', 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('$q$ [rad]');
legend('$q_1$','$q_2$','$q_3$','$q_4$','$q_5$','$q_6$','Location','best');
grid on;

subplot(2,1,2);
plot(t_rep2, Q_rep2', 'LineWidth', 1.2); hold on;
xline(T_end,   '--r', 'LineWidth', 1.5);
xline(2*T_end, '--r', 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('$q$ [rad]');
legend('$q_1$','$q_2$','$q_3$','$q_4$','$q_5$','$q_6$','Location','best');
grid on;

%% FIGURE 14: 3D Trajectory Comparison
figure('Name','Fig 14 - Method 3 3D Trajectory','NumberTitle','off');
subplot(1,2,1);
plot3(arrayfun(@(tk) x_d_nc(tk),t_vec), arrayfun(@(tk) y_d_nc(tk),t_vec), ...
      arrayfun(@(tk) z_d_nc(tk),t_vec), 'k--', 'LineWidth', 1.5); hold on;
plot3(P_actual_nc(1,:), P_actual_nc(2,:), P_actual_nc(3,:), ...
      'b', 'LineWidth', 1.5);
xlabel('$x$ [m]'); ylabel('$y$ [m]'); zlabel('$z$ [m]');
legend('$p_d$', lbl_nc, 'Location','best');
grid on; axis equal;

subplot(1,2,2);
plot3(xd_plot, yd_plot, zd_plot, 'k--', 'LineWidth', 1.5); hold on;
plot3(P_actual_c(1,:), P_actual_c(2,:), P_actual_c(3,:), ...
      'g', 'LineWidth', 1.5);
xlabel('$x$ [m]'); ylabel('$y$ [m]'); zlabel('$z$ [m]');
legend('$p_d$', lbl_c, 'Location','best');
grid on; axis equal;
sgtitle('Fig 14: Method 3 -- 3D Trajectory');

%% FIGURE 15: Primary Task Error
figure('Name','Fig 15 - Primary Task Error','NumberTitle','off');
subplot(2,1,1);
plot(t_vec, vecnorm(E1_aug_nc,2,1), 'b', 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('$\|e_1\|_2$ [m]');
grid on;

subplot(2,1,2);
plot(t_vec, vecnorm(E1_aug_c,2,1), 'g', 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('$\|e_1\|_2$ [m]');
grid on;
sgtitle('Fig 15: Primary Task Error');

%% FIGURE 16: Secondary Task Error
figure('Name','Fig 16 - Secondary Task Error','NumberTitle','off');
subplot(2,1,1);
plot(t_vec, E2_aug_nc, 'b', 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('$e_2$ [m]');
grid on;

subplot(2,1,2);
plot(t_vec, E2_aug_c, 'g', 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('$e_2$ [m]');
grid on;
sgtitle('Fig 16: Secondary Task Error');

%% FIGURE 17: Condition Number of Augmented Jacobian
figure('Name','Fig 17 - Condition Number','NumberTitle','off');
subplot(2,1,1);
semilogy(t_vec, COND_nc, 'b', 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('$\mathrm{cond}(J_a)$');
grid on;

subplot(2,1,2);
semilogy(t_vec, COND_c, 'g', 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('$\mathrm{cond}(J_a)$');
grid on;
sgtitle('Fig 17: Conflict Indicator -- cond($J_a$)');

%% FIGURE 18: Joint Trajectories
figure('Name','Fig 18 - Joint Trajectories Method 3','NumberTitle','off');
subplot(2,1,1);
plot(t_vec, Q_aug_nc', 'LineWidth', 1.2);
xlabel('Time [s]'); ylabel('$q$ [rad]');
legend('$q_1$','$q_2$','$q_3$','$q_4$','$q_5$','$q_6$','Location','best'); grid on;

subplot(2,1,2);
plot(t_vec, Q_aug_c', 'LineWidth', 1.2);
xlabel('Time [s]'); ylabel('$q$ [rad]');
legend('$q_1$','$q_2$','$q_3$','$q_4$','$q_5$','$q_6$','Location','best'); grid on;
sgtitle('Fig 18: Method 3 Joint Trajectories');

%% FIGURE 19: Z-Height Over Time
figure('Name','Fig 19 - Z Height','NumberTitle','off');
subplot(2,1,1);
plot(t_vec, P_actual_nc(3,:), 'b', 'LineWidth', 1.5); hold on;
yline(z_const, 'r--', 'LineWidth', 1.5, 'Label', 'z_{const}');
xlabel('Time [s]'); ylabel('$z$ [m]');
grid on;

subplot(2,1,2);
plot(t_vec, P_actual_c(3,:),  'g', 'LineWidth', 1.5); hold on;
plot(t_vec, arrayfun(@(tk) z_d_c(tk), t_vec), 'r--', 'LineWidth', 1.5);
yline(z_const, 'k:', 'LineWidth', 1.5, 'Label', 'z_{const}');
xlabel('Time [s]'); ylabel('$z$ [m]');
legend('Actual z','Desired z (Lissajous)','z const','Location','best');
grid on;
sgtitle('Fig 19: Z Height Tracking -- Conflict Demonstration');

%% FIGURE 20: Method Comparison Summary
figure('Name','Fig 20 - Full Method Comparison','NumberTitle','off');
plot(t_vec, vecnorm(E_dls,  2,1), 'b',  'LineWidth', 1.5); hold on;
plot(t_vec, vecnorm(E_ns,   2,1), 'g',  'LineWidth', 1.5);
plot(t_vec, vecnorm(E1_aug_nc,2,1),'m', 'LineWidth', 1.5);
plot(t_vec, vecnorm(E1_aug_c, 2,1),'r', 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('$\|e_p\|_2$ [m]');
legend(lbl_dls, lbl_ns, lbl_nc, lbl_c, 'Location','best');
grid on;

%% SAVE ALL FIGURES AS PDF
fprintf('Saving figures as PDF...\n');

output_folder = 'figures_pdf_sing_far';
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

fig_names = {
'Fig 1 - 3D Trajectory',              'Fig1_3D_Trajectory';
'Fig 2 - Position Error',             'Fig2_Position_Error';
'Fig 3 - Joint Trajectories',         'Fig3_Joint_Trajectories';
'Fig 4 - Manipulability',             'Fig4_Manipulability';
'Fig 5 - DLS Damping Factor',         'Fig5_DLS_Damping';
'Fig 6 - Repeatability',              'Fig6_Repeatability';
'Fig 7 - Error Components',           'Fig7_Error_Components';
'Fig 8 - Method 2 3D Trajectory',     'Fig8_Method2_3D';
'Fig 9 - Manipulability',             'Fig9_Manipulability_M1_M2';
'Fig 10 - Position Error',            'Fig10_Position_Error_M1_M2';
'Fig 11 - Joint Trajectories',        'Fig11_Joint_Trajectories_M1_M2';
'Fig 12 - Null-Space Magnitude',      'Fig12_NullSpace_Magnitude';
'Fig 13 - Repeatability Comparison',  'Fig13_Repeatability_Comparison';
'Fig 14 - Method 3 3D Trajectory',    'Fig14_Method3_3D';
'Fig 15 - Primary Task Error',        'Fig15_Primary_Error';
'Fig 16 - Secondary Task Error',      'Fig16_Secondary_Error';
'Fig 17 - Condition Number',          'Fig17_Condition_Number';
'Fig 18 - Joint Trajectories Method 3','Fig18_Joint_Traj_M3';
'Fig 19 - Z Height',                  'Fig19_Z_Height';
'Fig 20 - Full Method Comparison',    'Fig20_Full_Comparison';
};

for i = 1:size(fig_names, 1)
    fig_handles = findobj('Type', 'figure', 'Name', fig_names{i,1});
    if ~isempty(fig_handles)
        fig_handle = fig_handles(1);
        filepath = fullfile(output_folder, [fig_names{i,2} '.pdf']);
        saveFigureAsPDF(fig_handle, filepath);
        fprintf('  Saved: %s\n', filepath);
    else
        fprintf('  WARNING: Figure not found -- %s\n', fig_names{i,1});
    end
end
fprintf('All figures saved to folder: %s\n', output_folder);

%% ---- Local functions ---------------------------------------
function grad=grad_mani(q,J_func,fd_delta)
    n=length(q);
    grad=zeros(n,1);
    for i=1:n
        qp=q; qp(i)=qp(i)+fd_delta;
        qm=q; qm(i)=qm(i)-fd_delta;

        JP_p=J_func(qp(1),qp(2),qp(3),qp(4),qp(5),qp(6));
        JP_m=J_func(qm(1),qm(2),qm(3),qm(4),qm(5),qm(6));

        JP_p=JP_p(1:3,:);
        JP_m=JP_m(1:3,:);

        w_p=sqrt(abs(det(JP_p*(JP_p'))));
        w_m=sqrt(abs(det(JP_m*(JP_m'))));

        grad(i)=(w_p-w_m)/(2*fd_delta);
    end
end

function Ja_inv=aug_dls_pinv(Ja,lambda_max,epsilon)
    w_a=sqrt(abs(det(Ja*Ja')));
    if w_a>=epsilon
        lam2=0;
    else
        lam2=lambda_max^2*(1-w_a/epsilon)^2;
    end
    Ja_inv=Ja'/(Ja*Ja'+lam2*eye(4));
end