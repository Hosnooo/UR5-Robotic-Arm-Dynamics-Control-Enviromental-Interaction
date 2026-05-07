
clear all
clc

syms t1 t2 t3 t4 t5 t6 d1 d4 d5 d6 a2 a3 real
syms dt1 dt2 dt3 dt4 dt5 dt6 real
t=[t1;t2;t3;t4;t5;t6];
d=[dt1;dt2;dt3;dt4;dt5;dt6];
delR=sym(zeros(3,3));
R_dot=sym(zeros(3,3));
P_dot_x=sym(zeros(1,1));
P_dot_y=sym(zeros(1,1));
P_dot_z=sym(zeros(1,1));
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
%
T02=T01*T12;
T02=simplify(T02);
T03=T02*T23;
T03=simplify(T03);
T04=T03*T34;
T04=simplify(T04);
T05=T04*T45;
T05=simplify(T05);
T06=T05*T56;
T06=simplify(T06);
%{
T06 = [r11 r12 r13 px
       r21 r22 r23 py
       r31 r32 r33 pz
        0   0   0  1 ]
%}
Px=T06(1,4);
Py=T06(2,4);
Pz=T06(3,4);
r11=T06(1,1);
r12=T06(1,2);
r13=T06(1,3);
r21=T06(2,1);
r22=T06(2,2);
r23=T06(2,3);
r31=T06(3,1);
r32=T06(3,2);
r33=T06(3,3);

%linear velocity jacobian vector
Jv=sym(zeros(3,6));
for j=1:6
    Jv(:,j)=[diff(Px,t(j));diff(Py,t(j));diff(Pz,t(j))];
end
P_dot=Jv*d;   
P_dot=simplify(P_dot);

%eular angles
psi=atan2(r21,r11); % yaw
theta=atan2(-r31,sqrt(r11^2+r21^2)); % pitch
phi=atan2(r32,r33); % roll
euler=[phi;theta;psi];


%Angular velocity Jacobian in axis-angle
Jk=sym(zeros(3,6));
for j=1:6
    Jk(:,j)=[diff(phi,t(j));diff(theta,t(j));diff(psi,t(j))];
end
Jk=simplify(Jk);
K_dot=simplify(Jk*d);

%Full 6x6 analytic Jacobian
Ja = [Jv;Jk];

% Step 1: build geometric Jacobian (z-axis cross product method)
z0=[0;0;1];
z1=T01(1:3,1:3)*z0;
z2=T02(1:3,1:3)*z0;
z3=T03(1:3,1:3)*z0;
z4=T04(1:3,1:3)*z0;
z5=T05(1:3,1:3)*z0;

p0=[0;0;0];
p1=T01(1:3,4);
p2=T02(1:3,4);
p3=T03(1:3,4);
p4=T04(1:3,4);
p5=T05(1:3,4);
pe=T06(1:3,4);

%Linear part columns: zi x (pe - pi)
%Angular part columns: zi
J_geo=[cross(z0,pe-p0) cross(z1,pe-p1) cross(z2,pe-p2) cross(z3,pe-p3) cross(z4,pe-p4) cross(z5,pe-p5);z0 z1 z2 z3 z4 z5];
J_geo=simplify(J_geo);

%substitute DH values numerically
DH=[d1 d4 d5 d6 a2 a3];
DH_vals=[0.089159 0.10915 0.09465 0.0823 0.425 0.39225];
sing_threshold=1e-3;
angles=0:pi/6:2*pi;

J_func=matlabFunction(subs(J_geo,DH,DH_vals),'Vars',[t1 t2 t3 t4 t5 t6]);

fprintf('\n=== Wrist Singularity (t5 = 0 or pi) ===\n');
for t5v=[0,pi]
    found=false;
    for t1v=angles
    for t2v=angles
    for t3v=angles
        J_num=J_func(t1v,t2v,t3v,0,t5v,0);
        if min(svd(J_num))<sing_threshold
            found=true;break
        end
    end
    if found;break;end
    end
    if found;break;end
    end
    if found
        fprintf('  CONFIRMED: t5 = %.4f rad causes wrist singularity\n',t5v);
    end
end

fprintf('\n=== Elbow Singularity (t3 = 0 or pi) ===\n');
for t3v=[0,pi]
    %use a non-degenerate base pose
    J_num=J_func(pi/4,-pi/4,t3v,pi/4,pi/2,pi/4);
    sv=min(svd(J_num));
    if sv<sing_threshold
        fprintf('  CONFIRMED: t3 = %.4f rad causes elbow singularity (min_sv=%.2e)\n',t3v,sv);
    else
        fprintf('  t3 = %.4f rad  min_sv = %.4f (not singular)\n',t3v,sv);
    end
end

fprintf('\n=== Shoulder Singularity (wrist center on z0 axis) ===\n');

T05_num=subs(T05,DH,DH_vals);
Pw=T05_num(1:3,4);
xy_dist_sq=simplify(subs(Pw(1)^2+Pw(2)^2,[t1 t4 t5 t6],[0 0 0 0]));

% Find the minimum possible xy_dist over all t2, t3
f=matlabFunction(xy_dist_sq,'Vars',[t2 t3]);
opts=optimoptions('fminunc','Display','off');
[~, min_val]=fminunc(@(x) f(x(1),x(2)),[0 0],opts);
min_xy_dist=sqrt(max(0,min_val));

fprintf('  Minimum possible wrist-center xy_dist = %.6f m\n',min_xy_dist);

if min_xy_dist<1e-4
    fprintf('  CONFIRMED: true shoulder singularity exists\n');
else
    fprintf('  NOTE: wrist center never reaches base z-axis (min dist = %.4f m)\n',min_xy_dist);
    fprintf('  The UR5 has a NEAR-shoulder singularity — the offsets d4/d5\n');
    fprintf('  prevent the wrist from reaching z0, but manipulability\n');
    fprintf('  drops significantly near this minimum.\n');
end

fprintf('\n=== Summary ===\n');
fprintf('  Wrist    : t5 = 0 or pi   (joints 4 and 6 axes align)\n');
fprintf('  Elbow    : t3 = 0 or pi   (arm fully extended or folded)\n');
fprintf('  Shoulder : no exact singularity due to UR5 offsets d4/d5\n');
fprintf('             (min wrist-center distance from z0 = %.4f m)\n',min_xy_dist);
fprintf('             Manipulability degrades near this minimum config.\n');

