%% run_UR5_tree_vs_generated_test.m
% Compare rigidBodyTree against generated Maple/MATLAB UR5 files.
%
% Expected generated files:
%   UR5_fkine.m
%   UR5_fkall.m
%   UR5_jacobian_geometric.m
%   UR5_M.m
%   UR5_G.m
%   UR5_h.m
%   UR5_idyn.m
%
% This script:
%   1) builds rigidBodyTree from the SAME Pi vector
%   2) compares FK, FK-all, Jacobian, M, G, h, tau
%   3) also checks whether generated gravity has the opposite sign

clear; clc;

%% ----------------------- DRIVER -----------------------
addpath('maple_gen')

[p, Pi] = UR5_params();

nTests  = 50;
rngSeed = 1;

report = UR5_compare_tree_vs_generated(Pi, nTests, rngSeed);

disp(' ');
disp('Done.');
disp(report);

%% -------------------- LOCAL FUNCTIONS --------------------

function report = UR5_compare_tree_vs_generated(Pi, nTests, seed)

if nargin < 2 || isempty(nTests)
    nTests = 20;
end
if nargin < 3 || isempty(seed)
    seed = 1;
end

Pi = Pi(:);

requiredFiles = { ...
    'UR5_fkine.m', ...
    'UR5_fkall.m', ...
    'UR5_jacobian_geometric.m', ...
    'UR5_M.m', ...
    'UR5_G.m', ...
    'UR5_h.m', ...
    'UR5_idyn.m'};

for k = 1:numel(requiredFiles)
    assert(exist(requiredFiles{k}, 'file') == 2, ...
        'Missing required file: %s', requiredFiles{k});
end

rng(seed);
robot = UR5_rigidbodytree_from_Pi(Pi);

report = struct();
report.maxErr = struct( ...
    'fkine', 0, ...
    'fkall_T0e', 0, ...
    'fkall_each_body', zeros(6,1), ...
    'jacobian_linear', 0, ...
    'jacobian_angular', 0, ...
    'jacobian_full', 0, ...
    'M', 0, ...
    'G_raw', 0, ...
    'h', 0, ...
    'tau_raw', 0);

report.worstCase = struct();

for testIdx = 1:nTests
    q   = (2*rand(6,1) - 1) * pi;
    dq  = randn(6,1);
    ddq = randn(6,1);

    %% Generated model
    T_fkine_gen = UR5_fkine(q, Pi);
    [Tlist_gen, T0e_gen] = UR5_fkall(q, Pi);
    [Jv_gen, Jw_gen, Jg_gen] = UR5_jacobian_geometric(q, Pi);
    M_gen   = UR5_M(q, Pi);
    G_gen   = UR5_G(q, Pi);
    h_gen   = UR5_h(q, dq, Pi);
    tau_gen = UR5_idyn(q, dq, ddq, Pi);

    %% rigidBodyTree model
    T_tree = getTransform(robot, q, 'body6');

    % MATLAB returns geometricJacobian as [Jw; Jv]
    J_tree_raw = geometricJacobian(robot, q, 'body6');
    Jw_tree = J_tree_raw(1:3, :);
    Jv_tree = J_tree_raw(4:6, :);
    Jg_tree = [Jv_tree; Jw_tree];

    M_tree   = massMatrix(robot, q);
    G_tree   = gravityTorque(robot, q);
    h_tree   = velocityProduct(robot, q, dq);
    tau_tree = inverseDynamics(robot, q, dq, ddq);

    %% FK errors
    err_fkine = norm(T_tree - T_fkine_gen, 'fro');
    err_fkall = norm(T_tree - T0e_gen, 'fro');

    if err_fkine > report.maxErr.fkine
        report.maxErr.fkine = err_fkine;
        report.worstCase.fkine_q = q;
    end

    if err_fkall > report.maxErr.fkall_T0e
        report.maxErr.fkall_T0e = err_fkall;
        report.worstCase.fkall_q = q;
    end

    assert(iscell(Tlist_gen) && numel(Tlist_gen) >= 7, ...
        'Unexpected UR5_fkall output format.');

    for i = 1:6
        T_i_tree = getTransform(robot, q, sprintf('body%d', i));
        err_i = norm(T_i_tree - Tlist_gen{i+1}, 'fro');

        if err_i > report.maxErr.fkall_each_body(i)
            report.maxErr.fkall_each_body(i) = err_i;
            report.worstCase.(sprintf('body%d_q', i)) = q;
        end
    end

    %% Jacobian errors
    err_Jv = norm(Jv_tree - Jv_gen, 'fro');
    err_Jw = norm(Jw_tree - Jw_gen, 'fro');
    err_Jg = norm(Jg_tree - Jg_gen, 'fro');

    if err_Jv > report.maxErr.jacobian_linear
        report.maxErr.jacobian_linear = err_Jv;
        report.worstCase.Jv_q = q;
    end
    if err_Jw > report.maxErr.jacobian_angular
        report.maxErr.jacobian_angular = err_Jw;
        report.worstCase.Jw_q = q;
    end
    if err_Jg > report.maxErr.jacobian_full
        report.maxErr.jacobian_full = err_Jg;
        report.worstCase.Jg_q = q;
    end

    %% Dynamics errors
    err_M = norm(M_tree - M_gen, 'fro');
    err_h = norm(h_tree - h_gen, inf);

    % Raw generated gravity
    err_G_raw = norm(G_tree - G_gen, inf);

    % Raw generated inverse dynamics
    err_tau_raw = norm(tau_tree - tau_gen, inf);

    if err_M > report.maxErr.M
        report.maxErr.M = err_M;
        report.worstCase.M_q = q;
    end

    if err_h > report.maxErr.h
        report.maxErr.h = err_h;
        report.worstCase.h_qdq = [q; dq];
    end

    if err_G_raw > report.maxErr.G_raw
        report.maxErr.G_raw = err_G_raw;
        report.worstCase.G_raw_q = q;
    end

end

if err_tau_raw > report.maxErr.tau_raw
    report.maxErr.tau_raw = err_tau_raw;
    report.worstCase.tau_raw_qdqddq = [q; dq; ddq];
end

end

%% Summary
fprintf('\n=== rigidBodyTree vs generated files ===\n');
fprintf('Tests run                    : %d\n', nTests);
fprintf('max FK error (fkine)         : %.3e\n', report.maxErr.fkine);
fprintf('max FK error (fkall T0e)     : %.3e\n', report.maxErr.fkall_T0e);
for i = 1:6
    fprintf('max FK error (body%d)         : %.3e\n', i, report.maxErr.fkall_each_body(i));
end
fprintf('max Jv error                 : %.3e\n', report.maxErr.jacobian_linear);
fprintf('max Jw error                 : %.3e\n', report.maxErr.jacobian_angular);
fprintf('max Jg error                 : %.3e\n', report.maxErr.jacobian_full);
fprintf('max M error                  : %.3e\n', report.maxErr.M);
fprintf('max h error                  : %.3e\n', report.maxErr.h);
fprintf('max G error (raw)            : %.3e\n', report.maxErr.G_raw);
fprintf('max tau error (raw)          : %.3e\n', report.maxErr.tau_raw);

tol = 1e-10;

kinematicsGood = ...
    report.maxErr.fkine < tol && ...
    report.maxErr.fkall_T0e < tol && ...
    all(report.maxErr.fkall_each_body < tol) && ...
    report.maxErr.jacobian_linear < tol && ...
    report.maxErr.jacobian_angular < tol && ...
    report.maxErr.jacobian_full < tol;

if kinematicsGood
    fprintf('\nKinematics/Jacobian: PASS\n');
else
    fprintf('\nKinematics/Jacobian: MISMATCH\n');
end

if report.maxErr.M < tol && report.maxErr.h < tol
    fprintf('M and h: PASS\n');
else
    fprintf('M and h: MISMATCH\n');
end

if report.maxErr.G_raw < tol
    fprintf('G: PASS (raw)\n');
else
    fprintf('G: MISMATCH\n');
end

if report.maxErr.tau_raw < tol
    fprintf('tau: PASS (raw)\n');
else
    fprintf('tau: MISMATCH\n');
end


function robot = UR5_rigidbodytree_from_Pi(Pi)
% Build rigidBodyTree using the SAME Pi vector as the generated functions.
%
% IMPORTANT:
% - Pi inertia terms are assumed to be about link COM
% - rigidBody expects inertia about the body frame origin
% - so we convert using the parallel-axis theorem

Pi = Pi(:);
assert(numel(Pi) >= 67, 'Pi must have at least 67 elements.');

% Geometry
d1 = Pi(1);
a2 = Pi(2);
a3 = Pi(3);
d4 = Pi(4);
d5 = Pi(5);
d6 = Pi(6);

% Dynamics
m   = Pi(7:12);
rcx = Pi(13:18);
rcy = Pi(19:24);
rcz = Pi(25:30);
Ixx = Pi(31:36);
Iyy = Pi(37:42);
Izz = Pi(43:48);
Ixy = Pi(49:54);
Ixz = Pi(55:60);
Iyz = Pi(61:66);
g   = Pi(67);

robot = rigidBodyTree('DataFormat','column');
robot.Gravity = [0 0 -g];

% Maple model uses modified DH
% row i = [a_{i-1}, alpha_{i-1}, d_i, theta_i]
mdh = [ ...
    0    0      d1  0;
    0    pi/2   0   0;
    a2   0      0   0;
    a3   0      d4  0;
    0    pi/2   d5  0;
    0   -pi/2   d6  0];

parentName = robot.BaseName;

for i = 1:6
    body = rigidBody(sprintf('body%d', i));
    joint = rigidBodyJoint(sprintf('joint%d', i), 'revolute');
    setFixedTransform(joint, mdh(i,:), 'mdh');
    body.Joint = joint;

    r = [rcx(i); rcy(i); rcz(i)];

    % Inertia from Pi is about COM, expressed in body-frame orientation
    Ic_com = [ ...
        Ixx(i)  Ixy(i)  Ixz(i);
        Ixy(i)  Iyy(i)  Iyz(i);
        Ixz(i)  Iyz(i)  Izz(i)];

    % Parallel-axis theorem: inertia about body-frame origin
    I_body = Ic_com + m(i) * ((r.' * r) * eye(3) - (r * r.'));

    body.Mass = m(i);
    body.CenterOfMass = r.';

    % MATLAB order: [Ixx Iyy Izz Iyz Ixz Ixy]
    body.Inertia = [ ...
        I_body(1,1), ...
        I_body(2,2), ...
        I_body(3,3), ...
        I_body(2,3), ...
        I_body(1,3), ...
        I_body(1,2)];

    addBody(robot, body, parentName);
    parentName = body.Name;
end
end

function [p, Pi] = UR5_params()
%UR5_PARAMS Define the parameter set and pack it into the Pi vector.
%
% Pi order:
%   [d1; a2; a3; d4; d5; d6;
%    m1..m6;
%    rcx1..rcx6; rcy1..rcy6; rcz1..rcz6;
%    Ixx1..Ixx6; Iyy1..Iyy6; Izz1..Izz6;
%    Ixy1..Ixy6; Ixz1..Ixz6; Iyz1..Iyz6;
%    g]

% Geometry parameters
p.d1 = 0.0892;
p.a2 = -0.4250;
p.a3 = -0.39243;
p.d4 = 0.1090;
p.d5 = 0.0930;
p.d6 = 0.0820;

% Gravitational acceleration
p.g = 9.81;

% Link masses [kg]
p.m = [ ...
    3.7000;
    8.4000;
    2.3000;
    1.2000;
    1.2000;
    0.4000 ];

% Center-of-mass locations in link frames [m]
% Each row: [rcx rcy rcz]
p.rc = [ ...
    0.0000,  0.0000,  0.0300;
    -0.2100,  0.0000,  0.1200;
    -0.1800,  0.0000,  0.0200;
    0.0000,  0.0000,  0.0550;
    0.0000,  0.0000,  0.0400;
    0.0000,  0.0000,  0.0200 ];

% Inertia tensor diagonal terms about link COMs [kg*m^2]
p.Ixx = [0.0100; 0.0800; 0.0300; 0.0080; 0.0060; 0.0015];
p.Iyy = [0.0100; 0.0800; 0.0280; 0.0080; 0.0060; 0.0015];
p.Izz = [0.0060; 0.0200; 0.0150; 0.0040; 0.0030; 0.0010];

% Inertia tensor product terms about link COMs [kg*m^2]
p.Ixy = zeros(6,1);
p.Ixz = zeros(6,1);
p.Iyz = zeros(6,1);

% Pack parameter vector
Pi = [ ...
    p.d1; p.a2; p.a3; p.d4; p.d5; p.d6; ...
    p.m(:); ...
    p.rc(:,1); ...
    p.rc(:,2); ...
    p.rc(:,3); ...
    p.Ixx(:); ...
    p.Iyy(:); ...
    p.Izz(:); ...
    p.Ixy(:); ...
    p.Ixz(:); ...
    p.Iyz(:); ...
    p.g ];
end