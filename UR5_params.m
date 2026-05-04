function [p, Pi] = UR5_params()
%UR5_PARAMS Define the parameter set and pack it into the Pi vector.
%
% Outputs
%   p  : structure containing geometry and dynamic parameters
%   Pi : column vector matching the parameter order expected by the
%        generated UR5_* functions
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
% Each row corresponds to one link: [rcx rcy rcz]
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