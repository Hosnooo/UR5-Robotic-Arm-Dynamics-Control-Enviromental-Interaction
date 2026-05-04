function qdd = UR5_fdyn(q,dq,tau,Pi)
M = UR5_M(q,Pi);
h = UR5_h(q,dq,Pi);
G = UR5_G(q,Pi);
qdd = M \ (tau - h - G);
end
