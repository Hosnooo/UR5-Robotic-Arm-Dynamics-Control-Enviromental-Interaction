function tau = UR5_idyn(q,dq,ddq,Pi)
M = UR5_M(q,Pi);
h = UR5_h(q,dq,Pi);
G = UR5_G(q,Pi);
tau = M*ddq + h + G;
end
