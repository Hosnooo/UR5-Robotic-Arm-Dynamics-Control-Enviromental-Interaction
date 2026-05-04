function [M,C,G] = UR5_MCG(q,dq,Pi)
M = UR5_M(q,Pi);
C = UR5_C(q,dq,Pi);
G = UR5_G(q,Pi);
end
