function [error, dX] = ThickEquilibrium(Q0, Q0dot, F, DRX, dDRXdt, k1, k2, kF, Ntot, R, Rdot)

if nargin < 10
    R = 0;
end

[J1, J2] = ThickFilament_Dynamics(Q0, F, DRX, k1, k2, kF, Ntot, R);

% note: the term "Q0dot + Rdot" is the part of Q0dot due to regular
% attachment. explanation: if Rdot is great, it means that Q0 is losing to
% R. this implies that for a given Q0dot, the part due to regular
% attachment must be greater. thus, we need to add Rdot to Q0dot
dX = J1 - J2 - (Q0dot + Rdot);

error = dDRXdt - dX; 

end