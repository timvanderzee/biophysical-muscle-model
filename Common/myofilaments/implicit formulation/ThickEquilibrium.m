function [error, dX] = ThickEquilibrium(Q0, Q0dot, F, DRX, dDRXdt, k1, k2, kF, Ntot, R)

if nargin < 10
    R = 0;
end

[J1, J2] = ThickFilament_Dynamics(Q0, F, DRX, k1, k2, kF, Ntot, R);

dX = J1 - J2 - Q0dot;

error = dDRXdt - dX; 

end