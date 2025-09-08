function [error, dX] = ThickEquilibrium(Q0, Q0dot, F, DRX, dDRXdt, k1, k2, kF, Ntot)

[J1, J2] = ThickFilament_Dynamics(Q0, F, DRX, k1, k2, kF, Ntot);

dX = J1 - J2 - Q0dot;

error = dDRXdt - dX; 

end