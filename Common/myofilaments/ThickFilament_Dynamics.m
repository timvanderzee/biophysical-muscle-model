function[J1, J2] = ThickFilament_Dynamics(Q0, F, DRX, k1, k2, kF, Ntot, R)


SRX = 1 - Q0 - DRX - R;

J1 = k1 * (1 + kF * F/Ntot) .* SRX;
J2 = k2 .* DRX;
    
end