function[Q0dot, Q1dot, Q2dot, Rdot] = CrossBridge_Dynamics(Q0, p, q, f, w, k1, k2, IGef, Non, DRX, IG, b, k, R, dLcrit, ps2)

% attaching
beta = f .* [1 0 w^2];

% detaching
c1 = [Q0; p; 2*q];
phi0 = -IGef{1}(c1,k1) -IGef{1}(c1,k2);  
phi1 = -IGef{2}(c1,k1) -IGef{2}(c1,k2);  
phi2 = -IGef{3}(c1,k1) -IGef{3}(c1,k2);  
phi = [phi0; phi1; phi2];

% forcible detachment
% if b > 0
%     gamma = b .*  [1 ps2 w^2 + ps2.^2];
    gamma = b .*  [1 p w^2 + p.^2];
%     gamma = b .*  [1 ps2 w^2 + 0^2];
%     phiR0 = -k * (IG{1}(inf, c1) -IG{1}(dLcrit, c1)) + gamma(1) * R;
%     phiR1 = -k * (IG{2}(inf, c1) -IG{2}(dLcrit, c1)) + gamma(2) * R;
%     phiR2 = -k * (IG{3}(inf, c1) -IG{3}(dLcrit, c1)) + gamma(3) * R;
    
    phiR0 = -k * (IG{1}(10, c1) -IG{1}(dLcrit, c1)) + gamma(1) * R;
    phiR1 = -k * (IG{2}(10, c1) -IG{2}(dLcrit, c1)) + gamma(2) * R;
    phiR2 = -k * (IG{3}(10, c1) -IG{3}(dLcrit, c1)) + gamma(3) * R;
    
    Rdot = -phiR0;

% else
%     Rdot = 0;
%     phiR0 = 0;
%     phiR1 = 0;
%     phiR2 = 0;
% end

phiR = [phiR0; phiR1; phiR2];

% sum
phiT = phi + phiR;
% phiT = phi;

% cross-bridge dynamics
Q0dot = DRX .* (beta(1,1) .* (Non-Q0)) + phiT(1,:);
Q1dot = DRX .* (beta(1,2) .* (Non-Q0)) + phiT(2,:);
Q2dot = DRX .* (beta(1,3) .* (Non-Q0)) + phiT(3,:);

end