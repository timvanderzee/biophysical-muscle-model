function [error, Fdot] = MuscleEquilibrium(Q0, p, q, dQ0dt, dQ1dt, dQ2dt, f, w, k11, k12, k21, k22, Non, Ld, DRX)

% points where integrals is evaluated
k1 = [k11 k12];
k2 = [k21 -k22];
% w = 0.2;

% gamma = 108.3333; % length scaling

% without limit
% IGef{1} = @(c,k)(c(1,:)*k(1).*exp(c(3,:)*k(2)^2/4-c(2,:)*k(2)));
% IGef{2} = @(c,k)(c(1,:)*k(1).*exp(c(3,:)*k(2)^2/4-c(2,:)*k(2))).*(c(2,:)-c(3,:)*k(2)/2);
% IGef{3} = @(c,k)(c(1,:)*k(1).*exp(c(3,:)*k(2)^2/4-c(2,:)*k(2))).*((c(2,:)-c(3,:)*k(2)/2).^2+c(3,:)/2);

% with limit
IGef{1} = @(c,k)(c(1,:)*k(1).*exp(10*tanh((c(3,:)*k(2)^2/4-c(2,:)*k(2))/10)));
IGef{2} = @(c,k)(c(1,:)*k(1).*exp(10*tanh((c(3,:)*k(2)^2/4-c(2,:)*k(2))/10))).*(c(2,:)-c(3,:)*k(2)/2);
IGef{3} = @(c,k)(c(1,:)*k(1).*exp(10*tanh((c(3,:)*k(2)^2/4-c(2,:)*k(2))/10))).*((c(2,:)-c(3,:)*k(2)/2).^2+c(3,:)/2);

% Compute first-order moment
Q1 = p .* Q0;

% Compute Qdot
[Q0dot, Q1dot, Q2dot] = CrossBridge_Dynamics(Q0, p, q, f, w, k1, k2, IGef, Non, DRX);

% Length dynamics (ignoring tendon)
% k = 100;
% F = Q1 + Q0;
% F = log(1+exp(F*k))/k;
% dlse = parms.Lse_func(F, parms);
% dlse = log(F/kse0+1)/kse;

% stiffness
% kse = parms.kse_func(dlse, parms);
% kse = kse1 * (F + kse0);

% velocity - independent derivative
Fdot  = Q1dot + Q0dot;

% velocity from force constraint
% Ld  = (vMtilda .* gamma .* kse - Fdot) ./ (Q0 + kse);
% Ld = vMtilda * gamma;

error_Q0 = dQ0dt - Q0dot;
error_Q1 = dQ1dt - (Q1dot + 1 * Ld .* Q0);
error_Q2 = dQ2dt - (Q2dot + 2 * Ld .* Q1);

% error_Q1 = dQ1dt - Q1dot;
% error_Q2 = dQ2dt - Q2dot;

error = [error_Q0; error_Q1; error_Q2];

end