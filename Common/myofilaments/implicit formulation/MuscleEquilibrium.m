function [error_Q0, error_Q1, error_Q2, error_R, Fdot] = MuscleEquilibrium(Q0, Q1, p, q, dQ0dt, dQ1dt, dQ2dt, f, w, k11, k12, k21, k22, Non, Ld, DRX, dRdt, b, k, R, dLcrit)

% points where integrals is evaluated
k1 = [k11 k12];
k2 = [k21 -k22];

% without limit
% IGef{1} = @(c,k)(c(1,:)*k(1).*exp(c(3,:)*k(2)^2/4-c(2,:)*k(2)));
% IGef{2} = @(c,k)(c(1,:)*k(1).*exp(c(3,:)*k(2)^2/4-c(2,:)*k(2))).*(c(2,:)-c(3,:)*k(2)/2);
% IGef{3} = @(c,k)(c(1,:)*k(1).*exp(c(3,:)*k(2)^2/4-c(2,:)*k(2))).*((c(2,:)-c(3,:)*k(2)/2).^2+c(3,:)/2);

% with limit
IGef{1} = @(c,k)(c(1,:)*k(1).*exp(10*tanh((c(3,:)*k(2)^2/4-c(2,:)*k(2))/10)));
IGef{2} = @(c,k)(c(1,:)*k(1).*exp(10*tanh((c(3,:)*k(2)^2/4-c(2,:)*k(2))/10))).*(c(2,:)-c(3,:)*k(2)/2);
IGef{3} = @(c,k)(c(1,:)*k(1).*exp(10*tanh((c(3,:)*k(2)^2/4-c(2,:)*k(2))/10))).*((c(2,:)-c(3,:)*k(2)/2).^2+c(3,:)/2);

IG{1} =  @(x,c)1/2*c(1)*erf((x-c(2))/sqrt(c(3)));
IG{2} =  @(x,c)1/2*c(1)*erf((x-c(2))/sqrt(c(3)))*c(2)-1/2*sqrt(c(3))*c(1)/sqrt(pi)*exp(-(x-c(2)).^2/c(3));
IG{3} =  @(x,c)1/2*c(1)*erf((x-c(2))/sqrt(c(3)))*(c(2)^2+c(3)/2)-1/2*sqrt(c(3))*c(1)/sqrt(pi)*exp(-(x-c(2)).^2/c(3)).*(c(2)+min(x,1e4));

% Compute Qdot
[Q0dot, Q1dot, Q2dot, Rdot] = CrossBridge_Dynamics(Q0, p, q, f, w, k1, k2, IGef, Non, DRX, IG, b, k, R, dLcrit);

% velocity - independent derivative
Fdot  = Q1dot + Q0dot;

error_Q0 = dQ0dt - Q0dot;
error_Q1 = dQ1dt - (Q1dot + 1 * Ld .* Q0);
error_Q2 = dQ2dt - (Q2dot + 2 * Ld .* Q1);
error_R = dRdt - Rdot;

end