function [error_Q0, error_Q1, error_Q2, Fdot, Rdot] = MuscleEquilibrium(Q0, Q1, p, q, dQ0dt, dQ1dt, dQ2dt, f, w, k11, k12, k21, k22, Non, Ld, DRX, b, k, R, dLcrit, ps2, approx)

if nargin < 23
   approx = 1;
end

% points where integrals is evaluated
k1 = [k11 k12];
k2 = [k21 -k22];

% get functions
[IG, IGef] = get_IG_IGEf(approx);

% Compute Qdot
[Q0dot, Q1dot, Q2dot, Rdot] = CrossBridge_Dynamics(Q0, p, q, f, w, k1, k2, IGef, Non, DRX, IG, b, k, R, dLcrit, ps2);

% velocity - independent derivative
Fdot  = Q1dot + Q0dot;

error_Q0 = dQ0dt - Q0dot;
error_Q1 = dQ1dt - (Q1dot + 1 * Ld .* Q0);
error_Q2 = dQ2dt - (Q2dot + 2 * Ld .* Q1);

end