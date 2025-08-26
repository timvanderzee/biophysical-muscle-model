function[error_length] = LengthEquilibrium(Q0, Q1, Fdot, Ld, vMtilda, kse0, kse1)

gamma = 108.3333; % length scaling
% Fdot = dQ0dt + dQ1dt;

% k = 100;
F = Q1 + Q0;
% F = log(1+exp(F*k))/k;
kse = kse1 * (F + kse0);

error_length = Ld  .* (Q0 + kse) - (vMtilda .* gamma .* kse - Fdot);

end