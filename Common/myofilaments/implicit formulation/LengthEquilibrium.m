function[error_length] = LengthEquilibrium(Q0, F, F0dot, Ld, vMtilda, kse0, kse1)

gamma = 108.3333; % length scaling
% F = Q1 + Q0;
kse = kse1 * (F + kse0);

error_length = Ld  .* (Q0 + kse) - (vMtilda .* gamma .* kse - F0dot);

end