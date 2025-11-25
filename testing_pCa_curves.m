
close all
pCa = linspace(4.5, 9, 1000);
parms.n = 3;

pCa50 = 7;
parms.kappa = 10^-pCa50;

Ca = 10.^-pCa;

F = Ca.^parms.n./(parms.kappa^parms.n+Ca.^parms.n);

F2 = 1 ./ (1 + 10.^(parms.n * (pCa-pCa50)));

F3 = pCa.^parms.n ./ (pCa50^parms.n + pCa.^parms.n);

subplot(121)
semilogx(Ca, F); hold on

% subplot(122)
semilogx(Ca, F3, '--')