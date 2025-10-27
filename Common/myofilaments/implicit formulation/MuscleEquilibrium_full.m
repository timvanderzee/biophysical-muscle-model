function [error_n, Qdot] = MuscleEquilibrium_full(n, dndt, Q0, Non, DRX, f, w, xi, k11, k12, k21, k22, f_func, g_func)

% safety
n(n<0) = 0;

% points where integrals is evaluated
k1 = [k11 k12];
k2 = [k21 -k22];

% attachment and detachment at each strain
beta = f_func(xi, f, w);
phi = -(g_func(xi, k1(1), -k1(2)) + g_func(xi, k2(1), -k2(2))) .* n';   

% change in cross-bridge attachment
ndot = DRX * (beta * (Non - Q0)) + phi;

% first determine contraction velocity
Qdot = trapz(xi(:), [ndot(:) xi(:).*ndot(:)]);
    
error_n = dndt - ndot(:);

end