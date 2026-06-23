function[input] = isokinetic(Ca, T, dt, L0, v, T0)

    input.t  = 0:dt:T;

    input.v  = v*ones(size(input.t)) .* (input.t > T0);


    input.L  = L0 + cumtrapz(input.t, input.v);
    input.Cas = ones(size(input.t)) * Ca;

end