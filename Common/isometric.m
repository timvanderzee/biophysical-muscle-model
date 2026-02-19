function[input] = isometric(Ca, T, dt, L0)

    input.t  = 0:dt:T;
    input.L  = zeros(size(input.t))  + L0;
    input.v  = zeros(size(input.t));
    input.Cas = ones(size(input.t)) * Ca;

end