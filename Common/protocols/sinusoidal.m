function[input] = sinusoidal(Ca, T, dt, f, A, L0)

    input.t = 0:dt:T;
    input.L = A * (-.5 * cos(2*pi*f*input.t) + .5) + L0;
    input.v = A * 2*pi*f*.5*sin(2*pi*f*input.t);
    input.Cas = ones(size(input.t)) * Ca;

end