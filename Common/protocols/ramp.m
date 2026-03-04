function[input] = ramp(Ca, T, dt, ISI, A, v)

    t = 0:dt:T;
    N = length(t);

    dTt = .0383/v; % test stretch (= constant)
    dTc = A / v; % conditioning stretch

    [tis, Cas, Lis, vis] = create_input(T, dTt, dTc, ISI, Ca, N);

    input.t = tis;
    input.L = Lis;
    input.v = vis;
    input.Cas = Cas;
end