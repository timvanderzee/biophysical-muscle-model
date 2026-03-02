function[input] = ramp(Ca, T, dt, ISI, A)

    t = 0:dt:T;
    N = length(t);

    dTt = .0383/.4545; % test stretch (= constant)
    dTc = A / .4545; % conditioning stretch

    [tis, Cas, Lis, vis] = create_input(T, dTt, dTc, ISI, Ca, N);

    input.t = tis;
    input.L = Lis;
    input.v = vis;
    input.Cas = Cas;
end