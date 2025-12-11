function[IG] = get_initial_guess(tis, Cas, vts, Lts, parms)

% obtain initial guess
% intial guess is obtained through running a forward simulation with the
% first, simulate an isometric contraction


if parms.f > 0
    clc
    tic
    parms.vts = [0 0];
    parms.ti = [0 5];
    parms.Lts = [0 0];
    parms.Cas = Cas(1) * [1 1];
    
    parms.PE_isw_SE = 0;
           
    Q0 = 0;
    p0 = 0;
    q0 = .1;

    [Q00, Q20, lce0, Q10, Fse0] = find_steady_state(Q0, p0, q0, parms, 'adjusted');
    x0 = [Q00 Q20 Fse0 0 0 0]';
    xp0 = zeros(size(x0));
    
    if parms.J1 == 0
        x0(6) = 1;
    end
    
    odeopt = odeset('maxstep', 1e-2);
%     odeopt = [];
    
    % isometric contraction
    sol0 = ode15i(@(t,y,yp) fiber_dynamics_implicit_length_v2(t,y,yp, parms), [0 5], x0, xp0, odeopt);
    dx = fiber_dynamics_explicit_length_v2(sol0.x(end), sol0.y(:,end), parms);


    % dynamic contraction
    parms.vts = vts;
    parms.Lts = Lts;
%     parms.Lts = cumtrapz(tis, vts * parms.gamma);
    parms.ti = tis;
    parms.Cas = Cas;
    sol = ode15i(@(t,y,yp) fiber_dynamics_implicit_length_v2(t,y,yp,parms), [0 max(tis)], sol0.y(:,end), dx, odeopt);
%     sol = ode15s(@(t,y) fiber_dynamics_explicit_length(t,y,parms), [0 max(tis)], sol0.y(:,end), odeopt);
    [~, xdot] = deval(sol, sol.x);
    toc
    
    %% evaluate the errors in the derivative function
%     for i = 1:length(sol.x)
%         dxerrors(:,i) = fiber_dynamics_implicit_length(sol.x(i), sol.y(:,i), xdot(:,i), parms);
%     end
    
    %% extract the output and interpolate   
    % states
    Q0i     = interp1(sol.x, sol.y(1,:), tis); % zero-order moment
    Q2i     = interp1(sol.x, sol.y(2,:), tis); % second-order moment
    Fsei    = interp1(sol.x, sol.y(3,:), tis); % length
    Noni    = interp1(sol.x, sol.y(4,:), tis); % thin filament activation
    DRXi    = interp1(sol.x, sol.y(5,:), tis); % thick filament activation
    Ri      = interp1(sol.x, sol.y(6,:), tis); % 
%     Lti     = interp1(sol.x, sol.y(7,:), tis); %  
    
    % state derivatives
    dQ0dti  = interp1(sol.x, xdot(1,:), tis); % zero-order moment time derivative
    dQ2dti  = interp1(sol.x, xdot(2,:), tis); % second-order moment time derivative
    Fdi     = interp1(sol.x, xdot(3,:), tis); % velocity
    dNondti = interp1(sol.x, xdot(4,:), tis); % thin filament activation time derivative
    dDRXdti = interp1(sol.x, xdot(5,:), tis); % thick filament activation time derivative
    dRdti   = interp1(sol.x, xdot(6,:), tis);
%     vti     = interp1(sol.x, xdot(7,:), tis); %  
    
    % state-dependent variables
    dlse = parms.Lse_func(Fsei, parms);
    Li = Lts - dlse;
%     Fpe = parms.Fpe_func(Li, parms);
    Fpe = 0;
    Fce = Fsei - Fpe;
    Fce = max(Fce,0);
    kse = parms.kse_func(dlse, parms);
    kpe = parms.kpe_func(Li, parms);
%     kpe = parms.kpe .* (Li > 0);

    % mean and standard deviation
%     K = parms.K;
%     Q00i = log(1+exp(Q0i*K))/K;
    Q00i = max(Q0i, 1e-6);
    Q1i = Fce - Q00i;
    pi = Q1i./Q00i; 
    qi = Q2i./Q00i - pi.^2;  
    % q = log(1+exp(q*k))/k;
    qi = max(qi, 1e-6);

    
    dQ1dti = 0;
%     error_Q1i = 0;
    
    % get initial errors
%     error_thini      = ThinEquilibrium(parms.Cas, Q0i, Noni, dNondti, parms.kon, parms.koff, parms.koop, parms.Noverlap); % thin filament dynamics
%     error_thicki     = ThickEquilibrium(Q0i, dQ0dti, Fce, DRXi, dDRXdti, parms.J1, parms.J2, parms.JF, parms.Noverlap, Ri, dRdti); % thick filament dynamics
%     [error_Q0i, ~, error_Q2i, F0dot, Rdot] = MuscleEquilibrium(Q0i, Q1i, pi, qi, dQ0dti, dQ1dti, dQ2dti, parms.f, parms.w, parms.k11, parms.k12, parms.k21, parms.k22,  Noni, Ldi, DRXi, parms.b, parms.k, Ri, parms.dLcrit, parms.ps2, parms.approx); % cross-bridge dynamics
%     error_Ri = dRdti - Rdot;
%     error_lengthi    = LengthEquilibrium(Q0i, F0dot, Ldi, parms.vts, kse, kpe, parms.gamma);
        
    % errors
%     IG.error = [error_thini; error_thicki; error_Q0i; error_Q2i; error_Ri; error_lengthi];
%     
%     figure(1)
%     plot(IG.error')

    %% save to struct
    % states
    IG.Q0i = Q0i;
    IG.Q2i = Q2i;
    IG.Fsei = Fsei;
    IG.Noni = Noni;
    IG.DRXi = DRXi;
    IG.Ri = Ri;
%     IG.Lti = Lti;

%% compute Ldi
    kpe = 0;
    
    % points where integrals is evaluated
    k1 = [parms.k11 parms.k12];
    k2 = [parms.k21 -parms.k22];

    % get functions
    [IGf, IGef] = get_IG_IGEf(parms.approx);

    [Q0dot, Q1dot] = CrossBridge_Dynamics(Q0i, pi, qi, parms.f, parms.w, k1, k2, IGef, Noni, DRXi, IGf, parms.b, parms.k, Ri, parms.dLcrit, parms.ps2);

    % velocity - independent derivative
    F0dot  = Q1dot + Q0dot;
    
    kse = parms.kse * (Fsei + parms.kse0);
    Ldi = (vts .* parms.gamma .* kse - F0dot) ./  (Q0i + kse + kpe);

    %%
    % state derivatives
    IG.dQ0dti = dQ0dti;
    IG.dQ2dti = dQ2dti;
    IG.Ldi = Ldi;
    
    IG.dNondti = dNondti;
    IG.dDRXdti = dDRXdti;
    IG.dRdti = dRdti;
%     IG.vti = vti;
    
    % other variables
    IG.Q1i = Q1i;
    IG.Fcei = Fce;
%     IG.Fsei = Fsei;
    IG.Fpei = Fpe;
    IG.dQ1dti = dQ1dti;
%     IG.F0doti = F0dot;
    IG.pi = pi;
    IG.qi = qi;
    IG.Q00i = Q00i;

else
    parms.vts = 0;
    parms.ti = 0;
    parms.Cas = Cas(1);
    parms.Lts = 0;
    
    x0 =  0;
    xp0 = 0;
    odeopt = odeset('maxstep', 3e-3);
    
    % simulate
    sol0 = ode15i(@(t,y,yp) hill_type_implicit_v2(t,y,yp, parms), [0 max(tis)], x0, xp0, odeopt);
    
    % next, simulate response to specified velocity input vector
    parms.vts = vts;
    parms.ti = tis;
    parms.Cas = Cas;
    parms.Lts = Lts;
    
    % simulate
    sol = ode15i(@(t,y,yp) hill_type_implicit_v2(t,y,yp, parms), [0 max(tis)], sol0.y(end), xp0, odeopt);
    [~,xdot] = deval(sol, sol.x);
    
    %%
    % toc = sol.x;
    
    % interpolate solution to time nodes
    Li     = interp1(sol.x, sol.y(1,:), tis); % zero-order moment
    vi     = interp1(sol.x, xdot(1,:), tis); % zero-order moment time derivative
    
    Fce_rel = parms.vF_func(vi, parms);
    
    % activation from Ca
    a = parms.actfunc(Cas, parms);
    a(a<parms.amin) = parms.amin;
    
    Fce = a .* Fce_rel;
    
    % elastic elements
    Lse = parms.Lts - Li;
    
    Fse = parms.Fse_func(Lse, parms);
    Fse(Lse < 0) = 0;
    
    % error terms
    IG.error = Fse - Fce;
    IG.Fi = Fce;
    
    IG.vi = vi;
    IG.Li = Li;
    % IG.Firel = Fce_rel;
    
    %%
    % activation from Ca
    Act = parms.act_max * Cas.^parms.n ./ (parms.kappa^parms.n + Cas.^parms.n);
    
    Frel = Fce ./ Act;
    
    vt = parms.vmax/(2*parms.e(2))*(-exp(parms.e(4)/parms.e(1)-Frel/parms.e(1))+exp(Frel/parms.e(1)-parms.e(4)/parms.e(1))-2*parms.e(3));
    
    dL = log(Fce./parms.kse0 + 1) / parms.kse;
    
    Lt = parms.Lts - dL;
    
    % error terms
    error1 = vt - vi;
    error2 = Lt - Li;
    
    %     opti.subject_to((v(1:N-1) + v(2:N))*dt/2 + L(1:N-1) == L(2:N));
    
end

end
