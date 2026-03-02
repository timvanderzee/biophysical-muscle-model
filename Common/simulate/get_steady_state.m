function[Xss] = get_steady_state(model, odefunc, modelfunc, newparms, Cas)

X0 = get_initial_state(model, modelfunc, newparms);

newparms.vts = [0 0];
newparms.ti = [0 5];
newparms.Lts = [0 0];
newparms.Cas = Cas(1) * [1 1];

if contains(char(odefunc), 'i') % implicit
    sol = odefunc(@(t,y,yp) modelfunc(t,y,yp, newparms), [0 max(newparms.ti)], X0, zeros(size(X0)), []);
else % explicit
    sol = odefunc(@(t,y) modelfunc(t,y, newparms), [0 max(newparms.ti)], X0, []);
end

Xss = sol.y(:,end);

end

function[X0] = get_initial_state(model, modelfunc, newparms)

if contains(model, 'XB')
    
    % assumed initial cross-bridge (XB) state
    Q0 = 1e-3; % fraction of XBs bound
    p0 = 0; % mean strain of bound XBs (power-stroke centered and normalized)
    q0 = .1; % standard deviation strain of bound XBs (power-stroke normalized)
    
    % find the state in which Fse = Fce + Fpe, given the intial XB state
    [Q00, Q20, lce0, Q10, Fse0, Fpe0, Fce0] = find_steady_state(Q0, p0, q0, newparms, 'regular');
    
    if contains(char(modelfunc), 'discretized')
        costfunc = @(L, Fce, parms) (Fce + (newparms.kpe*(-L-parms.Lce0)) - (newparms.kse0*(exp(newparms.kse*L)-1))) .^2;
        dlse = fminsearch(@(L) costfunc(L, 0, newparms), 0);
        lce0 = -dlse;
        X0 = [Q0 * zeros(size(newparms.xi)) lce0 0 0 0];
    else
        X0 = [Q00 Q20 Fse0 0 0 0]';
    end
    
else
    X0 = 0;
    
end


if contains(model, '2-state')
    X0(end-1) = 1;
end

end