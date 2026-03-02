function[vMtilda, Lts, Act] = get_input_from_parms(t, parms)
    
% Velocity
if numel(parms.vts) == 1
    vMtilda = parms.vts;
else
    vMtilda = interp1(parms.ti, parms.vts, t);
end

% Length
if numel(parms.Lts) == 1
    Lts = parms.Lts;
else
    Lts = interp1(parms.ti, parms.Lts, t);
end

% Calcium
if numel(parms.Cas) == 1
    Ca = parms.Cas;
else
    Ca = interp1(parms.ti, parms.Cas, t);
end

% Activation
Act = parms.actfunc(Ca, parms);

end