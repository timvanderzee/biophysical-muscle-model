function[Fce, kpe, kse, Ntot] = get_force_and_stiffness(Fse, Lts, parms)

% SE
Fse(Fse<1e-6) = 1e-6;
dLse = max(parms.Lse_func(Fse, parms), 0); % can't be negative

% PE
Lce = Lts - dLse;
[Fpe, kpe] = get_PE(Lce, parms);

% CE
Fce = Fse - Fpe;
Fce = max(Fce, 0);

% SE stiffness
kse = parms.kse * (Fse + parms.kse0);

% overlap
if isfield(parms, 'FL_overlap')
    if parms.FL_overlap
        L = Lce/parms.gamma;
        Ntot = exp(-3*L.^2);
        
%         disp(Ntot)
%         Ntot = .5;
        
%         h = (.5 * parms.s / parms.gamma); % powerstroke size
%         L_hs = Lce * h * 1e9 + 1.3e3; % [nm]
%         Ntot = return_f_overlap(L_hs, parms);
    else
        Ntot = 1;
    end
else
    Ntot = 1;
end

end