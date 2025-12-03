function[Fpe, Lpe] = get_Fpe(t1, lce, t2, L, parms)

Lpe = L;

% if PE in series with SE, overwrite
if isfield(parms, 'PE_isw_SE')
    if parms.PE_isw_SE
        Lpe = interp1(t1, lce, t2);
    end
end

Fpe = parms.Fpe_func(Lpe, parms);

end