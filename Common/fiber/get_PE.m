function[Fpe, kpe] = get_PE(Lce, parms)

kpe = 0;
Fpe = 0;


% if isfield(parms, 'PE_isw_SE') % PE in series with SE
%     if parms.PE_isw_SE
%         kpe = parms.kpe_func(Lce, parms);
%         Fpe = parms.Fpe_func(Lce, parms);
%     end
% end

if isfield(parms, 'PE_isw_SE') % PE in series with SE
    dLce = Lce - parms.Lce0;

    if parms.K*dLce < 10
        kpe = parms.kpe .* (1 - 1./(exp(parms.K*dLce)+1));
        Fpe = parms.kpe * log(1+exp(dLce*parms.K))/parms.K;
    else
        kpe = parms.kpe;
        Fpe = parms.kpe * dLce;
    end
end
end