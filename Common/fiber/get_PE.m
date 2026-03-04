function[Fpe, kpe] = get_PE(Lce, parms)

dLce = Lce - parms.Lce0;

if parms.K*dLce < 10
    kpe = parms.kpe .* (1 - 1./(exp(parms.K*dLce)+1));
    Fpe = parms.kpe * log(1+exp(dLce*parms.K))/parms.K;
else
    kpe = parms.kpe;
    Fpe = parms.kpe * dLce;
end

end