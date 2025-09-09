function[texp, Lexp, Fexp, Tsrel] = get_data(data, sISIs, sAMPs, spCas)

% all ISIs and AMPs
AMPs = [0 12 38 121 216 288 383 682]/10000;
ISIs = [1 10 100 316 1000 3160 10000]/1000;
pCas = [4.5000    6.1000    6.2000    6.3000    6.4000    6.6000    9.0000];

% velocity 
v = .4545; % L0/s
tiso = 3;

% XData = nan(2e4, length(sISIs));
Fexp = nan(2e4, length(sISIs));
texp = nan(2e4, length(sISIs));
Lexp = nan(2e4, length(sISIs));
Tsrel = nan(length(sISIs), 7);

for kk = 1:length(sISIs)
    ISI = sISIs(kk);
    AMP = sAMPs(kk);
    pCa = spCas(kk);
    
    n = find(ISIs == round(ISI,4));
    m = find(AMPs == round(AMP,4));
    i = find(pCas == round(pCa,1));

    if isempty(data.Fexp(:,i,n,m))
        keyboard
    end
    
     % selected data
    Fexp(:,kk) = data.Fexp(:,i,n,m);
    texp(:,kk) = data.texp(:,i,n,m);
    Lexp(:,kk) = data.Lexp(:,i,n,m);
    
    T = [AMP AMP .0383]./v; % time of stretch (s)
    Ts = [tiso T(1) T(2) ISI T(3) tiso];
    
    Tsrel(kk,:) = [-10 cumsum(Ts(1,:)) - sum(Ts(1,1:4))];

end
end