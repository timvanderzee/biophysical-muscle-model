function[Data] = prep_data(data,n,m,Ks,tiso)

%% process data
% load short-range stiffness data (skinned rat soleus muscle fibers),
% cd('C:\Users\timvd\OneDrive - KU Leuven\9. Short-range stiffness\matlab\data')

Data.F = [];
Data.L = [];
Data.t = [];
Data.C = [];

Data.dTt = .0383/.4545; % test stretch (= constant)

for k = 1:length(Ks)
    
    F = data.Fexp(:,Ks(k),n,m);
    L = data.Lexp(:,Ks(k),n,m);
    t = data.texp(:,Ks(k),n,m);
    Ca = 10^(6-data.pCas(Ks(k))) * ones(size(F));
    
    Data.ISI = data.ISIs(n)/1000;
    Data.dTc = data.AMPs(m)/10000 / .4545; % conditioning stretch
    
    trange = [-.8 .2];
    id1 = t > trange(1) & t < trange(2);
    
    % save
    Data.t = [Data.t; t(id1) + k * tiso - 3*Data.dTt - .005]; % move to agree with idealized input
    Data.F = [Data.F; F(id1)];
    Data.L = [Data.L; L(id1)];
    Data.C = [Data.C; Ca(id1)];
    
end

Data.Ca = 10.^(6-data.pCas);

%% get velocity (through interpolating and filtering)
if ~isempty(Data.t)
    ti = linspace(0,max(Data.t), 10000);
    Li = interp1([0; Data.t], [0; Data.L], ti);
    
    % filter
    dt = mean(diff(ti),'omitnan');
    fs = 1/dt;
    Wn = 30 / (.5 * fs);
    [b,a] = butter(2, Wn);
    Lf = Li;
    Lf(isfinite(Li)) = filtfilt(b,a, Li(isfinite(Li)));
    
    % take time derivative
    vf = grad5(Lf(:), dt);
    
    % interpolate
    Data.v = interp1(ti, vf, Data.t);
else
    Data.v = [];
end



end
