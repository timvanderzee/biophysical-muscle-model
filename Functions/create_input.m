function[tis, Cas, Lis, vis, ts, Ts] = create_input(tiso, dTt, dTc, ISI, Ca, N)

if nargin < 6
    N = 500;
    N = 200;
end

%% create idealized input
ti = linspace(0, tiso, N);
% vis = [];
% Lis = [];

vis = [];
Cas = [];
Lis = [];

for i = 1:length(dTc)

    Ts = cumsum([tiso-(dTt*3+dTc(i)*2+ISI(i)); dTc(i); dTc(i); ISI(i); dTt; dTt; dTt]);
    tx = linspace(0, Ts(end), length(Ca)*N);

    % specify velocity
    vx = nan(size(tx));
    vx(tx <= Ts(1)) = 0;
    vx(tx > Ts(1) & tx <=  Ts(2)) = .4545;
    vx(tx > Ts(2) & tx <=  Ts(3)) = -.4545;
    vx(tx > Ts(3) & tx <=  Ts(4)) = 0;
    vx(tx > Ts(4) & tx <=  Ts(5)) = .4545;
    vx(tx > Ts(5) & tx <=  Ts(6)) = 0;
    vx(tx > Ts(6) & tx <=  Ts(7)) = -.4545;

    % integrate to get length
    Lx = cumtrapz(tx, vx);

    % set back to 0
    vx(tx > Ts(6) & tx <=  Ts(7)) = 0;

    % interpolate
    vi = interp1(tx, vx, ti, [], 'extrap');
    Li = interp1(tx, Lx, ti, [], 'extrap');
    
%     vis = [vis vi];
%     Lis = [Lis Li];
    
%     figure(1)
%     plot(ti, Li); hold on
   

    %% combine different pCa trials

    for k = 1:length(Ca)
        vis = [vis vi];
        Cas = [Cas Ca(k) * ones(size(vi))];
        Lis = [Lis Li];
    end

    ts = 0:tiso:(tiso*(length(Ca)-1));
end

    tis = linspace(0, tiso * length(Ca) * length(dTc), length(ti) * length(Ca)* length(dTc));



end

