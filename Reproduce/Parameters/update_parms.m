function[parms] = update_parms(newparms)

parms = newparms;
%         parms.xi = linspace(-20,20,1000);
parms.xi = linspace(-15,15,250);
parms.f_func = @(xi,f,w,mu)   f/sqrt((2*pi*w^2))*exp(-(xi-mu).^2./(2*w^2));
parms.g_func = @(xi,k1,k2) k1*exp(k2*xi);
%         parms.K = 1000;
% parms.approx = 0;
parms.amin = 1e-3;
%         parms.k = 500;
%         parms.b = 1e4;
%         parms.ps2 = 1;
%         parms.dLcrit = 1.5;

parms.myofilaments.thick_filament_length = 1530/2;
parms.myofilaments.thin_filament_length = 1120;
parms.myofilaments.bare_zone_length = 80;

    parms.Fpe_func = @(L, parms) parms.kpe*(L-parms.lmtc0).*(L>parms.lmtc0)+parms.Fpe0;
    parms.kpe_func = @(Lce, parms) parms.kpe .*(Lce>parms.lmtc0);

end