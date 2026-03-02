function[newparms] = complete_parms(redparms)

newparms = redparms;

% attachment functions
newparms.f_func = @(xi,f,w,mu) f/sqrt((2*pi*w^2))*exp(-(xi-mu).^2./(2*w^2));
newparms.g_func = @(xi,k1,k2) k1*exp(k2*xi);

% forcible detachment functions
newparms.k_func = @(parms)parms.k*(parms.xi>parms.dLcrit);
newparms.b_func = @(parms)parms.b*(parms.xi>(parms.ps-parms.w)&parms.xi<(parms.ps+parms.w));

% SE functions
newparms.kse_func = @(dlse,parms)parms.kse0*parms.kse*exp(parms.kse*dlse);
newparms.Fse_func = @(dlse,parms)parms.kse0*(exp(parms.kse*dlse)-1);
newparms.Lse_func = @(F,parms)log(F/parms.kse0+1)/parms.kse;

% force-velocity functions
newparms.Fv_func = @(Frel,parms)parms.vmax/(2*parms.e(2))*(-exp(parms.e(4)/parms.e(1)-Frel/parms.e(1))+exp(Frel/parms.e(1)-parms.e(4)/parms.e(1))-2*parms.e(3));
newparms.vF_func = @(vcerel,parms)parms.e(1)*log((parms.e(2)*vcerel./parms.vmax+parms.e(3))+sqrt((parms.e(2)*vcerel./parms.vmax+parms.e(3)).^2+1))+parms.e(4);

end