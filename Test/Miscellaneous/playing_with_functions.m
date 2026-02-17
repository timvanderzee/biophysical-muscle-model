close all
clc

x = linspace(-2,2,100);
c = [1 0 .1];

G = @(x,c) c(1)*exp(-(x-c(2)).^2/c(3));

plot(x,G(x,c))

IG = @(x,c) 1/2*sqrt(pi*c(3))*c(1) * erf((x-c(2))/sqrt(c(3)));

figure(2)
plot(x, IG(x,c) + IG(100,c)); hold on
plot(x, cumtrapz(x, G(x,c)),'--')

%% 
close all
G = @(x,c) c(1) / (sqrt(pi*c(3))) * exp(-(x-c(2)).^2/c(3));

% plot(x,G(x,c))

IGf = @(x,c) 1/2*c(1) * erf((x-c(2))/sqrt(c(3)));

figure(2)
plot(x, IGf(x,c) + IGf(100,c)); hold on
plot(x, cumtrapz(x, G(x,c)),'--')


%%
xi = linspace(-100,100,10000);
% G = @(x,c) c(1)/sqrt((2*pi*c(3))) * exp(-(x-c(2)).^2/(c(3)));

IG{1} =  @(x,c)1/2*c(1,:).*erf((x-c(2,:))./sqrt(c(3,:)));
IG{2} =  @(x,c)1/2*c(1,:).*erf((x-c(2,:))./sqrt(c(3,:))).*c(2,:)-1/2.*sqrt(c(3,:)).*c(1,:)/sqrt(pi).*exp(-(x-c(2,:)).^2./c(3,:));
IG{3} =  @(x,c)1/2*c(1,:).*erf((x-c(2,:))./sqrt(c(3,:))).*(c(2,:).^2+c(3,:)/2)-1/2*sqrt(c(3,:)).*c(1,:)/sqrt(pi).*exp(-(x-c(2,:)).^2./c(3,:)).*(c(2,:)+min(x,1e4));

Q0 = .5;
p = 0;
q = .3;
c1 = [Q0; p; q];

parms.b = 1000;
parms.ps2 = .1;

i = 3;

phi = -parms.k * (IG{i}(100, c1) -IG{i}(parms.dLcrit, c1))

phi2 = trapz(xi, -xi.^(i-1) .* parms.k .* (xi>parms.dLcrit) .* G(xi, c1))

gamma = parms.b .*  [1 parms.ps2 parms.w^2 + parms.ps2^2]

gamma2 = trapz(xi, xi.^(i-1) .* parms.f_func(xi, parms.b, parms.w, parms.ps2))
% gamma2 = trapz(xi, xi.^(i-1) .* G(xi, [parms.b parms.ps2 parms.w]))


% figure(1)
% % plot(xi, G()