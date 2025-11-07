close all

Q0 = 1;
Q1 = .1;
Q2 = .5;
p = Q1/Q0;
q = Q2./Q0 - p.^2; 
c = [Q0; p; 2*q; Q1; Q2];

x = linspace(-20,20,1000);

nfunc = @(x,c) c(1,:) / sqrt(pi*c(3,:)) * exp(-(x-(c(2,:))).^2 / c(3,:));


% as in code
IGa{1} =  @(x,c)1/2*c(1,:).*erf((x-c(2,:))./sqrt(c(3,:)));
IGa{2} =  @(x,c)1/2*c(1,:).*erf((x-c(2,:))./sqrt(c(3,:))).*c(2,:)-1/2.*sqrt(c(3,:)).*c(1,:)/sqrt(pi).*exp(-(x-c(2,:)).^2./c(3,:));
IGa{3} =  @(x,c)1/2*c(1,:).*erf((x-c(2,:))./sqrt(c(3,:))).*(c(2,:).^2+c(3,:)/2)-1/2*sqrt(c(3,:)).*c(1,:)/sqrt(pi).*exp(-(x-c(2,:)).^2./c(3,:)).*(c(2,:)+min(x,1e4));

% as in paper
IGb{1} =  @(x,c) 1/2*c(1,:).*erf((x-c(2,:))./sqrt(c(3,:)));
IGb{2} =  @(x,c) 1/2*c(4,:).*(erf((x-c(2,:))./sqrt(c(3,:))) + sqrt(c(3,:))/(sqrt(pi)*c(2,:)) .* exp(-(x-c(2,:)).^2./c(3,:)));
IGb{3} =  @(x,c) 1/2*c(5,:).*(erf((x-c(2,:))./sqrt(c(3,:))) + sqrt(c(3,:))/(sqrt(pi)*c(2,:)) .* exp(-(x-c(2,:)).^2./c(3,:)) .* (c(2,:)+x));

% translated from wolfram
IGa{1} =  @(x,c)1/2*c(1,:).*erf((x-c(2,:))./sqrt(c(3,:)));
IGa{2} =  @(x,c)1/2*c(1,:).*erf((x-c(2,:))./sqrt(c(3,:))).*c(2,:)-1/2.*sqrt(c(3,:)).*c(1,:)/sqrt(pi).*exp(-(x-c(2,:)).^2./c(3,:));
IGa{3} =  @(x,c)1/2*c(1,:).*erf((x-c(2,:))./sqrt(c(3,:))).*(c(2,:).^2+c(3,:)/2)-1/2*sqrt(c(3,:)).*c(1,:)/sqrt(pi).*exp(-(x-c(2,:)).^2./c(3,:)).*(c(2,:)+min(x,1e4));

figure(1)
for i = 1:3
    subplot(1,3,i)
    plot(x, IGa{i}(x, c), x, IGb{i}(x,c)); hold on

    plot(x, cumtrapz(x, nfunc(x, c).*x.^(i-1)) + IGa{i}(-20,c),'--')
end

% as in paper

% IGb{3} =  @(x,c) 1/2*c(5,:).*erf((x-c(2,:))./sqrt(c(3,:))) .* (c(2,:).^2+c(3,:)/2)-1/2*sqrt(c(3,:)).*c(1,:)/sqrt(pi).*exp(-(x-c(2,:)).^2./c(3,:)).*(c(2,:)+min(x,1e4));