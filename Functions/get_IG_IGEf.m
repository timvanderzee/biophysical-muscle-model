function[IG, IGef] = get_IG_IGEf(approx)

if approx
    % erf approximation
    erfap = @(x) (exp(x)-exp(-x)) ./ (exp(x)+exp(-x));
    IG{1} =  @(x,c)1/2*c(1,:).*erfap((x-c(2,:))./sqrt(c(3,:)));
    IG{2} =  @(x,c)1/2*c(1,:).*erfap((x-c(2,:))./sqrt(c(3,:))) .* c(2,:)-1/2.*sqrt(c(3,:)).*c(1,:)/sqrt(pi).*exp(-(x-c(2,:)).^2./c(3,:));
    IG{3} =  @(x,c)1/2*c(1,:).*erfap((x-c(2,:))./sqrt(c(3,:))) .* (c(2,:).^2+c(3,:)/2)-1/2*sqrt(c(3,:)).*c(1,:)/sqrt(pi).*exp(-(x-c(2,:)).^2./c(3,:)).*(c(2,:)+min(x,1e4));
    
    % exp approximation
    IGef{1} = @(c,k)(c(1,:)*k(1).*exp(10*tanh((c(3,:)*k(2)^2/4-c(2,:)*k(2))/10)));
    IGef{2} = @(c,k)(c(1,:)*k(1).*exp(10*tanh((c(3,:)*k(2)^2/4-c(2,:)*k(2))/10))).*(c(2,:)-c(3,:)*k(2)/2);
    IGef{3} = @(c,k)(c(1,:)*k(1).*exp(10*tanh((c(3,:)*k(2)^2/4-c(2,:)*k(2))/10))).*((c(2,:)-c(3,:)*k(2)/2).^2+c(3,:)/2);

else
    % exponential function
    IGef{1} = @(c,k)(c(1,:)*k(1).*exp(c(3,:)*k(2)^2/4-c(2,:)*k(2)));
    IGef{2} = @(c,k)(c(1,:)*k(1).*exp(c(3,:)*k(2)^2/4-c(2,:)*k(2))).*(c(2,:)-c(3,:)*k(2)/2);
    IGef{3} = @(c,k)(c(1,:)*k(1).*exp(c(3,:)*k(2)^2/4-c(2,:)*k(2))).*((c(2,:)-c(3,:)*k(2)/2).^2+c(3,:)/2);
    
    % erf function
    IG{1} =  @(x,c)1/2*c(1,:).*erf((x-c(2,:))./sqrt(c(3,:)));
    IG{2} =  @(x,c)1/2*c(1,:).*erf((x-c(2,:))./sqrt(c(3,:))).*c(2,:)-1/2.*sqrt(c(3,:)).*c(1,:)/sqrt(pi).*exp(-(x-c(2,:)).^2./c(3,:));
    IG{3} =  @(x,c)1/2*c(1,:).*erf((x-c(2,:))./sqrt(c(3,:))).*(c(2,:).^2+c(3,:)/2)-1/2*sqrt(c(3,:)).*c(1,:)/sqrt(pi).*exp(-(x-c(2,:)).^2./c(3,:)).*(c(2,:)+min(x,1e4));
end

end