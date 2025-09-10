clear all; close all; clc

ISIs = repmat(1:7,8,1);
AMPs = repmat((1:8)',1,7);

figure(1)
color = get(gca,'colororder');

f = 10;

for i = 1:7
    
    scolor = brighten(color(1,:) + [.3 0 0], .1*i);
    scolor2 = brighten([1 .2 .2], .1*i);
    
    Z = f * i * ones(size(AMPs));
    Z(1, 3:7,:) = nan;
    
    Z2 = nan(size(Z));
    Z2(6:7,3:4) = f * i;
    Z2(1:2,1:2) = f * i;
    
    surf(AMPs,ISIs, Z, 'FaceColor', scolor, 'Edgecolor', [.2 .2 .2]); hold on
    surf(AMPs,ISIs, Z2, 'FaceColor', scolor2, 'Edgecolor', [.2 .2 .2]); hold on
end

grid off
ylabel('Recovery time')
xlabel('Amplitude')
zlabel('Calcium activation')
axis([0 8 0 9 0 i*f+f])

set(gcf,'units','normalized')
set(gcf,'position', [ 0.1300    0.1100    0.4   0.8150]);

%%
figure(1)
view(70,8)