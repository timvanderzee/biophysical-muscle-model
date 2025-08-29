clear all; close all; clc

load('parms_v2.mat','sparms','pparms');

optparms = {'f', 'k11', 'k22', 'k21', 'JF', 'J1','J2', 'kon', 'koop', 'kse', 'kse0'};
iFs = [1 2 3, 5, 6, 7, 8, 10, 11];

figure(3)

Ys = nan(length(optparms), iFs(end));
Yp = nan(length(optparms), iFs(end));

for iF = iFs
    for i = 1:length(optparms)
        Ys(i,iF) = eval(['sparms(', num2str(iF),').',optparms{i}]);
        Yp(i,iF) = eval(['pparms(', num2str(iF),').',optparms{i}]);
    end
end

figure(1)
subplot(121);
bar(1:length(optparms), mean(Ys,2,'omitnan')); hold on
errorbar(1:length(optparms), mean(Ys,2,'omitnan'),  std(Ys,1,2,'omitnan'),'.')
% set(gca,'YScale','log')
ylim([0 2000])

subplot(122);
bar(1:length(optparms), mean(Yp,2,'omitnan')); hold on
errorbar(1:length(optparms), mean(Yp,2,'omitnan'),  std(Yp,1,2,'omitnan'),'.')
% set(gca,'YScale','log')
ylim([0 2000])
