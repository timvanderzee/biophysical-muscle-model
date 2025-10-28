clear all; close all; clc
savefig = 0;

[username, githubfolder] = get_paths();


filename = 'biophysical_full_regular_SRS';
versions = {'parms_v1d', 'parms_v4'};


for kk = 1:length(versions)
    
    cd(['C:\Users\u0167448\Documents\GitHub\biophysical-muscle-model\Model output\SRS\', versions{kk}])
    load(filename,'Stest', 'Scond', 'AMPs', 'pCas', 'ISIs', 'F0')
    
    
    data(kk).SRSrel = Stest ./ Scond(:,1,7,:);

    
end

%% plot 
% semilogx(ISIs, data(1).


%%

delta = data(1).SRSrel - data(2).SRSrel;

%%
aAMPs = repmat(AMPs, length(ISIs),1);
aISIs = repmat(ISIs(:), 1, length(AMPs));

close all
iFs = [2,3,5,6,7,8,11];

for iF = 5
% 
%     figure(1)
    for j = 1:2
        figure(j)
        
    for i = 1:size(delta,1)
        nexttile
%         surf(aAMPs, aISIs, squeeze(delta(i,:,:,iF)))
        
        surf(aAMPs, aISIs, squeeze(data(j).SRSrel(i,:,:,iF)))
        
        set(gca, 'yscale', 'log')
        zlim([0 2])
    end
    

    end
end

%% compare forces
% delta(3, 
