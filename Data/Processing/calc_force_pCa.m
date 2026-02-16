function[] = calc_force_pCa(datafolder, outputfolder)

fibers = {'12Dec2017a','13Dec2017a','13Dec2017b','14Dec2017a','14Dec2017b','18Dec2017a','18Dec2017b','19Dec2017a','6Aug2018a','6Aug2018b','7Aug2018a'};

pCas = [4.5 6.1 6.2 6.3 6.4 6.6 9];

F0 = nan(length(fibers), length(pCas));

for k = 1:length(fibers)
    % initialize
    Force   = nan(2e4, 30);
    Length  = nan(2e4, 30);
    Time    = nan(2e4, 30);
    
    
    cd(datafolder)
    filename = [fibers{k},'_activation_data.mat'];
    disp(filename)
    
    if exist(filename, 'file')
        load(filename)
        new_data = act_data;
        
        Fi = nan(1,length(new_data));
        pCa = nan(1, length(new_data));
        
        for j = 1:length(new_data)
            file_info = new_data(j).file_info_string;
            fgain = str2double(file_info(20:21));
            
            % overrule gain
            if k == 3 && j == 2
                fgain = 10;
            elseif k == 3 && j == 3
                fgain = 10;
            elseif k == 5 && j == 10
                fgain = 1;
            elseif k == 9 && (j == 8 || j == 9)
                fgain = 10;
            elseif k == 11 && j == 10
                fgain = 1;
            end
            
            % something went wrong
            if k == 4 && j == 8
                continue
            end
            
            N = length(new_data(j).force);
            Force(1:N,j) = new_data(j).force/1000 / fgain;
            Length(1:N,j) = new_data(j).fl;
            Time(1:N,j) = new_data(j).time;
            
            if k < 9
                tend(k) = .9;
            else
                tend(k) = .1;
            end
            
            id = new_data(j).time < tend(k); % first 100 ms
            Fi(j) = mean(Force(id,j));
            
            if ~isempty(new_data(j).pCa)
                pCa(j) = new_data(j).pCa;
            else
                pCa(j) = nan;
            end
            
            % overrule pCa
            if k == 5 && j == 12
                pCa(j) = 4.5;
            end
            
        end
    end
    
    %     Fmax(k) = max(Fi);
    
    figure(1)
    color = get(gca,'colororder');
    nexttile
    
    for j = 1:length(pCas)
        if sum(pCa == pCas(j)) > 0
            plot(Time(:, pCa == pCas(j)), Force(:, pCa == pCas(j)),'color',color(j,:)); hold on
            
            F0(k,j) = median(Fi(pCa == pCas(j)));
        end
        xline(tend(k),'k--')
        yline(max(F0(k,:)),'k--')
    end
    
    title(fibers{k})
    box off
    
    nexttile
    for j = 1:length(pCas)
        if sum(pCa == pCas(j)) > 0
            plot(pCa(pCa == pCas(j)), Fi(pCa == pCas(j)),'o','color',color(j,:),'markerfacecolor',color(j,:)); hold on
        end
        yline(max(F0(k,:)),'k--')
    end
    box off
    
end

%%
Fmax = max(F0, [], 2);
%
figure(2)
plot(pCas, F0./Fmax,'o')

cd(outputfolder)
save('force_pCa.mat','F0','pCas', 'Fmax')


end





