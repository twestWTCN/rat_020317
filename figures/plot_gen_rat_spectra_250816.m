function plot_gen_rat_spectra_250816(R)
close all
cmap = linspecer(3);
%% Generic plotter for bivariate connectivity measures.
spctralist = R.spectra.featspecs;
undirected = 1;

for feat = 1:length(spctralist) % run through features
    % get ft spectra names
    if strncmp(spctralist{feat},'power',8)
        spctrname = 'spectra'; dirna = 'Power';
        undirected = 1; variate = 1;
    elseif strncmp(spctralist{feat},'coherence',8)
        spctrname = 'cohspctrm'; dirna = 'Coherence';
        undirected = 1; variate = 2;
    elseif strncmp(spctralist{feat},'wpli',8)
        spctrname = 'wpli_debiasedspctrm'; dirna = 'WPLI';
        undirected = 1; variate = 2;
    elseif strncmp(spctralist{feat},'granger',8)
        spctrname = 'grangerspctrm'; dirna = 'Granger';
        undirected = 0; variate = 2;
    elseif strncmp(spctralist{feat},'npd',8)
        spctrname = 'npdspctrm'; dirna = 'NPD';
        undirected = 0; variate = 2;
    elseif strncmp(spctralist{feat},'dtf',8)
        spctrname = 'dtfspctrm'; dirna = 'DTF';
        undirected = 0; variate = 2;
    elseif strncmp(spctralist{feat},'nscoh',8)
        spctrname = 'nscohspctrm'; dirna = 'NS Coherence';
        undirected = 1; variate = 2;
    elseif strncmp(spctralist{feat},'icoherence',8)
        spctrname = 'cohspctrm'; dirna = 'Imaginary Coherence';
        undirected = 1; variate = 2; imagset = 1;
    end
    clear fyA
    for cond = 1:2
        for sub  = 1:length(R.subnames{cond})
            load([R.analysispath R.pipestamp '\data\processed\' R.subnames{cond}{sub} '_' R.condnames{cond} '_' R.pipestamp '.mat'])
            if variate == 1 % univariate (power)
                clear fxA
                for i = 1:length(R.sourcenames)
                    set1 = find(strncmp([R.sourcenames{i}], FTdata.freqPow.label,size(R.sourcenames{i},2)));
                    fxA(:,i) = FTdata.freqPow.freq;
                    fyA{i,cond,sub} = mean(FTdata.spectanaly.spectra{1}(set1,:),1);
                end
            elseif variate == 2 % bivariate (connectivity)
                if undirected == 1
                    
                    if exist('imagset','var')
                        spctralist{feat} = 'coherence';
                    end
                    
                    clear fxA
                    for i = 1:length(R.sourcenames)
                        set1 = find(strncmp([R.sourcenames{i}], FTdata.freqPow.label,size(R.sourcenames{i},2)));
                        for j = 1:length(R.sourcenames)
                            set2 = find(strncmp([R.sourcenames{j}], FTdata.freqPow.label,size(R.sourcenames{j},2)));
                            if exist('imagset','var')
                                fy = squeeze(nanmean(nanmean(imag(FTdata.(spctralist{feat}).(spctrname)(set1,set2,:)),1),2));
                            else
                                fy = squeeze(nanmean(nanmean(abs(FTdata.(spctralist{feat}).(spctrname)(set1,set2,:)),1),2));
                            end
                            if numel(fy)<2 || isempty(fy)
                                fy = NaN(length(FTdata.(spctralist{feat}).freq),1);
                            end
                            fyA{i,j,cond,sub} = fy';
                        end
                    end
                    fxA = FTdata.(spctralist{feat}).freq;
                elseif undirected == 0 % Group for directional measuremets
                    clear fxA
                    for i = 1:length(R.sourcenames)
                        set1 = find(strncmp([R.sourcenames{i}], FTdata.freqPow.label,size(R.sourcenames{i},2)));
                        for j = 1:length(R.sourcenames)
                            set2 = find(strncmp([R.sourcenames{j}], FTdata.freqPow.label,size(R.sourcenames{j},2)));
                            for dirc = 1:3
                                fy = squeeze(mean(mean(FTdata.(spctralist{feat}).(spctrname){dirc}(set1,set2,:),1),2));
                                fyA{i,j,dirc,cond,sub} = fy;
                            end
                        end
                    end
                    fxA{1} = FTdata.(spctralist{feat}).freq;
                    
                end % directed
            end % variate
        end % sub
    end % cond
    clear A
    for large = 0:1
        if large == 1
            FOI = [2 45];
        else
            FOI = R.FOI;
        end
        if undirected == 1 && variate == 1
            for logpl = 1:2
                for i = 1:length(R.sourcenames)
                    figure
                    for ex = 1:length(R.condnames)
                        A(ex,:) = nanmean(vertcat(fyA{i,ex,:}),1);
                        if logpl == 2
                            loglog(fxA(:,1)', A(ex,:)','Color',cmap(ex,:))
                        else
                            plot(fxA(:,1)', A(ex,:)','Color',cmap(ex,:))
                        end
                        hold on
                    end
                    ylim([0 1.2*max(A(:))]);xlim(FOI);% ylim([0 3e-5]);
                    title([R.sourcenames{i} ' Source Power Spectra'], 'Interpreter', 'none')
                    legend(R.condnames)
                    ylabel('Normalised Power'); xlabel('Frequency (Hz)');
                    if ~exist([R.analysispath R.pipestamp '\results\figures\spectra\' dirna], 'dir')
                        mkdir([R.analysispath R.pipestamp '\results\figures\spectra\' dirna]);
                    end
                    if logpl == 2
                        appender = 'logpl';
                    else
                        appender = [];
                    end
                    if large == 1
                        saveallfiguresFIL([R.analysispath R.pipestamp '\results\figures\spectra\' dirna '\' dirna '_' R.sourcenames{i} '_' appender '_large_spectra'],'-jpg'); close all
                    else
                        saveallfiguresFIL([R.analysispath R.pipestamp '\results\figures\spectra\' dirna '\' dirna '_' R.sourcenames{i} '_' appender '_spectra'],'-jpg'); close all
                    end
                end
            end
        elseif undirected == 1 && variate == 2
            for i = 1:length(R.sourcenames)
                for j = 1:length(R.sourcenames)
                    figure
                    for ex = 1:length(R.condnames)
                        A{ex} = nanmean(vertcat(fyA{i,j,ex,:}),1);
                        plot(fxA', A{ex}','Color',cmap(ex,:))
                        hold on
                    end
                    if sum(find(isnan([A{:}])))>1
                        ylim([0 1]);xlim(FOI)
                    elseif exist('imagset','var') && large == 1
                        ylim([-0.2 0.2]); xlim(FOI);
                    elseif exist('imagset','var')
                        ylim([-0.5 0.5]); xlim(FOI);
                    elseif strncmp(spctralist{feat},'coherence',8)
                        ylim([0 1]);xlim(FOI)
                        %                     else
                        %                         ylim([0 1.2*max(A(:))]);xlim(FOI)
                        
                    elseif large == 1
                        ylim([0 0.5]);xlim(FOI)
                        
                    else
                        ylim([0 1]); xlim(FOI);
                    end
                    title([R.sourcenames{i} ' <-> ' R.sourcenames{j} ' Source ' dirna ' Spectra'], 'Interpreter', 'none')
                    legend(R.condnames)
                    ylabel(spctralist{feat}); xlabel('Frequency (Hz)')
                    if ~exist([R.analysispath R.pipestamp '\results\figures\spectra\' dirna], 'dir')
                        mkdir([R.analysispath R.pipestamp '\results\figures\spectra\' dirna]);
                    end
                    if large
                        saveallfiguresFIL([R.analysispath R.pipestamp '\results\figures\spectra\' dirna '\' dirna '_' R.sourcenames{i} '_' R.sourcenames{j} '_large_spectra'],'-jpg'); close all
                    else
                        saveallfiguresFIL([R.analysispath R.pipestamp '\results\figures\spectra\' dirna '\' dirna '_' R.sourcenames{i} '_' R.sourcenames{j} '_spectra'],'-jpg'); close all
                    end
                end
            end
        elseif undirected == 0 && variate == 2 %% DIRECTIONAL PLOTS
            for i = 1:length(R.sourcenames)
                for j = 1:length(R.sourcenames)
                    figure
                    lst = {'-','--'};
                    for ex = 1:length(R.condnames)
                        for dirc = 1:3
                            figure(dirc)
                            A= nanmean([fyA{i,j,dirc,ex,:}],2);
                            plot(fxA{:}, A,'Color',cmap(dirc,:),'LineStyle',lst{ex})
                            hold on
                            if large == 1 && ~strncmp(spctralist{feat},'npd',8)
                                topcoh = 0.5;
                            elseif strncmp(spctralist{feat},'npd',8) && large == 1;
                                topcoh = 0.3;
                            else
                                topcoh = 1;
                                %                         legend({'LSN zero lag','LSN->','LSN <-','CNTRL zero lag','CNTRL ->','CNTRL <-'})
                            end
                            dirsign = {' <-> ',' -> ',' <- '};
                            ylim([0 topcoh]);xlim(FOI)
                            title([R.sourcenames{i} dirsign{dirc} R.sourcenames{j} ' Source ' dirna ' Spectra'], 'Interpreter', 'none')
                            ylabel(spctralist{feat}); xlabel('Frequency (Hz)')
                            legend(R.condnames)
                        end
                    end
                    if ~exist([R.analysispath R.pipestamp '\results\figures\spectra\' dirna], 'dir')
                        mkdir([R.analysispath R.pipestamp '\results\figures\spectra\' dirna]);
                    end
                    % save output
                    if large
                        saveallfiguresFIL([R.analysispath R.pipestamp '\results\figures\spectra\' dirna '\' dirna '_' R.sourcenames{i} '_' R.sourcenames{j} '_large_spectra'],'-jpg'); close all
                    else
                        saveallfiguresFIL([R.analysispath R.pipestamp '\results\figures\spectra\' dirna '\' dirna '_' R.sourcenames{i} '_' R.sourcenames{j} '_spectra'],'-jpg'); close all
                    end
                    
                end % j
            end % i
        end
    end % LARGE
    if exist('imagset','var')
        clear imagset
    end
end
