function plot_gen_rat_spectra_050816(R)
close all
try
    cmap = linspecer(3);
catch
    cmap = hsv(3);
end
%% Generic plotter for bivariate connectivity measures.
spctralist = R.spectra.featspecs;
undirected = 1;
for srcloc = 1:length(R.sourcenames)
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
            spctrname = 'npdspctrum'; dirna = 'NPD';
            undirected = 0; variate = 2;
        end
        for cond = 1:2
            for sub  = 1:length(R.subnames{cond})
                load([R.analysispath R.pipestamp '\data\processed\' R.subnames{cond}{sub} '_' R.condnames{cond} '_' R.pipestamp '.mat'])
                if variate == 1 % univariate (power)
                    inds = find(strncmp([R.sourcenames{srcloc}], FTdata.freqPow.label,size(R.sourcenames{srcloc},2)));
                    fxA = FTdata.freqPow.freq;
                    fyA{cond,sub} = FTdata.spectanaly.spectra{1}(inds,:);
                elseif variate == 2 % bivariate (connectivity)
                    if undirected == 1
                        clear fy
                        inds = find(strncmp([R.sourcenames{srcloc}], FTdata.freqPow.label,size(R.sourcenames{srcloc},2)));
                        for i = 1:length(inds)
                            set1 = inds(i);
                            set2 = find(strncmp(['STN'],FTdata.coherence.label,3));
                            fy(i,:) = mean(abs(squeeze(FTdata.(spctralist{feat}).(spctrname)(repmat(set1,1,size(set2,2)),set2,:))));
                            %                         if strncmp(spctralist{feat},'wpli',8)
                            %                             plot(FTdata.(spctralist{feat}).freq,fy(i,:))
                            %                             hold on
                            %                         end
                        end
                        if isempty(inds)
                            fy = [];
                        end
                        fxA = FTdata.(spctralist{feat}).freq;
                        fyA{cond,sub} = fy;
                        
                        
                    elseif undirected == 0 % Group for directional measuremets
                        clear fy
                        inds = find(strncmp([R.sourcenames{srcloc}], FTdata.freqPow.label,size(R.sourcenames{srcloc},2)));
                        for dirc = 1:3
                            for i = 1:length(inds)
                                set1 = inds(i);
                                set2 = find(strncmp(['STN'],FTdata.coherence.label,3));
                                fy(i,:) = mean(abs(squeeze(FTdata.(spctralist{feat}).(spctrname){1,dirc}(repmat(set1,1,size(set2,2)),set2,:))));
                            end
                            
                            if isempty(inds)
                                fy = [];
                            end
                            fxA = FTdata.(spctralist{feat}).freq;
                            fyA{dirc,cond,sub} = fy;
                        end
                    end
                end
            end
        end
        clear A
        for large = 0:1
            if large == 1
                FOI = [2 45];
            else
                FOI = R.FOI;
            end
            if undirected == 1 && variate == 1
                figure
                for i = 1:length(R.condnames)
                    A(i,:) = mean(vertcat(fyA{i,:}),1);
                    plot(fxA', A(i,:)','Color',cmap(i,:))
                    hold on
                end
                ylim([0 1.2*max(A(:))]);xlim(FOI)
                title([R.sourcenames{srcloc} ' Source Power Spectra'], 'Interpreter', 'none')
                legend(R.condnames)
                ylabel('Normalised Power'); xlabel('Frequency (Hz)');
                if ~exist([R.analysispath R.pipestamp '\results\figures\spectra\' dirna], 'dir')
                    mkdir([R.analysispath R.pipestamp '\results\figures\spectra\' dirna]);
                end
                if large == 1
                    saveallfiguresFIL([R.analysispath R.pipestamp '\results\figures\spectra\' dirna '\' dirna '_' R.sourcenames{srcloc} '_large_spectra'],'-jpg'); close all
                else
                    saveallfiguresFIL([R.analysispath R.pipestamp '\results\figures\spectra\' dirna '\' dirna '_' R.sourcenames{srcloc} '_spectra'],'-jpg'); close all
                end
                
            elseif undirected == 1 && variate == 2 % BIVARIATE UNDIRECTED
                
                figure
                for i = 1:length(R.condnames)
                    A(i,:) = mean(vertcat(fyA{i,:}),1);
                    plot(fxA', A(i,:)','Color',cmap(i,:))
                    hold on
                end
                
                if sum(find(isnan(A)))>1
                    ylim([0 1]);xlim(FOI)
                else
                    ylim([0 1.2*max(A(:))]);xlim(FOI)
                end
                title(['STN <-> ' R.sourcenames{srcloc} ' Source ' dirna ' Spectra'], 'Interpreter', 'none')
                legend(R.condnames)
                ylabel(spctralist{feat}); xlabel('Frequency (Hz)')
                if ~exist([R.analysispath R.pipestamp '\results\figures\spectra\' dirna], 'dir')
                    mkdir([R.analysispath R.pipestamp '\results\figures\spectra\' dirna]);
                end
                if large
                    saveallfiguresFIL([R.analysispath R.pipestamp '\results\figures\spectra\' dirna '\' dirna '_STN_' R.sourcenames{srcloc} '_large_spectra'],'-jpg'); close all
                else
                    saveallfiguresFIL([R.analysispath R.pipestamp '\results\figures\spectra\' dirna '\' dirna '_STN_' R.sourcenames{srcloc} '_spectra'],'-jpg'); close all
                end
            elseif undirected == 0 && variate == 2 %% DIRECTIONAL PLOTS
                figure
                for i = 1:length(R.condnames)
                    lst = {'-','--'};
                    for dirc = 1:3
                        A= mean(vertcat(fyA{dirc,i,:}),1);
                        plot(fxA{:}', A','Color',cmap(dirc,:),'LineStyle',lst{i})
                        hold on
                    end
                end
                if large == 1
                    topcoh = 0.2;
                else
                    topcoh = 1;
                    legend({'OFF zero lag','OFF->','OFF <-','ON zero lag','ON ->','ON <-'})
                end
                ylim([0 topcoh]);xlim(FOI)
                
                title(['STN <-> ' R.sourcenames{srcloc} ' Source ' dirna ' Spectra'], 'Interpreter', 'none')
                
                ylabel(spctralist{feat}); xlabel('Frequency (Hz)')
                if ~exist([R.analysispath R.pipestamp '\results\figures\spectra\' dirna], 'dir')
                    mkdir([R.analysispath R.pipestamp '\results\figures\spectra\' dirna]);
                end
                if large
                    saveallfiguresFIL([R.analysispath R.pipestamp '\results\figures\spectra\' dirna '\' dirna '_STN_' R.sourcenames{srcloc} '_large_spectra'],'-jpg'); close all
                else
                    saveallfiguresFIL([R.analysispath R.pipestamp '\results\figures\spectra\' dirna '\' dirna '_STN_' R.sourcenames{srcloc} '_spectra'],'-jpg'); close all
                end
                
            end
        end
    end
end