function plot_gen_rat_spectra_100816(R)
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
    end
    
    for cond = 1:2
        for sub  = 1:length(R.subnames{cond})
            load([R.analysispath R.pipestamp '\data\processed\' R.subnames{cond}{sub} '_' R.condnames{cond} '_' R.pipestamp '.mat'])
            if variate == 1 % univariate (power)
                for i = 1:length(R.sourcenames)
                    set1 = find(strncmp([R.sourcenames{i}], FTdata.freqPow.label,size(R.sourcenames{i},2)));
                    fxA(:,i) = FTdata.freqPow.freq;
                    fyA{i,cond,sub} = mean(FTdata.spectanaly.spectra{1}(set1,:),1);
                end
            elseif variate == 2 % bivariate (connectivity)
                clear fy
                for i = 1:length(R.sourcenames)
                    set1 = find(strncmp([R.sourcenames{i}], FTdata.freqPow.label,size(R.sourcenames{i},2)));
                    for j = 1:length(R.sourcenames)
                        set2 = find(strncmp([R.sourcenames{j}], FTdata.freqPow.label,size(R.sourcenames{j},2)));
                        
                        fy = squeeze(mean(mean(abs(FTdata.(spctralist{feat}).(spctrname)(set1,set2,:)),1),2));
                        fyA{i,j,cond,sub} = fy';
                    end
                end
                fxA = FTdata.(spctralist{feat}).freq;
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
            for i = 1:length(R.sourcenames)
                figure
                for ex = 1:length(R.condnames)
                    A(ex,:) = nanmean(vertcat(fyA{i,ex,:}),1);
                    plot(fxA', A(ex,:)','Color',cmap(ex,:))
                    hold on
                end
                ylim([0 1.2*max(A(:))]);xlim(FOI)
                title([R.sourcenames{i} ' Source Power Spectra'], 'Interpreter', 'none')
                legend(R.condnames)
                ylabel('Normalised Power'); xlabel('Frequency (Hz)');
                if ~exist([R.analysispath R.pipestamp '\results\figures\spectra\' dirna], 'dir')
                    mkdir([R.analysispath R.pipestamp '\results\figures\spectra\' dirna]);
                end
                if large == 1
                    saveallfiguresFIL([R.analysispath R.pipestamp '\results\figures\spectra\' dirna '\' dirna '_' R.sourcenames{i} '_large_spectra'],'-jpg'); close all
                else
                    saveallfiguresFIL([R.analysispath R.pipestamp '\results\figures\spectra\' dirna '\' dirna '_' R.sourcenames{i} '_spectra'],'-jpg'); close all
                end
            end
        elseif undirected == 1 && variate == 2
            for i = 1:length(R.sourcenames)
                for j = 1:length(R.sourcenames)
                    figure
                    for ex = 1:length(R.condnames)
                        A(ex,:) = nanmean(vertcat(fyA{i,j,ex,:}),1);
                        plot(fxA', A(ex,:)','Color',cmap(ex,:))
                        hold on
                    end
                    if sum(find(isnan(A)))>1
                        ylim([0 1]);xlim(FOI)
                    elseif strncmp(spctralist{feat},'coherence',8)
                        ylim([0 1]);xlim(FOI)
                    else
                        ylim([0 1.2*max(A(:))]);xlim(FOI)
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
        else % For more feature types (directed)
        end
    end
end
