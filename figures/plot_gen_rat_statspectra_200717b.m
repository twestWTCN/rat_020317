function plot_gen_rat_statspectra_200717b(R)
close all
cmap = linspecer(3);
%% Generic plotter for bivariate connectivity measures.
spctralist = R.spectra.featspecs;
undirected = 1; part = 0;

for feat = 1:length(spctralist) % run through features
    % get ft spectra names
    if strncmp(spctralist{feat},'nsIcoh',8)
        spctrname = 'nsIcohspctrm'; dirna = 'iCOH'; ylabd = 'iCOH';
        undirected = 1; variate = 2;
    elseif strncmp(spctralist{feat},'nsPow',8)
        spctrname = 'Powspctrm'; dirna = 'Power'; ylabd = 'log Power';
        undirected = 1; variate = 1;
    elseif strncmp(spctralist{feat},'power',8)
        spctrname = 'spectra'; dirna = 'Powera';
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
        spctrname = 'npdspctrm'; dirna = 'NPD'; savefl = 'True';
        undirected = 0; variate = 2; ylabd = 'NPD';
    elseif strncmp(spctralist{feat},'npdZ',8)
        spctrname = 'npdZspctrm'; dirna = 'NPDzstr';
        undirected = 0; variate = 2; part = 1; ylabd = 'NPD'; parSource = 'STR';
    elseif strncmp(spctralist{feat},'npdY',8)
        spctrname = 'npdYspctrm'; dirna = 'NPDxstn';
        undirected = 0; variate = 2; part = 1; ylabd = 'NPD'; parSource = 'STN';
    elseif strncmp(spctralist{feat},'npdX',8)
        spctrname = 'npdXspctrm'; dirna = 'NPDxgp';
        undirected = 0; variate = 2;part = 1; ylabd = 'NPD'; parSource = 'GPe';
    elseif strncmp(spctralist{feat},'npdW',8)
        spctrname = 'npdWspctrm'; dirna = 'NPDwm1';
        undirected = 0; variate = 2; part = 1; ylabd = 'NPD'; parSource = 'M1';
    elseif strncmp(spctralist{feat},'dtf',8)
        spctrname = 'dtfspctrm'; dirna = 'DTF';
        undirected = 0; variate = 2;
    elseif strncmp(spctralist{feat},'nscoh',8)
        spctrname = 'nscohspctrm'; dirna = 'NS Coherence';
        undirected = 1; variate = 2;
    elseif strncmp(spctralist{feat},'ncohXY',8)
        spctrname = 'ncohXYspctrm'; dirna = 'hoCOH on STN M1';
        undirected = 1; variate = 2;
    elseif strncmp(spctralist{feat},'icoherence',8)
        spctrname = 'cohspctrm'; dirna = 'Imaginary Coherence';
        undirected = 1; variate = 2; imagset = 1;
    elseif strncmp(spctralist{feat},'psi',8)
        spctrname = 'psispctrm'; dirna = 'Phase Slope Index';
        undirected = 0; variate = 2;
    elseif strncmp(spctralist{feat},'pdc',8)
        spctrname = 'pdcspctrm'; dirna = 'Partial Directed Coherence';
        undirected = 0; variate = 2;
    end
    if ~exist('savefl')
        savefl = 'False';
    end
    
    if ~exist('parSource')
        parSource = [];
    end
    clear fyA
    for cond = 1:2
        for sub  = 1:length(R.subnames{cond})
            load([R.analysispath R.pipestamp '\data\processed\' R.subnames{cond}{sub} '_' R.condnames{cond} '_' R.pipestamp '.mat'])
            if variate == 1 % univariate (power)
                clear fxA
                for i = 1:length(R.sourcenames)
                    set1 = find(strncmp([R.sourcenames{i}], FTdata.nsPow.label,size(R.sourcenames{i},2)));
                    fxA(:,i) = FTdata.nsPow.freq;
                    fyA{i,cond,sub} = mean(FTdata.nsPow.Powspctrm(set1,:),1);
                end
            elseif variate == 2 % bivariate (connectivity)
                if undirected == 1
                    
                    if exist('imagset','var')
                        spctralist{feat} = 'coherence';
                    end
                    
                    clear fxA
                    for i = 1:length(R.sourcenames)
                        set1 = find(strncmp([R.sourcenames{i}], FTdata.nsPow.label,size(R.sourcenames{i},2)));
                        for j = 1:length(R.sourcenames)
                            set2 = find(strncmp([R.sourcenames{j}], FTdata.nsPow.label,size(R.sourcenames{j},2)));
                            if exist('imagset','var')
                                fy = squeeze(nanmean(nanmean(abs(imag(FTdata.(spctralist{feat}).(spctrname)(set1,set2,:))),1),2));
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
                        set1 = find(strncmp([R.sourcenames{i}], FTdata.nsPow.label,size(R.sourcenames{i},2)));
                        for j = 1:length(R.sourcenames)
                            set2 = find(strncmp([R.sourcenames{j}], FTdata.nsPow.label,size(R.sourcenames{j},2)));
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
    if istrue(savefl)
        fxA_base = fxA; fyA_base = fyA;
        if ~exist([R.analysispath R.pipestamp '\statistics\save_base\'], 'dir')
            mkdir([R.analysispath R.pipestamp '\statistics\save_base\']);
        end
        save([R.analysispath R.pipestamp '\statistics\save_base\' dirna '_basesave'],'fxA_base','fyA_base');
    end
    if size(fxA,2)>size(fxA,1)
        fxA = fxA';
    end
    clear A
    for large = 1% 0:1;
        if large == 1
            appender = '_large';
            FOI = [2 45];
        else
            appender = [];
            FOI = R.FOI;
        end
        if undirected == 1 && variate == 1
            for i = 2:length(R.sourcenames)
                figure
                statspec =1;
                if statspec == 1
                    clf;
                    ON = vertcat(fyA{i,1,:});
                    OFF = vertcat(fyA{i,2,:});
                    ON = 10.^ON;
                    OFF = 10.^OFF;
                    alpha = 0.025;
                    [specstat] = spectralclusterstats060317(ON,OFF,fxA(:,1)',R.sourcenames{i},5000);
                    %                     writespecstat(R.sourcenames{i},specstat,R)
                    [ax clustat] = freqclusterplotlog_230517(OFF,ON,fxA(:,1)',specstat,alpha,[],[],[5e-5 2e-2],cmap);
                    if ~exist([R.analysispath R.pipestamp '\results\figures\spectra\' dirna '\clusterstat'], 'dir')
                        mkdir([R.analysispath R.pipestamp '\results\figures\spectra\' dirna '\clusterstat']);
                    end
                    csvwrite([R.analysispath R.pipestamp '\results\figures\spectra\' dirna '\clusterstat\' dirna '_' R.sourcenames{i} appender '_clusterstat.csv'],clustat)
                end
                
                % ylim([0 3e-5]);
               
                title([R.sourcenames{i} ' Power Spectra'], 'Interpreter', 'none','fontsize',20,'fontweight','bold')
                legend(ax,R.condnames)
                ylabel('log10 Normalised Power','fontsize',20); xlabel('Frequency (Hz)','fontsize',20);
                if ~exist([R.analysispath R.pipestamp '\results\figures\spectra\' dirna], 'dir')
                    mkdir([R.analysispath R.pipestamp '\results\figures\spectra\' dirna]);
                end
                savefigure_v2([R.analysispath R.pipestamp '\results\figures\spectra\' dirna '\'],[dirna '_' R.sourcenames{i} appender '_spectra'],[],'-tiff','-r300'); close all
%                 saveallfiguresFIL([R.analysispath R.pipestamp '\results\figures\spectra\' dirna '\' dirna '_' R.sourcenames{i} appender '_spectra'],'-jpg',0); close all
            end
        elseif undirected == 1 && variate == 2
            for i = 1:length(R.sourcenames)
                for j = 1:length(R.sourcenames)
                    figure
                    statspec =1;
                    
                    if exist('imagset','var') && large == 1
                        ylimz = [0 0.4]; xlim(FOI);
                    elseif exist('imagset','var')
                        ylimz = [0 0.6]; xlim(FOI);
                    elseif strncmp(spctralist{feat},'coherence',8)
                        ylimz =[0 1];xlim(FOI)
                        %                     else
                        %                         ylim([0 1.2*max(A(:))]);xlim(FOI)
                        
                    elseif large == 1
                        ylimz = [0 0.25]; xlim(FOI)
                        
                    else
                        ylimz = [0 0.5]; xlim(FOI);
                    end
                    
                    
                    if statspec == 1
                        clf;
                        ON =vertcat(fyA{i,j,1,:});
                        OFF = vertcat(fyA{i,j,2,:});
                        alpha = 0.05;
                        [specstat] = spectralclusterstats060317(ON,OFF,fxA(:,1)',R.sourcenames{i},5000);
                        [ax clustat] = freqclusterplot_200717(OFF,ON,fxA(:,1),specstat,alpha,[],[],ylimz,cmap);
                        if ~exist([R.analysispath R.pipestamp '\results\figures\spectra\' dirna '\clusterstat'], 'dir')
                            mkdir([R.analysispath R.pipestamp '\results\figures\spectra\' dirna '\clusterstat']);
                        end
                        csvwrite([R.analysispath R.pipestamp '\results\figures\spectra\' dirna '\clusterstat\' dirna '_' R.sourcenames{i} '_' R.sourcenames{j} appender '_clusterstat.csv'],clustat)
                        
                    end
                    %                     ylim(ylimz)
                    %                     title([R.sourcenames{i} ' <-> ' R.sourcenames{j} ' Source ' dirna ' Spectra'], 'Interpreter', 'none','fontsize',16,'fontweight','bold')
                    title([R.sourcenames{i} ' <-> ' R.sourcenames{j}], 'Interpreter', 'none','fontsize',16,'fontweight','bold')
                    legend(ax,R.condnames,'fontsize',12)
                    ylabel(ylabd,'fontsize',16); xlabel('Frequency (Hz)','fontsize',16)
                    if ~exist([R.analysispath R.pipestamp '\results\figures\spectra\' dirna], 'dir')
                        mkdir([R.analysispath R.pipestamp '\results\figures\spectra\' dirna]);
                    end
                    
                    saveallfiguresFIL([R.analysispath R.pipestamp '\results\figures\spectra\' dirna '\' dirna '_' R.sourcenames{i} '_' R.sourcenames{j} appender '_spectra'],'-jpg',0); close all
                end
            end
        elseif undirected == 0 && variate == 2 %% DIRECTIONAL PLOTS
            for i = 1:length(R.sourcenames)
                for j = 2:length(R.sourcenames)
                    if ~strcmp(parSource,R.sourcenames{j}) & ~strcmp(parSource,R.sourcenames{i})
                        figure
                        dirsign = {' <-> ',' -> ',' <- '};
                        
                        lst = {'-','--'};
                        for dirc = 2:3
                            figure(dirc)
                            if strncmp(spctralist{feat},'npd',3) && large == 1;
                                ylimz = [0 0.2];
                            elseif strncmp(spctralist{feat},'psi',8) && large == 0;
                                ylimz = [-.2 .2];
                            elseif strncmp(spctralist{feat},'psi',8) && large == 1;
                                ylimz = [-.5 .5];
                            elseif large == 1 && ~strncmp(spctralist{feat},'npd',3)
                                ylimz = [0 0.5];
                            else
                                ylimz = [0 .5];
                            end
                            if  dirc == 1
                                ylimz =[0 1];
                            end
                            statspec = 1;
                            if (strncmp(spctralist{feat},'dtf',8)||strncmp(spctralist{feat},'psi',8)||strncmp(spctralist{feat},'pdc',8))  && (dirc==1)
                                for ex = 1:2
                                    A= nanmean([fyA{i,j,dirc,ex,:}],2);
                                    plot(fxA{:}, A,'Color',cmap(ex,:),'LineStyle',lst{ex});
                                    hold on
                                end
                                ylim(ylimz);xlim(FOI)
                                legend(R.condnames)
                            elseif statspec == 1
                                clf;
                                ON =[fyA{i,j,dirc,1,:}];
                                OFF = [fyA{i,j,dirc,2,:}];
                                if dirc == 1 && i==j
                                    ON =[fyA{i,j,2,1,:}]+[fyA{i,j,3,1,:}];
                                    OFF = [fyA{i,j,2,2,:}]+[fyA{i,j,3,2,:}];
                                end
                                alpha = 0.05;
                                [specstat] = spectralclusterstats060317(ON',OFF',fxA{1}',R.sourcenames{i},5000);
                                [ax clustat] = freqclusterplot_200717(OFF',ON',fxA{1},specstat,alpha,[],ylabd,ylimz,cmap);
                                legend(ax,R.condnames)
                                %                             title([R.sourcenames{i} dirsign{dirc} R.sourcenames{j} ' Source ' dirna ' Spectra'], 'Interpreter', 'none','fontsize',20,'fontweight','bold')
                                title([R.sourcenames{i} dirsign{dirc} R.sourcenames{j}], 'Interpreter', 'none','fontsize',20,'fontweight','bold')
%                                 ylabel(ylabd,'fontsize',16); xlabel('Frequency (Hz)','fontsize',16)
                                
                                if ~exist([R.analysispath R.pipestamp '\results\figures\spectra\' dirna '\clusterstat'], 'dir')
                                    mkdir([R.analysispath R.pipestamp '\results\figures\spectra\' dirna '\clusterstat']);
                                end
                                csvwrite([R.analysispath R.pipestamp '\results\figures\spectra\' dirna '\clusterstat\' dirna '_' R.sourcenames{i} '_' R.sourcenames{j} '_' appender '_dirc' num2str(dirc) '_clusterstat.csv'],clustat)
                                
                                if part==1
                                    ON_part = ON; OFF_part = OFF;
                                    load([R.analysispath R.pipestamp '\statistics\save_base\' dirna(1:3) '_basesave'])
                                    ON =[fyA_base{i,j,dirc,1,:}];
                                    OFF = [fyA_base{i,j,dirc,2,:}];
                                    if dirc == 1 && i==j
                                        ON =[fyA_base{i,j,2,1,:}]+[fyA_base{i,j,3,1,:}];
                                        OFF = [fyA_base{i,j,2,2,:}]+[fyA_base{i,j,3,2,:}];
                                    end
                                    crmap = linspecer(5);
                                    [specstat] = spectralclusterstats060317(OFF',OFF_part',fxA{1}',R.sourcenames{i},5000);
                                    [clustat] = freqcluster_only_200717(OFF',ON',fxA{1},specstat,alpha,ylimz,crmap(1,:),0.007);
                                    csvwrite([R.analysispath R.pipestamp '\results\figures\spectra\' dirna '\clusterstat\' dirna '_' R.sourcenames{i} '_' R.sourcenames{j} '_' appender '_dirc' num2str(dirc) '_cluster_stat_baseline_ON.csv'],clustat)
                                    [specstat] = spectralclusterstats060317(ON',ON_part',fxA{1}',R.sourcenames{i},5000);
                                    [clustat] = freqcluster_only_200717(OFF',ON',fxA{1},specstat,alpha,ylimz,crmap(2,:),0.022);
                                    csvwrite([R.analysispath R.pipestamp '\results\figures\spectra\' dirna '\clusterstat\' dirna '_' R.sourcenames{i} '_' R.sourcenames{j} '_' appender '_dirc' num2str(dirc) '_cluster_stat_baseline_OFF.csv'],clustat)
                                end
                                
                            end
                            
                            
                            %                         legend({'LSN zero lag','LSN->','LSN <-','CNTRL zero lag','CNTRL ->','CNTRL <-'})
                            
                            clear ax
                        end
                        
                        if ~exist([R.analysispath R.pipestamp '\results\figures\spectra\' dirna], 'dir')
                            mkdir([R.analysispath R.pipestamp '\results\figures\spectra\' dirna]);
                        end
                        % save output
                        savefigure_v2([R.analysispath R.pipestamp '\results\figures\spectra\' dirna '\'],[dirna '_' R.sourcenames{i} '_' R.sourcenames{j} appender '_spectra'],1:3,'-tiff','-r300'); close all

%                         saveallfiguresFIL([R.analysispath R.pipestamp '\results\figures\spectra\' dirna '\' dirna '_' R.sourcenames{i} '_' R.sourcenames{j} appender '_spectra'],'-jpg',0); close all
                        close all
                    end
                end
            end % j        end % i
        end
    end % LARGE
    if exist('imagset','var')
        clear imagset
    end
end
