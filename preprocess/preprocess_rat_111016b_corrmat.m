function preprocess_rat_111016b_corrmat(R)
if R.clear.pp == 1
    delete([R.analysispath R.pipestamp '\data\processed\*.mat'])
end

for cond = 1:2
    for sub  = 1:length(R.subnames{cond})
        for L = 1:size(R.montname,2)
            load([R.analysispath R.pipestamp '\data\raw\' R.subnames{cond}{sub} '_' R.condnames{cond} '_' R.pipestamp '.mat'])
            
            R.montage = R.montname{L};
            %% Montage
            FTdata = montage_020317(FTdata,R);
            %         save([R.analysispath R.pipestamp '\data\processed\' R.subnames{cond}{sub} '_' R.condnames{cond} '_' R.pipestamp '.mat'],'FTdata')
            
            %% Preprocess Epoched
            cfg = [];
            cfg.hpfilter  = 'yes'; % highpass
            cfg.hpfreq = R.pp.hpfilt(2);
            cfg.lpfilter  = R.pp.lpfilt(2);  % low pass
            cfg.demean = 'yes';
            cfg.dftfilter = 'yes'; % line noise filter
            FTdata = ft_preprocessing(cfg,FTdata);
            
            % Epoched Data
            % First epoch into segments (1s)
            cfg = [];
            cfg.length  = 2;
            FTdataX = ft_redefinetrial(cfg, FTdata);
            
            
            trialstokeep = 1:numel(FTdataX.trial);
            trialstokeep([1,2 numel(FTdataX.trial)-1 numel(FTdataX.trial)]) = [];
            cfg = [];
            cfg.trials = trialstokeep;
            FTdataX = ft_preprocessing(cfg, FTdataX);
            
            % Muscle Artefact LFP
            cfg = [];
            % channel selection, cutoff and padding
            cfg.artfctdef.zvalue.channel = '*';
            cfg.artfctdef.zvalue.cutoff      = 4;
            cfg.artfctdef.zvalue.trlpadding  = 0;
            cfg.artfctdef.zvalue.fltpadding  = 0;
            cfg.artfctdef.zvalue.artpadding  = 0.1;
            
            % algorithmic parameters
            cfg.artfctdef.zvalue.bpfilter    = 'yes';
            cfg.artfctdef.zvalue.bpfreq      = [75 95];
            %         cfg.artfctdef.zvalue.bpfiltord   = 25;
            cfg.artfctdef.zvalue.bpfilttype  = 'but';
            cfg.artfctdef.zvalue.hilbert     = 'yes';
            cfg.artfctdef.zvalue.boxcar      = 0.2;
            
            [cfg, artifact_muscle] = ft_artifact_zvalue(cfg,FTdataX);
            
            
            
            %% Preprocess Continuous Data
            cfg = [];
            cfg.start = 2;
            cfg.end = 2;    % TRUNCATE
            ContData = preprocess_cont_rat_050816(cfg,FTdata);
            %         load([R.analysispath R.pipestamp '\data\processed\' R.subnames{cond}{sub} '_' R.condnames{cond} '_' R.pipestamp '.mat'])
            
            FTdata.EpochedData = FTdataX;
            FTdata.ContData = ContData;
            covarplot =1;
            if covarplot == 1
                figure
                imagesc(corr(ContData.trial{1}'))
                set(gca,'YDir','normal')
                set(gca,'XTick',[1:size(ContData.trial{1},1)]); set(gca,'XTickLabel',ContData.label); set(gca,'FontSize',7)
                set(gca,'YTick',[1:size(ContData.trial{1},1)]); set(gca,'YTickLabel',ContData.label); set(gca,'FontSize',7)
                title({['Correlation Matrix: ' R.montname{L}]; [R.condnames{cond} ' rat ' R.subnames{cond}{sub}]})
                set(gcf,'Position',[   593 49 1173 929])
                colorbar
                caxis([-1 1])
            end
            clear FTdata ContData EpochedData FTdataX
        end
        saveallfiguresSFLAP([R.analysispath R.pipestamp '\results\montaging\montagetest_' R.subnames{cond}{sub} '_' R.pipestamp '_' ],'-jpg')
        close all
        %         save([R.analysispath R.pipestamp '\data\processed\' R.subnames{cond}{sub} '_' R.condnames{cond} '_' R.pipestamp '.mat'],'FTdata')
    end
end
