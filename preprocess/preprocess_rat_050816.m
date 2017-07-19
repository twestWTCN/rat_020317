function preprocess_rat_050816(R)
if R.clear.pp == 1
    delete([R.analysispath R.pipestamp '\data\processed\*.mat'])
end

for cond = 1:2
    for sub  = 1:length(R.subnames{cond})
        load([R.analysispath R.pipestamp '\data\raw\' R.subnames{cond}{sub} '_' R.condnames{cond} '_' R.pipestamp '.mat'])
        
        % Preprocess
        cfg = [];
        cfg.hpfilter  = 'yes'; % highpass
        cfg.hpfreq = R.pp.hpfilt(2);
        cfg.lpfilter  = 'no';  % low pass
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
        
        %% Muscle Artefact LFP
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
        
        
        
        %% CONTINUOUS
        cfg = [];
        cfg.start = 2;
        cfg.end = 2;    % TRUNCATE
        ContData = preprocess_cont_rat_050816(cfg,FTdata);
        load([R.analysispath R.pipestamp '\data\raw\' R.subnames{cond}{sub} '_' R.condnames{cond} '_' R.pipestamp '.mat'])
        
        FTdata.EpochedData = FTdataX;
        FTdata.ContData = ContData;
        if ~exist([R.analysispath R.pipestamp '\data\processed\'], 'dir')
            mkdir([R.analysispath R.pipestamp '\data\processed\']);
        end
        save([R.analysispath R.pipestamp '\data\processed\' R.subnames{cond}{sub} '_' R.condnames{cond} '_' R.pipestamp '.mat'],'FTdata')
    end
end
