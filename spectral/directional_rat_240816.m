function directional_rat_240816(R)
granger = 1; npd = 0; dtf = 1; 
for cond = 1:2;
    for sub  = 1:length(R.subnames{cond})
        load([R.analysispath R.pipestamp '\data\processed\' R.subnames{cond}{sub} '_' R.condnames{cond} '_' R.pipestamp '.mat'])
        chcombs = combvec(1:length(FTdata.label),1:length(FTdata.label));
        cfg = [];
        cfg.method     = 'mtmfft';
        cfg.tapsmofrq  = 2.5;
        cfg.keeptrials = 'no';
        cfg.output     = 'fourier';
        cfg.keeptrials = 'yes';
        cfg.foilim    = [0 125];
        cfg.taper = 'dpss';
        cfg.pad = 2;
        cfg.padtype = 'zero';
        freqFour    = ft_freqanalysis(cfg,FTdata.EpochedData);
        
        %% Compute FT Non-parametric Granger
        if granger == 1
            cfg           = [];
            cfg.method    = 'granger';
            FTdata.granger       = ft_connectivityanalysis(cfg, freqFour);
            
            cfg           = [];
            cfg.method    = 'instantaneous_causality';
            FTdata.instcaus      = ft_connectivityanalysis(cfg,freqFour);
            for i = 1:length(chcombs)
                grangespctrm{1,1}(chcombs(1,i),chcombs(2,i),:) = FTdata.instcaus.instantspctrm(chcombs(1,i),chcombs(2,i),:);
                grangespctrm{1,2}(chcombs(1,i),chcombs(2,i),:) = FTdata.granger.grangerspctrm(chcombs(1,i),chcombs(2,i),:);
                grangespctrm{1,3}(chcombs(1,i),chcombs(2,i),:) = FTdata.granger.grangerspctrm(chcombs(2,i),chcombs(1,i),:);
            end
            FTdata.granger.grangerspctrm = grangespctrm;
        end
        
        %% DH Non-parametric directionality
        if npd == 1
            for i = 1:length(chcombs)
%                 ftspect = freqFour.fourierspctrm(:,[chcombs(1,i) chcombs(2,i)],:);
%                 [f13,t13,cl13] = sp2a2_R2_TW(FTdata.fsample,8,ftspect);
                 x = FTdata.ContData.trial{1}(chcombs(1,i),:);
                 y = FTdata.ContData.trial{1}(chcombs(2,i),:);
                [f13,~,~]=sp2a2_R2_mt(x',y',FTdata.fsample,8,'M1');
                
                npdspctrm{1,1}(chcombs(1,i),chcombs(2,i),:) = f13(:,10);
                npdspctrm{1,2}(chcombs(1,i),chcombs(2,i),:) = f13(:,12);
                npdspctrm{1,3}(chcombs(1,i),chcombs(2,i),:) = f13(:,11);
                nscohspctrm(chcombs(1,i),chcombs(2,i),:) = f13(:,4);
            end
            % neurospec coherence
            FTdata.nscoh.nscohspctrm = nscohspctrm;
            FTdata.nscoh.freq =  f13(:,1);
            FTdata.nscoh.label = FTdata.label;
            % neurospec NPD
            FTdata.npd.npdspctrm = npdspctrm;
            FTdata.npd.freq =  f13(:,1);
            FTdata.npd.dimord = 'chan_chan_freq';
            FTdata.npd.label = FTdata.label;
        end
        if dtf == 1
            cfg           = [];
            cfg.method    = 'dtf';
            FTdata.dtf      = ft_connectivityanalysis(cfg,freqFour);
            for i = 1:length(chcombs)
                dtfspctrm{1,1}(chcombs(1,i),chcombs(2,i),:) = NaN;
                dtfspctrm{1,2}(chcombs(1,i),chcombs(2,i),:) = FTdata.dtf.dtfspctrm(chcombs(1,i),chcombs(2,i),:);
                dtfspctrm{1,3}(chcombs(1,i),chcombs(2,i),:) = FTdata.dtf.dtfspctrm(chcombs(2,i),chcombs(1,i),:);
            end
             FTdata.dtf.dtfspctrm = dtfspctrm;
        end
        save([R.analysispath R.pipestamp '\data\processed\' R.subnames{cond}{sub} '_' R.condnames{cond} '_' R.pipestamp '.mat'],'FTdata')
    end
end
