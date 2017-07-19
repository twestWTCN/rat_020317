function directional_rat_revcomp_050916(R)
revcomp = 1;
for cond = 1:2;
    for sub  = 1:length(R.subnames{cond})
        load([R.analysispath R.pipestamp '\data\processed\' R.subnames{cond}{sub} '_' R.condnames{cond} '_' R.pipestamp '.mat'])
        for rev = 1:revcomp+1
            chcombs = combvec(1:length(FTdata.label),1:length(FTdata.label));
            if revcomp == 1
                EpochedData = FTdata.EpochedData;
                for tr = 1:length(FTdata.EpochedData.trial)
                    EpochedData.trial{tr} = FTdata.EpochedData.trial{tr}(:,size( FTdata.EpochedData.trial{tr},2):-1:1);
                end
                cfg = [];
                cfg.method     = 'mtmfft';
                cfg.tapsmofrq  = 1.5;
                cfg.keeptrials = 'no';
                cfg.output     = 'fourier';
                cfg.keeptrials = 'yes';
                cfg.foilim    = [0 125];
                cfg.taper = 'dpss';
                cfg.pad = 'nextpow2';
                cfg.padtype = 'zero';
                freqFour    = ft_freqanalysis(cfg,EpochedData);
            else
                freqFour = FTdata.freqFour;
            end
            %% Compute FT Non-parametric Granger
            if sum(strcmpi('granger',R.dirmeth))>0
                cfg           = [];
                cfg.method    = 'granger';
                granger       = ft_connectivityanalysis(cfg, freqFour);
                
                cfg           = [];
                cfg.method    = 'instantaneous_causality';
                instcaus      = ft_connectivityanalysis(cfg,freqFour);
                for i = 1:length(chcombs)
                    grangespctrm{1,1}(chcombs(1,i),chcombs(2,i),:) = instcaus.instantspctrm(chcombs(1,i),chcombs(2,i),:);
                    grangespctrm{1,2}(chcombs(1,i),chcombs(2,i),:) = granger.grangerspctrm(chcombs(1,i),chcombs(2,i),:);
                    grangespctrm{1,3}(chcombs(1,i),chcombs(2,i),:) = granger.grangerspctrm(chcombs(2,i),chcombs(1,i),:);
                end
                FTdata.granger(rev).grangerspctrm = grangespctrm;
                FTdata.granger(rev).freq = granger.freq;
            end
            %% DH Non-parametric directionality
            if sum(strcmpi('npd',R.dirmeth))>0
                for i = 1:length(chcombs)
                    %                 ftspect = freqFour.fourierspctrm(:,[chcombs(1,i) chcombs(2,i)],:);
                    %                 [f13,t13,cl13] = sp2a2_R2_TW(FTdata.fsample,8,ftspect);
                    x = FTdata.ContData.trial{1}(chcombs(1,i),:);
                    y = FTdata.ContData.trial{1}(chcombs(2,i),:);
                    [f13,~,~]=sp2a2_R2_mt(x',y',FTdata.fsample,8,'M1');
                    %                 ij = [chcombs(1,i) chcombs(2,i)];
                    %                 [f13,t,cl,sc] = npd_tw_010916(freqFour,ij,size(FTdata.EpochedData.trial{1},2));
                    
                    npdspctrm{1,1}(chcombs(1,i),chcombs(2,i),:) = f13(:,10);
                    npdspctrm{1,2}(chcombs(1,i),chcombs(2,i),:) = f13(:,12);
                    npdspctrm{1,3}(chcombs(1,i),chcombs(2,i),:) = f13(:,11);
                    nscohspctrm(chcombs(1,i),chcombs(2,i),:) = f13(:,4);
                end
                % neurospec coherence
                FTdata.nscoh(rev).nscohspctrm = nscohspctrm;
                FTdata.nscoh(rev).freq =  f13(:,1);
                FTdata.nscoh(rev).label = FTdata.label;
                % neurospec NPD
                FTdata.npd(rev).npdspctrm = npdspctrm;
                FTdata.npd(rev).freq =  f13(:,1);
                FTdata.npd(rev).dimord = 'chan_chan_freq';
                FTdata.npd(rev).label = FTdata.label;
            end
            if sum(strcmpi('dtf',R.dirmeth))>0
                cfg           = [];
                cfg.method    = 'dtf';
                FTdata.dtf      = ft_connectivityanalysis(cfg,freqFour);
                for i = 1:length(chcombs)
                    dtfspctrm{1,1}(chcombs(1,i),chcombs(2,i),:) = NaN;
                    dtfspctrm{1,2}(chcombs(1,i),chcombs(2,i),:) = FTdata.dtf.dtfspctrm(chcombs(1,i),chcombs(2,i),:);
                    dtfspctrm{1,3}(chcombs(1,i),chcombs(2,i),:) = FTdata.dtf.dtfspctrm(chcombs(2,i),chcombs(1,i),:);
                end
                FTdata.dtf(rev).dtfspctrm = dtfspctrm;
            end
        end
        save([R.analysispath R.pipestamp '\data\processed\' R.subnames{cond}{sub} '_' R.condnames{cond} '_' R.pipestamp '.mat'],'FTdata')
    end
end
