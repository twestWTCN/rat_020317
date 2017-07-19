function directional_rat_030517(R)
granger = 0; dtf = 0; psi = 0; pdc = 0; npdZ = 0; npdY = 0; npdX = 0; npd = 1;  npdW = 0;
revcheck = 1; groupstat = 0;
if groupstat == 1
    npd_xstore = [];
    npd_ystore = [];
end

for cond = 1:2;
    for sub  = 1:length(R.subnames{cond})
        load([R.analysispath R.pipestamp '\data\processed\' R.subnames{cond}{sub} '_' R.condnames{cond} '_' R.pipestamp '.mat'])
        chcombs = combvec(1:length(FTdata.label),1:length(FTdata.label));
        
        % MVAR COEFFICIENTS
        %         cfg         = [];
        %         cfg.order   = 25;
        %         cfg.toolbox = 'bsmart';
        %         mdata       = ft_mvaranalysis(cfg, FTdata.EpochedData);
        %         cfg = [];
        %         cfg.method     = 'mvar';
        %         freqFour    = ft_freqanalysis(cfg,mdata);
        
        %         cfg = [];
        %         cfg.method     = 'mtmfft';
        %         cfg.keeptrials = 'yes';
        %         cfg.output     = 'fourier';
        %         cfg.foilim    = [0 125];
        %         cfg.taper = 'dpss';
        
        %         cfg.tapsmofrq = 3;
        %         cfg.pad = 2;
        %         cfg.padtype = 'zero';
        %         freqFour    = ft_freqanalysis(cfg,FTdata.EpochedData);
        
        
        
        %% Compute FT Non-parametric Granger
        if granger == 1
            cfg           = [];
            cfg.method    = 'granger';
            FTdata.granger       = ft_connectivityanalysis(cfg, FTdata.freqFour);
            
            cfg           = [];
            cfg.method    = 'instantaneous_causality';
            FTdata.instcaus      = ft_connectivityanalysis(cfg,FTdata.freqFour);
            for i = 1:length(chcombs)
                grangespctrm{1,1}(chcombs(1,i),chcombs(2,i),:) = FTdata.instcaus.instantspctrm(chcombs(1,i),chcombs(2,i),:);
                grangespctrm{1,2}(chcombs(1,i),chcombs(2,i),:) = FTdata.granger.grangerspctrm(chcombs(1,i),chcombs(2,i),:);
                grangespctrm{1,3}(chcombs(1,i),chcombs(2,i),:) = FTdata.granger.grangerspctrm(chcombs(2,i),chcombs(1,i),:);
            end
            FTdata.granger.grangerspctrm = grangespctrm;
        end
        
        %% DH Non-parametric directionality
        if npd == 1
            for rev = 1:revcheck+1
                for i = 1:length(chcombs)
                    %                 ftspect = freqFour.fourierspctrm(:,[chcombs(1,i) chcombs(2,i)],:);
                    %                 [f13,t13,cl13] = sp2a2_R2_TW(FTdata.fsample,8,ftspect);
                    x = FTdata.ContData.trial{1}(chcombs(1,i),:);
                    y = FTdata.ContData.trial{1}(chcombs(2,i),:);
                    
                    if groupstat == 1
                        npd_xstore = [npd_xstore x];
                        npd_ystore = [npd_ystore x];
                    end
                    
                    if rev>1
                        x = fliplr(x);
                        y = fliplr(y);
                    end
                    %                 [f13,~,~]=sp2a2_R2_mt(x',y',FTdata.fsample,7,'M1');
                    [f13,~,~]=sp2a2_R2(x',y',FTdata.fsample,7);
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
        end
        
        if dtf == 1
            cfg           = [];
            cfg.method    = 'dtf';
            FTdata.dtf      = ft_connectivityanalysis(cfg,FTdata.freqFour);
            for i = 1:length(chcombs)
                dtfspctrm{1,1}(chcombs(1,i),chcombs(2,i),:) = NaN;
                dtfspctrm{1,2}(chcombs(1,i),chcombs(2,i),:) = FTdata.dtf.dtfspctrm(chcombs(1,i),chcombs(2,i),:);
                dtfspctrm{1,3}(chcombs(1,i),chcombs(2,i),:) = FTdata.dtf.dtfspctrm(chcombs(2,i),chcombs(1,i),:);
            end
            FTdata.dtf.dtfspctrm = dtfspctrm;
            clear dtfspctrm
        end
        
        if psi == 1
            cfg           = [];
            cfg.method    = 'psi';
            cfg.bandwidth = 2.5;
            FTdata.psi      = ft_connectivityanalysis(cfg,FTdata.freqFour);
            for i = 1:length(chcombs)
                psispctrm{1,1}(chcombs(1,i),chcombs(2,i),:) = NaN;
                psispctrm{1,2}(chcombs(1,i),chcombs(2,i),:) = FTdata.psi.psispctrm(chcombs(1,i),chcombs(2,i),:);
                psispctrm{1,3}(chcombs(1,i),chcombs(2,i),:) = FTdata.psi.psispctrm(chcombs(2,i),chcombs(1,i),:);
            end
            FTdata.psi.psispctrm = psispctrm;
            clear psispctrm
        end
        
        if npdZ == 1
            for rev = 1:revcheck+1
                for i = 1:length(chcombs)
                    x = FTdata.ContData.trial{1}(chcombs(1,i),:);
                    y = FTdata.ContData.trial{1}(chcombs(2,i),:);
                    if rev>1
                        x = fliplr(x);
                        y = fliplr(y);
                    end
                    z =  mean(FTdata.ContData.trial{1}(fix(mean(find(strncmp('STR', FTdata.label,3)))),:),1);
                    [f13,~,~]=sp2_R2a_pc1(x,y,z,FTdata.fsample,2^7);
                    npdZspctrm{1,1}(chcombs(1,i),chcombs(2,i),:) = f13(:,10);
                    npdZspctrm{1,2}(chcombs(1,i),chcombs(2,i),:) = f13(:,12);
                    npdZspctrm{1,3}(chcombs(1,i),chcombs(2,i),:) = f13(:,11);
                end
                
                % neurospec NPD
                FTdata.npdZ(rev).npdZspctrm = npdZspctrm;
                FTdata.npdZ(rev).freq =  f13(:,1);
                FTdata.npdZ(rev).dimord = 'chan_chan_freq';
                FTdata.npdZ(rev).label = FTdata.label;
            end
        end
        if npdY == 1
            for rev = 1:revcheck+1
                for i = 1:length(chcombs)
                    x = FTdata.ContData.trial{1}(chcombs(1,i),:);
                    y = FTdata.ContData.trial{1}(chcombs(2,i),:);
                    if rev>1
                        x = fliplr(x);
                        y = fliplr(y);
                    end
                    z =  mean(FTdata.ContData.trial{1}(fix(mean(find(strncmp('STN', FTdata.label,3)))),:),1);
                    [f13,~,~]=sp2_R2a_pc1(x,y,z,FTdata.fsample,2^7);
                    npdYspctrm{1,1}(chcombs(1,i),chcombs(2,i),:) = f13(:,10);
                    npdYspctrm{1,2}(chcombs(1,i),chcombs(2,i),:) = f13(:,12);
                    npdYspctrm{1,3}(chcombs(1,i),chcombs(2,i),:) = f13(:,11);
                end
                
                % neurospec NPD
                FTdata.npdY(rev).npdYspctrm = npdYspctrm;
                FTdata.npdY(rev).freq =  f13(:,1);
                FTdata.npdY(rev).dimord = 'chan_chan_freq';
                FTdata.npdY(rev).label = FTdata.label;
            end
        end
        if npdX == 1
            for rev = 1:revcheck+1
                for i = 1:length(chcombs)
                    x = FTdata.ContData.trial{1}(chcombs(1,i),:);
                    y = FTdata.ContData.trial{1}(chcombs(2,i),:);
                    if rev>1
                        x = fliplr(x);
                        y = fliplr(y);
                    end
                    z =  mean(FTdata.ContData.trial{1}(fix(mean(find(strncmp('GP', FTdata.label,2)))),:),1);
                    [f13,~,~]=sp2_R2a_pc1(x,y,z,FTdata.fsample,2^7);
                    npdXspctrm{1,1}(chcombs(1,i),chcombs(2,i),:) = f13(:,10);
                    npdXspctrm{1,2}(chcombs(1,i),chcombs(2,i),:) = f13(:,12);
                    npdXspctrm{1,3}(chcombs(1,i),chcombs(2,i),:) = f13(:,11);
                end
                
                % neurospec NPD
                FTdata.npdX(rev).npdXspctrm = npdXspctrm;
                FTdata.npdX(rev).freq =  f13(:,1);
                FTdata.npdX(rev).dimord = 'chan_chan_freq';
                FTdata.npdX(rev).label = FTdata.label;
            end
        end
        if npdW == 1
            for rev = 1:revcheck+1
                for i = 1:length(chcombs)
                    x = FTdata.ContData.trial{1}(chcombs(1,i),:);
                    if rev>1
                        x = fliplr(x);
                        y = fliplr(y);
                    end
                    y = FTdata.ContData.trial{1}(chcombs(2,i),:);
                    z =  mean(FTdata.ContData.trial{1}(fix(mean(find(strncmp('fEEG', FTdata.label,2)))),:),1);
                    [f13,~,~]=sp2_R2a_pc1(x,y,z,FTdata.fsample,2^7);
                    npdWspctrm{1,1}(chcombs(1,i),chcombs(2,i),:) = f13(:,10);
                    npdWspctrm{1,2}(chcombs(1,i),chcombs(2,i),:) = f13(:,12);
                    npdWspctrm{1,3}(chcombs(1,i),chcombs(2,i),:) = f13(:,11);
                end
                
                % neurospec NPD
                FTdata.npdW(rev).npdWspctrm = npdWspctrm;
                FTdata.npdW(rev).freq =  f13(:,1);
                FTdata.npdW(rev).dimord = 'chan_chan_freq';
                FTdata.npdW(rev).label = FTdata.label;
            end
        end
        
        if pdc == 1
            cfg           = [];
            cfg.method    = 'pdc';
            FTdata.pdc      = ft_connectivityanalysis(cfg,FTdata.freqFour);
            for i = 1:length(chcombs)
                pdcspctrm{1,1}(chcombs(1,i),chcombs(2,i),:) = NaN;
                pdcspctrm{1,2}(chcombs(1,i),chcombs(2,i),:) = FTdata.pdc.pdcspctrm(chcombs(1,i),chcombs(2,i),:);
                pdcspctrm{1,3}(chcombs(1,i),chcombs(2,i),:) = FTdata.pdc.pdcspctrm(chcombs(2,i),chcombs(1,i),:);
            end
            FTdata.pdc.pdcspctrm = pdcspctrm;
            FTdata.pdc = rmfield(FTdata.pdcbase,'pdcspctrm');
            clear pdcspctrm
        end
        %         cfg = [];
        %         cfg.parameter = 'pdcspctrm';
        %         ft_connectivityplot(cfg,  FTdata.pdc)
        save([R.analysispath R.pipestamp '\data\processed\' R.subnames{cond}{sub} '_' R.condnames{cond} '_' R.pipestamp '.mat'],'FTdata')
    end
    if groupstat == 1
        [f,t,cl] =sp2a2_R2(npd_xstore',npd_ystore',FTdata.fsample,7);
    end
end
