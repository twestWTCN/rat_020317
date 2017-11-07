function directional_rat_270217(R)
granger = 0; dtf = 0; psi = 0; pdc = 0; npdZ = 0; npdY = 0; npdX = 0; npd = 1;  npdW = 0; ncohXY = 0;

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
            for i = 1:length(chcombs)
                %                 ftspect = freqFour.fourierspctrm(:,[chcombs(1,i) chcombs(2,i)],:);
                %                 [f13,t13,cl13] = sp2a2_R2_TW(FTdata.fsample,8,ftspect);
                x = FTdata.ContData.trial{1}(chcombs(1,i),:);
                y = FTdata.ContData.trial{1}(chcombs(2,i),:);
                %                 [f13,~,~]=sp2a2_R2_mt(x',y',FTdata.fsample,7,'M1');
                [f13,~,~]=sp2a2_R2_mt(x',y',FTdata.fsample,7,'M1');
                npdspctrm{1,1}(chcombs(1,i),chcombs(2,i),:) = f13(:,10);
                npdspctrm{1,2}(chcombs(1,i),chcombs(2,i),:) = f13(:,12);
                npdspctrm{1,3}(chcombs(1,i),chcombs(2,i),:) = f13(:,11);
%                 nscohspctrm(chcombs(1,i),chcombs(2,i),:) = f13(:,4);
                nscohspctrm(chcombs(1,i),chcombs(2,i),:) = f13(:,11)+f13(:,12);
                nsicohspctrm(chcombs(1,i),chcombs(2,i),:) = f13(:,13); % Imaginary
                nsPowSpect(chcombs(1,i),:) = f13(:,2);
            end
            % neurospec coherence
            try
                FTdata = rmfield(FTdata,'nscoh');
            catch
            end
            FTdata.nscoh.nscohspctrm = nscohspctrm;
            FTdata.nscoh.freq =  f13(:,1);
            FTdata.nscoh.label = FTdata.label;
            % neurospec NPD
            try
                FTdata = rmfield(FTdata,'npd');
            catch
            end
            FTdata.npd.npdspctrm = npdspctrm;
            FTdata.npd.freq =  f13(:,1);
            FTdata.npd.dimord = 'chan_chan_freq';
            FTdata.npd.label = FTdata.label;
            % neurospec ICOH
            try
                FTdata = rmfield(FTdata,'nsIcoh');
            catch
            end
            FTdata.nsIcoh.nsIcohspctrm = nsicohspctrm;
            FTdata.nsIcoh.freq =  f13(:,1);
            FTdata.nsIcoh.dimord = 'chan_chan_freq';
            FTdata.nsIcoh.label = FTdata.label;
            
            % neurospec Pow
            try
                FTdata = rmfield(FTdata,'nsPow');
            catch
            end
            FTdata.nsPow.Powspctrm = nsPowSpect;
            FTdata.nsPow.freq =  f13(:,1);
            FTdata.nsPow.dimord = 'chan_freq';
            FTdata.nsPow.label = FTdata.label;
            
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
            for i = 1:length(chcombs)
                x = FTdata.ContData.trial{1}(chcombs(1,i),:);
                y = FTdata.ContData.trial{1}(chcombs(2,i),:);
                z =  mean(FTdata.ContData.trial{1}(fix(mean(find(strncmp('STR', FTdata.label,3)))),:),1);
                [f13,~,~]=sp2_R2a_pc1(x,y,z,FTdata.fsample,2^7);
                npdZspctrm{1,1}(chcombs(1,i),chcombs(2,i),:) = f13(:,10);
                npdZspctrm{1,2}(chcombs(1,i),chcombs(2,i),:) = f13(:,12);
                npdZspctrm{1,3}(chcombs(1,i),chcombs(2,i),:) = f13(:,11);
            end
            
            % neurospec NPD
            FTdata.npdZ.npdZspctrm = npdZspctrm;
            FTdata.npdZ.freq =  f13(:,1);
            FTdata.npdZ.dimord = 'chan_chan_freq';
            FTdata.npdZ.label = FTdata.label;
        end
        if npdY == 1
            for i = 1:length(chcombs)
                x = FTdata.ContData.trial{1}(chcombs(1,i),:);
                y = FTdata.ContData.trial{1}(chcombs(2,i),:);
                z =  mean(FTdata.ContData.trial{1}(fix(mean(find(strncmp('STN', FTdata.label,3)))),:),1);
                [f13,~,~]=sp2_R2a_pc1(x,y,z,FTdata.fsample,2^7);
                npdYspctrm{1,1}(chcombs(1,i),chcombs(2,i),:) = f13(:,10);
                npdYspctrm{1,2}(chcombs(1,i),chcombs(2,i),:) = f13(:,12);
                npdYspctrm{1,3}(chcombs(1,i),chcombs(2,i),:) = f13(:,11);
            end
            
            % neurospec NPD
            FTdata.npdY.npdYspctrm = npdYspctrm;
            FTdata.npdY.freq =  f13(:,1);
            FTdata.npdY.dimord = 'chan_chan_freq';
            FTdata.npdY.label = FTdata.label;
        end
        if npdX == 1
            for i = 1:length(chcombs)
                x = FTdata.ContData.trial{1}(chcombs(1,i),:);
                y = FTdata.ContData.trial{1}(chcombs(2,i),:);
                z =  mean(FTdata.ContData.trial{1}(fix(mean(find(strncmp('GPe', FTdata.label,2)))),:),1);
                [f13,~,~]=sp2_R2a_pc1(x,y,z,FTdata.fsample,2^7);
                npdXspctrm{1,1}(chcombs(1,i),chcombs(2,i),:) = f13(:,10);
                npdXspctrm{1,2}(chcombs(1,i),chcombs(2,i),:) = f13(:,12);
                npdXspctrm{1,3}(chcombs(1,i),chcombs(2,i),:) = f13(:,11);
            end
            
            % neurospec NPD
            FTdata.npdX.npdXspctrm = npdXspctrm;
            FTdata.npdX.freq =  f13(:,1);
            FTdata.npdX.dimord = 'chan_chan_freq';
            FTdata.npdX.label = FTdata.label;
        end
        if npdW == 1
            for i = 1:length(chcombs)
                x = FTdata.ContData.trial{1}(chcombs(1,i),:);
                y = FTdata.ContData.trial{1}(chcombs(2,i),:);
                z =  mean(FTdata.ContData.trial{1}(fix(mean(find(strncmp('M1', FTdata.label,2)))),:),1);
                [f13,~,~]=sp2_R2a_pc1(x,y,z,FTdata.fsample,2^7);
                npdWspctrm{1,1}(chcombs(1,i),chcombs(2,i),:) = f13(:,10);
                npdWspctrm{1,2}(chcombs(1,i),chcombs(2,i),:) = f13(:,12);
                npdWspctrm{1,3}(chcombs(1,i),chcombs(2,i),:) = f13(:,11);
            end
            
            % neurospec NPD
            FTdata.npdW.npdWspctrm = npdWspctrm;
            FTdata.npdW.freq =  f13(:,1);
            FTdata.npdW.dimord = 'chan_chan_freq';
            FTdata.npdW.label = FTdata.label;
        end
        
        if ncohXY == 1
            parfor i = 1:length(chcombs)
                if strcmp('M1',FTdata.ContData.label{chcombs(1,i)}) ||...
                        strncmp('STN',FTdata.ContData.label{chcombs(1,i)},3) || ...
                        strcmp('M1',FTdata.ContData.label{chcombs(2,i)}) ||...
                        strncmp('STN',FTdata.ContData.label{chcombs(2,i)},3)
                    npdspctrmZ2(i,:)= zeros(64,1);
                else
                    if chcombs(1,i) > 19
                        a = 2;
                    end
                    x = FTdata.ContData.trial{1}(chcombs(1,i),:);
                    y = FTdata.ContData.trial{1}(chcombs(2,i),:);
                    z1 =  FTdata.ContData.trial{1}(find(strncmp('M1', FTdata.label,3)),:);
                    z2 = FTdata.ContData.trial{1}(find(strncmp('STN', FTdata.label,3)),:);
                    z = [z1; z2];
                    [f1z2,t1z2] = HOpartcohtw_160517(x',y',z',FTdata.fsample,7,0);
                    npdspctrmZ2(i,:)= f1z2(:,4);
                    frqz(:,i) =  f1z2(:,1);
                end
                
                disp([sub cond i])
            end
            for i = 1:length(chcombs)
                npdspctrmZ2unf(chcombs(1,i),chcombs(2,i),:) = npdspctrmZ2(i,:);
            end
            % neurospec NPD
            FTdata.ncohXY.ncohXYspctrm = npdspctrmZ2unf;
            FTdata.ncohXY.freq =  squeeze(frqz(:,min(find(frqz(1,:)))));
            FTdata.ncohXY.dimord = 'chan_chan_freq';
            FTdata.ncohXY.label = FTdata.label;
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
        disp([sub cond])
    end
end
