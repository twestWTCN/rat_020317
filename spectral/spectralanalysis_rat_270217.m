function spectralanalysis_rat_270217(R)
for cond = 1:2
    for sub  = 1:length(R.subnames{cond})
        load([R.analysispath R.pipestamp '\data\processed\' R.subnames{cond}{sub} '_' R.condnames{cond} '_' R.pipestamp '.mat'])
        % Power
        cfg            = [];
        cfg.output     = 'powandcsd';
        cfg.method     = 'mtmfft';
        cfg.foilim     = R.FOI;
        cfg.taper = 'hanning';
%         cfg.tapsmofrq  = 3;
%         cfg.pad = 'nextpow2';
        cfg.keeptrials = 'no';
        freqPow    = ft_freqanalysis(cfg,FTdata.EpochedData);
        FTdata.freqPow = freqPow;
        
        % Fourier
        cfg = [];
        cfg.method     = 'mtmfft';
         cfg.tapsmofrq  = 2;
%         cfg.keeptrials = 'no';
        cfg.output     = 'fourier';
        cfg.keeptrials = 'yes';
%         cfg.taper = 'hanning';
         cfg.taper = 'dpss';
        freqFour    = ft_freqanalysis(cfg,FTdata.EpochedData);
        FTdata.freqFour = freqFour;
        
        % store to variables
        ft = freqPow.freq; fxx = freqPow.powspctrm;
        res= ft(2)-ft(1); % Frequency resolution
        
        % Compute confidence intervals
        CI = repmat(0.8512*sqrt(1/length(FTdata.EpochedData.trial)),size(fxx));    % Analytic (DH)
        
        %Correct powerline and normalise by full power
        for i = 1:size(fxx,1)
            fxxi = fxx(i,:);
            fxxi(round((48-R.FOI(1))/res+1):round((52-R.FOI(1))/res+1)) = mean(fxxi(round((46-R.FOI(1))/res+1):round((48-R.FOI(1))/res+1)));    % Correct powerline to mean power in prior 6 Hz
%             fxxi = fxxi/(res*trapz(fxxi(round((4-R.FOI(1))/res)+1:round((48-R.FOI(1))/res)+1)));  % Normalise by power between 4 and 48
            fxx(i,:) = fxxi;
            
            CIi = CI(i,:);  % Form confidence interval
            CIi = CIi/(res*trapz(CIi)); % normalise
            CI(i,:) = CIi;
        end

        FTdata.spectanaly.spectra = {fxx};
        FTdata.spectanaly.CI = CI;
        
        if R.spectanaly.cplot(1) == 1
            figure
            subplot(2,2,1)
            inds = find(strncmp(['fEEG'], freqPow.label,3)); % Source index
            plot(repmat(freqPow.freq,length(inds),1)',FTdata.spectanaly.spectra{1}(inds,:)');
            title('fEEG')
            subplot(2,2,2)
            inds = find(strncmp(['GP'], freqPow.label,2)); % Source index
            plot(repmat(freqPow.freq,length(inds),1)',FTdata.spectanaly.spectra{1}(inds,:)');
            title('GP')
            subplot(2,2,3)
            inds = find(strncmp(['STR'], freqPow.label,3)); % Source index
            plot(repmat(freqPow.freq,length(inds),1)',FTdata.spectanaly.spectra{1}(inds,:)');
            title('STR')
            subplot(2,2,4)
            inds = find(strncmp(['STN'], freqPow.label,3)); % Source index
            plot(repmat(freqPow.freq,length(inds),1)',FTdata.spectanaly.spectra{1}(inds,:)');
            title('STN')
        end
        
        % Compute Coherence
        cfg = [];
        cfg.method = 'coh';
        cfg.foilim     = R.FOI;
        cfg.complex = 'complex';
        coherence = ft_connectivityanalysis(cfg, freqFour);  % Compute coherence
        if R.spectanaly.cplot(2) == 1
            cfg = [];
            cfg.parameter = 'cohspctrm';
            figure;ft_connectivityplot(cfg, coherence);
        end
        FTdata.coherence = coherence;
        
        % Compute WPLI
        cfg = [];
        cfg.method = 'wpli_debiased';
        cfg.foilim     = R.FOI;
        wpli = ft_connectivityanalysis(cfg, freqFour);
        if R.spectanaly.cplot(3) == 1
            cfg = [];
            cfg.parameter = 'wpli_debiasedspctrm';
            figure;ft_connectivityplot(cfg, wpli);
        end
        FTdata.wpli = wpli;
        
        save([R.analysispath R.pipestamp '\data\processed\' R.subnames{cond}{sub} '_' R.condnames{cond} '_' R.pipestamp '.mat'],'FTdata')
    end
end