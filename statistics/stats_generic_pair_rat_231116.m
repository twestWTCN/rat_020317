function stats_generic_pair_rat_231116(R)
statlist = R.genpairstatlist;
% Generic compiler for stat structures
for srcloc = 2 %1:length(R.sourcenames)
    for cond = 2 %1:2
        for sub  = 1:length(R.subnames{cond})
            load([R.analysispath R.pipestamp '\data\processed\' R.subnames{cond}{sub} '_' R.condnames{cond} '_' R.pipestamp '.mat'])
            
            inds = find(strncmpi([R.sourcenames{srcloc}], FTdata.coherence.label,length(R.sourcenames{srcloc}))); % Source index
            bbounds = R.bbounds;
            %% Coherence
            if sum(strcmpi(statlist, 'COH'))>0
                if R.clear.genstat == 1 && srcloc == 1
                    FTdata.dirstats.coh = [];
                end
                clear fy
                intCoh = []; peakMag = []; peakFrq = [];
                if isempty(strmatch(R.sourcenames{srcloc},'fEEG'))
                    for i = 1:length(inds)
                        set1 = inds(i);
                        set2 = find(strncmpi(['fEEG'],FTdata.coherence.label,3)); % STN index
                        fy =squeeze(FTdata.coherence.cohspctrm(set1,set2,:)); % Combined Coherence spectra
                        if size(fy,1)>1 && size(fy,2)>1
                            fy = mean(fy,1);
                        end
                        fxA = FTdata.coherence.freq; % Frequency X vector
                        res = fxA(2) - fxA(1);  % Freq resolution
                        for band = 1:size(R.bbounds,1)
                            intCoh(band,i) = abs(res*trapz(fy(round((bbounds(band,1)-R.FOI(1))/res)+1:round((bbounds(band,2)-R.FOI(1))/res)+1)));
                            [peakMag(band,i) frq] = max(fy(round((bbounds(band,1)-R.FOI(1))/res)+1:round((bbounds(band,2)-R.FOI(1))/res)+1));
                            peakFrq(band,i) =  ((frq-1)*res)+bbounds(band,1);
                        end
                    end
                    FTdata.dirstats.coh.(R.sourcenames{srcloc}).banintCoh = intCoh;    % Integrated Coherence
                    FTdata.dirstats.coh.(R.sourcenames{srcloc}).peakMag = peakMag;  % peak magnitude within band
                    FTdata.dirstats.coh.(R.sourcenames{srcloc}).peakFrq = peakFrq;  % peak frequency (Hz)
                else
                    FTdata.dirstats.coh.(R.sourcenames{srcloc}).banintCoh = [];    % Integrated Coherence
                    FTdata.dirstats.coh.(R.sourcenames{srcloc}).peakMag = [];  % peak magnitude within band
                    FTdata.dirstats.coh.(R.sourcenames{srcloc}).peakFrq = [];  % peak frequency (Hz)
                end
                clear intCoh peakMag peakFrq
            end
            %% WPLI
            if sum(strcmpi(statlist, 'WPLI'))>0
                if R.clear.genstat == 1 && srcloc == 1
                    FTdata.dirstats.wpli = [];
                end
                intCoh = []; peakMag = []; peakFrq = [];
                if isempty(strmatch(R.sourcenames{srcloc},'fEEG'))
                    for i = 1:length(inds)
                        set1 = inds(i);
                        set2 = find(strncmpi(['fEEG'],FTdata.wpli.label,3)); % STN index
                        fy =squeeze(FTdata.wpli.wpli_debiasedspctrm(set1,set2,:)); % Combined Coherence spectra
                        if size(fy,1)>1 && size(fy,2)>1
                            fy = nanmean(fy,1);
                        end
                        fxA = FTdata.wpli.freq; % Frequency X vector
                        res = fxA(2) - fxA(1);  % Freq resolution
                        for band = 1:size(R.bbounds,1)
                            intCoh(band,i) = abs(res*trapz(fy(round((bbounds(band,1)-R.FOI(1))/res)+1:round((bbounds(band,2)-R.FOI(1))/res)+1)));
                            [peakMag(band,i) frq] = max(fy(round((bbounds(band,1)-R.FOI(1))/res)+1:round((bbounds(band,2)-R.FOI(1))/res)+1));
                            peakFrq(band,i) =  ((frq-1)*res)+bbounds(band,1);
                        end
                    end
                    FTdata.dirstats.wpli.(R.sourcenames{srcloc}).intCoh = intCoh;    % Integrated Coherence
                    FTdata.dirstats.wpli.(R.sourcenames{srcloc}).peakMag = peakMag;  % peak magnitude within band
                    FTdata.dirstats.wpli.(R.sourcenames{srcloc}).peakFrq = peakFrq;  % peak frequency (Hz)
                else
                    FTdata.dirstats.wpli.(R.sourcenames{srcloc}).intCoh = [];    % Integrated Coherence
                    FTdata.dirstats.wpli.(R.sourcenames{srcloc}).peakMag = [];   % peak magnitude within band
                    FTdata.dirstats.wpli.(R.sourcenames{srcloc}).peakFrq = [];   % peak frequency (Hz)
                end
                clear intCoh peakMag peakFrq
                
            end
            %% Mutual information
            if sum(strcmp(statlist, 'MI'))>0
                if R.clear.genstat == 1 && srcloc == 1
                    FTdata.dirstats.MI = [];
                    FTdata.dirstats.MIph = [];
                end
                
                clear fy
                for bands = 1:size(R.bbounds,1)
                    for i = 1:length(inds)
                        set1 = inds(i);
                        set2 = find(strncmpi(['fEEG'],FTdata.label,3));
                        fy(bands,i) = FTdata.MI.mi_norm_rec{bands}(set1,set2);
                        fy2(bands,i) = FTdata.MIphase.mi_norm_rec{bands}(set1,set2);
                    end
                end
                if isempty(inds)
                    fy = []; fy2 = [];
                end
                FTdata.dirstats.MI.(R.sourcenames{srcloc}).minorm = fy;
                FTdata.dirstats.MIph.(R.sourcenames{srcloc}).minorm = fy2;
            end
            %% DFA for Phase Sync
            if sum(strcmpi(statlist, 'DFAPS'))>0
                if R.clear.genstat == 1 && srcloc == 1
                    FTdata.dirstats.dfaps = [];
                end
                
                DFAPS = FTdata.DFAPS;
                if isfield(DFAPS,R.sourcenames{srcloc})
                    if ~isfield(DFAPS.(R.sourcenames{srcloc}),'DFAStore')
                        alpha = [];
                    else
                        alpha =  DFAPS.(R.sourcenames{srcloc}).DFAStore;   % alpha exponents
                    end
                    if ~isfield(DFAPS.(R.sourcenames{srcloc}),'peb')
                        peb    =  []; % peb scores
                    else
                        peb =  squeeze(DFAPS.(R.sourcenames{srcloc}).peb); % peb scores
                    end
                    
                    if ~isfield(DFAPS.(R.sourcenames{srcloc}),'ARCoeff')
                        ARCoeff    =  []; % peb scores
                    else
                        ARCoeff =  squeeze(DFAPS.(R.sourcenames{srcloc}).ARCoeff); % peb scores
                    end
                    if ~isfield(DFAPS.(R.sourcenames{srcloc}),'varfeat')
                        varfeat    =  []; % peb scores
                    else
                        varfeat =  squeeze(DFAPS.(R.sourcenames{srcloc}).varfeat); % peb scores
                    end
                    
                    alpha = alpha.*(peb>R.dfaps.BF_r);
                    rejrate = sum((peb>-2),1)./size(alpha,1)';
                    alpha(alpha==0) = NaN;
                    
                    pebstore{sub} = peb;
                    
                    FTdata.dirstats.dfaps.(R.sourcenames{srcloc}).alpha = alpha';
                    FTdata.dirstats.dfaps.(R.sourcenames{srcloc}).pebscore = peb';
                    FTdata.dirstats.dfaps.(R.sourcenames{srcloc}).rejrate = rejrate;
                    FTdata.dirstats.dfaps.(R.sourcenames{srcloc}).varfeat = varfeat';
                    FTdata.dirstats.dfaps.(R.sourcenames{srcloc}).ARCoeff = ARCoeff';
                else
                    FTdata.dirstats.dfaps.(R.sourcenames{srcloc}).alpha = [];
                    FTdata.dirstats.dfaps.(R.sourcenames{srcloc}).pebscore = [];
                    FTdata.dirstats.dfaps.(R.sourcenames{srcloc}).rejrate = [];
                    FTdata.dirstats.dfaps.(R.sourcenames{srcloc}).varfeat = [];
                    FTdata.dirstats.dfaps.(R.sourcenames{srcloc}).ARCoeff = [];
                end
            end
            if sum(strcmpi(statlist, 'DFAAE'))>0
                if R.clear.genstat == 1 && srcloc == 1
                    FTdata.dirstats.dfaae = [];
                end
                
                DFAAE = FTdata.DFAAE;
                if isfield(DFAAE,R.sourcenames{srcloc})
                    if ~isfield(DFAAE.(R.sourcenames{srcloc}),'DFAStore')
                        alpha = [];
                    else
                        alpha =  DFAAE.(R.sourcenames{srcloc}).DFAStore;   % alpha exponents
                    end
                    if ~isfield(DFAAE.(R.sourcenames{srcloc}),'peb')
                        peb    =  []; % peb scores
                    else
                        peb =  squeeze(DFAAE.(R.sourcenames{srcloc}).peb); % peb scores
                    end
                    if ~isfield(DFAAE.(R.sourcenames{srcloc}),'ARCoeff')
                        ARCoeff    =  []; % peb scores
                    else
                        ARCoeff =  squeeze(DFAAE.(R.sourcenames{srcloc}).ARCoeff); % peb scores
                    end
                    if ~isfield(DFAAE.(R.sourcenames{srcloc}),'varfeat')
                        varfeat    =  []; % peb scores
                    else
                        varfeat =  squeeze(DFAAE.(R.sourcenames{srcloc}).varfeat); % peb scores
                    end
                    alpha = alpha.*(peb>R.dfaae.BF_r);
                    rejrate = sum((peb>R.dfaps.BF_r),1)./size(alpha,1)';
                    alpha(alpha==0) = NaN;
                    
                    FTdata.dirstats.dfaae.(R.sourcenames{srcloc}).alpha = alpha';
                    FTdata.dirstats.dfaae.(R.sourcenames{srcloc}).pebscore = peb';
                    FTdata.dirstats.dfaae.(R.sourcenames{srcloc}).rejrate = rejrate';
                    FTdata.dirstats.dfaae.(R.sourcenames{srcloc}).varfeat = varfeat';
                    FTdata.dirstats.dfaae.(R.sourcenames{srcloc}).ARCoeff = ARCoeff';
                else
                    FTdata.dirstats.dfaae.(R.sourcenames{srcloc}).alpha = [];
                    FTdata.dirstats.dfaae.(R.sourcenames{srcloc}).pebscore = [];
                    FTdata.dirstats.dfaae.(R.sourcenames{srcloc}).rejrate = [];
                    FTdata.dirstats.dfaae.(R.sourcenames{srcloc}).varfeat = [];
                    FTdata.dirstats.dfaae.(R.sourcenames{srcloc}).ARCoeff = [];
                end
            end
            if sum(strcmpi(statlist, 'DFAIF'))>0
                if R.clear.genstat == 1 && srcloc == 1
                    FTdata.dirstats.dfaae = [];
                end
                
                DFAIF = FTdata.DFAIF;
                if isfield(DFAIF,R.sourcenames{srcloc})
                    if ~isfield(DFAIF.(R.sourcenames{srcloc}),'DFAStore')
                        alpha = [];
                    else
                        alpha =  DFAIF.(R.sourcenames{srcloc}).DFAStore;   % alpha exponents
                    end
                    if ~isfield(DFAIF.(R.sourcenames{srcloc}),'peb')
                        peb    =  []; % peb scores
                    else
                        peb =  squeeze(DFAIF.(R.sourcenames{srcloc}).peb); % peb scores
                    end
                    if ~isfield(DFAIF.(R.sourcenames{srcloc}),'ARCoeff')
                        ARCoeff    =  []; % peb scores
                    else
                        ARCoeff =  squeeze(DFAIF.(R.sourcenames{srcloc}).ARCoeff); % peb scores
                    end
                    if ~isfield(DFAIF.(R.sourcenames{srcloc}),'varfeat')
                        varfeat    =  []; % peb scores
                    else
                        varfeat =  squeeze(DFAIF.(R.sourcenames{srcloc}).varfeat); % peb scores
                    end
                    alpha = alpha.*(peb>R.dfaif.BF_r);
                    rejrate = sum((peb>R.dfaif.BF_r),1)./size(alpha,1)';
                    alpha(alpha==0) = NaN;
                    
                    FTdata.dirstats.dfaif.(R.sourcenames{srcloc}).alpha = alpha';
                    FTdata.dirstats.dfaif.(R.sourcenames{srcloc}).pebscore = peb';
                    FTdata.dirstats.dfaif.(R.sourcenames{srcloc}).rejrate = rejrate';
                    FTdata.dirstats.dfaif.(R.sourcenames{srcloc}).varfeat = varfeat';
                    FTdata.dirstats.dfaif.(R.sourcenames{srcloc}).ARCoeff = ARCoeff';
                else
                    FTdata.dirstats.dfaif.(R.sourcenames{srcloc}).alpha = [];
                    FTdata.dirstats.dfaif.(R.sourcenames{srcloc}).pebscore = [];
                    FTdata.dirstats.dfaif.(R.sourcenames{srcloc}).rejrate = [];
                    FTdata.dirstats.dfaif.(R.sourcenames{srcloc}).varfeat = [];
                    FTdata.dirstats.dfaif.(R.sourcenames{srcloc}).ARCoeff = [];
                end
            end
            
            if sum(strcmpi(statlist, 'npd'))>0
                if R.clear.genstat == 1 && srcloc == 1
                    FTdata.dirstats.npd = [];
                end
                intCoh = [];
                if isempty(strmatch(R.sourcenames{srcloc},'fEEG'))
                    for dir = 1:3
                        for i = 1:length(inds)
                            set1 = inds(i);
                            set2 = find(strncmpi(['fEEG'],FTdata.npd.label,3)); % STN index
                            fy = squeeze(FTdata.npd.npdspctrm{dir}(set1,set2,:)); % Combined Coherence spectra
                            if size(fy,1)>1 && size(fy,2)>1
                                fy = mean(fy,1);
                            end
                            
                            fxA = FTdata.npd.freq; % Frequency X vector
                            res = fxA(2) - fxA(1);  % Freq resolution
                            for band = 1:size(R.bbounds,1)
                                t_intCoh(band,i) = abs(res*trapz(fy(round((bbounds(band,1)-R.FOI(1))/res)+1:round((bbounds(band,2)-R.FOI(1))/res)+1)));
                                [t_peakMag(band,i) frq] = max(fy(round((bbounds(band,1)-R.FOI(1))/res)+1:round((bbounds(band,2)-R.FOI(1))/res)+1));
                                t_peakFrq(band,i) =  ((frq-1)*res)+bbounds(band,1);
                            end
                        end
                        intCoh{dir} = t_intCoh;
                        peakMag{dir} = t_peakMag;
                        peakFrq{dir} = t_peakFrq;
                    end
                    FTdata.dirstats.npd.(R.sourcenames{srcloc}).frwd_intCoh = intCoh{2};    % Integrated Coherence
                    FTdata.dirstats.npd.(R.sourcenames{srcloc}).frwd_peakMag = peakMag{2};  % peak magnitude within band
                    FTdata.dirstats.npd.(R.sourcenames{srcloc}).frwd_peakFrq = peakFrq{2};  % peak frequency (Hz)
                    FTdata.dirstats.npd.(R.sourcenames{srcloc}).rev_intCoh = intCoh{3};    % Integrated Coherence
                    FTdata.dirstats.npd.(R.sourcenames{srcloc}).rev_peakMag = peakMag{3};  % peak magnitude within band
                    FTdata.dirstats.npd.(R.sourcenames{srcloc}).rev_peakFrq = peakFrq{3};  % peak frequency (Hz)
                    
                else
                    FTdata.dirstats.npd.(R.sourcenames{srcloc}).frwd_intCoh = [];    % Integrated Coherence
                    FTdata.dirstats.npd.(R.sourcenames{srcloc}).rev_intCoh = [];   % peak magnitude within band
                    %                     FTdata.dirstats.npd.(R.sourcenames{srcloc}).peakFrq = [];   % peak frequency (Hz)
                end
            end
            
            
            %% Save Output
            save([R.analysispath R.pipestamp '\data\processed\' R.subnames{cond}{sub} '_' R.condnames{cond} '_' R.pipestamp '.mat'],'FTdata')
        end
    end
end
