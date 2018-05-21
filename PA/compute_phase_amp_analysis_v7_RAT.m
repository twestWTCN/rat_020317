function compute_phase_amp_analysis_v7_RAT(R)
if nargin<2
    R = buildheader_rat();
end
%%%
close all
for band = [2 3]
    for cond =1:2
        for sub  = 1:length(R.subnames{cond})
            load([R.analysispath R.pipestamp '\data\processed\' R.subnames{cond}{sub} '_' R.condnames{cond} '_' R.pipestamp '.mat'])
            stnlabs = find(strncmp('STN',FTdata.wpli.label,3));
            ctxlab = find(strncmp('M1',FTdata.wpli.label,3));
            stncoh = []; frqcoh = [];
            for stchi = 1:length(stnlabs)
                stnind = stnlabs(stchi);
                [stncoh, frqcoh] = max(FTdata.wpli.wpli_debiasedspctrm(1,stnind,FTdata.wpli.freq > R.bbounds(band,1) & FTdata.wpli.freq <= R.bbounds(band,2)));
                frq =frqcoh(stncoh == max(stncoh))+ R.bbounds(band,1);
                Xdata.time = FTdata.time;
                Xdata.trial{1} = FTdata.trial{1}([ctxlab stnind],:);
                Xdata.label = {'M2','STN'};
                Xdata.fsample = FTdata.fsample;
                % Find sub-band frequency
                [stncoh, frqcoh] = max(FTdata.wpli.wpli_debiasedspctrm(1,stnind,FTdata.wpli.freq > R.bbounds(band-1,1) & FTdata.wpli.freq <= R.bbounds(band-1,2)));
                stn_lb_frq =frqcoh(stncoh == max(stncoh))+ R.bbounds(band-1,1);
                
                % Create data structure
                % Get the overal signal amplitudes
                sig =  Xdata.trial{1};
                signalEnvAmp = median(abs(hilbert(sig)),2);
                
                % Compute optimal frequencies based on PLV
                [maxfrq,maxPLV] = compute_optimal_PhaseLockFrq(R,Xdata,band,stn_lb_frq); %
                [SRPeps Ampeps SNR_eps PLVeps] = phase_amp_surrComp(R,Xdata,band,maxfrq,stn_lb_frq,R.PA.AmpSurrN);
                SNR_eps_z(1) = 10*log10(SNR_eps(1)./signalEnvAmp(1));
                SNR_eps_z(2) =  10*log10(SNR_eps(2)./signalEnvAmp(2));
                SNR_eps_z(3) =  10*log10(SNR_eps(3)./signalEnvAmp(2));
                % Compute data transforms (Hilbert)
                [amp phi dphi_12 dphi_12_dt Xdata] = comp_instant_angle_phase(Xdata,maxfrq,stn_lb_frq,R.PA.bwid,band);
                % Sliding Window PLV
                if R.PA.SType == 1
                    fsamp = Xdata.fsample;
                    WinSize = floor(R.PA.slidingwindow*fsamp);
                    [PLV PLV_tvec] = slidingwindowPLV(WinSize,phi,R.PA.WinOver);
                    PLV_tvec = Xdata.time{1}(PLV_tvec);
                    amp_sw = cont2slidingwindow(amp,WinSize,floor(R.PA.WinOver*WinSize));
                    clear SNR_sw
                    SNR_sw(:,1) = 10.*log10(amp_sw(:,1)./signalEnvAmp(1));
                    SNR_sw(:,2) =  10.*log10(amp_sw(:,2)./signalEnvAmp(2));
                    SNR_sw(:,3) =  10.*log10(amp_sw(:,3)./signalEnvAmp(2));
                    dphi_12_sw = cont2slidingwindow(dphi_12,WinSize,round(R.PA.WinOver*WinSize));
                    
                    SW_sampr = max(diff(PLV_tvec));
                    tseries = dphi_12_sw; qstable = find(PLV>PLVeps);
                    mwid = R.PA.mwid;
                    period = (mwid/frq)/SW_sampr;
                    % Minimum number of cycles to consider sync
                    [phi_dist{stchi} amp_dist{stchi} seg_ddt1{stchi} segL_ddt{stchi} consecSegs H] = analysestablesegs(qstable,tseries,amp,period,mwid,1/SW_sampr,SNR_eps_z(1:2),[],[],Ampeps,SNR_sw);
                    %                  plot_example_phaseanalysis_SW(Xdata,amp,phi,PLV,seg_ddt1{i},PLVeps,PLV_tvec);
                elseif R.PA.SType == 2 % Sliding Window PhaseAng. Stability
                    clear SNR_sw
                    SNR_sw(:,1) = 10.*log10(amp(:,1)./signalEnvAmp(1));
                    SNR_sw(:,2) =  10.*log10(amp(:,2)./signalEnvAmp(2));
                    SNR_sw(:,3) =  10.*log10(amp(:,3)./signalEnvAmp(2));
                    % Plot SNR
                    %                 figure(1)
                    %                 SNR_Inspector(R,Xdata.time{1},SNR_sw,SNR_eps_z,band,{'M2','STN'})
                    
                    %%% Set Length Constraints
                    fsamp = Xdata.fsample;
                    mwid = R.PA.mwid; % Minimum number of cycles to consider sync
                    period = (mwid/maxfrq)*fsamp;
                    
                    %%% Find segments
                    %                 dphi_12_dt_sm = smooth(dphi_12_dt,period/6)';
                    qstable = find(abs(dphi_12_dt')<SRPeps); % The points that are below threshold
                    tseries = wrapToPi(dphi_12(2:end));
                    [phi_dist{stchi} amp_dist{stchi} seg_ddt1{stchi} segL_ddt{stchi} consecSegs H] = analysestablesegs(qstable,tseries,amp,period,mwid,fsamp,SNR_eps_z(1:2),[],[],Ampeps,SNR_sw);
                    %                 figure(2)
                    %                 plot_example_phaseanalysis_trace(Xdata,amp,phi,dphi_12_dt,seg_ddt1,SRPeps,fsamp);
                    %                 clf(1);
                end
            end
            % Save to data
            FTdata.PA(band).pA_pli_dist_save =phi_dist;
            FTdata.PA(band).amp_pli_dist_save = amp_dist;
            FTdata.PA(band).segL_pli_dist_save = segL_ddt;
            FTdata.PA(band).timevec = Xdata.time;
            %             FTdata.PA(band).H_dist_save{cond} = H;
            %             FTdata.history = [FTdata.history{:} {[mfilename '_' date]}];
            clear phi_dist amp_dist segL_ddt
            save([R.analysispath R.pipestamp '\data\processed\' R.subnames{cond}{sub} '_' R.condnames{cond} '_' R.pipestamp '.mat'],'FTdata')
        end
    end
end


% ! shutdown /h