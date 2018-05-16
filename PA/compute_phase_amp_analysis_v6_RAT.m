function compute_phase_amp_analysis_v6_RAT(R,idd)
clear
R = buildheader_rat;

for band = 2:3; %[1 3]
    for cond =1:2
        for sub  = 1:length(R.subnames{cond})
            load([R.analysispath R.pipestamp '\data\processed\' R.subnames{cond}{sub} '_' R.condnames{cond} '_' R.pipestamp '.mat'])
            stnlabs = find(strncmp('STN',FTdata.nsIcoh.label,3));
            % %             stncoh = []; frqcoh = [];
            % %             for i = 1:length(stnlabs)
            % %                [stncoh(i), frqcoh(i)] = max(FTdata.nsIcoh.nsIcohspctrm(1,stnlabs(i),FTdata.nsIcoh.freq > R.bbounds(band,1) & FTdata.nsIcoh.freq <= R.bbounds(band,2)));
            % %             end
            % %             stnind = find(stncoh == max(stncoh));
            % %             frq =  FTdata.nsIcoh.freq(FTdata.nsIcoh.freq > R.bbounds(band,1) & FTdata.nsIcoh.freq <= R.bbounds(band,2));
            % %             frq = frq(stnind);
            
            phi_dist_col = []; amp_dist_col = []; segL_ddt_col = []; tend = 0; H_col = [];
            for i = 1:length(stnlabs)
                stnind = stnlabs(i);
                stncoh = []; frqcoh = [];
                [stncoh, frqcoh] = max(FTdata.nsIcoh.nsIcohspctrm(1,stnind,FTdata.nsIcoh.freq > R.bbounds(band,1) & FTdata.nsIcoh.freq <= R.bbounds(band,2)));
                frq =  FTdata.nsIcoh.freq(FTdata.nsIcoh.freq > R.bbounds(band,1) & FTdata.nsIcoh.freq <= R.bbounds(band,2));
                frq = frq(frqcoh);
                
                % find lb
                [~, lbfrqcoh] = max(FTdata.nsIcoh.nsIcohspctrm(1,stnind,FTdata.nsIcoh.freq > R.bbounds(band-1,1) & FTdata.nsIcoh.freq <= R.bbounds(band-1,2)));
                lbfrq =  FTdata.nsIcoh.freq(FTdata.nsIcoh.freq > R.bbounds(band-1,1) & FTdata.nsIcoh.freq <= R.bbounds(band-1,2));
                lbfrq = lbfrq(lbfrqcoh);
                % Create data structure
                Xdata.fsample = FTdata.ContData.fsample;
                Xdata.label = {'M2','STN1','STN2'};
                Xdata.trial{1} = FTdata.ContData.trial{1}([1 stnind stnind],:);
                Xdata.time{1} = FTdata.ContData.time{1}; %{1}
                
                % Get the overal signal amplitudes
                signalEnvAmp = median(abs(hilbert(Xdata.trial{1})),2);
                maxfrq = frq;
                
                R.pp.cont.full.fs = FTdata.ContData.fsample; R.PA.stn_lb_frq = lbfrq;
                
                [SRPeps Ampeps SNR_eps] = phase_amp_surrComp(Xdata,R,maxfrq,100);
                SNR_eps_z(1) = 10*log10(SNR_eps(1)./signalEnvAmp(1));
                SNR_eps_z(2) =  10*log10(SNR_eps(2)./signalEnvAmp(2));
                SNR_eps_z(3) =  10*log10(SNR_eps(3)./signalEnvAmp(2));
                % Compute data transforms (Hilbert)
                [amp phi dphi_12 dphi_12_dt betaS] = comp_instant_angle_phase(Xdata,maxfrq,R.PA.stn_lb_frq,R.PA.bwid(band),Xdata.fsample,R.PA.LowAmpFix);
                
                % Sliding Window PLV
                if R.PA.SType == 1
                    WinSize = R.PA.slidingwindow*R.pp.cont.full.fs;
                    [PLV PLV_tvec] = slidingwindowPLV(WinSize,phi,R.PA.WinOver);
                    PLV_tvec = Xdata.time{1}(PLV_tvec);
                    amp_sw = cont2slidingwindow(amp,WinSize,floor(R.PA.WinOver*WinSize));
                    clear snr_sw
                    snr_sw(:,1) = log10(amp_sw(:,1)./signalEnvAmp(1));
                    snr_sw(:,2) =  log10(amp_sw(:,2)./signalEnvAmp(2));
                    snr_sw(:,3) =  log10(amp_sw(:,3)./signalEnvAmp(2));
                    dphi_12_sw = cont2slidingwindow(dphi_12,WinSize,round(R.PA.WinOver*WinSize));
                    
                    SW_sampr = max(diff(PLV_tvec));
                    tseries = dphi_12_sw; qstable = find(PLV>R.PA.PLVeps);
                    mwid = R.PA.mwid;
                    period = (mwid/frq)/SW_sampr;
                    % Minimum number of cycles to consider sync
                    [phi_dist amp_dist seg_ddt segL_ddt consecSegs H] = analysestablesegs(qstable,tseries,snr_sw,period,mwid,1/SW_sampr,R.PA.SNR(band),Ampeps);
                    
                elseif R.PA.SType == 2 % Sliding Window PhaseAng. Stability
                    clear snr_sw
                    snr_sw(:,1) = 10.*log10(amp(:,1)./signalEnvAmp(1));
                    snr_sw(:,2) =  10.*log10(amp(:,2)./signalEnvAmp(2));
                    snr_sw(:,3) =  10.*log10(amp(:,3)./signalEnvAmp(2));
                    mwid = R.PA.mwid; % Minimum number of cycles to consider sync
                    period = (mwid/maxfrq)*R.pp.cont.full.fs;
                    cycle = (1/maxfrq)*R.pp.cont.full.fs;
                    
                    %                             b = (1/floor(period))*ones(1,floor(period));
                    %                             a = 1;
                    %                             dphi_12_dt = filter(b,a,dphi_12_dt);
                    
                    %%% Find segments
                    %                 dphi_12_dt_sm = smooth(dphi_12_dt,period/6)';
                    qstable = find(abs(dphi_12_dt')<SRPeps); % The points that are below threshold
                    tseries = wrapToPi(dphi_12(2:end));
                    %                 tseries(:,2) = phi(:,2);
                    %                 tseries(:,1) = phi(:,1);
                    [phi_dist amp_dist seg_ddt1 segL_ddt consecSegs H] = analysestablesegs(qstable,tseries,amp,period,mwid,R.pp.cont.full.fs,SNR_eps_z(1:2),[],[],Ampeps,snr_sw);
                    %                                                 plot_example_phaseanalysis_trace(betaS,amp,phi,dphi_12_dt,seg_ddt1,SRPeps,R.pp.cont.full.fs);
                end
                tend = tend + (Xdata.time{1}(end)- Xdata.time{1}(1));
                %             phi_dist = wrapToPi(phi_dist-circ_mean(phi_dist(~isnan(phi_dist))'));
                phi_dist_col = [phi_dist_col phi_dist];
                amp_dist_col = [amp_dist_col amp_dist];
                segL_ddt_col = [segL_ddt_col segL_ddt];
                H_col = [H_col H];
            end
            
            %                 [gc_dist] = analysestablesegs_granger(qstable,tseries,amp,period,mwid,R.pp.cont.full.fs);
            FTdata.PA.gc_dist_save{1} = NaN; %gc_dist;
            FTdata.PA.H_dist_save{1} = H_col;
            FTdata.PA.pA_pli_dist_save{1} = phi_dist_col;
            FTdata.PA.amp_pli_dist_save{1} = amp_dist_col;
            FTdata.PA.segL_pli_dist_save{1} = segL_ddt_col;
            FTdata.PA.timevec{1} = tend; %{1}
            
            save([R.analysispath R.pipestamp '\data\processed\' R.subnames{cond}{sub} '_' R.condnames{cond} '_' R.pipestamp '_' R.bandnames{band} '.mat'],'FTdata')
            %                     figure
            %                     PLV_sw_plot(Xdata,betaS,amp,phi,snr_sw,seg_ddt,PLV,PLV_tvec,consecSegs,R)
            %                     %                 savefigure_v2([R.datapathr 'results\images\seganalysis\'],['example_seg_subject1_ON'],[],[],[]);
            %                                         close all
        end
    end
end

% ! shutdown /h