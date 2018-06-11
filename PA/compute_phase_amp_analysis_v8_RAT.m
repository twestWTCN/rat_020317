function compute_phase_amp_analysis_v8_RAT(R)
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
                [phi_dist{stchi} amp_dist{stchi} segL_ddt{stchi}] = compute_dynPhaseLocking(R,Xdata,band,stn_lb_frq)
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