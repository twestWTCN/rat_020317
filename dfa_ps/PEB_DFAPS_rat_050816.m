function [] = PEB_DFAPS_rat_050816(R)
numbands = length(R.bandnames);
for cond = 1:2
    for sub  = 1:length(R.subnames{cond})
        load([R.analysispath R.pipestamp '\data\processed\' R.subnames{cond}{sub} '_' R.condnames{cond} '_' R.pipestamp '.mat'])
        for srcloc = 2:length(R.sourcenames)
            DFAStore = []; boxStore = []; pebStore = [];
            
            xclean = FTdata.ContData.trial{1};
            inds = find(strncmpi([R.sourcenames{srcloc}], FTdata.freqPow.label,length(R.sourcenames{srcloc})));
            %             cohfreq = FTdata.dirstats.wpli.(R.sourcenames{srcloc}).peakFrq;   % if
            %             using WPLI coherent freq
            cohfreq = mean(R.bbounds,2);    % use centre of band
            for i = 1:length(inds)
                set1 = inds(i);
                set2 = find(strncmp(['STN'],FTdata.freqPow.label,3));
                for j = 1:length(set2)
                    for band = 1:numbands
                        x1 = xclean(set1,:); x2 = xclean(set2(j),:);
                        % Find PS-DFA
                        %                         cfreq = cohfreq(band,i);
                        cfreq = cohfreq(band);
                        lf = cfreq - R.psdfa.bwid(band); hf = cfreq + R.psdfa.bwid(band);
                        [DFAStore(i,j,band) boxStore(:,i,j,band) pebStore(:,i,j,band)] = dfaanaly(x1,x2,lf,hf,FTdata.fsample,R.psdfa.BF_r);
                    end
                end
            end
            
            FTdata.DFAPS.(R.sourcenames{srcloc}).DFAStore = reshape(DFAStore,[],3);
            FTdata.DFAPS.(R.sourcenames{srcloc}).boxStore = squeeze(mean(mean(boxStore,4),3));
            FTdata.DFAPS.(R.sourcenames{srcloc}).peb = reshape(pebStore(1,:,:,:),[],3);
            
            save([R.analysispath R.pipestamp '\data\processed\' R.subnames{cond}{sub} '_' R.condnames{cond} '_' R.pipestamp '.mat'],'FTdata')
            clear DFAStore boxStore pebStore
            disp([sub cond srcloc])
        end
    end
end
end
function [A B C] = dfaanaly(x1,x2,lf,hf,fsamp,BF_r)
x1_filt = filterEEG(x1,fsamp,lf,hf,round(6*fsamp/lf)); % bandpass filter
x2_filt = filterEEG(x2,fsamp,lf,hf,round(6*fsamp/lf)); % bandpass filter % Changes from 3 sample to 6 samples
dphdiff = psdfasigTrans(x1_filt,x2_filt);
maxFrac = 8; minBS = (1/lf)*12;
DFAP = [fsamp minBS (length(dphdiff)/maxFrac)/fsamp 100 0];
[bmod win evi alpha] = peb_dfa_cohproj_090616(dphdiff,DFAP,BF_r,0);
if isinf(evi)
    pause
end
if  evi < BF_r
    rej = 1;
    %     alpha = NaN;
else
    rej = 0;
end
A = [alpha]; B = [lf hf minBS]; C = [evi rej];
end