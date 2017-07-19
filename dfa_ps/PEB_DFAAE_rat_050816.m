function [] = PEB_DFAAE_rat_050816(R)
numbands = length(R.bandnames);
for cond = 1:2
    for sub  = 1:length(R.subnames{cond})
        load([R.analysispath R.pipestamp '\data\processed\' R.subnames{cond}{sub} '_' R.condnames{cond} '_' R.pipestamp '.mat'])
        for srcloc = 1:length(R.sourcenames)
            xclean = FTdata.ContData.trial{1};
            inds = find(strncmpi([R.sourcenames{srcloc}], FTdata.freqPow.label,size(R.sourcenames{srcloc},2)));
            powfreq = FTdata.dirstats.pow.(R.sourcenames{srcloc}).peakFrq;
            parfor i = 1:length(inds)
                for band = 1:numbands
                    x1 = xclean(inds(i),:);
                    % Find PS-DFA
                    cfreq = powfreq(band,i);
                    lf = cfreq - R.dfa.bwid(band); hf = cfreq + R.dfa.bwid(band);
                    [DFAStore(band,i) boxStore(:,band,i) pebStore(:,band,i)] = dfaanaly(x1,lf,hf,FTdata.fsample,R.dfa.BF_r);
                end
            end
            
            if isempty(inds)
                DFAStore = []; boxStore = []; pebStore = [];
            else
                for band = 1:numbands
                    rejProp(band) = sum(pebStore(3,band,:,:))/length(inds);
                    rejN(band) = sum(pebStore(3,band,:));
                end
            end
            
            FTdata.DFAAE.(R.sourcenames{srcloc}).DFAStore = DFAStore;
            FTdata.DFAAE.(R.sourcenames{srcloc}).boxStore = boxStore;
            FTdata.DFAAE.(R.sourcenames{srcloc}).peb = pebStore;
            FTdata.DFAAE.(R.sourcenames{srcloc}).pebrej = [rejProp; rejN];
            clear DFAStore boxStore pebStore
        end
        save([R.analysispath R.pipestamp '\data\processed\' R.subnames{cond}{sub} '_' R.condnames{cond} '_' R.pipestamp '.mat'],'FTdata')
        disp([sub cond srcloc])
    end
end
end
function [A B C] = dfaanaly(x1,lf,hf,fsamp,BF_r)
x1_filt = filterEEG(x1,fsamp,lf,hf,round(6*fsamp/lf)); % bandpass filter
x1AE = abs(hilbert(x1_filt));
maxFrac = 8; minBS = (1/lf)*12;
DFAP = [fsamp minBS (length(x1AE)/maxFrac)/fsamp 100 0];
[bmod win evi alpha] = peb_dfa_cohproj_090616(x1AE,DFAP,BF_r,0);
if  evi(2) < BF_r
    rej = 1;
    alpha = NaN;
else
    rej = 0;
    dfaparam = [fsamp minBS (length(x1AE)/maxFrac)/fsamp 100 1];
end
A = [alpha]; B = [lf hf minBS]; C = [evi rej];
end