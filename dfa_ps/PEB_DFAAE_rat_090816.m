function [] = PEB_DFAAE_rat_090816(R)
numbands = length(R.bandnames);
for cond = 1:2
    for sub  = 1:length(R.subnames{cond})
        load([R.analysispath R.pipestamp '\data\processed\' R.subnames{cond}{sub} '_' R.condnames{cond} '_' R.pipestamp '.mat'])
        for srcloc = 1:length(R.sourcenames)
            DFAStore = []; boxStore = []; pebStore = [];
            
            xclean = FTdata.ContData.trial{1};
            inds = find(strncmpi([R.sourcenames{srcloc}], FTdata.freqPow.label,length(R.sourcenames{srcloc})));
            powfreq = FTdata.dirstats.pow.(R.sourcenames{srcloc}).peakFrq;
            bwid = R.dfaae.bwid; BF_r = R.dfaae.BF_r; fsamp = FTdata.fsample;
            if ~isempty(inds)
                parfor i = 1:length(inds)
                    for band = 1:numbands
                        x1 = xclean(inds(i),:);
                        % Find PS-DFA
                        cfreq = powfreq(band,i);
                        lf = cfreq - bwid(band); hf = cfreq + bwid(band);
                        [DFAStore(i,band) boxStore(:,i,band) pebStore(:,i,band)] = dfaanaly(x1,lf,hf,fsamp,BF_r);
                    end
                end
                
                FTdata.DFAAE.(R.sourcenames{srcloc}).DFAStore = DFAStore;
                FTdata.DFAAE.(R.sourcenames{srcloc}).boxStore = squeeze(mean(mean(boxStore,3),2));
                FTdata.DFAAE.(R.sourcenames{srcloc}).peb = reshape(pebStore(1,:,:),[],3);
            end
            save([R.analysispath R.pipestamp '\data\processed\' R.subnames{cond}{sub} '_' R.condnames{cond} '_' R.pipestamp '.mat'],'FTdata')
            disp([sub cond srcloc])
            clear DFAStore boxStore pebStore
        end
    end
end
end
function [A B C] = dfaanaly(x1,lf,hf,fsamp,BF_r)
x1_filt = filterEEG(x1,fsamp,lf,hf,round(6*fsamp/lf)); % bandpass filter
x1AE = abs(hilbert(x1_filt));
maxFrac = 8; minBS = (1/lf)*12;
DFAP = [fsamp minBS (length(x1AE)/maxFrac)/fsamp 100 0];
[bmod win evi alpha] = peb_dfa_cohproj_090616(x1AE,DFAP,BF_r,0);
if  evi < BF_r
    rej = 1;
    %     alpha = NaN;
else
    rej = 0;
end
A = [alpha]; B = [lf hf minBS]; C = [evi rej];
end