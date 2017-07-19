function [] = PEB_DFAAE_rat_011116(R)
numbands = length(R.bandnames);
for cond = 1:2
    for sub  = 1:length(R.subnames{cond})
        load([R.analysispath R.pipestamp '\data\processed\' R.subnames{cond}{sub} '_' R.condnames{cond} '_' R.pipestamp '.mat'])
        for srcloc = 1:length(R.sourcenames)
            DFAStore = []; boxStore = []; pebStore = [];
            
            xclean = FTdata.ContData.trial{1};
            inds = find(strncmpi([R.sourcenames{srcloc}], FTdata.freqPow.label,length(R.sourcenames{srcloc})));
            bwid = R.dfaae.bwid; BF_r = R.dfaae.BF_r; fsamp = FTdata.fsample;
            if ~isempty(inds)
                powfreq = cenfreqcalc(FTdata,R.dfaae.method.cfreq,inds,R);
                for i = 1:length(inds)
                    parfor band = 1:numbands
                        x1 = xclean(inds(i),:);
                        % Find PS-DFA
                        cfreq = powfreq(i,band);
                        lf = cfreq - bwid(band); hf = cfreq + bwid(band);
                        [DFAStore(i,band) boxStore(:,i,band) pebStore(:,i,band) varfeat(i,band) ARCoeff(i,band)] = dfaanaly(x1,lf,hf,fsamp,BF_r);
                    end
                end
                
                FTdata.DFAAE.(R.sourcenames{srcloc}).DFAStore = DFAStore;
                FTdata.DFAAE.(R.sourcenames{srcloc}).boxStore = squeeze(mean(mean(boxStore,3),2));
                FTdata.DFAAE.(R.sourcenames{srcloc}).peb = reshape(pebStore(1,:,:),[],3);
                FTdata.DFAAE.(R.sourcenames{srcloc}).varfeat = varfeat;
                FTdata.DFAAE.(R.sourcenames{srcloc}).ARCoeff = ARCoeff;
                
            end
            clear DFAStore boxStore pebStore
        end
        save([R.analysispath R.pipestamp '\data\processed\' R.subnames{cond}{sub} '_' R.condnames{cond} '_' R.pipestamp '.mat'],'FTdata')
        disp([sub cond srcloc])
        
    end
end
end
function [A B C D E] = dfaanaly(x1,lf,hf,fsamp,BF_r)
% close all
x1_filt = filterEEG(x1,fsamp,lf,hf,6*fix(fsamp/lf)); % bandpass filter
x1AE = abs(hilbert(x1_filt));
maxFrac = 10; minBS = (1/lf)*6;
DFAP = [fsamp minBS (length(x1AE)/maxFrac)/fsamp 50 0];
[bmod win evi alpha] = peb_dfa_cohproj_090616(x1AE,DFAP,BF_r,0);
if  evi < BF_r
    rej = 1;
    %     alpha = NaN;
else
    rej = 0;
end
m = ar(x1AE,1);
A = [alpha]; B = [lf hf minBS]; C = [evi rej]; D = std(x1AE); E = getpvec(m);
end