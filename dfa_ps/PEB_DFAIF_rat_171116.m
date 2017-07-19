function [] = PEB_DFAIF_rat_171116(R)
numbands = length(R.bandnames);
for cond = 1:2
    for sub  = 1:length(R.subnames{cond})
        load([R.analysispath R.pipestamp '\data\processed\' R.subnames{cond}{sub} '_' R.condnames{cond} '_' R.pipestamp '.mat'])
        for srcloc = 1:length(R.sourcenames)
            DFAStore = []; boxStore = []; pebStore = [];
            
            xclean = FTdata.ContData.trial{1};
            inds = find(strncmpi([R.sourcenames{srcloc}], FTdata.freqPow.label,length(R.sourcenames{srcloc})));
            bwid = R.dfaif.bwid; BF_r = R.dfaif.BF_r; fsamp = FTdata.fsample;
            if ~isempty(inds)
                powfreq = cenfreqcalc(FTdata,R.dfaif.method.cfreq,inds,R);
                for i = 1:length(inds)
                    for band = 1:numbands
                        x1 = xclean(inds(i),:);
                        % Find PS-DFA
                        cfreq = powfreq(i,band);
                        lf = cfreq - bwid(band); hf = cfreq + bwid(band);
                        [DFAStore(i,band) boxStore(:,i,band) pebStore(:,i,band) varfeat(i,band) ARCoeff(i,band)] = dfaanaly(x1,lf,hf,fsamp,BF_r);
                    end
                end
                
                FTdata.DFAIF.(R.sourcenames{srcloc}).DFAStore = DFAStore;
                FTdata.DFAIF.(R.sourcenames{srcloc}).boxStore = squeeze(mean(mean(boxStore,3),2));
                FTdata.DFAIF.(R.sourcenames{srcloc}).peb = reshape(pebStore(1,:,:),[],3);
                FTdata.DFAIF.(R.sourcenames{srcloc}).varfeat = varfeat;
                FTdata.DFAIF.(R.sourcenames{srcloc}).ARCoeff = ARCoeff;
                
            end
            save([R.analysispath R.pipestamp '\data\processed\' R.subnames{cond}{sub} '_' R.condnames{cond} '_' R.pipestamp '.mat'],'FTdata')
            disp([sub cond srcloc])
            clear DFAStore boxStore pebStore
        end
    end
end
end
function [A B C D E] = dfaanaly(x1,lf,hf,fsamp,BF_r)
close all
x1_filt = filterEEG(x1,fsamp,lf,hf,round(6*fsamp/lf)); % bandpass filter
x1AE = normaliseV(abs(hilbert(x1_filt)));
x1IF = x1AE(2:end).*diff(unwrap(angle(hilbert(x1_filt))));
x2 = diff(unwrap(angle(hilbert(x1_filt))))
maxFrac = 10; minBS = (1/lf)*12;
DFAP = [fsamp minBS (length(x1IF)/maxFrac)/fsamp 100 1];
[bmod win evi alpha] = peb_dfa_cohproj_090616(x1IF,DFAP,BF_r,1);
figure
plot(x1IF,'r'); hold on; plot(x1AE); plot(x2,'g')
legend({'Amp Normalised dPhi/dt','Normalised Amplitude Envelope','dPhi/dt'})

if  evi < BF_r
    rej = 1;
    %     alpha = NaN;
else
    rej = 0;
end
m = ar(x1AE,1);
A = [alpha]; B = [lf hf minBS]; C = [evi rej]; D = std(x1AE); E = getpvec(m);
end