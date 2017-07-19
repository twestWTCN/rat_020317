function [] = PEB_DFAPS_rat_011116(R)
numbands = length(R.bandnames);
for cond = 1:2
    for sub  = 1:length(R.subnames{cond})
        load([R.analysispath R.pipestamp '\data\processed\' R.subnames{cond}{sub} '_' R.condnames{cond} '_' R.pipestamp '.mat'])
        for srcloc = 1:length(R.sourcenames)
            DFAStore = []; boxStore = []; pebStore = [];           
            xclean = FTdata.ContData.trial{1};
            
            set1 = find(strncmpi([R.sourcenames{srcloc}], FTdata.freqPow.label,length(R.sourcenames{srcloc})));
            set2 = find(strncmp(['fEEG'],FTdata.freqPow.label,3));
            [p,q] = meshgrid(set1, set2);
            pairs = [p(:) q(:)];
            pairs(pairs(:,2)==pairs(:,1),:) = []; % delete matching pairs
            bwid = R.dfaps.bwid; BF_r = R.dfaps.BF_r; fsamp = FTdata.fsample;
            if ~isempty(pairs)
                cohfreq = cenfreqcalc(FTdata,R.dfaps.method.cfreq,pairs,R);
                for i = 1:size(pairs,1)
                    parfor band = 1:numbands
                        x1 = xclean(pairs(i,1),:); x2 = xclean(pairs(i,2),:);
                        %                         cfreq = cohfreq(band,i);
                        cfreq = cohfreq(i,band);
                        [sub cond srcloc]
                        [i band]
                        lf = cfreq - bwid(band); hf = cfreq + bwid(band);
                        [DFAStore(i,band) boxStore(:,i,band) pebStore(:,i,band) varfeat(i,band) ARCoeff(i,band)] = dfaanaly(x1,x2,lf,hf,fsamp,BF_r);
                    end
                end
                if isempty(pebStore)
                    pebStore = [];
                    FTdata.DFAPS.(R.sourcenames{srcloc}).peb = [];
                else
                    FTdata.DFAPS.(R.sourcenames{srcloc}).peb = reshape(pebStore(1,:,:),[],3);
                end
                FTdata.DFAPS.(R.sourcenames{srcloc}).DFAStore = DFAStore;
                FTdata.DFAPS.(R.sourcenames{srcloc}).boxStore = squeeze(mean(mean(boxStore,4),3));
                FTdata.DFAPS.(R.sourcenames{srcloc}).varfeat = varfeat;
                 FTdata.DFAPS.(R.sourcenames{srcloc}).ARCoeff = ARCoeff;
            end
            clear DFAStore
        end
        if ~isfield(FTdata,'DFAPS')
            FTdata.DFAPS = 'empty';
        end
            save([R.analysispath R.pipestamp '\data\processed\' R.subnames{cond}{sub} '_' R.condnames{cond} '_' R.pipestamp '.mat'],'FTdata')
            disp([sub cond srcloc])
    end
end
end
function [A B C D E] = dfaanaly(x1,x2,lf,hf,fsamp,BF_r)
x1_filt = filterEEG(x1,fsamp,lf,hf,3*fix(fsamp/lf)); % bandpass filter
x2_filt = filterEEG(x2,fsamp,lf,hf,3*fix(fsamp/lf)); % bandpass filter % Changes from 3 sample to 6 samples
dphdiff = psdfasigTrans(x1_filt,x2_filt);
maxFrac = 8; minBS = (1/lf)*3;
DFAP = [fsamp minBS (length(dphdiff)/maxFrac)/fsamp 50 0];

[bmod win evi alpha] = peb_dfa_cohproj_090616(dphdiff,DFAP,BF_r,0);

if  evi < BF_r
    rej = 1;
    %     alpha = NaN;
else
    rej = 0;
end

m = ar(dphdiff,1);
A = [alpha]; B = [lf hf minBS]; C = [evi rej]; D = std(dphdiff); E = getpvec(m);
end