function [] = PEB_DFAPS_rat_090816(R)
numbands = length(R.bandnames);
for cond = 1:2
    for sub  = 1:length(R.subnames{cond})
        load([R.analysispath R.pipestamp '\data\processed\' R.subnames{cond}{sub} '_' R.condnames{cond} '_' R.pipestamp '.mat'])
        for srcloc = 1:length(R.sourcenames)
            DFAStore = []; boxStore = []; pebStore = [];
            
            xclean = FTdata.ContData.trial{1};
            cohfreq = mean(R.bbounds,2);    % use centre of band
            set1 = find(strncmpi([R.sourcenames{srcloc}], FTdata.freqPow.label,length(R.sourcenames{srcloc})));
            set2 = find(strncmp(['STN'],FTdata.freqPow.label,3));
            [p,q] = meshgrid(set1, set2);
            pairs = [p(:) q(:)];
            pairs(pairs(:,2)==pairs(:,1),:) = []; % delete matching pairs
            bwid = R.dfaps.bwid; BF_r = R.dfaps.BF_r; fsamp = FTdata.fsample;

            if ~isempty(set1)
                for i = 1:size(pairs,1)
                    parfor band = 1:numbands
                        x1 = xclean(pairs(i,1),:); x2 = xclean(pairs(i,2),:);
                        %                         cfreq = cohfreq(band,i);
                        cfreq = cohfreq(band);
                        [sub cond srcloc]
                        [i band]
                        lf = cfreq - bwid(band); hf = cfreq + bwid(band);
                        [DFAStore(i,band) boxStore(:,i,band) pebStore(:,i,band)] = dfaanaly(x1,x2,lf,hf,fsamp,BF_r);
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
                
            end
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

if  evi < BF_r
    rej = 1;
    %     alpha = NaN;
else
    rej = 0;
end
A = [alpha]; B = [lf hf minBS]; C = [evi rej];
end