%function clusterspec_020315(R)
%% Coherence Statistics- Integrated coherence, peak magnitude, WPLI
% if R.clear.specstat == 1
%     FTdata.dirstats.pow = [];
% end
bbounds = R.bbounds;
for cond = 1:2
    for sub  = 1:length(R.subnames{cond})
        load([R.analysispath R.pipestamp '\data\processed\' R.subnames{cond}{sub} '_' R.condnames{cond} '_' R.pipestamp '.mat'])
        for srcloc = 1:length(R.sourcenames)
            clear fy
            inds = find(strncmpi([R.sourcenames{srcloc} ' '], FTdata.freqPow.label,3)); % Source index
            intCoh = []; peakMag = []; peakFrq = [];
            for i = 1:length(inds)
                set1 = inds(i);
%                 fy(i,:) = squeeze(FTdata.freqPow.powspctrm(set1,:)); % Combined Coherence spectra
%                 fy(i,:) = squeeze(FTdata.spectanaly.spectra{1}(set1,:)); % Combined Coherence spectra
                fy(i,:) = squeeze(FTdata.wpli.wpli_debiasedspctrm(1,set1,:));
%                 fy(i,:) = squeeze(abs(FTdata.coherence.cohspctrm(1,set1,:)))';
%                 fxA = FTdata.nscoh.freq';
                fxA = FTdata.wpli.freq;
%                 fxA = FTdata.freqPow.freq; % Frequency X vector
                res = fxA(2) - fxA(1);  % Freq resolution
            end
            if ~isempty(inds)
                frqxsave{srcloc,sub,cond} = mean(fy,1);
%                 frqxsave{srcloc,sub,cond} = fy;
            end
        end
    end
end

%        label: {20x1 cell}
%        dimord: 'chan_freq'
%          freq: [1x145 double]
%     powspctrm: [20x145 double]
%      labelcmb: {190x2 cell}
%     crsspctrm: [190x145 double]
%           cfg: [1x1 struct]

close all

for srcloc = 2:4
ON = vertcat(frqxsave{srcloc,:,1});
OFF = vertcat(frqxsave{srcloc,:,2});
titular = 'wpli'; ylab = 'wpli'; alpha = 0.05;
[specstat] = spectralclusterstats060317(ON,OFF,fxA,R.sourcenames{srcloc})

figure(srcloc*10); clf
freqclusterplot(OFF,ON,fxA,specstat,alpha,titular,ylab)

clear p
for i = 1:size(fxA,2)
[h,p(i)] = ttest2(ON(:,i),OFF(:,i));
end
p(p>0.05) = 1;


% figure(srcloc)
% plot(fxA,p) 
% hold on
% plot(fxA,repmat(0.05,1,size(fxA,2)))
% xlim([4 70])
end