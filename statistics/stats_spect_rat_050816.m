function stats_spect_rat_050816(R)
%% Coherence Statistics- Integrated coherence, peak magnitude, WPLI
if R.clear.specstat == 1
    FTdata.dirstats.pow = [];
end
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
                fy = squeeze(FTdata.spectanaly.spectra{1}(set1,:)); % Combined Coherence spectra
                fxA = FTdata.freqPow.freq; % Frequency X vector
                res = fxA(2) - fxA(1);  % Freq resolution
                for band = 1:size(bbounds,1)
                    intCoh(band,i) = abs(res*trapz(fy(round((bbounds(band,1)-R.FOI(1))/res)+1:round((bbounds(band,2)-R.FOI(1))/res)+1)));
                    [peakMag(band,i) frq] = max(fy(round((bbounds(band,1)-R.FOI(1))/res)+1:round((bbounds(band,2)-R.FOI(1))/res)+1));
                    peakFrq(band,i) =  ((frq-1)*res)+bbounds(band,1);
                end
            end
            FTdata.dirstats.pow.(R.sourcenames{srcloc}).intCoh = intCoh;    % Integrated Coherence
            FTdata.dirstats.pow.(R.sourcenames{srcloc}).peakMag = peakMag;  % peak magnitude within band
            FTdata.dirstats.pow.(R.sourcenames{srcloc}).peakFrq = peakFrq;  % peak frequency (Hz)
            FTdata.dirstats.pow.(R.sourcenames{srcloc}).labels = {FTdata.freqPow.label{inds}};
            clear intCoh peakMag peakFrq
        end
        save([R.analysispath R.pipestamp '\data\processed\' R.subnames{cond}{sub} '_' R.condnames{cond} '_' R.pipestamp '.mat'],'FTdata')
    end
end