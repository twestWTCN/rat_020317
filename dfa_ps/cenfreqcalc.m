function cohfreq = cenfreqcalc(FTdata,method,pairs,R)
switch method
    case 'powcen'
        for i = 1:size(pairs,1)
            cohspec = FTdata.freqPow.powspctrm(pairs(i),:);
            frq = FTdata.freqPow.freq;
            for band = 1:3
                [dum k1] = find(frq>=R.bbounds(band,1) & frq<R.bbounds(band,2));
                frqb = frq(k1);
                [dum k2] = max(cohspec(k1)); cohfreq(i,band) = frqb(k2);
            end
        end
        
        
    case 'fboundcen' % centre of predefined bands
        cohfreq = repmat(mean(R.bbounds,2)',size(pairs,1),1);    % use centre of band
    case 'wpli'    % Peak of WPLI
        for i = 1:size(pairs,1)
            cohspec = squeeze(FTdata.wpli.wpli_debiasedspctrm(pairs(1,1),pairs(1,2),:));
            frq = FTdata.wpli.freq;
            for band = 1:3
                [dum k1] = find(frq>=R.bbounds(band,1) & frq<R.bbounds(band,2));
                frqb = frq(k1);
                [dum k2] = max(cohspec(k1)); cohfreq(i,band) = frqb(k2);
            end
        end
    case 'coh'  % Peak of Coh
        for i = 1:size(pairs,1)
            cohspec = squeeze(FTdata.coherence.cohspctrm(pairs(1,1),pairs(1,2),:));
            frq = FTdata.coherence.freq;
            for band = 1:3
                [dum k1] = find(frq>=R.bbounds(band,1) & frq<R.bbounds(band,2));
                frqb = frq(k1);
                [dum k2] = max(cohspec(k1)); cohfreq(i,band) = frqb(k2);
            end
        end
end
