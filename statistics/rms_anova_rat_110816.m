function  P = rms_anova_rat_110816(R,statcell,i,srcloc,band)
subn = 0; ZZ = [];
minreps = min([min(R.nreps{1}(:,srcloc)) min(R.nreps{2}(:,srcloc))]);
minreps(minreps<3) = 3;
% minreps = [min(minreps
for conder = 1:length(R.condnames)
    for suber = 1:length(R.subnames{conder})
        clear Z
        subn = subn+1;
        
        vec = statcell{i,conder,suber}(band,:); % load vector
        vec(isnan(vec)) = [];   % Remove missing data
        rn = length(vec);   % find length
        if ~isempty(vec)
        if rn == minreps(1)
            vec = vec;
        elseif rn>minreps(1)
            vec = vec(randsample(1:rn,minreps(1)));  % make random selection
        elseif rn<minreps(1) && numel(vec)>0
           def = minreps(1) - rn;
           vec(rn+1:minreps(1)) = normrnd(mean(vec(1:rn)),std(vec(1:rn)),[1 def]);
        end
        vec(isnan(vec(:,1))) = nanmean(vec);
        Z(:,1) = vec;
        Z(:,2) = conder;
        Z(:,3) = 1:length(vec);
        Z(:,4) = repmat(subn,length(vec),1);
        ZZ = [ZZ; Z];
        else
            subn = subn-1;
        end
    end
end

[SSQs, DFs, MSQs, Fs, Ps]=mixed_between_within_anova(ZZ,1);
P = Ps{1};
