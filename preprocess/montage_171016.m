function FTdata = montage_171016(FTdata,R)
% 2 bipolar neighbors, 3 is bipolar global average,
% 4 is bipolar distal global, 5 is bipolar distal local
switch R.montage
    case 'monopolar'
                %% Bipolar Fixed
    case 'bipolar_fixed'
        for i = 2:length(R.sourcenames)
            flabs = find(strncmp(R.sourcenames{i}, FTdata.label,2));
            flabs(flabs==2) = [];
            A = FTdata.trial{1}(flabs,:); B = FTdata.trial{1}(2,:);
            FTdata.trial{1}(flabs,:) = bsxfun(@minus,A,B);
        end
    case 'bipolar'
        trial = []; label = {};
        label = {R.sourcenames{1}}; trial = FTdata.trial{1}(1,:);
        for i = 2:length(R.sourcenames)
            flabs = find(strncmp(R.sourcenames{i}, FTdata.label,size(R.sourcenames{i},2))); % find indexs of source type
            flabsc = nchoosek(flabs,2); % find all permuations of indices
            mont = flabsc(find((flabsc(:,2)-flabsc(:,1))==1),:); % take only pairs 2 step apart
            if size(mont,1)<2
                mont = flabsc(find((flabsc(:,2)-flabsc(:,1))==1),:); % take only pairs 1 step apart
            end
            for j = 1:size(mont,1)
                trial = [trial; normaliseV(FTdata.trial{1}(mont(j,1),:)) -  normaliseV(FTdata.trial{1}(mont(j,2),:))]; % add to montage
                label = [label [R.sourcenames{i} ' ' sprintf('%d_%d',mont(j,1),mont(j,2))]]; % add label
            end
        end
        FTdata.trial{1} = trial;
        FTdata.label = label';
        %% Global Average
    case 'globalaverage'
        for i = 2:length(R.sourcenames)
            flabs = find(strncmp(R.sourcenames{i}, FTdata.label,2));
            A = FTdata.trial{1}(flabs,:); B = mean(FTdata.trial{1}(flabs,:),1);
            FTdata.trial{1}(flabs,:) = bsxfun(@minus,A,B);
        end
        %% Far Global
    case 'globaldistal'
        label = {R.sourcenames{1}}; trial = FTdata.trial{1}(1,:);
        for i = 2:length(R.sourcenames)
            iflabs = find(strncmp(R.sourcenames{i}, FTdata.label,size(R.sourcenames{i},2))); % find indexs of source type
            exflabs = setdiff(2:size(FTdata.label,1),iflabs)';
            flabsc = combvec(iflabs',exflabs')';
            mont = flabsc(find(abs((flabsc(:,2)-flabsc(:,1)))==8),:);
            % Now do the montaging
            for j = 1:size(mont,1)
                trial = [trial; normaliseV(FTdata.trial{1}(mont(j,1),:)) -  normaliseV(FTdata.trial{1}(mont(j,2),:))]; % add to montage
                label = [label [R.sourcenames{i} ' ' sprintf('%d_%d',mont(j,1),mont(j,2))]]; % add label
            end
        end
        FTdata.trial{1} = trial;
        FTdata.label = label';
        %% Far Local
    case 'localdistal'
        label = {R.sourcenames{1}}; trial = FTdata.trial{1}(1,:);
        for i = 2:length(R.sourcenames)
            clear mont
            flabs = find(strncmp(R.sourcenames{i}, FTdata.label,size(R.sourcenames{i},2))); % find indexs of source type
            for k = 1:size(flabs,1)
                abdifs = abs(flabs - flabs(k));
                [num ind] = max(abdifs);
                mont(k,:) = [flabs(k) flabs(ind)];
            end
            mont = unique(sort(mont,2),'rows'); % Use only unique combinations
            for j = 1:size(mont,1)
                trial = [trial; normaliseV(FTdata.trial{1}(mont(j,1),:)) -  normaliseV(FTdata.trial{1}(mont(j,2),:))]; % add to montage
                label = [label [R.sourcenames{i} ' ' sprintf('%d_%d',mont(j,1),mont(j,2))]]; % add label
            end
        end
        FTdata.trial{1} = trial;
        FTdata.label = label';
    case 'fastica'
        for i = 2:length(R.sourcenames)
            flabs = find(strncmp(R.sourcenames{i}, FTdata.label,2));
            A = FTdata.trial{1}(flabs,:); B = fastica(FTdata.trial{1}(2:end,:));
            FTdata.trial{1}(flabs,:) = bsxfun(@minus,A,B(1,:));
        end
%       case 'princomp'
%         for i = 2:length(R.sourcenames)
%             flabs = find(strncmp(R.sourcenames{i}, FTdata.label,2));
%             A = FTdata.trial{1}(flabs,:); B = pca(FTdata.trial{1}(2:end,:));
%             FTdata.trial{1}(flabs,:) = bsxfun(@minus,A,B(1,:));
%         end        
end