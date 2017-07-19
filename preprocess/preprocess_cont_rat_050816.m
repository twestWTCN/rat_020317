function cleandata = preprocess_cont_rat_050816(cfg,FTdata)
data = FTdata.trial{1};
% hpFilt = designfilt('highpassfir','StopbandFrequency',R.pp.hpfilt(1), ...
%     'PassbandFrequency',R.pp.hpfilt(2),'PassbandRipple',0.05, ...
%     'StopbandAttenuation',50,'DesignMethod','kaiserwin','SampleRate',FTdata.fsample);

for D = 1:size(data,1)
    X = data(D,:); % assign data
    % Truncate
    X = X(1,(cfg.start*FTdata.fsample):end-(cfg.end*FTdata.fsample)); %
    
    % Interpolate NaNs
    nanx = isnan(X);
    t    = 1:numel(X);
    X(nanx) = interp1(t(~nanx), X(~nanx), t(nanx));
    
    % De-mean
    X = X - mean(X);
    
%     % Filter
%     X = filtfilt(hpFilt,X);
    
    % Remove Suprathreshold Events
    Xdiff = abs(diff(X));
    xThresh = find(Xdiff>(median(Xdiff)+prctile(Xdiff,95)));
    for i = 1:length(xThresh)
        qI = [xThresh(i)-2,xThresh(i)-1];
        if numel(find(qI)) == 2
            rep = interp1(1:2,X(qI),linspace(1,2,3));
            X(xThresh(i)) = rep(2);
        end
    end
    data_clean(D,:) = X;
end
tvec = FTdata.time{1}(1,(cfg.start*FTdata.fsample):end-(cfg.end*FTdata.fsample));
cleandata.label = FTdata.label;
cleandata.fsample = FTdata.fsample;
cleandata.trial{1} = data_clean;
cleandata.time{1} = tvec;
