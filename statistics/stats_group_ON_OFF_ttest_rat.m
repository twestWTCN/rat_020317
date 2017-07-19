function statOut = stats_group_ON_OFF_ttest_rat(R,analysis,feat,yl,boxy,outlierset,colspec)
featlist = {'alpha','pebscore','varfeat','ARCoeff''banintCoh','peakMag','intCoh','peakFrq','minorm','frwd_intCoh','rev_intCoh'};
featname = {'exponents from','model evidences for','mean variance for signals from','AR coeff for signals from','band ','peak magnitudes for','integrated','peak frequency of','mutual info','Forward Integrated Coh.','Reverse Integrated Coh.'};
analylist = {'dfaps','coh','wpli','pow','MI','MIph','dfaae','npd'};
analyname = {'DFA-PS','Coherence','WPLI','Power','Mutual Information','Mutual Phase Information','DFA-AE','Non-parametric directionality'};
if isempty(yl)
    yl = 'auto';
end
for srcloc = 1:length(R.sourcenames)
    if srcloc == 1
        statmeth = 1;
    else
        statmeth = 2;
    end
    statcell = [];
    for cond = 1:2
        for sub  = 1:length(R.subnames{cond})
            load([R.analysispath R.pipestamp '\data\processed\' R.subnames{cond}{sub} '_' R.condnames{cond} '_' R.pipestamp '.mat'])
            %             if (strncmp(analysis,'dfaps',5)|strncmp(analysis,'dfaae',5)) & sub == 7 %& srcloc == 1
            if sub == 7
                for i = 1:length(feat)
                    statcell{i,cond,sub} = repmat(NaN,3,1);
                end
            else
                if strmatch('pow.source', analysis)
                    istat = FTdata.dirstats.pow.source;
                elseif strmatch('pow.STN', analysis)
                    istat = FTdata.dirstats.pow.STN;
                else
                    istat = FTdata.dirstats.(analysis);
                end
                %         statcell{1,cond,sub} = mean(istat.intCoh,2);
                %         statcell{2,cond,sub} = mean(istat.peakMag,2);
                %         statcell{3,cond,sub} = mean(istat.peakFrq,2);
                for i = 1:length(feat)
                    try
                        insrt = istat.(R.sourcenames{srcloc}).(feat{i});
                    catch
                        insrt = repmat(NaN,3,1);
                    end
                    if size(insrt,1) ~= 3
                        insrt = repmat(NaN,3,1);
                    end
                    if size(insrt) == [1 3]
                        insrt = [];
                    end
                    statcell{i,cond,sub} = insrt;
                end
            end
        end
    end
    
    for i = 1:length(feat)
        boxgroups = []; flag = 0;
        clear p ci stat sumstat
        for band = 1:size(R.bbounds,1)
            lesion = [statcell{i,2,:}]; control = [statcell{i,1,:}];
            X = lesion(band,:); Y = control(band,:);
            N1 = numel(X); N2 = numel(Y);
            if outlierset == 1
                X = rmvoutliers(X,'STD','NaN',1,2.5);
                Y = rmvoutliers(Y,'STD','NaN',1,2.5);
            end
            n1 = numel(X); n2 = numel(Y);
            [chi2stat chipval] = chi2testTW(n1,N1,n2,N2);
            
            
            if numel(X)>3 && numel(Y)>3
                if range(X)>0 && range(Y)>0
                    if statmeth == 1
                        if swtest(X) || swtest(Y)
                            [p(band) ~] = ranksum(X,Y); % non-parametric
                            test{band} = 'ranksum';
                            CI = confint_ttest([nanmean(X) nanmean(Y)],[nanstd(X) nanstd(Y)],numel(~isnan(Y)));
                            L = 1;
                            sumstat(:,band) = [median(X) iqr(X) numel(~isnan(X)) median(Y) iqr(Y) numel(~isnan(Y)) median(X)-median(Y) CI p(band) L];
                        else
                            [dummy p(band)] = ttest2(X,Y);  % parametric
                            test{band} = '2sampTtest';
                            L = 0;
                            CI = confint_ttest([nanmean(X) nanmean(Y)],[nanstd(X) nanstd(Y)],size(control,2)-1);
                            sumstat(:,band) = [mean(X) std(X) numel(~isnan(X)) mean(Y) std(Y) numel(~isnan(Y)) mean(X)-mean(Y) CI p(band) L];
                        end
                    else % Use repeated measures anova
                        p(band) = rms_anova_rat_110816(R,statcell,i,srcloc,band);
                        L = 1;
                        CI = confint_ttest([nanmean(X) nanmean(Y)],[nanstd(X) nanstd(Y)],size(control,2)-1);
                        sumstat(:,band) = [mean(X) std(X) numel(~isnan(X)) mean(Y) std(Y) numel(~isnan(Y)) mean(X)-mean(Y) CI p(band) L];
                    end
                else
                    warning('There is no variance in one or more of your feature vectors')
                    flag = 1;
                    sumstat = NaN(10,size(R.bbounds,1));
                end
                boxgroups = [boxgroups band];
            else
                sumstat(:,band) = repmat(NaN,10,1);
            end
            if strncmp(analysis,'dfaae',5)|| strncmp(analysis,'dfaps',5) % if DFA then use CHI squared count stat
                sumstat(end,band) = chipval;
            end
            
        end
        % Now for plotting
        if boxy == 1 && numel(boxgroups) ~= 0 && flag == 0
            figure
            adata = {control(boxgroups,:)' lesion(boxgroups,:)'};
            [centre height] = aboxplotTW(adata,'labels',{R.bandnames{boxgroups}},'Colorgrad',colspec);
            %             legend({'control','lesion'},'Location','SE')
            title({['Comparison of mean ' featname{find(strcmp(featlist,feat{i}))}];[ analyname{find(strcmp(analylist,analysis))} ' in ' R.sourcenames{srcloc} ' signals']})
            %             ylim([0 .01])
            if numel(boxgroups) == 3
                groups={[centre(1)-.1 centre(1)+.1],[centre(2)-.1
                    centre(2)+.1],[centre(3)-.1 centre(3)+.1]}; % IF BBOUNDS IS 3
            elseif numel(boxgroups) == 2
                groups={[centre(1)-.1 centre(1)+.1],[centre(2)-.1 centre(2)+.1]};
            elseif numel(boxgroups) == 1
                groups={[centre(1)-.1 centre(1)+.1]};
            end
            H = sigstarTW(groups,[p(boxgroups)],0,max(height,[],2)+(.1.*(max(height,[],2))) );
            %             ylim(yl);
            if find(strcmp({'wpli','coh','npd'},analysis))
                pos = 'NorthWest';
            elseif find(strcmp({'dfaps'},analysis))
                pos = 'SouthEast';
            else
                pos = 'NorthEast';
            end
            legend(R.condnames,'Location',pos)
        end
        statOut{srcloc,i} = sumstat';
    end
    
end
end
