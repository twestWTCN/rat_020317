function plot_phase_amp_analysis_PLIs_SegStatsPlot_RATV8(R)
if nargin<1
    R = buildheader_rat();
end
R.condname = R.condnames;
close all
% WIP
% Percentage Change over what? Baseline or baseline of segments?
QX = 8 ; % Bin Size
for band = 3 %[2 3]
    cmapint = linspecer(3);
    cmap = linspecer(5);
    cmapint(4,:) = cmap(5,:);
    for cond =1:2
        for sub  = 2%:length(R.subnames{cond})
            load([R.analysispath R.pipestamp '\data\processed\' R.subnames{cond}{sub} '_' R.condnames{cond} '_' R.pipestamp '.mat'])
            stnlabs = find(strncmp('STN',FTdata.wpli.label,3));
            for stchi = 1:length(stnlabs)
                phi = FTdata.PA(band).pA_pli_dist_save{stchi}';
                relativePhiCol = phi'; %wrapToPi(phi-circ_mean(phi(~isnan(phi'))) )'; %circ_mean(pA_pli_dist_save{2,nrOFF}(amp_pli_dist_save{2,nrOFF}(3,:)>95)')); %wrapToPi(pA_pli_dist_save{cond,nr}); %-circ_mean([pA_pli_dist_save{2,nrOFF}]'));
                segLCol =  FTdata.PA(band).segL_pli_dist_save{stchi}; %((segL_pli_dist_save{cond,nr} - mean([segL_pli_dist_save{cond,nr}],2))./mean([segL_pli_dist_save{cond,nr}],2)  )*100;
                ampSegCol = FTdata.PA(band).amp_pli_dist_save{stchi}; %((amp_pli_dist_save{cond,nr}  - mean([amp_pli_dist_save{cond,nr}],2))./mean([amp_pli_dist_save{cond,nr}],2) )*100;
                %                 HdistSegCol = FTdata.PA(band).H_dist_save(1,:);
                tendtot = FTdata.PA(band).timevec{1}(end)-FTdata.PA(band).timevec{1}(1);
                
                phiBin = linspace(-pi,pi,QX);
                for i = 1:3
                    [shiftPhiCol phipeak(i) binind]= findAmpPhi(R,ampSegCol(i,:),relativePhiCol,phiBin);
                    selind{i} = find(phi>=phiBin(binind) & phi<=phiBin(binind+1));
                end
                [shiftPhiCol phipeak(4) bind]= findAmpPhi(R,segLCol,relativePhiCol,phiBin);
                 selind{4} = find(phi>=phiBin(binind) & phi<=phiBin(binind+1));
                
                % Length/Amp Correlations
                for i = 1:3
                    x = segLCol(selind{i})'; y = ampSegCol(i,selind{i})';
                    [x,y] = remnan(x,y);
                    % %                         nleg = 50; nseg  = fix(length(x)/nleg);
                    % %                         x = reshape(x(1:nleg*nseg),nseg,nleg);
                    % %                         y = reshape(y(1:nleg*nseg),nseg,nleg);
                    % %                         for ip =1:nseg
                    % %                             [r(ip) p(ip)] = corr(log10(x(ip,:))',y(ip,:)','Type','Spearman');
                    % %                         end
                    % %                         rGroup(i,sub,side,cond) = mean(r);
                    % %                         pGroup(i,sub,side,cond) = mean(p);
                    try
                        [r p] = corr(log10(x),y,'Type','Spearman');
                    catch
                        r = 0; p = 1;
                        disp('Correlation Not Computed!!!')
                    end
                    rGroup{i,1,sub,cond} = r;
                    pGroup{i,1,sub,cond} = p;
                end
                
                % Relative Phase
                dp = pi/12;
                phibins = -pi:dp:pi;
                phibinmid = phibins(1:end-1)+(dp/2);
                c = histcounts(relativePhiCol,phibins);
                phiGroup{1,sub,cond} = c./tendtot;
                % Segment Lengths
                db = 0.1;
                segbins = -1.5:db:1;
                segbinmid = segbins(1:end-1)+(db/2);
                c = histcounts(log10(segLCol),segbins); %,'Normalization','probability');
                segLGroup{1,sub,cond} = c./tendtot;
                %
            end
        end
    end
    
    figure(1+(10*band))
    PhiPlot(phiGroup,phibinmid,cmap,R);
    figure(2+(10*band))
    segLPlot(segLGroup,segbinmid,cmap,R);
    figure(3+(10*band))
    obs = {[R.sourcenames{1} ' ' R.bandinits{band}],['STN ' R.bandinits{band}],['STN ' R.bandinits{band-1}],[R.bandinits{band} ' Seg. Length']};
    for i = 1:3
        r1 = squeeze(vertcat(rGroup{i,:,:,1}));
        r2 = squeeze(vertcat(rGroup{i,:,:,2}));
        p1 = squeeze(vertcat(pGroup{i,:,:,1})); p1 = p1(:);
        p2 = squeeze(vertcat(pGroup{i,:,:,2})); p2 = p2(:);
        s1 = sum(p1<0.05 & p1>0); t1 = sum(p1>0);
        r1s =r1(p1<0.05 & p1>0)
        mean(r1s)
        s2 = sum(p2<0.05 & p2>0); t2 = sum(p2>0); r2s =r2(p1<0.05 & p1>0)
        
        [h,p, chi2stat,df] = prop_test([s1 s2], [t1 t2],'false')
        
        
        subplot(1,3,i)
        boxploter(r1,r2,cmap(i,:),R,obs{i})
    end
    clear segLGroup phiGroup rGroup
    close all
end
%% FUNCTION ENDS HERE

