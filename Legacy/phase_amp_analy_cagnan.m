% function plot_phase_amp_analysis_PLIs_V3_RAT(R)
% if nargin<1
%     R = makeHeader_SubCort_Cort_Networks();
% end
clear
R = buildheader_rat;

% Collate the points and then do paired ttest between bins? A lot of
% information is lost by cutting off in histograms. Either that or somehow
% increase the number of samples - perhaps include all of the STN channels?
% Maybe combine left and right into one histogram? Collapse all of one
% subjects sync segs into one histogram?? All should be sampling from the
% same distribution so build up better resolved image of underlying 2D PDF.
% Or look at doing a hotelling test between the collated data??

close all
analynames = {'Segment Length','CTX High Beta Amp','STN High Beta Amp','STN Low Beta Amp','STN/CTX High Beta Amp Correlation','Causal Density'};

QX = 8 ; QY = 6; % 6 6
% QX = [-pi -pi/8 pi/8 pi]';
% QX = linspace(-pi,pi,QX);
% load([R.datapathr 'subject_hbWPLI075'])
% subscreen = squeeze(sum(subject_hbcohscreen>R.PA.WPLIscreen)==2);
a = logspace(0.5,log10(150),4); logscalez = [-100 -50 -15 15 50 100]; %-40:10:40; %linspace(-50,50,QY); %logscalez = [-3.^(4:-1:2) 3.^(2:4)];
for band = [2 3] %[1 3]
    for cond =1:2
        for sub  = 1:length(R.subnames{cond})
            load([R.analysispath R.pipestamp '\data\processed\' R.subnames{cond}{sub} '_' R.condnames{cond} '_' R.pipestamp '_' R.bandnames{band} '.mat'],'FTdata')
            pA_pli_dist_save = FTdata.PA.pA_pli_dist_save;
            segL_pli_dist_save = FTdata.PA.segL_pli_dist_save;
            amp_pli_dist_save = FTdata.PA.amp_pli_dist_save;
            H_dist_save = FTdata.PA.H_dist_save;
            timevec = FTdata.PA.timevec{1};
            
            relativePhi{1} = pA_pli_dist_save{1}; %wrapToPi(pA_pli_dist_save{1}-circ_mean(pA_pli_dist_save{1}') ); %circ_mean(pA_pli_dist_save{2,nrOFF}(amp_pli_dist_save{2,nrOFF}(3,:)>95)')); %wrapToPi(pA_pli_dist_save{cond,nr}); %-circ_mean([pA_pli_dist_save{2,nrOFF}]'));
            segL{1} =  segL_pli_dist_save{1}; %((segL_pli_dist_save{cond,nr} - mean([segL_pli_dist_save{cond,nr}],2))./mean([segL_pli_dist_save{cond,nr}],2)  )*100;
            ampSeg{1} = amp_pli_dist_save{1}; %((amp_pli_dist_save{cond,nr}  - mean([amp_pli_dist_save{cond,nr}],2))./mean([amp_pli_dist_save{cond,nr}],2) )*100;
            HdistSeg{1} = H_dist_save{1}(1,:);
            %             tends(nr) = timevec{1}(end)- timevec{1}(1);
            tendtot = timevec;
            relativePhiCol = [relativePhi{:,:}];
            segLCol = [segL{:,:}];
            ampSegCol = [ampSeg{:,:}];
            HdistSegCol = [HdistSeg{:,:}];
            
            phiBin = linspace(-pi,pi,QX);
            phiBinMid = phiBin(1:end-1)+((phiBin(2)-phiBin(1))/2);
            ampBin = [];
            for i = 1:length(phiBin)-1
                binDat = ampSegCol(:,relativePhiCol>=phiBin(i) & relativePhiCol<=phiBin(i+1))';
                ampBinMu(:,i) = nanmean(binDat);
                ampBinSEM(:,i) = nanstd(binDat)/size(binDat,1);
                
                binSeg = segLCol(:,relativePhiCol>=phiBin(i) & relativePhiCol<=phiBin(i+1))';
                segBinMu(:,i) = nanmean(binSeg);
                segBinSEM(:,i) =nanstd(binSeg)/size(binSeg,1);
            end
            
            
            ampBinGroup{cond,sub} = ampBinMu;
            segBinGroup{cond,sub} = segBinMu;
            figure(1+(10*band))
            cmap = linspecer(3);
            panlist = [1 3 5 ; 2 4 6];
            obs = {['M2 ' R.bandnames{band}],['STN ' R.bandnames{band}],['STN ' R.bandnames{band-1}]};
            for i = 1:3
                subplot(3,2,panlist(cond,i))
                hl = scatter(phiBinMid', ampBinMu(i,:)','MarkerEdgeColor',cmap(i,:)); hold on
                a =plot(phiBinMid', ampBinMu(i,:)','color',cmap(i,:));
                if cond == 1; a.LineStyle = '--'; end
                %                 [hl, hp] = boundedline(phiBinMid', ampBinMu(i,:)',ampBinSEM(i,:)','cmap',cmap(i,:));
                %                 hp.FaceAlpha = 1;
                if cond == 1; hl.Marker = '*'; end
                ylim([-75 115])
                title(['STN-M2 Phase vs ' obs{i}])
                ylabel(['Amplification of ' obs{i}]); xlabel('Relative Phase')
                grid on; hold on
                
                [xq yq R2 flag] = sinfit(phiBinMid,ampBinMu(i,:),50,[25; 2*pi; 0; -500],[300;2*pi; 2*pi; 500],0);
                if flag
                    a =plot(xq,yq,'color',cmap(i,:),'MarkerSize',0.05)
                    if cond == 1; a.LineStyle = '--'; end
                end
            end
            
            figure(2+(10*band))
            cmap = linspecer(5);
            subplot(1,2,cond)
            hl = scatter(phiBinMid', segBinMu(1,:)','MarkerEdgeColor',cmap(5,:)); hold on
            a =plot(phiBinMid', segBinMu(1,:)','color',cmap(5,:));
            if cond == 1; a.LineStyle = '--'; end
            %             [hl, hp] = boundedline(phiBinMid', segBinMu(1,:)',segBinSEM(1,:)','cmap',cmap(5,:));
            hp.FaceAlpha = 0.4;
            if cond == 1; hl.Marker = '*'; end
            ylim([0 1.1]); grid on
            title(['STN-M2 Phase vs ' R.bandnames{band} ' Frame Length'])
            ylabel(['LB Segment Length (s)']); xlabel('Relative Phase')
            
            [xq yq R2 flag] = sinfit(phiBinMid,segBinMu(1,:),50,[0.05; 2*pi; 0; 0],[0.3; 2*pi; 2*pi; 5],0);
            if flag
                a = plot(xq,yq,'color',cmap(5,:),'MarkerSize',0.05);
                if cond == 1; a.LineStyle = '--'; end
                
            end
            
        end
        %         save([R.analysispath R.pipestamp '\results\PAdata\PAgroup'],'segLsave','amp_pli_dist_save','densesave','hdistsave','phaseAng_dist_save'); %,'gc_dist_sub_save','GC_stat_table')
    end
    figure(1+(10*band))
    set(gcf,'Position',[680          94         864        1004])
    for cond=1:2
        y =  horzcat(ampBinGroup{cond,:});
        for i = 1:3
            subplot(3,2,panlist(cond,i)); hold on
            a = plot(phiBinMid,nanmean(reshape(y(i,:),size(phiBinMid,2),[])'),'color',cmap(i,:),'LineWidth',2);
            if cond == 1; a.LineStyle = '--'; end
        end
    end
    annotation(gcf,'textbox',...
        [0.187 0.956 0.205 0.035],...
        'String',{'Control'},...
        'LineStyle','none',...
        'HorizontalAlignment','center',...
        'FontSize',18,...
        'FitBoxToText','off');
    
    annotation(gcf,'textbox',...
        [0.624 0.959 0.205 0.0349],...
        'String','Lesion',...
        'LineStyle','none',...
        'HorizontalAlignment','center',...
        'FontSize',18,...
        'FitBoxToText','off');
    
    figure(2+(10*band))
    set(gcf,'Position',[681   761   862   336])
    cmap = linspecer(5);
    for cond=1:2
        y =  vertcat(segBinGroup{cond,:});
        subplot(1,2,cond); hold on
        a = plot(phiBinMid,nanmean(y),'color',cmap(5,:),'LineWidth',2);
        if cond == 1; a.LineStyle = '--'; end
    end
    for bandr = 1:3
        clear ON OFF
        for cond =1:2
            for sub  = 1:length(R.subnames{cond})
                if cond ==1
                    ON(:,sub) = ampBinGroup{cond,sub}(bandr,:);
                else
                    OFF(:,sub) = ampBinGroup{cond,sub}(bandr,:);
                end
            end
        end
        %         ON(isnan(ON)) = 0; OFF(isnan(OFF)) = 0;
        [h pvar] = vartest2(ON',OFF')
        [h ptt] = ttest2(ON',OFF')
        figure(1+(10*band))
        for cond=1:2
            subplot(3,2,panlist(cond,bandr)); hold on
            x = phiBinMid(pvar<0.05);
            scatter(x,repmat(80,size(x)),50,'d','filled',...
                'MarkerFaceColor',cmap(bandr,:),'MarkerEdgeColor',cmap(bandr,:),'LineWidth',2)
            x = phiBinMid(ptt<0.05);
            scatter(x,repmat(105,size(x)),65,'s','filled',...
                'MarkerFaceColor',cmap(bandr,:),'MarkerEdgeColor',cmap(bandr,:),'LineWidth',2)
            
            %             if cond==1
            %                 ONphi =repmat(phiBinMid,size(ON,2),1);
            %                 ONphi =ONphi(~isnan(ON)); ON = ON(~isnan(ON));
            %                 [xq yq R2 flag] = sinfit(ONphi(:),ON(:),50,[25; 2*pi; 0; -500],[300;2*pi; 2*pi; 500],0);
            %                 if flag == 1
            %                     plot(xq,yq,'-','color',cmap(bandr,:),'MarkerSize',0.05,'LineWidth',3)
            %                 end
            %             else
            %                 OFFphi =repmat(phiBinMid,size(OFF,2),1);
            %                 OFFphi =OFFphi(~isnan(OFF)); OFF = OFF(~isnan(OFF));
            %                 [xq yq R2 flag] = sinfit(OFFphi(:),OFF(:),50,[25; 2*pi; 0; -500],[300;2*pi; 2*pi; 500],0);
            %                 if flag == 1
            %                     plot(xq,yq,'-','color',cmap(bandr,:),'MarkerSize',0.05,'LineWidth',3)
            %                 end
            %             end
            
            
        end
    end % End of stars for amp
    
    %% Now for segL Group Stats
    clear ON OFF
    for cond =1:2
        for sub  = 1:length(R.subnames{cond})
            if cond ==1
                ON(:,sub) = segBinGroup{cond,sub};
            else
                OFF(:,sub) = segBinGroup{cond,sub};
            end
            
        end
    end
    
    figure(2+(10*band))
    [h pvar] = vartest2(ON',OFF')
    [h ptt] = ttest2(ON',OFF')
    for cond=1:2
        subplot(1,2,cond); hold on
        x = phiBinMid(pvar<0.05);
        scatter(x,repmat(0.9,size(x)),50,'d','filled',...
            'MarkerFaceColor',cmap(5,:),'MarkerEdgeColor',cmap(5,:),'LineWidth',2)
        x = phiBinMid(ptt<0.05);
        scatter(x,repmat(0.97,size(x)),65,'s','filled',...
            'MarkerFaceColor',cmap(5,:),'MarkerEdgeColor',cmap(5,:),'LineWidth',2)
        %         if cond==1
        %             ONphi =repmat(phiBinMid,size(ON,2),1);
        %             ONphi =ONphi(~isnan(ON)); ON = ON(~isnan(ON));
        %             [xq yq R2 flag] = sinfit(ONphi(:),ON(:),50,[0.05; 2*pi; 0; 0],[0.3; 2*pi; 2*pi; 5],0);
        %             if flag == 1
        %                 plot(xq,yq,'-','color',cmap(5,:),'MarkerSize',0.05,'LineWidth',3)
        %             end
        %         else
        %             OFFphi =repmat(phiBinMid,size(OFF,2),1);
        %             OFFphi =OFFphi(~isnan(OFF)); OFF = OFF(~isnan(OFF));
        %             [xq yq R2 flag] = sinfit(OFFphi(:),OFF(:),50,[0.05; 2*pi; 0; 0],[0.3; 2*pi; 2*pi; 5],0);
        %             if flag == 1
        %                 plot(xq,yq,'-','color',cmap(5,:),'MarkerSize',0.05,'LineWidth',3)
        %             end
        %         end
    end % End of Frame Length Group stat/plot
    
end






