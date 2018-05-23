function plot_phase_amp_analysis_PLIs_RATV7_Cagnan(R)
if nargin<1
R = buildheader_rat();
end
close all
analynames = {'Segment Length','CTX High Beta Amp','STN High Beta Amp','STN Low Beta Amp','STN/CTX High Beta Amp Correlation','Causal Density'};

QX = 8 ; % Bin Size
for band = [2 3]
    for cond =1:2
        for sub  = 1:length(R.subnames{cond})
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
                phiBinMid = phiBin(1:end-1)+((phiBin(2)-phiBin(1))/2);
                [sub cond]
                clear ampBinMu ampBinSEM
                for i = 1:3
                    shiftPhiCol = findAmpPhi(R,ampSegCol(i,:),relativePhiCol,phiBin);
                    [ampBinMu(i,:) ampBinSEM(i,:)] = binstats(shiftPhiCol,ampSegCol(i,:),phiBin);
                end
                
                shiftPhiCol = findAmpPhi(R,segLCol,relativePhiCol,phiBin);
                [segBinMu segBinSEM] = binstats(shiftPhiCol,segLCol,phiBin);
                segBinSEM(isnan(segBinSEM)) = 0;
                
                
                ampBinGroup{cond,sub,stchi} = ampBinMu;
                segBinGroup{cond,sub,stchi} = segBinMu;
                figure(1)
                cmap = linspecer(3);
                panlist = [1 3 5 ; 2 4 6]; ylimlist = {{[0 0],[0 0],[0 0]},{[-50 75];[-50 75];[-50 75]},{[-50 75];[-50 75];[-50 75]}};
                obs = {[R.sourcenames{1} ' ' R.bandnames{band}],['STN ' R.bandnames{band}],['STN ' R.bandnames{band-1}]};
                for i = 1:3
                    subplot(4,2,panlist(cond,i))
                    hl = plot(phiBinMid', ampBinMu(i,:)','--','color',cmap(i,:)); hold on
                    % %                         [hl, hp] = boundedline(phiBinMid', ampBinMu(i,:)',ampBinSEM(i,:)','cmap',cmap(i,:)); hold on
                    % %                         if cond == 1; hl.LineStyle = '--'; end
                    % %                         hp.FaceAlpha = 0.4;
                    % %                 [xq yq R2 exitflag s(:,sub,cond,i)] = VMfit(phiBinMid',ampBinMu(i,:)',20,[0,-2*pi,-200,-100],[6*pi,2*pi,200,100],0);
                    % %                 if exitflag == 1
                    % %                     hold on; plot(xq,yq,'color',cmap(i,:));
                    % %                 end
                    % %                         [xq yq R2 exitflag] = sinfit(phiBinMid', ampBinMu(i,:)',20,[25; 2*pi; 0; -500],[300;2*pi; 2*pi; 500],0);
                    % %                         plot(xq,yq,'color',cmap(i,:));
                    ylim(ylimlist{band}{i})
                    title(['STN-' R.sourcenames{1} ' ' R.bandnames{band} ' Phase vs ' obs{i} ' Power'])
                    ylabel(['% Change in ' obs{i}]); xlabel('Relative Phase')
                    grid on
                end
                
                cmap = linspecer(5);
                panlist(:,4) = [7; 8];
                subplot(4,2,   panlist(cond,4))
                hl = plot(phiBinMid', segBinMu(1,:)','--','color',cmap(5,:)); hold on
                % %                     [hl, hp] = boundedline(phiBinMid', segBinMu(1,:)',segBinSEM(1,:)','cmap',cmap(5,:));
                % %                     if cond == 1; hl.LineStyle = '--'; end
                % %                     hp.FaceAlpha = 0.4;
                
                % %             [xq yq R2 exitflag s(:,sub,cond,4)] = VMfit(phiBinMid',segBinMu(1,:)',20,[0,-1,-5,-5],[6*pi,1,5,5],0);
                %             [xq yq R2 exitflag] = sinfit(phiBinMid', segBinMu(1,:)',20,[0.01; 2*pi; 0; -1],[0.5; 2*pi; 2*pi; 1],0);
                %             if exitflag == 1
                %                 hold on; plot(xq,yq,'color',cmap(5,:));
                %             end
                ylimlist = {[0 0],[0 1.5],[0 1.5]};
                ylim(ylimlist{band}); grid on
                title(['STN-' R.sourcenames{1} ' ' R.bandnames{band} ' Phase vs ' R.bandnames{band} ' Frame Length'])
                ylabel(['LB Segment Length (s)']); xlabel('Relative Phase')
            end
        end
    end
    
    figure(1)
    set(gcf,'Position',[680    -3   864   999])
    for cond=1:2
        y =  horzcat(ampBinGroup{cond,:,:});
        for i = 1:3
            subplot(4,2,panlist(cond,i)); hold on
            y1 = reshape(y(i,:),size(phiBinMid,2),[])';
            ysave{cond,i} = y1;
            [hl hp] = boundedline(phiBinMid', nanmean(y1)',nanstd(y1)','cmap',cmap(i,:)); %
            hl.LineWidth = 2;
            L = get(gca,'children'); L(2) = L(end); L(end) = hp;set(gca,'children',L)
            %             a = plot(phiBinMid,nanmean(reshape(y(i,:),size(phiBinMid,2),[])'),'color',cmap(i,:),'LineWidth',2);
            % %             if cond == 1; a.LineStyle = '--'; end
        end
    end
    cmapA = linspecer(3)
    cmap = linspecer(5);
    cmapA(4,:) = cmap(5);
    for cond=1:2
        y =  vertcat(segBinGroup{cond,:,:});
        ysave{cond,4} = y;
        subplot(4,2, panlist(cond,4)); hold on
        [hl hp] = boundedline(phiBinMid', nanmean(y)',nanstd(y),'cmap',cmap(5,:));
        hl.LineWidth = 2;
        L = get(gca,'children'); L(2) = L(end); L(end) = hp; set(gca,'children',L)
        %         a = plot(phiBinMid,nanmean(y),'color',cmap(5,:),'LineWidth',2);
        % %         if cond == 1; a.LineStyle = '--'; end
    end
    
    
    [ptt,pvart] = cagnan_stats_test(phiBinMid,ysave,1,panlist,cmapA);
    
    annotation(gcf,'textbox',...
        [0.187 0.956 0.205 0.035],...
        'String',R.condnames{1},...
        'LineStyle','none',...
        'HorizontalAlignment','center',...
        'FontSize',18,...
        'FitBoxToText','off');
    
    annotation(gcf,'textbox',...
        [0.624 0.959 0.205 0.0349],...
        'String',R.condnames{2},...
        'LineStyle','none',...
        'HorizontalAlignment','center',...
        'FontSize',18,...
        'FitBoxToText','off');
    savefigure_v2([R.analysispath R.pipestamp '\results\'],...
        [R.bandnames{band} '_Group_CagnanAnalysis'],[],[],'-r100'); close all
    close all
    clear ampBinGroup segBinGroup
end




