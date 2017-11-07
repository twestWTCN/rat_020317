function [clustat] = freqcluster_only_200717(OFF,ON,fxA,stat,alpha,ylimd,cmap,adj)
if nargin<5
    alpha = 0.05;
end

if size(fxA,1)>size(fxA,2)
    fxA = fxA';
end

lineind = find(fxA>48.5 & fxA<51.5);
OFF(:,lineind) = NaN(size(OFF,1),size(lineind,2));
ON(:,lineind) = NaN(size(ON,1),size(lineind,2));

% if dir == 1
%     plot(repmat(fxA,size(OFF,1),1)',OFF','Color',cmap(1,:),'LineStyle','--','linewidth',1);
%     hold on
%     ax(2) = plot(fxA,nanmean(OFF,1),'color',cmap(1,:),'linewidth',3);
%
%     plot(repmat(fxA,size(ON,1),1)',ON','Color',cmap(2,:),'LineStyle','--','linewidth',1);
%     hold on
%     ax(1) = plot(fxA,nanmean(ON,1),'color',cmap(1,:),'linewidth',3);
% else
%     plot(repmat(fxA,size(OFF,1),1)',OFF','Color',cmap(1,:),'LineStyle','--','linewidth',1);
%     hold on
%     ax(2) = plot(fxA,nanmean(OFF,1),'color',cmap(1,:),'linewidth',3);
%
%     plot(repmat(fxA,size(ON,1),1)',ON','Color',cmap(2,:),'LineStyle','--','linewidth',1);
%     hold on
%     ax(1) = plot(fxA,nanmean(ON,1),'color',cmap(2,:),'linewidth',3);
% end
stat.mask(stat.freq> 48 & stat.freq < 52) = 0;
sigfreq = stat.mask.*stat.freq;
sigfreq(sigfreq==0) = NaN;
sigpow = (stat.mask.*nanmean(nanmean(ON,1),2).*7.5)+(adj);
if ~isempty(ylimd)
    sigpow = (stat.mask.*ylimd(2)*0.75)+(adj);
end
% sigpow = stat.mask.*1; %(nanmean(nanmean(ON,1),2).*7.5);

sigpow(sigpow==0) = NaN;
if ~isempty(ylimd)
    ylim(ylimd);
end
xlim([4 75]);
clustat = [];
if sum(stat.mask)>0
    for i = 1:size(stat.posclusters,2)
        if stat.posclusters(i).prob<(alpha/2)
            labs = stat.posclusterslabelmat;
            group = find(labs==i);
            if length(group)<2
                scatter(sigfreq(group),sigpow(group),'ks','filled','MarkerFaceColor',cmap);
            else
                plot(sigfreq(group),sigpow(group),'k','linewidth',6,'color',cmap);
            end
            freqcen = mean(stat.freq(1,group));
            if min(sigfreq(group))<6
                shift = 3;
            else
                shift = 5;
            end
            [figx figy] = dsxy2figxy(gca, freqcen-shift, (max(sigpow)*1.052));
            h = annotation('textbox',[figx figy 0.01 .01],'String',{['(+) P = ' num2str(stat.posclusters(i).prob,'%.3f')]},'FitBoxToText','on','LineStyle','none','fontsize',9.5,'fontweight','bold');
            clustat = [clustat; min(stat.freq(1,group)) max(stat.freq(1,group)) stat.posclusters(i).prob];
        end
    end
    for i = 1:size(stat.negclusters,2)
        if stat.negclusters(i).prob<(alpha/2)
            labs = stat.negclusterslabelmat;
            group = find(labs==i);
            if length(group)<2
                scatter(sigfreq(group),sigpow(group),'ks','filled','MarkerFaceColor',cmap);
            else
                plot(sigfreq(group),sigpow(group),'k','linewidth',6,'color',cmap);
            end
            freqcen = mean(stat.freq(1,group));
            if min(sigfreq(group))<6
                shift = 3;
            else
                shift = 5;
            end
            [figx figy] = dsxy2figxy(gca, freqcen-shift, (max(sigpow)*1.052));
            h = annotation('textbox',[figx figy .01 .01],'String',{['(-) P = ' num2str(stat.negclusters(i).prob,'%.3f')]},'FitBoxToText','on','LineStyle','none','fontsize',9.5,'fontweight','bold');
            clustat = [clustat; min(stat.freq(1,group)) max(stat.freq(1,group)) stat.negclusters(i).prob];
        end
    end
end
grid on