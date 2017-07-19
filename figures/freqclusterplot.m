function [ax clustat] = freqclusterplot(OFF,ON,fxA,stat,alpha,titular,ylab,ylimd,cmap,dir)
if nargin<5
    alpha = 0.05;
end
if nargin<6
    titular = 'spectral stats for unknown analysis';
end
if nargin<7
    ylab = 'scaler';
end
if nargin <10
    dir = 0;
end
if size(fxA,1)>size(fxA,2)
    fxA = fxA';
end

lineind = find(fxA>48.5 & fxA<51.5);
OFF(:,lineind) = NaN(size(OFF,1),size(lineind,2));
ON(:,lineind) = NaN(size(ON,1),size(lineind,2));

if dir == 1
    plot(repmat(fxA,size(OFF,1),1)',OFF','Color',cmap(1,:),'LineStyle','--','linewidth',1);
    hold on
    ax(2) = plot(fxA,nanmean(OFF,1),'color',cmap(1,:),'linewidth',3);
    
    plot(repmat(fxA,size(ON,1),1)',ON','Color',cmap(2,:),'LineStyle','--','linewidth',1);
    hold on
    ax(1) = plot(fxA,nanmean(ON,1),'color',cmap(1,:),'linewidth',3);
else
    plot(repmat(fxA,size(OFF,1),1)',OFF','Color',cmap(1,:),'LineStyle','--','linewidth',1);
    hold on
    ax(2) = plot(fxA,nanmean(OFF,1),'color',cmap(1,:),'linewidth',3);
    
    plot(repmat(fxA,size(ON,1),1)',ON','Color',cmap(2,:),'LineStyle','--','linewidth',1);
    hold on
    ax(1) = plot(fxA,nanmean(ON,1),'color',cmap(2,:),'linewidth',3);
end
% emsure mask at 49-51 is zero
stat.mask(stat.freq> 49 & stat.freq < 51) = 0;
sigfreq = stat.mask.*stat.freq;
sigfreq(sigfreq==0) = NaN;
sigpow = stat.mask.*(nanmean(nanmean(ON,1),2).*10);
if ~isempty(ylimd)
    sigpow = stat.mask.*ylimd(2)*0.85;
end
% sigpow = stat.mask.*1; %(nanmean(nanmean(ON,1),2).*7.5);

sigpow(sigpow==0) = NaN;
% if sum(~isnan(sigpow))<2
%     scatter(sigfreq,sigpow,'ks','filled');
% else
%     plot(sigfreq,sigpow,'k','linewidth',6);
% end
if ~isempty(ylimd)
    ylim(ylimd);
end
xlim([4 100]);
clustat = [];
if sum(stat.mask)>0
    for i = 1:size(stat.posclusters,2)
        if stat.posclusters(i).prob<alpha
            labs = stat.posclusterslabelmat;
            group = find(labs==i);
            if length(group)<2
                scatter(sigfreq(group),sigpow(group),'ks','filled');
            else
                plot(sigfreq(group),sigpow(group),'k','linewidth',6);
            end
            freqcen = mean(group);
            [figx figy] = dsxy2figxy(gca, freqcen-5, ylimd(2)*0.9);
            h = annotation('textbox',[figx figy .01 .01],'String',{['P = ' num2str(stat.posclusters(i).prob,'%.4f')]},'FitBoxToText','on','LineStyle','none','fontsize',8);
            clustat = [clustat; min(stat.freq(1,group)) max(stat.freq(1,group)) stat.posclusters(i).clusterstat stat.posclusters(i).prob];
        end
    end
    for i = 1:size(stat.negclusters,2)
        if stat.negclusters(i).prob<alpha
            labs = stat.negclusterslabelmat;
            group = find(labs==i);
            if length(group)<2
                scatter(sigfreq(group),sigpow(group),'ks','filled');
            else
                plot(sigfreq(group),sigpow(group),'k','linewidth',6);
            end
            freqcen = mean(group);
            [figx figy] = dsxy2figxy(gca, freqcen-5, ylimd(2)*0.9);
            h = annotation('textbox',[figx figy .01 .01],'String',{['P = ' num2str(stat.negclusters(i).prob,'%.3f')]},'FitBoxToText','on','LineStyle','none','fontsize',8);
            clustat = [clustat; min(stat.freq(1,group)) max(stat.freq(1,group)) stat.negclusters(i).clusterstat stat.negclusters(i).prob];
        end
    end
end
xlabel('Frequency (Hz)','fontsize',16); ylabel(ylab); title(titular);
grid on