function [ax] = freqclusterplot_dir(OFF,ON,fxA,stat,alpha,titular,ylab,ylimd,cmap)
if nargin<5
    alpha = 0.05;
end

if nargin<6
    titular = 'spectral stats for unknown analysis';
end
if nargin<6
    ylab = 'scaler';
end

plot(repmat(fxA,size(OFF,1),1)',OFF','Color',cmap(1,:),'LineStyle','-','linewidth',1);
hold on
ax(1) = plot(fxA,nanmean(OFF,1),'color',cmap(1,:),'linewidth',4);

plot(repmat(fxA,size(ON,1),1)',ON','Color',cmap(2,:),'LineStyle','--','linewidth',1);
hold on
ax(2) = plot(fxA,nanmean(ON,1),'color',cmap(1,:),'LineStyle','--','linewidth',4);

sigfreq = stat.mask.*stat.freq;
sigfreq(sigfreq==0) = NaN;
sigpow = stat.mask.*(nanmean(nanmean(ON,1),2).*10);
if ~isempty(ylimd)
sigpow = stat.mask.*ylimd(2)*0.85;
end
% sigpow = stat.mask.*1; %(nanmean(nanmean(ON,1),2).*7.5);

sigpow(sigpow==0) = NaN;
 plot(sigfreq,sigpow,'k','linewidth',6);
if ~isempty(ylimd)
    ylim(ylimd);
end
xlim([4 70]);
if sum(stat.mask)>0
    for i = 1:size(stat.posclusters,2)
        if stat.posclusters(i).prob<alpha
        labs = stat.posclusterslabelmat;
        dx = diff(diff(labs));
                if labs(1) == 1
            dx(1) = -1;
                end
        y = reshape(find(dx==-1),[],2) + 1;
        
            freqcen = mean(stat.freq(1,y(i,1):y(i,2)));
            [figx figy] = dsxy2figxy(gca, freqcen-5, ylimd(2)*0.9);
            h = annotation('textbox',[figx figy .01 .01],'String',{['P = ' num2str(stat.posclusters(i).prob,'%.3f')]},'FitBoxToText','on','LineStyle','none','fontsize',8);
        end
    end
    for i = 1:size(stat.negclusters,2)
        if stat.negclusters(i).prob<alpha
                    labs = stat.negclusterslabelmat;
        dx = diff(diff(labs));
        y = reshape(find(dx==-1),[],2) + 1;

            freqcen = mean(stat.freq(1,y(i,1):y(i,2)));
            [figx figy] = dsxy2figxy(gca, freqcen-5, ylimd(2)*0.9);
            h = annotation('textbox',[figx figy .01 .01],'String',{['P = ' num2str(stat.negclusters(i).prob,'%.3f')]},'FitBoxToText','on','LineStyle','none','fontsize',8);
        end
    end   
end

xlabel('Frequency (Hz)'); ylabel(ylab); title(titular);
