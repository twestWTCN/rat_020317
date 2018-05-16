function [ax clustat] = freqclusterplotlog(OFF,ON,fxA,stat,alpha,titular,ylab,ylimd,cmap,dir)
if nargin<5
    alpha = 0.05;
end
if nargin<6
    titular = 'spectral stats for unknown analysis';
end
if nargin<6
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


% Now plot
fxA =log10(fxA);
figure
[ax(2), hp] = boundedline(fxA,-mean(log10(OFF)),std(log10(OFF))/sqrt(log10(size(OFF,1))),'cmap',cmap(1,:),'alpha','transparency',0.45);
hold on
[xCalc,yCalc] = linregress(fxA',-mean(log10(OFF))');
plot(xCalc(:,2),yCalc,'--','color',cmap(1,:),'linewidth',2);

%     hl(1).LineWidth = 2;
%     hout = outlinebounds(hl, hp);
hold on
% plot(repmat(fxA,size(OFF,1),1)',-log10(OFF)','Color',cmap(1,:),'LineStyle','--','linewidth',1);

[ax(1), hp] = boundedline(fxA,-mean(log10(ON)),std(log10(ON))/sqrt(log10(size(ON,1))),'cmap',cmap(2,:),'alpha','transparency',0.45);
hold on
[xCalc,yCalc] = linregress(fxA',-mean(log10(ON))');
plot(xCalc(:,2),yCalc,'--','color',cmap(2,:),'linewidth',2);

%     hl(1).LineWidth = 2;
%     hout = outlinebounds(hl, hp);
hold on
% plot(repmat(fxA,size(ON,1),1)',-log10(ON)','Color',cmap(2,:),'LineStyle','--','linewidth',1);

% emsure mask at 49-51 is zero
stat.mask = stat.prob<alpha
stat.mask(stat.freq> 49 & stat.freq < 51) = 0;

sigfreq = stat.mask.*stat.freq;
sigfreq(sigfreq==0) = NaN;
sigpow = stat.mask.*(nanmean(nanmean(ON,1),2));
if ~isempty(ylimd)
    sigpow = stat.mask.*(ylimd(2))*1.42;
end
% sigpow = stat.mask.*1; %(nanmean(nanmean(ON,1),2).*7.5);

sigpow(sigpow==0) = NaN;
if sum(~isnan(sigpow))<2
    scatter(log10(sigfreq),sigpow,'ks','filled');
else
    plot(log10(sigfreq),sigpow,'k','linewidth',6);
end
if ~isempty(ylimd)
    ylim(ylimd);
end
xlim(log10([4 100]));
clustat = [];
if sum(stat.mask)>0
    for i = 1:size(stat.posclusters,2)
        if stat.posclusters(i).prob<alpha
            labs = stat.posclusterslabelmat;
            freqcen = mean(stat.freq(1,find(labs==i)));
            freqcen = log10(freqcen);
            [figx figy] = dsxy2figxy(gca, freqcen-0.12, (ylimd(2)*1.25)); %-log10(5)
            h = annotation('textbox',[figx figy .01 .01],'String',{['P = ' num2str(stat.posclusters(i).prob,'%.3f')]},'FitBoxToText','on','LineStyle','none','fontsize',12,'fontweight','bold');
            clustat = [clustat; min(stat.freq(1,find(labs==i))) max(stat.freq(1,find(labs==i))) stat.posclusters(i).clusterstat stat.posclusters(i).prob];
            
        end
    end
    for i = 1:size(stat.negclusters,2)
        if stat.negclusters(i).prob<alpha
            labs = stat.negclusterslabelmat;
            freqcen = mean(stat.freq(1,find(labs==i)));
            freqcen = log10(freqcen);
            [figx figy] = dsxy2figxy(gca, freqcen-0.12, (ylimd(2)*1.25));
            h = annotation('textbox',[figx figy .01 .01],'String',{['P = ' num2str(stat.negclusters(i).prob,'%.3f')]},'FitBoxToText','on','LineStyle','none','fontsize',12,'fontweight','bold');
            clustat = [clustat; min(stat.freq(1,find(labs==i))) max(stat.freq(1,find(labs==i))) stat.posclusters(i).clusterstat stat.negclusters(i).prob];
        end
    end
end

xlabel('log10 Frequency (Hz)'); ylabel(ylab); title(titular);
if  strcmp(getenv('COMPUTERNAME'), 'SFLAP-2') == 1
    x = gca;
    for i = 1:size(x.XTick,2)
        xtlab{i} = sprintf('%1.f',10^(x.XTick(i)));
    end
    x.XTickLabel = xtlab;
else
    xti = get(gca,'XTick')
    for i = 1:size(xti,2)
        xtlab{i} = num2str(10^(xti(i)),2);
        %                      xtlab(i) = 10^(xti(i));
    end
    set(gca,'XTickLabel',xtlab)
end
grid on