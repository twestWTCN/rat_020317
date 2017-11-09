function [ax] = freqplotlog(fyA,fxA,titular,ylab,ylimd,cmap)

if size(fxA,1)>size(fxA,2)
    fxA = fxA';
end

lineind = find(fxA>48.5 & fxA<51.5);
fyA(:,lineind) = NaN(size(fyA,1),size(lineind,2));


%
fxA =log10(fxA);
% fyA= log10(fyA);
for i = 1:size(fyA,1)
    ax(i) = plot(fxA,fyA(i,:),'color',cmap(i,:),'linewidth',3);
    hold on
    [xCalc,yCalc] = linregress(fxA',fyA(i,:)');
    plot(xCalc(:,2),yCalc,'--','color',cmap(i,:),'linewidth',1);    
end

xlabel('log10 Frequency (Hz)'); ylabel(ylab); title(titular);
if  strcmp(getenv('COMPUTERNAME'), 'SFLAP-2') == 1
    x = gca;
    for i = 1:size(x.XTick,2)
        xtlab{i} = num2str(10^(x.XTick(i)),2);
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