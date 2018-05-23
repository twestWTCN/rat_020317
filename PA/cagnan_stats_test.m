function [ptt,pvart] = cagnan_stats_test(x,ysave,plotop,panlist,cmap)

for i = 1:size(ysave,2)
    x1 = ysave{1,i};
    x2 = ysave{2,i};
    
    [h ptt(:,i)] =ttest2(x1,x2);
    [h pvart(:,i)] =vartest2(x1,x2);
end

if plotop == 1
    for cond=1:size(ysave,1)
        for i = 1:size(ysave,2)
            alpha = 0.05/size(ptt,1);
            a =  subplot(4,2,panlist(cond,i));
            sigx = x'.*(ptt(:,i)<=alpha);
            sigx(sigx==0) = NaN;
            scatter(sigx,repmat(a.YLim(2)*0.8,size(sigx)),'x','MarkerEdgeColor',cmap(i,:))
            
            sigx = x'.*(pvart(:,i)<=alpha);
            sigx(sigx==0) = [];
            plot(sigx,repmat(a.YLim(2)*0.9,size(sigx)),'*','MarkerEdgeColor',cmap(i,:))
        end
    end
end