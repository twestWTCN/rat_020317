function ax = revcheckplot(OFF,ON,fxA,cmap,ax)
if size(fxA,1)>size(fxA,2)
    fxA = fxA';
end

ax(3) = plot(fxA,nanmean(OFF,1),'Color',cmap(1,:),'LineStyle','--','linewidth',3)
hold on
ax(4) = plot(fxA,nanmean(ON,1),'Color',cmap(2,:),'LineStyle','--','linewidth',3)
