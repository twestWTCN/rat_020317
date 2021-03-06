function plot_example_phaseanalysis_SW(Xdata,amp,phi,dphi_12_dt,seg_ddt,ddphi_ci,timevec) %,PLI_tvec,PLI,consecSegs)     
        ax(1) = subplot(4,1,1);
        cmap = linspecer(2);
        plot(Xdata.time{1},normaliseV(Xdata.trial{1}(1,:)),'color',cmap(1,:));hold on
        plot(Xdata.time{1},normaliseV(Xdata.trial{1}(2,:)),'color',cmap(2,:))
        Amp1 = normaliseV(amp(:,1)'); Amp1 = Amp1-min(Amp1);
        Amp2 = normaliseV(amp(:,2)'); Amp2 = Amp2-min(Amp2);
        plot(Xdata.time{1},Amp1,'color',cmap(1,:)); plot(Xdata.time{1},Amp2,'color',cmap(2,:));
%         TAmp = Amp1.*Amp2;
%         TAmpNorm = (TAmp./mean(TAmp));
%         TAmpNormNeg = TAmpNorm;
%         TAmpNormNeg(TAmpNorm>1) = 1;
        
%         plot(Xdata.time{1},TAmpNormNeg,'k-')
%         xmed1 = median(amp(:,1)); xmed2 = median(amp(:,2));
%         plot(Xdata.time{1},repmat(xmed1,1,size(Xdata.time{1},2)),'--','color',cmap(1,:));
%         plot(Xdata.time{1},repmat(xmed2,1,size(Xdata.time{1},2)),'--','color',cmap(2,:));
        
        %xlim([60 70])
        ylabel('\beta activity','FontSize',14,'FontWeight','bold'); title('Example Analysis')


        ax(2) = subplot(4,1,2);
        plot(Xdata.time{1}(1:length(phi)),phi(:,1),Xdata.time{1}(1:length(phi)),phi(:,2),'linestyle','none','marker','.');%xlim([60 70])
        ylabel('\phi_{1/2}','FontSize',14,'FontWeight','bold');


        % ddt Time Series
        ax(3) = subplot(4,1,3);
        plot(timevec,dphi_12_dt); %xlim([60 70]);
         yvec = nan(size(timevec)); yvec(seg_ddt) = dphi_12_dt(seg_ddt);
         xvec = nan(size(timevec)); xvec(seg_ddt) = timevec(seg_ddt);
        hold on; plot(xvec,yvec,'LineWidth',2)
        hold on; plot([0 timevec(end)],[ddphi_ci ddphi_ci],'k--');
        linkaxes(ax,'x');
        
        set(gcf,'Position',[287 72 1446 932]);
%         xlim([80 85])
