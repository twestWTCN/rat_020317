%function dwell = getHists_phase_amp_RAT(R)
% consider saving the difference of ON and OFF - recentre that to some
% peak? Or the mean of the ON and OFF together?
clear
R = buildheader_rat;

close all
QX = 6 ; QY = 8;
% load([R.datapathr 'subject_hbWPLI075'])
% subscreen = squeeze(sum(subject_hbcohscreen>R.PA.WPLIscreen)==2);
a = logspace(0.5,log10(150),4); logscalez = -80:20:80; %linspace(-50,50,QY); %logscalez = [-3.^(4:-1:2) 3.^(2:4)];
condcr = {'r','b'};
for band = 3; %[1 3]
    for cond =1:2
        for sub  = 1:length(R.subnames{cond})
            nr = 1;
            load([R.analysispath R.pipestamp '\data\processed\' R.subnames{cond}{sub} '_' R.condnames{cond} '_' R.pipestamp '.mat'],'FTdata')
            pA_pli_dist_save = FTdata.PA.pA_pli_dist_save;
            segL_pli_dist_save = FTdata.PA.segL_pli_dist_save;
            amp_pli_dist_save = FTdata.PA.amp_pli_dist_save;
            H_dist_save = FTdata.PA.H_dist_save;
            timevec = FTdata.PA.timevec{1};
            
            
            
            segL = segL_pli_dist_save{1};
            x = timevec{1}(end);
            dwellrat(cond,sub) = sum(segL)/x;
            
            figure(10)
            subplot(1,2,cond)
            phase_ang = pA_pli_dist_save{1};
            histogram(phase_ang,linspace(-pi,pi,6),'FaceColor',condcr{cond},'FaceAlpha',0.1,'Normalization','Probability'); hold on
            xlabel('Phase Angle');ylabel('P(X)'); ylim([0 0.5]); grid on; title(R.condnames{cond})
            [N B] = histcounts(phase_ang,linspace(-pi,pi,6));
            phaseang_dist(cond,sub).N = N;
            phaseang_dist(cond,sub).B = B;
            segL = segL_pli_dist_save{1};
            
            figure(20)
            subplot(1,2,cond)
            histogram(segL,logspace(log10(0.05),log10(1.5),16),'FaceColor',condcr{cond},'FaceAlpha',0.1,'Normalization','Probability'); hold on
            set(gca,'xscale','log')
            xlabel('Segment Length');ylabel('P(X)');  grid on; title(R.condnames{cond});ylim([0 0.5]); xlim([-1 1.5]);
            [N B] = histcounts(segL,linspace(0,3.5,12));
            segL_dist(nr,cond,sub).N = N;
            segL_dist(nr,cond,sub).B = B;
            
            ampname = {'HB CTX','HB STN','LB STN'};
            for i = 1:3
                amp = amp_pli_dist_save{1}(i,:);
                figure(30+i)
                subplot(1,2,cond)
                histogram(amp,linspace(-100,250,12),'FaceColor',condcr{cond},'FaceAlpha',0.1,'Normalization','Probability'); hold on
                xlabel(ampname{i});ylabel('P(X)'); ylim([0 0.7]); grid on; title(R.condnames{cond})
                [N B] = histcounts(amp,linspace(-100,100,14));
                amp_dist(i,nr,cond,sub).N  = N;
                amp_dist(i,nr,cond,sub).B  = B;
            end
            
            H = H_dist_save{1}(1,:);
            figure(40)
            subplot(1,2,cond)
            histogram(H,linspace(0,1,12),'FaceColor',condcr{cond},'FaceAlpha',0.1,'Normalization','Probability'); hold on
            xlabel('Amp Env Correlation');ylabel('P(X)'); ylim([0 0.3]); grid on; title(R.condnames{cond})
            [N B] = histcounts(segL,linspace(0,1,14));
            H_dist(cond,sub).N = N;
            H_dist(cond,sub).B = B;
            Hcorr(nr,cond,sub) = mean(H(1,:))
            
            figure(252)
            subplot(1,2,cond)
            y =  H_dist_save{1}; x = segL_pli_dist_save{1};
            scatter(x,y(1,:),'filled','MarkerFaceColor',condcr{cond},'MarkerFaceAlpha',0.1);hold on;[r1 p1] = corrcoef(x,y(1,:)); length_amp_corr(:,cond,sub,nr) = [r1(2) p1(2)];
            xlabel('Segment Length'); ylabel('Amplitude Correlation'); grid on; title(R.condnames{cond})
            set(gca,'xscale','log'); xlim([0.2 1.5])
            for i=1:3
                figure(260+i)
                subplot(1,2,cond)
                y =  H_dist_save{1}; x = amp_pli_dist_save{1}(i,:); %segL_pli_dist_save{cond,nr};
                scatter(x,y(1,:),'filled','MarkerFaceColor',condcr{cond},'MarkerFaceAlpha',0.1);hold on;[r1 p1] = corrcoef(x,y(1,:)); amp_amp_corr(:,i,cond,sub,nr) = [r1(2) p1(2)];
                xlabel(ampname{i}); ylabel('Amplitude Correlation'); grid on; title(R.condnames{cond});xlim([-100 250]); ylim([-1 1])
            end
        end
    end
    % save([R.datapathr '\results\seganalysis\groupseganaly'],'ampsave','densesave','hdistsave','phaseAng_dist_save'); %,'gc_dist_sub_save','GC_stat_table')
    
    
    x = (squeeze(dwellrat(1,:))); x(x==0) = [];
    y = (squeeze(dwellrat(2,:))); y(y==0) = [];
    
    figure(1)
    subplot(1,2,1)
    histogram(x,linspace(0.1,0.3,4),'FaceColor',condcr{1},'Normalization','Probability'); xlim([0 1])
    xlabel('Dwell/Escape'); ylabel('P(X)'); grid on; ylim([0 0.7])
    subplot(1,2,2);
    histogram(y,linspace(0.1,0.3,4),'FaceColor',condcr{2},'Normalization','Probability'); xlim([0 1]);
    xlabel('Dwell/Escape'); ylabel('P(X)'); grid on; ylim([0 0.7])
    [h p] = ttest2(x,y)
    title(num2str(p))
    dwell = {x,y,p};
    
    x = (squeeze(Hcorr(:,1,:))); x(x==0) = [];
    y = (squeeze(Hcorr(:,2,:))); y(y==0) = [];
    
    figure(2)
    subplot(1,2,1)
    histogram(x,linspace(0.2,0.5,5),'FaceColor',condcr{1},'Normalization','Probability'); xlim([0 1])
    xlabel('STN/M1 Amp Corr'); ylabel('P(X)'); grid on; ylim([0 0.7])
    subplot(1,2,2);
    histogram(y,linspace(0.2,0.5,5),'FaceColor',condcr{2},'Normalization','Probability'); xlim([0 1]);
    xlabel('STN/M1 Amp Corr'); ylabel('P(X)'); grid on; ylim([0 0.7])
    [h p] = ttest2(x,y)
    title(num2str(p))
    dwell = {x,y,p};
    savefigure_v2([R.datapathr '\results\seganalysis\PLI\partests\'],[ num2str(idd) '_PLI_surrogate_bwid_' num2str(R.PA.bwid) '_PLVeps_' num2str(R.PA.PLVeps)],[],[],[]); close all
end
% figure(3)
% x = sum(squeeze(length_amp_corr(1,1,:,:)),3);x(x==0) = [];
% y = sum(squeeze(length_amp_corr(1,2,:,:)),3);y(y==0) = [];
% subplot(1,2,1)
% histogram(x,linspace(-1,1,12),'FaceColor',condcr{1},'Normalization','Probability');% xlim([0 1])
% xlabel('Frame Length'); ylabel('Amplitude Correlation'); grid on;% ylim([0 0.7])
% subplot(1,2,2);
% histogram(y,linspace(-1,1,12),'FaceColor',condcr{2},'Normalization','Probability'); %xlim([0 1]);
% xlabel('Frame Length'); ylabel('Amplitude Correlation'); grid on; %ylim([0 0.7])
% [h p] = ttest2(x,y)
% title(num2str(p))

%
% for i = 1:3
% figure(400+i)
% x = sum(squeeze(amp_amp_corr(1,i,1,:,:)),3);x(x==0) = [];
%     y = sum(squeeze(amp_amp_corr(1,i,2,:,:)),3);y(y==0) = [];
%     subplot(1,2,1)
%     histogram(x,linspace(-1,1,12),'FaceColor',condcr{1},'Normalization','Probability');% xlim([0 1])
%     xlabel(ampname{i}); ylabel('Amplitude Correlation'); grid on;% ylim([0 0.7])
%     subplot(1,2,2);
%     histogram(y,linspace(-1,1,12),'FaceColor',condcr{2},'Normalization','Probability'); %xlim([0 1]);
%     xlabel(ampname{i}); ylabel('Amplitude Correlation'); grid on; %ylim([0 0.7])
%     [h p] = ttest2(x,y)
%     title(num2str(p))
% end
% savefigure_v2([R.datapathr '\results\seganalysis\PLI\partests\'],[ num2str(idd) '_PLI_surrogate_.tiff'],[],[],[]); close all

