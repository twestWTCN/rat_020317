figure
plot(linspace(0,size(x1AE,2)/fsamp,size(x1AE,2)),x1AE,'r','LineWidth',2);
hold on
plot(linspace(0,size(x1_filt,2)/fsamp,size(x1_filt,2)),x1_filt,'b','LineWidth',2);

 xlim([50 60]); xlabel('Time (s)'); ylabel('Amplitude (mV)')
 set(gcf,'Position',[680   714   787   384])
 set(gca,'FontSize',12)
 figure
 DFAP = [fsamp minBS (length(x1AE)/maxFrac)/fsamp 100 1];
 [bmod win evi alpha] = peb_dfa_cohproj_090616(x1AE,DFAP,BF_r,0);
 set(gcf,'Position',[680   714   787   384])
 set(gca,'FontSize',12)