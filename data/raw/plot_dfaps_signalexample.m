figure
plot(linspace(0,size(dphdiff,2)/fsamp,size(dphdiff,2)),dphdiff,'LineWidth',2);
 xlim([50 60]); xlabel('Time (s)'); ylabel('d(\phi_1 - \phi_2)/dt (rads)')
 set(gcf,'Position',[680   714   787   384])
 set(gca,'FontSize',12)
 figure
 [bmod win evi alpha] = peb_dfa_cohproj_090616(dphdiff,DFAP,BF_r,0);
 set(gcf,'Position',[680   714   787   384])
 set(gca,'FontSize',12)