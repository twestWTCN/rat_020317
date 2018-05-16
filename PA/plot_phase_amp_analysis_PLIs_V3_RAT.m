% function plot_phase_amp_analysis_PLIs_V3_RAT(R)
% if nargin<1
%     R = makeHeader_SubCort_Cort_Networks();
% end
clear
R = buildheader_rat;

% Collate the points and then do paired ttest between bins? A lot of
% information is lost by cutting off in histograms. Either that or somehow
% increase the number of samples - perhaps include all of the STN channels?
% Maybe combine left and right into one histogram? Collapse all of one
% subjects sync segs into one histogram?? All should be sampling from the
% same distribution so build up better resolved image of underlying 2D PDF.
% Or look at doing a hotelling test between the collated data??

close all
analynames = {'Segment Length','CTX High Beta Amp','STN High Beta Amp','STN Low Beta Amp','STN/CTX High Beta Amp Correlation','Causal Density'};

QX = 8 ; QY = 8; % 6 6
% QX = [-pi -pi/8 pi/8 pi]';
QX = linspace(-pi,pi,QX);
% load([R.datapathr 'subject_hbWPLI075'])
% subscreen = squeeze(sum(subject_hbcohscreen>R.PA.WPLIscreen)==2);
a = logspace(0.5,log10(150),4); logscalez = [-100 -50 -15 15 50 100]; %-40:10:40; %linspace(-50,50,QY); %logscalez = [-3.^(4:-1:2) 3.^(2:4)];
for band = 2; %[1 3]
    for cond =1:2
        for sub  = 1:length(R.subnames{cond})
            nr = 1;
            load([R.analysispath R.pipestamp '\data\processed\' R.subnames{cond}{sub} '_' R.condnames{cond} '_' R.pipestamp '.mat'],'FTdata')
            pA_pli_dist_save = FTdata.PA.pA_pli_dist_save;
            segL_pli_dist_save = FTdata.PA.segL_pli_dist_save;
            amp_pli_dist_save = FTdata.PA.amp_pli_dist_save;
            H_dist_save = FTdata.PA.H_dist_save;
            timevec = FTdata.PA.timevec{1};
            
            relativePhi{nr} = pA_pli_dist_save{1}; %wrapToPi(pA_pli_dist_save{1}-circ_mean(pA_pli_dist_save{1}') ); %circ_mean(pA_pli_dist_save{2,nrOFF}(amp_pli_dist_save{2,nrOFF}(3,:)>95)')); %wrapToPi(pA_pli_dist_save{cond,nr}); %-circ_mean([pA_pli_dist_save{2,nrOFF}]'));
            segL{nr} =  segL_pli_dist_save{1}; %((segL_pli_dist_save{cond,nr} - mean([segL_pli_dist_save{cond,nr}],2))./mean([segL_pli_dist_save{cond,nr}],2)  )*100;
            ampSeg{nr} = amp_pli_dist_save{1}; %((amp_pli_dist_save{cond,nr}  - mean([amp_pli_dist_save{cond,nr}],2))./mean([amp_pli_dist_save{cond,nr}],2) )*100;
            HdistSeg{nr} = H_dist_save{1}(1,:);
%             tends(nr) = timevec{1}(end)- timevec{1}(1);
            tendtot = timevec;
            relativePhiCol = [relativePhi{:,:}];
            segLCol = [segL{:,:}];
            ampSegCol = [ampSeg{:,:}];
            HdistSegCol = [HdistSeg{:,:}];
            
            % Plots of rel phase vs seg length
            f(1) = figure(1*10 + 4); % PA vs SegL
            f(2) = figure((2)*10 +1); set(gcf,'Position',[263 517 1403 352])
            %             subplot(1,3,cond);
            [L4(cond) segLbin binmid pdense(1,cond)] = plot_PA_Dep_relation(cond,relativePhiCol,segLCol,...
                'Segment Length',QX,linspace(0.1,1,QY),f,tendtot); % ,logspace(log10(0.5),log10(1.2),6)
            xlim([-3.5 3.5]); %ylim([0 1.75]);
            figure(f(2)); set(gcf,'Position',[263   343   630   526]);  title(R.condnames{cond})
            caxis([0 0.1]);
            
            % Plots of rel phase vs M1 HB Amp
            f(1) = figure(1*10 + 5); % PA vs M1 Amp
            f(2) = figure((2)*10 +2);
            %             subplot(1,3,cond);
            [L5(cond) ampbin(:,nr,cond,1) binmid pdense(2,cond)] = plot_PA_Dep_relation(cond,relativePhiCol,ampSegCol(1,:),...
                'M1 High Beta Amp',QX,linspace(-50,250,QY),f,tendtot); %linspace(0,7,8) ,linspace(-100,200,QY)
            xlim([-3.5 3.5]); %ylim([0 8])
            figure(f(2)); set(gcf,'Position',[263   343   630   526]);  title(R.condnames{cond})
            caxis([0 0.1]);
            
            % Plots of rel phase vs STN HB Amp
            f(1) = figure(1*10 + 6); % PA vs SegL
            f(2) = figure((2)*10 +3);
            %             subplot(1,3,cond)
            [L6(cond) ampbin(:,nr,cond,2) binmid pdense(3,cond)] = plot_PA_Dep_relation(cond,relativePhiCol,ampSegCol(2,:),...
                'STN High Beta Amp',QX,linspace(-50,300,QY),f,tendtot); %linspace(0,7,8) linspace(-100,200,QY)
            xlim([-3.5 3.5]); %ylim([0 8]);
            figure(f(2)); set(gcf,'Position',[263   343   630   526]);  title(R.condnames{cond})
            caxis([0 0.1]);
            
            % Plots of rel phase vs STN LB Amp
            f(1) = figure(1*10 + 7); % PA vs SegL
            f(2) = figure((2)*10 +4); set(gcf,'Position',[263 517 1403 352])
            %             subplot(1,3,cond)
            [L7(cond) ampbin(:,nr,cond,3) binmid pdense(4,cond)] = plot_PA_Dep_relation(cond,relativePhiCol,ampSegCol(3,:),...
                'STN Low Beta Amp',QX,linspace(-50,200,QY),f,tendtot); % linspace(0,4,8) linspace(-100,200,QY)
            xlim([-3.5 3.5]); %ylim([0 5]);
            figure(f(2)); set(gcf,'Position',[263   343   630   526]);  title(R.condnames{cond})
            caxis([0 0.2]);
            % Plots of rel phase vs Segment Length DDT
            figure(100)
            [H1(cond) hdist(:,cond,1)] = plot_segL_histogram(relativePhiCol,segLCol,'Segment Length',cond);
            
            %Plots of rel phase vs Segment Length DDT
            f(1) = figure(1*10 + 8); % PA vs SegL
            f(2) = figure((2)*10 +5);
            %             subplot(1,3,cond)
            set(gcf,'Position',[263 517 1403 352])
            %                 [H2(cond) hdist(:,cond,2)] = plot_segL_histogram(pA_pli_dist_save{cond,nr},segL_pli_dist_save{cond,nr},'PLI Segment Length',cond);
            [L6(cond) hbin(:,nr,cond,1) binmid pdense(5,cond)] = plot_PA_Dep_relation(cond,relativePhiCol,HdistSegCol,...
                'STN CTX H Beta Amp Corr',QX,linspace(-1,1,QY),f,tendtot); % linspace(0,4,8)
            figure(f(2)); set(gcf,'Position',[263   343   630   526]);  title(R.condnames{cond})
            caxis([0 0.1]);
            %% Phase Angle Rose DDT
            figure(200)
            pA_dist = [pA_pli_dist_save{1}];
            pA_dist(isnan(pA_dist)) = [];
            R_1(cond) = polarhistogram(pA_dist,18,'FaceAlpha',0.75); hold on
            phaseAng_dist(cond).circ_mean_std = [circ_mean(pA_dist') circ_var(pA_dist')];
            [h mu] =circ_mtest(pA_dist',circ_mean(pA_dist'));
            phaseAng_dist(cond).circ_meantest = [h mu];
            %                     [pval, z] = circ_rtest(pA_dist');
            
            phaseAng_dist(cond).circ_RayTest = [1 1]; %[pval, z];
            title('Phase Angle Distribution Phase Lock','FontSize',18)
            %% correlation
            %                     figure;
            %                     x =  amp_dist_save{1,nr}; y = segL_dist_save{1,nr};
            %                     scatter(x(3,:),y);hold on
            %                     x =  amp_dist_save{2,nr}; y = segL_dist_save{2,nr};
            %                     scatter(x(3,:),y);hold on
            figure(252)
            y =  HdistSegCol; x = segLCol;
            
            scatter(x,y(1,:));hold on;[r1 p1] = corrcoef(x,y(1,:)); length_amp_corr(:,cond,sub) = [r1(2) p1(2)];
            xlabel('Segment Length'); ylabel('Amplitude Correlation')
            
%             savefigure_v2([R.analysispath R.pipestamp '\results\figures\PA\subject\'],[R.subnames{cond}{sub} '_PA_analysis'],[],[],[]);
            close all
            % Now shift relative to OFF condition
            segLsave{sub,cond} = segLbin;

            densesave{sub,cond} = pdense;
            hdistsave{sub,cond} = hdist;
            phaseAng_dist_save{sub,cond} = phaseAng_dist;
            %         gc_dist_sub_save{sub,side} = gc_dist_c;
            %         gc_cd_dist_sub_save{sub,side} =gc_dist_cd;
            clear segLbin ampbin pdense hdist phaseAng_dist gc_dist_c gc_dist_cd
        end
    end
    save([R.analysispath R.pipestamp '\results\PAdata\PAgroup'],'segLsave','amp_pli_dist_save','densesave','hdistsave','phaseAng_dist_save'); %,'gc_dist_sub_save','GC_stat_table')
end
    
    
    
    
    
