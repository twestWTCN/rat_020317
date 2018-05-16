function dataviewer_rat_071117(R)
for scale = 1%:2
    for cond = 1:2
        for sub  = 1:length(R.subnames{cond})
            set(0,'DefaultAxesFontSize',16)
            set(0,'defaultlinelinewidth',1.5)
            set(0,'DefaultLineMarkerSize',5)
            set(0, 'DefaultFigurePosition', [408.5000  396.5000  954.5000  584.0000]);
            cmap = linspecer(2,'qualitative');
            
            % Waterfall Plot raw data
            figure
            set(gcf,'Position',[408.5000  362.0000  965.5000  636.0000])
            
            %         subplot(1,2,cond)
            load([R.analysispath R.pipestamp '\data\processed\' R.subnames{cond}{sub} '_' R.condnames{cond} '_' R.pipestamp '.mat'])
            raw = FTdata;
            nmelist = {[]}; ilist = [];
            for i = 1:numel(FTdata.label)
                nme = FTdata.label{i};
                if ~any(strncmp(nmelist,nme,3)) | i == 1
                    if length(nme)>3; nme = nme(1:3); end
                    nmelist{end+1} = nme;
                    ilist = [ilist i];
                end
            end
            nmelist = {nmelist{2:end}};
            processed = FTdata.ContData;
            cmap = linspecer(size(ilist,2));
            
            for i = 1:size(ilist,2)
                ofset = 0.5*(i-1);
                plot(raw.time{1},raw.trial{1}(ilist(i),:)+ofset,'color',cmap(2,:),'linewidth',2.5)
                hold on
                %                 plot(processed.time{1},processed.trial{1}(ilist(i),:)+ofset,'color',cmap(2,:),'linewidth',2.5)
            end
            ylim([0-0.5 size(ilist,2)*0.5])
            if scale == 2
                xlim([30 35])
            end
            %             h = legend({'raw da','preprocessed'},'FontSize',18,'FontWeight','bold');
            %             set(h,'Position',[0.0118    0.0204    0.1272    0.1212])
           
            set(gca,'YTick',[0:0.5:(size(ilist,2)-1)*0.5])
            set(gca,'YTickLabel',nmelist)
            a = get(gca,'YTickLabel');
            set(gca,'YTickLabel',a,'fontsize',14);
             set(gca,'Position',[0.1834    0.1108    0.7855    0.8090])
            xlabel('time (s)','fontsize',24,'FontWeight','bold'); ylabel('amplitude (uV)','FontSize',24,'FontWeight','bold'); 

            title([R.condnames{cond} ' animal recording'],'FontSize',28,'FontWeight','bold');
            
            grid on
            
            % time statistics
            tlength(cond,sub) = max(raw.time{1});
            % Save Figures
%             if scale == 2
%                 saveallfiguresFIL([R.analysispath R.pipestamp '\results\preprocessing\' FTdata.subject '_short'],'-jpg',0,'-r150'); close all
%             else
%                 saveallfiguresFIL([R.analysispath R.pipestamp '\results\preprocessing\' FTdata.subject '_long'],'-jpg',0,'-r150'); close all
%             end
            %             close all
            if scale == 1
                %% Plot Spectra
                set(0,'DefaultAxesFontSize',14)
                set(0,'defaultlinelinewidth',2)
                set(0,'DefaultLineMarkerSize',9)
                set(0, 'DefaultFigurePosition', [669         189        1131         705]);
                cmap = linspecer(size(ilist,2));
                figure; set(gcf,'Position',[1038 363 753 636])
                [ax] = freqplotlog(FTdata.nsPow.Powspctrm(ilist,:),FTdata.nsPow.freq,'test',[],[],cmap)
                xlabel('frequency (Hz) (log scale)','fontsize',24,'FontWeight','bold'); ylabel('log normalised power','fontsize',24,'FontWeight','bold')
                title([R.condnames{cond} ' animal spectra'],'FontSize',28,'FontWeight','bold');
                h = legend(ax,nmelist,'FontWeight')
                set(h,'FontSize',18)
                xlim(log10([4 100])); ylim([-5.5 -2])
                %                 l = axes;
                %                 for i = 1:size(raw.trial{1},1)
                %                     fxA = FTdata.nsPow.freq;
                %                     fyA = FTdata.nsPow.Powspctrm(i,:);
                %                     plot(fxA,fyA,'color',cmap(i,:)); hold on
                %                 end
                %                 xlim([3 30])
                xti = get(gca,'XTick')
                for i = 1:size(xti,2)
                    xtlab{i} = num2str(10^(xti(i)),2);
                    %                      xtlab(i) = 10^(xti(i));
                end
                set(gca,'XTickLabel',xtlab)
                grid on
                
%                 saveallfiguresFIL([R.analysispath R.pipestamp '\results\preprocessing\' FTdata.subject '_spectra'],'-jpg',0,'-r150'); close all
            end
            close all
        end
    end
end
a = 1;