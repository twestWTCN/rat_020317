function dataviewer_rat_100816(R)
for scale = 1:2
    for cond = 1:2
        for sub  = 1:length(R.subnames{cond})
            set(0,'DefaultAxesFontSize',16)
            set(0,'defaultlinelinewidth',1.5)
            set(0,'DefaultLineMarkerSize',5)
            set(0, 'DefaultFigurePosition', [296         318        1494         678]);
            cmap = linspecer(2,'qualitative');
            
            % Waterfall Plot raw data
            figure
            %         subplot(1,2,cond)
            load([R.analysispath R.pipestamp '\data\processed\' R.subnames{cond}{sub} '_' R.condnames{cond} '_' R.pipestamp '.mat'])
            raw = FTdata;
            processed = FTdata.ContData;
            for i = 1:size(raw.trial{1},1)
                ofset = 0.5*(i-1);
                plot(raw.time{1},raw.trial{1}(i,:)+ofset,'color',cmap(1,:),'linewidth',2.5)
                hold on
                plot(processed.time{1},processed.trial{1}(i,:)+ofset,'color',cmap(2,:),'linewidth',2.5)
            end
            xlabel('Time (s)','FontSize',18); ylabel('Amplitude uV','FontSize',18); ylim([-.5 size(raw.trial{1},1)*0.5])
            if scale == 2
                xlim([30 35])
            end
            title(['Rat ' FTdata.subject],'FontSize',18,'FontWeight','bold');
            h = legend({'raw','preprocessed'});
            set(h,'Position',[0.0118    0.0204    0.1272    0.1212])
            set(gca,'Position',[0.1834    0.0625    0.7855    0.8778])
            set(gca,'YTick',[0:0.5:(size(raw.trial{1},1)-1)*0.5])
            set(gca,'YTickLabel',raw.label)
            a = get(gca,'YTickLabel');
            set(gca,'YTickLabel',a,'fontsize',14);
            set(gcf,'Position',[229         151        1654         845])
           % Save Figures
                        if scale == 2
                            saveallfiguresFIL([R.analysispath R.pipestamp '\results\preprocessing\' FTdata.subject '_short'],'-jpg',0,'-r150'); close all
                        else
                            saveallfiguresFIL([R.analysispath R.pipestamp '\results\preprocessing\' FTdata.subject '_long'],'-jpg',0,'-r150'); close all
                        end
%             close all
            if scale == 1
                %% Plot Spectra
                set(0,'DefaultAxesFontSize',14)
                set(0,'defaultlinelinewidth',2)
                set(0,'DefaultLineMarkerSize',9)
                set(0, 'DefaultFigurePosition', [669         189        1131         705]);
                cmap = linspecer(size(raw.trial{1},1),'qualitative');
                figure
                for i = 1:size(raw.trial{1},1)
                    fxA = log10(FTdata.nsPow.freq);
                    fyA = FTdata.nsPow.Powspctrm(i,:);
                    plot(fxA,fyA,'color',cmap(i,:)); hold on
                end
                xlabel('Frequency (Hz)','fontsize',18); ylabel('log Normalised Magnitude','fontsize',18)
                title(['Rat ' FTdata.subject ' Spectral Estimate'],'FontSize',18,'FontWeight','bold');
                h = legend(raw.label);
                set(h,'FontSize',10)
                xlim(log10([4 100]))
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
                set(gcf,'Position',[1038 363 753 636])
                saveallfiguresFIL([R.analysispath R.pipestamp '\results\preprocessing\' FTdata.subject '_spectra'],'-jpg',0,'-r150'); close all
            end
        end
    end
end
