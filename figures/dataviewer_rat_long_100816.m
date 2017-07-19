function dataviewer_rat_long_100816(R)
set(0,'DefaultAxesFontSize',16)
set(0,'defaultlinelinewidth',1.5)
set(0,'DefaultLineMarkerSize',5)
set(0, 'DefaultFigurePosition', [296         318        1494         678]);
for cond = 1:2
    for sub  = 1:length(R.subnames{cond})
        % Plot raw data
        figure(sub+cond*10)
        %         subplot(1,2,cond)
        load([R.analysispath R.pipestamp '\data\raw\' R.subnames{cond}{sub} '_' R.condnames{cond} '_' R.pipestamp '.mat'])
        raw = FTdata;
        load([R.analysispath R.pipestamp '\data\processed\' R.subnames{cond}{sub} '_' R.condnames{cond} '_' R.pipestamp '.mat'])
        processed = FTdata.ContData;
        for i = 1:size(raw.trial{1},1)
            ofset = 0.5*(i-1);
            plot(raw.time{1},raw.trial{1}(i,:)+ofset,'b')
            hold on
            plot(processed.time{1},processed.trial{1}(i,:)+ofset,'r')
            xlabel('Time (s)'); ylabel('Amplitude uV'); ylim([-.5 size(raw.trial{1},1)*0.5])
            title(['Rat ' FTdata.subject]);
            h = legend({'raw','preprocessed'});
            set(h,'Position',[0.0118    0.0204    0.1272    0.1212])
            set(gca,'Position',[0.1834    0.0625    0.7855    0.8778])
            set(gca,'YTick',[0:0.5:(size(raw.trial{1},1)-1)*0.5])
            set(gca,'YTickLabel',raw.label)
        end
        
    end
end