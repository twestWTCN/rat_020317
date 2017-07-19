function boxplot_ttest_rat_050816(R)
boxy = 1;
close all
% R.boxplot.featspecs = {'psdfa','coh','wpli','pow','powSTN','MI','MIph'};
specs = R.boxplot.featspecs;
for i = 1:length(specs)
    if strncmp(specs{i},'dfaps',7)
        stat_psdfa = stats_group_ON_OFF_ttest_rat(R,'dfaps',{'alpha' 'pebscore','varfeat'},[0.45 1],boxy,1,'blue_up'); %,'pebscore','rejrate'}
        mean_diffstat_tabsave(stat_psdfa,specs{i},R)
        saveallfiguresFIL([R.analysispath R.pipestamp '\results\figures\boxplots\DFAPS\DFAPS_boxplots'],'-jpg'); close all
    elseif strncmp(specs{i},'coh',7)
        stat_coh = stats_group_ON_OFF_ttest_rat(R,'coh',{'banintCoh' 'peakMag'},[],boxy,1,'red_up');
        mean_diffstat_tabsave({stat_coh{:,2}},specs{i},R)
        saveallfiguresFIL([R.analysispath R.pipestamp '\results\figures\boxplots\Coherence\Coherence_boxplots'],'-jpg'); close all
    elseif strncmp(specs{i},'wpli',7)
        stat_wpli = stats_group_ON_OFF_ttest_rat(R,'wpli',{'intCoh' 'peakMag'},[],boxy,1,'red_up');
        mean_diffstat_tabsave({stat_wpli{:,2}},specs{i},R)
        saveallfiguresFIL([R.analysispath R.pipestamp '\results\figures\boxplots\WPLI\WPLI_boxplots'],'-jpg'); close all
    elseif strncmp(specs{i},'pow',7)
        stat_powalphasource = stats_group_ON_OFF_ttest_rat(R,'pow',{'intCoh' 'peakMag'},[],boxy,1,'green_up');
        mean_diffstat_tabsave({stat_powalphasource{:,1}},specs{i},R)
        saveallfiguresFIL([R.analysispath R.pipestamp '\results\figures\boxplots\Power\powsource_boxplots'],'-jpg'); close all
    elseif strncmp(specs{i},'powSTN',7)
        stat_powalphasource = stats_group_ON_OFF_ttest_rat(R,'pow.STN',{'intCoh' 'peakMag' 'peakFrq'},[],boxy,1);
        saveallfiguresFIL([R.analysispath R.pipestamp '\results\figures\boxplots\Power\powSTN_boxplots'],'-jpg'); close all
    elseif strncmp(specs{i},'MI',7)
        stat_MIalphasource = stats_group_ON_OFF_ttest_rat(R,'MI',{'minorm'},[],boxy,1);
        saveallfiguresFIL([R.analysispath R.pipestamp '\results\figures\boxplots\MI\MIalphasource_boxplots'],'-jpg'); close all
    elseif strncmp(specs{i},'MI',7)
        stat_MIalphasource = stats_group_ON_OFF_ttest_rat(R,'MIph',{'minorm'},[],boxy,1);
        saveallfiguresFIL([R.analysispath R.pipestamp '\results\figures\boxplots\MIph\MIphalphasource_boxplots'],'-jpg'); close all
    elseif strncmp(specs{i},'dfaae',7)
        stat_dfaae = stats_group_ON_OFF_ttest_rat(R,'dfaae',{'alpha' 'pebscore','varfeat'},[0.45 1],boxy,1,'blue_down'); %,'pebscore','rejrate'}
        mean_diffstat_tabsave(stat_dfaae,specs{i},R)
        saveallfiguresFIL([R.analysispath R.pipestamp '\results\figures\boxplots\DFAAE\DFAAE_boxplots'],'-jpg'); close all
    elseif strncmp(specs{i},'dfaif',7)
        stat_dfaif = stats_group_ON_OFF_ttest_rat(R,'dfaif',{'alpha' 'pebscore','varfeat'},[0.45 1],boxy,1,'blue_down'); %,'pebscore','rejrate'}
        mean_diffstat_tabsave(stat_dfaif,specs{i},R)
        saveallfiguresFIL([R.analysispath R.pipestamp '\results\figures\boxplots\DFAIF\DFAIF_boxplots'],'-jpg'); close all
    elseif strncmp(specs{i},'npd',7)
        stat_npd = stats_group_ON_OFF_ttest_rat(R,'npd',{'frwd_intCoh','rev_intCoh'},[],boxy,1,'orange_down'); %,'pebscore','rejrate'}
        mean_diffstat_tabsave(stat_npd,specs{i},R)
        saveallfiguresFIL([R.analysispath R.pipestamp '\results\figures\boxplots\NPD\NPD_boxplots'],'-jpg'); close all
    end
end
