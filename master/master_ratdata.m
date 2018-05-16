%% RAT DATA ANALYSIS
% clear all;
close all
dbstop if error
% COMMAND PROGRAM - SETS ANALYSIS PROCEDURES

%% Build header to be used for analysis programs
R = buildheader_rat;
% Select methods
section = [3];    % Sections
procs = {[2,3,4,7]...           % Analysis
    [1:2]...           % Statistics
    [2]}...     % Plotting
    ;
for o = 1:length(section)
    switch section(o)
        %% EXTRACTION, PREPROCESSING and ANALYSIS
        case 1
            for i = 1:length(procs{section(o)})
                switch procs{section(o)}(i)
                    %% Data Preparation
                    case 1
                        R.clear.raw = 1;
                        % Now read and convert to FT format
                        convertFT_rat_190717(R)
                    case 2
                        % Preprocess Data
                        R.clear.pp = 1;
                        preprocess_rat_111016b(R)
                    case 3 % Preprocesing Long Epoch
                        %                         R.clear.ppe = 0;
                        %                         R.longE_length = 10;
                        %                         preprocess_longE_rat_050816(R)
                        %% Spectral Analyses
                    case 4 % Spectral analysis including coherence and WPLI
                        R.spectanaly.cplot = [0 0 0]; % Plotting options
                        spectralanalysis_rat_270217(R)
                        R.clear.specstat = 1;
                        stats_spect_rat_050816(R)
                        % - - - - - - - -
                        % Qs:
                    case 5 % DFA-PS
                        R.dfaps.bwid = [1.5 1.5 1.5];
                        R.dfaps.BF_r = -6;
                        R.dfaps.method.cfreq = 'wpli';
                        PEB_DFAPS_rat_011116(R)
                    case 6 % DFA-AE
                        R.dfaae.bwid = [2.5 2.5 2.5];
                        R.dfaae.BF_r = -6;
                        R.dfaae.method.cfreq = 'powcen';
                        PEB_DFAAE_rat_011116(R)
                    case 7 % Directional
                        directional_rat_270217(R)
                end
            end
            %% STATISTICS
        case 2
            for i = 1:length(procs{section(o)})
                switch procs{section(o)}(i)
                    case 1 % Power Stats
                        R.clear.specstat = 1;
                        stats_spect_rat_050816(R)
                    case 2 % Generic Stats
                        R.clear.genstat = 1;
                        R.dfaae.BF_r = -6; R.dfaps.BF_r = -6; R.dfaif.BF_r = -6;
                        R.genpairstatlist = {'COH','ICOH','WPLI','NPD'}; %,'MI'}; %,'MIph'};,
                        stats_generic_pair_rat_270217(R)
                end
            end
            
        case 3
            %% PLOTTING
            % set plot defaults
            set(0,'DefaultAxesFontSize',16)
            set(0,'defaultlinelinewidth',3)
            set(0,'DefaultLineMarkerSize',9)
            set(0, 'DefaultFigurePosition', [12    57   605   550]);
%             figure_clear(R)
            for i = 1:length(procs{section(o)})
                switch procs{section(o)}(i)
                    case 1 % Outcomes of preprocessing
                        dataviewer_rat_071117(R)
                    case 2 % Plot generic connectivity spectra
                        set(0,'DefaultAxesFontSize',16)
                        set(0,'defaultlinelinewidth',3)
                        set(0,'DefaultLineMarkerSize',9)
                        set(0, 'DefaultFigurePosition', [12    57   605   550]);
                        R.spectra.featspecs = {'nsPow'}; %'npdX','npdY','npdZ','npdW',} %'nsPow','nsIcoh','npd','npdX','npdY',,'nsIcoh','npd'}; %,'npdW','npdZ'}; %,'npdX','npdW','npdZ','npdZ'} %,'npdY','npdZ'}; %'nsPow''nsIcoh','npd','nsicoh','npd','npdX','npdY','npdZ','npdW'}; %'ncohXY'}; 'nsPow', %power','coherence','wpli','npd','npdX','npdY','npdZ','npdW'}; %'npd','dtf','power','coherence','wpli','granger','icoherence','npd'};%,'power','coherence','wpli'
%                         plot_gen_rat_statspectra_060317(R)
                        plot_gen_rat_statspectra_200717b(R)
                    case 3 % Boxplots
                        R.boxplot.featspecs =  {'npd'};%,'dfaae','power','coh','wpli','npd','pow'};
                        boxplot_ttest_rat_050816(R);
                end
            end
            close all
            
    end
end
