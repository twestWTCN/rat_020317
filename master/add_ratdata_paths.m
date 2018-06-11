function [] = add_ratdata_paths()
if strcmp(getenv('COMPUTERNAME'),'SFLAP-2')
    addpath(genpath('C:\Users\twest\Documents\Work\PhD\LitvakProject\rat_data\pipeline\rat_020317'))
    rmpath(genpath('C:\Users\twest\Documents\Work\PhD\LitvakProject\rat_data\pipeline\rat_020317\Legacy'))
    pathCell = regexp(path, pathsep, 'split'); onPath = any(strcmpi('C:\Users\Tim\Documents\MATLAB_ADDONS\mvgc_v1.0', pathCell));
    
    % Grangerpath
    if ~onPath; addpath('C:\Users\Tim\Documents\MATLAB_ADDONS\mvgc_v1.0'); run startup.m; end
    addpath('C:\Users\Tim\Documents\MATLAB_ADDONS\gramm-master\@gramm')
    
elseif strcmp(getenv('COMPUTERNAME'),'FREE')
    addpath('C:\shared')
    addpath(genpath('C:\Users\twest\Documents\Work\PhD\LitvakProject\rat_data\pipeline\rat_020317'))
    rmpath(genpath('C:\Users\twest\Documents\Work\PhD\LitvakProject\rat_data\pipeline\rat_020317\Legacy'))
    addpath(genpath('C:\Users\twest\Documents\Work\Github\DynPhaseConnectivity'))
    addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\mvgc_v1.0')
    addpath('C:\shared')
    run startup.m
    addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\gramm-master\@gramm')
    pathCell = regexp(path, pathsep, 'split'); onPath = any(strcmpi('C:\spm12', pathCell));
    
    if ~onPath; addpath('C:\spm12'); spm eeg; close all; end
    addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\Circular_Statistics_Toolbox')
    addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\SplitVec')
    addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\FMINSEARCHBND')
    addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\aboxplot')
    addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\ATvDFA-package')
    addpath(genpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\boundedline-pkg'))
    addpath(genpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\DylanMuir-ParforProgMon-9a1c257'))
    addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\ATvDFA-package')
    
    addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\linspecer')
    addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\TWtools')
    addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\sigstar-master')
    addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\Neurospec\neurospec21')
    addpath(genpath('C:\spm12\toolbox\xjview96\xjview'))
end