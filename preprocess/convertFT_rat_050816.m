function convertFT_rat_050816(R)
if R.clear.raw == 1
  delete([R.analysispath R.pipestamp '\data\raw\*.mat'])
end
for cond = 1:2
    for sub  = 1:length(R.subnames{cond})
        % Load Data
        load([R.datapath R.subnames{cond}{sub} '.mat'])
%         data = ImportSMR([R.datapath R.subnames{cond}{sub} '.smr']);
%         % HAYRIYES FILTER OUT MULTIUNIT
%         
%         for i = 1:size(SmrData.WvData,1)
%             MuData = makemua_hayriye3(SmrData.WvData(i,:),0.001,0.003,SmrData.SR,SmrData.SR,4)
%         end
% hpsig = unit (already high pass filtered);
% msback = 0.001;
% msfoward = 0.003;
% muasr = 16000;
% uthresh = 4;
% Highpass filter signal

        
        % Create Fieldtrip File
        FTdata.label = SmrData.WvTits;
        FTdata.trial = {SmrData.WvData};
        FTdata.time = {linspace(0,size(SmrData.WvData,2)/SmrData.SR,size(SmrData.WvData,2))};
        FTdata.fsample = SmrData.SR;
        
        % Downsample
        cfg = [];
        cfg.resamplefs = R.pp.ds;
        cfg.detrend  = 'no';
        FTdata = ft_resampledata(cfg, FTdata);
        
        FTdata.subject = R.subnames{cond}{sub};
        FTdata.cond = R.condnames{cond};
        FTdata.crdate = date;
        
        save([R.analysispath R.pipestamp '\data\raw\' R.subnames{cond}{sub} '_' R.condnames{cond} '_' R.pipestamp '.mat'],'FTdata')     
    end
end
