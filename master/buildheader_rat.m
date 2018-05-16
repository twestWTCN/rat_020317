function R = buildheader_rat()
R.pipestamp = 'rat_020317';
R.subnames = {{'C1','C2','C3','C4','C5','C6','C8'}...
    {'L4','L6','L13','L18','L19','L20','L21','L22','L23'}};
% EXCLUDED RAT C7 - IpsiEEG
R.bbounds = [5 10; 14 24; 25 40]; % Bands of interest
R.bandnames = {'Alpha','Low Beta','High beta'};
R.bandnames_nospace = {'alpha','low_beta','high_beta'};
R.source_inport_names = {'fEEG','GP','STR','STN'};
R.sourcenames = {'M2','GPe','STR','STN'};
R.condnames = {'control','lesion'};
R.FOI = [4 100]; % Frequencies of interest
R.pp.ds = 250; % sample rate to downsample to
R.pp.hpfilt = [1 2]; % highpass [stop pass]
R.pp.lpfilt = [100 95]; % lowpass [stop pass]
R.montage = 'monopolar';
R.montname = {'monopolar','bipolar','globalaverage','globaldistal','localdistal','bipolar_fixed','bipolarEEG'}; %,'princomp'};
% number of channels per source per subject nreps{cond}(sub srcloc)
R.nreps{1} = [
    1	8	8	4
    1	9	6	4
    1	9	6	5
    1	8	7	5
    1	9	5	4
    1	7	6	3
    1	9	3	4
    1	9	5	4
    ];
R.nreps{2} = [
    1	9	6	4
    1	9	6	4
    1	10	4	2
    1	9	4	4
    1	11	4	2
    1	10	5	5
    1	9	6	2
    1	8	7	3
    1	9	6	3
    ];
if  strcmp(getenv('COMPUTERNAME'), 'FREE') == 1
    R.datapath = 'C:\home\data\ratdata_050816\';
    R.analysispath = 'C:\Users\twest\Documents\Work\PhD\LitvakProject\rat_data\pipeline\';
else
    R.datapath = 'C:\home\data\ratdata_050816\';
    R.analysispath = 'C:\Users\Tim\Documents\Work\GIT\';
    
end

R.PA.SType = 2; % 1 = sliding window PLI and 2 = SRP
R.PA.bwid = [0.5 0.5 0.5];
R.PA.mwid = 3; % minimum SRP length (cycles)
R.PA.SRPeps = 0.006; %0.006;
R.PA.LowAmpFix = 0; % 1 if SRP is adjusted to account for low amplitude

R.PA.slidingwindow = 1;
R.PA.PLVeps =  0.50;
R.PA.WinOver = 0.98;
% R.PA.stn_lb_frq = 14;
R.PA.SNR_eps = -1;
R.PA.SNReps_prctile = 50;
R.PA.SRPeps_prctile = 10;

R.PA.frqrange{1} = R.bbounds(1,1):0.5: R.bbounds(1,2);
R.PA.frqrange{2} =  R.bbounds(2,1):0.5: R.bbounds(2,2);
R.PA.frqrange{3} =  R.bbounds(3,1):0.5: R.bbounds(3,2);
R.PA.SNR = [-2 -1.5 -1.25]; %-1.5;

