clear all; close all
% load('C:\Users\twest\Documents\Work\PhD\LitvakProject\rat_data\pipeline\rat_121016\data\processed\L4_lesion_rat_121016.mat')

load('C:\Users\twest\Documents\Work\PhD\LitvakProject\Analysis_Pipe\Current Pipe\SN_SF\Data\STN_MEGbetasources_JB_Rest1_offPP_FT_SN_SF.mat')
cfg = [];
cfg.length  = 0.5;
FTdataX = ft_redefinetrial(cfg, FTdata);
cfg_red = cfg;

% Preprocess
cfg = [];
cfg.hpfilter  = 'yes';
cfg.hpfreq = 2;
cfg.lpfilter  = 'yes';
cfg.lpfreq = 120;
cfg.dftfilter  = 'yes';
cfg.demean = 'yes';
%         cfg.detrend = 'yes';
FTdataX = ft_preprocessing(cfg,FTdataX);

%         % Frequency
%         cfg            = [];
%         cfg.output     = 'pow';
%         cfg.method     = 'mtmfft';
%         cfg.foilim     = [4 145];
%         cfg.tapsmofrq  = 1;
%         cfg.keeptrials = 'no';
%         cfg.channel    = FTdata.label(1:6);
%         cfg.pad = [];
%         freqPow    = ft_freqanalysis(cfg, FTdataX);
%         cfg_freq = cfg;


W = 4;
L=2*W-1;    % # orthogonal tapers

% xfreqcohere(x(16,:),x(16,:),L)
% FTdataX = FTdata.EpochedData; % for rat
for i = 1
for t = 1:size(FTdataX.trial,2)
x(t,:) = FTdataX.trial{t}(3,:);
end
[ampsny2] = xfreqcohere(x',L);
y(:,:,:,i) = ampsny2;
end
ampsny2m = mean(ampsny2,3);
N=size(ampsny2,1);
freq=linspace(-150,150,N); %
figure
imagesc(freq,freq,abs(ampsny2m))
xlim([-50 50]); ylim([-50 50])
colorbar; caxis([0 0.3])
shg

[trd,rd1,g] = svd(reshape(ampsny2,N^2,size(ampsny2,3)),'econ');
trabs = reshape(trd,N,N,size(ampsny2,3));
figure
for i = 1:6
subplot(2,3,i)
imagesc(freq,freq,squeeze(abs(trabs(:,:,i))))
end