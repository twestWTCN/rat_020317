clear all
cmap = linspecer(3);
R = buildheader_rat;

load('C:\Users\twest\Documents\Work\PhD\temp-\fxA.mat')
load('C:\Users\twest\Documents\Work\PhD\temp-\OFF_GPe_STN.mat')
GPe_STN = ON;
load('C:\Users\twest\Documents\Work\PhD\temp-\OFF_STN_GPe.mat')
STN_GPe = ON;

[specstat] = spectralclusterstats060317(GPe_STN',STN_GPe',fxA{1}',R.sourcenames{2},5000);
[ax clustat] = freqclusterplot_200717(STN_GPe',GPe_STN',fxA{1},specstat,0.05,[],'NPD',[0 .5],cmap);
