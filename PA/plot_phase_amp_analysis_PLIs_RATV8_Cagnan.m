function plot_phase_amp_analysis_PLIs_RATV8_Cagnan(R)
if nargin<1
    R = buildheader_rat();
end
close all

QX = 8 ; % Bin Size
for band = 3 %[2 3]
    cmapint = linspecer(3);
    cmap = linspecer(5);
    cmapint(4,:) = cmap(5,:);
    ampBinGroup = []; segBinGroup = []; phipeakGroup = []; RcoeffGroup = [];
    for cond =1:2
        for sub  = 1:length(R.subnames{cond})
            load([R.analysispath R.pipestamp '\data\processed\' R.subnames{cond}{sub} '_' R.condnames{cond} '_' R.pipestamp '.mat'])
            stnlabs = find(strncmp('STN',FTdata.wpli.label,3));
            for stchi = 1:length(stnlabs)
                phi = FTdata.PA(band).pA_pli_dist_save{stchi}';
                relativePhiCol = phi'; %wrapToPi(phi-circ_mean(phi(~isnan(phi'))) )'; %circ_mean(pA_pli_dist_save{2,nrOFF}(amp_pli_dist_save{2,nrOFF}(3,:)>95)')); %wrapToPi(pA_pli_dist_save{cond,nr}); %-circ_mean([pA_pli_dist_save{2,nrOFF}]'));
                segLCol =  FTdata.PA(band).segL_pli_dist_save{stchi}; %((segL_pli_dist_save{cond,nr} - mean([segL_pli_dist_save{cond,nr}],2))./mean([segL_pli_dist_save{cond,nr}],2)  )*100;
                ampSegCol = FTdata.PA(band).amp_pli_dist_save{stchi}; %((amp_pli_dist_save{cond,nr}  - mean([amp_pli_dist_save{cond,nr}],2))./mean([amp_pli_dist_save{cond,nr}],2) )*100;
                %                 HdistSegCol = FTdata.PA(band).H_dist_save(1,:);
                tendtot = FTdata.PA(band).timevec{1}(end)-FTdata.PA(band).timevec{1}(1);
                
                
                ylimlistS{1} = {{[0 0],[0 0],[0 0]},{[-25 75];[-25 75];[-25 75]},{[-50 100];[-50 50];[-50 150]}};
                ylimlistS{2} = {{[0 0]},{[0 2]},{[0 2]}};
                panlist = [1 3 5 ; 2 4 6];
                phiBin = linspace(-pi,pi,QX);
                obs = {[R.sourcenames{1} ' ' R.bandinits{band}],['STN ' R.bandinits{band}],['STN ' R.bandinits{band-1}],[R.bandinits{band} ' Seg. Length']};
                for i = 1:3
                    titleR{1,i} = ['STN-' R.sourcenames{1} ' ' R.bandinits{band} ' Phase vs ' obs{i} ' Power'];
                    titleR{2,i} = ['STN-' R.sourcenames{1} ' ' R.bandinits{band} ' Phase vs ' obs{i} ' Power'];
                    titleR{3,i} = ['STN-' R.sourcenames{1} ' ' R.bandinits{band} ' Phase vs ' R.bandinits{band} ' Frame Length'];
                end
                titleR{4,4} = ['STN-' R.sourcenames{1} ' ' R.bandinits{band} ' Phase vs ' obs{4}];

                [ampBinGroup segBinGroup phipeakGroup RcoeffGroup] =  plot_phase_amp_analysis_generic(R,ylimlistS,obs,band,cond,stchi,sub,phiBin,...
                    cmapint,segLCol,relativePhiCol,ampSegCol,ampBinGroup,segBinGroup,phipeakGroup,RcoeffGroup,...
                    panlist,titleR);
            end
        end
    end
    R.condname = R.condnames;
    plot_phase_amp_analysis_generic_group(R,panlist,obs,phiBin,ampBinGroup,segBinGroup,phipeakGroup,RcoeffGroup)
    close all
end




