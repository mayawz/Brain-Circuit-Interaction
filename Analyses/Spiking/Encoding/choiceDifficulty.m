% OFC-PCC paper reply to reviewer analysis
% last used: 2021/11/03

% beta correlation for easy vs. difficulty choices for cho 12 in OFCin
% take the rho to predict choLR
% see if easier choice 12 drives larger choice lr signal

clear all; close all; clc
% addpath(genpath('/'));

dpath = '/Users/wangz37/OneDrive - National Institutes of Health/_paper_submission/OFC_PCC/2021_R2R/wrapped/PScombined/';
spath = '/Users/wangz37/OneDrive - National Institutes of Health/_paper_submission/OFC_PCC/2021_R2R/R2R_result/';
fpath = '/Users/wangz37/OneDrive - National Institutes of Health/_paper_submission/OFC_PCC/2021_R2R/R2R_code/';
addpath(genpath('/Users/wangz37/OneDrive - National Institutes of Health/01Code/General_functions/'));
addpath(genpath('/Users/wangz37/OneDrive - National Institutes of Health/_paper_submission/OFC_PCC/2021_R2R/00 GeneralFunctions/'))
addpath(fpath);
cd(dpath);
% atOFFER1: opt1=301:400;opt2=401:500; chocie: 500:614; outcome: 615:730 (according to trial length)
% atCHOICE: choice=200:300; outcome=300:400; ITI=400:500
% 300 is when choice strobe is sent

%%
ofcin = load([dpath 'PS.StagOps.OFCin_atOFFER1.mat']);
pcc = load([dpath 'PS.StagOps.PCC_atOFFER1.mat']);
%% data

clc
md_in = fitRegs_Difficulty(ofcin.data);

md_pcc = fitRegs_Difficulty(pcc.data);


%%
save('ChoiceDifficulty.mat','md_in', 'md_pcc')

%% ###########    OFC    ############
close all; clc
% load('ChoiceDifficulty.mat')

% 'Kendall'
% 'Spearman'
% 'Pearson'
for tp = 1:size(md_in(1).hard.bev1,2)
[r_ofc.hard(tp),p_ofc.hard(tp)] = corr(md_in(1).hard.bev1(:,tp),md_in(1).hard.bev2(:,tp),'Type','Spearman');
[r_ofc.easy(tp),p_ofc.easy(tp)] = corr(md_in(1).easy.bev1(:,tp),md_in(1).easy.bev2(:,tp),'Type','Spearman');
end


%%
close all
histogram(r_ofc.hard)
hold on
histogram(r_ofc.easy)
legend('hard','easy')
%% #############    PCC    ##########
clear r p
close all; clc
% load('ChoiceDifficulty.mat')

% 'Kendall'
% 'Spearman'
% 'Pearson'
for tp = 1:size(md_pcc(1).hard.bev1,2)
[r.hard(tp),p.hard(tp)] = corr(md_pcc(1).hard.bev1(:,tp),md_pcc(1).hard.bev2(:,tp),'Type','Spearman');
[r.easy(tp),p.easy(tp)] = corr(md_pcc(1).easy.bev1(:,tp),md_pcc(1).easy.bev2(:,tp),'Type','Spearman');
end


%%
close all
histogram(r.hard)
hold on
histogram(r.easy)
legend('hard','easy')


%% #######################
clear r p
close all; clc
% load('ChoiceDifficulty.mat')

% 'Kendall'
% 'Spearman'
% 'Pearson'
for tp = 1:size(md_pcc(1).hard.bev1,2)
[r_pcc.hard(tp),p_pcc.hard(tp)] = corr(md_pcc(2).hard.blev(:,tp),md_pcc(2).hard.brev(:,tp),'Type','Spearman');
[r_pcc.easy(tp),p_pcc.easy(tp)] = corr(md_pcc(2).easy.blev(:,tp),md_pcc(2).easy.brev(:,tp),'Type','Spearman');
end


%%
close all
histogram(r_pcc.hard)
hold on
histogram(r_pcc.easy)
legend('hard','easy')


%% granger causality
clc
% for nlag=1:48
for nlag=1:99

%     [~,pvalue1(nlag),stat1(nlag),~]=gctest(In_rEV1EV2',Pcc_rEVLEVR','NumLags',nlag,'Trend',1);
%     
%     [~,pvalue,stat2(nlag),~]=gctest(Pcc_rEVLEVR',In_rEV1EV2','NumLags',nlag,'Trend',1);
%     
%     [~,pvalue,stat3(nlag),~]=gctest(Out_rEV1EV2',Pcc_rEVLEVR','NumLags',nlag,'Trend',1);
    
    

%     [h,pvalue,stat,cvalue]=gctest(Out_rEV1EV2',Pcc_rEVLEVR',...
%         [In_rEVLEVR',Pcc_rEV1EV2',],'NumLags',nlag,'Trend',1);
%     [h,pvalue,stat,cvalue]=gctest(In_rEV1EV2',Pcc_rEVLEVR',...
%         [In_rEVLEVR',Pcc_rEV1EV2',],'NumLags',nlag,'Trend',1);
%     [h,pvalue,stat,cvalue]=gctest(Pcc_rEVLEVR',In_rEV1EV2',...
%         [In_rEVLEVR',Pcc_rEV1EV2',],'NumLags',nlag,'Trend',1);



%   [h,pvalue,stat,cvalue]=gctest(varChoOpt.in_t',varChoLR.pcc_t',...
%         [varChoLR.in_t', varChoOpt.pcc_t',varChoLR.out_t',varChoOpt.out_t'],...
%         'NumLags',nlag,'Trend',1);
    
%     [h,pvalue,stat,cvalue]=gctest(varChoLR.pcc_t',varChoOpt.in_t',...
%         [varChoLR.in_t', varChoOpt.pcc_t',varChoLR.out_t',varChoOpt.out_t'],...
%         'NumLags',nlag,'Trend',1);
    
%     [h,pvalue,stat,cvalue]=gctest(varChoLR.pcc_t',varChoOpt.in_t',...
%         'NumLags',nlag,'Trend',1);

% [h,pvalue,stat,cvalue]=gctest(r_ofc.hard',r_pcc.hard','NumLags',nlag,'Trend',1);    
    
   [h,pvalue,stat,cvalue]=gctest(r_ofc.easy',r_pcc.easy','NumLags',nlag,'Trend',1);   
    
    if pvalue<0.025
        disp('sig')
        nlag
    end
end

disp('##################### DONE #####################')



%%
[h,pvalue,stat,cvalue]=gctest(r_ofc.hard',r_pcc.hard','NumLags',93,'Trend',1) 

disp('########')
[h,pvalue,stat,cvalue]=gctest(r_ofc.easy',r_pcc.easy','NumLags',79,'Trend',1)








%% 200-neuron EV1-EV2
close all; clc
cd(spath)
% Calculate the 42nd percentile. Y = prctile(x,42);
sigh_in = prctile(r_in200_ev12,2.5);
sigl_in = prctile(r_in200_ev12,97.5);
histogram(r_in200_ev12,'BinWidth', 0.02, 'faceColor',[0.6350 0.0780 0.1840],'lineStyle','none');
xlim([-1,1]);
hold on
vline(-0.36,[1 0 0]);
vline(sigh_in);
vline(sigl_in);
sigh_out = prctile(r_out200_ev12,2.5);
sigl_out = prctile(r_out200_ev12,97.5);
histogram(r_out200_ev12,'BinWidth', 0.02, 'faceColor',[0.8500 0.3250 0.0980],'lineStyle','none');
vline(-0.18,[255,140,0]/255);
vline(sigh_out,[0.5,0.5,0.5]);
vline(sigl_out,[0.5,0.5,0.5]);
legend('cOFCm','cOFCm data', 'cOFCm-sig-cutoff','cOFCm-sig-cutoff',...
    'cOFCl','cOFCl data', 'cOFCl-sig-cutoff','cOFCl-sig-cutoff');
xlabel('Spearman Correlation Coefficient','fontSize',16);
title('EV1 vs. EV2 - MCresampled200','fontSize',26);
hold off
saveas(gcf,[spath 'EV12_200neuron'],'epsc')

[h,p,ks2stat] = kstest2(r_in200_ev12,r_out200_ev12)

%% Jackknife EV1-EV2
close all; clc
cd(spath)
% Calculate the 42nd percentile. Y = prctile(x,42);
sigh_in = prctile(r_inJK_ev12,2.5);
sigl_in = prctile(r_inJK_ev12,97.5);
histogram(r_inJK_ev12,'BinWidth', 0.02, 'faceColor',[0.6350 0.0780 0.1840],'lineStyle','none');
xlim([-1,1]);
hold on
vline(-0.36,[1 0 0]);
vline(sigh_in);
vline(sigl_in);
sigh_out = prctile(r_outJK_ev12,2.5);
sigl_out = prctile(r_outJK_ev12,97.5);
histogram(r_outJK_ev12,'BinWidth', 0.02, 'faceColor',[0.8500 0.3250 0.0980],'lineStyle','none');
vline(-0.18,[255,140,0]/255);
vline(sigh_out,[0.5,0.5,0.5]);
vline(sigl_out,[0.5,0.5,0.5]);
legend('cOFCm','cOFCm data', 'cOFCm-sig-cutoff','cOFCm-sig-cutoff',...
    'cOFCl','cOFCl data', 'cOFCl-sig-cutoff','cOFCl-sig-cutoff');
xlabel('Spearman Correlation Coefficient','fontSize',16);
title('EV1 vs. EV2 - Jackknife','fontSize',26);
hold off
saveas(gcf,[spath 'EV12_Jackknife'],'epsc')

[h,p,ks2stat] = kstest2(r_inJK_ev12,r_outJK_ev12)

%% 200-neuron EVL-EVR
close all; clc
cd(spath)
% Calculate the 42nd percentile. Y = prctile(x,42);
sigh_in = prctile(r_in200_evlr,2.5);
sigl_in = prctile(r_in200_evlr,97.5);
histogram(r_in200_evlr,'BinWidth', 0.02, 'faceColor',[0 0.4470 0.7410],'lineStyle','none');
xlim([-1,1]);
hold on
sigh_out = prctile(r_out200_evlr,2.5);
sigl_out = prctile(r_out200_evlr,97.5);
histogram(r_out200_evlr,'BinWidth', 0.02, 'faceColor',[0.3010 0.7450 0.9330],'lineStyle','none');
vline(-0.16,[0 0 1]);
vline(sigh_in);
vline(sigl_in);

vline(0.10,[30,144,255]/255);
vline(sigh_out,[0.5,0.5,0.5]);
vline(sigl_out,[0.5,0.5,0.5]);
legend('cOFCm','cOFCl','cOFCm data', 'cOFCm-sig-cutoff','cOFCm-sig-cutoff',...
    'cOFCl data', 'cOFCl-sig-cutoff','cOFCl-sig-cutoff');
xlabel('Spearman Correlation Coefficient','fontSize',16);
title('EVl vs. EVr - MCresampled200','fontSize',26);
hold off
saveas(gcf,[spath 'EVlr_200neuron'],'epsc')

[h,p,ks2stat] = kstest2(r_in200_evlr,r_out200_evlr)

%% Jackknife EVL-EVR
close all; clc
cd(spath)
% Calculate the 42nd percentile. Y = prctile(x,42);
sigh_in = prctile(r_inJK_evlr,5);
sigl_in = prctile(r_inJK_evlr,95);
histogram(r_inJK_evlr,'BinWidth', 0.02, 'faceColor',[0 0.4470 0.7410],'lineStyle','none');
xlim([-1,1]);
hold on
sigh_out = prctile(r_outJK_evlr,2.5);
sigl_out = prctile(r_outJK_evlr,97.5);
histogram(r_outJK_evlr,'BinWidth', 0.02, 'faceColor',[0.3010 0.7450 0.9330],'lineStyle','none');
vline(-0.16,[0 0 1]);
vline(sigh_in);
vline(sigl_in);

vline(0.10,[30,144,255]/255);
vline(sigh_out,[0.5,0.5,0.5]);
vline(sigl_out,[0.5,0.5,0.5]);

legend('cOFCm','cOFCl','cOFCm data', 'cOFCm-sig-cutoff','cOFCm-sig-cutoff',...
    'cOFCl data', 'cOFCl-sig-cutoff','cOFCl-sig-cutoff');
xlabel('Spearman Correlation Coefficient','fontSize',16);
title('EVl vs. EVr - Jackknife','fontSize',26);
hold off
saveas(gcf,[spath 'EVlr_Jackknife'],'epsc')

[h,p,ks2stat] = kstest2(r_inJK_evlr,r_outJK_evlr)


%%
median([0.68, 0.33, 0.41, 0.31, 0.27, 0.36, 0.2])







