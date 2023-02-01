% OFC-PCC paper reply to reviewer analysis
% last used: 2021/10/14

% this addresses "sample too small"
% Monte-Carlo resampling
% Jackknife -- take one out.
% 200-neuron resampling with replacement 1000 times

% reversed power analysis -- in R

clear all; close all; clc
% addpath(genpath('/'));

dpath = '/Users/wangz37/OneDrive - National Institutes of Health/_paper_submission/OFC_PCC/2021_R2R/wrapped/PScombined/';
spath = '/Users/wangz37/OneDrive - National Institutes of Health/_paper_submission/OFC_PCC/2021_R2R/R2R_result/';
fpath = '/Users/wangz37/OneDrive - National Institutes of Health/_paper_submission/OFC_PCC/2021_R2R/R2R_code/';
addpath(genpath('/Users/wangz37/OneDrive - National Institutes of Health/01Code/General_functions/'));
addpath(fpath);
cd(dpath);
% atOFFER1: opt1=301:400;opt2=401:500;
% atCHOICE: choice=200:300; outcome=300:400; ITI=400:500
% 300 is when choice strobe is sent

%%
ofcin = load([dpath 'PS.StagOps.OFCin_atOFFER1.mat']);
ofcout = load([dpath 'PS.StagOps.OFCout_atOFFER1.mat']);
%% data

md = fitRegs(ofcin.data);
% md = fitRegs(ofcout.data);

close all; clc
% mdn = 1;
bev1 = md(1).bev1;
bev2 = md(1).bev2;
bevl = md(2).blev;
bevr = md(2).brev;

% 'Kendall'
% 'Spearman'
% 'Pearson'
[r,p] = corr(bev1,bev2,'Type','Spearman')
[r,p] = corr(bevl,bevr,'Type','Spearman')


%% jackknife - ofcin
for xx = 1:1000
    clear ofcinJK data md
    ofcinJK = randperm(size(ofcin.data,2),size(ofcin.data,2)-1);
    data = ofcin.data(:,ofcinJK);
    
    md = fitRegs(data);
    bev1 = md(1).bev1;
    bev2 = md(1).bev2;
    bevl = md(2).blev;
    bevr = md(2).brev;
    
    [r_inJK_ev12(xx),~] = corr(bev1,bev2,'Type','Spearman');
    [r_inJK_evlr(xx),~] = corr(bevl,bevr,'Type','Spearman');
    
    disp(['Done with OFCin JK repetition ' num2str(xx)]);
end

save([spath 'inJK.mat'],'r_inJK_ev12','r_inJK_evlr')

% jackknife - ofcout
for xx = 1:1000
    clear ofcoutJK data md
    ofcoutJK = randperm(size(ofcout.data,2),size(ofcout.data,2)-1);
    data = ofcout.data(:,ofcoutJK);
    
    md = fitRegs(data);
    bev1 = md(1).bev1;
    bev2 = md(1).bev2;
    bevl = md(2).blev;
    bevr = md(2).brev;
    
    [r_outJK_ev12(xx),~] = corr(bev1,bev2,'Type','Spearman');
    [r_outJK_evlr(xx),~] = corr(bevl,bevr,'Type','Spearman');
    
    disp(['Done with OFout JK repetition ' num2str(xx)]);
end
save([spath 'outJK.mat'],'r_outJK_ev12','r_outJK_evlr')


% 200-neuron - ofcin

for xx = 1:1000
    clear ofcin200 data md
    ofcin200 = randsample(size(ofcin.data,2),200,true);
    data = ofcin.data(:,ofcin200);
    
    md = fitRegs(data);
    bev1 = md(1).bev1;
    bev2 = md(1).bev2;
    bevl = md(2).blev;
    bevr = md(2).brev;
    
    [r_in200_ev12(xx),~] = corr(bev1,bev2,'Type','Spearman');
    [r_in200_evlr(xx),~] = corr(bevl,bevr,'Type','Spearman');
    
    disp(['Done with OFCin repetition ' num2str(xx)]);
end
save([spath 'in200.mat'],'r_in200_ev12','r_in200_evlr')

% 200-neuron - ofcout
for xx = 1:1000
    clear ofcout200 data md
    ofcout200 = randsample(size(ofcout.data,2),200,true);
    data = ofcout.data(:,ofcout200);
    
    md = fitRegs(data);
    bev1 = md(1).bev1;
    bev2 = md(1).bev2;
    bevl = md(2).blev;
    bevr = md(2).brev;
    
    [r_out200_ev12(xx),~] = corr(bev1,bev2,'Type','Spearman');
    [r_out200_evlr(xx),~] = corr(bevl,bevr,'Type','Spearman');
    
    disp(['Done with OFCout repetition ' num2str(xx)]);
end

save([spath 'out200.mat'],'r_out200_ev12','r_out200_evlr')
cd(spath)



%%

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







