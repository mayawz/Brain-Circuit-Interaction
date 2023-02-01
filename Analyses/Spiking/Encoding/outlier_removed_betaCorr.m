% OFC-PCC paper reply to reviewer analysis
% last used: 2021/10/14

% this addresses "PCCgurys vs. PCCsulcus"

clear all; close all; clc
% addpath(genpath('/'));

dpath = '/Volumes/GoogleDrive/My Drive/01Data/00 w_Anatomy/2020final/wrapped/PScombined/';
spath = '/Volumes/GoogleDrive/My Drive/01Data/00 w_Anatomy/2020final/R2R_result/';
fpath = '/Volumes/GoogleDrive/My Drive/01Data/00 w_Anatomy/2020final/R2R_code/';
addpath(genpath('/Volumes/GoogleDrive/My Drive/01Data/00 GeneralFunctions/'));
addpath(fpath);
cd(dpath);
% atOFFER1: opt1=301:400;opt2=401:500;
% atCHOICE: choice=200:300; outcome=300:400; ITI=400:500
% 300 is when choice strobe is sent

%%
tmp = load([dpath 'PS.StagOps.OFCin_atOFFER1.mat'],'data');
ofcin = tmp.data([1:40,42:end]);
clear tmp

tmp = load([dpath 'PS.StagOps.PCC_atOFFER1.mat'],'data');
pcc = tmp.data([1:117,119:end]);
clear tmp


%% PCC w/o outlier, EV1-EV2
% md = fitRegs(data,t,binSize)
md1 = fitRegs(pcc,436,20);

close all; clc

% 'Kendall'
% 'Spearman'
% 'Pearson'
[r1,p1] = corr(md1(1).bev1,md1(1).bev2,'Type','Spearman')

%%
% [zval, p]=FishersTransformation(r1, r2, n1, n2, tailOption, sameSample)
[zval, p]=FishersTransformation(r1, 0.02, size(pcc,2), 213, 'both', 0)

%% OFCin w/o outlier, EV1-EV2
% md = fitRegs(data,t,binSize)
md1 = fitRegs(ofcin,436,20);

close all; clc

% 'Kendall'
% 'Spearman'
% 'Pearson'
[r1,p1] = corr(md1(1).bev1,md1(1).bev2,'Type','Spearman')

%%
% [zval, p]=FishersTransformation(r1, r2, n1, n2, tailOption, sameSample)
[zval, p]=FishersTransformation(r1, -0.36, size(ofcin,2), 44, 'both', 0)




