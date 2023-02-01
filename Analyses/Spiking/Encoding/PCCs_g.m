% OFC-PCC paper reply to reviewer analysis
% last used: 2021/10/14

% this addresses "PCCgurys vs. PCCsulcus"

clear all; close all; clc
% addpath(genpath('/'));

dpath = '/Users/wangz37/OneDrive - National Institutes of Health/_paper_submission/OFC_PCC/2021_R2R/wrapped/PScombined/';
spath = '/Users/wangz37/OneDrive - National Institutes of Health/_paper_submission/OFC_PCC/2021_R2R/R2R_result/';
fpath = '/Users/wangz37/OneDrive - National Institutes of Health/_paper_submission/OFC_PCC/2021_R2R/R2R_code/';
addpath(genpath('/Users/wangz37/OneDrive - National Institutes of Health/_paper_submission/OFC_PCC/2021_R2R/00 GeneralFunctions/'));
addpath(fpath);
cd(dpath);
% atOFFER1: opt1=301:400;opt2=401:500;
% atCHOICE: choice=200:300; outcome=300:400; ITI=400:500
% 300 is when choice strobe is sent

%%
tmp = load([dpath 'PCCg_atOFFER1.mat'],'data');
PCCg_off = tmp.data;
clear tmp

tmp = load([dpath 'PCCg_atCHOICE.mat'],'data');
PCCg_cho = tmp.data;
clear tmp

tmp = load([dpath 'PCCs_atOFFER1.mat'],'data');
PCCs_off = tmp.data;
clear tmp

tmp = load([dpath 'PCCs_atCHOICE.mat'],'data');
PCCs_cho = tmp.data;
clear tmp

%% PCCg vs. PCCs, EV1-EV2
% md = fitRegs(data,t,binSize)
md1 = fitRegs(PCCg_off,435,25);
md2 = fitRegs(PCCs_off,435,25);

close all; clc
% mdn = 1;
bev1 = md(1).bev1;
bev2 = md(1).bev2;


% 'Kendall'
% 'Spearman'
% 'Pearson'
[r1,p1] = corr(md1(1).bev1,md1(1).bev2,'Type','Spearman')
[r2,p2] = corr(md2(1).bev1,md2(1).bev2,'Type','Spearman')

%%
% [zval, p]=FishersTransformation(r1, r2, n1, n2, tailOption, sameSample)
[zval, p]=FishersTransformation(r1, r2, size(PCCg_off,2), size(PCCs_off,2), 'both', 0)
[zval, p]=FishersTransformation(r1, 0.02, size(PCCg_off,2), size(PCCg_off,2)+size(PCCs_off,2), 'both', 0)
[zval, p]=FishersTransformation(r2, 0.02, size(PCCs_off,2), size(PCCg_off,2)+size(PCCs_off,2), 'both', 0)
%% PCCg vs. PCCs, EVl-EVr
% md = fitRegs(data,t,binSize)
% md3 = fitRegs(PCCg_cho,150,50);
% md4 = fitRegs(PCCs_cho,200,50);
% 
% md3 = fitRegs(PCCg_cho,330,50);
% md4 = fitRegs(PCCs_cho,340,50);

md3 = fitRegs(PCCg_off,453,40); % matching original code

md4 = fitRegs(PCCs_off,453,40);


close all; clc

% 'Kendall'
% 'Spearman'
% 'Pearson'
[r3,p3] = corr(md3(2).blev,md3(2).brev,'Type','Spearman')
[r4,p4] = corr(md4(2).blev,md4(2).brev,'Type','Spearman')

%%
[zval, p]=FishersTransformation(r3, r4, size(PCCg_cho,2), size(PCCs_cho,2), 'both', 0)
[zval, p]=FishersTransformation(r3, -0.24, size(PCCg_cho,2), size(PCCg_cho,2)+size(PCCs_cho,2), 'both', 0)
[zval, p]=FishersTransformation(r4, -0.24, size(PCCs_cho,2), size(PCCg_cho,2)+size(PCCs_cho,2), 'both', 0)
