% AnalyzeSlideRegress -- Proportion of cells
% last used: Mar 6 2020

% 1. proportion of cell
% 2. latency
% 3. encoding strength abs(t-stat): % variance explained by each factor in population

% 4. LSTM model of predictive relations between factors.

close all; clear all; clc
dpath='/Users/mayazwang/Google Drive/01Data/00 w_Anatomy/2020final/results/slideRegress/';
fpath='/Users/mayazwang/Google Drive/01Data/00 w_Anatomy/2020final/code/SlideRegress/';
spath='/Users/mayazwang/Google Drive/01Data/00 w_Anatomy/2020final/results/slideRegress//';
cd(dpath)

addpath(genpath(fpath))

% at OFFER
List=dir('*atOFFER*');      
List.name
OFFERin=load(List(1).name);
OFFERout=load(List(2).name);
OFFERpcc=load(List(3).name);

clear List

% at CHOICE
List=dir('*atCHOICE*');
List.name
CHOICEin=load(List(1).name);
CHOICEout=load(List(2).name);
CHOICEpcc=load(List(3).name);

%% ####################### OFFER VALUE #######################
%% OFCin proportion of cells tuning for EV1 during epoch 1 & 2
clc
bin=301:500;
[prop, peakBin]=max(mean(OFFERin.md(1).mdlPval(:,bin)<0.05))
prop*size(OFFERin.md(1).mdlPval,1)
[prop, peakBin]=max(mean(OFFERin.md(16).pev1(:,bin)<0.05))
prop*size(OFFERin.md(16).pev1,1)


%% OFCin proportion of cells tuning for EV2 during epoch 1 & 2
clc
[prop, peakBin]=max(mean(OFFERin.md(2).mdlPval(:,bin)<0.05))
prop*size(OFFERin.md(2).mdlPval,1)
[prop, peakBin]=max(mean(OFFERin.md(16).pev2(:,bin)<0.05))
prop*size(OFFERin.md(16).pev2,1)
%% OFCout proportion of cells tuning for EV1 during epoch 1 & 2
clc
[prop, peakBin]=max(mean(OFFERout.md(1).mdlPval(:,bin)<0.05))
prop*size(OFFERout.md(1).mdlPval,1)
[prop, peakBin]=max(mean(OFFERout.md(16).pev1(:,bin)<0.05))
prop*size(OFFERout.md(16).pev1,1)


%% OFCout proportion of cells tuning for EV2 during epoch 1 & 2
clc
[prop, peakBin]=max(mean(OFFERout.md(2).mdlPval(:,bin)<0.05))
prop*size(OFFERout.md(1).mdlPval,1)
[prop, peakBin]=max(mean(OFFERout.md(16).pev2(:,bin)<0.05))
prop*size(OFFERout.md(16).pev2,1)
%% PCC proportion of cells tuning for EV1 during epoch 1 & 2
clc
[prop, peakBin]=max(mean(OFFERpcc.md(1).mdlPval(:,bin)<0.05))
prop*size(OFFERpcc.md(1).mdlPval,1)
[prop, peakBin]=max(mean(OFFERpcc.md(16).pev1(:,bin)<0.05))
prop*size(OFFERpcc.md(16).pev1,1)

%% PCC
clc
[prop, peakBin]=max(mean(OFFERpcc.md(2).mdlPval(:,bin)<0.05))
prop*size(OFFERpcc.md(2).mdlPval,1)
[prop, peakBin]=max(mean(OFFERpcc.md(16).pev2(:,bin)<0.05))
prop*size(OFFERpcc.md(16).pev2,1)


%% EV1 --> wasn't super duper pretty
close all
subplot(1,2,1)
plot(smooth(mean(OFFERin.md(1).pev1(:,OFFERepoch)<0.05),smo),'b-')
hold on
plot(smooth(mean(OFFERout.md(1).pev1(:,OFFERepoch)<0.05),smo),'r-')
plot(smooth(mean(OFFERpcc.md(1).pev1(:,OFFERepoch)<0.05),smo),'k-')
vline(20)
vline(120)
vline(220)
hold off

subplot(1,2,2)
plot(smooth(mean(CHOICEin.md(1).pev1(:,CHOICEepoch)<0.05),smo),'b-')
hold on
plot(smooth(mean(CHOICEout.md(1).pev1(:,CHOICEepoch)<0.05),smo),'r-')
plot(smooth(mean(CHOICEpcc.md(1).pev1(:,CHOICEepoch)<0.05),smo),'k-')
vline(100)
vline(200)
hold off
%% ####################### REWARD VALUE #######################
%% OFCin during outcome epoch
clc
bin=301:400;
% mdn=11 12 17 18 19
[prop, peakBin]=max(mean(CHOICEin.md(11).prwdOtc(:,bin)<0.05))
prop*size(CHOICEin.md(11).prwdOtc,1)

%% OFCout during outcome epoch
clc
bin=301:400;
% mdn=11 12 17 18 19
[prop, peakBin]=max(mean(CHOICEout.md(11).prwdOtc(:,bin)<0.05))
prop*size(CHOICEout.md(11).prwdOtc,1)

%% PCC proportion of cells tuning for rewardOutcome during outcome epoch
clc
bin=301:400;
% mdn=11 12 17 18 19
[prop, peakBin]=max(mean(CHOICEpcc.md(11).prwdOtc(:,bin)<0.05))
prop*size(CHOICEpcc.md(11).prwdOtc,1)


%% ####################### CHOICE OPTION #######################
bin1=426:600;
bin2=151:325;
% model 4 11 12 17 18 19 
mdn=11;
%% OFCin 
close all; clc

[prop, peakBin]=max(mean(OFFERin.md(mdn).pchoOpt(:,bin1)<0.05))
prop*size(OFFERin.md(mdn).pchoOpt,1)

[prop, peakBin]=max(mean(CHOICEin.md(mdn).pchoOpt(:,bin2)<0.05))
prop*size(CHOICEin.md(mdn).pchoOpt,1)
% 2273; 14/16; 20/14; 14/25; 14/25; 16/14

plot(mean(OFFERin.md(mdn).pchoOpt(:,bin1)<0.05),'r-')
hold on
plot(mean(CHOICEin.md(mdn).pchoOpt(:,bin2)<0.05),'b-')
%% OFCout 
close all; clc

[prop, peakBin]=max(mean(OFFERout.md(mdn).pchoOpt(:,bin1)<0.05))
prop*size(OFFERout.md(mdn).pchoOpt,1)

[prop, peakBin]=max(mean(CHOICEout.md(mdn).pchoOpt(:,bin2)<0.05))
prop*size(CHOICEout.md(mdn).pchoOpt,1)
% 2273; 14/16; 20/14; 14/25; 14/25; 16/14
% bins: offer 448 = choice 115???

plot(mean(OFFERout.md(mdn).pchoOpt(:,bin1)<0.05),'r-')
hold on
plot(mean(CHOICEout.md(mdn).pchoOpt(:,bin2)<0.05),'b-')

%% PCC 
close all; clc

[prop, peakBin]=max(mean(OFFERpcc.md(mdn).pchoOpt(:,bin1)<0.05))
prop*size(OFFERpcc.md(mdn).pchoOpt,1)

[prop, peakBin]=max(mean(CHOICEpcc.md(mdn).pchoOpt(:,bin2)<0.05))
prop*size(CHOICEpcc.md(mdn).pchoOpt,1)
% 2273; 14/16; 20/14; 14/25; 14/25; 16/14
% bins: offer 448 = choice 115???

plot(mean(OFFERpcc.md(mdn).pchoOpt(:,bin1)<0.05),'r-')
hold on
plot(mean(CHOICEpcc.md(mdn).pchoOpt(:,bin2)<0.05),'b-')

%% ####################### CHOICE LR #######################
bin1=420:600;
bin2=120:300;
% model 3 11 12 17 18 19 
mdn=11;

%% OFCin 
close all; clc

[prop, peakBin]=max(mean(OFFERin.md(mdn).pchoLR(:,bin1)<0.05))
prop*size(OFFERin.md(mdn).pchoLR,1)

[prop, peakBin]=max(mean(CHOICEin.md(mdn).pchoLR(:,bin2)<0.05))
prop*size(CHOICEin.md(mdn).pchoLR,1)
% 2273; 14/16; 20/14; 14/25; 14/25; 16/14
% bins: offer 448 = choice 115???

plot(mean(OFFERin.md(mdn).pchoLR(:,bin1)<0.05),'r-')
hold on
plot(mean(CHOICEin.md(mdn).pchoLR(:,bin2)<0.05),'b-')
%% OFCout 
close all; clc

[prop, peakBin]=max(mean(OFFERout.md(mdn).pchoLR(:,bin1)<0.05))
prop*size(OFFERout.md(mdn).pchoLR,1)

[prop, peakBin]=max(mean(CHOICEout.md(mdn).pchoLR(:,bin2)<0.05))
prop*size(CHOICEout.md(mdn).pchoLR,1)
% 2273; 14/16; 20/14; 14/25; 14/25; 16/14
% bins: offer 448 = choice 115???

plot(mean(OFFERout.md(mdn).pchoLR(:,bin1)<0.05),'r-')
hold on
plot(mean(CHOICEout.md(mdn).pchoLR(:,bin2)<0.05),'b-')

%% PCC 
close all; clc

[prop, peakBin]=max(mean(OFFERpcc.md(mdn).pchoLR(:,bin1)<0.05))
prop*size(OFFERpcc.md(mdn).pchoLR,1)

[prop, peakBin]=max(mean(CHOICEpcc.md(mdn).pchoLR(:,bin2)<0.05))
prop*size(CHOICEpcc.md(mdn).pchoLR,1)
% 2273; 14/16; 20/14; 14/25; 14/25; 16/14
% bins: offer 448 = choice 115???

plot(mean(OFFERpcc.md(mdn).pchoLR(:,bin1)<0.05),'r-')
hold on
plot(mean(CHOICEpcc.md(mdn).pchoLR(:,bin2)<0.05),'b-')









