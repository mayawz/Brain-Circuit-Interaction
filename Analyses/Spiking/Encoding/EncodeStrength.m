% AnalyzeSlideRegress -- encoding strength
% last used: Mar 6 2020

% 1. proportion of cell
% 2. latency
% 3. encoding strength abs(t-stat): % variance explained by each factor in population

% 4. LSTM model of predictive relations between factors.

close all; clear all; clc
dpath='/Users/mwang/Google Drive/01Data/00 w_Anatomy/2020final/results/slideRegress/';
fpath='/Users/mwang/Google Drive/01Data/00 w_Anatomy/2020final/code/SlideRegress/';
spath='/Users/mwang/Google Drive/01Data/00 w_Anatomy/2020final/results/slideRegress/crossPredict/';
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

%% Offer1
clc
bin=301:500;
stgPCC=nanmean(abs(OFFERpcc.md(16).tev1(:,bin)),2);
stgIn=nanmean(abs(OFFERin.md(16).tev1(:,bin)),2);
stgOut=nanmean(abs(OFFERout.md(16).tev1(:,bin)),2);
stgOut(find(isinf(stgOut)==1))=0;

X=nan(length(stgPCC),3);
X(1:length(stgIn),1)=stgIn;
X(1:length(stgOut),2)=stgOut;
X(1:length(stgPCC),3)=stgPCC;

[p,~,stats] = kruskalwallis(X)

%% Reward
close all; clc
bin=301:450;
stgPCC=nanmean(abs(CHOICEpcc.md(11).trwdOtc(:,bin)),2);
stgIn=nanmean(abs(CHOICEin.md(11).trwdOtc(:,bin)),2);
stgOut=nanmean(abs(CHOICEout.md(11).trwdOtc(:,bin)),2);
% stgOut(find(isinf(stgOut)==1))=0;

X=nan(length(stgPCC),3);
X(1:length(stgIn),1)=stgIn;
X(1:length(stgOut),2)=stgOut;
X(1:length(stgPCC),3)=stgPCC;

[p,~,stats] = kruskalwallis(X)

multcompare(stats,'CType','tukey-kramer')

%% Choice Option -- use OFFER
close all; clc
bin=401:550;
stgPCC=nanmean(abs(OFFERpcc.md(11).tchoOpt(:,bin)),2);
stgIn=nanmean(abs(OFFERin.md(11).tchoOpt(:,bin)),2);
stgOut=nanmean(abs(OFFERout.md(11).tchoOpt(:,bin)),2);
% stgOut(find(isinf(stgOut)==1))=0;

X=nan(length(stgPCC),3);
X(1:length(stgIn),1)=stgIn;
X(1:length(stgOut),2)=stgOut;
X(1:length(stgPCC),3)=stgPCC;

[p,~,stats] = kruskalwallis(X)

% multcompare(stats,'CType','tukey-kramer')

%% save for cross prediction
varChoOpt.in=stgIn;
varChoOpt.out=stgOut;
varChoOpt.pcc=stgPCC;

varChoOpt.in_t=nanmean(abs(OFFERin.md(11).tchoOpt(:,bin)),1);
varChoOpt.out_t=nanmean(abs(OFFERout.md(11).tchoOpt(:,bin)),1);
varChoOpt.pcc_t=nanmean(abs(OFFERpcc.md(11).tchoOpt(:,bin)),1);

save([spath 'varChoOpt.mat'], 'varChoOpt')
cd(spath)
%% Choice Location -- use OFFER

close all; clc
bin=401:550;
stgPCC=nanmean(abs(OFFERpcc.md(11).tchoLR(:,bin)),2);
stgIn=nanmean(abs(OFFERin.md(11).tchoLR(:,bin)),2);
stgOut=nanmean(abs(OFFERout.md(11).tchoLR(:,bin)),2);
% stgOut(find(isinf(stgOut)==1))=0;

X=nan(length(stgPCC),3);
X(1:length(stgIn),1)=stgIn;
X(1:length(stgOut),2)=stgOut;
X(1:length(stgPCC),3)=stgPCC;

[p,~,stats] = kruskalwallis(X)

% multcompare(stats,'CType','tukey-kramer')

%% save for cross prediction
varChoLR.in=stgIn;
varChoLR.out=stgOut;
varChoLR.pcc=stgPCC;

varChoLR.in_t=nanmean(abs(OFFERin.md(11).tchoLR(:,bin)),1);
varChoLR.out_t=nanmean(abs(OFFERout.md(11).tchoLR(:,bin)),1);
varChoLR.pcc_t=nanmean(abs(OFFERpcc.md(11).tchoLR(:,bin)),1);

save([spath 'varChoLR.mat'], 'varChoLR')
cd(spath)

%% Choice Option -- use CHOICE

close all; clc
bin=151:300;

stgPCC=nanmean(abs(CHOICEpcc.md(11).tchoOpt(:,bin)),2);
stgIn=nanmean(abs(CHOICEin.md(11).tchoOpt(:,bin)),2);
stgOut=nanmean(abs(CHOICEout.md(11).tchoOpt(:,bin)),2);
% stgOut(find(isinf(stgOut)==1))=0;

X=nan(length(stgPCC),3);
X(1:length(stgIn),1)=stgIn;
X(1:length(stgOut),2)=stgOut;
X(1:length(stgPCC),3)=stgPCC;

[p,~,stats] = kruskalwallis(X)

% multcompare(stats,'CType','tukey-kramer')

%% Choice Option -- use CHOICE

close all; clc
bin=151:300;

stgPCC=nanmean(abs(CHOICEpcc.md(11).tchoLR(:,bin)),2);
stgIn=nanmean(abs(CHOICEin.md(11).tchoLR(:,bin)),2);
stgOut=nanmean(abs(CHOICEout.md(11).tchoLR(:,bin)),2);
% stgOut(find(isinf(stgOut)==1))=0;

X=nan(length(stgPCC),3);
X(1:length(stgIn),1)=stgIn;
X(1:length(stgOut),2)=stgOut;
X(1:length(stgPCC),3)=stgPCC;

[p,~,stats] = kruskalwallis(X)

% multcompare(stats,'CType','tukey-kramer')






