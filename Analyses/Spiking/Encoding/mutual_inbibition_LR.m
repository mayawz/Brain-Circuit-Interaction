% % AnalyzeSlideRegress -- mutual inbibition offer LR
% last used: Mar 7 2020

close all; clear all; clc
dpath='/Users/mwang/Google Drive/01Data/00 w_Anatomy/2020final/results/slideRegress/';
fpath='/Users/mwang/Google Drive/01Data/00 w_Anatomy/2020final/code/SlideRegress/';
spath='/Users/mwang/Google Drive/01Data/00 w_Anatomy/2020final/results/slideRegress/Mutual_Inhibition_OfferLR/';
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

smo=1;

%% ####################### 	OFCin  #######################
clc; close all
close all; clc
sBin=401:550;
mdn=17;
blev=OFFERin.md(mdn).blev(:,sBin);
brev=OFFERin.md(mdn).brev(:,sBin);
for tt=1:length(sBin)
    [rO(tt),pO(tt)] = corr(blev(:,tt),brev(:,tt),'Type','Spearman');
end
% 'Kendall'
% 'Spearman'
% 'Pearson'
sBinC=151:300;

blevC=CHOICEin.md(mdn).blev(:,sBinC);
brevC=CHOICEin.md(mdn).brev(:,sBinC);
for tt=1:length(sBinC)
    [rC(tt),pC(tt)] = corr(blevC(:,tt),brevC(:,tt),'Type','Spearman');
end


subplot(1,2,1)
plot(smooth(rO,smo))
hold on
plot(find(pO<0.05),zeros(1,length(find(pO<0.05))),'ro')
hline(0)
hold off

subplot(1,2,2)
plot(smooth(rC,smo))
hold on
plot(find(pC<0.05),zeros(1,length(find(pC<0.05))),'ro')
hline(0)
hold off

%% save for cross correlation
In_rEVLEVR=rO;

save([spath 'In_rEVLEVR.mat'],'In_rEVLEVR')


%%
[val,indx]=min(rO)
sigslides=find(pO<0.05)
% 116   117   118   119   123   125   126   127   128   129
%% OFFER
close all; clc
% indx=53;
indx=116
sBin=401:550;
blev=OFFERin.md(mdn).blev(:,sBin);
brev=OFFERin.md(mdn).brev(:,sBin);

[r,p] = corr(blev(:,indx),brev(:,indx),'Type','Spearman')

ax_x=1;
ax_y=1;
figBetaCorr(blev(:,indx),brev(:,indx),nan,nan,r,ax_x,ax_y)

% @ 53 r = -0.1617 p =0.2934
% @ 127 r = -0.3869 p =0.0099
%% CHOICE
close all; clc
indx=110;
sBinC=151:300;
mdn
blevC=CHOICEin.md(mdn).blev(:,sBinC);
brevC=CHOICEin.md(mdn).brev(:,sBinC);

[r,p] = corr(blevC(:,indx),brevC(:,indx),'Type','Spearman')

ax_x=1;
ax_y=1;
figBetaCorr(blevC(:,indx),brevC(:,indx),nan,nan,r,ax_x,ax_y)
% r = 0.0032 p =0.9836

%% ####################### 	OFCout  #######################
clc; close all
close all; clc
sBin=401:550;
mdn=17;
blev=OFFERout.md(mdn).blev(:,sBin);
brev=OFFERout.md(mdn).brev(:,sBin);
for tt=1:length(sBin)
    [rO(tt),pO(tt)] = corr(blev(:,tt),brev(:,tt),'Type','Spearman');
end
% 'Kendall'
% 'Spearman'
% 'Pearson'
sBinC=151:300;

blevC=CHOICEout.md(mdn).blev(:,sBinC);
brevC=CHOICEout.md(mdn).brev(:,sBinC);
for tt=1:length(sBinC)
    [rC(tt),pC(tt)] = corr(blevC(:,tt),brevC(:,tt),'Type','Spearman');
end


subplot(1,2,1)
plot(smooth(rO,smo))
hold on
plot(find(pO<0.05),zeros(1,length(find(pO<0.05))),'ro')
hline(0)
hold off

subplot(1,2,2)
plot(smooth(rC,smo))
hold on
plot(find(pC<0.05),zeros(1,length(find(pC<0.05))),'ro')
hline(0)
hold off

%% save for cross correlation
Out_rEVLEVR=rO;

save([spath 'Out_rEVLEVR.mat'],'Out_rEVLEVR')

%%
[val,indx]=min(rO)
sigslides=find(pO<0.05)

% 1    64    74    76   142
%%
close all; clc
% indx=53;
indx=74
sBin=401:550;
blev=OFFERout.md(mdn).blev(:,sBin);
brev=OFFERout.md(mdn).brev(:,sBin);

[r,p] = corr(blev(:,indx),brev(:,indx),'Type','Spearman')

ax_x=1;
ax_y=1;
figBetaCorr(blev(:,indx),brev(:,indx),nan,nan,r,ax_x,ax_y)

% @53 r = 0.1787 p =0.1954
% @53 r = 0.1787 p =0.1954
%% CHOICE
close all; clc
indx=110;
sBinC=151:300;
mdn
blevC=CHOICEout.md(mdn).blev(:,sBinC);
brevC=CHOICEout.md(mdn).brev(:,sBinC);

[r,p] = corr(blevC(:,indx),brevC(:,indx),'Type','Spearman')

ax_x=1;
ax_y=1;
figBetaCorr(blevC(:,indx),brevC(:,indx),nan,nan,r,ax_x,ax_y)
% r = -0.1427 p =0.3022

%% ####################### 	PCC  #######################
clc; close all
close all; clc
sBin=401:550;
mdn=17;
blev=OFFERpcc.md(mdn).blev(:,sBin);
brev=OFFERpcc.md(mdn).brev(:,sBin);
for tt=1:length(sBin)
    [rO(tt),pO(tt)] = corr(blev(:,tt),brev(:,tt),'Type','Spearman');
end
% 'Kendall'
% 'Spearman'
% 'Pearson'
sBinC=151:300;

blevC=CHOICEpcc.md(mdn).blev(:,sBinC);
brevC=CHOICEpcc.md(mdn).brev(:,sBinC);
for tt=1:length(sBinC)
    [rC(tt),pC(tt)] = corr(blevC(:,tt),brevC(:,tt),'Type','Spearman');
end


subplot(1,2,1)
plot(smooth(rO,smo))
hold on
% ylim([-0.5 0.3])
plot(find(pO<0.05),zeros(1,length(find(pO<0.05))),'ro')
hline(0)
hold off

subplot(1,2,2)
plot(smooth(rC,smo))
hold on
% ylim([-0.5 0.3])
plot(find(pC<0.05),zeros(1,length(find(pC<0.05))),'ro')
hline(0)
hold off
%% save for cross correlation
Pcc_rEVLEVR=rO;

save([spath 'Pcc_rEVLEVR.mat'],'Pcc_rEVLEVR')
%%
[val,indx]=min(rO)
sigslides=find(pO<0.05)

%     53     65     73    74    75  
%%
close all; clc
% indx=53;
indx=75;
sBin=401:550;
mdn
blev=OFFERpcc.md(mdn).blev(:,sBin);
brev=OFFERpcc.md(mdn).brev(:,sBin);

[r,p] = corr(blev(:,indx),brev(:,indx),'Type','Spearman')

ax_x=1;
ax_y=1;
figBetaCorr(blev(:,indx),brev(:,indx),nan,nan,r,ax_x,ax_y)
% @53:r = -0.2360  p = 5.3230e-04
% @75: r = -0.1864  p = 0.0064
%% CHOICE
close all; clc
indx=110;
sBinC=151:300;
mdn
blevC=CHOICEpcc.md(mdn).blev(:,sBinC);
brevC=CHOICEpcc.md(mdn).brev(:,sBinC);

[r,p] = corr(blevC(:,indx),brevC(:,indx),'Type','Spearman')

ax_x=1;
ax_y=1;
figBetaCorr(blevC(:,indx),brevC(:,indx),nan,nan,r,ax_x,ax_y)
% r = 0.0023 p =0.9731




