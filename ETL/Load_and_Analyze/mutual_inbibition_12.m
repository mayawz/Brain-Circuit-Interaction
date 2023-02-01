% % AnalyzeSlideRegress -- mutual inbibition offer 12
% last used: Mar 7 2020

close all; clear all; clc
dpath='/Volumes/GoogleDrive/My Drive/01Data/00 w_Anatomy/2020final/results/slideRegress/';
fpath='/Volumes/GoogleDrive/My Drive/01Data/00 w_Anatomy/2020final/results/slideRegress/';
spath='/Volumes/GoogleDrive/My Drive/01Data/00 w_Anatomy/2020final/results/slideRegress/Mutual_Inhibition_Offer12/';
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


%% save betas for outlier test
mdN = 16;
binN = 401+36-1;
clear betas
betas_in = [OFFERin.md(mdN).bev1(:,binN),...
    OFFERin.md(mdN).bev2(:,binN),...
    OFFERin.md(17).blev(:,453),...
    OFFERin.md(17).brev(:,453)];
betas_out = [OFFERout.md(mdN).bev1(:,binN),...
    OFFERout.md(mdN).bev2(:,binN),...
    OFFERout.md(17).blev(:,453),...
    OFFERout.md(17).brev(:,453)];
betas_pcc = [OFFERpcc.md(mdN).bev1(:,binN),...
    OFFERpcc.md(mdN).bev2(:,binN),...
    OFFERpcc.md(17).blev(:,453),...
    OFFERpcc.md(17).brev(:,453)];

writematrix(betas_in,'betas_in.csv', 'FileType', 'text');
writematrix(betas_out,'betas_out.csv', 'FileType', 'text');
writematrix(betas_pcc,'betas_pcc.csv', 'FileType', 'text');


%% ####################### 	OFCin  #######################
clc; close all
close all; clc
sBin=401:550;
mdn=16;
bev1=OFFERin.md(mdn).bev1(:,sBin);
bev2=OFFERin.md(mdn).bev2(:,sBin);
for tt=1:length(sBin)
    [rO(tt),pO(tt)] = corr(bev1(:,tt),bev2(:,tt),'Type','Spearman');
end
% 'Kendall'
% 'Spearman'
% 'Pearson'
sBinC=151:300;

bev1C=CHOICEin.md(mdn).bev1(:,sBinC);
bev2C=CHOICEin.md(mdn).bev2(:,sBinC);
for tt=1:length(sBinC)
    [rC(tt),pC(tt)] = corr(bev1C(:,tt),bev2C(:,tt),'Type','Spearman');
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
In_rEV1EV2=rO;

save([spath 'In_EV1EV2.mat'],'In_rEV1EV2')


%%
[val,indx]=min(rO)
sigslides=find(pO<0.05)

%% OFFER
close all; clc
% indx=36;
sBin=401:550;
bev1=OFFERin.md(mdn).bev1(:,sBin);
bev2=OFFERin.md(mdn).bev2(:,sBin);

[r,p] = corr(bev1(:,indx),bev2(:,indx),'Type','Spearman')

ax_x=1;
ax_y=1;
figBetaCorr(bev1(:,indx),bev2(:,indx),nan,nan,r,ax_x,ax_y)


%% CHOICE
close all; clc
indx=110;
sBinC=151:300;
mdn
bev1C=CHOICEin.md(mdn).bev1(:,sBinC);
bev2C=CHOICEin.md(mdn).bev2(:,sBinC);

[r,p] = corr(bev1C(:,indx),bev2C(:,indx),'Type','Spearman')

ax_x=1;
ax_y=1;
figBetaCorr(bev1C(:,indx),bev2C(:,indx),nan,nan,r,ax_x,ax_y)


%% ####################### 	OFCout  #######################
clc; close all
close all; clc
sBin=401:550;
mdn=16;
bev1=OFFERout.md(mdn).bev1(:,sBin);
bev2=OFFERout.md(mdn).bev2(:,sBin);
for tt=1:length(sBin)
    [rO(tt),pO(tt)] = corr(bev1(:,tt),bev2(:,tt),'Type','Spearman');
end
% 'Kendall'
% 'Spearman'
% 'Pearson'
sBinC=151:300;

bev1C=CHOICEout.md(mdn).bev1(:,sBinC);
bev2C=CHOICEout.md(mdn).bev2(:,sBinC);
for tt=1:length(sBinC)
    [rC(tt),pC(tt)] = corr(bev1C(:,tt),bev2C(:,tt),'Type','Spearman');
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
Out_rEV1EV2=rO;

save([spath 'Out_rEV1EV2.mat'],'Out_rEV1EV2')

%%
close all; clc
indx=36;
sBin=401:550;
bev1=OFFERout.md(mdn).bev1(:,sBin);
bev2=OFFERout.md(mdn).bev2(:,sBin);

[r,p] = corr(bev1(:,indx),bev2(:,indx),'Type','Spearman')

ax_x=1;
ax_y=1;
figBetaCorr(bev1(:,indx),bev2(:,indx),nan,nan,r,ax_x,ax_y)

%% CHOICE
close all; clc
indx=110;
sBinC=151:300;
mdn
bev1C=CHOICEout.md(mdn).bev1(:,sBinC);
bev2C=CHOICEout.md(mdn).bev2(:,sBinC);

[r,p] = corr(bev1C(:,indx),bev2C(:,indx),'Type','Spearman')

ax_x=1;
ax_y=1;
figBetaCorr(bev1C(:,indx),bev2C(:,indx),nan,nan,r,ax_x,ax_y)

%% ####################### 	PCC  #######################
clc; close all
close all; clc
sBin=401:550;
mdn=16;
bev1=OFFERpcc.md(mdn).bev1(:,sBin);
bev2=OFFERpcc.md(mdn).bev2(:,sBin);
for tt=1:length(sBin)
    [rO(tt),pO(tt)] = corr(bev1(:,tt),bev2(:,tt),'Type','Spearman');
end
% 'Kendall'
% 'Spearman'
% 'Pearson'
sBinC=151:300;

bev1C=CHOICEpcc.md(mdn).bev1(:,sBinC);
bev2C=CHOICEpcc.md(mdn).bev2(:,sBinC);
for tt=1:length(sBinC)
    [rC(tt),pC(tt)] = corr(bev1C(:,tt),bev2C(:,tt),'Type','Spearman');
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
Pcc_rEV1EV2=rO;

save([spath 'Pcc_rEV1EV2.mat'],'Pcc_rEV1EV2')
%%
close all; clc
indx=36;
sBin=401:550;
bev1=OFFERpcc.md(mdn).bev1(:,sBin);
bev2=OFFERpcc.md(mdn).bev2(:,sBin);

[r,p] = corr(bev1(:,indx),bev2(:,indx),'Type','Spearman')

ax_x=1;
ax_y=1;
figBetaCorr(bev1(:,indx),bev2(:,indx),nan,nan,r,ax_x,ax_y)


%% CHOICE
close all; clc
indx=110;
sBinC=151:300;
mdn
bev1C=CHOICEpcc.md(mdn).bev1(:,sBinC);
bev2C=CHOICEpcc.md(mdn).bev2(:,sBinC);

[r,p] = corr(bev1C(:,indx),bev2C(:,indx),'Type','Spearman')

ax_x=1;
ax_y=1;
figBetaCorr(bev1C(:,indx),bev2C(:,indx),nan,nan,r,ax_x,ax_y)




