% AnalyzeSlideRegress -- Latency
% last used: Mar 6 2020



close all; clear all; clc
dpath='/Users/mayazwang/Google Drive/01Data/00 w_Anatomy/2020final/results/slideRegress/';
% dpath='/Users/mayazwang/Google Drive/01Data/00 w_Anatomy/2020final/results/slideRegress/Nov2019SlideRegress/';
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

%% ####################### CHOICE LR #######################
% sorted peak t-stat plot
close all; clc

mdn=11;

bin1=415:550;
bin2=251:400;

vline1=80;
vline2=70;
colorBarLim=[0,6];

subplot(2,3,1)
tmp=abs(OFFERin.md(mdn).tchoLR(:,bin1));
[~,IND]=max(tmp,[],2);
[~,ind] = sort(IND);
imagesc(tmp(ind,:))
colormap hot
colorbar
caxis(colorBarLim)
hold on
vline(vline1)
title('OFCin - Offer')
clear tmp IND ind

subplot(2,3,2)
tmp=abs(OFFERout.md(mdn).tchoLR(:,bin1));
[~,IND]=max(tmp,[],2);
[~,ind] = sort(IND);
imagesc(tmp(ind,:))
colormap hot
colormap hot
colorbar
caxis(colorBarLim)
hold on
vline(vline1)
title('OFCout - Offer')
clear tmp IND ind

subplot(2,3,3)
tmp=abs(OFFERpcc.md(mdn).tchoLR(:,bin1));
[~,IND]=max(tmp,[],2);
[~,ind] = sort(IND);
imagesc(tmp(ind,:))
colormap hot
colorbar
caxis(colorBarLim)
hold on
vline(vline1)
title('PCC - Offer')
clear tmp IND ind

subplot(2,3,4)
tmp=abs(CHOICEin.md(mdn).tchoLR(:,bin2));
[~,IND]=max(tmp,[],2);
[~,ind] = sort(IND);
imagesc(tmp(ind,:))
colormap hot
colorbar
caxis(colorBarLim)
hold on
vline(vline2)
title('OFCin_Choice')
clear tmp IND ind

subplot(2,3,5)
tmp=abs(CHOICEout.md(mdn).tchoLR(:,bin2));
[~,IND]=max(tmp,[],2);
[~,ind] = sort(IND);
imagesc(tmp(ind,:))
colormap hot
colorbar
caxis(colorBarLim)
hold on
vline(vline2)
title('OFCout - Choice')
clear tmp IND ind

subplot(2,3,6)
tmp=abs(CHOICEpcc.md(mdn).tchoLR(:,bin2));
[~,IND]=max(tmp,[],2);
[~,ind] = sort(IND);
imagesc(tmp(ind,:))
colormap hot
colorbar
caxis(colorBarLim)
hold on
vline(vline2)
title('PCC - Choice')
clear tmp IND ind
% conclusion: 530 in OFFER =?= 300 in CHOICE
%% sliding procedure
close all; clc 
Bwidth=30;
starts=421:5:551-Bwidth;
mdn=11;

for tt=1:length(starts)
    %%
    tt=5
    sBin=starts(tt):starts(tt)-1+Bwidth;
    tmp1=abs(OFFERpcc.md(mdn).tchoLR(:,sBin));
    tmp2=abs(OFFERin.md(mdn).tchoLR(:,sBin));
    tmp3=abs(OFFERout.md(mdn).tchoLR(:,sBin));
    [~, latPCC]=max(tmp1,[],2);
    [~, latIn]=max(tmp2,[],2);
    [~, latOut]=max(tmp3,[],2);
    grpVar=zeros(length(latPCC)+length(latIn)+length(latOut),1);
    grpVar(1:length(latPCC),1)=1;
    grpVar(length(latPCC)+1:length(latPCC)+length(latIn),1)=2;
    grpVar(length(latPCC)+length(latIn)+1:length(latPCC)+length(latIn)+length(latOut),1)=3;
    
    yy=[latPCC;latIn;latOut];
    xx=categorical(grpVar);
    mdl = fitglm(xx,yy,'Distribution','gamma');
    mdl
    %%
    p=double(mdl.devianceTest.pValue);
    mdltsig(tt)=p(2);
    clear mdl
    
    %%
    
    sigPCC=OFFERpcc.md(mdn).pchoLR(:,sBin)<0.05;
    sigIn=OFFERin.md(mdn).pchoLR(:,sBin)<0.05;
    sigOut=OFFERout.md(mdn).pchoLR(:,sBin)<0.05;
    % find first sig in a row
    [~, latPCC]=max(sigPCC,[],2);
    [~, latIn]=max(sigIn,[],2);
    [~, latOut]=max(sigOut,[],2);
    grpVar=zeros(length(latPCC)+length(latIn)+length(latOut),1);
    grpVar(1:length(latPCC),1)=1;
    grpVar(length(latPCC)+1:length(latPCC)+length(latIn),1)=2;
    grpVar(length(latPCC)+length(latIn)+1:length(latPCC)+length(latIn)+length(latOut),1)=3;
    
    yy=[latPCC;latIn;latOut];
    xx=categorical(grpVar);
    mdl = fitglm(xx,yy,'Distribution','gamma');
    p=double(mdl.devianceTest.pValue);
    mdlpsig(tt)=p(2);
    clear mdl
end
%%
clc
find(mdltsig<0.05)
useTStart=starts(mdltsig<0.05)
find(mdlpsig<0.05)
usePStart=starts(mdlpsig<0.05)
% used mdltsig (5)= 441

%%
clc
inPCC=median(latPCC) % 14
inIn=median(latIn) % 15
inOut=median(latOut) % 23
% 441-420+23
%% compare 2
clc
grpVar=zeros(length(latPCC)+length(latIn)+length(latOut),1);
grpVar(1:length(latOut),1)=1;
grpVar(length(latOut)+1:length(latOut)+length(latIn),1)=2;
grpVar(length(latOut)+length(latIn)+1:length(latOut)+length(latIn)+length(latPCC),1)=3;
yy=[latOut;latIn;latPCC];
xx=categorical(grpVar);
mdl = fitglm(xx,yy,'Distribution','gamma')

%% compare 3
clc
grpVar=zeros(length(latPCC)+length(latIn)+length(latOut),1);
grpVar(1:length(latIn),1)=1;
grpVar(length(latIn)+1:length(latIn)+length(latOut),1)=2;
grpVar(length(latIn)+length(latOut)+1:length(latIn)+length(latOut)+length(latPCC),1)=3;
yy=[latIn;latOut;latPCC];
xx=categorical(grpVar);
mdl = fitglm(xx,yy,'Distribution','gamma')

%% Kruskal-Wallis
X=nan(length(latPCC),3);
X(1:length(latIn),1)=latIn;
X(1:length(latOut),2)=latOut;
X(1:length(latPCC),3)=latPCC;

[p,~,stats] = kruskalwallis(X)
multcompare(stats,'CType','tukey-kramer')

%% ####################### CHOICE 1 2 #######################
% sorted peak t-stat plot
% sorted peak t-stat plot
close all; clc

mdn=11;

bin1=415:550;
bin2=251:400;

vline1=80;
vline2=70;
colorBarLim=[0,6];

subplot(2,3,1)
tmp=abs(OFFERin.md(mdn).tchoOpt(:,bin1));
[~,IND]=max(tmp,[],2);
[~,ind] = sort(IND);
imagesc(tmp(ind,:))
colormap hot
colorbar
caxis(colorBarLim)
hold on
vline(vline1)
title('OFCin - Offer')
clear tmp IND ind

subplot(2,3,2)
tmp=abs(OFFERout.md(mdn).tchoOpt(:,bin1));
[~,IND]=max(tmp,[],2);
[~,ind] = sort(IND);
imagesc(tmp(ind,:))
colormap hot
colormap hot
colorbar
caxis(colorBarLim)
hold on
vline(vline1)
title('OFCout - Offer')
clear tmp IND ind

subplot(2,3,3)
tmp=abs(OFFERpcc.md(mdn).tchoOpt(:,bin1));
[~,IND]=max(tmp,[],2);
[~,ind] = sort(IND);
imagesc(tmp(ind,:))
colormap hot
colorbar
caxis(colorBarLim)
hold on
vline(vline1)
title('PCC - Offer')
clear tmp IND ind

subplot(2,3,4)
tmp=abs(CHOICEin.md(mdn).tchoOpt(:,bin2));
[~,IND]=max(tmp,[],2);
[~,ind] = sort(IND);
imagesc(tmp(ind,:))
colormap hot
colorbar
caxis(colorBarLim)
hold on
vline(vline2)
title('OFCin_Choice')
clear tmp IND ind

subplot(2,3,5)
tmp=abs(CHOICEout.md(mdn).tchoOpt(:,bin2));
[~,IND]=max(tmp,[],2);
[~,ind] = sort(IND);
imagesc(tmp(ind,:))
colormap hot
colorbar
caxis(colorBarLim)
hold on
vline(vline2)
title('OFCout - Choice')
clear tmp IND ind

subplot(2,3,6)
tmp=abs(CHOICEpcc.md(mdn).tchoOpt(:,bin2));
[~,IND]=max(tmp,[],2);
[~,ind] = sort(IND);
imagesc(tmp(ind,:))
colormap hot
colorbar
caxis(colorBarLim)
hold on
vline(vline2)
title('PCC - Choice')
clear tmp IND ind

%% sliding procedure
close all; clc 
Bwidth=30;
starts=421:5:551-Bwidth;
mdn=11;

for tt=1:length(starts)
    %%
    tt=9
    sBin=starts(tt):starts(tt)-1+30;
    tmp1=abs(OFFERpcc.md(mdn).tchoOpt(:,sBin));
    tmp2=abs(OFFERin.md(mdn).tchoOpt(:,sBin));
    tmp3=abs(OFFERout.md(mdn).tchoOpt(:,sBin));
    [~, latPCC]=max(tmp1,[],2);
    [~, latIn]=max(tmp2,[],2);
    [~, latOut]=max(tmp3,[],2);
    grpVar=zeros(length(latPCC)+length(latIn)+length(latOut),1);
    grpVar(1:length(latPCC),1)=1;
    grpVar(length(latPCC)+1:length(latPCC)+length(latIn),1)=2;
    grpVar(length(latPCC)+length(latIn)+1:length(latPCC)+length(latIn)+length(latOut),1)=3;
    
    yy=[latPCC;latIn;latOut];
    xx=categorical(grpVar);
    mdl = fitglm(xx,yy,'Distribution','gamma');
    mdl
    %%
    p=double(mdl.devianceTest.pValue);
    mdltsig(tt)=p(2);
    clear mdl
    
    sigPCC=OFFERpcc.md(mdn).pchoOpt(:,sBin)<0.05;
    sigIn=OFFERin.md(mdn).pchoOpt(:,sBin)<0.05;
    sigOut=OFFERout.md(mdn).pchoOpt(:,sBin)<0.05;
    % find first sig in a row
    [~, latPCC]=max(sigPCC,[],2);
    [~, latIn]=max(sigIn,[],2);
    [~, latOut]=max(sigOut,[],2);
    grpVar=zeros(length(latPCC)+length(latIn)+length(latOut),1);
    grpVar(1:length(latPCC),1)=1;
    grpVar(length(latPCC)+1:length(latPCC)+length(latIn),1)=2;
    grpVar(length(latPCC)+length(latIn)+1:length(latPCC)+length(latIn)+length(latOut),1)=3;
    
    yy=[latPCC;latIn;latOut];
    xx=categorical(grpVar);
    mdl = fitglm(xx,yy,'Distribution','gamma');
    p=double(mdl.devianceTest.pValue);
    mdlpsig(tt)=p(2);
    clear mdl
end
%%
find(mdltsig<0.05)
useTStart=starts(mdltsig<0.05)
find(mdlpsig<0.05)
usePStart=starts(mdlpsig<0.05)
% used mdltsig (9)= 461

%%
clc
inPCC=median(latPCC) % 15
inIn=median(latIn) % 9
inOut=median(latOut) % 17

% 461-420+17

%% ####################### OFFER 1 #######################
sBin=300:450;
tmp1=abs(OFFERpcc.md(mdn).tchoLR(:,sBin));
    tmp2=abs(OFFERin.md(mdn).tchoLR(:,sBin));
    tmp3=abs(OFFERout.md(mdn).tchoLR(:,sBin));
    [~, latPCC]=max(tmp1,[],2);
    [~, latIn]=max(tmp2,[],2);
    [~, latOut]=max(tmp3,[],2);
    grpVar=zeros(length(latPCC)+length(latIn)+length(latOut),1);
    grpVar(1:length(latPCC),1)=1;
    grpVar(length(latPCC)+1:length(latPCC)+length(latIn),1)=2;
    grpVar(length(latPCC)+length(latIn)+1:length(latPCC)+length(latIn)+length(latOut),1)=3;
    
    yy=[latPCC;latIn;latOut];
    xx=categorical(grpVar);
    mdl = fitglm(xx,yy,'Distribution','gamma')

%% ####################### OFFER 1 #######################
sBin=300:450;
tmp1=abs(OFFERpcc.md(mdn).tchoLR(:,sBin));
    tmp2=abs(OFFERin.md(mdn).tchoLR(:,sBin));
    tmp3=abs(OFFERout.md(mdn).tchoLR(:,sBin));
    [~, latPCC]=max(tmp1,[],2);
    [~, latIn]=max(tmp2,[],2);
    [~, latOut]=max(tmp3,[],2);
    grpVar=zeros(length(latPCC)+length(latIn)+length(latOut),1);
    grpVar(1:length(latPCC),1)=1;
    grpVar(length(latPCC)+1:length(latPCC)+length(latIn),1)=2;
    grpVar(length(latPCC)+length(latIn)+1:length(latPCC)+length(latIn)+length(latOut),1)=3;
    
    yy=[latPCC;latIn;latOut];
    xx=categorical(grpVar);
    mdl = fitglm(xx,yy,'Distribution','gamma')
