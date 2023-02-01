% PCAt5
% last used May 8 2020
% This pandemic is really weighing me down
% plots the trial-averaged distance against shuffled trial-averaged dist

clear all; close all; clc

dpath1='/Users/mwang/Google Drive/01Data/00 w_Anatomy/2020final/results/PCAt/';
dpath2='/Users/mwang/Google Drive/01Data/00 w_Anatomy/2020final/results/PCAt/ShuffledTrlAvgDist/';
fpath='/Users/mwang/Google Drive/01Data/00 w_Anatomy/2020final/code/Trajectory/';
spath='/Users/mwang/Google Drive/01Data/00 w_Anatomy/2020final/results/PCAt/Figures/';
addpath(genpath(fpath))

cd(dpath2)

% List=dir('*P2017*DtrlAvgShuffle*');
List=dir('*S2018*DtrlAvgShuffle*');

% Subject 1=pumbaa; 2=spock;
% subn=1;
subn=2;

load(List(1).name);
rIn=rD;
clear rD
load(List(2).name);
rOut=rD;
clear rD
load(List(3).name);
rPcc=rD;
clear rD
clear List

cd(dpath1)
load('OFCinDistAvg.mat');
in=D;
clear D
load('OFCoutDistAvg.mat');
out=D;
clear D
load('PCCDistAvg.mat');
pcc=D;
clear D

% load('OFCinDistTbT.mat');
% in=td(subn);
% clear td
% load('OFCoutDistTbT.mat');
% out=td(subn);
% clear td
% load('PCCDistTbT.mat');
% pcc=td(subn);
% clear td

stepSize=5;
binStarts=251:stepSize:750-stepSize+1;

% ########################################

useDr=rIn.c12;
% useDr=rIn.clr;
% useDr=rIn.cev1;
inL25=prctile(useDr,2.5,2);
inH25=prctile(useDr,97.5,2);
clear useDr
inD=in.c12(subn,:);  
% inD=in.clr(subn,:);
% inD=in.cev1(subn,:);

useDr=rPcc.c12;
% useDr=rPcc.clr;
% useDr=rPcc.cev1;
pccL25=prctile(useDr,2.5,2);
pccH25=prctile(useDr,97.5,2);
clear useDr
pccD=pcc.c12(subn,:);
% pccD=pcc.clr(subn,:);
% pccD=pcc.cev1(subn,:);

useDr=rOut.c12;
% useDr=rOut.clr;
% useDr=rOut.cev1;
outL25=prctile(useDr,2.5,2);
outH25=prctile(useDr,97.5,2);
clear useDr
outD=out.c12(subn,:);
% outD=out.clr(subn,:);
% outD=out.cev1(subn,:);

titleContent='Choice 1 vs 2, Correct';
% titleContent='Choice L vs R, Correct';
% titleContent='EV1 H vs L, Correct';


% useD=mean([in.c12.Dst.aD12;in.c12.Dst.aD21],1)';
% useD=mean([pcc.clr.Dst.aD12;pcc.clr.Dst.aD21],1)';
%% 
clc;close all

subplot(1,3,1)
plot(inD,'r-','lineWidth',1.5)
hold on
plot(inH25,'r:','lineWidth',1.5)
plot(inL25,'r:','lineWidth',1.5)
ylim([0 7]) % 3.2
vline(11);
vline(31);
vline(51);
vline(71);
vline(91);
legend('Mean','upper 2.5%','lower 2.5%',...
    num2str(binStarts(11)),num2str(binStarts(31)),...
    num2str(binStarts(51)),...
    num2str(binStarts(71)),num2str(binStarts(91)));

hold off

subplot(1,3,2)
plot(pccD,'b-','lineWidth',1.5)
hold on
plot(pccH25,'b:','lineWidth',1.5)
plot(pccL25,'b:','lineWidth',1.5)
ylim([0 7]) % 3.2
vline(11);
vline(31);
vline(51);
vline(71);
vline(91);
legend('Mean','upper 2.5%','lower 2.5%',...
    num2str(binStarts(11)),num2str(binStarts(31)),...
    num2str(binStarts(51)),...
    num2str(binStarts(71)),num2str(binStarts(91)));
title(titleContent)
hold off

subplot(1,3,3)
plot(outD,'k-','lineWidth',1.5)
hold on
plot(outH25,'k:','lineWidth',1.5)
plot(outL25,'k:','lineWidth',1.5)
ylim([0 7]) %3.2
vline(11);
vline(31);
vline(51);
vline(71);
vline(91);
legend('Mean','upper 2.5%','lower 2.5%',...
    num2str(binStarts(11)),num2str(binStarts(31)),...
    num2str(binStarts(51)),...
    num2str(binStarts(71)),num2str(binStarts(91)));
xlabel('Time');
ylabel('Distance');

%% for error trials

clear all; close all; clc

dpath1='/Users/mwang/Google Drive/01Data/00 w_Anatomy/2020final/results/PCAt/';
dpath2='/Users/mwang/Google Drive/01Data/00 w_Anatomy/2020final/results/PCAt/ShuffledTrlAvgDist/';
fpath='/Users/mwang/Google Drive/01Data/00 w_Anatomy/2020final/code/Trajectory/';
spath='/Users/mwang/Google Drive/01Data/00 w_Anatomy/2020final/results/PCAt/Figures/';
addpath(genpath(fpath))

cd(dpath2)

% List=dir('*P2017*DtrlAvgShuffle*');
List=dir('*S2018*DtrlAvgShuffle*');

% Subject 1=pumbaa; 2=spock;
% subn=1;
subn=2;

load(List(1).name);
rIn=rD;
clear rD
load(List(2).name);
rOut=rD;
clear rD
load(List(3).name);
rPcc=rD;
clear rD
clear List

cd(dpath1)
load('OFCinDistAvg.mat');
in=D;
clear D
load('OFCoutDistAvg.mat');
out=D;
clear D
load('PCCDistAvg.mat');
pcc=D;
clear D

% load('OFCinDistTbT.mat');
% in=td(subn);
% clear td
% load('OFCoutDistTbT.mat');
% out=td(subn);
% clear td
% load('PCCDistTbT.mat');
% pcc=td(subn);
% clear td

stepSize=5;
binStarts=251:stepSize:750-stepSize+1;

% ########################################

useDr=rIn.e12;
% useDr=rIn.elr;
% useDr=rIn.eev1;
inL25=prctile(useDr,2.5,2);
inH25=prctile(useDr,97.5,2);
clear useDr
inD=in.e12(subn,:); 
% inD=mean([in.e12.Dst.aD12;in.e12.Dst.aD21],1)'; % doesn't work too large
% inD=in.elr(subn,:);
% inD=in.eev1(subn,:);

useDr=rPcc.e12;
% useDr=rPcc.elr;
% useDr=rPcc.eev1;
pccL25=prctile(useDr,2.5,2);
pccH25=prctile(useDr,97.5,2);
clear useDr
pccD=pcc.e12(subn,:);
% pccD=pcc.elr(subn,:);
% useD=mean([pcc.elr.Dst.aD12;pcc.elr.Dst.aD21],1)';
% pccD=pcc.eev1(subn,:);

useDr=rOut.e12;
% useDr=rOut.elr;
% useDr=rOut.eev1;
outL25=prctile(useDr,2.5,2);
outH25=prctile(useDr,97.5,2);
clear useDr
outD=out.e12(subn,:);
% outD=out.elr(subn,:);
% outD=out.eev1(subn,:);

% titleContent='Choice 1 vs 2, Error';
% titleContent='Choice L vs R, Error';
titleContent='EV1 H vs L, Error';


%% 
clc;close all

subplot(1,3,1)
plot(inD,'r-','lineWidth',1.5)
hold on
plot(inH25,'r:','lineWidth',1.5)
plot(inL25,'r:','lineWidth',1.5)
ylim([-0.5 15]) % 3.2
vline(11);
vline(31);
vline(51);
vline(71);
vline(91);
legend('Mean','upper 2.5%','lower 2.5%',...
    num2str(binStarts(11)),num2str(binStarts(31)),...
    num2str(binStarts(51)),...
    num2str(binStarts(71)),num2str(binStarts(91)));

hold off

subplot(1,3,2)
plot(pccD,'b-','lineWidth',1.5)
hold on
plot(pccH25,'b:','lineWidth',1.5)
plot(pccL25,'b:','lineWidth',1.5)
ylim([-0.5 15]) % 3.2
vline(11);
vline(31);
vline(51);
vline(71);
vline(91);
legend('Mean','upper 2.5%','lower 2.5%',...
    num2str(binStarts(11)),num2str(binStarts(31)),...
    num2str(binStarts(51)),...
    num2str(binStarts(71)),num2str(binStarts(91)));
title(titleContent)
hold off

subplot(1,3,3)
plot(outD,'k-','lineWidth',1.5)
hold on
plot(outH25,'k:','lineWidth',1.5)
plot(outL25,'k:','lineWidth',1.5)
ylim([-0.5 15]) %3.2
vline(11);
vline(31);
vline(51);
vline(71);
vline(91);
legend('Mean','upper 2.5%','lower 2.5%',...
    num2str(binStarts(11)),num2str(binStarts(31)),...
    num2str(binStarts(51)),...
    num2str(binStarts(71)),num2str(binStarts(91)));
xlabel('Time');
ylabel('Distance');



