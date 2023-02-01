% posteriorPrediction
% last used: Mar 12 2020
% use the decoding posterior from each region to predict that of another
% region

% for figures and chi-square tests use 500ms step ones

% for Granger causality use 50ms step ones on server.

% Data info
% 4 fold cross validation -- 4 columns
% 6 time window
% at OFFER: 251:5:550
% at CHOICE:  201:5:500

clear all; close all;clc
% file path
fpath='/Users/mwang/Google Drive/01Data/00 w_Anatomy/2020final/code/Decode/';
addpath(genpath(fpath))
dpath='/Users/mwang/Google Drive/01Data/00 w_Anatomy/2020final/results/Decoding/results/step500ms/';
spath='/Users/mwang/Google Drive/01Data/00 w_Anatomy/2020final/results/Decoding/findings/';
cd(dpath)
List=dir('*decodePS*');

for fn=1:length(List)
    
    List(fn).name
    data(fn)=load(List(fn).name);
    
    
end

cd(spath)
%% decoding accuracy choOptC
clc;close all
Colors=[178,34,34;... %fire
    255,127,80;... %coral
    255,215,0]/255; % gold
Xs=1:6;

subplot(1,2,1)
e=errorbar(Xs,data(2).A.choOptC.lda,data(2).Asem.choOptC.lda,data(2).Asem.choOptC.lda);
e.Color=Colors(1,:);
xlim([0 7]);
ylim([0.4 0.7])
hold on
e=errorbar(Xs,data(4).A.choOptC.lda,data(4).Asem.choOptC.lda,data(4).Asem.choOptC.lda);
e.Color=Colors(2,:);
e=errorbar(Xs,data(6).A.choOptC.lda,data(6).Asem.choOptC.lda,data(6).Asem.choOptC.lda);
e.Color=Colors(3,:);
hline(0.5)
vline(1.5)
vline(3.5)
vline(5.5)
ylabel('Decoding Accuracy')
xlabel('OFFER: 251   301   351   401   451   501')
title('choOptC')

subplot(1,2,2)
e=errorbar(Xs,data(1).A.choOptC.lda,data(1).Asem.choOptC.lda,data(1).Asem.choOptC.lda);
e.Color=Colors(1,:);
xlim([0 7]);
ylim([0.4 0.7])
hold on
e=errorbar(Xs,data(3).A.choOptC.lda,data(3).Asem.choOptC.lda,data(3).Asem.choOptC.lda);
e.Color=Colors(2,:);
e=errorbar(Xs,data(5).A.choOptC.lda,data(5).Asem.choOptC.lda,data(5).Asem.choOptC.lda);
e.Color=Colors(3,:);
hline(0.5)
vline(2.5)
vline(4.5)
legend('In','Out', 'PCC')
ylabel('Decoding Accuracy')
xlabel('CHOICE: 201   251   301   351   401   451')

%% sig test
clc
data(2).A.choOptC.lda
data(2).A.choOptC.Size
data(2).A.choOptC.lda.*data(2).A.choOptC.Size
%% sig test
useBin=5;
clc
disp('OFCin')
data(2).A.choOptC.lda(useBin)
data(2).A.choOptC.Size(useBin)
data(2).A.choOptC.lda(useBin)*data(2).A.choOptC.Size(useBin)
disp('OFCout')
data(4).A.choOptC.lda(useBin)
data(4).A.choOptC.Size(useBin)
data(4).A.choOptC.lda(useBin)*data(4).A.choOptC.Size(useBin)
disp('PCC')
data(6).A.choOptC.lda(useBin)
data(6).A.choOptC.Size(useBin)
data(6).A.choOptC.lda(useBin)*data(6).A.choOptC.Size(useBin)
%% decoding accuracy choOptE
clc;close all
Colors=[178,34,34;... %fire
    255,127,80;... %coral
    255,215,0]/255; % gold
Xs=1:6;

subplot(1,2,1)
e=errorbar(Xs,data(2).A.choOptE.lda,data(2).Asem.choOptE.lda,data(2).Asem.choOptE.lda,':');
e.Color=Colors(1,:);
xlim([0 7]);
ylim([0.4 0.7])
hold on
e=errorbar(Xs,data(4).A.choOptE.lda,data(4).Asem.choOptE.lda,data(4).Asem.choOptE.lda,':');
e.Color=Colors(2,:);
e=errorbar(Xs,data(6).A.choOptE.lda,data(6).Asem.choOptE.lda,data(6).Asem.choOptE.lda,':');
e.Color=Colors(3,:);
hline(0.5)
vline(1.5)
vline(3.5)
vline(5.5)
ylabel('Decoding Accuracy')
xlabel('OFFER: 251   301   351   401   451   501')
title('choOptE')

subplot(1,2,2)
e=errorbar(Xs,data(1).A.choOptE.lda,data(1).Asem.choOptE.lda,data(1).Asem.choOptE.lda,':');
e.Color=Colors(1,:);
xlim([0 7]);
ylim([0.4 0.7])
hold on
e=errorbar(Xs,data(3).A.choOptE.lda,data(3).Asem.choOptE.lda,data(3).Asem.choOptE.lda,':');
e.Color=Colors(2,:);
e=errorbar(Xs,data(5).A.choOptE.lda,data(5).Asem.choOptE.lda,data(5).Asem.choOptE.lda,':');
e.Color=Colors(3,:);
hline(0.5)
vline(2.5)
vline(4.5)
legend('In','Out', 'PCC')
ylabel('Decoding Accuracy')
xlabel('CHOICE: 201   251   301   351   401   451')

%% sig test
useBin=5;
clc
disp('OFCin')
data(2).A.choOptE.lda(useBin)
data(2).A.choOptE.Size(useBin)
data(2).A.choOptE.lda(useBin)*data(2).A.choOptE.Size(useBin)
disp('OFCout')
data(4).A.choOptE.lda(useBin)
data(4).A.choOptE.Size(useBin)
data(4).A.choOptE.lda(useBin)*data(4).A.choOptE.Size(useBin)
disp('PCC')
data(6).A.choOptE.lda(useBin)
data(6).A.choOptE.Size(useBin)
data(6).A.choOptE.lda(useBin)*data(6).A.choOptE.Size(useBin)
%% decoding accuracy choLRC
clc;close all
Colors=[0,0,128;... % navy
    70,130,180;... % blue
    32,178,170]/255; % turquoise

Xs=1:6;

subplot(1,2,1)
e=errorbar(Xs,data(2).A.choLRC.lda,data(2).Asem.choLRC.lda,data(2).Asem.choLRC.lda);
e.Color=Colors(1,:);
xlim([0 7]);
ylim([0.4 0.7])
hold on
e=errorbar(Xs,data(4).A.choLRC.lda,data(4).Asem.choLRC.lda,data(4).Asem.choLRC.lda);
e.Color=Colors(2,:);
e=errorbar(Xs,data(6).A.choLRC.lda,data(6).Asem.choLRC.lda,data(6).Asem.choLRC.lda);
e.Color=Colors(3,:);
hline(0.5)
vline(1.5)
vline(3.5)
vline(5.5)
ylabel('Decoding Accuracy')
xlabel('OFFER: 251   301   351   401   451   501')
title('choLRC')

subplot(1,2,2)
e=errorbar(Xs,data(1).A.choLRC.lda,data(1).Asem.choLRC.lda,data(1).Asem.choLRC.lda);
e.Color=Colors(1,:);
xlim([0 7]);
ylim([0.4 0.7])
hold on
e=errorbar(Xs,data(3).A.choLRC.lda,data(3).Asem.choLRC.lda,data(3).Asem.choLRC.lda);
e.Color=Colors(2,:);
e=errorbar(Xs,data(5).A.choLRC.lda,data(5).Asem.choLRC.lda,data(5).Asem.choLRC.lda);
e.Color=Colors(3,:);
hline(0.5)
vline(2.5)
vline(4.5)
legend('In','Out', 'PCC')
ylabel('Decoding Accuracy')
xlabel('CHOICE: 201   251   301   351   401   451')
%% sig test at offer
useBin=5;
clc
disp('OFCin')
data(2).A.choLRC.lda(useBin)
data(2).A.choLRC.Size(useBin)
data(2).A.choLRC.lda(useBin)*data(2).A.choLRC.Size(useBin)
disp('OFCout')
data(4).A.choLRC.lda(useBin)
data(4).A.choLRC.Size(useBin)
data(4).A.choLRC.lda(useBin)*data(4).A.choLRC.Size(useBin)
disp('PCC')
data(6).A.choLRC.lda(useBin)
data(6).A.choLRC.Size(useBin)
data(6).A.choLRC.lda(useBin)*data(6).A.choLRC.Size(useBin)
%% sig test at choice
useBin=2;
clc
disp('OFCin')
data(1).A.choLRC.lda(useBin)
data(1).A.choLRC.Size(useBin)
data(1).A.choLRC.lda(useBin)*data(1).A.choLRC.Size(useBin)
disp('OFCout')
data(3).A.choLRC.lda(useBin)
data(3).A.choLRC.Size(useBin)
data(3).A.choLRC.lda(useBin)*data(3).A.choLRC.Size(useBin)
disp('PCC')
data(5).A.choLRC.lda(useBin)
data(5).A.choLRC.Size(useBin)
data(5).A.choLRC.lda(useBin)*data(5).A.choLRC.Size(useBin)
%% decoding accuracy choLRE
clc;close all
Colors=[0,0,128;... % navy
    70,130,180;... % blue
    32,178,170]/255; % turquoise

Xs=1:6;

subplot(1,2,1)
e=errorbar(Xs,data(2).A.choLRE.lda,data(2).Asem.choLRE.lda,data(2).Asem.choLRE.lda,':');
e.Color=Colors(1,:);
xlim([0 7]);
ylim([0.4 0.7])
hold on
e=errorbar(Xs,data(4).A.choLRE.lda,data(4).Asem.choLRE.lda,data(4).Asem.choLRE.lda,':');
e.Color=Colors(2,:);
e=errorbar(Xs,data(6).A.choLRE.lda,data(6).Asem.choLRE.lda,data(6).Asem.choLRE.lda,':');
e.Color=Colors(3,:);
hline(0.5)
vline(1.5)
vline(3.5)
vline(5.5)
ylabel('Decoding Accuracy')
xlabel('OFFER: 251   301   351   401   451   501')
title('choLRE')

subplot(1,2,2)
e=errorbar(Xs,data(1).A.choLRE.lda,data(1).Asem.choLRE.lda,data(1).Asem.choLRE.lda,':');
e.Color=Colors(1,:);
xlim([0 7]);
ylim([0.4 0.7])
hold on
e=errorbar(Xs,data(3).A.choLRE.lda,data(3).Asem.choLRE.lda,data(3).Asem.choLRE.lda,':');
e.Color=Colors(2,:);
e=errorbar(Xs,data(5).A.choLRE.lda,data(5).Asem.choLRE.lda,data(5).Asem.choLRE.lda,':');
e.Color=Colors(3,:);
hline(0.5)
vline(2.5)
vline(4.5)
legend('In','Out', 'PCC')
ylabel('Decoding Accuracy')
xlabel('CHOICE: 201   251   301   351   401   451')

%% sig test at choice
useBin=2;
clc
disp('OFCin')
data(1).A.choLRE.lda(useBin)
data(1).A.choLRE.Size(useBin)
data(1).A.choLRE.lda(useBin)*data(1).A.choLRE.Size(useBin)
disp('OFCout')
data(3).A.choLRE.lda(useBin)
data(3).A.choLRE.Size(useBin)
data(3).A.choLRE.lda(useBin)*data(3).A.choLRE.Size(useBin)
disp('PCC')
data(5).A.choLRE.lda(useBin)
data(5).A.choLRE.Size(useBin)
data(5).A.choLRE.lda(useBin)*data(5).A.choLRE.Size(useBin)
%% decoding accuracy ev1highC
clc;close all
Colors=[75,0,130;... % indigo
    147,112,219;... % purple
    199,21,133]/255; % turquoise

Xs=1:6;

subplot(1,2,1)
e=errorbar(Xs,data(2).A.ev1highC.lda,data(2).Asem.ev1highC.lda,data(2).Asem.ev1highC.lda);
e.Color=Colors(1,:);
xlim([0 7]);
ylim([0.4 0.7])
hold on
e=errorbar(Xs,data(4).A.ev1highC.lda,data(4).Asem.ev1highC.lda,data(4).Asem.ev1highC.lda);
e.Color=Colors(2,:);
e=errorbar(Xs,data(6).A.ev1highC.lda,data(6).Asem.ev1highC.lda,data(6).Asem.ev1highC.lda);
e.Color=Colors(3,:);
hline(0.5)
vline(1.5)
vline(3.5)
vline(5.5)
ylabel('Decoding Accuracy')
xlabel('OFFER: 251   301   351   401   451   501')
title('ev1highC')

subplot(1,2,2)
e=errorbar(Xs,data(1).A.ev1highC.lda,data(1).Asem.ev1highC.lda,data(1).Asem.ev1highC.lda);
e.Color=Colors(1,:);
xlim([0 7]);
ylim([0.4 0.7])
hold on
e=errorbar(Xs,data(3).A.ev1highC.lda,data(3).Asem.ev1highC.lda,data(3).Asem.ev1highC.lda);
e.Color=Colors(2,:);
e=errorbar(Xs,data(5).A.ev1highC.lda,data(5).Asem.ev1highC.lda,data(5).Asem.ev1highC.lda);
e.Color=Colors(3,:);
hline(0.5)
vline(2.5)
vline(4.5)
legend('In','Out', 'PCC')
ylabel('Decoding Accuracy')
xlabel('CHOICE: 201   251   301   351   401   451')

%% sig test at offer
useBin=3;
clc
disp('OFCin')
data(2).A.ev1highC.lda(useBin)
data(2).A.ev1highC.Size(useBin)
data(2).A.ev1highC.lda(useBin)*data(2).A.ev1highC.Size(useBin)
disp('OFCout')
data(4).A.ev1highC.lda(useBin)
data(4).A.ev1highC.Size(useBin)
data(4).A.ev1highC.lda(useBin)*data(4).A.ev1highC.Size(useBin)
disp('PCC')
data(6).A.ev1highC.lda(useBin)
data(6).A.ev1highC.Size(useBin)
data(6).A.ev1highC.lda(useBin)*data(6).A.ev1highC.Size(useBin)
%% decoding accuracy ev1highE
clc;close all
Colors=[75,0,130;... % indigo
    147,112,219;... % purple
    199,21,133]/255; % turquoise

Xs=1:6;

subplot(1,2,1)
e=errorbar(Xs,data(2).A.ev1highE.lda,data(2).Asem.ev1highE.lda,data(2).Asem.ev1highE.lda,':');
e.Color=Colors(1,:);
xlim([0 7]);
ylim([0.4 0.7])
hold on
e=errorbar(Xs,data(4).A.ev1highE.lda,data(4).Asem.ev1highE.lda,data(4).Asem.ev1highE.lda,':');
e.Color=Colors(2,:);
e=errorbar(Xs,data(6).A.ev1highE.lda,data(6).Asem.ev1highE.lda,data(6).Asem.ev1highE.lda,':');
e.Color=Colors(3,:);
hline(0.5)
vline(1.5)
vline(3.5)
vline(5.5)
ylabel('Decoding Accuracy')
xlabel('OFFER: 251   301   351   401   451   501')
title('ev1highE')

subplot(1,2,2)
e=errorbar(Xs,data(1).A.ev1highE.lda,data(1).Asem.ev1highE.lda,data(1).Asem.ev1highE.lda,':');
e.Color=Colors(1,:);
xlim([0 7]);
ylim([0.4 0.7])
hold on
e=errorbar(Xs,data(3).A.ev1highE.lda,data(3).Asem.ev1highE.lda,data(3).Asem.ev1highE.lda,':');
e.Color=Colors(2,:);
e=errorbar(Xs,data(5).A.ev1highE.lda,data(5).Asem.ev1highE.lda,data(5).Asem.ev1highE.lda,':');
e.Color=Colors(3,:);
hline(0.5)
vline(2.5)
vline(4.5)
legend('In','Out', 'PCC')
ylabel('Decoding Accuracy')
xlabel('CHOICE: 201   251   301   351   401   451')

%% decoding accuracy ev2highC
clc;close all
Colors=[75,0,130;... % indigo
    147,112,219;... % purple
    199,21,133]/255; % turquoise

Xs=1:6;

subplot(1,2,1)
e=errorbar(Xs,data(2).A.ev2highC.lda,data(2).Asem.ev2highC.lda,data(2).Asem.ev2highC.lda);
e.Color=Colors(1,:);
xlim([0 7]);
ylim([0.4 0.7])
hold on
e=errorbar(Xs,data(4).A.ev2highC.lda,data(4).Asem.ev2highC.lda,data(4).Asem.ev2highC.lda);
e.Color=Colors(2,:);
e=errorbar(Xs,data(6).A.ev2highC.lda,data(6).Asem.ev2highC.lda,data(6).Asem.ev2highC.lda);
e.Color=Colors(3,:);
hline(0.5)
vline(1.5)
vline(3.5)
vline(5.5)
ylabel('Decoding Accuracy')
xlabel('OFFER: 251   301   351   401   451   501')
title('ev2highC')

subplot(1,2,2)
e=errorbar(Xs,data(1).A.ev2highC.lda,data(1).Asem.ev2highC.lda,data(1).Asem.ev2highC.lda);
e.Color=Colors(1,:);
xlim([0 7]);
ylim([0.4 0.7])
hold on
e=errorbar(Xs,data(3).A.ev2highC.lda,data(3).Asem.ev2highC.lda,data(3).Asem.ev2highC.lda);
e.Color=Colors(2,:);
e=errorbar(Xs,data(5).A.ev2highC.lda,data(5).Asem.ev2highC.lda,data(5).Asem.ev2highC.lda);
e.Color=Colors(3,:);
hline(0.5)
vline(2.5)
vline(4.5)
legend('In','Out', 'PCC')
ylabel('Decoding Accuracy')
xlabel('CHOICE: 201   251   301   351   401   451')

%% decoding accuracy ev2highE
clc;close all
Colors=[75,0,130;... % indigo
    147,112,219;... % purple
    199,21,133]/255; % turquoise

Xs=1:6;

subplot(1,2,1)
e=errorbar(Xs,data(2).A.ev2highE.lda,data(2).Asem.ev2highE.lda,data(2).Asem.ev2highE.lda,':');
e.Color=Colors(1,:);
xlim([0 7]);
ylim([0.4 0.7])
hold on
e=errorbar(Xs,data(4).A.ev2highE.lda,data(4).Asem.ev2highE.lda,data(4).Asem.ev2highE.lda,':');
e.Color=Colors(2,:);
e=errorbar(Xs,data(6).A.ev2highE.lda,data(6).Asem.ev2highE.lda,data(6).Asem.ev2highE.lda,':');
e.Color=Colors(3,:);
hline(0.5)
vline(1.5)
vline(3.5)
vline(5.5)
ylabel('Decoding Accuracy')
xlabel('OFFER: 251   301   351   401   451   501')
title('ev2highE')

subplot(1,2,2)
e=errorbar(Xs,data(1).A.ev2highE.lda,data(1).Asem.ev2highE.lda,data(1).Asem.ev2highE.lda,':');
e.Color=Colors(1,:);
xlim([0 7]);
ylim([0.4 0.7])
hold on
e=errorbar(Xs,data(3).A.ev2highE.lda,data(3).Asem.ev2highE.lda,data(3).Asem.ev2highE.lda,':');
e.Color=Colors(2,:);
e=errorbar(Xs,data(5).A.ev2highE.lda,data(5).Asem.ev2highE.lda,data(5).Asem.ev2highE.lda,':');
e.Color=Colors(3,:);
hline(0.5)
vline(2.5)
vline(4.5)
legend('In','Out', 'PCC')
ylabel('Decoding Accuracy')
xlabel('CHOICE: 201   251   301   351   401   451')

%%
% similar was done to posterior --> not as pretty









