% PCAt2
% last used: May 14 2020
% use data from PrepPCAt2020 & PCAt1
% what this does:

% calculate distance on trial-averaged trajectories
% SEPARATED BY SUBJECT
% Instead of averaged across subjects
% Bc the time dynamics were different for each

% granger causality eventually used in the paper at near bottom

clear all; close all; clc
dpath='/Users/mwang/Google Drive/01Data/00 w_Anatomy/2020final/wrapped/PCAt_states/';
fpath='/Users/mwang/Google Drive/01Data/00 w_Anatomy/2020final/code/';
spath='/Users/mwang/Google Drive/01Data/00 w_Anatomy/2020final/results/PCAt/';
addpath(genpath(fpath))

cd(dpath)

%%
useNames={'OFCin','OFCout','PCC'};
for nn=1:length(useNames)
    clear List
    cd(dpath)
    List=dir(['*' useNames{nn} '*PCA*']);
    for fn=1:length(List)
        load(List(fn).name);
        for tp=1:size(C12.top70PCs,1)/2
            clear p1 p2
            p1=C12.top70PCs(tp,:);
            p2=C12.top70PCs(tp+100,:);
            D.c12(fn,tp)=adjDist1(p1,p2);
            
            clear p1 p2
            p1=E12.top70PCs(tp,:);
            p2=E12.top70PCs(tp+100,:);
            D.e12(fn,tp)=adjDist1(p1,p2);          
            
            clear p1 p2
            p1=cLR.top70PCs(tp,:);
            p2=cLR.top70PCs(tp+100,:);
            D.clr(fn,tp)=adjDist1(p1,p2);
            
            clear p1 p2
            p1=eLR.top70PCs(tp,:);
            p2=eLR.top70PCs(tp+100,:);
            D.elr(fn,tp)=adjDist1(p1,p2);
            
            clear p1 p2
            p1=cEV1.top70PCs(tp,:);
            p2=cEV1.top70PCs(tp+100,:);
            D.cev1(fn,tp)=adjDist1(p1,p2);
            
            clear p1 p2
            p1=eEV1.top70PCs(tp,:);
            p2=eEV1.top70PCs(tp+100,:);
            D.eev1(fn,tp)=adjDist1(p1,p2);
            
        end % tp time point
    end % fn flie in a list
    
    D.memo={'D.condition(row=subject, 1=P,2=S;',...
        'column=timepoint 251:5:750)'};

    save([spath useNames{nn} 'DistAvg.mat'], 'D')
    
    cd(spath)
    
end % nn next list

clear List1 List2

%%
clear all; close all; clc
dpath='/Users/mwang/Google Drive/01Data/00 w_Anatomy/2020final/wrapped/PCAt_states/';
fpath='/Users/mwang/Google Drive/01Data/00 w_Anatomy/2020final/code/';
spath='/Users/mwang/Google Drive/01Data/00 w_Anatomy/2020final/results/PCAt/';

cd(spath)
load('OFCinDistAvg.mat');
in=D;
clear D
load('OFCoutDistAvg.mat');
out=D;
clear D
load('PCCDistAvg.mat');
pcc=D;
clear D
%% pick a subject
% pumbaa
subn=1; 
% spock
% subn=2;

%% just plot the distance
subplot(1,3,1)
plot(in.c12(subn,:),'r-','lineWidth',1.5)
hold on
plot(in.e12(subn,:),'r:','lineWidth',1.5)
plot(pcc.c12(subn,:),'b-','lineWidth',1.5)
plot(pcc.e12(subn,:),'b:','lineWidth',1.5)
plot(out.c12(subn,:),'k-','lineWidth',1.5)
plot(out.e12(subn,:),'k:','lineWidth',1.5)
ylim([-1,5])
title('1 vs 2')
hold off

subplot(1,3,2)
plot(in.clr(subn,:),'r-','lineWidth',1.5)
hold on
plot(in.elr(subn,:),'r:','lineWidth',1.5)
plot(pcc.clr(subn,:),'b-','lineWidth',1.5)
plot(pcc.elr(subn,:),'b:','lineWidth',1.5)
plot(out.clr(subn,:),'k-','lineWidth',1.5)
plot(out.elr(subn,:),'k:','lineWidth',1.5)
ylim([-1,5])
title('L vs R')
hold off

subplot(1,3,3)
plot(in.cev1(subn,:),'r-','lineWidth',1.5)
hold on
plot(in.eev1(subn,:),'r:','lineWidth',1.5)
plot(pcc.cev1(subn,:),'b-','lineWidth',1.5)
plot(pcc.eev1(subn,:),'b:','lineWidth',1.5)
plot(out.cev1(subn,:),'k-','lineWidth',1.5)
plot(out.eev1(subn,:),'k:','lineWidth',1.5)
ylim([-1,5])
title('EV1 high vs low')
hold off
%%
subplot(1,2,1)
hist(pcc.clr)
subplot(1,2,2)
hist(out.clr)
%% comparing correct and error 
clc; clear dd cData eData

% Subject: pumbaa=1; spock=2
subn=1; 

% Area
dd=in;
% dd=out;
% dd=pcc;

cData=dd.cev1(subn,:)';
eData=dd.eev1(subn,:)';

% cData=dd.c12(subn,:)';
% eData=dd.e12(subn,:)';

% cData=dd.clr(subn,:)';
% eData=dd.elr(subn,:)';
%% Kruskal Wallis - C vs E 
[p,tbl,stats]=kruskalwallis([cData,eData])

%% Two-sample Kolmogorov-Smirnov
[h,p,ks2stat] = kstest2(cData,eData,'Tail','unequal') 
figure(1)
histogram(cData,ceil(length(cData)*0.4));
hold on
histogram(eData,ceil(length(eData)*0.4));
legend('cData','eData')

figure(2)
ecdf(cData,'Bounds','on')
hold on
cdfplot(eData)
grid on
legend('cData','','','eData')
%% KS alternative Hypo: 'smaller' 'larger'
[h,p,ks2stat] = kstest2(cData,eData,'Tail','smaller') 

%% Wilcoxon signed rank test
[p,h,stats] = signrank(cData,eData)
%% comparing across 3 areas
clc; close all; clear Din Dout Dpcc subn

% Subject: pumbaa=1; spock=2
subn=1; 

% just correct
% Din=in.cev1(subn,:)';
% Dout=out.cev1(subn,:)';
% Dpcc=pcc.cev1(subn,:)';

% Din=in.c12(subn,:)';
% Dout=out.c12(subn,:)';
% Dpcc=pcc.c12(subn,:)';

Din=in.clr(subn,:)';
Dout=out.clr(subn,:)';
Dpcc=pcc.clr(subn,:)';

% just error

% Din=in.eev1(subn,:)';
% Dout=out.eev1(subn,:)';
% Dpcc=pcc.eev1(subn,:)';

% Din=in.e12(subn,:)';
% Dout=out.e12(subn,:)';
% Dpcc=pcc.e12(subn,:)';

% Din=in.elr(subn,:)';
% Dout=out.elr(subn,:)';
% Dpcc=pcc.elr(subn,:)';

%%
%% comparing across 3 areas
clc; close all; clear Din Dout Dpcc 

% just correct
% Din=[in.cev1(1,:)';in.cev1(2,:)'];
% Dout=[out.cev1(1,:)';out.cev1(2,:)'];
% Dpcc=[pcc.cev1(1,:)';pcc.cev1(2,:)'];

% Din=[in.c12(1,:)';in.c12(2,:)'];
% Dout=[out.c12(1,:)';out.c12(2,:)'];
% Dpcc=[pcc.c12(1,:)';pcc.c12(2,:)'];

% Din=[in.clr(1,:)';in.clr(2,:)'];
% Dout=[out.clr(1,:)';out.clr(2,:)'];
% Dpcc=[pcc.clr(1,:)';pcc.clr(2,:)'];

% just error
% Din=[in.eev1(1,:)';in.eev1(2,:)'];
% Dout=[out.eev1(1,:)';out.eev1(2,:)'];
% Dpcc=[pcc.eev1(1,:)';pcc.eev1(2,:)'];

% Din=[in.e12(1,:)';in.e12(2,:)'];
% Dout=[out.e12(1,:)';out.e12(2,:)'];
% Dpcc=[pcc.e12(1,:)';pcc.e12(2,:)'];

% Din=[in.elr(1,:)';in.elr(2,:)'];
% Dout=[out.elr(1,:)';out.elr(2,:)'];
% Dpcc=[pcc.elr(1,:)';pcc.elr(2,:)'];

%% Kruskal Wallis

[p,tbl,stats]=kruskalwallis([Din,Dout,Dpcc])
figure(3)
multcompare(stats)

%% ######################## Granger Causality ########################
% #################################both###################################
clear all; close all; clc
dpath='/Users/mwang/Google Drive/01Data/00 w_Anatomy/2020final/wrapped/PCAt_states/';
fpath='/Users/mwang/Google Drive/01Data/00 w_Anatomy/2020final/code/';
spath='/Users/mwang/Google Drive/01Data/00 w_Anatomy/2020final/results/PCAt/';

cd(spath)
load('OFCinDistAvg.mat');
in=D;
clear D
load('OFCoutDistAvg.mat');
out=D;
clear D
load('PCCDistAvg.mat');
pcc=D;
clear D


% average
in12c=nanmean(in.c12,1)';
in12e=nanmean(in.e12,1)';
inlrc=nanmean(in.clr,1)';
inlre=nanmean(in.elr,1)';

pcc12c=nanmean(pcc.c12,1)';
pcc12e=nanmean(pcc.e12,1)';
pcclrc=nanmean(pcc.clr,1)';
pcclre=nanmean(pcc.elr,1)';

out12c=nanmean(out.c12,1)';
out12e=nanmean(out.e12,1)';
outlrc=nanmean(out.clr,1)';
outlre=nanmean(out.elr,1)';

%% stationary test
clc
% test stationary: augmented Dickey-Fuller test, h = adftest(var,'Model',"ard")
h=adftest(in12c,'Model','ard') % non-stationary  
% detrend: detrend(sdata); % detrended out data is stationary
% adftest(detrend(in12c),'Model','ard') 
if h==0 % non-stationary
    trendSwitch=1;
else
    trendSwitch=0;
end
trendSwitch
%% VAR model check num of lags to use
clc;

X1 = in12c;
X2 = pcclrc;
T = numel(X1);

numseries = 2;
numlags = (1:14)'; 
nummdls = numel(numlags);

% Partition time base.
maxp = max(numlags); % Maximum number of required presample responses
idxpre = 1:maxp;
idxest = (maxp + 1):T;

% Preallocation
EstMdl(nummdls) = varm(numseries,0);
aic = zeros(nummdls,1);

% Fit VAR models to data.
Y0 = [X1(idxpre) X2(idxpre)]; % Presample
Y = [X1(idxest) X2(idxest)];  % Estimation sample
for j = 1:numel(numlags)
    Mdl = varm(numseries,numlags(j));
    Mdl.SeriesNames = ["in12" "pcclr"];
    EstMdl(j) = estimate(Mdl,Y,'Y0',Y0);
    results = summarize(EstMdl(j));
    aic(j) = results.AIC;
end

%%
clc
NLAG = numlags(aic == min(aic))
trendSwitch
[h,pvalue,stat] = gctest(X1,X2,'NumLags', NLAG, 'Trend', trendSwitch)
%%
clc
[h,pvalue,stat] = gctest(in12c,pcclrc,[pcc12c,inlrc,out12c,outlrc],....
    'NumLags', NLAG, 'Trend',trendSwitch, 'Alpha',0.025)
%%
clc
[h,pvalue,stat] = gctest(pcclrc,in12c,[pcc12c,inlrc,out12c,outlrc],....
    'NumLags', 7, 'Trend',trendSwitch, 'Alpha',0.025)
    
%%
clc; close all
trendSwitch=1;
for nlags=1:11 
    
    % dist 12 --> lr 
    [~,pvalue(nlags),stat(nlags),cvalue(nlags)] = gctest(in12c,pcclrc,...
        [pcc12c,inlrc,out12c,outlrc],....
        'NumLags', nlags, 'Trend',trendSwitch, 'Alpha',0.025);
    
    [~,pvalue2(nlags),stat2(nlags),cvalue2(nlags)] = gctest(pcclrc,in12c,...
        [pcc12c,inlrc,out12c,outlrc],....
        'NumLags', nlags, 'Trend',trendSwitch, 'Alpha',0.025);
%     
%     [~,pvalue3(nlags),stat3(nlags),cvalue3(nlags)] = gctest(pcc12c,pcclrc,...
%         [in12c,inlrc,out12c,outlrc],....
%         'NumLags', nlags, 'Trend',trendSwitch, 'Alpha',0.025);

    [~,pvalue3(nlags),stat3(nlags),cvalue3(nlags)] = gctest(out12c,pcclrc,...
        [in12c,inlrc,pcc12c,outlrc],....
        'NumLags', nlags, 'Trend',trendSwitch, 'Alpha',0.025);

%     [h4(nlags),pvalue4(nlags),stat4(nlags),cvalue4(nlags)] = gctest(inlrc,pcclrc,...
%         [in12c,pcc12c,out12c,outlrc],....
%         'NumLags', nlags, 'Trend',1, 'Alpha',0.025);
    
end
%%
% 3*50=150
% 7*50=350
% 9*50=450
find(pvalue3<0.025)
stat2(7)
pvalue2(7)
%%
stat(9)
pvalue(9)
%%
close all
plot(stat,'r-','lineWidth',2)
hold on
plot(find(pvalue<0.025), stat(pvalue<0.025),'ro','lineWidth',2)
% plot(cvalue,'r:','lineWidth',2)

plot(stat2,'b-','lineWidth',2)
plot(find(pvalue2<0.025), stat2(pvalue2<0.025),'bo','lineWidth',2)
% plot(cvalue2,'b:','lineWidth',2)

plot(stat3,'k-','lineWidth',1.5)
plot(find(pvalue3<0.025), stat3(pvalue3<0.025),'ko','lineWidth',2)
% plot(cvalue3,'k:','lineWidth',2)

legend('in.c12->pcc.clr','sig',...
    'pcc.clr->in.c12','sig',...
    'pcc.c12->pcc.clr','sig')
% 
% plot(stat4,'b:','lineWidth',1.5)
% plot(find(pvalue4<0.025), stat4(pvalue4<0.025),'bo','lineWidth',1.5)
xlim([0 12])
%%
clc
nn=9 % 3 9
stat(nn)
pvalue(nn)

%%
clc
nn=7
stat2(nn)
pvalue2(nn)

%%
clc
nn=3
stat3(nn)
pvalue3(nn)



%% Granger Causality: single subject
clc
% Subject: pumbaa=1; spock=2
subn=1; 
% when check "feedback" granger causality, reduce significance to half
% 'Alpha', 0.025 then test A-->B and then B-->A

% test stationary: augmented Dickey-Fuller test, h = adftest(var,'Model',"ard")
% adftest(in.c12(subn,:)','Model','ard') % stationary 
% adftest(pcc.clr(subn,:)','Model','ard') % stationary 
adftest(detrend(out.c12(subn,:)'),'Model','ard') % non-stationary 
% adftest(out.c12(subn,:)','Model','ard')
% detrend: detrend(sdata); % detrended out data is stationary
%% VAR model check num of lags to use
clc;

X1 = in.c12(subn,:)';
X2 = pcc.clr(subn,:)';
T = numel(X1);

numseries = 2;
numlags = (1:12)'; % bc that's 2 cycles of 1 Hz?
nummdls = numel(numlags);

% Partition time base.
maxp = max(numlags); % Maximum number of required presample responses
idxpre = 1:maxp;
idxest = (maxp + 1):T;

% Preallocation
EstMdl(nummdls) = varm(numseries,0);
aic = zeros(nummdls,1);

% Fit VAR models to data.
Y0 = [X1(idxpre) X2(idxpre)]; % Presample
Y = [X1(idxest) X2(idxest)];  % Estimation sample
for j = 1:numel(numlags)
    Mdl = varm(numseries,numlags(j));
    Mdl.SeriesNames = ["in12" "pcclr"];
    EstMdl(j) = estimate(Mdl,Y,'Y0',Y0);
    results = summarize(EstMdl(j));
    aic(j) = results.AIC;
end

nlag = numlags(aic == min(aic))

[h,pvalue,stat,cvalue] = gctest(X1,X2,'NumLags', nlag, 'Trend',0)

%%
clc; close all



for nlags=1:12 
    
    % dist 12 --> lr 
    [h(nlags),pvalue(nlags),stat(nlags),cvalue(nlags)] = gctest(in.c12(subn,:)',pcc.clr(subn,:)',...
        [pcc.c12(subn,:)',in.clr(subn,:)',out.c12(subn,:)',out.clr(subn,:)'],....
        'NumLags', nlags, 'Trend',1, 'Alpha',0.025);
    
    [h2(nlags),pvalue2(nlags),stat2(nlags),cvalue2(nlags)] = gctest(pcc.clr(subn,:)',in.c12(subn,:)',...
        [pcc.c12(subn,:)',in.clr(subn,:)',out.c12(subn,:)',out.clr(subn,:)'],....
        'NumLags', nlags, 'Trend',1, 'Alpha',0.025);
    
%     [h3(nlags),pvalue3(nlags),stat3(nlags),cvalue3(nlags)] = gctest(in.e12(subn,:)',pcc.elr(subn,:)',...
%         'NumLags', nlags, 'Trend',1, 'Alpha',0.025);
%          %[pcc.e12(subn,:)',in.elr(subn,:)',out.e12(subn,:)',out.elr(subn,:)'],....
% 
%     [h4(nlags),pvalue4(nlags),stat4(nlags),cvalue4(nlags)] = gctest(pcc.elr(subn,:)',in.e12(subn,:)',...
%         'NumLags', nlags, 'Trend',1, 'Alpha',0.025);
%             %[pcc.e12(subn,:)',in.elr(subn,:)',out.e12(subn,:)',out.elr(subn,:)'],....
% 
%     
%     
    
end

plot(stat,'r-','lineWidth',2)
hold on
plot(find(pvalue<0.025), stat(pvalue<0.025),'ro','lineWidth',2)

plot(stat2,'b-','lineWidth',2)
plot(find(pvalue2<0.025), stat(pvalue2<0.025),'bo','lineWidth',2)

legend('in.c12->pcc.clr','sig','pcc.clr->in.c12','sig')
% plot(stat3,'r:','lineWidth',1.5)
% hold on
% plot(find(pvalue3<0.025), stat(pvalue3<0.025),'ro','lineWidth',1.5)
% 
% plot(stat4,'b:','lineWidth',1.5)
% plot(find(pvalue4<0.025), stat(pvalue4<0.025),'bo','lineWidth',1.5)
