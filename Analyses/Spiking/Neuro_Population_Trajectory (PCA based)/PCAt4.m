% PCAt4
% last used May 4 2020
% may the 4th be with you
% calculate distance for trial by trial data
% use data from .m files PrepPCAt2020 & PCAt1


%% ################## Across 3 Regions ##################
clear all; close all; clc

dpath='/Users/mwang/Google Drive/01Data/00 w_Anatomy/2020final/results/PCAt/';
fpath='/Users/mwang/Google Drive/01Data/00 w_Anatomy/2020final/code/Trajectory/';
spath=dpath;
addpath(genpath(fpath))
cd(dpath)

stepSize=5;
binStarts=251:stepSize:750-stepSize+1;

% subject pumbaa=1 spock=2
subn=2;

load('OFCinDistTbT.mat');
%in=td(subn);
in=td;
clear td
load('OFCoutDistTbT.mat');
% out=td(subn);
out=td;
clear td
load('PCCDistTbT.mat');
% pcc=td(subn);
pcc=td;
clear td

%% ############## distance 3 areas -- correct ##############
% combine subjects
close all; clc; clear disIn disPcc disOut

% disIn=[in(1).c12.Dst.aD12;in(1).c12.Dst.aD21;...
%     in(2).c12.Dst.aD12;in(2).c12.Dst.aD21];
% disPcc=[pcc(1).c12.Dst.aD12;pcc(1).c12.Dst.aD21;...
%     pcc(2).c12.Dst.aD12;pcc(2).c12.Dst.aD21];
% disOut=[out(1).c12.Dst.aD12;out(1).c12.Dst.aD21;...
%     out(2).c12.Dst.aD12;out(2).c12.Dst.aD21];
% 
% disIn=[in(1).clr.Dst.aD12;in(1).clr.Dst.aD21;...
%     in(2).clr.Dst.aD12;in(2).clr.Dst.aD21];
% disPcc=[pcc(1).clr.Dst.aD12;pcc(1).clr.Dst.aD21;...
%     pcc(2).clr.Dst.aD12;pcc(2).clr.Dst.aD21];
% disOut=[out(1).clr.Dst.aD12;out(1).clr.Dst.aD21;...
%     out(2).clr.Dst.aD12;out(2).clr.Dst.aD21];

disIn=[in(1).cev1.Dst.aD12;in(1).cev1.Dst.aD21;...
    in(2).cev1.Dst.aD12;in(2).cev1.Dst.aD21];
disPcc=[pcc(1).cev1.Dst.aD12;pcc(1).cev1.Dst.aD21;...
    pcc(2).cev1.Dst.aD12;pcc(2).cev1.Dst.aD21];
disOut=[out(1).cev1.Dst.aD12;out(1).cev1.Dst.aD21;...
    out(2).cev1.Dst.aD12;out(2).cev1.Dst.aD21];

% Over all
close all; clear dIn dOut dPcc
dIn=nanmean(disIn,1)';
dPcc=nanmean(disPcc,1)';
dOut=nanmean(disOut,1)';

smo=5; yLims=[2 16];
DistByTimePlot3_e2(disIn,disPcc,disOut,smo,yLims);

%% distance 3 areas -- error
% combine subjects
close all; clc; clear disIn disPcc disOut

% disIn=[in(1).e12.Dst.aD12;in(1).e12.Dst.aD21;...
%     in(2).e12.Dst.aD12;in(2).e12.Dst.aD21];
% disPcc=[pcc(1).e12.Dst.aD12;pcc(1).e12.Dst.aD21;...
%     pcc(2).e12.Dst.aD12;pcc(2).e12.Dst.aD21];
% disOut=[out(1).e12.Dst.aD12;out(1).e12.Dst.aD21;...
%     out(2).e12.Dst.aD12;out(2).e12.Dst.aD21];
% 
% % disIn=[in(1).elr.Dst.aD12;in(1).elr.Dst.aD21;...
% %     in(2).elr.Dst.aD12;in(2).elr.Dst.aD21];
% % disPcc=[pcc(1).elr.Dst.aD12;pcc(1).elr.Dst.aD21;...
% %     pcc(2).elr.Dst.aD12;pcc(2).elr.Dst.aD21];
% % disOut=[out(1).elr.Dst.aD12;out(1).elr.Dst.aD21;...
% %     out(2).elr.Dst.aD12;out(2).elr.Dst.aD21];

disIn=[in(1).eev1.Dst.aD12;in(1).eev1.Dst.aD21;...
    in(2).eev1.Dst.aD12;in(2).eev1.Dst.aD21];
disPcc=[pcc(1).eev1.Dst.aD12;pcc(1).eev1.Dst.aD21;...
    pcc(2).eev1.Dst.aD12;pcc(2).eev1.Dst.aD21];
disOut=[out(1).eev1.Dst.aD12;out(1).eev1.Dst.aD21;...
    out(2).eev1.Dst.aD12;out(2).eev1.Dst.aD21];

% Over all
close all; clear dIn dOut dPcc
dIn=nanmean(disIn,1)';
dPcc=nanmean(disPcc,1)';
dOut=nanmean(disOut,1)';

smo=5; yLims=[2 16];
DistByTimePlot3_e2(disIn,disPcc,disOut,smo,yLims);

%% ############## sig test 3 areas ##############
%% Kruskal Wallis
clc; close all
usebins=1:100;
[p,tbl,stats]=kruskalwallis([dIn(usebins),dPcc(usebins),dOut(usebins)]);

figure(3)
multcompare(stats)
%% Kolmogorov-Smirnov
clc;close all
[h,p,ks2stat] = kstest2(dIn,dPcc,'Tail','unequal') 
figure(1)
histogram(dPcc,ceil(length(dPcc)*0.4));
hold on
histogram(dIn,ceil(length(dIn)*0.4));
legend('dPcc','dIn')
figure(2)
ecdf(dIn,'Bounds','on')
hold on
cdfplot(dPcc)
grid on
legend('dIn','','','dPcc')

%% KS alternative Hypo: 'smaller' 'larger'
[h,p,ks2stat] = kstest2(dIn,dPcc,'Tail','larger') 
%% Wilcoxon signed rank test
clc;close all
[p,h,stats] = signrank(dIn,dPcc)

%% dispersion 3 areas  -- correct
clc; clear tin tout tpcc
tin=in.c12.Dsp;
tout=out.c12.Dsp;
tpcc=pcc.c12.Dsp;
smo=1;
DispByTimePlot3(tin,tpcc,tout,smo);

%% dispersion 3 areas  -- error
clc; clear tin tout tpcc
tin=in.e12.Dsp;
tout=out.e12.Dsp;
tpcc=pcc.e12.Dsp;
smo=1;
DispByTimePlot3(tin,tpcc,tout,smo);


%% ################## correct vs error ##################
clear all; close all; clc

dpath='/Users/mwang/Google Drive/01Data/00 w_Anatomy/2020final/results/PCAt/';
fpath='/Users/mwang/Google Drive/01Data/00 w_Anatomy/2020final/code/Trajectory/';
spath=dpath;
addpath(genpath(fpath))
cd(dpath)

stepSize=5;
binStarts=251:stepSize:750-stepSize+1;

% subject pumbaa=1 spock=2
subn=1;

load('OFCinDistTbT.mat');
%in=td(subn);
in=td;
clear td
load('OFCoutDistTbT.mat');
% out=td(subn);
out=td;
clear td
load('PCCDistTbT.mat');
% pcc=td(subn);
pcc=td;
clear td

clear tmpd
tmpd=pcc;
% cdd=[tmpd(1).c12.Dst.aD12;tmpd(1).c12.Dst.aD21;...
%     tmpd(2).c12.Dst.aD12;tmpd(2).c12.Dst.aD21];
% edd=[tmpd(1).e12.Dst.aD12;tmpd(1).e12.Dst.aD21;...
%     tmpd(2).e12.Dst.aD12;tmpd(2).e12.Dst.aD21];

% cdd=[tmpd(1).clr.Dst.aD12;tmpd(1).clr.Dst.aD21;...
%     tmpd(2).clr.Dst.aD12;tmpd(2).clr.Dst.aD21];
% edd=[tmpd(1).elr.Dst.aD12;tmpd(1).elr.Dst.aD21;...
%     tmpd(2).elr.Dst.aD12;tmpd(2).elr.Dst.aD21];

cdd=[tmpd(1).cev1.Dst.aD12;tmpd(1).cev1.Dst.aD21;...
    tmpd(2).cev1.Dst.aD12;tmpd(2).cev1.Dst.aD21];
edd=[tmpd(1).eev1.Dst.aD12;tmpd(1).eev1.Dst.aD21;...
    tmpd(2).eev1.Dst.aD12;tmpd(2).eev1.Dst.aD21];
% Over all
clear cD eD
cD=nanmean(cdd,1)';
eD=nanmean(edd,1)';

%% ############## sig test correct vs error ##############
%% Kruskal Wallis
clc;close all
[kw.p,kw.tbl,kw.stats]=kruskalwallis([cD,eD])

%% Wilcoxon signed rank test
clc; close all
[wsr.p,~,wsr.stats] = signrank(cD,eD)



%% Two-sample Kolmogorov-Smirnov
clc
[h,p,ks2stat] = kstest2(cD,eD,'Tail','unequal') 
figure(1)
histogram(cD,ceil(length(cD)*0.4));
hold on
histogram(eD,ceil(length(eD)*0.4));
legend('cD','eD')

figure(2)
ecdf(cD,'Bounds','on')
hold on
cdfplot(eD)
grid on
legend('cD','','','eD')
%% KS alternative Hypo: 'smaller' 'larger'
[h,p,ks2stat] = kstest2(cD,eD,'Tail','smaller') 

%% dispersion correct vs error
clc; clear tmpd cdd edd
tmpd=pcc;

% cdd=[tmpd(1).c12.Dsp.cnd1;tmpd(1).c12.Dsp.cnd2;...
%     tmpd(2).c12.Dsp.cnd1;tmpd(2).c12.Dsp.cnd2];
% edd=[tmpd(1).e12.Dsp.cnd1;tmpd(1).e12.Dsp.cnd2;...
%     tmpd(2).e12.Dsp.cnd1;tmpd(2).e12.Dsp.cnd2];

% cdd=[tmpd(1).clr.Dsp.cnd1;tmpd(1).clr.Dsp.cnd2;...
%     tmpd(2).clr.Dsp.cnd1;tmpd(2).clr.Dsp.cnd2];
% edd=[tmpd(1).elr.Dsp.cnd1;tmpd(1).elr.Dsp.cnd2;...
%     tmpd(2).elr.Dsp.cnd1;tmpd(2).elr.Dsp.cnd2];


cdd=[tmpd(1).cev1.Dsp.cnd1;tmpd(1).cev1.Dsp.cnd2;...
    tmpd(2).cev1.Dsp.cnd1;tmpd(2).cev1.Dsp.cnd2];
edd=[tmpd(1).eev1.Dsp.cnd1;tmpd(1).eev1.Dsp.cnd2;...
    tmpd(2).eev1.Dsp.cnd1;tmpd(2).eev1.Dsp.cnd2];

% Over all
clear cD eD
cD=nanmean(cdd,1)';
eD=nanmean(edd,1)';

clc;close all
[kw.p,kw.tbl,kw.stats]=kruskalwallis([cD,eD])

%% ######################## Granger Causality ########################
% ############################ both subs ################################
% when check "feedback" granger causality, reduce significance to half
% 'Alpha', 0.025 then test A-->B and then B-->A

clear all; close all; clc

dpath='/Users/mwang/Google Drive/01Data/00 w_Anatomy/2020final/results/PCAt/';
fpath='/Users/mwang/Google Drive/01Data/00 w_Anatomy/2020final/code/Trajectory/';
spath=dpath;
addpath(genpath(fpath))
cd(dpath)

stepSize=5;
binStarts=251:stepSize:750-stepSize+1;

% subject pumbaa=1 spock=2
% subn=1;

load('OFCinDistTbT.mat');
%in=td(subn);
in=td;
clear td
load('OFCoutDistTbT.mat');
% out=td(subn);
out=td;
clear td
load('PCCDistTbT.mat');
% pcc=td(subn);
pcc=td;
clear td

% merge
in12ct=[in(1).c12.Dst.aD12;in(1).c12.Dst.aD21;...
    in(2).c12.Dst.aD12;in(2).c12.Dst.aD21];
in12et=[in(1).e12.Dst.aD12;in(1).e12.Dst.aD21;...
    in(2).e12.Dst.aD12;in(2).e12.Dst.aD21];

inlrct=[in(1).clr.Dst.aD12;in(1).clr.Dst.aD21;...
    in(2).clr.Dst.aD12;in(2).clr.Dst.aD21];
inlret=[in(1).elr.Dst.aD12;in(1).elr.Dst.aD21;...
    in(2).elr.Dst.aD12;in(2).elr.Dst.aD21];

pcc12ct=[pcc(1).c12.Dst.aD12;pcc(1).c12.Dst.aD21;...
    pcc(2).c12.Dst.aD12;pcc(2).c12.Dst.aD21];
pcc12et=[pcc(1).e12.Dst.aD12;pcc(1).e12.Dst.aD21;...
    pcc(2).e12.Dst.aD12;pcc(2).e12.Dst.aD21];

pcclrct=[pcc(1).elr.Dst.aD12;pcc(1).elr.Dst.aD21;...
    pcc(2).elr.Dst.aD12;pcc(2).elr.Dst.aD21];
pcclret=[pcc(1).clr.Dst.aD12;pcc(1).clr.Dst.aD21;...
    pcc(2).clr.Dst.aD12;pcc(2).clr.Dst.aD21];

out12ct=[out(1).c12.Dst.aD12;out(1).c12.Dst.aD21;...
    out(2).c12.Dst.aD12;out(2).c12.Dst.aD21];
out12et=[out(1).e12.Dst.aD12;out(1).e12.Dst.aD21;...
    out(2).e12.Dst.aD12;out(2).e12.Dst.aD21];

outlrct=[out(1).clr.Dst.aD12;out(1).clr.Dst.aD21;...
    out(2).clr.Dst.aD12;out(2).clr.Dst.aD21];
outlret=[out(1).elr.Dst.aD12;out(1).elr.Dst.aD21;...
    out(2).elr.Dst.aD12;out(2).elr.Dst.aD21];

% average
in12c=nanmean(in12ct,1)';
in12e=nanmean(in12et,1)';
inlrc=nanmean(inlrct,1)';
inlre=nanmean(inlret,1)';

pcc12c=nanmean(pcc12ct,1)';
pcc12e=nanmean(pcc12et,1)';
pcclrc=nanmean(pcclrct,1)';
pcclre=nanmean(pcclret,1)';

out12c=nanmean(out12ct,1)';
out12e=nanmean(out12et,1)';
outlrc=nanmean(outlrct,1)';
outlre=nanmean(outlret,1)';
% test stationary: augmented Dickey-Fuller test, h = adftest(var,'Model',"ard")
% adftest(in.c12(subn,:)','Model','ard') % stationary 
% adftest(pcc.clr(subn,:)','Model','ard') % stationary 
% adftest(detrend(out.c12(subn,:)'),'Model','ard') % non-stationary 
% detrend: detrend(sdata); % detrended out data is stationary

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

NLAG = numlags(aic == min(aic))

[h,pvalue,stat,cvalue] = gctest(X1,X2,'NumLags', NLAG, 'Trend', 0)

%%
clc; close all
for nlags=1:13 
    
    % dist 12 --> lr 
    [h(nlags),pvalue(nlags),stat(nlags),cvalue(nlags)] = gctest(in12c,pcclrc,...
        [pcc12c,inlrc,out12c,outlrc],....
        'NumLags', nlags, 'Trend',1, 'Alpha',0.025);
    
    [h2(nlags),pvalue2(nlags),stat2(nlags),cvalue2(nlags)] = gctest(in12e,pcclre,...
        [pcc12e,inlre,out12e,outlre],....
        'NumLags', nlags, 'Trend',1, 'Alpha',0.025);
    
    [h3(nlags),pvalue3(nlags),stat3(nlags),cvalue3(nlags)] = gctest(pcclrc,in12c,...
        [pcc12c,inlrc,out12c,outlrc],....
        'NumLags', nlags, 'Trend',1, 'Alpha',0.025);
% 
%     [h4(nlags),pvalue4(nlags),stat4(nlags),cvalue4(nlags)] = gctest(pcc.elr(subn,:)',in.e12(subn,:)',...
%         'NumLags', nlags, 'Trend',1, 'Alpha',0.025);
%             %[pcc.e12(subn,:)',in.elr(subn,:)',out.e12(subn,:)',out.elr(subn,:)'],....
% 
%     
%     
    
end

%%
close all
plot(stat,'r-','lineWidth',2)
hold on
plot(find(pvalue<0.025), stat(pvalue<0.025),'ro','lineWidth',2)

plot(stat2,'b-','lineWidth',2)
plot(find(pvalue2<0.025), stat2(pvalue2<0.025),'bo','lineWidth',2)

legend('in.c12->pcc.clr','sig','pcc.clr->in.c12','sig')
plot(stat3,'r:','lineWidth',1.5)
hold on
plot(find(pvalue3<0.025), stat3(pvalue3<0.025),'ro','lineWidth',1.5)
% 
% plot(stat4,'b:','lineWidth',1.5)
% plot(find(pvalue4<0.025), stat(pvalue4<0.025),'bo','lineWidth',1.5)


%%
stepSize=5;
binStarts=251:stepSize:750-stepSize+1;
binStarts(11)


%% ######################## Granger Causality ########################
% ############################ single subs ################################
% when check "feedback" granger causality, reduce significance to half
% 'Alpha', 0.025 then test A-->B and then B-->A

clear all; close all; clc

dpath='/Users/mwang/Google Drive/01Data/00 w_Anatomy/2020final/results/PCAt/';
fpath='/Users/mwang/Google Drive/01Data/00 w_Anatomy/2020final/code/Trajectory/';
spath=dpath;
addpath(genpath(fpath))
cd(dpath)

stepSize=5;
binStarts=251:stepSize:750-stepSize+1;

% subject pumbaa=1 spock=2
subn=2;

load('OFCinDistTbT.mat');
in=td(subn);
% in=td;
clear td
load('OFCoutDistTbT.mat');
out=td(subn);
% out=td;
clear td
load('PCCDistTbT.mat');
pcc=td(subn);
% pcc=td;
clear td

% merge
in12ct=[in.c12.Dst.aD12;in.c12.Dst.aD21];
in12et=[in.e12.Dst.aD12;in.e12.Dst.aD21];

inlrct=[in.clr.Dst.aD12;in.clr.Dst.aD21];
inlret=[in.elr.Dst.aD12;in.elr.Dst.aD21];

pcc12ct=[pcc.c12.Dst.aD12;pcc.c12.Dst.aD21];
pcc12et=[pcc.e12.Dst.aD12;pcc.e12.Dst.aD21];

pcclrct=[pcc.elr.Dst.aD12;pcc.elr.Dst.aD21];
pcclret=[pcc.clr.Dst.aD12;pcc.clr.Dst.aD21];

out12ct=[out.c12.Dst.aD12;out.c12.Dst.aD21;];
out12et=[out.e12.Dst.aD12;out.e12.Dst.aD21];

outlrct=[out.clr.Dst.aD12;out.clr.Dst.aD21];
outlret=[out.elr.Dst.aD12;out.elr.Dst.aD21];

% average
in12c=nanmean(in12ct,1)';
in12e=nanmean(in12et,1)';
inlrc=nanmean(inlrct,1)';
inlre=nanmean(inlret,1)';

pcc12c=nanmean(pcc12ct,1)';
pcc12e=nanmean(pcc12et,1)';
pcclrc=nanmean(pcclrct,1)';
pcclre=nanmean(pcclret,1)';

out12c=nanmean(out12ct,1)';
out12e=nanmean(out12et,1)';
outlrc=nanmean(outlrct,1)';
outlre=nanmean(outlret,1)';
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

NLAG = numlags(aic == min(aic))

[h,pvalue,stat,cvalue] = gctest(X1,X2,'NumLags', NLAG, 'Trend', 0)

%%
clc; close all
for nlags=1:13 
    
    % dist 12 --> lr 
    [h(nlags),pvalue(nlags),stat(nlags),cvalue(nlags)] = gctest(in12c,pcclrc,...
        [pcc12c,inlrc,out12c,outlrc],....
        'NumLags', nlags, 'Trend',1, 'Alpha',0.025);
    
    [h2(nlags),pvalue2(nlags),stat2(nlags),cvalue2(nlags)] = gctest(in12e,pcclre,...
        [pcc12e,inlre,out12e,outlre],....
        'NumLags', nlags, 'Trend',1, 'Alpha',0.025);
    
    [h3(nlags),pvalue3(nlags),stat3(nlags),cvalue3(nlags)] = gctest(pcclrc,in12c,...
        [pcc12c,inlrc,out12c,outlrc],....
        'NumLags', nlags, 'Trend',1, 'Alpha',0.025);
% 
%     [h4(nlags),pvalue4(nlags),stat4(nlags),cvalue4(nlags)] = gctest(pcc.elr(subn,:)',in.e12(subn,:)',...
%         'NumLags', nlags, 'Trend',1, 'Alpha',0.025);
%             %[pcc.e12(subn,:)',in.elr(subn,:)',out.e12(subn,:)',out.elr(subn,:)'],....
% 
%     
%     
    
end

%%
close all
plot(stat,'r-','lineWidth',2)
hold on
plot(find(pvalue<0.025), stat(pvalue<0.025),'ro','lineWidth',2)

plot(stat2,'b-','lineWidth',2)
plot(find(pvalue2<0.025), stat2(pvalue2<0.025),'bo','lineWidth',2)

legend('in.c12->pcc.clr','sig','pcc.clr->in.c12','sig')
plot(stat3,'r:','lineWidth',1.5)
hold on
plot(find(pvalue3<0.025), stat3(pvalue3<0.025),'ro','lineWidth',1.5)
% 
% plot(stat4,'b:','lineWidth',1.5)
% plot(find(pvalue4<0.025), stat(pvalue4<0.025),'bo','lineWidth',1.5)


%%
stepSize=5;
binStarts=251:stepSize:750-stepSize+1;
binStarts(11)

%% dist 3 areas - Correct - single subject

clc; clear tin tout tpcc

% subject pumbaa=1 spock=2
subn=2;

% tin=in(subn).c12.Dst;
% tout=out(subn).c12.Dst;
% tpcc=pcc(subn).c12.Dst;

tin=in(subn).clr.Dst;
tout=out(subn).clr.Dst;
tpcc=pcc(subn).clr.Dst;

% tin=in(subn).cev1.Dst;
% tout=out(subn).cev1.Dst;
% tpcc=pcc(subn).cev1.Dst;

smo=5;yLims=[1 9]; % 9 % 30
DistByTimePlot3(tin,tpcc,tout,smo);

% disIn=[tin.aD12;tin.aD21];
% disPcc=[tpcc.aD12;tpcc.aD21];
% disOut=[tout.aD12;tout.aD21];

%% dist 3 areas - Error - single subject

close all;clc; clear tin tout tpcc

tin=in(subn).e12.Dst;
tout=out(subn).e12.Dst;
tpcc=pcc(subn).e12.Dst;

% tin=in(subn).elr.Dst;
% tout=out(subn).elr.Dst;
% tpcc=pcc(subn).elr.Dst;

% tin=in(subn).eev1.Dst;
% tout=out(subn).eev1.Dst;
% tpcc=pcc(subn).eev1.Dst;

smo=5;yLims=[1 9];
DistByTimePlot3(tin,tpcc,tout,smo);

disIn=[tin.aD12;tin.aD21];
disPcc=[tpcc.aD12;tpcc.aD21];
disOut=[tout.aD12;tout.aD21];

% Over all
clear dIn dOut dPcc
dIn=nanmean(disIn,1)';
dPcc=nanmean(disPcc,1)';
dOut=nanmean(disOut,1)';

%% distance correct vs error
clc; clear cdata edata cdd edd cD eD
% cdata=in.c12.Dst;
% edata=in.e12.Dst;
% cdata=in.clr.Dst;
% edata=in.elr.Dst;
% cdata=in.cev1.Dst;
% edata=in.eev1.Dst;

% cdata=pcc.c12.Dst;
% edata=pcc.e12.Dst;
% cdata=pcc.clr.Dst;
% edata=pcc.elr.Dst;
cdata=pcc.cev1.Dst;
edata=pcc.eev1.Dst;

% cdata=out.c12.Dst;
% edata=out.e12.Dst;
% cdata=out.clr.Dst;
% edata=out.elr.Dst;
% cdata=out.cev1.Dst;
% edata=out.eev1.Dst;
smo=1;
DistByTimePlot2(cdata,edata,smo);

% ############ sig test for 2 things ############
cdd=[cdata.aD12;cdata.aD21];
edd=[edata.aD12;edata.aD21];
% Over all
clear cD eD
cD=nanmean(cdd,1)';
eD=nanmean(edd,1)';
% By time point

