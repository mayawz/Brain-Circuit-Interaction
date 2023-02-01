% stagops_behavior
% Maya Zhe Wang
% last used: nov 12 2019

% double check all var export and manipulation
% combine behavioral data from both subjects
% include all trials not cutting safe ones
% only trials paired with the recording (1 session from P and 1 from S)
% wil add RSC-OFC11 session from P later

clear all; close all; clc
fpath='/Users/mwang/Google Drive/01Data/00 w_Anatomy/Gambling/Nov2019Wrapped/psth/combined/';
spath='/Users/mwang/Google Drive/01Data/00 w_Anatomy/Gambling/Nov2019Results/';
cd(fpath);
fname='PS.StagOps.PCC_atOFFER1.mat';
load([fpath fname])
d1=data;
clear data
cutSafe=0;
%%
for sesNum=1:2
    
    switch sesNum
        case 1
            i=19;
            
            c=d1(i);
            
            if cutSafe
                temp=find(c.vars(:,9)~=0);% 0:Safe 1:Lose 2:Win
                psth=c.psth(temp,:);
                vars=c.vars(temp,:);
            else
                psth=c.psth;
                vars=c.vars;
            end
            
            % NORMALIZE
            psth=psth-mean(psth(:));
            psth=psth./var(psth(:));
            
            lp=vars(:,3);% prob of left opt
            temp=vars(:,6);% magnitude of left opt
            templ=temp+1;
            rwdsize=[165 240 125];
            lm=rwdsize(templ)';
            clear temp;
            
            rp=vars(:,4);% prob of right opt
            temp=vars(:,7);% magnitude of right opt
            tempr=temp+1;
            rm=rwdsize(tempr)';
            clear temp;
            
            firstAppear=vars(:,5);% L=1; R=2
            choice=vars(:,8); % 1:Left 2:Right
            choLR=choice; % 1:Left 2:Right
            choLR(choLR==2)=-1;% 1=left -1=right
            outcome=vars(:,9);% 0:Safe 1:Lose 2:Win
            outcome(outcome==1)=0;
            outcome(outcome==2)=1;
            
            lFirstTrl=find(firstAppear==1);
            rFirstTrl=find(firstAppear==2);
            
            
            p1=zeros([length(vars) 1]);
            p2=zeros([length(vars) 1]);
            m1=zeros([length(vars) 1]);
            m2=zeros([length(vars) 1]);
            
            
            p1(lFirstTrl)=lp(lFirstTrl);
            p1(rFirstTrl)=rp(rFirstTrl);
            % p1=normalize(p1,'center','mean');
            
            
            p2(lFirstTrl)=rp(lFirstTrl);
            p2(rFirstTrl)=lp(rFirstTrl);
            
            m1(lFirstTrl)=lm(lFirstTrl);
            m1(rFirstTrl)=rm(rFirstTrl);
            
            m2(lFirstTrl)=rm(lFirstTrl);
            m2(rFirstTrl)=lm(rFirstTrl);
            
            ev1=p1.*m1;
            ev2=p2.*m2;
            
            choOpt=[];
            choOpt(1:length(ev1),1)=2;
            choOpt(choice==firstAppear)=1;
            
            lev=lp.*lm;
            rev=rp.*rm;
            lCorr=intersect(find(lev>rev),find(choice==1));
            rCorr=intersect(find(lev<rev),find(choice==2));
            correctTrls=sort([lCorr; rCorr],'ascend');
            numCorrect=length(correctTrls);
            totalTrl=length(vars(:,1));
            
            receivedRwd=outcome;
            receivedRwd(choOpt==1)=m1(choOpt==1).*outcome(choOpt==1);
            receivedRwd(choOpt==2)=m2(choOpt==2).*outcome(choOpt==2);
            
            
            NumCorrect=numCorrect;
            EV1=ev1;
            EV2=ev2;
            ChoOpt=choOpt;
            ChoLR=choLR;
            RwdOtc=receivedRwd;
            
        case 2
            clearvars -except fpath spath fname cutSafe sesNum d1 NumCorrect...
                EV1 EV2 ChoOpt RwdOtc ChoOpt ChoLR
            
            
            c=d1(end);
            
            if cutSafe
                temp=find(c.vars(:,9)~=0);% 0:Safe 1:Lose 2:Win
                psth=c.psth(temp,:);
                vars=c.vars(temp,:);
            else
                psth=c.psth;
                vars=c.vars;
            end
            
            % NORMALIZE
            psth=psth-mean(psth(:));
            psth=psth./var(psth(:));
            
            lp=vars(:,3);% prob of left opt
            temp=vars(:,6);% magnitude of left opt
            templ=temp+1;
            rwdsize=[165 240 125];
            lm=rwdsize(templ)';
            clear temp;
            
            rp=vars(:,4);% prob of right opt
            temp=vars(:,7);% magnitude of right opt
            tempr=temp+1;
            rm=rwdsize(tempr)';
            clear temp;
            
            firstAppear=vars(:,5);% L=1; R=2
            choice=vars(:,8); % 1:Left 2:Right
            choLR=choice; % 1:Left 2:Right
            choLR(choLR==2)=-1;% 1=left -1=right
            outcome=vars(:,9);% 0:Safe 1:Lose 2:Win
            outcome(outcome==1)=0;
            outcome(outcome==2)=1;
            
            lFirstTrl=find(firstAppear==1);
            rFirstTrl=find(firstAppear==2);
            
            
            p1=zeros([length(vars) 1]);
            p2=zeros([length(vars) 1]);
            m1=zeros([length(vars) 1]);
            m2=zeros([length(vars) 1]);
            
            
            p1(lFirstTrl)=lp(lFirstTrl);
            p1(rFirstTrl)=rp(rFirstTrl);
            % p1=normalize(p1,'center','mean');
            
            
            p2(lFirstTrl)=rp(lFirstTrl);
            p2(rFirstTrl)=lp(rFirstTrl);
            
            m1(lFirstTrl)=lm(lFirstTrl);
            m1(rFirstTrl)=rm(rFirstTrl);
            
            m2(lFirstTrl)=rm(lFirstTrl);
            m2(rFirstTrl)=lm(rFirstTrl);
            
            ev1=p1.*m1;
            ev2=p2.*m2;
            
            choOpt=[];
            choOpt(1:length(ev1),1)=2;
            choOpt(choice==firstAppear)=1;
            
            lev=lp.*lm;
            rev=rp.*rm;
            lCorr=intersect(find(lev>=rev),find(choice==1));
            rCorr=intersect(find(lev<=rev),find(choice==2));
            correctTrls=sort([lCorr; rCorr],'ascend');
            numCorrect=length(correctTrls);
            totalTrl=length(vars(:,1));
            
            receivedRwd=outcome;
            receivedRwd(choOpt==1)=m1(choOpt==1).*outcome(choOpt==1);
            receivedRwd(choOpt==2)=m2(choOpt==2).*outcome(choOpt==2);
            
            
            NumCorrect=[NumCorrect; numCorrect];
            EV1=[EV1; ev1];
            EV2=[EV2; ev2];
            ChoOpt=[ChoOpt; choOpt];
            ChoLR=[ChoLR; choLR];
            RwdOtc=[RwdOtc; receivedRwd];
            
            
    end
    size(ChoOpt)
end

%%
save([fpath fname(1:11) 'Bhv.mat'], 'NumCorrect','EV1','EV2','ChoOpt','ChoLR','RwdOtc');
cd(fpath)

%% 
clear all; close all; clc
fpath='/Users/mwang/Google Drive/01Data/00 w_Anatomy/Gambling/Nov2019Wrapped/psth/combined/';
spath='/Users/mwang/Google Drive/01Data/00 w_Anatomy/Gambling/Nov2019Results/';
cd(fpath);
load('PS.StagOps.Bhv.mat')
%% prob of choosing the larger reward

clc
numCorrect=sum(NumCorrect) %all:848
totalTrl=length(ChoLR) % all: 1160
prop_correct=numCorrect/totalTrl % all: 


%% prob of choosing the left offer

clc
numL=sum(ChoLR==1) % 285
prop_L=numL/totalTrl % 0.4364


%% prob of choosing the first offer

clc
numCho1st=sum(ChoOpt==1) %544 
prop_cho1st=numCho1st/totalTrl % 46.90%

%% prob of choosing offer1 as a value difference of offer1 and offer2

valDiff=EV1-EV2;
nChoOpt1=sum(ChoOpt==1);
ChoOpt1=(ChoOpt==1);

[logitCoef,dev,stats1] = glmfit(valDiff,ChoOpt1,'binomial','logit');
logitFit = glmval(logitCoef,valDiff,'logit');
% [yhat,dylo,dyhi] = glmval(logitCoef,0,'logit',stats1);
% plot(valDiff,logitFit*100,'bo');
% hold on
%
[failedPred,dlo,dhi] = glmval(logitCoef,valDiff,'logit',stats1,.95,100);
errorbar(valDiff,failedPred,dlo,dhi,'o');

xlabel('EV1 minus EV2')
ylabel('Prob Choose Offer1')
hline(50)
vline(0)

%% ################ P ################ trial num 653
clear all; close all; clc
fpath='/Users/mwang/Google Drive/01Data/00 w_Anatomy/Gambling/Nov2019Wrapped/psth/combined/';
spath='/Users/mwang/Google Drive/01Data/00 w_Anatomy/Gambling/Nov2019Results/';
cd(fpath);

% fname1='PS.StagOps.OFCin_atOFFER1.mat';
% fname2='PS.StagOps.OFCin_atCHOICE.mat';

% fname1='PS.StagOps.OFCout_atOFFER1.mat';
% fname2='PS.StagOps.OFCout_atCHOICE.mat';

fname='PS.StagOps.PCC_atOFFER1.mat';
% fname2='PS.StagOps.PCC_atCHOICE.mat';

% fpath='/Users/mwang/Google Drive/01Data/05 StagOps1113/';
% fname='P.StagOps13.mat';
% i=68;

load([fpath fname])
d1=data;
clear data
cutSafe=0;
% opt1=301:400;opt2=401:500; choice=501:650; outcome=651:730;

%% 

i=19;

c=d1(i);

if cutSafe
    temp=find(c.vars(:,9)~=0);% 0:Safe 1:Lose 2:Win
    psth=c.psth(temp,:);
    vars=c.vars(temp,:);
else
    psth=c.psth;
    vars=c.vars;
end

% NORMALIZE
psth=psth-mean(psth(:));
psth=psth./var(psth(:));

lp=vars(:,3);% prob of left opt
temp=vars(:,6);% magnitude of left opt
templ=temp+1;
rwdsize=[165 240 125];
lm=rwdsize(templ)';
clear temp;

rp=vars(:,4);% prob of right opt
temp=vars(:,7);% magnitude of right opt
tempr=temp+1;
rm=rwdsize(tempr)';
clear temp;

firstAppear=vars(:,5);% L=1; R=2
choice=vars(:,8); % 1:Left 2:Right
choLR=choice; % 1:Left 2:Right
choLR(choLR==2)=-1;% 1=left -1=right
outcome=vars(:,9);% 0:Safe 1:Lose 2:Win
outcome(outcome==1)=0;
outcome(outcome==2)=1;

lFirstTrl=find(firstAppear==1);
rFirstTrl=find(firstAppear==2);


p1=zeros([length(vars) 1]);
p2=zeros([length(vars) 1]);
m1=zeros([length(vars) 1]);
m2=zeros([length(vars) 1]);


p1(lFirstTrl)=lp(lFirstTrl);
p1(rFirstTrl)=rp(rFirstTrl);
% p1=normalize(p1,'center','mean');


p2(lFirstTrl)=rp(lFirstTrl);
p2(rFirstTrl)=lp(rFirstTrl);

m1(lFirstTrl)=lm(lFirstTrl);
m1(rFirstTrl)=rm(rFirstTrl);

m2(lFirstTrl)=rm(lFirstTrl);
m2(rFirstTrl)=lm(rFirstTrl);

ev1=p1.*m1;
ev2=p2.*m2;

choOpt=[];
choOpt(1:length(ev1),1)=2;
choOpt(choice==firstAppear)=1;

winTrl=find(outcome~=1);
cho1=find(firstAppear==choice);
rej1=find(firstAppear~=choice);

% checks
%     mean(firstAppear==2)
%     length(cho1Trl)/length(firstAppear)
%     intersect(cho1,rej1) %should be empty

lev=lp.*lm;
rev=rp.*rm;
lCorr=intersect(find(lev>rev),find(choice==1));
rCorr=intersect(find(lev<rev),find(choice==2));
correctTrls=sort([lCorr; rCorr],'ascend');
numCorrect=length(correctTrls);
totalTrl=length(vars(:,1));

receivedRwd=outcome;
receivedRwd(choOpt==1)=m1(choOpt==1).*outcome(choOpt==1);
receivedRwd(choOpt==2)=m2(choOpt==2).*outcome(choOpt==2);

cho1Trl=intersect(cho1,winTrl);
rej1Trl=intersect(rej1,winTrl);

%   normalize(A,2,'center','mean')
%   normalize each row across the 2nd dim col

% keyboard
% m1=zscore(m1);
% m2=zscore(m2);
% p1=zscore(p1);
% p2=zscore(p2);
% ev1=zscore(ev1);
% ev2=zscore(ev2);
% choOpt=zscore(choOpt);
% rwdOtc=zscore(receivedRwd);

%% prob of choosing the larger reward

numCorrect_P=numCorrect %all:479
totalTrl_P=totalTrl % all: 653
prop_correct_P=numCorrect_P/totalTrl_P %=0.7360 (cutSafe) % 0.7335 (all)



%% prob of choosing the left offer

numL_P=sum(choLR==1) % 285
prop_L_P=numL_P/totalTrl_P % 0.4364


%% prob of choosing the first offer

numCho1st_P=sum(choOpt==1) %295
prop_cho1st_P=numCho1st_P/totalTrl_P % 0.4518

%% prob of choosing offer1 as a value difference of offer1 and offer2

valDiff=ev1-ev2;
nChoOpt1=sum(choOpt==1);
ChoOpt1=(choOpt==1);

[logitCoef,dev,stats1] = glmfit(valDiff,ChoOpt1,'binomial','logit');
logitFit = glmval(logitCoef,valDiff,'logit');
% [yhat,dylo,dyhi] = glmval(logitCoef,0,'logit',stats1);
% plot(valDiff,logitFit*100,'bo');
% hold on
%
[failedPred,dlo,dhi] = glmval(logitCoef,valDiff,'logit',stats1,.95,100);
errorbar(valDiff,failedPred,dlo,dhi,'o');

xlabel('EV1 minus EV2')
ylabel('Prob Choose Offer1')
hline(50)
vline(0)

%% ################ S ################


clear all; close all; clc
fpath='/Users/mwang/Google Drive/01Data/00 w_Anatomy/Gambling/Nov2019Wrapped/psth/combined/';
spath='/Users/mwang/Google Drive/01Data/00 w_Anatomy/Gambling/Nov2019Results/';
cd(fpath);

% fname1='PS.StagOps.OFCin_atOFFER1.mat';
% fname2='PS.StagOps.OFCin_atCHOICE.mat';

% fname1='PS.StagOps.OFCout_atOFFER1.mat';
% fname2='PS.StagOps.OFCout_atCHOICE.mat';

fname='PS.StagOps.PCC_atOFFER1.mat';
% fname2='PS.StagOps.PCC_atCHOICE.mat';

% fpath='/Users/mwang/Google Drive/01Data/05 StagOps1113/';
% fname='P.StagOps13.mat';
% i=68;

load([fpath fname])
d1=data;
clear data
cutSafe=0;
% opt1=301:400;opt2=401:500; choice=501:650; outcome=651:730;

%% subject S trial num 507

i=100;
c=d1(i);

if cutSafe
    temp=find(c.vars(:,9)~=0);% 0:Safe 1:Lose 2:Win
    psth=c.psth(temp,:);
    vars=c.vars(temp,:);
else
    psth=c.psth;
    vars=c.vars;
end

% NORMALIZE
psth=psth-mean(psth(:));
psth=psth./var(psth(:));

lp=vars(:,3);% prob of left opt
temp=vars(:,6);% magnitude of left opt
templ=temp+1;
rwdsize=[165 240 125];
lm=rwdsize(templ)';
clear temp;

rp=vars(:,4);% prob of right opt
temp=vars(:,7);% magnitude of right opt
tempr=temp+1;
rm=rwdsize(tempr)';
clear temp;

firstAppear=vars(:,5);% L=1; R=2
choice=vars(:,8); % 1:Left 2:Right
choLR=choice; % 1:Left 2:Right
choLR(choLR==2)=-1;% 1=left -1=right
outcome=vars(:,9);% 0:Safe 1:Lose 2:Win
outcome(outcome==1)=0;
outcome(outcome==2)=1;

lFirstTrl=find(firstAppear==1);
rFirstTrl=find(firstAppear==2);


p1=zeros([length(vars) 1]);
p2=zeros([length(vars) 1]);
m1=zeros([length(vars) 1]);
m2=zeros([length(vars) 1]);


p1(lFirstTrl)=lp(lFirstTrl);
p1(rFirstTrl)=rp(rFirstTrl);
% p1=normalize(p1,'center','mean');


p2(lFirstTrl)=rp(lFirstTrl);
p2(rFirstTrl)=lp(rFirstTrl);

m1(lFirstTrl)=lm(lFirstTrl);
m1(rFirstTrl)=rm(rFirstTrl);

m2(lFirstTrl)=rm(lFirstTrl);
m2(rFirstTrl)=lm(rFirstTrl);

ev1=p1.*m1;
ev2=p2.*m2;

choOpt=[];
choOpt(1:length(ev1),1)=2;
choOpt(choice==firstAppear)=1;

winTrl=find(outcome~=1);
cho1=find(firstAppear==choice);
rej1=find(firstAppear~=choice);

% checks
%     mean(firstAppear==2)
%     length(cho1Trl)/length(firstAppear)
%     intersect(cho1,rej1) %should be empty

lev=lp.*lm;
rev=rp.*rm;
lCorr=intersect(find(lev>rev),find(choice==1));
rCorr=intersect(find(lev<rev),find(choice==2));
correctTrls=sort([lCorr; rCorr],'ascend');
numCorrect=length(correctTrls);
totalTrl=length(vars(:,1));

receivedRwd=outcome;
receivedRwd(choOpt==1)=m1(choOpt==1).*outcome(choOpt==1);
receivedRwd(choOpt==2)=m2(choOpt==2).*outcome(choOpt==2);

cho1Trl=intersect(cho1,winTrl);
rej1Trl=intersect(rej1,winTrl);

%   normalize(A,2,'center','mean')
%   normalize each row across the 2nd dim col

% keyboard
% m1=zscore(m1);
% m2=zscore(m2);
% p1=zscore(p1);
% p2=zscore(p2);
% ev1=zscore(ev1);
% ev2=zscore(ev2);
% choOpt=zscore(choOpt);
% rwdOtc=zscore(receivedRwd);

%% prob of choosing the larger reward

numCorrect_S=numCorrect % 345 (cutSafe) 367
totalTrl_S=totalTrl % 467 (cutSafe) 507
prop_correct_S=numCorrect_S/totalTrl_S %0.7388 (cutSafe) % 0.7239 (all)

%% prob of choosing the left offer

numL_S=sum(choLR==1) % 224
prop_L_S=numL_S/totalTrl_S % 0.4418


%% prob of choosing the first offer

numCho1st_S=sum(choOpt==1) % 249
prop_cho1st_S=numCho1st_S/totalTrl_S % 0.4911

%% prob of choosing offer1 as a value difference of offer1 and offer2

valDiff=ev1-ev2;
nChoOpt1=sum(choOpt==1);
ChoOpt1=(choOpt==1);

[logitCoef,dev,stats1] = glmfit(valDiff,ChoOpt1,'binomial','logit');
logitFit = glmval(logitCoef,valDiff,'logit');
% [yhat,dylo,dyhi] = glmval(logitCoef,0,'logit',stats1);
% plot(valDiff,logitFit*100,'bo');
% hold on
%
[failedPred,dlo,dhi] = glmval(logitCoef,valDiff,'logit',stats1,.95,100);
errorbar(valDiff,failedPred,dlo,dhi,'o');

xlabel('EV1 minus EV2')
ylabel('Prob Choose Offer1')
hline(50)
vline(0)
















