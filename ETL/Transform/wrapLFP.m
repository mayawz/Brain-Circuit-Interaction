% wrap LFP
% last used: Sep 8 2019
clear all; close all; clc

% sampling rate 1Ksamples/s = bin size 1ms


% cd('/Users/mwang/Google Drive/01Data/00 w_Anatomy/LFP/OFCin/')
% cd('/Users/mwang/Google Drive/01Data/00 w_Anatomy/LFP/OFCout/')
cd('/Users/mwang/Google Drive/01Data/00 w_Anatomy/LFP/PCC/')

spath='/Users/mwang/Google Drive/01Data/00 w_Anatomy/Gambling/Nov2019Wrapped/forCoherence/LFP/';
%%
% d=load(['/Users/mwang/Google Drive/01Data/00 w_Anatomy/Gambling/Sorted/OFCin/' ...
%     'P20171207_OFC13_Array6_sorted.mat']);
% d=load(['/Users/mwang/Google Drive/01Data/00 w_Anatomy/Gambling/Sorted/OFCin/' ...
%     'S20180731_OFC13_Array4_sorted.mat']);

% d=load(['/Users/mwang/Google Drive/01Data/00 w_Anatomy/Gambling/Sorted/OFCout/' ...
%     'P20171207_OFC13_Array5_sorted.mat']);
% d=load(['/Users/mwang/Google Drive/01Data/00 w_Anatomy/Gambling/Sorted/OFCout/' ...
%     'S20180731_OFC13_Array6_sorted.mat']);
% 
% d=load(['/Users/mwang/Google Drive/01Data/00 w_Anatomy/Gambling/Sorted/PCC/' ...
%     'P20171207_PCC_Array2_sorted.mat']);
d=load(['/Users/mwang/Google Drive/01Data/00 w_Anatomy/Gambling/Sorted/PCC/' ...
    'S20180731_PCC_Array2_sorted.mat']);

Strobed=d.Strobed;
clear d

% fname='P20171207_OFC13_Array6_LFP.mat';
% fname='S20180731_OFC13_Array45_LFP.mat';

% fname='P20171207_OFC13_Array5_LFP.mat';
% fname='S20180731_OFC13_Array67_LFP.mat';


% fname='P20171207_PCC_LFP.mat';
fname='S20180731_PCC_LFP.mat';

load(fname);
% 
% startStrobe=4001;
% align='OFFER1';


startStrobe=4007
align='CHOICE';

% area='OFCin';
% area='OFCout';
area='PCC';


%%
% memo for all the strobes:
% task
% 4001 1st opt appears
% 4002 1st opt disappears
% 4003 2st opt appears
% 4035 2st opt disappears
% 4004 Fixation dot appears
% 4005 3rd op appears
% 4006 3rd op disappears
% 4007 choice onset "go signal"
% 4008 feedback
% 4009 ITI
% 4051 firstly looked at option 1... 4052 firstly looked option 2
% 4061 fixation
% 4062 fixation lost
% 4073 fixed on right

% other vars

%	1 - 2000: trial numbers
%	8201 - 8203 : Number of options

%	30XX: gamble prob for left XX = prob*100
%	33XX: gamble prob for right XX = prob*100
%	35XX: gamble prob for center XX = prob*100

%	12ABC: Option order, ABC=(order(1)*100)+( order(2)*10 )+ order(3) 1:Left 2:Right 3:Center

%	13ABC: ABC=(notBlueOps(1)*100)+(notBlueOps(2)*10)+notBlueOps(3)
%   notBlueOps ABC [position] corresponds to A=left B=right C=central
%   the value on each [notBlueOps ABC position]
%   2:safe gamble, medium reward
%   1 && gamble win: Green, huge reward
%   0 && gamble win: Blue, large reward
%   gamble not win: nothing


%	8001 - 8003 : Choice 1:Left 2:Right 3:Center
%	10001 - 10003: Gamble outcome 0:Safe 1:Lose 2:Win

%   reward size: medium=150ul, large=180ul, huge=210ul, +10500
%   --> medium=10650, large=10680, huge= 10710

%   chance3op=0 no third option
%	7000, 7200, 7400: chance3op*100+7000, chance of a third option, of safe (grey), large (blue) huge (green)
%	20000 - 20004: Location of three options. Not used.


strobeName(:,1)=Strobed(:,2);

temp=find(strobeName==20000);
theCpltEnd=temp(end); % the end of all complete trials
clear temp;

temp=find(strobeName==startStrobe);
theCpltBeginning=temp(1); % the beginning of all complete trials
clear temp;

Strobed=Strobed(theCpltBeginning:theCpltEnd,:);

% get vars of one cell out of one day, one session

% vars = extractVars2(Strobed);
[vars strobesFromVars]= extractVars3(Strobed);

%% align lfp
trialStartStrobe=find(Strobed(:,2)==startStrobe);% define start strobes of each trial
% ########### comment this part out if aligned for choice ###########
% temp=[trialStartStrobe(1);strobesFromVars];
% clear trialStartStrobe
% trialStartStrobe=temp;
%####################################################################
trialStartTimes=Strobed(trialStartStrobe,1);% find the time of each trial starts
lfp.uV=[];
for trialNum=1:length(trialStartTimes)
    clear uV startat stopat
    trialStartTime=trialStartTimes(trialNum);% start time of each trial
    startat=trialStartTime-0.5; % 500ms before actual start time
    stopat=trialStartTime+9.5; % 950ms after actual start time
    LFPtimes=intersect(find(LFP.timeStamp>startat),find(LFP.timeStamp<stopat));% SPKs are actually time stamps
    uV=LFP.uV(:,LFPtimes)';
    if size(uV,1)~=10000
        tmean=mean(uV,1);
        uV(end+1:10000,:)=repmat(tmean,10000-size(uV,1),1);
        clear tmean;
    end
    if ~isempty(lfp.uV)
        lfp.uV=[lfp.uV, uV];
    else
        lfp.uV=uV;
    end

%     keyboard;
%     lfp.trl(trialNum).uV=uV;
    
    disp(trialNum);
end

lfp.vars=vars;

% save([fname(1:9) '_OFCin_LFP.mat'], 'lfp', '-v7.3')
% save([fname(1:9) '_OFCout_LFP.mat'], 'lfp', '-v7.3')
% save([fname(1:9) '_PCCwrapped_LFP.mat'], 'lfp', '-v7.3')

save([spath fname(1:9) '_' area '_LFPat' align '.mat'], 'lfp', '-v7.3')
cd(spath)


