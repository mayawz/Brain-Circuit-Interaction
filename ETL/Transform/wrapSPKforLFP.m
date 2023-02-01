% wrap SPK for coherence with LFP
% Last used on:Feb 8 2020
% aligned at offer1 onset and then at choice



close all;
clear all;
clc;

fpath='/Users/mwang/Google Drive/01Data/00 w_Anatomy/Gambling/Sorted/OFCin/';
% fpath='/Users/mwang/Google Drive/01Data/00 w_Anatomy/Gambling/Sorted/OFCout/';
% fpath='/Users/mwang/Google Drive/01Data/00 w_Anatomy/Gambling/Sorted/PCC/';
cd(fpath)

spath='/Users/mwang/Google Drive/01Data/00 w_Anatomy/2020final/wrapped/forCoherence/';
% cd(spath)

List=dir('*P2017*.mat');
% List=dir('*S2018*.mat');

area='OFCin';
% area='OFCout';
% area='PCC';

% startStrobe=4001
% align='OFFER1';


startStrobe=4007
align='CHOICE';

% %% get vars and trl num
% tmp=load(List(1).name);
% fname=List(1).name(1:9);
% clear Strobed strobeName
% 
% Strobed=tmp.Strobed;
% strobeName(:,1)=Strobed(:,2);
% 
% temp=find(strobeName==20000);
% theCpltEnd=temp(end); % the end of all complete trials
% clear temp;
% 
% temp=find(strobeName==4001);
% theCpltBeginning=temp(1); % the beginning of all complete trials
% clear temp;
% 
% Strobed=Strobed(theCpltBeginning:theCpltEnd,:);
% 
% [vars strobesFromVars]= extractVars3(Strobed);
% clear tmp
% 
% spk.vars=vars;
%%

spk.psth=[];
for fn=1:length(List)
    clear vars Strobed strobeName d
    %% get vars and trl num
    tmp=load(List(fn).name);
    fname=List(fn).name(1:9);
    clear Strobed strobeName
    
    Strobed=tmp.Strobed;
    strobeName(:,1)=Strobed(:,2);
    
    temp=find(strobeName==20000);
    theCpltEnd=temp(end); % the end of all complete trials
    clear temp;
    
    temp=find(strobeName==4001);
    theCpltBeginning=temp(1); % the beginning of all complete trials
    clear temp;
    
    Strobed=Strobed(theCpltBeginning:theCpltEnd,:);
    
    [vars strobesFromVars]= extractVars3(Strobed);
    
    spk.vars=vars;
    %%
    
    D(fn).data=load(List(fn).name);
    d=struct2cell(D(fn).data);
    
    for i=1:length(d)-1 % loop cell num
        %%
        clear SPK psth
        
        disp(['Start Cell # ' num2str(i)])
        
        SPK=d{i+1};
        cellNum=i;
        
        % this gets 500ms before and 950ms after the strobe
        [psth,meanTrialDuration] = extractPSTH_forLFP(Strobed,SPK,startStrobe,0.001,strobesFromVars);
        
%         keyboard
        %%
        if isempty(spk.psth)
            spk.psth=psth';
        else
            spk.psth=[spk.psth, psth'];
        end
        
        %% by trial
%         for tn=1:size(vars,1)
%             if fn==1
%                 spk.trl(tn).smpBycell(:,i)=psth(tn,:)';
%             else
%                 if i==1
%                     tmp=size(spk.trl(tn).smpBycell,2);
%                 end
%                 spk.trl(tn).smpBycell(:,tmp+i)=psth(tn,:)';
%                 
%             end
%         end % loop trl num
    end
end

%%
% save(['/Users/mwang/Google Drive/01Data/00 w_Anatomy/Gambling/Sorted/' fname '_' area '_SPKatOFFER1.mat'],'spk', '-v7.3')

%%
save([spath fname '_' area '_SPKat' align '.mat'],'spk', '-v7.3')
cd(spath)

%%
size(spk.psth,2)/size(vars,1)
   


%%
%     names = {'TrialNum';'CumulativeTrialNum';'LeftProb';'RightProb';...
%         'LorRAppearedFirst';'LeftRwdMag';'RightRwdMag';...
%         'Choice'; 'GambleOutcome'; 'Aligned at 1st opt on'; fileName; num2str(simultaneousDay); Area};
%  

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
