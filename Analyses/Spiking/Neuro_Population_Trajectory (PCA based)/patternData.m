% patternData
% last used Feb 10 2020

% We first defined a 125-dimensional neuronal space, with each neuron
% taking up one dimension. Then we computed the activation state for each
% of five 300 ms epochs (offer 1 cue, offer 1 value, offer 2, choice and
% outcome) in each trial type, by averaging firing rates for each neuron
% across all trials and across time bins in each epoch. We subsequently
% conducted a principal component analysis on the 125-dimensional, 5-epoch,
% 2-trial-type, population responses.

% offer1 H - offer2 L - cho1 - L/R correct
% offer1 L - offer2 G - cho2 - L/R correct
% offer1 H - offer2 L - cho2 - L/R error
% offer1 L - offer2 H - cho1 - L/R error
% there are 4 correct condition and 4 error conditon. 

% PCA trajectories will be on correct trials only
% then project error trial activities onto the top PC space

% two kinds of PCA will be conducted
% 1. trial averaged
% columns: # of neurons
% rows: conditon X time point
% cho1LC all time points, meanFR(all trls in condition)
% cho1RC all time points
% etc etc

% 2. trial by trial
% columns: # of neurons
% rows: conditon X time point
% cho1LC all time points, FR(just 1 trial)
% cho1RC all time points
% etc etc

clear all; close all; clc
% data path
dpath='/Users/mwang/Google Drive/01Data/00 w_Anatomy/2020final/wrapped/';
% file path
% fpath='/Users/mwang/Google Drive/01Data/00 w_Anatomy/Gambling/Nov2019Analyses/';
fpath='/Users/mwang/Google Drive/01Data/00 w_Anatomy/2020final/code/';
% save path
spath='/Users/mwang/Google Drive/01Data/00 w_Anatomy/2020final/wrapped/PCAt_states/';

cutSafe=1;

cd(dpath)

for segments=1:2
    if segments==1
        List=dir('*atOFFER1*');
        bins=251:550;
        stepSize=50;
        binStarts=251:stepSize:550-stepSize+1;
        
    else
        List=dir('*atCHOICE*');
        bins=201:500;
        stepSize=50;
        binStarts=201:stepSize:500-stepSize+1;
    end
    for fn=1:length(List)
        clear psth vars data lm rm lp rp firstAppear choice choLR choOpt ...
            outcome trls condTrls
        
        fname=List(fn).name(1:end-4)
        load([dpath List(fn).name]);
        
        
        
        D=data;
        clear data
        for cn=1:length(D)
            
            cn
            if cutSafe
                temp=find(D(cn).vars(:,9)~=0);% 0:Safe 1:Lose 2:Win
                psth=D(cn).psth(temp,:);
                vars=D(cn).vars(temp,:);
            else
                psth=D(cn).psth;
                vars=D(cn).vars;
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
            choLR=vars(:,8); % 1:Left 2:Right
            choLR(choLR==2)=-1; % 1=left; -1=right
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
            
            corr1=intersect(find(ev1>=ev2),find(choOpt==1));
            corr2=intersect(find(ev2>=ev1),find(choOpt==2));
            correctTrls=sort([corr1; corr2],'ascend');
            
            erro1=intersect(find(ev1<ev2),find(choOpt==1));
            erro2=intersect(find(ev2<ev1),find(choOpt==2));
            errorTrls=sort([erro1; erro2],'ascend');
            
            % length(correctTrls)+length(errorTrls)
            
            cho1Corr=intersect(find(choOpt==1),correctTrls);
            cho2Corr=intersect(find(choOpt==2),correctTrls);
            cho1Err=intersect(find(choOpt==1),errorTrls);
            cho2Err=intersect(find(choOpt==2),errorTrls);
            
            cTrls{1}={intersect(cho1Corr,find(choLR==1))};
            cTrls{2}={intersect(cho1Corr,find(choLR==-1))};
            cTrls{3}={intersect(cho2Corr,find(choLR==1))};
            cTrls{4}={intersect(cho2Corr,find(choLR==-1))};
            
            eTrls{1}={intersect(cho1Err,find(choLR==1))};
            eTrls{2}={intersect(cho1Err,find(choLR==-1))};
            eTrls{3}={intersect(cho2Err,find(choLR==1))};
            eTrls{4}={intersect(cho2Err,find(choLR==-1))};
            for condn=1:length(cTrls)
                for tp=1:length(binStarts)
                    binSlide=binStarts(tp):binStarts(tp)+stepSize-1;
                    for tn=1:length(cell2mat(cTrls{condn}))
                        clear trlN
                        trlN=cell2mat(cTrls{condn});
                        cTrlN(condn)=length(trlN);
                        rowN=tp+length(binStarts)*(condn-1);
                        
                        correct.M{rowN,cn,tn}=mean(psth(trlN(tn),binSlide));
                        
                        correct.avg(rowN,cn)=mean(mean(psth(trlN,binSlide)));
                        % keyboard;
                        % unique(mean(psth(trlN,binSlide)))
                        % unique(median(psth(trlN,binSlide)))
                    end % tn= correct trial number
                    
                    for tn=1:length(cell2mat(eTrls{condn}))
                        
                        clear trlN
                        trlN=cell2mat(eTrls{condn});
                        eTrlN(condn)=length(trlN);
                        rowN=tp+length(binStarts)*(condn-1);
                        
                        error.M{rowN,cn,tn}=mean(psth(trlN(tn),binSlide));
                        error.avg(rowN,cn)=mean(mean(psth(trlN,binSlide)));
                    end % tn=error trial number
                end % tp=time point
            end % condn=condition number
            
        end %  cn=cell number
        
        
        %%
        % keyboard
        % cell2mat(correct.M(1,:,1))
        correct.correctOrder={'cho1LC','cho1RC','cho2LC','cho2RC'};
        correct.cTrlN=cTrlN;
        error.errorOrder={'cho1LE','cho1RE','cho2LE','cho2RE'};
        error.eTrlN=eTrlN;
        save([spath List(fn).name(1:end-4) '_states.mat'],'correct','error')
        clear D cho
        
        
        
    end % fn=file count in the list
    
    
    
end % end of segments=atOFFERS or atCHOICE

cd(spath)

% %% double check the data
% clear all; close all;clc
% 
% spath='/Users/mwang/Google Drive/01Data/00 w_Anatomy/2020final/wrapped/patternRaw/';
% 
% cd(spath)
% 
% load('P20171207_OFCin_atOFFER1_states.mat')
% correct.cTrlN % = 88   123   111   141
% correct.M(:,:,89)
% cell2mat(correct.M(:,:,89))
% 
% correct.M(:,:,123)
% 
% mean(cell2mat(correct.M(:,:,88)),3)
% 
% unique(mean(cell2mat(correct.M(:,:,88)),3))











