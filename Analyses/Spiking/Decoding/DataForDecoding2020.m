% DataForDecoding
% last used: Feb 10 2020

% Here organize data into
% columns: # of neuron dimentional
% rows: sliding time of 200ms with 50ms steps by deviding dimension of
% interest
% deviding by: (1) choOpt: offer1 offer2;
%              (2) choLR: choLeft choRight
%              (3) EV1>=mean(EV1)
%              (4) EV2>=mean(EV2)

clear all; close all; clc
% data path
dpath='/Users/mwang/Google Drive/01Data/00 w_Anatomy/2020final/wrapped/';
% file path
% fpath='/Users/mwang/Google Drive/01Data/00 w_Anatomy/Gambling/Nov2019Analyses/';
fpath='/Users/mwang/Google Drive/01Data/00 w_Anatomy/2020final/code/';
addpath(genpath(fpath))
% save path
spath='/Users/mwang/Google Drive/01Data/00 w_Anatomy/2020final/wrapped/decodePatterns/';

cutSafe=1;

cd(dpath)
%%
for segments=1:2
    if segments==1
        List=dir('*atOFFER1*');
        bins=251:550;
        stepSize=5;
        binStarts=251:stepSize:550-stepSize+1;
        
    else
        List=dir('*atCHOICE*');
        bins=201:500;
        stepSize=5;
        binStarts=201:stepSize:500-stepSize+1;
    end
    for fn=1:length(List)
        clear psth vars data lm rm lp rp firstAppear choice choLR choOpt ...
            outcome trls errorTrls correctTrls cho1Corr cho2Corr cho1Err ...
            cho2Err choLCorr choLErr choRCorr choRErr
        
%         if fn==4
%             keyboard
%         end
        
        fname=List(fn).name(1:end-4);
        load([dpath List(fn).name]);
        
        D=data;
        clear data
        
        for cn=1:length(D)
            disp(fname)
            disp(['Cell: ' num2str(cn)])
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
            choLR=choice; % 1:Left 2:Right
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
            
            lev=lp.*lm;
            rev=rp.*rm;
            lCorr=intersect(find(lev>rev),find(choice==1));
            rCorr=intersect(find(lev<rev),find(choice==2));
            correctTrls=sort([lCorr; rCorr],'ascend');
            
            lErr=intersect(find(lev<rev),find(choice==1));
            rErr=intersect(find(lev>rev),find(choice==2));
            errorTrls=sort([lErr; rErr],'ascend');
            
            cho1Corr=intersect(find(choOpt==1),correctTrls);
            cho2Corr=intersect(find(choOpt==2),correctTrls);
            cho1Err=intersect(find(choOpt==1),errorTrls);
            cho2Err=intersect(find(choOpt==2),errorTrls);
            
            choLCorr=intersect(find(choLR==1),correctTrls);
            choRCorr=intersect(find(choLR==-1),correctTrls);
            choLErr=intersect(find(choLR==1),errorTrls);
            choRErr=intersect(find(choLR==-1),errorTrls);
            
            
            
            for bn=1:length(binStarts)
                
                binSlide=binStarts(bn):binStarts(bn)+stepSize-1;
                
                % activity pattern for EV1 decode for high/low
                for tn=1:length(ev1)
                    patt(bn).ev1fr(tn,cn)=mean(psth(tn,binSlide));
                    if cn==1
                        patt(bn).EV1_high=ev1>=mean((ev1+ev2)/2);
                    end
                end
                
                % activity pattern for ev1 Correct
                for tn=1:length(correctTrls)
                    patt(bn).ev1fr_corr(tn,cn)=mean(psth(correctTrls(tn),binSlide));
                    if cn==1
                        tmp=ev1>=mean((ev1+ev2)/2);
%                         keyboard
                        patt(bn).EV1_high_corr=tmp(correctTrls);
                        clear tmp
                    end
                end
                
                % activity pattern for ev1 Error
                 for tn=1:length(errorTrls)
                    patt(bn).ev1fr_erro(tn,cn)=mean(psth(errorTrls(tn),binSlide));
                    if cn==1
                        tmp=ev1>=mean((ev1+ev2)/2);
                        patt(bn).EV1_high_erro=tmp(errorTrls);
                        clear tmp
                    end
                 end
              
                % activity pattern for EV2 decode for high/low
                
                for tn=1:length(ev1)
                    patt(bn).ev2fr(tn,cn)=mean(psth(tn,binSlide));
                    if cn==1
                        patt(bn).EV2_high=ev2>=mean((ev1+ev2)/2);
                    end
                end
                
                % activity pattern for ev2 Correct
                 for tn=1:length(correctTrls)
                    patt(bn).ev2fr_corr(tn,cn)=mean(psth(correctTrls(tn),binSlide));
                    if cn==1
                        tmp=ev2>=mean((ev1+ev2)/2);
                        patt(bn).EV2_high_corr=tmp(correctTrls);
                        clear tmp
                    end
                 end
                 
                % activity pattern for ev2 Error
                for tn=1:length(errorTrls)
                    patt(bn).ev2fr_erro(tn,cn)=mean(psth(errorTrls(tn),binSlide));
                    if cn==1
                        tmp=ev2>=mean((ev1+ev2)/2);
                        patt(bn).EV2_high_erro=tmp(errorTrls);
                        clear tmp
                    end
                end
               
                % activity pattern for choOpt Correct
                for tn=1:length(correctTrls)
                    patt(bn).choOptCfr(tn,cn)=mean(psth(correctTrls(tn),binSlide));
                    if cn==1
                        patt(bn).choOptC=choOpt(correctTrls);
                    end
                end
                
                % activity pattern for choOpt Error
                for tn=1:length(errorTrls)
                    patt(bn).choOptEfr(tn,cn)=mean(psth(errorTrls(tn),binSlide));
                    if cn==1
                        patt(bn).choOptE=choOpt(errorTrls);
                    end
                end
                
                % activity pattern for choLR Correct
                for tn=1:length(correctTrls)
                    patt(bn).choLRCfr(tn,cn)=mean(psth(correctTrls(tn),binSlide));
                    if cn==1
                        patt(bn).choLRC=choLR(correctTrls);
                    end
                end
                
                % activity pattern for choLR Error
                for tn=1:length(errorTrls)
                    patt(bn).choLREfr(tn,cn)=mean(psth(errorTrls(tn),binSlide));
                    if cn==1
                        patt(bn).choLRE=choLR(errorTrls);
                    end
                end
                
            end % bn= time point number / moving window
            
        end %  cell number
        
        %%
        save([spath fname '_patterns.mat'],'patt')
        clear D cho patt
    end
end % end of segments=atOFFERS or atCHOICE

%%
cd(spath)












