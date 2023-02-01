% decodeGC
% last used: Mar 13 2020
clear all; close all;clc
% fpath='/home/jeeves-raid2/benh-data/Maze/AR_Sorted/w_Anatomy/2020final/code/Decode/';
fpath='/Users/mwang/Google Drive/01Data/00 w_Anatomy/2020final/code/Decode/Decode_fromServer/';
addpath(genpath(fpath))
% dpath='/home/jeeves-raid2/benh-data/Maze/AR_Sorted/w_Anatomy/2020final/results/decode/';
% dpath='/Users/mwang/Google Drive/01Data/00 w_Anatomy/2020final/results/Decoding/fromServer/';
dpath='/Users/mwang/Google Drive/01Data/00 w_Anatomy/2020final/results/Decoding/results/';

% spath=dpath;
cd(dpath)
List=dir('*decodePS*');

%%
for fn=1:length(List)
%     List(fn).name
%     keyboard
    data(fn)=load([dpath List(fn).name]);
end

%% 
clc
in_lrC=[data(2).A.choLRC.lda, data(1).A.choLRC.lda]';
out_lrC=[data(4).A.choLRC.lda, data(3).A.choLRC.lda]';
pcc_lrC=[data(6).A.choLRC.lda, data(5).A.choLRC.lda]';

in_lrE=[data(2).A.choLRE.lda, data(1).A.choLRE.lda]';
out_lrE=[data(4).A.choLRE.lda, data(3).A.choLRE.lda]';
pcc_lrE=[data(6).A.choLRE.lda, data(5).A.choLRE.lda]';

in_12C=[data(2).A.choOptC.lda, data(1).A.choOptC.lda]';
out_12C=[data(4).A.choOptC.lda, data(3).A.choOptC.lda]';
pcc_12C=[data(6).A.choOptC.lda, data(5).A.choOptC.lda]';

in_12E=[data(2).A.choOptE.lda, data(1).A.choOptE.lda]';
out_12E=[data(4).A.choOptE.lda, data(3).A.choOptE.lda]';
pcc_12E=[data(6).A.choOptE.lda, data(5).A.choOptE.lda]';

in_ev1C=[data(2).A.ev1highC.lda, data(1).A.ev1highC.lda]';
out_ev1C=[data(4).A.ev1highC.lda, data(3).A.ev1highC.lda]';
pcc_ev1C=[data(6).A.ev1highC.lda, data(5).A.ev1highC.lda]';

in_ev1E=[data(2).A.ev1highE.lda, data(1).A.ev1highE.lda]';
out_ev1E=[data(4).A.ev1highE.lda, data(3).A.ev1highE.lda]';
pcc_ev1E=[data(6).A.ev1highE.lda, data(5).A.ev1highE.lda]';

in_ev2C=[data(2).A.ev2highC.lda, data(1).A.ev2highC.lda]';
out_ev2C=[data(4).A.ev2highC.lda, data(3).A.ev2highC.lda]';
pcc_ev2C=[data(6).A.ev2highC.lda, data(5).A.ev2highC.lda]';

in_ev2E=[data(2).A.ev2highE.lda, data(1).A.ev2highE.lda]';
out_ev2E=[data(4).A.ev2highE.lda, data(3).A.ev2highE.lda]';
pcc_ev2E=[data(6).A.ev2highE.lda, data(5).A.ev2highE.lda]';

%%
clc

for nlag=1:13 %13

    %  ################## CORRECT TRIALs ##################
    % hypothesis test direction 1
    % sig!!! at lag 4 
%       [h,pvalue,stat,cvalue]=gctest(in_12C,pcc_lrC,...
%           [in_lrC,pcc_12C,out_12C,out_lrC],'NumLags',nlag,'Trend',1);

     % hypothesis test direction 2
%       [h,pvalue,stat,cvalue]=gctest(pcc_lrC,in_12C,...
%           [in_lrC,pcc_12C,out_12C,out_lrC],'NumLags',nlag,'Trend',1);

    % control test XX--> PCC lr
      [h,pvalue,stat,cvalue]=gctest(out_12C,pcc_lrC,...
          [in_lrC,pcc_12C,in_12C,out_lrC],'NumLags',nlag,'Trend',1);

%       [h,pvalue,stat,cvalue]=gctest(pcc_12C,pcc_lrC,...
%           [in_12C,in_lrC,out_12C,out_lrC],'NumLags',nlag,'Trend',1);

%       [h,pvalue,stat,cvalue]=gctest(in_lrC,pcc_lrC,...
%           [in_12C,pcc_12C,out_12C,out_lrC],'NumLags',nlag,'Trend',1);
   
%       [h,pvalue,stat,cvalue]=gctest(out_lrC,pcc_lrC,...
%           [in_lrC,pcc_12C,in_12C,out_12C],'NumLags',nlag,'Trend',1);


      % ################## EVs --> pcc lr C ##################
      % not sig 
%       [h,pvalue,stat,cvalue]=gctest(in_ev1C, pcc_lrC, ...
%           [in_lrC,pcc_12C,in_12C,out_12C,out_lrC],'NumLags',nlag,'Trend',1);
      
      % not sig 
%       [h,pvalue,stat,cvalue]=gctest(in_ev2C, pcc_lrC, ...
%           [in_lrC,pcc_12C,in_12C,out_12C,out_lrC],'NumLags',nlag,'Trend',1);
%         
%       % not sig 
%       [h,pvalue,stat,cvalue]=gctest(pcc_ev1C, pcc_lrC, ...
%           [in_lrC,pcc_12C,in_12C,out_12C,out_lrC],'NumLags',nlag,'Trend',1);
%       
%       % sig  3 4 8
%       [h,pvalue,stat,cvalue]=gctest(pcc_ev2C, pcc_lrC, ...
%           [in_lrC,pcc_12C,in_12C,out_12C,out_lrC],'NumLags',nlag,'Trend',1);
%       
      
      % ################## EVs --> EVs ##################
%       % not sig 
%       [h,pvalue,stat,cvalue]=gctest(in_ev1C, in_ev2C, ...
%           'NumLags',nlag,'Trend',1);
      
        % sig 9
%       [h,pvalue,stat,cvalue]=gctest(in_ev1C, in_12C, ...
%           [in_ev2C],'NumLags',nlag,'Trend',1);

        % no sig
%       [h,pvalue,stat,cvalue]=gctest(in_ev2C, in_12C, ...
%           [in_ev1C],'NumLags',nlag,'Trend',1);
      
      
      % no sig
%       [h,pvalue,stat,cvalue]=gctest(in_ev1C, pcc_12C, ...
%           [in_ev2C],'NumLags',nlag,'Trend',1);
      
      % no sig
%       [h,pvalue,stat,cvalue]=gctest(in_ev2C, pcc_12C, ...
%           [in_ev1C],'NumLags',nlag,'Trend',1);
%       
%       % sig 9
%       [h,pvalue,stat,cvalue]=gctest(in_ev1C, pcc_ev1C, ...
%           [in_ev2C,pcc_ev2C],'NumLags',nlag,'Trend',1);

        % no sig 
%       [h,pvalue,stat,cvalue]=gctest(in_ev1E, pcc_ev1E, ...
%           [in_ev2E,pcc_ev2E],'NumLags',nlag,'Trend',1);

%       % no sig
%       [h,pvalue,stat,cvalue]=gctest(in_ev2C,pcc_ev2C, ...
%           [in_ev1C, pcc_ev1C,],'NumLags',nlag,'Trend',1);
      
%       % no sig
%       [h,pvalue,stat,cvalue]=gctest(out_ev1C, pcc_ev1C, ...
%           [out_ev2C, pcc_ev2C,],'NumLags',nlag,'Trend',1);
      
      
      % ################## in 12 --> pcc lr E ##################
    % not sig!!
%       [h,pvalue,stat,cvalue]=gctest(in_12E,pcc_lrE,...
%           [in_lrE,pcc_12E,out_12E,out_lrE],'NumLags',nlag,'Trend',1);
      
    % not sig!!! 
%       [h,pvalue,stat,cvalue]=gctest(in_lrE,pcc_lrE,...
%           [in_12E,pcc_12E,out_12E,out_lrE],'NumLags',nlag,'Trend',1);

    % sig 15
%       [h,pvalue,stat,cvalue]=gctest(in_12E,pcc_12E,...
%           [in_lrE,pcc_lrE,out_12E,out_lrE],'NumLags',nlag,'Trend',1);

%       % not sig 
%       [h,pvalue,stat,cvalue]=gctest(pcc_12E, pcc_lrE,...
%           [in_lrE,in_12E,out_12E,out_lrE],'NumLags',nlag,'Trend',1);

      % sig 15 
%       [h,pvalue,stat,cvalue]=gctest(pcc_lrE,pcc_12E, ...
%           [in_lrE,in_12E,out_12E,out_lrE],'NumLags',nlag,'Trend',1);

    % ################## what give rise to in_lr C ##################
    
    % not sig 
%       [h,pvalue,stat,cvalue]=gctest(in_12C, in_lrC,...
%           [pcc_12C,pcc_lrC,out_12C,out_lrC],'NumLags',nlag,'Trend',1);

    % sig 2 !!!! 
%       [h,pvalue,stat,cvalue]=gctest(pcc_12C, in_lrC,...
%           [in_12C,pcc_lrC,out_12C,out_lrC],'NumLags',nlag,'Trend',1);

    % NOT SIG for E !!!! 
%           [h,pvalue,stat,cvalue]=gctest(pcc_12E, in_lrE,...
%               [in_12E,pcc_lrE,out_12E,out_lrE],'NumLags',nlag,'Trend',1);
      
    % not sig 
%       [h,pvalue,stat,cvalue]=gctest(pcc_12C, in_lrC,...
%           [in_12C,pcc_lrC,out_12C,out_lrC],'NumLags',nlag,'Trend',1);

    %  not sig
%       [h,pvalue,stat,cvalue]=gctest(out_12C, in_lrC,...
%           [pcc_12C,in_12C,pcc_lrC,out_lrC],'NumLags',nlag,'Trend',1);  
      
    % not sig 
%       [h,pvalue,stat,cvalue]=gctest(out_lrC, in_lrC,...
%           [pcc_12C,in_12C,out_12C,pcc_lrC],'NumLags',nlag,'Trend',1); 
    
    % ################## what give rise to pcc_12 C ##################
    % sig lag 15:16
%       [h,pvalue,stat,cvalue]=gctest(in_lrC,pcc_12C,...
%           [in_12C,pcc_lrC,out_12C,out_lrC],'NumLags',nlag,'Trend',1);
          
      % sig lag 8!!!
%       [h,pvalue,stat,cvalue]=gctest(pcc_lrC,in_12C,...
%           [in_lrC,pcc_12C,out_12C,out_lrC],'NumLags',nlag,'Trend',1);

      % sig 2
%       [h,pvalue,stat,cvalue]=gctest(pcc_12C, in_lrC,...
%           [pcc_lrC,in_12C,out_12C,out_lrC],'NumLags',nlag,'Trend',1);
      
      % sig 
%       [h,pvalue,stat,cvalue]=gctest(pcc_lrC, pcc_lrC,...
%           [in_lrC,in_12C,out_12C,out_lrC],'NumLags',nlag,'Trend',1);
%       
%       % sig 
%       [h,pvalue,stat,cvalue]=gctest(pcc_12C, pcc_lrC,...
%           [in_lrC,in_12C,out_12C,out_lrC],'NumLags',nlag,'Trend',1);
%       
%       % sig 
%       [h,pvalue,stat,cvalue]=gctest(pcc_12C, pcc_lrC,...
%           [in_lrC,in_12C,out_12C,out_lrC],'NumLags',nlag,'Trend',1);
%       
%       % sig 
%       [h,pvalue,stat,cvalue]=gctest(pcc_12C, pcc_lrC,...
%           [in_lrC,in_12C,out_12C,out_lrC],'NumLags',nlag,'Trend',1);
%       
%       % sig 
%       [h,pvalue,stat,cvalue]=gctest(pcc_12C, pcc_lrC,...
%           [in_lrC,in_12C,out_12C,out_lrC],'NumLags',nlag,'Trend',1);
%       
      

      if pvalue<0.05
        disp('sig')
        nlag
      end
end


disp('##################### DONE #####################')


%%
clc
[Beta, sigma, E, covb,logL]=mvregress([in_12C,in_lrC],[pcc_12C,pcc_lrC])

%%
clc
% md=fitlm([in_12C,in_lrC],pcc_lrC)
md=fitlm([in_12C(1:6),in_lrC(1:6)],pcc_lrC(7:12))

%%
clc
md=fitlm(in_12C(1:6),pcc_lrC(7:12))





