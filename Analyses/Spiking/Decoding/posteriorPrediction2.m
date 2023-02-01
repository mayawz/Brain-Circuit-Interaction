% posteriorPrediction
% last used: Mar 12 2020
% use the decoding posterior from each region to predict that of another
% region

% Data info
% 4 fold cross validation -- 4 columns
% 6 time window
% at OFFER: 251:50:550
% 251   301   351   401   451   501
% at CHOICE:  201:50:500
% 201   251   301   351   401   451

clear all; close all;clc
% file path
fpath='/Users/mwang/Google Drive/01Data/00 w_Anatomy/2020final/code/Decode/';
addpath(genpath(fpath))
dpath='/Users/mwang/Google Drive/01Data/00 w_Anatomy/2020final/results/Decoding/results/';
spath='/Users/mwang/Google Drive/01Data/00 w_Anatomy/2020final/results/Decoding/results/';
cd(dpath)
ftnames={'OFCin_atOFF', 'OFCout_atOFF','PCC_atOFF',...
    'OFCin_atCHO', 'OFCout_atCHO','PCC_atCHO'};
for ftn=1:length(ftnames)
    List(ftn,1:2)=dir(['*' ftnames{ftn} '*']);
end
%

%%

for ln=1:size(List,1)
    
    for fn=1:size(List,2)
        
        data(fn)=load([dpath List(ln,fn).name]);
        
    end
%%    
    for tp=1:size(data(fn).accuracy,2)
        tmp=[data(1).accuracy(tp).choOptC.lda;...
            data(2).accuracy(tp).choOptC.lda];
        A.choOptC.lda(tp)=nanmean(tmp);
        Asem.choOptC.lda(tp)=confidenceInterval(tmp,1);
        A.choOptC.Size(tp)=length(tmp);
        clear tmp
        tmp=[data(1).accuracy(tp).choOptC.nbl;...
            data(2).accuracy(tp).choOptC.nbl];
        A.choOptC.nbl(tp)=nanmean(tmp);
        Asem.choOptC.nbl(tp)=confidenceInterval(tmp,1);
        clear tmp
        
        tmp=[data(1).accuracy(tp).choOptE.lda;...
            data(2).accuracy(tp).choOptE.lda];
        A.choOptE.lda(tp)=nanmean(tmp);
        Asem.choOptE.lda(tp)=confidenceInterval(tmp,1);
        A.choOptE.Size(tp)=length(tmp);
        clear tmp
        tmp=[data(1).accuracy(tp).choOptE.nbl;...
            data(2).accuracy(tp).choOptE.nbl];
        A.choOptE.nbl(tp)=nanmean(tmp);
        Asem.choOptE.nbl(tp)=confidenceInterval(tmp,1);
        clear tmp
        
        % ###########################################
        tmp=[data(1).accuracy(tp).choLRC.lda;...
            data(2).accuracy(tp).choLRC.lda];
        A.choLRC.lda(tp)=nanmean(tmp);
        Asem.choLRC.lda(tp)=confidenceInterval(tmp,1);
        A.choLRC.Size(tp)=length(tmp);
        clear tmp
        tmp=[data(1).accuracy(tp).choLRC.nbl;...
            data(2).accuracy(tp).choLRC.nbl];
        A.choLRC.nbl(tp)=nanmean(tmp);
        Asem.choLRC.nbl(tp)=confidenceInterval(tmp,1);
        clear tmp
        
        tmp=[data(1).accuracy(tp).choLRE.lda;...
            data(2).accuracy(tp).choLRE.lda];
        A.choLRE.lda(tp)=nanmean(tmp);
        Asem.choLRE.lda(tp)=confidenceInterval(tmp,1);
        A.choLRE.Size(tp)=length(tmp);
        clear tmp
        tmp=[data(1).accuracy(tp).choLRE.nbl;...
            data(2).accuracy(tp).choLRE.nbl];
        A.choLRE.nbl(tp)=nanmean(tmp);
        Asem.choLRE.nbl(tp)=confidenceInterval(tmp,1);
        clear tmp
        
        % ###########################################
        tmp=[data(1).accuracy(tp).ev1highC.lda;...
            data(2).accuracy(tp).ev1highC.lda];
        A.ev1highC.lda(tp)=nanmean(tmp);
        Asem.ev1highC.lda(tp)=confidenceInterval(tmp,1);
        A.ev1highC.Size(tp)=length(tmp);
        clear tmp
        tmp=[data(1).accuracy(tp).ev1highC.nbl;...
            data(2).accuracy(tp).ev1highC.nbl];
        A.ev1highC.nbl(tp)=nanmean(tmp);
        Asem.ev1highC.nbl(tp)=confidenceInterval(tmp,1);
        clear tmp
        
        tmp=[data(1).accuracy(tp).ev1highE.lda;...
            data(2).accuracy(tp).ev1highE.lda];
        A.ev1highE.lda(tp)=nanmean(tmp);
        Asem.ev1highE.lda(tp)=confidenceInterval(tmp,1);
        A.ev1highE.Size(tp)=length(tmp);
        clear tmp
        tmp=[data(1).accuracy(tp).ev1highE.nbl;...
            data(2).accuracy(tp).ev1highE.nbl];
        A.ev1highE.nbl(tp)=nanmean(tmp);
        Asem.ev1highE.nbl(tp)=confidenceInterval(tmp,1);
        clear tmp
        
        % ###########################################
        tmp=[data(1).accuracy(tp).ev2highC.lda;...
            data(2).accuracy(tp).ev2highC.lda];
        A.ev2highC.lda(tp)=nanmean(tmp);
        Asem.ev2highC.lda(tp)=confidenceInterval(tmp,1);
        A.ev2highC.Size(tp)=length(tmp);
        clear tmp
        tmp=[data(1).accuracy(tp).ev2highC.nbl;...
            data(2).accuracy(tp).ev2highC.nbl];
        A.ev2highC.nbl(tp)=nanmean(tmp);
        Asem.ev2highC.nbl(tp)=confidenceInterval(tmp,1);
        clear tmp
        
        tmp=[data(1).accuracy(tp).ev2highE.lda;...
            data(2).accuracy(tp).ev2highE.lda];
        A.ev2highE.lda(tp)=nanmean(tmp);
        Asem.ev2highE.lda(tp)=confidenceInterval(tmp,1);
        A.ev2highE.Size(tp)=length(tmp);
        clear tmp
        tmp=[data(1).accuracy(tp).ev2highE.nbl;...
            data(2).accuracy(tp).ev2highE.nbl];
        A.ev2highE.nbl(tp)=nanmean(tmp);
        Asem.ev2highE.nbl(tp)=confidenceInterval(tmp,1);
        clear tmp
        
        % ###########################################
        tmp=[data(1).accuracy(tp).EV1highCombo.lda;...
            data(2).accuracy(tp).EV1highCombo.lda];
        A.EV1highCombo.lda(tp)=nanmean(tmp);
        Asem.EV1highCombo.lda(tp)=confidenceInterval(tmp,1);
        A.EV1highCombo.Size(tp)=length(tmp);
        clear tmp
        tmp=[data(1).accuracy(tp).EV1highCombo.nbl;...
            data(2).accuracy(tp).EV1highCombo.nbl];
        A.EV1highCombo.nbl(tp)=nanmean(tmp);
        Asem.EV1highCombo.nbl(tp)=confidenceInterval(tmp,1);
        clear tmp
        
        % ###########################################
        tmp=[data(1).accuracy(tp).EV2highCombo.lda;...
            data(2).accuracy(tp).EV2highCombo.lda];
        A.EV2highCombo.lda(tp)=nanmean(tmp);
        Asem.EV2highCombo.lda(tp)=confidenceInterval(tmp,1);
        A.EV2highCombo.Size(tp)=length(tmp);
        clear tmp
        tmp=[data(1).accuracy(tp).EV2highCombo.nbl;...
            data(2).accuracy(tp).EV2highCombo.nbl];
        A.EV2highCombo.nbl(tp)=nanmean(tmp);
        Asem.EV2highCombo.nbl(tp)=confidenceInterval(tmp,1);
        clear tmp
        
        tmp=[data(1).post(tp).choOptC.lda;...
            data(2).post(tp).choOptC.lda];
        P.choOptC.lda(tp)=nanmean(tmp);
        Psem.choOptC.lda(tp)=confidenceInterval(tmp,1);
        P.choOptC.Size(tp)=length(tmp);
        clear tmp
        tmp=[data(1).post(tp).choOptC.nbl;...
            data(2).post(tp).choOptC.nbl];
        P.choOptC.nbl(tp)=nanmean(tmp);
        Psem.choOptC.nbl(tp)=confidenceInterval(tmp,1);
        clear tmp
        
        tmp=[data(1).post(tp).choOptE.lda;...
            data(2).post(tp).choOptE.lda];
        P.choOptE.lda(tp)=nanmean(tmp);
        Psem.choOptE.lda(tp)=confidenceInterval(tmp,1);
        P.choOptE.Size(tp)=length(tmp);
        clear tmp
        tmp=[data(1).post(tp).choOptE.nbl;...
            data(2).post(tp).choOptE.nbl];
        P.choOptE.nbl(tp)=nanmean(tmp);
        Psem.choOptE.nbl(tp)=confidenceInterval(tmp,1);
        clear tmp
        
        % ###########################################
        tmp=[data(1).post(tp).choLRC.lda;...
            data(2).post(tp).choLRC.lda];
        P.choLRC.lda(tp)=nanmean(tmp);
        Psem.choLRC.lda(tp)=confidenceInterval(tmp,1);
        P.choLRC.Size(tp)=length(tmp);
        clear tmp
        tmp=[data(1).post(tp).choLRC.nbl;...
            data(2).post(tp).choLRC.nbl];
        P.choLRC.nbl(tp)=nanmean(tmp);
        Psem.choLRC.nbl(tp)=confidenceInterval(tmp,1);
        clear tmp
        
        tmp=[data(1).post(tp).choLRE.lda;...
            data(2).post(tp).choLRE.lda];
        P.choLRE.lda(tp)=nanmean(tmp);
        Psem.choLRE.lda(tp)=confidenceInterval(tmp,1);
        P.choLRE.Size(tp)=length(tmp);
        clear tmp
        tmp=[data(1).post(tp).choLRE.nbl;...
            data(2).post(tp).choLRE.nbl];
        P.choLRE.nbl(tp)=nanmean(tmp);
        Psem.choLRE.nbl(tp)=confidenceInterval(tmp,1);
        clear tmp
        
        % ###########################################
        tmp=[data(1).post(tp).ev1highC.lda;...
            data(2).post(tp).ev1highC.lda];
        P.ev1highC.lda(tp)=nanmean(tmp);
        Psem.ev1highC.lda(tp)=confidenceInterval(tmp,1);
        P.ev1highC.Size(tp)=length(tmp);
        clear tmp
        tmp=[data(1).post(tp).ev1highC.nbl;...
            data(2).post(tp).ev1highC.nbl];
        P.ev1highC.nbl(tp)=nanmean(tmp);
        Psem.ev1highC.nbl(tp)=confidenceInterval(tmp,1);
        clear tmp
        
        tmp=[data(1).post(tp).ev1highE.lda;...
            data(2).post(tp).ev1highE.lda];
        P.ev1highE.lda(tp)=nanmean(tmp);
        Psem.ev1highE.lda(tp)=confidenceInterval(tmp,1);
        P.ev1highE.Size(tp)=length(tmp);
        clear tmp
        tmp=[data(1).post(tp).ev1highE.nbl;...
            data(2).post(tp).ev1highE.nbl];
        P.ev1highE.nbl(tp)=nanmean(tmp);
        Psem.ev1highE.nbl(tp)=confidenceInterval(tmp,1);
        clear tmp
        
        % ###########################################
        tmp=[data(1).post(tp).ev2highC.lda;...
            data(2).post(tp).ev2highC.lda];
        P.ev2highC.lda(tp)=nanmean(tmp);
        Psem.ev2highC.lda(tp)=confidenceInterval(tmp,1);
        P.ev2highC.Size(tp)=length(tmp);
        clear tmp
        tmp=[data(1).post(tp).ev2highC.nbl;...
            data(2).post(tp).ev2highC.nbl];
        P.ev2highC.nbl(tp)=nanmean(tmp);
        Psem.ev2highC.nbl(tp)=confidenceInterval(tmp,1);
        clear tmp
        
        tmp=[data(1).post(tp).ev2highE.lda;...
            data(2).post(tp).ev2highE.lda];
        P.ev2highE.lda(tp)=nanmean(tmp);
        Psem.ev2highE.lda(tp)=confidenceInterval(tmp,1);
        P.ev2highE.Size(tp)=length(tmp);
        clear tmp
        tmp=[data(1).post(tp).ev2highE.nbl;...
            data(2).post(tp).ev2highE.nbl];
        P.ev2highE.nbl(tp)=nanmean(tmp);
        Psem.ev2highE.nbl(tp)=confidenceInterval(tmp,1);
        clear tmp
        
        % ###########################################
        tmp=[data(1).post(tp).EV1highCombo.lda;...
            data(2).post(tp).EV1highCombo.lda];
        P.EV1highCombo.lda(tp)=nanmean(tmp);
        Psem.EV1highCombo.lda(tp)=confidenceInterval(tmp,1);
        P.EV1highCombo.Size(tp)=length(tmp);
        clear tmp
        tmp=[data(1).post(tp).EV1highCombo.nbl;...
            data(2).post(tp).EV1highCombo.nbl];
        P.EV1highCombo.nbl(tp)=nanmean(tmp);
        Psem.EV1highCombo.nbl(tp)=confidenceInterval(tmp,1);
        clear tmp
        
        % ###########################################
        tmp=[data(1).post(tp).EV2highCombo.lda;...
            data(2).post(tp).EV2highCombo.lda];
        P.EV2highCombo.lda(tp)=nanmean(tmp);
        Psem.EV2highCombo.lda(tp)=confidenceInterval(tmp,1);
        P.EV2highCombo.Size(tp)=length(tmp);
        clear tmp
        tmp=[data(1).post(tp).EV2highCombo.nbl;...
            data(2).post(tp).EV2highCombo.nbl];
        P.EV2highCombo.nbl(tp)=nanmean(tmp);
        Psem.EV2highCombo.nbl(tp)=confidenceInterval(tmp,1);
        clear tmp
    end
    
    
    %%
    save([spath 'decodePS.' List(ln,fn).name(11:21) '.mat'])
    %%
    clear data A
    
    
end
