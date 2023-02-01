% posteriorPrediction
% last used: Mar 12 2020
% use the decoding posterior from each region to predict that of another
% region
clear all; close all;clc


% file path
fpath='/Users/mwang/Google Drive/01Data/00 w_Anatomy/2020final/code/Decode/';
addpath(genpath(fpath))
dpath='/Users/mwang/Google Drive/01Data/00 w_Anatomy/2020final/results/Decoding/';
spath='/Users/mwang/Google Drive/01Data/00 w_Anatomy/2020final/results/Decoding/results/';
cd(dpath)
List=dir('*mat*');
%% Data info
% 4 fold cross validation -- 4 columns
% 6 time window
% at OFFER: 251:50:550
% 251   301   351   401   451   501
% at CHOICE:  201:50:500
% 201   251   301   351   401   451

%%
for fn=1:length(List)
    load([dpath List(fn).name])
    % time point by cross validation
    for tp=1:60
        for cv=1:4
            if cv==1
                
                % accuracy(tp)
                
                tmp=lda.choOptC(tp,cv).predict(1).predictedLabel==lda.choOptC(tp,cv).predict(1).trueLabel;
                accuracy(tp).choOptC.lda=tmp;
                clear tmp
                tmp=nbl.choOptC(tp,cv).predict(1).predictedLabel==nbl.choOptC(tp,cv).predict(1).trueLabel;
                accuracy(tp).choOptC.nbl=tmp;
                clear tmp
                tmp=lda.choOptC(tp,cv).predict(2).predictedLabel==lda.choOptC(tp,cv).predict(2).trueLabel;
                accuracy(tp).choOptE.lda=tmp;
                clear tmp
                tmp=nbl.choOptC(tp,cv).predict(2).predictedLabel==nbl.choOptC(tp,cv).predict(2).trueLabel;
                accuracy(tp).choOptE.nbl=tmp;
                clear tmp
                
                tmp=lda.choLRC(tp,cv).predict(1).predictedLabel==lda.choLRC(tp,cv).predict(1).trueLabel;
                accuracy(tp).choLRC.lda=tmp;
                clear tmp
                tmp=nbl.choLRC(tp,cv).predict(1).predictedLabel==nbl.choLRC(tp,cv).predict(1).trueLabel;
                accuracy(tp).choLRC.nbl=tmp;
                clear tmp
                tmp=lda.choLRC(tp,cv).predict(2).predictedLabel==lda.choLRC(tp,cv).predict(2).trueLabel;
                accuracy(tp).choLRE.lda=tmp;
                clear tmp
                tmp=nbl.choLRC(tp,cv).predict(2).predictedLabel==nbl.choLRC(tp,cv).predict(2).trueLabel;
                accuracy(tp).choLRE.nbl=tmp;
                clear tmp
                
                tmp=lda.EV1_high(tp,cv).predict.predictedLabel==lda.EV1_high(tp,cv).predict.trueLabel;
                accuracy(tp).EV1highCombo.lda=tmp;
                clear tmp
                tmp=nbl.EV1_high(tp,cv).predict.predictedLabel==nbl.EV1_high(tp,cv).predict.trueLabel;
                accuracy(tp).EV1highCombo.nbl=tmp;
                clear tmp
                tmp=lda.EV2_high(tp,cv).predict.predictedLabel==lda.EV2_high(tp,cv).predict.trueLabel;
                accuracy(tp).EV2highCombo.lda=tmp;
                clear tmp
                tmp=nbl.EV2_high(tp,cv).predict.predictedLabel==nbl.EV2_high(tp,cv).predict.trueLabel;
                accuracy(tp).EV2highCombo.nbl=tmp;
                clear tmp
                
                tmp=lda.ev1fr_corr(tp,cv).predict(1).predictedLabel==lda.ev1fr_corr(tp,cv).predict(1).trueLabel;
                accuracy(tp).ev1highC.lda=tmp;
                clear tmp
                tmp=nbl.ev1fr_corr(tp,cv).predict(1).predictedLabel==nbl.ev1fr_corr(tp,cv).predict(1).trueLabel;
                accuracy(tp).ev1highC.nbl=tmp;
                clear tmp
                tmp=lda.ev1fr_corr(tp,cv).predict(2).predictedLabel==lda.ev1fr_corr(tp,cv).predict(2).trueLabel;
                accuracy(tp).ev1highE.lda=tmp;
                clear tmp
                tmp=nbl.ev1fr_corr(tp,cv).predict(2).predictedLabel==nbl.ev1fr_corr(tp,cv).predict(2).trueLabel;
                accuracy(tp).ev1highE.nbl=tmp;
                clear tmp
                
                tmp=lda.ev2fr_corr(tp,cv).predict(1).predictedLabel==lda.ev2fr_corr(tp,cv).predict(1).trueLabel;
                accuracy(tp).ev2highC.lda=tmp;
                clear tmp
                tmp=nbl.ev2fr_corr(tp,cv).predict(1).predictedLabel==nbl.ev2fr_corr(tp,cv).predict(1).trueLabel;
                accuracy(tp).ev2highC.nbl=tmp;
                clear tmp
                tmp=lda.ev2fr_corr(tp,cv).predict(2).predictedLabel==lda.ev2fr_corr(tp,cv).predict(2).trueLabel;
                accuracy(tp).ev2highE.lda=tmp;
                clear tmp
                tmp=nbl.ev2fr_corr(tp,cv).predict(2).predictedLabel==nbl.ev2fr_corr(tp,cv).predict(2).trueLabel;
                accuracy(tp).ev2highE.nbl=tmp;
                clear tmp
                
                % posterior over trueLabel
                use1=find(lda.choOptC(tp,cv).predict(1).trueLabel=='1');
                use2=find(lda.choOptC(tp,cv).predict(1).trueLabel=='2');
                post(tp).choOptC.lda(use1,1)=lda.choOptC(tp,cv).predict(1).posterior(use1,1);
                post(tp).choOptC.lda(use2,1)=lda.choOptC(tp,cv).predict(1).posterior(use2,2);
                clear use1 use2
                use1=find(lda.choOptC(tp,cv).predict(2).trueLabel=='1');
                use2=find(lda.choOptC(tp,cv).predict(2).trueLabel=='2');
                post(tp).choOptE.lda(use1,1)=lda.choOptC(tp,cv).predict(2).posterior(use1,1);
                post(tp).choOptE.lda(use2,1)=lda.choOptC(tp,cv).predict(2).posterior(use2,2);
                clear use1 use2
                
                use1=find(nbl.choOptC(tp,cv).predict(1).trueLabel=='1');
                use2=find(nbl.choOptC(tp,cv).predict(1).trueLabel=='2');
                post(tp).choOptC.nbl(use1,1)=nbl.choOptC(tp,cv).predict(1).posterior(use1,1);
                post(tp).choOptC.nbl(use2,1)=nbl.choOptC(tp,cv).predict(1).posterior(use2,2);
                clear use1 use2
                use1=find(nbl.choOptC(tp,cv).predict(2).trueLabel=='1');
                use2=find(nbl.choOptC(tp,cv).predict(2).trueLabel=='2');
                post(tp).choOptE.nbl(use1,1)=nbl.choOptC(tp,cv).predict(2).posterior(use1,1);
                post(tp).choOptE.nbl(use2,1)=nbl.choOptC(tp,cv).predict(2).posterior(use2,2);
                clear use1 use2
                
                use1=find(lda.choLRC(tp,cv).predict(1).trueLabel=='1');
                use2=find(lda.choLRC(tp,cv).predict(1).trueLabel=='2');
                post(tp).choLRC.lda(use1,1)=lda.choLRC(tp,cv).predict(1).posterior(use1,1);
                post(tp).choLRC.lda(use2,1)=lda.choLRC(tp,cv).predict(1).posterior(use2,2);
                clear use1 use2
                use1=find(lda.choLRC(tp,cv).predict(2).trueLabel=='1');
                use2=find(lda.choLRC(tp,cv).predict(2).trueLabel=='2');
                post(tp).choLRE.lda(use1,1)=lda.choLRC(tp,cv).predict(2).posterior(use1,1);
                post(tp).choLRE.lda(use2,1)=lda.choLRC(tp,cv).predict(2).posterior(use2,2);
                clear use1 use2
                
                use1=find(nbl.choLRC(tp,cv).predict(1).trueLabel=='1');
                use2=find(nbl.choLRC(tp,cv).predict(1).trueLabel=='2');
                post(tp).choLRC.nbl(use1,1)=nbl.choLRC(tp,cv).predict(1).posterior(use1,1);
                post(tp).choLRC.nbl(use2,1)=nbl.choLRC(tp,cv).predict(1).posterior(use2,2);
                clear use1 use2
                use1=find(nbl.choLRC(tp,cv).predict(2).trueLabel=='1');
                use2=find(nbl.choLRC(tp,cv).predict(2).trueLabel=='2');
                post(tp).choLRE.nbl(use1,1)=nbl.choLRC(tp,cv).predict(2).posterior(use1,1);
                post(tp).choLRE.nbl(use2,1)=nbl.choLRC(tp,cv).predict(2).posterior(use2,2);
                clear use1 use2
                
                use1=find(lda.EV1_high(tp,cv).predict.trueLabel=='false');
                use2=find(lda.EV1_high(tp,cv).predict.trueLabel=='true');
                post(tp).EV1highCombo.lda(use1,1)=lda.EV1_high(tp,cv).predict.posterior(use1,1);
                post(tp).EV1highCombo.lda(use2,1)=lda.EV1_high(tp,cv).predict.posterior(use2,2);
                clear use1 use2
                use1=find(nbl.EV1_high(tp,cv).predict.trueLabel=='false');
                use2=find(nbl.EV1_high(tp,cv).predict.trueLabel=='true');
                post(tp).EV1highCombo.nbl(use1,1)=nbl.EV1_high(tp,cv).predict.posterior(use1,1);
                post(tp).EV1highCombo.nbl(use2,1)=nbl.EV1_high(tp,cv).predict.posterior(use2,2);
                clear use1 use2
            
                use1=find(lda.EV2_high(tp,cv).predict.trueLabel=='false');
                use2=find(lda.EV2_high(tp,cv).predict.trueLabel=='true');
                post(tp).EV2highCombo.lda(use1,1)=lda.EV2_high(tp,cv).predict.posterior(use1,1);
                post(tp).EV2highCombo.lda(use2,1)=lda.EV2_high(tp,cv).predict.posterior(use2,2);
                clear use1 use2
                use1=find(nbl.EV2_high(tp,cv).predict.trueLabel=='false');
                use2=find(nbl.EV2_high(tp,cv).predict.trueLabel=='true');
                post(tp).EV2highCombo.nbl(use1,1)=nbl.EV2_high(tp,cv).predict.posterior(use1,1);
                post(tp).EV2highCombo.nbl(use2,1)=nbl.EV2_high(tp,cv).predict.posterior(use2,2);
                clear use1 use2
                
                use1=find(lda.ev1fr_corr(tp,cv).predict(1).trueLabel=='false');
                use2=find(lda.ev1fr_corr(tp,cv).predict(1).trueLabel=='true');
                post(tp).ev1highC.lda(use1,1)=lda.ev1fr_corr(tp,cv).predict(1).posterior(use1,1);
                post(tp).ev1highC.lda(use2,1)=lda.ev1fr_corr(tp,cv).predict(1).posterior(use2,2);
                clear use1 use2
                use1=find(lda.ev1fr_corr(tp,cv).predict(2).trueLabel=='false');
                use2=find(lda.ev1fr_corr(tp,cv).predict(2).trueLabel=='true');
                post(tp).ev1highE.lda(use1,1)=lda.ev1fr_corr(tp,cv).predict(2).posterior(use1,1);
                post(tp).ev1highE.lda(use2,1)=lda.ev1fr_corr(tp,cv).predict(2).posterior(use2,2);
                clear use1 use2
                
                use1=find(nbl.ev1fr_corr(tp,cv).predict(1).trueLabel=='false');
                use2=find(nbl.ev1fr_corr(tp,cv).predict(1).trueLabel=='true');
                post(tp).ev1highC.nbl(use1,1)=nbl.ev1fr_corr(tp,cv).predict(1).posterior(use1,1);
                post(tp).ev1highC.nbl(use2,1)=nbl.ev1fr_corr(tp,cv).predict(1).posterior(use2,2);
                clear use1 use2
                use1=find(nbl.ev1fr_corr(tp,cv).predict(2).trueLabel=='false');
                use2=find(nbl.ev1fr_corr(tp,cv).predict(2).trueLabel=='true');
                post(tp).ev1highE.nbl(use1,1)=nbl.ev1fr_corr(tp,cv).predict(2).posterior(use1,1);
                post(tp).ev1highE.nbl(use2,1)=nbl.ev1fr_corr(tp,cv).predict(2).posterior(use2,2);
                clear use1 use2
                
                use1=find(lda.ev2fr_corr(tp,cv).predict(1).trueLabel=='false');
                use2=find(lda.ev2fr_corr(tp,cv).predict(1).trueLabel=='true');
                post(tp).ev2highC.lda(use1,1)=lda.ev2fr_corr(tp,cv).predict(1).posterior(use1,1);
                post(tp).ev2highC.lda(use2,1)=lda.ev2fr_corr(tp,cv).predict(1).posterior(use2,2);
                clear use1 use2
                use1=find(lda.ev2fr_corr(tp,cv).predict(2).trueLabel=='false');
                use2=find(lda.ev2fr_corr(tp,cv).predict(2).trueLabel=='true');
                post(tp).ev2highE.lda(use1,1)=lda.ev2fr_corr(tp,cv).predict(2).posterior(use1,1);
                post(tp).ev2highE.lda(use2,1)=lda.ev2fr_corr(tp,cv).predict(2).posterior(use2,2);
                clear use1 use2
                
                use1=find(nbl.ev2fr_corr(tp,cv).predict(1).trueLabel=='false');
                use2=find(nbl.ev2fr_corr(tp,cv).predict(1).trueLabel=='true');
                post(tp).ev2highC.nbl(use1,1)=nbl.ev2fr_corr(tp,cv).predict(1).posterior(use1,1);
                post(tp).ev2highC.nbl(use2,1)=nbl.ev2fr_corr(tp,cv).predict(1).posterior(use2,2);
                clear use1 use2
                use1=find(nbl.ev2fr_corr(tp,cv).predict(2).trueLabel=='false');
                use2=find(nbl.ev2fr_corr(tp,cv).predict(2).trueLabel=='true');
                post(tp).ev2highE.nbl(use1,1)=nbl.ev2fr_corr(tp,cv).predict(2).posterior(use1,1);
                post(tp).ev2highE.nbl(use2,1)=nbl.ev2fr_corr(tp,cv).predict(2).posterior(use2,2);
                clear use1 use2
            else
                
                % accuracy(tp)
                
                tmp=lda.choOptC(tp,cv).predict(1).predictedLabel==lda.choOptC(tp,cv).predict(1).trueLabel;
                accuracy(tp).choOptC.lda=[accuracy(tp).choOptC.lda;tmp];
                clear tmp
                tmp=nbl.choOptC(tp,cv).predict(1).predictedLabel==nbl.choOptC(tp,cv).predict(1).trueLabel;
                accuracy(tp).choOptC.nbl=[accuracy(tp).choOptC.nbl;tmp];
                clear tmp
                tmp=lda.choOptC(tp,cv).predict(2).predictedLabel==lda.choOptC(tp,cv).predict(2).trueLabel;
                accuracy(tp).choOptE.lda=[accuracy(tp).choOptE.lda;tmp];
                clear tmp
                tmp=nbl.choOptC(tp,cv).predict(2).predictedLabel==nbl.choOptC(tp,cv).predict(2).trueLabel;
                accuracy(tp).choOptE.nbl=[accuracy(tp).choOptE.nbl;tmp];
                clear tmp
                
                tmp=lda.choLRC(tp,cv).predict(1).predictedLabel==lda.choLRC(tp,cv).predict(1).trueLabel;
                accuracy(tp).choLRC.lda=[accuracy(tp).choLRC.lda;tmp];
                clear tmp
                tmp=nbl.choLRC(tp,cv).predict(1).predictedLabel==nbl.choLRC(tp,cv).predict(1).trueLabel;
                accuracy(tp).choLRC.nbl=[accuracy(tp).choLRC.nbl;tmp];
                clear tmp
                tmp=lda.choLRC(tp,cv).predict(2).predictedLabel==lda.choLRC(tp,cv).predict(2).trueLabel;
                accuracy(tp).choLRE.lda=[accuracy(tp).choLRE.lda;tmp];
                clear tmp
                tmp=nbl.choLRC(tp,cv).predict(2).predictedLabel==nbl.choLRC(tp,cv).predict(2).trueLabel;
                accuracy(tp).choLRE.nbl=[accuracy(tp).choLRE.nbl;tmp];
                clear tmp
                
                tmp=lda.EV1_high(tp,cv).predict.predictedLabel==lda.EV1_high(tp,cv).predict.trueLabel;
                accuracy(tp).EV1highCombo.lda=[accuracy(tp).EV1highCombo.lda;tmp];
                clear tmp
                tmp=nbl.EV1_high(tp,cv).predict.predictedLabel==nbl.EV1_high(tp,cv).predict.trueLabel;
                accuracy(tp).EV1highCombo.nbl=[accuracy(tp).EV1highCombo.nbl;tmp];
                clear tmp
                tmp=lda.EV2_high(tp,cv).predict.predictedLabel==lda.EV2_high(tp,cv).predict.trueLabel;
                accuracy(tp).EV2highCombo.lda=[accuracy(tp).EV2highCombo.lda;tmp];
                clear tmp
                tmp=nbl.EV2_high(tp,cv).predict.predictedLabel==nbl.EV2_high(tp,cv).predict.trueLabel;
                accuracy(tp).EV2highCombo.nbl=[accuracy(tp).EV2highCombo.nbl;tmp];
                clear tmp
                
                tmp=lda.ev1fr_corr(tp,cv).predict(1).predictedLabel==lda.ev1fr_corr(tp,cv).predict(1).trueLabel;
                accuracy(tp).ev1highC.lda=[accuracy(tp).ev1highC.lda;tmp];
                clear tmp
                tmp=nbl.ev1fr_corr(tp,cv).predict(1).predictedLabel==nbl.ev1fr_corr(tp,cv).predict(1).trueLabel;
                accuracy(tp).ev1highC.nbl=[accuracy(tp).ev1highC.nbl;tmp];
                clear tmp
                tmp=lda.ev1fr_corr(tp,cv).predict(2).predictedLabel==lda.ev1fr_corr(tp,cv).predict(2).trueLabel;
                accuracy(tp).ev1highE.lda=[accuracy(tp).ev1highE.lda;tmp];
                clear tmp
                tmp=nbl.ev1fr_corr(tp,cv).predict(2).predictedLabel==nbl.ev1fr_corr(tp,cv).predict(2).trueLabel;
                accuracy(tp).ev1highE.nbl=[accuracy(tp).ev1highE.nbl;tmp];
                clear tmp
                
                tmp=lda.ev2fr_corr(tp,cv).predict(1).predictedLabel==lda.ev2fr_corr(tp,cv).predict(1).trueLabel;
                accuracy(tp).ev2highC.lda=[accuracy(tp).ev2highC.lda;tmp];
                clear tmp
                tmp=nbl.ev2fr_corr(tp,cv).predict(1).predictedLabel==nbl.ev2fr_corr(tp,cv).predict(1).trueLabel;
                accuracy(tp).ev2highC.nbl=[accuracy(tp).ev2highC.nbl;tmp];
                clear tmp
                tmp=lda.ev2fr_corr(tp,cv).predict(2).predictedLabel==lda.ev2fr_corr(tp,cv).predict(2).trueLabel;
                accuracy(tp).ev2highE.lda=[accuracy(tp).ev2highE.lda;tmp];
                clear tmp
                tmp=nbl.ev2fr_corr(tp,cv).predict(2).predictedLabel==nbl.ev2fr_corr(tp,cv).predict(2).trueLabel;
                accuracy(tp).ev2highE.nbl=[accuracy(tp).ev2highE.nbl;tmp];
                clear tmp
                
                % posterior over true label
                
                use1=find(lda.choOptC(tp,cv).predict(1).trueLabel=='1');
                use2=find(lda.choOptC(tp,cv).predict(1).trueLabel=='2');
                tmp(use1,1)=lda.choOptC(tp,cv).predict(1).posterior(use1,1);
                tmp(use2,1)=lda.choOptC(tp,cv).predict(1).posterior(use2,2);
                post(tp).choOptC.lda=[post(tp).choOptC.lda;tmp];
                clear use1 use2 tmp
                use1=find(lda.choOptC(tp,cv).predict(2).trueLabel=='1');
                use2=find(lda.choOptC(tp,cv).predict(2).trueLabel=='2');
                tmp(use1,1)=lda.choOptC(tp,cv).predict(2).posterior(use1,1);
                tmp(use2,1)=lda.choOptC(tp,cv).predict(2).posterior(use2,2);
                post(tp).choOptE.lda=[post(tp).choOptE.lda;tmp];
                clear use1 use2 tmp
                
                use1=find(nbl.choOptC(tp,cv).predict(1).trueLabel=='1');
                use2=find(nbl.choOptC(tp,cv).predict(1).trueLabel=='2');
                tmp(use1,1)=nbl.choOptC(tp,cv).predict(1).posterior(use1,1);
                tmp(use2,1)=nbl.choOptC(tp,cv).predict(1).posterior(use2,2);
                post(tp).choOptC.nbl=[post(tp).choOptC.nbl;tmp];
                clear use1 use2 tmp
                use1=find(nbl.choOptC(tp,cv).predict(2).trueLabel=='1');
                use2=find(nbl.choOptC(tp,cv).predict(2).trueLabel=='2');
                tmp(use1,1)=nbl.choOptC(tp,cv).predict(2).posterior(use1,1);
                tmp(use2,1)=nbl.choOptC(tp,cv).predict(2).posterior(use2,2);
                post(tp).choOptE.nbl=[post(tp).choOptC.nbl;tmp];
                clear use1 use2 tmp
                
                use1=find(lda.choLRC(tp,cv).predict(1).trueLabel=='1');
                use2=find(lda.choLRC(tp,cv).predict(1).trueLabel=='2');
                tmp(use1,1)=lda.choLRC(tp,cv).predict(1).posterior(use1,1);
                tmp(use2,1)=lda.choLRC(tp,cv).predict(1).posterior(use2,2);
                post(tp).choLRC.lda=[post(tp).choLRC.lda;tmp];
                clear use1 use2 tmp
                use1=find(lda.choLRC(tp,cv).predict(2).trueLabel=='1');
                use2=find(lda.choLRC(tp,cv).predict(2).trueLabel=='2');
                tmp(use1,1)=lda.choLRC(tp,cv).predict(2).posterior(use1,1);
                tmp(use2,1)=lda.choLRC(tp,cv).predict(2).posterior(use2,2);
                post(tp).choLRE.lda=[post(tp).choLRE.lda;tmp];
                clear use1 use2 tmp
                
                use1=find(nbl.choLRC(tp,cv).predict(1).trueLabel=='1');
                use2=find(nbl.choLRC(tp,cv).predict(1).trueLabel=='2');
                tmp(use1,1)=nbl.choLRC(tp,cv).predict(1).posterior(use1,1);
                tmp(use2,1)=nbl.choLRC(tp,cv).predict(1).posterior(use2,2);
                post(tp).choLRC.nbl=[post(tp).choLRC.nbl;tmp];
                clear use1 use2 tmp
                use1=find(nbl.choLRC(tp,cv).predict(2).trueLabel=='1');
                use2=find(nbl.choLRC(tp,cv).predict(2).trueLabel=='2');
                tmp(use1,1)=nbl.choLRC(tp,cv).predict(2).posterior(use1,1);
                tmp(use2,1)=nbl.choLRC(tp,cv).predict(2).posterior(use2,2);
                post(tp).choLRE.nbl=[post(tp).choLRE.nbl;tmp];
                clear use1 use2 tmp
                
                use1=find(lda.EV1_high(tp,cv).predict.trueLabel=='false');
                use2=find(lda.EV1_high(tp,cv).predict.trueLabel=='true');
                tmp(use1,1)=lda.EV1_high(tp,cv).predict.posterior(use1,1);
                tmp(use2,1)=lda.EV1_high(tp,cv).predict.posterior(use2,2);
                post(tp).EV1highCombo.lda=[post(tp).EV1highCombo.lda;tmp];
                clear use1 use2 tmp
                use1=find(nbl.EV1_high(tp,cv).predict.trueLabel=='false');
                use2=find(nbl.EV1_high(tp,cv).predict.trueLabel=='true');
                tmp(use1,1)=nbl.EV1_high(tp,cv).predict.posterior(use1,1);
                tmp(use2,1)=nbl.EV1_high(tp,cv).predict.posterior(use2,2);
                post(tp).EV1highCombo.nbl=[post(tp).EV1highCombo.nbl;tmp];
                clear use1 use2 tmp
            
                use1=find(lda.EV2_high(tp,cv).predict.trueLabel=='false');
                use2=find(lda.EV2_high(tp,cv).predict.trueLabel=='true');
                tmp(use1,1)=lda.EV2_high(tp,cv).predict.posterior(use1,1);
                tmp(use2,1)=lda.EV2_high(tp,cv).predict.posterior(use2,2);
                post(tp).EV2highCombo.lda=[post(tp).EV2highCombo.lda;tmp];
                clear use1 use2 tmp
                use1=find(nbl.EV2_high(tp,cv).predict.trueLabel=='false');
                use2=find(nbl.EV2_high(tp,cv).predict.trueLabel=='true');
                tmp(use1,1)=nbl.EV2_high(tp,cv).predict.posterior(use1,1);
                tmp(use2,1)=nbl.EV2_high(tp,cv).predict.posterior(use2,2);
                post(tp).EV2highCombo.nbl=[post(tp).EV2highCombo.nbl;tmp];
                clear use1 use2 tmp
                
                use1=find(lda.ev1fr_corr(tp,cv).predict(1).trueLabel=='false');
                use2=find(lda.ev1fr_corr(tp,cv).predict(1).trueLabel=='true');
                tmp(use1,1)=lda.ev1fr_corr(tp,cv).predict(1).posterior(use1,1);
                tmp(use2,1)=lda.ev1fr_corr(tp,cv).predict(1).posterior(use2,2);
                post(tp).ev1highC.lda=[post(tp).ev1highC.lda;tmp];
                clear use1 use2 tmp
                use1=find(lda.ev1fr_corr(tp,cv).predict(2).trueLabel=='false');
                use2=find(lda.ev1fr_corr(tp,cv).predict(2).trueLabel=='true');
                tmp(use1,1)=lda.ev1fr_corr(tp,cv).predict(2).posterior(use1,1);
                tmp(use2,1)=lda.ev1fr_corr(tp,cv).predict(2).posterior(use2,2);
                post(tp).ev1highE.lda=[post(tp).ev1highE.lda;tmp];
                clear use1 use2 tmp
                
                use1=find(nbl.ev1fr_corr(tp,cv).predict(1).trueLabel=='false');
                use2=find(nbl.ev1fr_corr(tp,cv).predict(1).trueLabel=='true');
                tmp(use1,1)=nbl.ev1fr_corr(tp,cv).predict(1).posterior(use1,1);
                tmp(use2,1)=nbl.ev1fr_corr(tp,cv).predict(1).posterior(use2,2);
                post(tp).ev1highC.nbl=[post(tp).ev1highC.nbl;tmp];
                clear use1 use2 tmp
                use1=find(nbl.ev1fr_corr(tp,cv).predict(2).trueLabel=='false');
                use2=find(nbl.ev1fr_corr(tp,cv).predict(2).trueLabel=='true');
                tmp(use1,1)=nbl.ev1fr_corr(tp,cv).predict(2).posterior(use1,1);
                tmp(use2,1)=nbl.ev1fr_corr(tp,cv).predict(2).posterior(use2,2);
                post(tp).ev1highE.nbl=[post(tp).ev1highE.nbl;tmp];
                clear use1 use2 tmp
                
                use1=find(lda.ev2fr_corr(tp,cv).predict(1).trueLabel=='false');
                use2=find(lda.ev2fr_corr(tp,cv).predict(1).trueLabel=='true');
                tmp(use1,1)=lda.ev2fr_corr(tp,cv).predict(1).posterior(use1,1);
                tmp(use2,1)=lda.ev2fr_corr(tp,cv).predict(1).posterior(use2,2);
                post(tp).ev2highC.lda=[post(tp).ev2highC.lda;tmp];
                clear use1 use2 tmp
                use1=find(lda.ev2fr_corr(tp,cv).predict(2).trueLabel=='false');
                use2=find(lda.ev2fr_corr(tp,cv).predict(2).trueLabel=='true');
                tmp(use1,1)=lda.ev2fr_corr(tp,cv).predict(2).posterior(use1,1);
                tmp(use2,1)=lda.ev2fr_corr(tp,cv).predict(2).posterior(use2,2);
                post(tp).ev2highE.lda=[post(tp).ev2highE.lda;tmp];
                clear use1 use2 tmp
                
                use1=find(nbl.ev2fr_corr(tp,cv).predict(1).trueLabel=='false');
                use2=find(nbl.ev2fr_corr(tp,cv).predict(1).trueLabel=='true');
                tmp(use1,1)=nbl.ev2fr_corr(tp,cv).predict(1).posterior(use1,1);
                tmp(use2,1)=nbl.ev2fr_corr(tp,cv).predict(1).posterior(use2,2);
                post(tp).ev2highC.nbl=[post(tp).ev2highC.nbl;tmp];
                clear use1 use2 tmp
                use1=find(nbl.ev2fr_corr(tp,cv).predict(2).trueLabel=='false');
                use2=find(nbl.ev2fr_corr(tp,cv).predict(2).trueLabel=='true');
                tmp(use1,1)=nbl.ev2fr_corr(tp,cv).predict(2).posterior(use1,1);
                tmp(use2,1)=nbl.ev2fr_corr(tp,cv).predict(2).posterior(use2,2);
                post(tp).ev2highE.nbl=[post(tp).ev2highE.nbl;tmp];
                clear use1 use2 tmp
                
            end
        end
    end
    clear nbl lda
%     keyboard
    save([spath List(fn).name(1:end-13) 'decode.mat'], 'accuracy','post')
end


cd(spath)













