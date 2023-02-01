% DecodeForPosterior
% last used: Feb 10 2020

% get decoding posterior for
%              (1) choOpt: offer1 offer2;
%              (2) choLR: choLeft choRight
%              (3) EV1_high: EV1>=mean(EV1)
%              (4) EV2_high: EV2>=mean(EV2)
% separating correct and error trials
% training decoding with only correct trials except for EVs


clear all; close all; clc
% data path
dpath='/Users/mwang/Google Drive/01Data/00 w_Anatomy/2020final/wrapped/decodePatterns/';
% file path
fpath='/Users/mwang/Google Drive/01Data/00 w_Anatomy/2020final/code/Decode/';
addpath(genpath(fpath))
addpath(genpath('/Users/mwang/Google Drive/04Matlab/Functions'))
% save path
spath='/Users/mwang/Google Drive/01Data/00 w_Anatomy/2020final/results/Decoding/';

cd(dpath)
List=dir([dpath '*patterns*']);
for fn=1:length(List)
    load(List(fn).name)
    fname=List(fn).name(1:end-12)

    for tp=1:length(patt)
        
        % ###### EV1_high #########
        Xs=patt(tp).ev1fr;
        Ylabels=categorical(patt(tp).EV1_high);
        
        foldN=4;
        % indices = crossvalind('HoldOut',Ylabels,0.25);
        indices = crossvalind('KFold',length(Ylabels),4);
        for cv = 1:foldN
            test = (indices == cv);
            train = ~test;
            % ##### LDA #####
            % [predictStruct]=decoderRoutine(decoderName,...
            %    Xs_train,Ylabels_train,...
            %    Xs_test,Ylabels_test,...
            %    varargin)
            % decoderNameSet={'LDA','QDA','NBL','NBQ','SVMl','SVMn'};
            [predictStruct]=decoderRoutine('LDA',...
                Xs(train,:),Ylabels(train,:),...
                Xs(test,:),Ylabels(test,:));
            lda.EV1_high(tp,cv).predict=predictStruct;
            clear predictStruct
            
            [predictStruct]=decoderRoutine('NBL',...
                Xs(train,:),Ylabels(train,:),...
                Xs(test,:),Ylabels(test,:));
            nbl.EV1_high(tp,cv).predict=predictStruct;
            clear predictStruct 
        end
        clear Xs Ylabels
        
        % ###### EV2_high #########
        Xs=patt(tp).ev2fr;
        Ylabels=categorical(patt(tp).EV2_high);
    
        foldN=4;
        indices = crossvalind('KFold',length(Ylabels),4);
        
        for cv = 1:foldN
            test = (indices == cv);
            train = ~test;
            [predictStruct]=decoderRoutine('LDA',...
                Xs(train,:),Ylabels(train,:),...
                Xs(test,:),Ylabels(test,:));
            lda.EV2_high(tp,cv).predict=predictStruct;
            clear predictStruct
            
            [predictStruct]=decoderRoutine('NBL',...
                Xs(train,:),Ylabels(train,:),...
                Xs(test,:),Ylabels(test,:));
            nbl.EV2_high(tp,cv).predict=predictStruct;
            clear predictStruct
        end
        clear Xs Ylabels
        
        % ###### choOptC #########
        Xs=patt(tp).choOptCfr;
        Ylabels=categorical(patt(tp).choOptC);
        
        Xs_test2=patt(tp).choOptEfr;
        Ylabels_test2=categorical(patt(tp).choOptE);
        
        foldN=4;
        indices = crossvalind('KFold',length(Ylabels),4);
        
        for cv = 1:foldN
            test = (indices == cv);
            train = ~test;
            [predictStruct]=decoderRoutine('LDA',...
                Xs(train,:),Ylabels(train,:),...
                Xs(test,:),Ylabels(test,:),Xs_test2,Ylabels_test2);
            lda.choOptC(tp,cv).predict=predictStruct;
            clear predictStruct
            
            [predictStruct]=decoderRoutine('NBL',...
                Xs(train,:),Ylabels(train,:),...
                Xs(test,:),Ylabels(test,:),Xs_test2,Ylabels_test2);
            nbl.choOptC(tp,cv).predict=predictStruct;
            clear predictStruct
        end
        clear Xs Ylabels Xs_test2 Ylabels_test2

        % ###### choLRC #########
        Xs=patt(tp).choLRCfr;
        Ylabels=categorical(patt(tp).choLRC);
        
        Xs_test2=patt(tp).choLREfr;
        Ylabels_test2=categorical(patt(tp).choLRE);
        
        foldN=4;
        indices = crossvalind('KFold',length(Ylabels),4);
        
        for cv = 1:foldN
            test = (indices == cv);
            train = ~test;
            [predictStruct]=decoderRoutine('LDA',...
                Xs(train,:),Ylabels(train,:),...
                Xs(test,:),Ylabels(test,:),Xs_test2,Ylabels_test2);
            lda.choLRC(tp,cv).predict=predictStruct;
            clear predictStruct
            
            [predictStruct]=decoderRoutine('NBL',...
                Xs(train,:),Ylabels(train,:),...
                Xs(test,:),Ylabels(test,:),Xs_test2,Ylabels_test2);
            nbl.choLRC(tp,cv).predict=predictStruct;
            clear predictStruct
        end
        clear Xs Ylabels Xs_test2 Ylabels_test2
        
        % ###### ev1fr_corr #########
        Xs=patt(tp).ev1fr_corr;
        Ylabels=categorical(patt(tp).EV1_high_corr);
        
        Xs_test2=patt(tp).ev1fr_erro;
        Ylabels_test2=categorical(patt(tp).EV1_high_erro);
        
        foldN=4;
        indices = crossvalind('KFold',length(Ylabels),4);
        
        for cv = 1:foldN
            test = (indices == cv);
            train = ~test;
            [predictStruct]=decoderRoutine('LDA',...
                Xs(train,:),Ylabels(train,:),...
                Xs(test,:),Ylabels(test,:),Xs_test2,Ylabels_test2);
            lda.ev1fr_corr(tp,cv).predict=predictStruct;
            clear predictStruct
            
            [predictStruct]=decoderRoutine('NBL',...
                Xs(train,:),Ylabels(train,:),...
                Xs(test,:),Ylabels(test,:),Xs_test2,Ylabels_test2);
            nbl.ev1fr_corr(tp,cv).predict=predictStruct;
            clear predictStruct
        end
        clear Xs Ylabels Xs_test2 Ylabels_test2
        
        % ###### ev2fr_corr #########
        Xs=patt(tp).ev2fr_corr;
        Ylabels=categorical(patt(tp).EV2_high_corr);
        
        Xs_test2=patt(tp).ev2fr_erro;
        Ylabels_test2=categorical(patt(tp).EV2_high_erro);
        
        foldN=4;
        indices = crossvalind('KFold',length(Ylabels),4);
        
        for cv = 1:foldN
            test = (indices == cv);
            train = ~test;
            [predictStruct]=decoderRoutine('LDA',...
                Xs(train,:),Ylabels(train,:),...
                Xs(test,:),Ylabels(test,:),Xs_test2,Ylabels_test2);
            lda.ev2fr_corr(tp,cv).predict=predictStruct;
            clear predictStruct
            
            [predictStruct]=decoderRoutine('NBL',...
                Xs(train,:),Ylabels(train,:),...
                Xs(test,:),Ylabels(test,:),Xs_test2,Ylabels_test2);
            nbl.ev2fr_corr(tp,cv).predict=predictStruct;
            clear predictStruct
        end
        clear Xs Ylabels Xs_test2 Ylabels_test2
    end % tp=time point
    
    save([spath fname 'posterior.mat'],'lda','nbl')
    
end % fn=file number from List

%%
cd(spath)










