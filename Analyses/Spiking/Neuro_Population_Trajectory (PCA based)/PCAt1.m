% PCAt1
% last used: Apr 28 2020
% use data from PrepPCAt2020

% what this does:
% conduct PCA on trial-averaged activity in correct trials only
% project trial-by-trial activity onto the top X PC space
% to counterbalance for sample/cell# size difference
% the number of top PCs is what captures >=70% of the variance
% calculate the distance indx
% dist(across condition)/dist(within condition)


% distIndx in OFCin and OFCout predict that in PCC and vise versa

% separete correct and incorrect trials and see whether that differentiates
% the predictions



clear all; close all; clc
dpath='/Volumes/GoogleDrive/My Drive/01Data/00 w_Anatomy/2020final/wrapped/PCAt_states/';
fpath='/Volumes/GoogleDrive/My Drive/01Data/00 w_Anatomy/2020final/code/';
spath='/Volumes/GoogleDrive/My Drive/01Data/00 w_Anatomy/2020final/results/PCAt/';
addpath(genpath(fpath))

cd(dpath)

List=dir('*states.mat*');

% bins=251:750; %choice at potentially 530
stepSize=5;
binStarts=251:stepSize:750-stepSize+1;
% Bwin=20;
% find(binStarts==301)%11
% find(binStarts==401)%31
% find(binStarts==501)%51
% find(binStarts==601)%71
% find(binStarts==701)%91
% find(binStarts==351)%21
% printNum=[11 31 51 61 71 91];
% mappedNum=[301 401 501 551 601 701];
plotSections=[1 find(binStarts==301) find(binStarts==401)...
    find(binStarts==501) find(binStarts==601) find(binStarts==701) 100];
plotColors1=[192,192,192;... % silver
    255,228,225;... % ll pink
    255,105,180;... % l pink
    255,0,0;... %  red
    178,34,34;... % dark red
    0,0,0]/255;% black
plotColors2=[192,192,192;... % silver
    176,196,222;... % ll blue
    70,130,180;... % l blue
    0,0,255;... % blue
    0,0,128;... % dark blue
    0,0,0]/255;% black
%%
for fn=1:length(List)
    %%
    clear clr c12 cLR C12 cev1 ev1hl
    fn=1
    cd(dpath)
    %%
    fname=List(fn).name
    load(fname);
    
    % ############    cholr       ############
    clr=[nanmean(cholr.cl,3);nanmean(cholr.cr,3)];
    [cLR.loadings,cLR.pcs,cLR.variance,~,cLR.explained,~] = pca(clr);
    cLR.cumExp=cumsum(cLR.explained);
    tmp=find(cLR.cumExp>70);
    cLR.top70=tmp(1); clear tmp
    if cLR.top70<3
        cLR.top70=3;
    end
    
%     figure(1)
%     plot(cLR.explained,'-s')
%     xlabel('Component')
%     ylabel('% Explained')
%     %     axis([0 10 0 50 ])
%     vline(cLR.top70)
%     txt=['    Cumulative%Explained = ' num2str(cLR.cumExp(cLR.top70))];
%     text(cLR.top70,cLR.explained(cLR.top70),txt)
    
% %     see PCA_Examples.m to see the reasoning
% %     distance from the training/averaged data correct trial
   
    cLR.top70PCs=centering(clr,1)*cLR.loadings(:,1:cLR.top70);

%     figure(2)
%     clear tmp1 tmp2
%     tmp1=cLR.top70PCs(1:100,1:3);
%     tmp2=cLR.top70PCs(101:200,1:3);
%     for sn=1:length(plotSections)-1
%         psec=plotSections(sn):plotSections(sn+1);
%         trajectoryPlot(tmp1(psec,:,:),plotColors1(sn,:),'o-')
%         hold on
%         trajectoryPlot(tmp2(psec,:,:),plotColors2(sn,:),'s-')
%     end
%     legend('ITI','ITI','Offer1','Offer1','Offer2','Offer2',...
%         'Choice','Choice','Reward','Reward','ITI','ITI')
%     xlabel('PC1')
%     ylabel('PC2')
%     zlabel('PC3')
%     title('Choice L vs R trl-avg')
%     hold off
%     cd(spath)

    elr=[nanmean(cholr.el,3);nanmean(cholr.er,3)];
    eLR.top70PCs=centering(elr,1)*cLR.loadings(:,1:cLR.top70);

%     figure(3) % just the error trials
%     clear tmp3 tmp4
%     tmp3=eLR.top70PCs(1:100,1:3);
%     tmp4=eLR.top70PCs(101:200,1:3);
%     
%     for sn=1:length(plotSections)-1
%         psec=plotSections(sn):plotSections(sn+1);
%         trajectoryPlot(tmp3(psec,:,:),plotColors1(sn,:),'o:')
%         hold on
%         trajectoryPlot(tmp4(psec,:,:),plotColors2(sn,:),'s:')
%     end
%     legend('ITI','ITI','Offer1','Offer1','Offer2','Offer2',...
%         'Choice','Choice','Reward','Reward','ITI','ITI')
%     xlabel('PC1')
%     ylabel('PC2')
%     zlabel('PC3')
%     title('Choice L vs R trl-avg ERROR trls')
%     hold off
%     cd(spath)
    
    % ############    cho12       ############
    c12=[nanmean(cho12.c1,3);nanmean(cho12.c2,3)];
    [C12.loadings,C12.pcs,C12.variance,~,C12.explained,~] = pca(c12);
    C12.cumExp=cumsum(C12.explained);
    tmp=find(C12.cumExp>70);
    C12.top70=tmp(1); clear tmp
    if C12.top70<3
        C12.top70=3;
    end
    
%     figure(1)
%     plot(C12.explained,'-s')
%     xlabel('Component')
%     ylabel('% Explained')
%     %     axis([0 10 0 50 ])
%     vline(C12.top70)
%     txt=['    Cumulative%Explained = ' num2str(C12.cumExp(C12.top70))];
%     text(C12.top70,C12.explained(C12.top70),txt)
    
    C12.top70PCs=centering(c12,1)*C12.loadings(:,1:C12.top70);    
    
    figure(2)
    clear tmp1 tmp2
    tmp1=C12.top70PCs(1:100,1:3);
    tmp2=C12.top70PCs(101:200,1:3);
    for sn=1:length(plotSections)-1
        psec=plotSections(sn):plotSections(sn+1);
        trajectoryPlot(tmp1(psec,:,:),plotColors1(sn,:),'o-')
        hold on
        trajectoryPlot(tmp2(psec,:,:),plotColors2(sn,:),'s-')
    end
    legend('ITI','ITI','Offer1','Offer1','Offer2','Offer2',...
        'Choice','Choice','Reward','Reward','ITI','ITI')
    xlabel('PC1')
    ylabel('PC2')
    zlabel('PC3')
    title('Choice 1 vs 2 trl-avg')
    hold off
    cd(spath)
    
    e12=[nanmean(cho12.e1,3);nanmean(cho12.e1,3)];
    E12.top70PCs=centering(e12,1)*C12.loadings(:,1:C12.top70);
    
%     figure(3) % just the error trials
%     clear tmp3 tmp4
%     tmp3=E12.top70PCs(1:100,1:3);
%     tmp4=E12.top70PCs(101:200,1:3);
%     
%     for sn=1:length(plotSections)-1
%         psec=plotSections(sn):plotSections(sn+1);
%         trajectoryPlot(tmp3(psec,:,:),plotColors1(sn,:),'o:')
%         hold on
%         trajectoryPlot(tmp4(psec,:,:),plotColors2(sn,:),'s:')
%     end
%     legend('ITI','ITI','Offer1','Offer1','Offer2','Offer2',...
%         'Choice','Choice','Reward','Reward','ITI','ITI')
%     xlabel('PC1')
%     ylabel('PC2')
%     zlabel('PC3')
%     title('Choice 1 vs 2 trl-avg ERROR trls')
%     hold off
%     cd(spath)
%     
    
    % ############    ev1hl       ############
    cev1=[nanmean(ev1hl.ch,3);nanmean(ev1hl.cl,3)];
    [cEV1.loadings,cEV1.pcs,cEV1.variance,~,cEV1.explained,~] = pca(cev1);
    cEV1.cumExp=cumsum(cEV1.explained);
    tmp=find(cEV1.cumExp>70);
    cEV1.top70=tmp(1); clear tmp
    if cEV1.top70<3
        cEV1.top70=3;
    end
    
%     figure(1)
%     plot(cEV1.explained,'-s')
%     xlabel('Component')
%     ylabel('% Explained')
%     %     axis([0 10 0 50 ])
%     vline(cEV1.top70)
%     txt=['    Cumulative%Explained = ' num2str(cEV1.cumExp(cEV1.top70))];
%     text(cEV1.top70,cEV1.explained(cEV1.top70),txt)
    
    cEV1.top70PCs=centering(cev1,1)*cEV1.loadings(:,1:cEV1.top70);
    
%     figure(2)
%     clear tmp1 tmp2
%     tmp1=cEV1.top70PCs(1:100,1:3);
%     tmp2=cEV1.top70PCs(101:200,1:3);
%     for sn=1:length(plotSections)-1
%         psec=plotSections(sn):plotSections(sn+1);
%         trajectoryPlot(tmp1(psec,:,:),plotColors1(sn,:),'o-')
%         hold on
%         trajectoryPlot(tmp2(psec,:,:),plotColors2(sn,:),'s-')
%     end
%     legend('ITI','ITI','Offer1','Offer1','Offer2','Offer2',...
%         'Choice','Choice','Reward','Reward','ITI','ITI')
%     xlabel('PC1')
%     ylabel('PC2')
%     zlabel('PC3')
%     title('EV1 h vs l trl-avg')
%     hold off
%     cd(spath)
    
    eev1=[nanmean(cho12.e1,3);nanmean(cho12.e1,3)];
    eEV1.top70PCs=centering(eev1,1)*cEV1.loadings(:,1:cEV1.top70);

%     figure(3) % just the error trials
%     clear tmp3 tmp4
%     tmp3=eEV1.top70PCs(1:100,1:3);
%     tmp4=eEV1.top70PCs(101:200,1:3);   
%     for sn=1:length(plotSections)-1
%         psec=plotSections(sn):plotSections(sn+1);
%         trajectoryPlot(tmp3(psec,:,:),plotColors1(sn,:),'o:')
%         hold on
%         trajectoryPlot(tmp4(psec,:,:),plotColors2(sn,:),'s:')
%     end
%     legend('ITI','ITI','Offer1','Offer1','Offer2','Offer2',...
%         'Choice','Choice','Reward','Reward','ITI','ITI')
%     xlabel('PC1')
%     ylabel('PC2')
%     zlabel('PC3')
%     title('EV1 h vs l trl-avg ERROR trls')
%     hold off
%     cd(spath)


save([dpath fname(1:end-11) '_PCAtrlAvg.mat'], 'cLR','eLR','C12','E12','cEV1','eEV1')
cd(dpath)

end % files fn



%%

cd(dpath)

