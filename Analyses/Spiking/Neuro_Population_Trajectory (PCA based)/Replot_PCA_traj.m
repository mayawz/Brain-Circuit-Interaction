% R2R replot the population trajectory
% Last used Oct 2021

clear all; close all; clc
dpath='/Volumes/GoogleDrive/My Drive/01Data/00 w_Anatomy/2020final/wrapped/PCAt_states/';
fpath='/Volumes/GoogleDrive/My Drive/01Data/00 w_Anatomy/2020final/code/';
spath='/Volumes/GoogleDrive/My Drive/01Data/00 w_Anatomy/2020final/results/PCAt/R2R_oct2021/';
addpath(genpath(fpath))

cd(dpath)

% cd(spath)

%%
% load('P20171207_OFCin_PCAtrlAvg.mat');
% load('S20180731_OFCin_PCAtrlAvg.mat');

% load('P20171207_OFCout_PCAtrlAvg.mat');
% load('S20180731_OFCout_PCAtrlAvg.mat');

load('P20171207_PCC_PCAtrlAvg.mat');
% load('S20180731_PCC_PCAtrlAvg.mat');

% %%
% p = load('P20171207_OFCin_PCAtrlAvg.mat');
% s = load('S20180731_OFCin_PCAtrlAvg.mat');
% 
% in.cCho1 = [p.C12.top70PCs(1:100,1:3); s.C12.top70PCs(1:100,1:3)];
% in.cCho2 = [p.C12.top70PCs(101:200,1:3); s.C12.top70PCs(101:200,1:3)];
% in.cChoL = [p.cLR.top70PCs(1:100,1:3); s.C12.top70PCs(1:100,1:3)];
% in.cChoR = [p.cLR.top70PCs(101:200,1:3); s.C12.top70PCs(101:200,1:3)];
% 
% clear p s
% 
% p = load('P20171207_OFCout_PCAtrlAvg.mat');
% s = load('S20180731_OFCout_PCAtrlAvg.mat');
% 
% out.cCho1 = [p.C12.top70PCs(1:100,1:3); s.C12.top70PCs(1:100,1:3)];
% out.cCho2 = [p.C12.top70PCs(101:200,1:3); s.C12.top70PCs(101:200,1:3)];
% out.cChoL = [p.cLR.top70PCs(1:100,1:3); s.C12.top70PCs(1:100,1:3)];
% out.cChoR = [p.cLR.top70PCs(101:200,1:3); s.C12.top70PCs(101:200,1:3)];
% 
% clear p s
% 
% p = load('P20171207_PCC_PCAtrlAvg.mat');
% s = load('S20180731_PCC_PCAtrlAvg.mat');
% 
% pcc.cCho1 = [p.C12.top70PCs(1:100,1:3); s.C12.top70PCs(1:100,1:3)];
% pcc.cCho2 = [p.C12.top70PCs(101:200,1:3); s.C12.top70PCs(101:200,1:3)];
% pcc.cChoL = [p.cLR.top70PCs(1:100,1:3); s.C12.top70PCs(1:100,1:3)];
% pcc.cChoR = [p.cLR.top70PCs(101:200,1:3); s.C12.top70PCs(101:200,1:3)];

%%

stepSize=5;
binStarts=251:stepSize:750-stepSize+1;

Labels = {'Offer1','Offer2','Choice','Reward'};
% nums = [21, 41, 61, 81];

nums =[find(binStarts==331) find(binStarts==451)...
    find(binStarts==551) find(binStarts==681)];

% middle of 'Offer1','Offer2','Choice','Reward'
% 21 41 61 81

clear tmp1 tmp2
% ########### for each subject ###################
% tmp1=C12.top70PCs(1:100,1:3);
% tmp2=C12.top70PCs(101:200,1:3);

tmp1=cLR.top70PCs(1:100,1:3);
tmp2=cLR.top70PCs(101:200,1:3);

% ########### for 2 subjects combined ###################

% tmp1 = in.cCho1;
% tmp2 = in.cCho2;


% tmp1 = pcc.cChoL;
% tmp2 = pcc.cChoR;

%% CORRECT trials
close all

smo = 22; % smo = 25;

point11 = smooth(tmp1(:,1),smo);
point12 = smooth(tmp1(:,2),smo);
point13 = smooth(tmp1(:,3),smo);

point21 = smooth(tmp2(:,1),smo);
point22 = smooth(tmp2(:,2),smo);
point23 = smooth(tmp2(:,3),smo);

plot3(point11,point12,point13,'-','color','r','lineWidth',2.5);
hold on
plot3(point21,point22,point23,'-','color','b','lineWidth',2.5);
plot3(point11(nums),point12(nums),point13(nums),'s','color','r','MarkerFaceColor','r','MarkerSize',14);
plot3(point21(nums),point22(nums),point23(nums),'s','color','b','MarkerFaceColor','b','MarkerSize',14);



for tp = 1:length(Labels)
    text(point11(nums(tp)),point12(nums(tp)),point13(nums(tp)),['   ' Labels{tp}],'Color','red','FontSize',14);
    text(point21(nums(tp)),point22(nums(tp)),point23(nums(tp)),['   ' Labels{tp}],'Color','blue','FontSize',14);
    
end
% 
% xlim([-1.5 1.5])
% ylim([-1.5 1.5])
% zlim([-1.5 1.5])


xlabel('PC1')
ylabel('PC2')
zlabel('PC3')

subj = 'P_';
% subj = 'S_';

% filename = 'OFCin c12';
% filename = 'OFCin cLR';
% filename = 'OFCout c12';
% filename = 'OFCout cLR';
% filename = 'PCC c12';
filename = 'PCC cLR';

title(filename,'FontSize',16);
hold off

%%
% saveas(gcf,[spath subj filename],'tif')
print([spath subj filename],'-depsc','-painters');
cd(spath)
% 
    
 
 %% ERROR trials
 
clear tmp1 tmp2
% ########### for each subject ###################

tmp1=E12.top70PCs(1:100,1:3);
tmp2=E12.top70PCs(101:200,1:3);

% tmp1=eLR.top70PCs(1:100,1:3);
% tmp2=eLR.top70PCs(101:200,1:3);
% ########### for 2 subjects combined ###################


close all

smo = 22; % smo = 25;

point11 = smooth(tmp1(:,1),smo);
point12 = smooth(tmp1(:,2),smo);
point13 = smooth(tmp1(:,3),smo);

point21 = smooth(tmp2(:,1),smo);
point22 = smooth(tmp2(:,2),smo);
point23 = smooth(tmp2(:,3),smo);

plot3(point11,point12,point13,':','color','r','lineWidth',2.5);
hold on
plot3(point21,point22,point23,':','color','b','lineWidth',2.5);
plot3(point11(nums),point12(nums),point13(nums),'s','color','r','MarkerFaceColor','r','MarkerSize',14);
plot3(point21(nums),point22(nums),point23(nums),'s','color','b','MarkerFaceColor','b','MarkerSize',14);



for tp = 1:length(Labels)
    text(point11(nums(tp)),point12(nums(tp)),point13(nums(tp)),['   ' Labels{tp}],'Color','red','FontSize',14);
    text(point21(nums(tp)),point22(nums(tp)),point23(nums(tp)),['   ' Labels{tp}],'Color','blue','FontSize',14);
    
end

% xlim([-1.5 1.5])
% ylim([-1.5 1.5])
% zlim([-1.5 1.5])
xlabel('PC1')
ylabel('PC2')
zlabel('PC3')


% subj = 'P_';
subj = 'S_';
% filename = 'OFCin e12';
% filename = 'OFCin eLR';
% filename = 'OFCout e12';
% filename = 'OFCout eLR';
filename = 'PCC e12';
% filename = 'PCC eLR';

title(filename,'FontSize',16);
hold off

%%
saveas(gcf,[spath subj filename],'tif')
cd(spath)







