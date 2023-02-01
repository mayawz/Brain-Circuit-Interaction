% PCAt3
% last used: Apr 28 2020
% use data from PrepPCAt2020 & PCAt1
% project trial-by-trial activity onto the top X PC space
% to counterbalance for sample/cell# size difference
% the number of top PCs is what captures >=70% of the variance
% calculate the distance indx
% dist(fn)(across condition)/dist(fn)(within condition)

% dist(fn)Indx in OFCin and OFCout predict that in PCC and vise versa

% separete correct and incorrect trials and see whether that differentiates
% the predictions
clear all; close all; clc

dpath='/Users/mwang/Google Drive/01Data/00 w_Anatomy/2020final/wrapped/PCAt_states/';
fpath='/Users/mwang/Google Drive/01Data/00 w_Anatomy/2020final/code/Trajectory/';
spath='/Users/mwang/Google Drive/01Data/00 w_Anatomy/2020final/results/PCAt/';
addpath(genpath(fpath))
cd(dpath)

stepSize=5;
binStarts=251:stepSize:750-stepSize+1;

%%
cd(dpath)
useNames={'OFCin','OFCout','PCC'};
for nn=1:length(useNames)
    clear List1 List2 dist
    List1=dir(['*' useNames{nn} '*PCA*']);
    List2=dir(['*' useNames{nn} '*states*']);
    for fn=1:length(List1)
        clear ps ds
        ps=load(List1(fn).name);
        ds=load(List2(fn).name);
        % #################### 12 C ###############################
        t1=ds.cho12.c1Trl;
        t2=ds.cho12.c2Trl;
        d1=ds.cho12.c1;
        d2=ds.cho12.c2;
        pdata=ps.C12;
        
        % keyboard
        [Dst,Dsp]=projDist(t1,t2,d1,d2,pdata);
        td(fn).c12.Dst=Dst;
        td(fn).c12.Dsp=Dsp;
        clear t2 t2 d1 d2 pdata Dst Dsp
        
        % #################### 12 E ###############################
        t1=ds.cho12.e1Trl;
        t2=ds.cho12.e2Trl;
        d1=ds.cho12.e1;
        d2=ds.cho12.e2;
        pdata=ps.C12;
        
        [Dst,Dsp]=projDist(t1,t2,d1,d2,pdata);
        td(fn).e12.Dst=Dst;
        td(fn).e12.Dsp=Dsp;
        clear t2 t2 d1 d2 pdata Dst Dsp
        
        % #################### lr  C ###############################
        t1=ds.cholr.clTrl;
        t2=ds.cholr.crTrl;
        d1=ds.cholr.cl;
        d2=ds.cholr.cr;
        pdata=ps.cLR;
        
        [Dst,Dsp]=projDist(t1,t2,d1,d2,pdata);
        td(fn).clr.Dst=Dst;
        td(fn).clr.Dsp=Dsp;
        clear t2 t2 d1 d2 pdata Dst Dsp
        
        % #################### lr  E ###############################
        t1=ds.cholr.elTrl;
        t2=ds.cholr.erTrl;
        d1=ds.cholr.el;
        d2=ds.cholr.er;
        pdata=ps.cLR;
        
        [Dst,Dsp]=projDist(t1,t2,d1,d2,pdata);
        td(fn).elr.Dst=Dst;
        td(fn).elr.Dsp=Dsp;
        clear t2 t2 d1 d2 pdata Dst Dsp
        
        % #################### ev C ###############################
        t1=ds.ev1hl.chTrl;
        t2=ds.ev1hl.clTrl;
        d1=ds.ev1hl.ch;
        d2=ds.ev1hl.cl;
        pdata=ps.cEV1;
        
        [Dst,Dsp]=projDist(t1,t2,d1,d2,pdata);
        td(fn).cev1.Dst=Dst;
        td(fn).cev1.Dsp=Dsp;
        clear t2 t2 d1 d2 pdata Dst Dsp
        
        % #################### ev E ###############################
        t1=ds.ev1hl.ehTrl;
        t2=ds.ev1hl.elTrl;
        d1=ds.ev1hl.eh;
        d2=ds.ev1hl.el;
        pdata=ps.cEV1;
        
        [Dst,Dsp]=projDist(t1,t2,d1,d2,pdata);
        td(fn).eev1.Dst=Dst;
        td(fn).eev1.Dsp=Dsp;
        clear t2 t2 d1 d2 pdata Dst Dsp
        
        td(fn).memo={'trl_dist, td(row=subj P=1 S=2).condition.Dst OR .Dsp'};
    end
    
    
    save([spath useNames{nn} 'DistTbT.mat'], 'td')
    
end
%%
cd(spath)

