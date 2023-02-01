% Coherence1Par
% last used: Feb 9 2020
% can use [tapers,eigs]=dpsschk(tapers,N,Fs) to check tapers
% research gate users say default number [3 5] is ok
% this is used to calculate spectral gram, coherence, ect.
% from the actual data
close all; clear all; clc

% /opt/local/matlab2017b/bin/matlab

addpath(genpath('####################################################/chronux_2_12/'))

dpath='####################################################/wrapped/';
fpath='####################################################/Coherence/';
addpath(genpath(fpath))
cd(dpath)


%%
% numWorkers=12;
% parallel setup
% if isempty(gcp('nocreate'))
%     poolobj=parpool('local',numWorkers);
% end

lfpdpath='###################################################/wrapped/LFPforCoherence/';
cd(lfpdpath)
List_lfp=dir( '*LFP*');
%%
% parfor combinations=1:length(List_lfp)
% for combinations=1:length(List_lfp)  
for combinations=9:10 
% combinations=2
    
    spath='/###################################################/LFPforCohCleaned/';

    cd(lfpdpath)
    List_lfp=dir( '*LFP*');
    tmp=load(List_lfp(combinations).name);
    lfp=tmp.lfp;
    
    
    %% remove low frequency movement artifacts
    
    disp('dlfp')
    
    dlfp=locdetrend(lfp.uV,1000,[.1 .05]);
    
    
    
    
    %% remove 60Hz harmonics
    disp('60Hz')
    
    params=struct;
    params.tapers=[5 9]; % taper parameters
    params.Fs=1000; % sampling frequency Hz: cycles per second
    params.pad=0; % pad factor for fft
    plt='y';
    LFP=rmlinesc(dlfp,params,plt);
    
%%
disp('save cleaned lfp')

fname=[spath List_lfp(combinations).name(1:end-4) '_cleaned.mat'];

save(fname, 'LFP','-v7.3')

%%
cd(spath)   
end



%% shuts down workers
% delete(poolobj)

%%
% function parsave(fname, var1)
%     save(fname, var1,'-v7.3')
% end




