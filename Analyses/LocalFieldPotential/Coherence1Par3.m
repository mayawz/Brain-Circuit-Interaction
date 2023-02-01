% Coherence1Par3
% last used: Feb 25 2020

% code from Feb 09 only saved the high gamma band
% because cohgram=struct; declared a new struct with looping of every band

% can use [tapers,eigs]=dpsschk(tapers,N,Fs) to check tapers
% research gate users say default number [3 5] is ok
% this is used to calculate spectral gram, coherence, ect.
% from the actual data
close all; clear all; clc

% not run successfully

% /opt/local/matlab2017b/bin/matlab

addpath(genpath('/###################################################/chronux_2_12/'))

dpath='/###################################################/wrapped/';
fpath='/###################################################/Coherence/';
addpath(genpath(fpath))
cd(dpath)


%%

numWorkers=3;
% parallel setup
if isempty(gcp('nocreate'))
    poolobj=parpool('local',numWorkers);
end

%%
% parfor combinations=1:3
parfor tt=1:3
    combinations=10:12
    
    %%
    % there are 24 combinations of spk-field coherence in total
    % in the orders of
    % Pumbaa atOFFER all 6 combo, P@CHOICE, S@OFFER, S@CHOICE
    % condN={'SPKin-LFPpcc','SPKout-LFPpcc','SPKpcc-LFPin',...
    % 'SPKpcc-LFPout', 'SPKin-LFPout','SPKout-LFPin'};
    
    spath='/###################################################/results/Coherence/';
    
    spkorders=[2 4 6 6 2 4,...
        1 3 5 5 1 3 ,...
        8    10    12    12     8    10,...
        7     9    11    11     7     9];
    
    lfporders=[6 6 2 4 4 2,...
        5 5 1 3 3 1 ,...
        12    12     8    10    10     8,...
        11    11     7     9     9     7];
    
    
    condN={'SPKin-LFPpcc','SPKout-LFPpcc','SPKpcc-LFPin',...
        'SPKpcc-LFPout', 'SPKin-LFPout','SPKout-LFPin',...
        'SPKin-LFPpcc','SPKout-LFPpcc','SPKpcc-LFPin',...
        'SPKpcc-LFPout', 'SPKin-LFPout','SPKout-LFPin',...
        'SPKin-LFPpcc','SPKout-LFPpcc','SPKpcc-LFPin',...
        'SPKpcc-LFPout', 'SPKin-LFPout','SPKout-LFPin',...
        'SPKin-LFPpcc','SPKout-LFPpcc','SPKpcc-LFPin',...
        'SPKpcc-LFPout', 'SPKin-LFPout','SPKout-LFPin'};
    
    spkdpath='/###################################################/wrapped/SPKforCoherence/';
    
    cd(spkdpath)
    List_spk=dir( '*SPK*');
    tmp=load(List_spk(spkorders(combinations(tt))).name);
    spk=tmp.spk;
    
    lfpdpath='/###################################################/wrapped/LFPforCohCleaned/';
    
    cd(lfpdpath)
    List_lfp=dir( '*LFP*');
    tmp=load(List_lfp(lfporders(combinations(tt))).name);
    LFP=tmp.LFP;
    
    subj=List_spk(spkorders(combinations(tt))).name(1);
    condName=condN{combinations(tt)}
    aligned=List_spk(spkorders(combinations(tt))).name(end-11:end-4);
    
    %%
        
    % Instead of bandpass use frequency resolved FFT
    % NHP Theta usually is wider than rodent theta: 5-10 Hz, peaks around 8 Hz,
    % 3-5 cycles short-lived
    % Lower frequency longer bin, bin size should be based on center frequency
    % *3-5 cycle 10hz 500ms etc
    
    % Brain waves
    % Delta wave ? (0.5 ? 3 Hz)
    % Theta wave ? (4 ? 7 Hz)--> NHP(5-10 Hz)
    % Alpha wave ? (7 ? 15 Hz)
    % Mu wave ? (7.5 ? 12.5 Hz)
    % SMR wave ? (12.5 ? 15.5 Hz)
    % Beta wave ? (15 ? 30 Hz)
    % Gamma wave ? (>30 Hz)
    
    maxHz=[5;10;15;30;100];
    minHz=[0.5;6;11;16;31];
    % use 4 cycles
    mvgWinSizes=(1./maxHz)*4;
    
    cohgram=struct;
    
    for bandn=1:length(maxHz)
        
        disp(['Bandwidths ' num2str(bandn) '/' num2str(length(maxHz))])
        params=struct;
        params.Fs=1000; % sampling frequency Hz: cycles per second
        params.fpass=[0.5 100]; % band of frequencies to be kept
        % taper parameters default [3 5] for LFP
        % (higher the tapers lower the variance?)
        % example for S-L coherence was [10 19] strong smoothing of spk
        % example used [5 9] also looked good
        params.tapers=[5 9];
        params.pad=1; % pad factor for fft
        params.err=[1 0.05];
        params.trialave=1;
        
        movingwin=[mvgWinSizes(bandn) 0.05];
        % movingwin (in the form [window winstep] -- required for cohgramcpb
        % (in the form [window winstep] i.e length of moving window and step size)
        % Note units have to be consistent. Thus, if movingwin is in seconds, Fs
        % has to be in Hz. see chronux.m for more information.
        
        
        
        
        
        
        % since LFP and spikes are not of the same size, we use leave M out
        % cross-validation to equally subsample the same sizes of samples
        
        %         if size(LFP,2)>size(spk.psth,2)
        %             N=size(LFP,2);
        %             M=size(LFP,2)-size(spk.psth,2);
        %             [train,test] = crossvalind('LeaveMOut',N,1)
        %         elseif size(LFP,2)<size(spk.psth,2)
        %         else
        %         end
        % data1 (continuous data in form samples x trials) -- required
        % data2 (binned point process data in form samples x trials) -- required
        
        % ######## cohgramcpb for time*frequency*coherence heatmap #######
        if size(LFP,2)>= size(spk.psth,2)
            disp(['cohgramcpb: ' condName ' ' subj aligned ' 1/2'])
            [C1,phi1,S12_1,S1_1,S2_1,t1,f1,zerosp1,confC1,phierr1]=cohgramcpb(LFP(:,1:size(spk.psth,2)),spk.psth,movingwin,params);
            
            disp(['cohgramcpb: ' condName ' ' subj aligned ' 2/2'])
            [C2,phi2,S12_2,S1_2,S2_2,t2,f2,zerosp2,confC2,phierr2]=cohgramcpb(LFP(:,size(LFP,2)-size(spk.psth,2)+1:end),spk.psth,movingwin, params);
        else
            disp(['cohgramcpb: ' condName ' ' subj aligned ' 1/2'])
            [C1,phi1,S12_1,S1_1,S2_1,t1,f1,zerosp1,confC1,phierr1]=cohgramcpb(LFP,spk.psth(:,1:size(LFP,2)),movingwin,params);
            
            disp(['cohgramcpb: ' condName ' ' subj aligned ' 2/2'])
            [C2,phi2,S12_2,S1_2,S2_2,t2,f2,zerosp2,confC2,phierr2]=cohgramcpb(LFP,spk.psth(:,size(spk.psth,2)-size(LFP,2)+1:end),movingwin, params);
        end
        
        
        cohgram(bandn).C=(C1+C2)/2;
        cohgram(bandn).phi=(phi1+phi2)/2;
        cohgram(bandn).S12=(S12_1+S12_2)/2;
        cohgram(bandn).S1=(S1_1+S1_2)/2;
        cohgram(bandn).S2=(S2_1+S2_2)/2;
        cohgram(bandn).t=(t1+t2)/2;
        cohgram(bandn).f=(f1+f2)/2;
        cohgram(bandn).zerosp=(zerosp1+zerosp2)/2;
        cohgram(bandn).confC=(confC1+confC2)/2;
        cohgram(bandn).phierr=(phierr1+phierr2)/2;
        
    end

    %%
    % save([spath 'Coh_' subject '_' condName '.mat'], 'cohgram', 'coh', 'spectra', '-v7.3')
    fname=[spath 'Cohgram2_'  condName '_' subj aligned '.mat'];
    parsave(fname, cohgram)
    
    %     save(fname, 'cohgram', 'spectra','-v7.3')
    
    
end

%% shuts down workers


delete(poolobj)

%
function parsave(fname, cohgram)
save(fname, 'cohgram','-v7.3')
end




