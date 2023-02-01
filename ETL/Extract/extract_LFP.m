% extract_LFP

clear all; close all; clc; fclose('all');

cd('/home/jeeves-raid2/benh-data/Maze/Raw/Spock/Spock/PCC_OFC13/')



% openNSx('S2018073110002.ns2')
%
% openNEV('P201712070002.nev')
% %%
% strobes=NEV.Data.SerialDigitalIO.UnparsedData;
%
% %%
% NS2.ElectrodesInfo(31:end).ElectrodeID

%%
% This is a script template for accessesing different data from a Trellis
% saved data set.  If you are using a data set that utilizes pausing, the
% time between data will be padded with ZEROS.

dataType = 'LFP'; % dataType can be:
% Neural - 'Spike','LFP', 'Hi-Res', or 'Raw';
% Stim - 'Stim'
% Analog I/O - 'Analog 30k' or 'Analog 1k'
% Digital I/O - 'Digital'

% dataChannel = 129;   % A-1 starts at 1, A-2 starts at 33,
% B-1 starts at 129
% Analog I/O starts at 10241
% Digital is 1 to 4 or 'parallel'
% dataChannel = 'parallel';

plotStatus = true; % Whether or not to plot the data


% completeFilePath = 'Y:\Isaac\SampleGrapevineData\SampleGrapevineData.nev';
completeFilePath = 'S2018073110002.ns2';

% Open the file and extract some basic information
[ns_status, hFile] = ns_OpenFile(completeFilePath)

%%
for d=1:192
    % Determine correct entityID for desired datastream
    dataChannel=d+128;
    
    switch dataType
        
        case {'Spike' 'LFP' 'Hi-Res' 'Raw' 'Analog 30k' 'Analog 1k'}
            EntityIndices = find([hFile.Entity(:).ElectrodeID] == dataChannel);
            for i = 1:length(EntityIndices)
                fileTypeNum = hFile.Entity(EntityIndices(i)).FileType;
                fileType = hFile.FileInfo(fileTypeNum).Type;
                switch dataType
                    case 'Spike'
                        if strcmp('nev', fileType); entityID = EntityIndices(i); break; end
                    case {'LFP' 'Analog 1k'}
                        if strcmp('ns2', fileType); entityID = EntityIndices(i); break; end
                    case 'Hi-Res'
                        if strcmp('nf3', fileType); entityID = EntityIndices(i); break; end
                    case {'Raw' 'Analog 30k'}
                        if strcmp('ns5', fileType); entityID = EntityIndices(i); break; end
                end
            end
            
        case 'Stim'
            entityID = find([hFile.Entity(:).ElectrodeID] == dataChannel + 5120);
            
        case 'Digital'
            if isnumeric(dataChannel)
                entityID = find(cellfun(@strcmpi, {hFile.Entity.Reason},...
                    repmat({['Input Ch ' num2str(dataChannel)]},...
                    size({hFile.Entity.Reason}))));
            else
                %             entityID = find(cellfun(@strcmpi, {hFile.Entity.Reason},...
                %                 repmat({'Digital Input'}, size({hFile.Entity.Reason}))));
                entityID = find(cellfun(@strcmpi, {hFile.Entity.Reason},...
                    repmat({'Parallel Input'}, size({hFile.Entity.Reason}))));
            end
    end
    
    % Extract channel info
    [ns_RESULT, entityInfo] = ns_GetEntityInfo(hFile, entityID(end));
    
    % Extract data
    switch dataType
        case 'Spike'
            for i = 1:entityInfo.ItemCount
                if i == 1; [ns_RESULT, nsSegmentSourceInfo] = ns_GetSegmentSourceInfo(hFile, entityID, i); end;
                [ns_RESULT, spikeEventTime_s(i), spikeWindowData(:,i), sample_count(i), unit_id(i)] = ns_GetSegmentData(hFile, entityID, i);
            end
            % Sorted spikes
            uniqueUnits = unique(unit_id(unit_id~=0));
            for i = 1:length(uniqueUnits)
                SortedSpikeData{i} = spikeWindowData(:,unit_id==uniqueUnits(i));
            end
            
            if plotStatus
                figure(); hold on;
                if isempty(uniqueUnits)
                    for i = 1:size(spikeWindowData,1); plot(spikeWindowData(:,i),'b'); end
                    xlabel('Sample'); ylabel('Spike Data (\muV)');
                else
                    for i = 1:length(uniqueUnits); plot(SortedSpikeData{i}); end
                end
                hold off
            end
            
            
        case {'LFP' 'Hi-Res' 'Raw' 'Analog 30k' 'Analog 1k'}
            [ns_RESULT, analogInfo] = ns_GetAnalogInfo(hFile, entityID(end));     % analog info contains things like range and sampling rate
            
            TimeStamps = hFile.FileInfo(hFile.Entity(entityID).FileType).TimeStamps;
            numSamples = sum(TimeStamps(:,end));
            analogInputData = zeros(1,numSamples);
            startIndex = 1;
            indexCount = TimeStamps(2,1);
            for i = 1:size(TimeStamps,2)
                [~, ~, tempData] = ns_GetAnalogData(hFile, entityID, startIndex, indexCount);
                dataRange = TimeStamps(1,i) + (1:TimeStamps(2,i));
                analogInputData(dataRange) = tempData';
                clear tempData
                if i ~= size(TimeStamps,2)
                    startIndex = startIndex + TimeStamps(2,i);
                    indexCount = TimeStamps(2,i+1);
                end
            end
            analogInputDataTime_s = (0:numSamples-1)' ./ analogInfo.SampleRate;
            
            if strcmp(dataType, 'Analog')
                analogInputData_mV = analogInputData; clear analogInputData;
                if plotStatus; plot(analogInputDataTime_s,analogInputData_mV); xlabel('Time (s)'); ylabel ([dataType ' (mV)']); end;
            else
                analogInputData_uV = analogInputData; clear analogInputData;
                if plotStatus; plot(analogInputDataTime_s,analogInputData_uV); xlabel('Time (s)'); ylabel ([dataType ' (uV)']); end;
            end
            
        case 'Stim'
            % work through it backwards to setup variables first
            for i = hFile.Entity(entityID).Count:-1:1
                [ns_RESULT, TimeStamp(i), Data(i,:), SampleCount] = ns_GetSegmentData(hFile, entityID, i);
                if i == hFile.Entity(entityID).Count
                    time = 0:1/3e4:TimeStamp(hFile.Entity(entityID).Count)+51/3e4;
                    stimData = zeros(size(time));
                end
                Ts = find(time < TimeStamp(i),1,'last');
                stimData(Ts+1:Ts+52) = Data(i,:);
            end
            if plotStatus; figure(); plot(time,stimData); end
            
        case 'Digital'
            % Get events and time stamps
            numCount = entityInfo.ItemCount;
            data = NaN(1, numCount); timeStamps = NaN(1, numCount); sz = NaN(1, numCount);
            for i = 1:numCount
                [~, timeStamps(i), data(i), dataSize(i)] = ns_GetEventData(hFile, entityID, i);
            end
    end
    
    %%
    if d==1
        LFP.time_Stamp=analogInputDataTime_s;
        LFP.uV(d,:)=analogInputData_uV;
    else
        LFP.uV(d,:)=analogInputData_uV;
    end
    
    disp(d)
end

save([completeFilePath(1:9) 'LFP.mat'],'LFP','-v7.3');

%%
% 1:16:192 = 1    17    33    49    65    81    97   113   129   145   161   177
tmp1= LFP.time_Stamp;
tmp2 = LFP.uV;

clearvars -except tmp1 tmp2 

%%
clear LFP
LFP.timeStamp=tmp1;
% 49:64 81:96 113:128 129:160 161:192

% LFP.uV=tmp2([49:64 81:96 113:128],:);
LFP.uV=tmp2(81:96,:);


save('S20180731_PCC_Array2_LFP.mat','LFP')

% save('S20180731_PCC_LFP.mat','LFP','-v7.3')






