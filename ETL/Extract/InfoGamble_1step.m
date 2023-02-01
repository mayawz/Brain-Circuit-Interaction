function InfoGamble_1step(initial)
%ZMW 20170124


global eye;
global control;
global step;
global vars;
global color;
global visual;
global time;
global toplexon;
global counter;
global data;
rng('shuffle')
KbName('UnifyKeyNames'); 
Screen('Preference', 'VisualDebugLevel', 0);% warning('off','MATLAB:dispatcher:InexactCaseMatch')
Screen('Preference', 'Verbosity', 0);%Hides PTB Warnings

try
    [eye, control, step, vars, color, visual, time, toplexon, counter,data]=deal([]);
    
    % ###################### vars might need changes ######################
    %Initialize strobes -- old program
    control.strobesOn=0;
    
    % eye: eyetracking parameters
    eye.side = 2;%Tracked eye (L = 1, R = 2)
    % rewards
    vars.magnitude=[125 165 250]; % small=125; medium=165; large=250
    vars.rwdDuration=[0.090; 0.150; 0.230]*1.8; % juicer open duration
    color.loss= [255  0    0];     %Red - no reward
    color.rwd = [255 255 0; % yellow - small
        0    0    255;   %Blue - medium
        0    255  0];    %Green - large
    color.info=[0    255  255];    %cyan
    
    control.marketInfluence=1;
    % 0=no market, prob uniformly distributed
    % 1= good/bad/neutral
    % good market, prob=normrnd(0.7,0.08)
    % bad market, prob=normrnd(0.3,0.08)
    % neutral market, uniform
    if control.marketInfluence==1
        vars.blkSZ=randi([1,2],1,200);
        control.cumSumBlkSz=cumsum(vars.blkSZ);
        control.mktBlkChange=zeros(1,4000);
        control.mktBlkChange(control.cumSumBlkSz+1)=1;
        control.mktBlkChange(1)=1;
    end
    % ##################################################################
    
%     toplexon = NewStrobe;
    eye.fixating=0;%initialize the fixation 0=not 1=is fixated
    
    % counter: for sending strobes during eyetracking
    counter.fix1=0;
    counter.fix2=0;
    counter.choice=0;
    counter.type=0;
    
    % time: define event duration
    time.ITI=1; % ITI 1 second
    time.minFix=0.1; % minimum fixation duration 200ms
    time.trlTypeCue=0.5; % time showing/choose trl type
    time.offer=0.5; % 800ms offer presentation
    time.delay=0.5; % 200ms offer off/ delay
    time.feedback=0.15; % 150ms chosen option outlined
    time.outcome=1; % 300ms reveal outcome + reward delivery
    
    % vars: variables to be saved
    vars.trial=[];%trial number for stamping all other varialbles
    vars.daysTrials = 0;
    vars.offersPrst=[];% offers presented on this trial 1=safe, 2=gamble with small size,3 gamble with large size
    vars.firstLeft=[];% position of first offer, 1=left 0=right
    
    vars.size1=[];
    vars.size2=[];
    vars.prob1=[]; % probability for 1st offer rand(1)
    vars.prob2=[]; % probability for 2nd offer
    vars.ev1=[];
    vars.ev2=[];
    vars.win=[]; % win gamble=1; lose gamble=0; chose safe =2;
    vars.unchosenWin=999;
    vars.choseLeft=[];% chose which side: 1=left 0=right
    vars.outcome=[]; % actually received reward size 0=lost 125=choseSafe 165=winSmall 240=winLarge
    vars.correct=[]; % whether monkey chose the option with bigger EV, 1=correct; 0=incorrect
    vars.ProjName='E:\Matlab code\Data\InfoGamble';
    [vars.filename, vars.foldername] = createFile(vars.ProjName, 'InfoGamble_1step', initial);
    vars.daysTrials = countDayTrials(vars.ProjName, vars.foldername); %Count day's cumulative trials
    
    % color: initialize colors for visua stimuli
    color.white = WhiteIndex(0);
    color.black = BlackIndex(0);
    color.gray  = color.white/2;
    color.choiceOutline = color.white;  %Chosen option outline color white
    
    visual.baseRect = [0 0 120 340];
    % make a fixation square where monkeys can fixate on it with some
    % wiggle room around the dot
    visual.fixRect = [0 0 300 300];
    visual.infoRect= [0 0 45 45];
    % outline square for choice feedback
    visual.outline = [0 0 200 400];
    % trial Type rect
	visual.trlTypeRect = [0 0 200 150];
	%Connect to Eyelink
	%     visual.screen = setupEyelink;
	if ~Eyelink('IsConnected')
		Eyelink('initialize');%connects to eyelink computer
	end
	if ~Eyelink('IsConnected'), Eyelink('initialize');end % Connects to eyelink computer
    Eyelink('command','screen_pixel_coords = %ld %ld %ld %ld', 0, 0, 1920, 1080); %3440 1440        VR
    Eyelink('startrecording'); % Turns on the recording of eye position
    Eyelink('Command', 'calibration_type = HV5');
    Eyelink('Command', 'randomize_calibration_order = NO');
    Eyelink('Command', 'force_manual_accept = YES');
    Eyelink('StartSetup');
    
    % Send the keypress 'o' to put Eyelink in output mode
    Eyelink('SendKeyButton',double('o'),0,10);
    Eyelink('SendKeyButton',double('o'),0,10);
    
	% Open an on screen window
    [visual.window, windowRect] = Screen('OpenWindow', 1, color.black);
    % Get the size of the on screen window
    [screenXpixels, screenYpixels] = Screen('WindowSize', visual.window);
    % Get the centre coordinate of the window
    [visual.xCenter, visual.yCenter] = RectCenter(windowRect);
    % X coordinate for options
    visual.squareXpos = [screenXpixels * 0.25 screenXpixels * 0.75];
%   visual.squareXpos = [480 1440];
    % Query the frame duration
    visual.ifi = Screen('GetFlipInterval', visual.window);
    % Numer of frames to wait when specifying good timing
    visual.waitframes = 1;
    
    
    % Do dummy calls to GetSecs, WaitSecs, KbCheck to make sure they are
    % loaded and ready when we need them - without delays in the wrong
    % moment:
    KbCheck;
    WaitSecs(0.05);
    GetSecs;
    ShowCursor;
	
    % Maximum priority level
    topPriorityLevel = MaxPriority(visual.window);
    Priority(topPriorityLevel);
    
    
    %***** Ask to start *******
    go = 0;
    step=1;
    disp('Right Arrow to start');
    vars.StartTime=GetSecs;
    gokey=KbName('RightArrow');
    nokey=KbName('ESCAPE');
    while(go == 0)
        [keyIsDown,~,keyCode] = KbCheck;
        if keyCode(gokey)
            go = 1;
        elseif keyCode(nokey)
            go = -1;
            
        end
    end
    while keyIsDown
        [keyIsDown,~,~] = KbCheck;
    end
    home
    
    step=1;
    
    %***** Run trials *******
    while(go == 1)
        
        switch step
            case 1, step_ITI;
                
            case 2, step_fixaion1;
                
            case 3, step_trlTypeCue;
                
            case 4, step_chooseTrlType;
                
            case 5, step_offer1;
                
            case 6, step_delay1;
                
            case 7, step_offer2;
                
            case 8, step_delay2;
                
            case 9, step_fixation2;
                
            case 10, step_choice;
                
            case 11, step_computeOutcome;
                
            case 12, step_chosenReveal;
                
            case 13, step_allReveal;
                
        end
        
        go = keyCapture;
        
    end % of while-go
    sca
    
catch
    Screen('CloseAll');
    ShowCursor;
    fclose('all');
    Priority(0);
    close all;
    disp('end of session')
    % output the error message that describe the error
    
%         psychrethrow(psychlasterror);
end

end

% step 1
function step_ITI
global eye;global control;global step;global vars;global color;global visual;global time;
global toplexon; global counter;

Screen('FillRect', visual.window, color.black);
% Flip to the screen
Screen('Flip', visual.window);

if control.strobesOn==1, toplexon(6010);end %strobe: ITI -- start of trial
WaitSecs(time.ITI);


% ######## update the trial number for this trial ########
if isempty(vars.trial)
    vars.trial=1;
else
    vars.trial=vars.trial+1;
end

if(vars.trial ~= (vars.trial + vars.daysTrials))
    disp(['Current # ' num2str(vars.trial) '/' 'Cumulative #' num2str(vars.trial + vars.daysTrials)]);
    disp([])
else
    disp(['Trial #' num2str(vars.trial)]);
    disp([])
end

% ################ initialize paras for current trial ################
if rand(1)<=0.5
    vars.trlType=randi([1 2],1); %1=forced non-info 2=forced
else
    vars.trlType=3;
end

disp('Trial Type [1=forced non-info 2=forced info 3=freeChoice]')
disp(vars.trlType)
disp([])
% %%
% x = [0:.001:1];
% norm = betapdf(x,4,6.5);
% plot(x,norm,'b-')
% hold on
% x1 = [0:.001:1];
% norm1 = betapdf(x,6.5,4);
% plot(x1,norm1,'r-')
% %%
% for i=1:1000
%  y(i)=betarnd(3,6);
% end
% 
% hist(y)

if control.marketInfluence==0
    vars.prob1=rand(1);
    vars.prob2=rand(1);
    vars.market=0;
    disp('Market control off -- Neutral market')
    disp([])
    disp([])
else
    if control.mktBlkChange(vars.trial)==1
        vars.market=randi([-1 1],1);
        % -1=bad market; 0=neutral market; 3=good market;
        switch vars.market
            case -1
                vars.prob1=betarnd(4,6.5);
                vars.prob2=betarnd(4,6.5);
            case 0
                vars.prob1=rand(1);
                vars.prob2=rand(1);
            case 1
                vars.prob1=betarnd(6.5,4);
                vars.prob2=betarnd(6.5,4);
        end
    end
    disp('Market:-1=bad; 0=neutral; 1=good')
    disp(vars.market)
    disp([])
    disp([])
end


% strobe counter for fixation dot and choice where eyetracking is used
counter.fix1=0;
counter.fix2=0;
counter.choice=0;
counter.type=0;

% decide whether first offer is on the left or on the right
% position of first offer, 1=left 0=right
if rand(1)<=.5
    vars.firstLeft=1;
    visual.offerX1=visual.squareXpos(1);
    visual.offerX2=visual.squareXpos(2);
else
    vars.firstLeft=0;
    visual.offerX1=visual.squareXpos(2);
    visual.offerX2=visual.squareXpos(1);
end

% decide magnitude and probability based on vars.offersPrst
% vars.magnitude=[125 165 240]; % safe=125; small=165; large=240
temp=randperm(3);
vars.offersPrst=temp(1:2);
clear temp
vars.size1=vars.magnitude(vars.offersPrst(1));
vars.size2=vars.magnitude(vars.offersPrst(2));

color.offer1=color.rwd(vars.offersPrst(1),:);
color.offer2=color.rwd(vars.offersPrst(2),:);


step=2; %update step to proceed with the task
end

% step 2
function step_fixaion1
global eye;global control;global step;global vars;global color;global visual;
global time;global toplexon; global counter;


counter.fix1=counter.fix1+1;

% initialize a fixation square as an object with the same color of background
visual.fixObj=CenterRectOnPointd(visual.fixRect, visual.xCenter, visual.yCenter);
Screen('FillRect', visual.window, color.black, visual.fixObj);

visual.fixsqa=CenterRectOnPointd([0 0 20 20], visual.xCenter, visual.yCenter);

switch vars.market
    case -1
        % Screen('FillPoly',visual.window, color.white, [xPosVector; yPosVector], 1);
        xposvec=[visual.xCenter-20 visual.xCenter visual.xCenter+20  visual.xCenter visual.xCenter-20]';
        yposvec=[visual.yCenter visual.yCenter+30 visual.yCenter visual.yCenter-30 visual.yCenter]';
        Screen('FillPoly',visual.window, color.white, [xposvec yposvec], 1);
    case 0
        % Draw a white fixation dot where the mouse is
        Screen('DrawDots', visual.window, [visual.xCenter, visual.yCenter], 20, color.white, [], 2);
        
    case 1
        Screen('FillRect', visual.window, color.white, visual.fixsqa);
end


% Screen('DrawDots', visual.window, [visual.xCenter, visual.yCenter], 15, color.white, [], 2);


% Flip to the screen
Screen('Flip', visual.window);

if counter.fix1==1 && control.strobesOn==1, toplexon(6020);end % strobe: fixation 1 on screen

if rand(1)<=.5
    visual.typeX1=visual.squareXpos(1);
    visual.typeX2=visual.squareXpos(2);
else
    visual.typeX1=visual.squareXpos(2);
    visual.typeX2=visual.squareXpos(1);
end


% ######### Check eye position #########
e = Eyelink('newestfloatsample');

% e.gx(eye.side)
% e.gy(eye.side)
inside=IsInRect(e.gx(eye.side), e.gy(eye.side), visual.fixObj);
if inside==1
    if eye.fixating ~= 1 % this was initialized as 0
        eye.fixtime = GetSecs;
        eye.fixating = 1;
    elseif GetSecs >= (time.minFix + eye.fixtime)
        step=5;
        eye.fixating = 0;
    end
elseif eye.fixating == 1
    eye.fixating = 0;
end

end

% step 3
function step_trlTypeCue
global eye;global control;global step;global vars;global color;global visual;
global time;global toplexon; global counter;

% Sync us and get a time stamp
vbl = Screen('Flip', visual.window);
% Length of time and number of frames we will use for each drawing test
numSecs = time.trlTypeCue;
numFrames = round(numSecs / visual.ifi);

visual.cueObj=CenterRectOnPointd(visual.trlTypeRect, visual.xCenter, visual.yCenter);



if vars.trlType==1
    %     disp('Forced non-info trial')
    for frame = 1:numFrames
        
        % Draw the rect to the screen
        % Screen('FillRect', windowPtr [,color] [,rect] )
        Screen('DrawTexture', visual.window, visual.texture1, [], visual.cueObj, 0);
        
        % Flip to the screen
        vbl  = Screen('Flip', visual.window, vbl + (visual.waitframes - 0.5) * visual.ifi);
        
        % keycapture function in case we want to stop the task when still in this loop
        stopInWhile;
    end
    
elseif vars.trlType==2
    %     disp('Forced info trial')
    for frame = 1:numFrames
        
        % Draw the rect to the screen
        % Screen('FillRect', windowPtr [,color] [,rect] )
        Screen('DrawTexture', visual.window, visual.texture2, [], visual.cueObj, 0);
        % Flip to the screen
        vbl  = Screen('Flip', visual.window, vbl + (visual.waitframes - 0.5) * visual.ifi);
        
        % keycapture function in case we want to stop the task when still in this loop
        stopInWhile;
    end
    
end

step=5;
end

% step 5
function step_offer1
global eye;global control;global step;global vars;global color;global visual;
global time;global toplexon; global counter;

visual.Base1 = CenterRectOnPointd(visual.baseRect, visual.offerX1, visual.yCenter);
opt1Y=visual.yCenter+visual.baseRect(4)* (1-vars.prob1)/2;
visual.Opt1=CenterRectOnPointd(visual.baseRect.*[1 1 1 vars.prob1], visual.offerX1, opt1Y);
% Sync us and get a time stamp
vbl = Screen('Flip', visual.window);
% Length of time and number of frames we will use for each drawing test
numSecs = time.offer;
numFrames = round(numSecs / visual.ifi);
frame=1;

switch vars.trlType %1=forced non-info 2=forced info 3=freeChoice
    case 1
        vars.off1info=0;
        while frame <= numFrames
            Screen('FillRect', visual.window, color.loss, visual.Base1);
            Screen('FillRect', visual.window, color.offer1, visual.Opt1);
            
            stopInWhile;
            
            vbl  = Screen('Flip', visual.window, vbl + (visual.waitframes - 0.5) * visual.ifi);
            
            frame=frame+1;
        end
        
    case 2
        vars.off1info=1;
        visual.infoObj1=CenterRectOnPointd(visual.infoRect, visual.offerX1, visual.yCenter);
        while frame <= numFrames
            Screen('FillRect', visual.window, color.loss, visual.Base1);
            Screen('FillRect', visual.window, color.offer1, visual.Opt1);
            Screen('FillOval', visual.window,color.info, visual.infoObj1);
            
            stopInWhile;
            
            vbl  = Screen('Flip', visual.window, vbl + (visual.waitframes - 0.5) * visual.ifi);
            
            frame=frame+1;
        end
    case 3
        if rand(1)>0.5
            vars.off1info=1;
            visual.infoObj1=CenterRectOnPointd(visual.infoRect, visual.offerX1, visual.yCenter);
            while frame <= numFrames
                Screen('FillRect', visual.window, color.loss, visual.Base1);
                Screen('FillRect', visual.window, color.offer1, visual.Opt1);
                Screen('FillOval', visual.window,color.info, visual.infoObj1);
                
                stopInWhile;
                
                vbl  = Screen('Flip', visual.window, vbl + (visual.waitframes - 0.5) * visual.ifi);
                
                frame=frame+1;
            end
        else
            vars.off1info=0;
            while frame <= numFrames
                Screen('FillRect', visual.window, color.loss, visual.Base1);
                Screen('FillRect', visual.window, color.offer1, visual.Opt1);
                
                stopInWhile;
                
                vbl  = Screen('Flip', visual.window, vbl + (visual.waitframes - 0.5) * visual.ifi);
                
                frame=frame+1;
            end
        end
end

step =6;
end

% step 6
function step_delay1
global eye;global control;global step;global vars;global color;global visual;
global time;global toplexon; global counter;


% draw background screen, here we use black as background
Screen('FillRect', visual.window, color.black);
% Flip to the screen
Screen('Flip', visual.window);

WaitSecs(time.delay);

step=7;
end

% step 7
function step_offer2
global eye;global control;global step;global vars;global color;global visual;
global time;global toplexon; global counter;

visual.Base2 = CenterRectOnPointd(visual.baseRect, visual.offerX2, visual.yCenter);
opt2Y=visual.yCenter+visual.baseRect(4)* (1-vars.prob2)/2;
visual.Opt2=CenterRectOnPointd(visual.baseRect.*[1 1 1 vars.prob2], visual.offerX2, opt2Y);

% Sync us and get a time stamp
vbl = Screen('Flip', visual.window);
% Length of time and number of frames we will use for each drawing test
numSecs = time.offer;
numFrames = round(numSecs / visual.ifi);
frame=1;

switch vars.trlType %1=forced non-info 2=forced info 3=freeChoice
    case 1
        vars.off2info=0;
        while frame <= numFrames
            Screen('FillRect', visual.window, color.loss, visual.Base2);
            Screen('FillRect', visual.window, color.offer2, visual.Opt2);
            
            stopInWhile;
            
            vbl  = Screen('Flip', visual.window, vbl + (visual.waitframes - 0.5) * visual.ifi);
            
            frame=frame+1;
        end
        
    case 2
        vars.off2info=1;
        visual.infoObj2=CenterRectOnPointd(visual.infoRect, visual.offerX2, visual.yCenter);
        while frame <= numFrames
            Screen('FillRect', visual.window, color.loss, visual.Base2);
            Screen('FillRect', visual.window, color.offer2, visual.Opt2);
            Screen('FillOval', visual.window,color.info, visual.infoObj2);
            
            stopInWhile;
            
            vbl  = Screen('Flip', visual.window, vbl + (visual.waitframes - 0.5) * visual.ifi);
            
            frame=frame+1;
        end
    case 3
        if vars.off1info==0
           vars.off2info=1;
            visual.infoObj2=CenterRectOnPointd(visual.infoRect, visual.offerX2, visual.yCenter);
            while frame <= numFrames
                Screen('FillRect', visual.window, color.loss, visual.Base2);
                Screen('FillRect', visual.window, color.offer2, visual.Opt2);
                Screen('FillOval', visual.window,color.info, visual.infoObj2);
                
                stopInWhile;
                
                vbl  = Screen('Flip', visual.window, vbl + (visual.waitframes - 0.5) * visual.ifi);
                
                frame=frame+1;
            end
        else
            vars.off2info=0;
            while frame <= numFrames
                Screen('FillRect', visual.window, color.loss, visual.Base2);
                Screen('FillRect', visual.window, color.offer2, visual.Opt2);
                
                stopInWhile;
                
                vbl  = Screen('Flip', visual.window, vbl + (visual.waitframes - 0.5) * visual.ifi);
                
                frame=frame+1;
            end
        end
end

step =8;

end

% step 8
function step_delay2
global eye;global control;global step;global vars;global color;global visual;
global time;global toplexon; global counter;

% draw background screen, here we use black as background
Screen('FillRect', visual.window, color.black);
% Flip to the screen
Screen('Flip', visual.window);

WaitSecs(time.delay);

step=10;

end

% step 10
function step_choice
global eye;global control;global step;global vars;global color;global visual;
global time;global toplexon; global counter;

counter.choice=counter.choice+1;

% make a choice fixation rectangle as an object with back ground color
visual.choice1 = CenterRectOnPointd(visual.outline, visual.offerX1, visual.yCenter);
visual.choice2 = CenterRectOnPointd(visual.outline, visual.offerX2, visual.yCenter);
% draw both options on the screen
Screen('FillRect', visual.window, color.black, visual.choice1);
Screen('FillRect', visual.window, color.black, visual.choice2);
Screen('FillRect', visual.window, color.loss, visual.Base1);
Screen('FillRect', visual.window, color.offer1, visual.Opt1);
Screen('FillRect', visual.window, color.loss, visual.Base2);
Screen('FillRect', visual.window, color.offer2, visual.Opt2);
switch vars.trlType %1=forced non-info 2=forced info 3=freeChoice
    case 2
        Screen('FillOval', visual.window,color.info, visual.infoObj1);
        Screen('FillOval', visual.window,color.info, visual.infoObj2);
    case 3
        if vars.off1info==1
            Screen('FillOval', visual.window,color.info, visual.infoObj1);
        else
            Screen('FillOval', visual.window,color.info, visual.infoObj2);
        end
end
% Flip to the screen
Screen('Flip', visual.window);

if counter.choice==1
    if control.strobesOn==1, toplexon(6070);end %strobe: choice
end

t0=GetSecs;

%***** Check eye position *****
e = Eyelink('newestfloatsample');
inside1 = IsInRect(e.gx(eye.side), e.gy(eye.side), visual.choice1);
% inside1 = 0;
inside2 = IsInRect(e.gx(eye.side), e.gy(eye.side), visual.choice2);

if inside2==1
    if eye.fixating ~= 2
        eye.fixtime = GetSecs;
        eye.fixating = 2;
    elseif GetSecs >= (time.minFix + eye.fixtime)
        tEnd=GetSecs;
        vars.rt=tEnd-t0;
        vars.choice = 2;
        if vars.firstLeft==1
            vars.choseLeft=0;
        else
            vars.choseLeft=1;
        end
        
        % draw chosen outline and both options on the screen
        Screen('FillRect', visual.window, color.black, visual.choice1);
        Screen('FillRect', visual.window, color.choiceOutline, visual.choice2);
        Screen('FillRect', visual.window, color.loss, visual.Base1);
        Screen('FillRect', visual.window, color.offer1, visual.Opt1);
        Screen('FillRect', visual.window, color.loss, visual.Base2);
        Screen('FillRect', visual.window, color.offer2, visual.Opt2);
        switch vars.trlType %1=forced non-info 2=forced info 3=freeChoice
            case 2
                Screen('FillOval', visual.window,color.info, visual.infoObj1);
                Screen('FillOval', visual.window,color.info, visual.infoObj2);
            case 3
                if vars.off1info==1
                    Screen('FillOval', visual.window,color.info, visual.infoObj1);
                    vars.choseType=0;
                else
                    Screen('FillOval', visual.window,color.info, visual.infoObj2);
                    vars.choseType=1;
                end
        end
        % Flip to the screen
        Screen('Flip', visual.window);
        WaitSecs(time.feedback);
        
        step = 11;
        eye.fixating = 0;
        
    end
elseif eye.fixating == 2
    eye.fixating = 0;
    
end

if inside1==1
    if eye.fixating ~= 1
        eye.fixtime = GetSecs;
        eye.fixating = 1;
    elseif GetSecs >= (time.minFix + eye.fixtime)
        tEnd=GetSecs;
        vars.rt=tEnd-t0;
        vars.choice = 1;
        if vars.firstLeft==1
            vars.choseLeft=1;
        else
            vars.choseLeft=0;
        end
        
        % draw chosen outline and both options on the screen
        Screen('FillRect', visual.window, color.choiceOutline, visual.choice1);
        Screen('FillRect', visual.window, color.black, visual.choice2);
        Screen('FillRect', visual.window, color.loss, visual.Base1);
        Screen('FillRect', visual.window, color.offer1, visual.Opt1);
        Screen('FillRect', visual.window, color.loss, visual.Base2);
        Screen('FillRect', visual.window, color.offer2, visual.Opt2);
        switch vars.trlType %1=forced non-info 2=forced info 3=freeChoice
            case 2
                Screen('FillOval', visual.window,color.info, visual.infoObj1);
                Screen('FillOval', visual.window,color.info, visual.infoObj2);
            case 3
                if vars.off1info==1
                    Screen('FillOval', visual.window,color.info, visual.infoObj1);
                    vars.choseType=1;
                else
                    Screen('FillOval', visual.window,color.info, visual.infoObj2);
                    vars.choseType=0;
                end
        end
        % Flip to the screen
        Screen('Flip', visual.window);
        WaitSecs(time.feedback);
        
        step = 11;
        eye.fixating = 0;
        
    end
elseif eye.fixating == 1
    eye.fixating = 0;
end






end

% step 11
function step_computeOutcome
global eye;global control;global step;global vars;global color;global visual;
global time;global toplexon; global counter;

% vars.outcome=[];
% actually received reward size 0=lost 125=small
% 165=medium 240=large
% vars.offersPrst
% offers presented on this trial 1=small, 2=medium,3=large
temp=rand(1);
if vars.choice==1
    if temp>vars.prob1
        vars.win=0;
        vars.outcome=0;
        time.juicer=0;
    else
        vars.win=1;
        vars.outcome=vars.magnitude(vars.offersPrst(1));
        time.juicer=vars.rwdDuration(vars.offersPrst(1));
    end
else
    if temp>vars.prob2
        vars.win=0;
        vars.outcome=0;
        time.juicer=0;
    else
        vars.win=1;
        vars.outcome=vars.magnitude(vars.offersPrst(2));
        time.juicer=vars.rwdDuration(vars.offersPrst(2));
    end
end

vars.ev1=vars.prob1*vars.size1;
vars.ev2=vars.prob2*vars.size2;

if vars.choice==1
    if vars.ev1>vars.ev2
        vars.correct=1;
        disp ('Monkey is correct')
    else
        vars.correct=0;
        disp ('Monkey is incorrect')
    end
    
else
    if vars.ev1<vars.ev2
        vars.correct=1;
        disp ('Monkey is correct')
    else
        vars.correct=0;
        disp ('Monkey is incorrect')
    end
end

switch vars.trlType
    case 1 % forced non-info
        step=12;
    case 2 % forced info
        step=13;
    case 3
        if vars.choseType==1
            step=13;
        else
            step=12;
        end
end

end

% step 12
function step_chosenReveal
global eye;global control;global step;global vars;global color;global visual;
global time;global toplexon; global counter; global data;global juiceAO



% draw revealed chosen option on the screen
if vars.choice==1
    Screen('FillRect', visual.window, color.choiceOutline, visual.choice1);
    switch vars.win
        case 1
            Screen('FillRect', visual.window, color.offer1, visual.Base1);
        case 0
            Screen('FillRect', visual.window, color.loss, visual.Base1);
    end
    
else
    Screen('FillRect', visual.window, color.choiceOutline, visual.choice2);
    switch vars.win
        case 1
            Screen('FillRect', visual.window, color.offer2, visual.Base2);
        case 0
            Screen('FillRect', visual.window, color.loss, visual.Base2);
    end
    
end

% Flip to the screen
Screen('Flip', visual.window);

reward_digital_Juicer1(time.juicer);

WaitSecs(time.outcome);


data{vars.trial} = vars;

eval(['save ' vars.filename ' data']);




step=1;


end

% step 13
function step_allReveal
global eye;global control;global step;global vars;global color;global visual;
global time;global toplexon; global counter; global data;global juiceAO



% vars.unchosenWin: the result of unchosen gamble
temp=rand(1);
if vars.choice==1
    if temp>vars.prob2
        vars.unchosenWin=0;
    else
        vars.unchosenWin=1;
    end
else
    if temp>vars.prob1
        vars.unchosenWin=0;
    else
        vars.unchosenWin=1;
    end
end
% draw revealed both options on the screen
if vars.choice==1
    Screen('FillRect', visual.window, color.choiceOutline, visual.choice1);
    switch vars.win
        case 1
            Screen('FillRect', visual.window, color.offer1, visual.Base1);
        case 0
            Screen('FillRect', visual.window, color.loss, visual.Base1);
    end
    
    switch vars.unchosenWin
        case 1
            Screen('FillRect', visual.window, color.offer2, visual.Base2);
        case 0
            Screen('FillRect', visual.window, color.loss, visual.Base2);
    end
    
else
    Screen('FillRect', visual.window, color.choiceOutline, visual.choice2);
    switch vars.win
        case 1
            Screen('FillRect', visual.window, color.offer2, visual.Base2);
        case 0
            Screen('FillRect', visual.window, color.loss, visual.Base2);
    end
    
    switch vars.unchosenWin
        case 1
            Screen('FillRect', visual.window, color.offer1, visual.Base1);
        case 0
            Screen('FillRect', visual.window, color.loss, visual.Base1);
    end
    
end

% Flip to the screen
Screen('Flip', visual.window);

reward_digital_Juicer1(time.juicer);

WaitSecs(time.outcome);



data{vars.trial} = vars;

eval(['save ' vars.filename ' data']);


step=1;


end



function go = keyCapture
global step; global vars;
go = 1;
stopkey=KbName('ESCAPE');
pause=KbName('LeftArrow');
[keyIsDown,~,keyCode] = KbCheck;
if keyCode(stopkey)
    if control.strobesOn==1, toplexon(8003);end
    vars.EndTime=GetSecs;
    vars.TotalRunTime=(vars.EndTime-vars.StartTime)/60;
    disp('Total Time of Running in minutes')
    disp(vars.TotalRunTime)
    go = 0;
    Eyelink('Stoprecording');
elseif keyCode(pause) && step ~= 8
    step = 8;
elseif keyCode(pause) && step == 8
    step = 1;
    
end
while keyIsDown
    [keyIsDown,~,~] = KbCheck;
end
end

function reward_digital_Juicer1(rewardDuration)%MAM 20150706
%Changed by MAM 20160707 to use with Sesison-based Interface
%This function is to be used with the NI USB 6501 card.
%Pin 17/P0.0 (Juicer 1) and 25/GND(ground)
%Pin 18/P0.1 (Juicer 2) and 26/GND(ground)
warning('off','all');
% rewardDuration = 1;
%%%Move this into the running code so it initializes when you start the
%%%program.  The addline command will also turn on strobing capability.
%%%
s = daq.createSession('ni');
addDigitalChannel(s,'Dev3','Port0/line0:1','OutputOnly');

% outputSingleScan(s,[1 0])= juicer1
% outputSingleScan(s,[0 1])= juicer2
outputSingleScan(s,[1 0]);
tic;
while toc < rewardDuration;
end
outputSingleScan(s,[0 0]);

pause(0.001);


end

function [filename, foldername] = createFile(projFolder, projInits, initial) %CES 9/27/2013
% cd 'E:\Matlab code\Data'
cd(projFolder);
if exist(projFolder(7:end), 'dir')==0, mkdir(projFolder(7:end)); end
cd(projFolder);
dateS = datestr(now, 'yymmdd');
filename = [initial dateS '.1.' projInits '.mat'];
foldername = [initial dateS];
if exist(foldername, 'dir')==0, mkdir(foldername); end
cd(foldername)
trynum = 1;
while(exist(filename, 'file')~=0)
    trynum = trynum + 1;
    filename = [initial dateS '.' num2str(trynum) '.' projInits '.mat'];
end
home
end
  
function daysTrials = countDayTrials(projFolder, foldername) %CES 5/6/2013
cd(projFolder);
thesefiles = dir(foldername);
cd(foldername);
fileIndex = find(~[thesefiles.isdir]);
daysTrials = 0;
for i = 1:length(fileIndex)
    thisfile = thesefiles(fileIndex(i)).name;
    thisdata = importdata(thisfile);
    daysTrials = daysTrials + length(thisdata);
end
end

function stopInWhile
global vars; global onset; global step;
[keyIsDown,~,keyCode] = KbCheck;

if keyCode(KbName('ESCAPE'))
    vars.EndTime=GetSecs;
    vars.TotalRunTime=(vars.EndTime-vars.StartTime)/60;
    disp('Total Time of Running in minutes')
    disp(vars.TotalRunTime)
    onset = 0;
    step = 1;
    Screen('Preference', 'Verbosity', 0);%Hides PTB Warnings
    Gamepad('Unplug');
    sca;
    
    
end
end






