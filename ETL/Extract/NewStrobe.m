function NewStrobe( n )
%%% 20170217 MAM, PPB
%%% Uses NI USB6501 device connected to Plexon Digital input
global InitPL
% make sure to initialize global InitPl = []; at the start of the task code
% New edition done by Seng Bum Michael Yoo, 2016.01.15
if isempty(InitPL)
    InitPL = PL_DOInitDevice(1, 0); %(2,0)
%     InitPL = nan;
else
    
end

PL_DOGetDigitalOutputInfo;
if InitPL ~= 0;
    InitPL = PL_DOInitDevice(1, 0);
    
    [numDOCards, deviceNumbers, numBits, numLines] = PL_DOGetDigitalOutputInfo;
    [getDeviceStringResult, deviceString] = PL_DOGetDeviceString(1);
end
% DOInitDeviceResult = PL_DOInitDevice(2, 0);%This can only be called one time; has to be called before the strobe will
%be sent; if called another time the strobe will not work and you will have
%to restart Matlab and reinitialize the DIO

%% So I want to be able to call this again and again, but I only want to initialize it one time.
PL_DOSetWord( 1, 1, 15, n );% This is the event flag number you want to send over.
PL_DOPulseBit(1, 16, 0 );%This pin is connected to Pin 22(White/Black/Orange)and to Pin3 on the NI card, and functions as a Strobe
% memory; % Memory command is purposed for just checking.
end