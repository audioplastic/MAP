close all; clear all; clc

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
channelBFs= [250; 500; 1000; 2000; 4000;]; %MUST BE A COLUMN!!!!
% channelBFs= 250*2.^(-1:0.5:3.5)';
% channelBFs= [400; 800; 1600; 3200];
% channelBFs= [200; 400; 800; 1600; 3200];

TAspl = 50 * ones(size(channelBFs)); %Thresholds actual (in dB SPL)
TDspl = 10 * ones(size(channelBFs)); %Thresholds desired (in dB SPL)
TCspl = TDspl + 30;      %Compression thresholds (in dB SPL)
TMspl = TDspl + 10 ;     %MOC thresholds (in dB SPL)

ARtau  = 0.03;            % decay time constant
ARthresholddB = 80;      % dB SPL (input signal level) =>200 to disable

MOCtau = 0.1;  %0.06 
MOCfactor = 3e6;       % now that the conversions between velocity and displacement have been removed, internal values are 100 times greater (10kHz / 100 Hz) and so defaults need to be at least 100 times less than the old value of 2e9

DRNLc = 0.2 * ones(size(channelBFs));
bwOct =    1 * ones(size(channelBFs)); %Octaves

filterOrder  = 4; %This sounds better than 2nd order


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONSTANTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stapesScalar =	6e-008; % useful for calibrating the system
DRNLa1to1 = 1e4; %Default DRNLa calibrated to give unity gain

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DERRIVED CONSTANTS (DO NOT EDIT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numChannels = numel(channelBFs);
ARthresholdPa=  20e-6*10^(ARthresholddB/20);% Pa thresh for triggering AR

DRNLaBaseline = DRNLa1to1 * 10.^((TAspl-TDspl)/20);  
DRNLb    =   DRNLaBaseline  .*  (stapesScalar * 2e-5 .* 10.^(TCspl/20)) .^ (1-DRNLc)  ;
MOCthreshold =  min([...
    DRNLaBaseline *  stapesScalar *  2e-5 .* 10.^(TMspl/20)...
    DRNLb        .* (stapesScalar .* 2e-5 .* 10.^(TMspl/20)  ) .^ (DRNLc)...
    ], [], 2); 

opScaling_dB = 20*log10(  1  /  (DRNLa1to1*stapesScalar)  );


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ENUMERATIONS USED IN THE FRAME PROCESSOR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fake enumeration - must be kept up to date with juce enum
enumC_ARb = 0;
enumC_ARa = 2;
enumC_MOCb = 4;
enumC_MOCa = 6;

% enumC_BPb1 = 8;
% enumC_BPa1 = 13;
% enumC_BPb2 = 18;
% enumC_BPa2 = 23;
% enumC_BPb3 = 28;
% enumC_BPa3 = 33;
% enumC_BPb4 = 38;
% enumC_BPa4 = 43;

enumS_AR   = 0;

% enumS_MOC1  = 1;
% enumS_BPin_1_1 = 2;
% enumS_BPin_2_1 = 6;
% enumS_BPout_1_1 = 10;
% enumS_BPout_2_1 = 14;
% 
% enumS_MOC2 = 18;
% enumS_BPin_1_2 = 19;
% enumS_BPin_2_2 = 23;
% enumS_BPout_1_2 = 27;
% enumS_BPout_2_2 = 31;
% ...


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EMULATION OF THE GUI PARAMETER CONVERSIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numSamples = 14;
biggestNumSamples = numSamples; %6192 %Value comes form the maximum that juce can support


currentSpeechLevel = 75; %dB
[speech, sampleRate] = wavread('demo.wav');
speech = speech./sqrt(mean(speech.^2)); %Normalize RMS to 1
speech = speech * 20e-6*10^(currentSpeechLevel/20); %Convert RMS to pascals at desired level




filterStates = (zeros(300,1));
filterCoeffs = (zeros(500,1));

ipScalar = 1;
opScalar = 10^(opScaling_dB/20);



%filter coefficients
ARcutOff=1/(2*pi*ARtau);
[b,a] = butter(1,ARcutOff/(sampleRate/2));
filterCoeffs(enumC_ARb+1:enumC_ARb+2) = b;
filterCoeffs(enumC_ARa+1:enumC_ARa+2) = a;

MOCcutOff=1/(2*pi*MOCtau);
[bMOC,aMOC] = butter(1,MOCcutOff/(sampleRate/2));
filterCoeffs(enumC_MOCb+1:enumC_MOCb+2) = bMOC;
filterCoeffs(enumC_MOCa+1:enumC_MOCa+2) = aMOC;



for filterCount = 1:numChannels
    %-----------------------------------
    % nonlinear path - filter bws
    %-----------------------------------
    %Now defined in terms of octaves
    lowerCutOff=channelBFs(filterCount)*2^(-bwOct(filterCount)/2);
    upperCutOff=channelBFs(filterCount)*2^( bwOct(filterCount)/2);
    
    [b_DRNL,a_DRNL] = butter(2,[lowerCutOff upperCutOff]/(sampleRate/2));
    filterCoeffs(10*(filterCount-1)+9 :10*(filterCount-1)+13) = b_DRNL;
    filterCoeffs(10*(filterCount-1)+14:10*(filterCount-1)+18) = a_DRNL;
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EMULATION OF THE JUCE IO CALLBACK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
frameBuffer = buffer(speech, numSamples);
nFrames = size(frameBuffer,2);

pad = zeros(1,biggestNumSamples-numSamples);
ARamp=ones(1,biggestNumSamples);
MOCcontrol = ones(numChannels, biggestNumSamples);


op = zeros(1,biggestNumSamples);
tic
for nn = 1:nFrames
    frameBufferPadded = [frameBuffer(:,nn)' pad];
    
[ outBuffer, filterStates, ARamp, MOCcontrol ] = EssexAidProcess14Frame( ...    
    frameBufferPadded,...
    filterStates,...
    filterCoeffs,...
    numChannels,...
    numSamples,...
    ARamp,...
    MOCcontrol,...
    ARthresholdPa,...
    filterOrder,...
    DRNLaBaseline,...
    DRNLb,...
    DRNLc,...
    MOCthreshold,...
    MOCfactor,...
    stapesScalar);

    outBuffer = (outBuffer(1:numSamples));
    op(1+(nn-1)*numSamples:nn*numSamples) = outBuffer;
end;
toc

soundsc([speech' op], sampleRate)
% soundsc(op, sampleRate)

%% For emlc MEX
z = { ...    
    frameBufferPadded,...
    filterStates,...
    filterCoeffs,...
    numChannels,...
    numSamples,...
    ARamp,...
    MOCcontrol,...
    ARthresholdPa,...
    filterOrder,...
    DRNLaBaseline,...
    DRNLb,...
    DRNLc,...
    MOCthreshold,...
    MOCfactor,...
    stapesScalar};

% 6912 is the maximum number of frame samples accepted by JUCE Audio Editor
% 10 is a guess of the maximum number of channels we'll ever use

% USE THE FOLLOWING CODE TO COMPILE TO MEX
% emlc EssexAidProcessFrame -eg z -report

%% For emlc RTW:LIB
% 
% z = { ...    
%     emlcoder.egs(single(0),[1 6912]),...%frameBuffer(:,nn)',...
%     single(filterStates),...
%     single(filterCoeffs),...
%     int32(numChannels),...
%     emlcoder.egs(single(0),[1 6912]),...%ARamp,...
%     emlcoder.egs(single(0),[10 6912]),...%MOCcontrol,...
%     single(ARthresholdPa),...
%     int32(filterOrder),...
%     emlcoder.egs(single(0),10),...%DRNLaBaseline,...
%     emlcoder.egs(single(0),10),...%DRNLb,...
%     emlcoder.egs(single(0),10),...%DRNLc,...
%     emlcoder.egs(single(0),10),...%MOCthreshold,...
%     single(MOCfactor),...
%     single(stapesScalar),...
%     single(ipScalar),...
%     single(opScalar)};
% 
rtw_config = emlcoder.RTWConfig;
% rtw_config.TargetFunctionLibrary = 'C89/C90 (ANSI)'; %DO NOT USE prevents having to declare all the "extern" crap in C++
rtw_config.FilePartitionMethod = 'SingleFile';
rtw_config.GenCodeOnly = true;

% % % USE THE FOLLOWING TO COMPILE LIB
% % emlc EssexAidProcessFrame -eg z -launchreport -T rtw -s rtw_config


%% Filter initial coeffs
% clc
% for nn = 1:70
%     disp(  ['eml_filterCoeffs[' num2str(nn-1) '] = ' num2str(filterCoeffs(nn),'%0.10f') 'f;']  )
% end
    
    