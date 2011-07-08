function [ op ] = EssexAid_wrapper( inputSignal, sampleRate, param )
%ESSEXAID_WRAPPER Summary of this function goes here
%   Detailed explanation goes here


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONSTANTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stapesScalar =	6e-008; % useful for calibrating the system
DRNLa1to1 = 1e4; %Default DRNLa calibrated to give unity gain

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DERRIVED CONSTANTS (DO NOT EDIT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numChannels = numel(param.channelBFs);
ARthresholdPa=  20e-6*10^(param.ARthresholddB/20);% Pa thresh for triggering AR

DRNLaBaseline = DRNLa1to1 * 10.^((param.TAspl-param.TDspl)/20);  
DRNLb    =   DRNLaBaseline  .*  (stapesScalar * 2e-5 .* 10.^(param.TCspl/20)) .^ (1-param.DRNLc)  ;
MOCthreshold =  min([...
    DRNLaBaseline *  stapesScalar *  2e-5 .* 10.^(param.TMspl/20)...
    DRNLb        .* (stapesScalar .* 2e-5 .* 10.^(param.TMspl/20)  ) .^ (param.DRNLc)...
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

% enumS_AR   = 0;

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
numSamples = 1024;

% currentSpeechLevel = 75; %dB
% [speech, sampleRate] = wavread('demo.wav');
% speech = speech./sqrt(mean(speech.^2)); %Normalize RMS to 1
% speech = speech * 20e-6*10^(currentSpeechLevel/20); %Convert RMS to pascals at desired level



ARamp=ones(1,numSamples);

filterStates = (zeros(300,1));
filterCoeffs = (zeros(500,1));

ipScalar = 1;
opScalar = 10^(opScaling_dB/20);

MOCcontrol = ones(numChannels, numSamples);


%filter coefficients
ARcutOff=1/(2*pi*param.ARtau);
[b,a] = butter(1,ARcutOff/(sampleRate/2));
filterCoeffs(enumC_ARb+1:enumC_ARb+2) = b;
filterCoeffs(enumC_ARa+1:enumC_ARa+2) = a;

MOCcutOff=1/(2*pi*param.MOCtau);
[bMOC,aMOC] = butter(1,MOCcutOff/(sampleRate/2));
filterCoeffs(enumC_MOCb+1:enumC_MOCb+2) = bMOC;
filterCoeffs(enumC_MOCa+1:enumC_MOCa+2) = aMOC;



for filterCount = 1:numChannels
    %-----------------------------------
    % nonlinear path - filter bws
    %-----------------------------------
    %Now defined in terms of octaves
    lowerCutOff=param.channelBFs(filterCount)*2^(-param.bwOct(filterCount)/2);
    upperCutOff=param.channelBFs(filterCount)*2^( param.bwOct(filterCount)/2);
    
    [b_DRNL,a_DRNL] = butter(2,[lowerCutOff upperCutOff]/(sampleRate/2));
    filterCoeffs(10*(filterCount-1)+9 :10*(filterCount-1)+13) = b_DRNL;
    filterCoeffs(10*(filterCount-1)+14:10*(filterCount-1)+18) = a_DRNL;
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EMULATION OF THE JUCE IO CALLBACK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

frameBuffer = buffer(inputSignal, numSamples);
nFrames = size(frameBuffer,2);

biggestNumSamples = numSamples;
pad = zeros(1,biggestNumSamples-numSamples);

op = zeros(size(inputSignal));

for nn = 1:nFrames
    frameBufferPadded = [frameBuffer(:,nn)' pad];
    
    [ outBuffer, filterStates, ARamp, MOCcontrol ] =...
    EssexAidProcessVFrame( ...
    frameBufferPadded,...
    filterStates,...
    filterCoeffs,...
    numChannels,...
    numSamples,...
    ARamp,...
    MOCcontrol,...
    ARthresholdPa,...
    param.filterOrder,...
    DRNLaBaseline,...
    DRNLb,...
    param.DRNLc,...
    MOCthreshold,...
    param.MOCfactor,...
    stapesScalar);
    
%     op = (op(1:numSamples));
    op(1+(nn-1)*numSamples:nn*numSamples) = outBuffer;
%     op = [op outBuffer];
end

op = (op(1:length(inputSignal)));


end

