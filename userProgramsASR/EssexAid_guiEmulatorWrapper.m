function [ op ] = EssexAid_guiEmulatorWrapper( speech, sampleRate, params )
%ESSEXAID_GUIEMULATORWRAPPER Summary of this function goes here
%   Detailed explanation goes here

% HOW TO GENERATE THE ALGO FOR THE JUCE IMPLEMENTATION
% 1) Make sure the the buffer size used corresponds to the name of the c
% files that you want to generate. For the VFrame version, ensure that the
% buffer size is max (6912)
% 2) Copy and paste the appropriate emlc line into the command window
% 3) Change the MOCcontrol memory allocation line in the appropriate
% function to reflect the buffer size.



% close all; clear all; clc


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DO NOT EDIT %
fNmax = 4/params.bwOct;
channelBFs = 250 * 2.^((0:fNmax)'*params.bwOct);
% END OF DO NOT EDIT %

mainGain = params.mainGain;

TCdBO = params.TCdBO;      %Compression thresholds (in dB OUTPUT from 2nd filt)
TMdBO = params.TMdBO;      %MOC thresholds (in dB OUTPUT from 2nd filt)

ARtau  = params.ARtau;            % decay time constant
ARthresholddB = params.ARthresholddB;      % dB SPL (input signal level) =>200 to disable

MOCtau = params.MOCtau; 
MOCfactor = params.MOCfactor;   %dB per dB OUTPUT

DRNLc = params.DRNLc;
bwOct = params.bwOct * ones(size(channelBFs)); %Octaves

filterOrder  = params.filterOrder; 

numSamples = params.numSamples; %MAX=6912


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DERRIVED CONSTANTS (DO NOT EDIT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numChannels = numel(channelBFs);
ARthresholdPa=  20e-6*10^(ARthresholddB/20);% Pa thresh for triggering AR
DRNLb    =     ( 2e-5 .* 10.^(TCdBO/20)) .^ (1-DRNLc)  ;

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
biggestNumSamples = numSamples; %6192 %Value comes form the maximum that juce can support


filterStatesL = (zeros(3000,1));
filterStatesR = filterStatesL;
filterCoeffs = (zeros(5000,1));

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
frameBufferL = buffer(speech(:,1), numSamples);
frameBufferR = buffer(speech(:,2), numSamples);
nFrames = size(frameBufferL,2);

pad = zeros(1,biggestNumSamples-numSamples);
ARampL=ones(1,biggestNumSamples);
ARampR = ARampL;
MOCcontrol = ones(numChannels, biggestNumSamples);

peakIPL = zeros(5,1);
peakOPL = peakIPL;
rmsIPL  = peakIPL;
rmsOPL  = peakIPL;

peakIPR = peakIPL;
peakOPR = peakIPL;
rmsIPR  = peakIPL;
rmsOPR  = peakIPL;



MOCend = zeros(numChannels,1);

op = [];
for nn = 1:nFrames
    frameBufferPadL = [frameBufferL(:,nn)' pad];
    frameBufferPadR = [frameBufferR(:,nn)' pad];

  [ outBufferL, outBufferR, filterStatesL, filterStatesR,  ARampL, ARampR, MOCend, peakIPL, peakOPL, rmsIPL, rmsOPL, peakIPR, peakOPR, rmsIPR, rmsOPR, MOCcontrol ] =...
    EssexAidProcessVFrameFBack( ...
    frameBufferPadL,...
    frameBufferPadR,...
    filterStatesL,...
    filterStatesR,...
    filterCoeffs,...
    numChannels,...
    numSamples,...
    ARampL,...
    ARampR,...
    ARthresholdPa,...
    filterOrder,...
    DRNLb,...
    DRNLc,...
    TMdBO,...
    MOCfactor,...
    peakIPL,...
    peakOPL,...
    rmsIPL,...
    rmsOPL,...
    peakIPR,...
    peakOPR,...
    rmsIPR,...
    rmsOPR,...
    MOCend,...
    MOCcontrol,...
    mainGain);
    

    outBuffer = ( [outBufferL(:, 1:numSamples); outBufferR(:, 1:numSamples)] );

    op = [op outBuffer];

end;


end

