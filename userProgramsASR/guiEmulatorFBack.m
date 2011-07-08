% HOW TO GENERATE THE ALGO FOR THE JUCE IMPLEMENTATION
% 1) Make sure the the buffer size used corresponds to the name of the c
% files that you want to generate. For the VFrame version, ensure that the
% buffer size is max (6912)
% 2) Copy and paste the appropriate emlc line into the command window
% 3) Change the MOCcontrol memory allocation line in the appropriate
% function to reflect the buffer size.



close all; clear all; clc


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bwOct = 1/1;

% DO NOT EDIT %
fNmax = 4/bwOct;
channelBFs = 250 * 2.^((0:fNmax)'*bwOct);
% END OF DO NOT EDIT %

mainGain = 1 * ones(size(channelBFs));

TCdBO = 2150 * ones(size(channelBFs));      %Compression thresholds (in dB OUTPUT from 2nd filt)
TMdBO = 30 * ones(size(channelBFs));      %MOC thresholds (in dB OUTPUT from 2nd filt)

ARtau  = 0.03;            % decay time constant
ARthresholddB = 150;      % dB SPL (input signal level) =>200 to disable

MOCtau = 0.3; 
MOCfactor = 0.5;   %dB per dB OUTPUT

DRNLc = 0.2 * ones(size(channelBFs));
bwOct = bwOct * ones(size(channelBFs)); %Octaves

filterOrder  = 2; 


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
numSamples = 6912; %MAX=6912
biggestNumSamples = numSamples; %6192 %Value comes form the maximum that juce can support


% currentSpeechLevel = 90; %dB
% [speech, sampleRate] = wavread('demoM.wav');
% 
% speech = speech./sqrt(mean(speech.^2)); %Normalize RMS to 1
% noiseL = 0.2 * randn(size(speech));
% noiseR = 0.2 * randn(size(speech));
% 
% speech = [speech speech];
% 
% speech = speech * 20e-6*10^(currentSpeechLevel/20); %Convert RMS to pascals at desired level

%%%%% Sine pulse code %%%%%
sampleRate = 48e3;
pulseDur = 0.1;
dt = 1/sampleRate;
tAxis = dt:dt:pulseDur;
freq = 2000;
sPulse = sin(2*pi*freq*tAxis);
sPulse = sPulse./sqrt(mean(sPulse.^2));
rms2dBspl = @(dBspl)20e-6*10^(dBspl/20);
zPadDur = 0.5;
zPad = zeros(1,ceil(sampleRate*zPadDur));
speech = [  sPulse*rms2dBspl(20)  zPad... 
            sPulse*rms2dBspl(40)  zPad... 
            sPulse*rms2dBspl(60)  zPad... 
            sPulse*rms2dBspl(80)  zPad... 
            sPulse*rms2dBspl(100) zPad  ];
% figure; plot(speech);
speech = [speech' speech'];
%%%%% END of sine pulse code %%%%%



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


MOCobs = [];
MOCend = zeros(numChannels,1);

op = [];
tic
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
    

    MOCobs = [MOCobs MOCend];
    outBuffer = ( [outBufferL(:, 1:numSamples); outBufferR(:, 1:numSamples)] );
%     op(:, 1+(nn-1)*numSamples:nn*numSamples) = outBuffer;
    op = [op outBuffer];
%     peakOPsaved = [peakOPsaved peakOP];
%     peakIPsaved = [peakIPsaved peakIP];
%     rmsOPsaved = [rmsOPsaved rmsOP];
%     rmsIPsaved = [rmsIPsaved rmsIP];
end;
toc


figure; plot(20*log10(MOCobs'))
20*log10(sqrt(mean(op(1,:).^2)) / 2e-5);

ipdb = 20*log10(abs(speech/20e-6)+(1/(2^32)));
opdb = 20*log10(abs(op/20e-6)+(1/(2^32)));
figure; plot(ipdb(:,1),'k'); hold on; plot(opdb(1,:),'r');
ylim([0 100])

% figure
% subplot(2,1,1)
% plot(peakIPsaved')
% title('INPUT')
% subplot(2,1,2)
% plot(peakOPsaved'-80)
% title('OUTPUT')
% 
% figure
% subplot(2,1,1)
% plot(rmsIPsaved')
% title('INPUT')
% subplot(2,1,2)
% plot(rmsOPsaved'-80)
% title('OUTPUT')


speechPly = speech(:,1)/sqrt(mean(speech(:,1).^2));
opPly = op(1,:)/sqrt(mean(op(1,:).^2));
soundsc([speechPly; opPly'], sampleRate)
% soundsc(op, sampleRate)

%% For emlc MEX
z = { ...    
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
    mainGain};

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
% % emlc EssexAidProcessVFrameFBack -eg z -launchreport -T rtw -s rtw_config


%% Filter initial coeffs
% clc
% for nn = 1:70
%     disp(  ['eml_filterCoeffs[' num2str(nn-1) '] = ' num2str(filterCoeffs(nn),'%0.10f') 'f;']  )
% end
    
    