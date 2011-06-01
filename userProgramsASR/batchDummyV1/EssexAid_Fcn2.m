% Function version of EssexAid_COPY2
% Nov 2010

function output = EssexAid_Fcn2(inputSignal, sampleRate, param)
% close all; clear all; clc

% Another major hack in Feb 2011
% The whole thing converted to worki in terms of velocity. I spent a long
% time formulating parameters around the conversions between displacement
% and velocity to find in the end that there were no real advantages of
% using a displacement measure in the aid. Hahahaha - live and learn.

%MOC threshold now needs to be different in each channel if we want to
%specify it as a meaningful dB value

% % % Fixed a bug ARthresholdPa=  28e-6*10^(ARthresholddB/20); ->20e-6

% EssexAid after Nick has hacked away at it even more:
% * Noticed that GUI analytical I/O functions and the I/O functions here are
% not in agreement -> RMS was calculated over whole sample for each input
% level (including silence) in original EssexAid.m. Now rms is only
% calculated in the energetic region. Results from both GUI and this
% function now agree. YIPPEE!
% * I also changed the (potentially hazardous) "if MOC>0" section. This was
% a hangover from the days when the MOC control was more simple. It has
% been converted to play nicely with the time varying MOC function. If any
% MOC <=0 then DRNLa is set to zero for the whole segment. I think this was
% the original intention

% EssexAid after Nick has hacked away at it:
% * Now has hp filter at op to convert back to pressure
% * Single gain stage at the end: opScaling
%       - hearingAidAmp parameter removed
% * Filter BWs now specified in octaves
% * DRNLaBaseline can be tweaked in each chan re. DRNLaBaseline
% * DRNLb can be tweaked in each chan re. DRNLb

% EssexAid
% Complete hearing aid simulation test bed.
% Sounds are either locally generated (signaltype='tone')
%   or read from file (signaltype='file')
%  For simplicity all sound is sampled at 44100
%  All sounds are mixed with background noise read from file specified
%   in *backgroundFileName* (also 44100 Hz).
%  Speech and noise signalLevels are in dB SPL
%
% Aid filters can be either gammatone or Butterworth (useGammatone= 1/0)
%
% Outputs are saved as normalised vectors in 'outputSummary.wav'
%  these can be played back (playSounds=1)
%  at a volume set in hearingAidGain (dB)
% The output is effectively a hypothesis as to the response of a
%  normal compressive basilar membrane suitable for input to a
%  the ear of a patient with no compression, acoustic reflex
%  or MOC reflex regulation.

% parameters
% param.opScaling_dB=110;         % 70;
%
% param.ARthresholddB=60;         % 60
% param.ARtau=0.03;               % 0.03
%
% param.DRNLaBaseline=1e4;        % 1e4
% param.DRNLb=0.5e-5;             % 0.5e-5 ?? --> YES but misleading. This = 5e-6 as in the readme (watch out fot the 0.5 vs 5)
% param.DRNLc=0.2;
%
% param.MOCfactor=2e9;            % 2e9
% param.MOCtau=0.06;              % 0.06;
%
% param.doHpFilt                    = 1;    %TRUE
% param.bwOct                       = [4/3 4/3 4/3 4/3]; %Octaves
% param.DRNLaBaselineRelative_dB    = [0 0 0 0]; % [0 0 0 0];
% param.DRNLbRelative_dB            = [0 0 0 0]; % [0 0 0 0];
% param.filterOrder                 = 3; %3
%
%
% numChannels=4;
% minBF=250; 	maxBF= 4000;
% param.channelBFs=round(logspace(log10(minBF), log10(maxBF), numChannels));

% **  **  **  **  **  **  **  **  **  **  **  **   aid parameters
opScaling_dB=param.opScaling_dB;      % dB (BM displacement -> Pa), applied to output

%                             **stapes**
stapesScalar=	6e-008; % useful for calibrating the system

%                          **Acoustic reflex**
% This monitors recent sound levels and applies attenuation to
%   continuously pulls signal input level down to the AR threshold level
ARtau=param.ARtau;            % decay time constant
ARthresholddB=param.ARthresholddB;       % dB SPL (input signal level) =1 to disable
ARthresholdPa=  20e-6*10^(ARthresholddB/20);% Pa thresh for triggering AR

DRNLaBaseline=param.DRNLaBaseline;% * 10.^(param.DRNLaBaselineRelative_dB/20)';
DRNLb=param.DRNLb; %* 10.^(param.DRNLbRelative_dB/20);%repmat(DRNLb,1, numChannels); % DRNLb is a vector: one per channel

c=param.DRNLc;                  % compression exponent (fixed)

%                              **MOC**
% Within-channel attenuation operates as AGC at relatively low levels
MOCtau=param.MOCtau;            % decay time constant
MOCfactor=param.MOCfactor;          % arbitrary strength factor
efferentDelay= 0.010;   % feedback delay (also controls segment duration)

dt=1/sampleRate;
duration=dt*length(inputSignal);

%                                     segment realated information
% The signal is divided into 10 ms segments to allow for
% feedback signals to be updated.
% 10 ms is approximately the mandatory delay imposed
%  by delays in the nervous system
%
%  earObject length should be a multiple of  segment length
%   (trim any excess signal)
segmentDuration=efferentDelay;
segmentLength=round(segmentDuration/dt);
nSignalSegments=floor(duration/(segmentDuration));
trimmedLength=nSignalSegments*segmentLength;
inputSignal=inputSignal (1:trimmedLength);

%                                             pre-allocation
numChannels = numel(param.channelBFs);
ARrmsBndry=[];
NLg1bndry=cell(1,numChannels);
NLg2bndry=cell(1,numChannels);
smoothedRMSbndry=cell(1,numChannels);
filteredSignal=zeros(numChannels, length(inputSignal));
MOCcontrol=repmat(DRNLaBaseline,1,segmentLength);

%% signal processing begins here
% NB all units are international units and are kept close to
% realistic values for the purpose of comparing aid response to
% BM response. This will not be necessesary in the hardware
for segmentStart=1:segmentLength:length(inputSignal)
    segmentEnd=segmentStart+segmentLength-1;
    prevSegStart=segmentStart-segmentLength;
    prevSegEnd=segmentEnd-segmentLength;
    inputThisSeg=inputSignal(segmentStart:segmentEnd);
    
    %                                                 Acoustic Reflex
    % find  rms of smoothed input signal
    %  this will be used to trigger the AR reflex
    cutOff=1/(2*pi*ARtau);
    % square and smooth
    [y  ARrmsBndry]=  ...
        UTIL_lowPassFilter(inputThisSeg.^2, cutOff, 1, dt, ARrmsBndry);
    % restore Pa scale
    smoothedARrms(segmentStart:segmentEnd)= y.^0.5;
    
    if prevSegStart>0   % not possible for first segment
        % compare levels in the previous segment with AR threshold
        ARamp = ...
            smoothedARrms(prevSegStart:prevSegEnd)/ARthresholdPa;
        % all sub-treshold values are set to 1
        ARamp(ARamp<1)=1;
    else
        % first segment: no AR applicable
        ARamp=ones(size(inputThisSeg));
    end
    % attenuate input (NB cross product used)
    inputThisSeg=inputThisSeg./ARamp;
    
    % tympanic membrane response in meters
    stapesAR=inputThisSeg*stapesScalar;
    
    %                      DRNL
    filterCount=0;
    for cf=param.channelBFs
        filterCount=filterCount+1;
        %      get stapesAR for all channels
        y=stapesAR;        
        
        %-----------------------------------
        % nonlinear path - filter bws
        %-----------------------------------
        %Now defined in terms of octaves
        lowerCutOff=cf*2^(-param.bwOct(filterCount)/2);
        upperCutOff=cf*2^( param.bwOct(filterCount)/2);
        
        [y  NLg1bndry{filterCount}] = ...
            UTIL_bandPassFilter(y, param.filterOrder, lowerCutOff, upperCutOff, dt,NLg1bndry{filterCount});
        
        
        %                            NonLin compression function
        % first identify MOC attenuation to be applied (see below)
        %  this will be different for each channel
        MOC=MOCcontrol(filterCount,:);
        if all(MOC>0)
            DRNLa = DRNLaBaseline(filterCount)./MOC;
        else
            % prevention of denormals that may occur as a result of DIV/0 
            % - this should never happen as . . 
            % the MOC = x = 1+x*MOCfactor - where x is always +ve         
            DRNLa=zeros(size(MOC));
        end
        
        %                            broken stick compression function
        y = DRNL_brokenstick_nl (y, DRNLa, DRNLb(filterCount), c);
        
        [DRNLoutputThisSeg  NLg2bndry{filterCount}] = ...
            UTIL_bandPassFilter(y, param.filterOrder, lowerCutOff,upperCutOff, dt,NLg2bndry{filterCount});                        
        
        filteredSignal(filterCount, segmentStart:segmentEnd) = DRNLoutputThisSeg;
        
        % compute MOC control (smoothed squared DRNLoutput)
        MOCcutOff=1/(2*pi*MOCtau);
        [x  smoothedRMSbndry{filterCount}]=...
            UTIL_lowPassFilter(DRNLoutputThisSeg.^2,MOCcutOff, 1, dt, smoothedRMSbndry{filterCount});
        x=x.^0.5;    % restore to meaningful scale (meters)
        
        % compare MOC control signal to specified threshold
        x= x -param.MOCthreshold(filterCount);
        % and zero all values below threshold
        x(x<0)=0;
        
        % default attenuation is 1 ( i.e. no MOC gives unit gain)
        x= 1+ x *MOCfactor ;
        
        % save as control signal to attenuate output in next segment
        MOCcontrol(filterCount,:)= x;
        
    end  % BF channel
end  % segments

% OUTPUT=> combine across all channels
output=sum(filteredSignal,1);
hearingAidAmp = 10^(opScaling_dB/20);
output = output*hearingAidAmp;
