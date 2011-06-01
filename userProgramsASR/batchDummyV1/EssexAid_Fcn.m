% Function version of EssexAid_COPY2
% Nov 2010

function output = EssexAid_Fcn(inputSignal, sampleRate, param)
% close all; clear all; clc

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








%%
% global signaltype speechFileName sampleRate
% global toneFrequency toneDuration toneRampDuration
% global backgroundFileName noiseLevel backgroundRampDuration 
% global precedingBackgroundDuration trailingBackgrouneDuration
% global param
% 
% signaltype= 'tone';         % 'tone' or 'file'
% backgroundFileName='babble';% 'babble', 'PinkNoise2sec'
% 
% signalLevels=0:10:100;      % dB SPL
% noiseLevel= -80;            % dB SPL
% 
% toneDuration= .2;           % (s)
% toneRampDuration=0.02;      % (s)
% 
% precedingBackgroundDuration= 0.25;   % (s)
% trailingBackgrouneDuration= 0.05;    % (s)
% backgroundRampDuration=.005;         % onset and offset ramps for noise
% 
% 
% switch signaltype
%     case'tone'
%         % choose a frequency corresponding to one of the filters
%         toneFrequency=param.channelBFs(toneFrequencyIDX); 
%         speechFileName=['tone: ' num2str(toneFrequency) ' Hz'];
%     case 'speech'
%         speechFileName='twister_44kHz';
%         % speechFileName='1z67931a_44kHz';
%         % speechFileName='CallasShort';
%         % speechFileName='ieee01s01';
%         % speechFileName='13secspeech';
% end
% 
% % dbstop if error
% % figure(11), clf, drawnow
% % tic

% **  **  **  **  **  **  **  **  **  **  **  **   aid parameters
opScaling_dB=param.opScaling_dB;      % dB (BM displacement -> Pa), applied to output

%                             **stapes**
stapesScalar=	6e-008; % useful for calibrating the system

%                          **Acoustic reflex**
% This monitors recent sound levels and applies attenuation to
%   continuously pulls signal input level down to the AR threshold level
ARtau=param.ARtau;            % decay time constant
ARthresholddB=param.ARthresholddB;       % dB SPL (input signal level) =1 to disable
ARthresholdPa=  28e-6*10^(ARthresholddB/20);% Pa thresh for triggering AR

%                          **DRNL nonlinear path**
% useGammatone=0;         % 1=gammatone 0=butterworth
% filterOrder=3;          % used for both gammatone and Butterworth

% bandwidths (BW) are greater at higher frequencies
% NOW DEFINED IN THE PROCESSING LOOP IN OCTAVES
% p=.000246;   q=1.410546;
% nonlinBWs=param.channelBFs./(p*param.channelBFs+q); 
% this equation is ad hoc (better would be channelBF/p +q)

% DRNLbase is DRNL nonlin gain parameter, 'a' before MOC attenuation 
% DRNLaBaseline=param.DRNLaBaseline;      % DRNLa=MOCcontrol*DRNLaBaseline
DRNLaBaseline=param.DRNLaBaseline * 10.^(param.DRNLaBaselineRelative_dB/20)';

% DRNLb sets the level at which nonlin compression begins
%  normally about 20 dB above threshold
%  this parameter is best set visually using the tone I/O facility above
%  Disable the MOC function for clearer guidance
%DRNLb=param.DRNLb;            % to disable compression, set b to 1.
DRNLb=param.DRNLb * 10.^(param.DRNLbRelative_dB/20);%repmat(DRNLb,1, numChannels); % DRNLb is a vector: one per channel

c=param.DRNLc;                  % compression exponent (fixed)

%                       **DRNL linear path not used
% to activate this path uncomment the linear path code below
% LINbw=100;            % Hz fixed bandwidth (needs changing)
% LINgain=400;      
% LINgain=0;            % set to 0 to disable linear component
% 

%                              **MOC**
% Within-channel attenuation operates as AGC at relatively low levels
MOCtau=param.MOCtau;            % decay time constant
% MOCtau=0.4;             % decay time constant
MOCfactor=param.MOCfactor;          % arbitrary strength factor
efferentDelay= 0.010;   % feedback delay (also controls segment duration)

%                                                   **  signal generation
% outputSummary=[];   inputSummary=[]; % waveforms
% ARsummary=[];MOCsummary=[];          % control values
% inputLevelSummary=[]; outputLevelSummary=[]; % signalLevels
% 
% disp(' ')
% disp('start')
% disp(signaltype)
% disp(['numChannels: ' num2str(numChannels)])
% disp(['useGammatone: ' num2str(useGammatone)])
% % disp(['lin gain=' num2str(LINgain)])
% disp(['DRNLb= ' num2str(DRNLb(1))])
% disp(['  MOCtau= ' num2str(MOCtau)])
% disp(['  MOCfactor= ' num2str(MOCfactor,'%6.0e')])

maxLINEAR=0;
MAXnonlin=0;
% levelNumber=0; nLevels=length(signalLevels);
% for level=signalLevels
%     probeLevel=level;
%     SNR=probeLevel-noiseLevel;
%     levelNumber=levelNumber+1;
%     disp(' ')
%     disp(['level= ' num2str(level)])

%     [inputSignal, inputdBSPL]...
%         =aidStimulusPreparation(probeLevel);
    
    dt=1/sampleRate;
    duration=dt*length(inputSignal);
%     probeTime=dt:dt:duration;
    
    %                                     segment realated information
    % The signal is divided into 10 ms segments to allow for
    % feedback signals to be updated.
    % 10 ms is approximately the mandatory delay imposed 
    %  by delays in the nervous system
    %
    %  earObject length should be a multiple of  segment length
    %   (trim any excess signal)
    segmentDuration=efferentDelay;
%     [nEarObjectRows earObjectLength]=size(inputSignal);
    segmentLength=round(segmentDuration/dt);
    nSignalSegments=floor(duration/(segmentDuration));
    trimmedLength=nSignalSegments*segmentLength;
    inputSignal=inputSignal (1:trimmedLength);
%     [nEarObjectRows earObjectLength]=size(inputSignal);

    %                                             pre-allocation
    numChannels = numel(param.channelBFs);
    stapesBndry=[]; ARrmsBndry=[];
    NLg1bndry=cell(1,numChannels);
    NLg2bndry=cell(1,numChannels);
%     LINbndry=cell(1,numChannels);
    smoothedRMSbndry=cell(1,numChannels);
    stapesRecord=zeros(1, length(inputSignal));
    filteredSignal=zeros(numChannels, length(inputSignal));
    MOCreport=zeros(numChannels, length(inputSignal));
%     smoothedDRNL=zeros(numChannels, segmentLength);
    ARrecord=zeros(1, segmentLength);
%     MOCcontrol=repmat(DRNLaBaseline,numChannels,segmentLength);
    MOCcontrol=repmat(DRNLaBaseline,1,segmentLength);

    %% signal processing begins here
    % NB all units are international units and are kept close to 
    % realistic values for the purpose of comparing aid response to
    % BM response. This will not be necessesary in the hardware
    

    for segmentStart=1:segmentLength:length(inputSignal)
        segmentEnd=segmentStart+segmentLength-1;
        prevSegStart=segmentStart-segmentLength;
        prevSegEnd=segmentEnd-segmentLength;
%         nextSegmentStart=segmentStart+segmentLength;
%         nextSegmentEnd=segmentEnd+segmentLength;
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
        % record of AR control for monitoring purposes only
        ARrecord(segmentStart:segmentEnd)=  ARamp;

        % convert air velocity to displacement. Use lowpass filter
        [inputThisSeg stapesBndry]=...
            UTIL_lowPassFilter(inputThisSeg, 100, 1, dt, stapesBndry);
        
        % tympanic membrane response in meters
        stapesAR=inputThisSeg*stapesScalar;

        % record of stapesAR displacement for monitoring purposes only
        stapesRecord(segmentStart:segmentEnd)=stapesAR;

        %                                                             DRNL
        filterCount=0;
        for cf=param.channelBFs
            filterCount=filterCount+1;
            %      get stapesAR for all channels
            y=stapesAR;

            %                                               nonlinear path
%             if useGammatone
%                 %                            first gammatone
%                 [y  NLg1bndry{filterCount}] = ...
%                     UTIL_gammatone(y, filterOrder, cf,  ...
%                     nonlinBWs(filterCount), dt,NLg1bndry{filterCount});
%             else
                %-----------------------------------
                % nonlinear path - filter bws
                %-----------------------------------
                %Now defined in terms of octaves
                lowerCutOff=cf*2^(-param.bwOct(filterCount)/2);
                upperCutOff=cf*2^( param.bwOct(filterCount)/2);
%                 lowerCutOff=cf-nonlinBWs(filterCount)/2;
%                 upperCutOff=cf+nonlinBWs(filterCount)/2;
                [y  NLg1bndry{filterCount}] = ...
                    UTIL_bandPassFilter(y, param.filterOrder, ...
                    lowerCutOff,upperCutOff, dt,NLg1bndry{filterCount});
%             end

            %                            NonLin compression function
            % first identify MOC attenuation to be applied (see below)
            %  this will be different for each channel
            MOC=MOCcontrol(filterCount,:);
            if all(MOC>0)
                DRNLa = DRNLaBaseline(filterCount)./MOC;
            else
                % catch negative values of MOC
                DRNLa=zeros(size(MOC));
            end

            %                            broken stick compression function
            y = DRNL_brokenstick_nl (y, DRNLa, DRNLb(filterCount), c);

            %                            second gammatone
%             if useGammatone
%                 [nonlinOutput  NLg2bndry{filterCount}] = ...
%                     UTIL_gammatone(y, filterOrder, cf,  ...
%                     nonlinBWs(filterCount), dt,NLg2bndry{filterCount});
%             else
                [nonlinOutput  NLg2bndry{filterCount}] = ...
                    UTIL_bandPassFilter(y, param.filterOrder, ...
                    lowerCutOff,upperCutOff, dt,NLg2bndry{filterCount});
%             end


            % linear path. To activate, uncomment and set parameters above
            linOutput=zeros(size(nonlinOutput));
%             if LINgain>0
%                 y=stapesAR*LINgain;
%                 %  linearpath single  gammatone
%                 if useGammatone
%                     %                            first gammatone
%                     [linOutput  LINbndry{filterCount}] = ...
%                         UTIL_gammatone(y, filterOrder, cf,   ...
%                         LINbw, dt,LINbndry{filterCount});
%                 else
%                     [linOutput  LINbndry{filterCount}] = ...
%                         UTIL_bandPassFilter(y, filterOrder, ...
%                         cf-LINbw/2,cf+LINbw/2, dt,LINbndry{filterCount});
%                 end
%             end

            % combine linear and nonlinear path outputs
            DRNLoutputThisSeg=linOutput + nonlinOutput;

            % report peak output (across segments)
            maxLINEAR=max(maxLINEAR, max(linOutput));
            MAXnonlin=max(MAXnonlin,max(nonlinOutput));

            filteredSignal(filterCount, segmentStart:segmentEnd)=...
                DRNLoutputThisSeg;

            % compute MOC control (smoothed squared DRNLoutput)
            MOCcutOff=1/(2*pi*MOCtau);
            [x  smoothedRMSbndry{filterCount}]=...
                UTIL_lowPassFilter(DRNLoutputThisSeg.^2,MOCcutOff,...
                1, dt, smoothedRMSbndry{filterCount});
            x=x.^0.5;    % restore to meaningful scale (meters)
            
            % compare MOC control signal to specified threshold
            x= x -param.MOCthreshold;
            % and zero all values below threshold
            x(x<0)=0;
            
            % default attenuation is 1 ( i.e. no MOC gives unit gain)
            x= 1+ x *MOCfactor ;
            
            % save as control signal to attenuate output in next segment
            MOCcontrol(filterCount,:)= x;
            % report later the control signal
            MOCreport(filterCount, segmentStart:segmentEnd)=   ...
                MOCcontrol(filterCount,:);
        end  % BF channel
    end  % segments

    % monitor feedback control signals as dB
%     ARsummary=[ARsummary 20*log10(max(ARrecord))];
%     MOCsummary=[MOCsummary 20*log10(max(max(MOCreport)))];

%     disp(['maxLINEAR= ' num2str(20*log10(maxLINEAR/8e-11), '%5.0f')])
%     disp(['MAXnonlin= ' num2str(20*log10(MAXnonlin/1.7e-8), '%5.0f')])
%     disp(['maxAR= ' num2str(20*log10(max(ARrecord)))])
%     disp(['maxMOC= ' num2str(20*log10(max(max(MOCreport))))])

    % OUTPUT=> combine across all channels
    output=sum(filteredSignal,1);
    
    % filter to restore displacement to velocity (pressure)
    if param.doHpFilt
%         output =UTIL_highPassFilter(output, 10000, 1, dt, []);
        [b,a]=butter(1,10000/sampleRate*2,'high');
        output = filter(b,a,output);
    end
    hearingAidAmp = 10^(opScaling_dB/20);
    output = output*hearingAidAmp;

%     figure; plot(  20*log10(1./MOCreport')  ); legend(num2str(param.channelBFs')); ylabel('Atten (dB)')
%     disp(['inputdBSPL= ' num2str(inputdBSPL, '%5.0f')])  % found earlier
    
%     rms=(mean(output.^2)).^0.5;
%     rms=(mean(output(1.3e4:1.9e4).^2)).^0.5; %dont include silence in RMS
%     outputdBSPL=20*log10(rms/20e-6);
%     disp(['outputdBSPL= ' num2str(outputdBSPL, '%5.0f')])

    % collect data for I/O calculations
%     inputLevelSummary=[inputLevelSummary level ];
%     outputLevelSummary=[outputLevelSummary outputdBSPL ];
%     inputSummary=[inputSummary inputSignal];
%     outputSummary=[outputSummary output];
%     
%     % figure 1
%     if length(signalLevels)==1
%         aidStageDisplay(inputSignal,output,ARrecord,MOCreport,probeLevel)
%     end
%     toc
%     soundsc([inputSignal output],sampleRate)
% end  % levels

% % figure 11
% normalizedInput=inputSummary/(1.01*max(abs(inputSummary)));
% % time=dt:dt:dt*length(inputSummary);
% if length(signalLevels)>1
% level= min(signalLevels):(max(signalLevels)-min(signalLevels))/...
%     (length(normalizedInput)-1):max(signalLevels);
% else
%     level=dt:dt:dt*length(inputSummary);
% end
% figure(11), subplot(3,1,1)
% plot(level, normalizedInput, 'k')
% xlim([min(level) max(level)])
% title(['input: ' num2str(toneFrequency,'%6.0f') ' Hz.  '...
%     speechFileName '  output gain='  num2str(opScaling_dB) ' dB'])
% 
% normalizedOutput=outputSummary/(1.01*max(abs(outputSummary)));
% figure(11), subplot(3,1,2)
% set(gcf,'Position', [26   134   560   954])
% plot(level,normalizedOutput,'k')
% xlim([min(level) max(level)])
% title(['output: AR (tau thresholddB)=' ...
%     num2str([ARtau ARthresholddB],'%5.3f%4.0f') ...
%     ',  MOCtau (tau factor)=' num2str([MOCtau MOCfactor],'%5.3f %7g')])

% % save both input and output signals to file
% try
%     % the output level is arbitrary and the wavfile is normalised
%     wavwrite(normalizedInput, sampleRate, 'inputSummary')
%     wavwrite(normalizedOutput, sampleRate, 'outputSummary')
%     if playSounds
%         %         wavplay(normalizedInput, sampleRate)
%         switch signaltype
%             case 'file';
%                 % x = UTIL_lowPassFilter(normalizedOutput, 100, 1, dt);
%                 % wavplay(x/max(abs(x)), sampleRate)
%                 wavplay(normalizedOutput, sampleRate, 'async')
%         end
%     end
% catch
%     % error is caused by the file being open in Adobe audition
%     wavwrite(normalizedInput, sampleRate, 'reserveInput')
%     wavwrite(normalizedOutput, sampleRate, 'reserveOutput')
%     % all is not lost
%     [ding FS]=  wavread(['..\wavFileStore\ding']);
%     wavplay(ding, FS)
% end

% figure(11), subplot(3,1,3)
% plot(inputLevelSummary,outputLevelSummary,'r:', 'linewidth', 3)
% hold on
% plot([min(inputLevelSummary) max(inputLevelSummary)], ...
%     [ min(inputLevelSummary)  max(inputLevelSummary)], '-') % linear ref
% plot(signalLevels, ARsummary,'g')
% plot(signalLevels, MOCsummary,'k')
% hold off
% grid on
% title([ 'I/O:   DRNLa=' num2str(DRNLaBaseline', '%1g')...
%     '  DRNL.b=' num2str(DRNLb,' %6g')] )
% ylabel('dB (rms)')
% ylim([-10 100])
% xlabel( 'speech level (dB SPL)')
% if min(inputLevelSummary(:,1)) < max(inputLevelSummary(:,1))
%     xlim([min(inputLevelSummary(:,1)) max(inputLevelSummary(:,1))])
% end
% legend({'output','lin ref','AR','MOC'},'location','eastoutside')
% 
