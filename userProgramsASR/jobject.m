classdef jobject
    %JOBJECT Summary of this class goes here
    %   Detailed explanation goes here
    
    %%  *********************************************************
    %  properties                      _   _
    %                                 | | (_)
    %  _ __  _ __ ___  _ __   ___ _ __| |_ _  ___  ___
    % | '_ \| '__/ _ \| '_ \ / _ \ '__| __| |/ _ \/ __|
    % | |_) | | | (_) | |_) |  __/ |  | |_| |  __/\__ \
    % | .__/|_|  \___/| .__/ \___|_|   \__|_|\___||___/
    % | |             | |
    % |_|             |_|
    %************************************************************
    
    %% **********************************************************
    % Public properties - can be set by user
    %************************************************************
    properties(Access = public)
        wavFolder
        opFolder
        noiseFolder
        
        wavList
        todoStatus
        
        
        participant        = 'Normal';%'DEMO2_multiSpont';
        noiseName          = '8TalkerBabble';
        numWavs            =  5;
        noiseLevToUse      = -200;
        speechLevToUse     =  50;
        
        speechLevStd       = 0;
        noiseLevStd        = 0;
        freezeNoise        = 0;
        meanSNR            = 20;
        speechDist         = 'None';
        noiseDist          = 'None';
        
        noisePreDur        = 0.0;
        noisePostDur       = 0.0;
        truncateDur        = 0.0; %Dr. RF used 0.550
        
        currentSpeechLevel
        currentNoiseLevel
        
        
        useSpectrogram = false;
        numCoeff = 9;
        
        %************************************************************
        % plotting handles
        %************************************************************
        probHaxes   = [];
        probHaxesSM = [];
        sacfHaxes   = [];
        sacfHaxesSM = [];
        featHaxes   = [];
        reconHaxes  = [];
        
        %************************************************************
        % SACF params
        %************************************************************
        useSACF             = false;
        
        
        SACFacfTau          = 2; % > 1 = Wiegrebe mode
        SACFnBins           = 128;
        SACFminLag          = 1 / 4000;
        SACFmaxLag          = 1 / 50;
        SACFlambda          = 10e-3;
        
        %************************************************************
        % MAP params
        %************************************************************
        MAProot                 = fullfile('..');
%         mapVer                  = 'MAP1_12_experimental';
%         MAPuseEfferent          = 1;
        MAPplotGraphs           = 0;
        MAPDRNLSave             = 0;
        
        MAPopLSR                = 0;
        MAPopMSR                = 0;
        MAPopHSR                = 1;
        
        MAPparamChanges = {};
        
        %************************************************************
        % HTK stuff - writing to HTK recogniser format
        %************************************************************
        frameshift              = 10;       % shift between frames (ms)
        sampPeriodFromMsFactor  = 10000;    % appropriate for 10ms frame rate
        paramKind               = 9;        % HTK USER format for user-defined features (p73 of HTK book)
        
        removeEnergyStatic      = false;
        
        %************************************************************
        % Portable EssexAid params
        %************************************************************        
%         sampleRate = 48e3;        
        bwOct = 1/1;
        filterOrder  = 2;  
                
        mainGain = [ 1;    1;    1;    1;    1];     % gain in linear units
        TCdBO    = [40;   40;   40;   40;   40];      %Compression thresholds (in dB OUTPUT from 2nd filt)
        TMdBO    = [10;   10;   10;   10;   10];      %MOC thresholds (in dB OUTPUT from 2nd filt)
        DRNLc    = [ 0.2;  0.2;  0.2;  0.2;  0.2;]
        
        ARtau  = 0.03;            % decay time constant
        ARthresholddB = 85;       % dB SPL (input signal level) =>200 to disable        
        MOCtau = 0.3;
        MOCfactor = 0.5;   %dB per dB OUTPUT
                      
        numSamples = 1024; %MAX=6912        
        useAid = 0;
    end
        
    %%  *********************************************************
    % Protected properties - inter/ra class communication
    %************************************************************
    properties(Access = protected)
        jobLockFid               
        
        %Nick C comment on this:
        %OK. The big-endian thing works because in the config file
        %'config_tr_zcpa12' there is a flag called NATURALREADORDER that is set to
        %FALSE and thus appears to override x86 standard
        %little-endian. Endianess has **!nothing!** to do with win vs *nix
        byteOrder = 'le';  % byte order is big endian
    end
    properties(Dependent = true)
        jobLockTxtFile
    end
    
    %%  *********************************************************
    % methods        _   _               _
    %               | | | |             | |
    % _ __ ___   ___| |_| |__   ___   __| |___
    %| '_ ` _ \ / _ \ __| '_ \ / _ \ / _` / __|
    %| | | | | |  __/ |_| | | | (_) | (_| \__ \
    %|_| |_| |_|\___|\__|_| |_|\___/ \__,_|___/
    %************************************************************
    
    methods
        %% **********************************************************
        % Constructor
        %************************************************************
        function obj = jobject(LearnOrRecWavs, jobFolder)
            if nargin < 1
                warning('job:ambiguouscorpus',...
                    'We need to know whether the (L)earning or (R)ecognition corpus is being used, assuming ''L''');
                LearnOrRecWavs = 'L';
            end
            if nargin > 1
                obj.opFolder = jobFolder;
            else
                if isunix
                    if ismac
                        obj.opFolder = '~/ASR/exps/_foo';
                    else
                        obj.opFolder = '/scratch/nrclark/exps/_foo';
                    end
                else
                    obj.opFolder = 'D:\exps\_foo';
                end                    
            end

            obj = obj.assignWavPaths(LearnOrRecWavs);
%             mkdir(obj.opFolder);
            obj = obj.initMAP;
            
            
        end % ------ OF CONSTRUCTOR
        
        %% **********************************************************
        % Set Wav Paths
        %************************************************************
        function obj = assignWavPaths(obj, LearnOrRecWavs)
            if isunix
                if ismac
                    lWAVpath = '~/ASR/reducedAURORA/TrainingData-Clean/';
                    rWAVpath  = '~/ASR/reducedAURORA/TripletTestData/';
                    obj.noiseFolder = '~/ASR/noises';
                else                    
                    lWAVpath = '/scratch/nrclark/corpora/AURORA digits (wav)/TrainingData-Clean/';
                    rWAVpath  = '/scratch/nrclark/corpora/AURORA digits (wav)/TripletTestData/';
                    obj.noiseFolder = '/scratch/nrclark/corpora/noises';
                end
            else
                %lWAVpath = 'D:\ASRexperiments\Stimuli\AURORA digits (wav)\TrainingData-Clean';
                %rWAVpath  = 'D:\ASRexperiments\Stimuli\AURORA digits (wav)\TripletTestData';
                lWAVpath = 'C:\corpora\AURORA digits (wav)\TrainingData-Clean';
                rWAVpath = 'C:\corpora\AURORA digits (wav)\TripletTestData';
                %obj.noiseFolder = 'D:\ASRexperiments\Stimuli\noises';
                obj.noiseFolder = 'C:\corpora\noises';
            end
            
            if strcmpi(LearnOrRecWavs, 'l')
                obj.wavFolder = lWAVpath;
            elseif strcmpi(LearnOrRecWavs, 'r')
                obj.wavFolder = rWAVpath;
            else
                error('First argument to constructor must be ''L'' or ''R''');
            end
        end
        
        %% **********************************************************
        % mutex related    _
        %                 | |
        %  _ __ ___  _   _| |_ _____  __
        % | '_ ` _ \| | | | __/ _ \ \/ /
        % | | | | | | |_| | ||  __/>  <
        % |_| |_| |_|\__,_|\__\___/_/\_\
        %************************************************************
        
        %% **********************************************************
        % lockJobList - File mutex protecting from race conditions
        %************************************************************
        function obj = lockJobList(obj)
            lockON = false;
            while (~lockON)
                if numel(dir(obj.jobLockTxtFile)) %Check to see if lock already in place
                    wTime = randi(750)+250; %3,000-10,000 ms
                    disp(['File mutex in place. Retrying in ' num2str(wTime) ' ms'])
                    pause(wTime/1000);
                else
                    obj.jobLockFid = fopen(obj.jobLockTxtFile,'w');
                    disp('Job locked');
                    pause(50/1000);
                    lockON = true;
                end
            end
            fclose(obj.jobLockFid);
        end% ------ OF LOCKJOBLIST
        
        %% **********************************************************
        % unlockJobList - unlocks for other workers
        %************************************************************
        function obj = unlockJobList(obj)
            lockON = logical(numel(dir(obj.jobLockTxtFile)));
            while(lockON)
                try
                    delete(obj.jobLockTxtFile);
                    disp('Job unlocked');
                    pause(50/1000);
                    lockON = false;
                catch %#ok<CTCH>
                    disp('Unjamming lock - retrying immediately')
                    pause(200/1000)
                end
            end
        end% ------ OF UNLOCKJOBLIST
        
        %% **********************************************************
        % storeSelf - save a copy of this object in opFolder
        %************************************************************
        function storeSelf(obj)
            doneSaving = false;
            while(~doneSaving)
                try
                    save(fullfile(obj.opFolder, 'jobObject'), 'obj');
                    doneSaving = true;
                catch  %#ok<CTCH>
                    wTime = randi(3800)+200; %200-4000 ms
                    disp(['Save collision (THIS IS VERY VERY BAD - have you not used the mutex?). Retrying in ' num2str(wTime) ' ms']);
                    pause(wTime/1000);
                end
            end
        end % ------ OF STORESELF
        
        %% **********************************************************
        % loadSelf - load a copy of this object from opFolder
        %************************************************************
        function obj = loadSelf(obj)
            doneLoading = false;
            while(~doneLoading)
                try
                    load(fullfile(obj.opFolder,'jobObject'));
                    doneLoading = true;
                catch  %#ok<CTCH>
                    wTime = randi(3800)+200; %200-4000 ms
                    disp(['Load collision (THIS IS VERY VERY BAD - have you not used the mutex?). Retrying in ' num2str(wTime) ' ms'])
                    pause(wTime/1000);
                end
            end
        end% ------ OF LOADSELF
        
        
        %%  *********************************************************
        %  _                          _                   _
        % | |                        | |                 (_)
        % | |__   ___  _   _ ___  ___| | _____  ___ _ __  _ _ __   __ _
        % | '_ \ / _ \| | | / __|/ _ \ |/ / _ \/ _ \ '_ \| | '_ \ / _` |
        % | | | | (_) | |_| \__ \  __/   <  __/  __/ |_) | | | | | (_| |
        % |_| |_|\___/ \__,_|___/\___|_|\_\___|\___| .__/|_|_| |_|\__, |
        %                                          | |             __/ |
        %                                          |_|            |___/
        %************************************************************
        
        %% **********************************************************
        % get method for dependtent var jobLockTxtFile
        %************************************************************        
        function value = get.jobLockTxtFile(obj)
            %Name the MUTEX file here
            value = [fullfile(obj.opFolder, 'jobLock') '.txt'];
        end
        
        %% **********************************************************
        % checkStatus - see how much of the current job is complete
        %************************************************************
        function checkStatus(obj)
            NOWopen = sum(obj.todoStatus==0); %Starting from R2010b, Matlab supports enumerations. For now we use horrible integers for compatibility.
            NOWpend = sum(obj.todoStatus==1);
            NOWdone = sum(obj.todoStatus==2);
            
            zz = clock;
            disp([num2str(zz(3)) '-' num2str(zz(2)) '-' num2str(zz(1)) '   ' num2str(zz(4)) ':' num2str(zz(5))])
            disp(['CURRENT JOB:' obj.opFolder]);
            disp(' ')
            disp(['open - ' num2str(NOWopen) ' || pending - ' num2str(NOWpend) ' || complete - ' num2str(NOWdone)])
            
            % ---- PROGRESS BAR ----
            pcDone = 100*NOWdone/obj.numWavs;
            progBarLength = 40;
            charBars = repmat('=',1,floor(pcDone/100 * progBarLength));
            charWhiteSpace = repmat(' ',1, progBarLength - numel(charBars));
            disp(' ')
            disp([' -[' charBars charWhiteSpace ']-  ' num2str(pcDone, '%0.1f') '%'])
            disp(' ')
            disp(' ')
            % -- END PROGRESS BAR ---
        end% ------ OF CHECKSTATUS
        
        %% **********************************************************
        % initMAP - add MAP stuff to path
        %************************************************************
        function obj = initMAP(obj)
            addpath(...fullfile(obj.MAProot, 'modules'),...
                fullfile(obj.MAProot, 'utilities'),...
                fullfile(obj.MAProot, 'MAP'),...
                fullfile(obj.MAProot, 'parameterStore'),...
                fullfile('ASR files'));
        end % ------ OF INIT MAP
        
        %% **********************************************************
        % assign files to testing and training sets
        %************************************************************
        function obj = assignFiles(obj)
            speechWavs  = dir(fullfile(obj.wavFolder, '*.wav'));
            assert(obj.numWavs <= size(speechWavs, 1) ,...
                'not enough files available in the corpus');  % make sure we have enough wavs
            
            randomWavs = rand(1, size(speechWavs, 1));
            [~,b] = sort(randomWavs);
            trainFileIdx = b(1:obj.numWavs);
            
            obj.wavList  = speechWavs(trainFileIdx); %This is a record of all of the wavs that should be done
            
            %Starting from R2010b, Matlab supports enumerations. For now we
            %use horrible integers for compatibility.
            %0=open, 1=processing, 2=done
            obj.todoStatus = zeros(numel(obj.wavList), 1); 
            
            
        end % ------ OF ASSIGN FILES
        
        %% **********************************************************
        % generate  feature
        %************************************************************
        function obj = genFeat(obj, currentWav)                        
            fprintf(1,'Processing: %s \n', currentWav);
            if strcmpi(obj.speechDist,'Gaussian')
                tempSpeechLev = obj.speechLevToUse + obj.speechLevStd*randn;
            elseif strcmpi(obj.speechDist,'Uniform')
                % for a uniform distribution, the standard deviation is
                % range/sqrt(12)
                % http://mathforum.org/library/drmath/view/52066.html
                tempSpeechLev = obj.speechLevToUse - obj.speechLevStd*sqrt(12)/2 + obj.speechLevStd*sqrt(12)*rand;
            elseif strcmpi(obj.speechDist,'None')
                tempSpeechLev = obj.speechLevToUse;
            end
            
            if strcmpi(obj.noiseDist,'Gaussian')
                tempNoiseLev  = speechLev - obj.meanSNR  + obj.noiseLevStd*randn;
            elseif strcmpi(obj.noiseDist,'Uniform')
                % for a uniform distribution, the standard deviation is
                % range/sqrt(12)
                % http://mathforum.org/library/drmath/view/52066.html
                tempNoiseLev  = tempSpeechLev - obj.meanSNR - obj.noiseLevStd*sqrt(12)/2 + obj.noiseLevStd*sqrt(12)*rand;
            elseif strcmpi(obj.noiseDist,'None')
                tempNoiseLev  = obj.noiseLevToUse;
            end
            
            disp(['Current speech level = ' num2str(tempSpeechLev)]);
            disp(['Current noise level  = ' num2str(tempNoiseLev)]);
            
            obj.currentSpeechLevel = tempSpeechLev;
            obj.currentNoiseLevel = tempNoiseLev;
            %             obj.currentFeatureFile = currentWav;
            [finalFeatures, ~] = processWavs(obj, currentWav); %discard the output from ANprobabilityResponse and method using ~                                    
            opForHTK(obj, currentWav, finalFeatures);                        
        end % ------ OF GENFEAT
        
        %% **********************************************************
        % write o/p in HTK friendly format
        %************************************************************
        function obj = opForHTK(obj, currentWav, featureData)
            
            featureName = strrep(currentWav, '.wav','.map');
            targetFilename = fullfile(obj.opFolder, featureName);
            
            % write in a format HTK compliant for the recogniser to use
            writeHTK(...
                targetFilename,...
                featureData,...
                size(featureData,2),...
                obj.frameshift*obj.sampPeriodFromMsFactor,...
                size(featureData,1)*4,...
                obj.paramKind,...
                obj.byteOrder);
            
            % make lists of training data/files - SAVE FOR HMM STAGE
            %             createSCP(obj, targetPath); % filename list (training)
            %             createMLF(obj, targetPath); %, 0); % training data (without short pauses)
        end % ------ OF opForHTK
        
        
        %% **********************************************************
        %      _                   _                                     _
        %     (_)                 | |                                   (_)
        %  ___ _  __ _ _ __   __ _| |  _ __  _ __ ___   ___ ___  ___ ___ _ _ __   __ _
        % / __| |/ _` | '_ \ / _` | | | '_ \| '__/ _ \ / __/ _ \/ __/ __| | '_ \ / _` |
        % \__ \ | (_| | | | | (_| | | | |_) | | | (_) | (_|  __/\__ \__ \ | | | | (_| |
        % |___/_|\__, |_| |_|\__,_|_| | .__/|_|  \___/ \___\___||___/___/_|_| |_|\__, |
        %         __/ |               | |                                         __/ |
        %        |___/                |_|                                        |___/
        %************************************************************
        
        %% **********************************************************
        % getStimulus - what it says on the tin
        %************************************************************
        function [stimulus, sampleRate] = getStimulus(obj, currentWav)
            
            % NRCgetStimulus.m - N.R.Clark Aug 2010
            % Modified version of RTF script to include:
            %
            % 1)Signal and noise samples with different sample rate. The component with
            % lower sample rate is upsampled to the rate of that with the higher rate.
            %
            % 2) Neater level setting
            %
            % 3) Parameter to change noise intro duration
            %
            % 4) Added some silence on the tail (read babble) same duration as noise
            % pre dur
            %
            % ORIGINAL HEADER:
            % getStimulus.m
            %
            % Robert T. Ferry
            % 13th May 2009
            
            
            % Set levels
            [speech speechSampleRate] = wavread(fullfile(obj.wavFolder, currentWav ));
            speech = speech./sqrt(mean(speech.^2)); %Normalize RMS to 1
            speech = speech * 20e-6 * 10^(obj.currentSpeechLevel/20); %Convert RMS to pascals at desired level
            %disp(20*log10(sqrt(mean(speech.^2))/20e-6))
            
            [noise noiseSampleRate] = wavread(fullfile(obj.noiseFolder, obj.noiseName ));
            noise = noise./sqrt(mean(noise.^2)); %Normalize RMS to 1
            noise = noise * 20e-6*10^(obj.currentNoiseLevel/20); %Convert RMS to pascals at desired level
            %disp(20*log10(sqrt(mean(noise.^2))/20e-6))
            
            
            % Do sample rate conversion if needed
            % Will always convert stimulus component with lower rate up to that with
            % higher rate.
            if speechSampleRate > noiseSampleRate
                %     disp('S>N')
                [p,q] = rat(speechSampleRate/noiseSampleRate,0.0001);
                noise = resample(noise, p, q);
                noiseSampleRate = speechSampleRate;
            elseif noiseSampleRate > speechSampleRate
                %     disp('N>S')
                [p,q] = rat(noiseSampleRate/speechSampleRate,0.0001);
                speech = resample(speech, p, q);
                speechSampleRate = noiseSampleRate; %#ok<NASGU>
            else
                %DO NOTHING BUT ASSERT
                assert(speechSampleRate == noiseSampleRate);
            end
            %disp(20*log10(sqrt(mean(speech.^2))/20e-6))
            sampleRate = noiseSampleRate;
            dt = 1/sampleRate;
            
            % mix stimuli
            % Everything from here down (mostly) is RTF's original
            % add 1/2  second prior to stimulus
            % and 0 s following it
            silenceStart = floor(obj.noisePreDur*sampleRate);
            silenceEnd = floor(obj.noisePostDur*sampleRate);
            
            silencePointsStart = zeros(silenceStart,1);
            silencePointsEnd = zeros(silenceEnd,1);
            
            speech = [silencePointsStart; speech; silencePointsEnd];
            
            stimLength = length(speech);
            noiseLength = length(noise);
            if obj.freezeNoise
                idx = 1;
            else
                idx = ceil(rand*(noiseLength-stimLength));
            end
            noise = noise(idx:idx+stimLength-1);
            
            stimulus = speech+noise;
            %disp(20*log10(sqrt(mean(stimulus.^2))/20e-6))
            
            % add ramps to noise
            stimInNoiseTime = dt:dt:dt*length(stimulus);
            rampDuration = 0.100;
            rampTime = dt : dt : rampDuration;
            ramp = [0.5*(1+cos(2*pi*rampTime/(2*rampDuration)+pi)) ones(1,length(stimInNoiseTime)-length(rampTime))];
            stimulus = stimulus'.*ramp;
            ramp = fliplr(ramp);
            stimulus = stimulus.*ramp;
            stimulus = stimulus';
            %disp(20*log10(sqrt(mean(stimulus.^2))/20e-6))
            
            % add additional silent period to start of stimulus for model to 'settle down'
            additionalSilenceLength = round(0.050*sampleRate);
            additionalSilencePoints = zeros(additionalSilenceLength,1);
            stimulus = [additionalSilencePoints; stimulus]'; %now rotated.
            %disp(20*log10(sqrt(mean(stimulus.^2))/20e-6))
        end% ------ OF GETSTIMULUS
        
        
        
        %% **********************************************************
        % processWavs - do all of the MAP signal processing
        %************************************************************
        function [finalFeatures, ANprobabilityResponse] = processWavs(obj, currentWav)
            
            %**********************************************************
            % FIRST GET THE STIMULUS
            %**********************************************************
            [stimulus, sampleRate] = obj.getStimulus(currentWav);
%             disp(20*log10(sqrt(mean(stimulus.^2))/20e-6))

            
            %**********************************************************
            % NOW TO LOAD IN THE HEARING AID
            %**********************************************************
            if obj.useAid
%                 disp(size(stimulus))
                stimulus = [stimulus; stimulus]'; %EsxAid requires stereo stim
%                 disp(size(stimulus))
                stimulus = EssexAid_guiEmulatorWrapper(stimulus, sampleRate, obj);
                stimulus = stimulus(1,:); %convert back to mono
%                 disp(size(stimulus))
            end
            
%             %**********************************************************
%             % Map related bits
%             %**********************************************************
%             %             % global variables - necessary for changing model pars
%             global AN_IHCsynapseParams inputStimulusParams IHCpreSynapseParams DRNLParams; %#ok<NUSED>
%             
%             paramsToLoad = ['method=MAPparams', obj.participant, '();'];
%             eval(paramsToLoad);
%             %             DRNLParams.fixedMOCAttenuation= -20;
%             
%             % now that the correct parameters have been loaded, change
%             % them!
%             inputStimulusParams.sampleRate = sampleRate;
%OBSOLETE              method.plotGraphs =	obj.MAPplotGraphs;
%             method.DRNLSave = obj.MAPDRNLSave;
            
            %**********************************************************
            % If using efferent change some params pre processing
            %**********************************************************
%             if obj.MAPuseEfferent
%                 method.saveANprobability = 1; %new
%                 method.useEfferent = 1;
%                 MAPmoduleSequence       = [1 2 3 4 5 6 7 8 9];
%                 AN_IHCsynapseParams.mode= 'spikes'; %now define here and hide from class interface (the global)
%             else
%                 method.saveANprobability = 0; %new
%                 method.useEfferent = 0;
%                 MAPmoduleSequence       = [1 2 3 4 5 6]; %= MAP 1_13, [1 2 3 4 5 6 7];= MAP1_12
%                 AN_IHCsynapseParams.mode= 'probability';
%                 
%                 %Below is some stuff related to HSR / LSR
%                 %This could be done lower down but it makes little
%                 %difference. Probably more efficient to just not calculate
%                 %it for probability mode
%                 
%                 if numel(IHCpreSynapseParams.tauCa)==2
%                     %                     LSRtauCa=IHCpreSynapseParams.tauCa(1);
%                     HSRtauCa=IHCpreSynapseParams.tauCa(end);
%                 elseif numel(IHCpreSynapseParams.tauCa)==3
%                     LSRtauCa=IHCpreSynapseParams.tauCa(1);
%                     MSRtauCa=IHCpreSynapseParams.tauCa(2);
%                     HSRtauCa=IHCpreSynapseParams.tauCa(end);
%                 else
%                     LSRtauCa=[];
%                     HSRtauCa=IHCpreSynapseParams.tauCa; %incase the model used only has one fibre type
%                 end
%                 
%                 IHCpreSynapseParams.tauCa= []; %wipe it
%                 if obj.MAPopLSR
%                     assert(~isempty(LSRtauCa), 'No LSR data available for this parameter file')
%                     IHCpreSynapseParams.tauCa = LSRtauCa;
%                 end
%                 if obj.MAPopMSR
%                     assert(~isempty(MSRtauCa), 'No MSR data available for this parameter file')
%                     IHCpreSynapseParams.tauCa = [IHCpreSynapseParams.tauCa MSRtauCa];
%                 end
%                 if obj.MAPopHSR
%                     IHCpreSynapseParams.tauCa = [IHCpreSynapseParams.tauCa HSRtauCa];
%                 end
%             end
%             
%             %             disp(AN_IHCsynapseParams);
%             % run the model
%             %             method.saveANprobability = 1;
%             %             method.returnSynapseContents = 0;
%             %             method.saveSynapseContents = 0;
% %             [ANprobabilityResponse, method, ~] = MAPsequenceSeg(stimulus, method, MAPmoduleSequence);

%             paramChanges = {};
%             if ~obj.MAPuseEfferent
%                 %disp('efferent off')
%                 paramChanges{numel(paramChanges)+1}= 'OMEParams.rateToAttenuationFactor=0;';        % 0= off
%                 paramChanges{numel(paramChanges)+1}= 'OMEParams.rateToAttenuationFactorProb=0;';    % 0= off
%                 paramChanges{numel(paramChanges)+1}= 'DRNLParams.rateToAttenuationFactor = 0;';     % 0 = MOC off (probability)
%                 paramChanges{numel(paramChanges)+1}= 'DRNLParams.rateToAttenuationFactorProb = 0;'; % 0 = MOC off (spikes)
%             end
            AN_spikesOrProbability = 'probability';
            
            if obj.useSpectrogram
                lowestBF=100; 	highestBF= 4500; 	numChannels=30;
                F=round(logspace(log10(lowestBF),log10(highestBF),numChannels));
                
                nfft = 1024;
                hopSamples = 64;
                noverlap = nfft - hopSamples;
                dt = hopSamples/sampleRate;
                %method.dt = ANdt;
                
                [~,~,~,P] = spectrogram(stimulus,nfft,noverlap,F,sampleRate);    
                
                ANprobabilityResponse = 10*log10(  abs(P) /  ((20e-6)^2)  ); %now correct [(a^2)/(b^2) = (a/b)^2]
                
                
            else
                [ANprobabilityResponse, dt, myBFlist] = MAPwrap(stimulus, sampleRate, -1, obj.participant, AN_spikesOrProbability, obj.MAPparamChanges);
            end
            nChannels = numel(myBFlist);
            
            
            
           
            
            %**********************************************************
            % If using efferent do some different post processing
            %**********************************************************
            
                
%             if obj.MAPuseEfferent
%                 ANprobabilityResponse = method.ANprobabilities;
%             else
%                 %do nothing - in the conditional statement a few lines
%                 %above, the output mode is set to probabilities if the
%                 %efferent system is switched off. The probabilities do not
%                 %need to be extracted from something stored.
%             end
            
            % remove the additional silence periods from start
            time_ANresponse = dt:dt:dt*size(ANprobabilityResponse,2);
            idx = time_ANresponse > obj.truncateDur; %RTF had this @ 0.550
            ANprobabilityResponse = ANprobabilityResponse(:, idx);
            
             
            
%             %**********************************************************
%             % Quick method to see whether we have both LSR and HSR fibres
%             %**********************************************************
%             if obj.MAPuseEfferent %only test for MacGrates if we've calulated them (module #9)
%                 if numel(method.MacGratesAllCells) == 2
%                     nChannels = size(ANprobabilityResponse,1) / 2;
%                 elseif numel(method.MacGratesAllCells) == 3
%                     nChannels = size(ANprobabilityResponse,1) / 3;
%                 else
%                     nChannels = size(ANprobabilityResponse,1);
%                 end
%             else
%                 nChannels = size(ANprobabilityResponse,1)/sum([obj.MAPopHSR obj.MAPopMSR obj.MAPopLSR]);
%             end
            
%             % The following code is dependent on wheter efferent is used or
%             % not. When efferent is used then both HSR and LSR responses
%             % are generated bacause the AR efferent process is driven by
%             % Low spont o/p and MOC process is drivven by high spont o/p.
%             % When the effernt process is disabled then only the spont
%             % rates desired need to be calculated.
%             if obj.MAPuseEfferent
%                 if obj.MAPopLSR && ~obj.MAPopHSR && ~obj.MAPopMSR
%                     ANprobabilityResponse = ANprobabilityResponse(1:nChannels, :); %use LSR
%                 elseif ~obj.MAPopLSR && obj.MAPopHSR && ~obj.MAPopMSR
%                     ANprobabilityResponse = ANprobabilityResponse(end-nChannels+1:end, :); %or use HSR
%                 elseif ~obj.MAPopLSR && ~obj.MAPopHSR && obj.MAPopMSR
%                     ANprobabilityResponse = ANprobabilityResponse(nChannels+1:end-nChannels, :); %or use MSR
%                 else %User hasn't specified either - assertion check anyway
%                     assert(sum([obj.MAPopLSR obj.MAPopHSR obj.MAPopMSR] ) == 1, 'You nead output from ONLY one fibre type!!');
%                 end
%             else
%                 if sum([obj.MAPopLSR obj.MAPopHSR obj.MAPopMSR] ) == 1 %only one or the other is generated in this instance
%                     ANprobabilityResponse = ANprobabilityResponse(:, :); %use LSR or HSR
%                 else %User hasn't specified either - assertion check anyway
%                     assert(sum([obj.MAPopLSR obj.MAPopHSR obj.MAPopMSR] ) == 1, 'You nead output from ONLY one fibre type!!');
%                 end
%             end
            
%             nChannels = numel(method.nonlinCF); %nonlin-cf is one of the few remaining method fields in 1_14
            if ~obj.useSpectrogram
                if obj.MAPopLSR && ~obj.MAPopHSR
                    ANprobabilityResponse = ANprobabilityResponse(1:nChannels, :); %use LSR
                end
                if ~obj.MAPopLSR && obj.MAPopHSR
                    ANprobabilityResponse = ANprobabilityResponse(end-nChannels+1:end, :); %or use HSR
                end
                if obj.MAPopMSR
                    assert(0,'Not working with MSR at the mo')
                end
            end
            
            % OPTIONAL PLOTTING
            YTickIdx = 1:floor(numel(myBFlist)/6):numel(myBFlist);
            YTickIdxRev = numel(myBFlist)+1-YTickIdx;
            if ~isempty(obj.probHaxes)
                axes(obj.probHaxes);  %#ok<MAXES>
                imagesc(ANprobabilityResponse); colorbar('peer', obj.probHaxes)                
                set(obj.probHaxes, 'YTick', YTickIdx);                
                set(obj.probHaxes, 'YTickLabel', num2str(    myBFlist(YTickIdxRev)'     ));
                ylabel('cf in Hz')                
            end
            
            % OPTIONAL PLOTTING SMOOTHED
            if ~isempty(obj.probHaxesSM)
                axes(obj.probHaxesSM); %#ok<MAXES>
                imagesc([],1:size(ANprobabilityResponse,1),flipud(obj.makeANsmooth(ANprobabilityResponse, 1/dt)));
                set(obj.probHaxesSM, 'YTick', YTickIdx);
                set(obj.probHaxesSM, 'YTickLabel', num2str(    myBFlist(YTickIdxRev)'     ));
                shading(obj.probHaxesSM, 'flat'); colorbar('peer', obj.probHaxesSM)
                ylabel('cf in Hz')
            end
         
            
            %**********************************************************
            % optional SACF stage
            %**********************************************************
            if obj.useSACF
                
                % A slightly ugly copying is needed 
%                 SACFparams.minLag = obj.SACFminLag;
%                 SACFparams.maxLag = obj.SACFmaxLag;       
                SACFparams.lambda = obj.SACFlambda;
                SACFparams.acfTau = obj.SACFacfTau;                
                SACFparams.lags = logspace(log10(obj.SACFminLag),log10(obj.SACFmaxLag),obj.SACFnBins);
                SACFparams.lags = linspace(obj.SACFminLag, obj.SACFmaxLag,obj.SACFnBins );
                
                SACFmethod.dt = dt;
                SACFmethod.nonlinCF = myBFlist;
                
                %This is slightly misleading as the ANprob is now a SACF
                [ANprobabilityResponse, ~, ~, ~] = filteredSACF(ANprobabilityResponse, SACFmethod, SACFparams);
                
                % OPTIONAL PLOTTING
                YTickIdx = 1:floor(obj.SACFnBins/6):obj.SACFnBins;
                %YTickIdxRev = obj.SACFnBins+1-YTickIdx;
                if ~isempty(obj.sacfHaxes)
                    axes(obj.sacfHaxes);  %#ok<MAXES>
                    imagesc(flipud(ANprobabilityResponse)); shading(obj.sacfHaxes, 'flat'); colorbar('peer', obj.sacfHaxes)
                    set(obj.sacfHaxes, 'YTick', YTickIdx);
                    set(obj.sacfHaxes, 'YTickLabel', num2str(    1./SACFparams.lags(YTickIdx)'    ,'%0.1f' ));
                    ylabel('Pitch in Hz')                    
                end
                
                % OPTIONAL PLOTTING SMOOTHED
                if ~isempty(obj.sacfHaxesSM)
                    axes(obj.sacfHaxesSM);  %#ok<MAXES>
                    imagesc(flipud(obj.makeANsmooth(ANprobabilityResponse, 1/dt))); shading(obj.sacfHaxesSM, 'flat'); colorbar('peer', obj.sacfHaxesSM)
                    set(obj.sacfHaxesSM, 'YTick', YTickIdx);
                    set(obj.sacfHaxesSM, 'YTickLabel', num2str(    1./SACFparams.lags(YTickIdx)'    ,'%0.1f' ));
                    ylabel('Pitch in Hz') 
                end

                
            end
            
                
            finalFeatures = obj.makeANfeatures(  ...
                obj.makeANsmooth(ANprobabilityResponse, 1/dt), obj.numCoeff  );            
            
            if obj.removeEnergyStatic
                finalFeatures = finalFeatures(2:end,:);
                % disp(size(finalFeatures))
            end
            
            % OPTIONAL PLOTTING
            if ~isempty(obj.featHaxes)
                pcolor(obj.featHaxes, finalFeatures); shading(obj.featHaxes, 'flat'); colorbar('peer', obj.featHaxes)
            end
            if ~isempty(obj.reconHaxes)
                 reconsData = idct(finalFeatures,obj.SACFnBins);
                 axes(obj.reconHaxes);  %#ok<MAXES>
                 imagesc(flipud( reconsData ));
            end
            
            opForHTK(obj, currentWav, finalFeatures);
        end % ------ OF PROCESSWAVS        
        
    end % ------ OF METHODS
    
    %%  *********************************************************
    %      _        _   _                       _   _               _
    %     | |      | | (_)                     | | | |             | |
    %  ___| |_ __ _| |_ _  ___   _ __ ___   ___| |_| |__   ___   __| |___
    % / __| __/ _` | __| |/ __| | '_ ` _ \ / _ \ __| '_ \ / _ \ / _` / __|
    % \__ \ || (_| | |_| | (__  | | | | | |  __/ |_| | | | (_) | (_| \__ \
    % |___/\__\__,_|\__|_|\___| |_| |_| |_|\___|\__|_| |_|\___/ \__,_|___/
    %************************************************************
    
    methods(Static)        
        %% ********************************************************
        % makeANsmooth - smooth the AN response into hanning windowed chunks
        %**********************************************************        
        function ANsmooth = makeANsmooth(ANresponse, sampleRate, winSize, hopSize)
            if nargin < 3
                winSize = 25; %default 25 ms window
            end
            if nargin < 4
                hopSize = 10; %default 10 ms jump between windows
            end
            
            winSizeSamples = round(winSize*sampleRate/1000);
            hopSizeSamples = round(hopSize*sampleRate/1000);
            
            % smooth
            hann = GJB_hanning(winSizeSamples);
            
            ANsmooth = [];%Cannot pre-allocate a size as it is unknown until the enframing
            for chan = 1:size(ANresponse,1)
                f = enframe(ANresponse(chan,:), hann, hopSizeSamples);
                ANsmooth(chan,:) = mean(f,2)'; %#ok<AGROW> see above comment
            end
        end% ------ OF makeANsmooth
        
        %% ********************************************************
        % makeANfeatures - dct wizardry
        %********************************************************** 
        function ANfeatures = makeANfeatures(ANrate, numCoeff)
            % make feature vectors
            features = GJB_dct(ANrate);
            ANfeatures = features(1:numCoeff,:);            
        end % ------ OF makeANfeatures
        
    end % ------ OF STATIC METHODS
    
end % ------ OF CLASS
