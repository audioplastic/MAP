% function efferentTest
% this program is used to explore the acoustic reflex system
% see Moss PhD for details and Hung and Dallos (1972)
%
% If the level is set to 110 dB SPL, an additional figure is shown to
% compare the time course of the effect with Hung and Dallos. for this to
% work, the end silence should be substantial

global  inputStimulusParams outerMiddleEarParams IHCpreSynapseParams
global AN_IHCsynapseParams  DRNLParams
global MacGregorMultiParams  MacGregorParams

global ANprobRateOutput  ANdt

paramsFileName='Normal';

% stimulus details
signalType='tone';
sampleRate=44100;
toneFrequency=1000;
% toneFrequency=[ 500];
% toneFrequency=[ 250];
% toneFrequency=[ 4000];

% % or
% signalType='noise';
% frequencies=[0 ];

% signalLevels=[0:10:100];
% signalLevels=[80:10:100];

signalLevels=[90 100];

beginSilenceDuration= 0.05;
stimulusDuration= 0.2;
endSilenceDuration=0.05;

rampDuration=0.005;

% %single shot:  NB 101 dB triggers plot of AR compared with data
% % at 101 dB SPL compare onset and offset reactance trajectories of
% % model with data from Hung and Dallos (1972) Fig. 1.
% signalLevels=101;
% signalType='noise';
% beginSilenceDuration= 0.2;
% stimulusDuration= 0.55;
% endSilenceDuration=0.5;
% rampDuration=0.02;

% %single shot:  NB 61 dB triggers plot of MOC compared with data
% % at 61 dB SPL compare onset and offset reactance trajectories of
% % model with data from Hung and Dallos (1972) Fig. 1.
% signalLevels=61;
% signalType='noise';
% beginSilenceDuration= 0.2;
% stimulusDuration= 0.55;
% endSilenceDuration=0.5;
% rampDuration=0.02;
% segmentDuration=0.005;

% default sequence of modules to IC
moduleSequence= [1 2 3 4 5 6 7 8 ];  	

lowestBF=250; 	highestBF= 8000; 	numChannels=21;
% includes BFs at 250 500 1000 2000 4000 8000
BFlist=round(logspace(log10(lowestBF), log10(highestBF), numChannels));

% identify closest channel to tone frequency
% [euclid bestCh]= min((BFlist-toneFrequency).^2);
[x bestCh]=min((BFlist-toneFrequency).^2);
BFused=BFlist(bestCh);

% create access to all MAP 1_8 facilities
addpath ('..\modules', '..\utilities', ...
    '..\parameterStore',  '..\wavFileStore' , '..\testPrograms')
figure(14), clf, set (gcf,'name',paramsFileName, ...
    'position', [19  277  560  420] )
figure(15), clf, set (gcf,'name','efferent cross-channel',...
    'position', [1030    38   333   278])
figure(16), clf, set (gcf,'name','AR MOC control',...
    'position', [1043   455   316   242])

doPlot=1; % MAP plotting (normally off)

ARattenSummary=[];earObject=[];allLevels=[];
nFrequencies=1;
nLevels=length(signalLevels);

rateMacGHSR=[]; rateMacGLSR=[]; rateANLSR=[]; rateANHSR=[]; 
signalLevelsSoFar=[];     mocATTdBs=[]; ARdBs=[];

for leveldB=signalLevels
    disp(['leveldB=' num2str(leveldB)])
    signalLevelsSoFar=[signalLevelsSoFar; leveldB];
    
    % create stimulus
    plotStimulus=0;
    stim=[];
    stim.type=signalType;
    stim.phases='sin';
    stim.toneDuration=stimulusDuration;
    stim.frequencies=toneFrequency;
    stim.amplitudesdB=leveldB;
    stim.beginSilence=beginSilenceDuration;
    stim.endSilence=endSilenceDuration;
    stim.rampOnDur=rampDuration;
    stim.rampOffDur=rampDuration;
    switch signalType
        case 'noise'
            stim.filter=[100 10000 2]; %bandwidth of noise
    end
    globalStimParams.FS=sampleRate;
    overallDuration=...
        stimulusDuration+beginSilenceDuration+endSilenceDuration;  % s
    globalStimParams.overallDuration=overallDuration;
    
    inputSignal=stimulusCreate(globalStimParams, stim, plotStimulus);
    inputSignal=inputSignal(:,1)';  % mono, row vector
    
    %Overall level (dB SPL) = Spectrum level (dB SPL/Hz) + 10*log10(Noise
    %bandwidth)
    switch signalType
        case 'noise'
            overallToSpectrumLevel=10*log10(stim.filter(2)-stim.filter(1));
            amp=10^(overallToSpectrumLevel/20);
            inputSignal=inputSignal*amp;
    end
    
    
    % parameters

cmd=['method=MAPparams' paramsFileName '();'];
eval(cmd)

%     inputStimulusParams.sampleRate=sampleRate;
    method.dt=1/sampleRate;
    
    method.plotGraphs=	doPlot;
    method.showSummaryStatistics=0;
    method.AN_IHCsynapseSave=1;  % special save for summary
    method.MacGregorSave=1;  % special save for summary
    
    % run the model
%     [earObjectThisSegment, method]=...
%         MAPsequenceSeg(inputSignal, method, moduleSequence);
    
paramChanges = {};
MAP1_14(inputSignal, sampleRate, BFlist, paramsFileName, 'probability', paramChanges)
    
    
        
    % AN LSR summary
%     [totalNfibers nDataPoints]=size(method.AN_IHCsynapseData);
    [totalNfibers nDataPoints]=size(ANprobRateOutput);
    fibersPerChannel=method.numANfibers;
    nLSRfibers=round(totalNfibers/length(method.tauCa));
    nHSRfibers=nLSRfibers;
    % use all LSR fibers i.e. top half of AN spike matrix
    % mean spike rate across all channels
    spikeRate=sum(method.AN_IHCsynapseData(1:nLSRfibers,:),2);
    spikeRate=sum(spikeRate)/ (nLSRfibers)/globalStimParams.overallDuration;   
    % or only spikeRate at BF
    %  BFrate=method.AN_IHCsynapseData((bestCh-1)*...
    %     fibersPerChannel+1:bestCh*fibersPerChannel,:);
    %  spikeRate=mean(sum(BFrate,2)/ globalStimParams.overallDuration);
    rateANLSR=[rateANLSR; spikeRate];
    disp(['ANLSR spikeRate= ' num2str(spikeRate)])
    
    % AN HSR summary        use lower half of matrix
    % mean spike rate across all channels
    spikeRate=sum(method.AN_IHCsynapseData(nLSRfibers+1:end,:));
    spikeRate=sum(spikeRate)/ (nLSRfibers)/ globalStimParams.overallDuration;
    % % or spikeRate at BF
    %  BFrate=method.AN_IHCsynapseData(nLSRfibers+(bestCh-1)*...
    %   fibersPerChannel+1:nLSRfibers+bestCh*fibersPerChannel,:);
    %  spikeRate=mean(sum(BFrate,2)/globalStimParams.overallDuration);
    rateANHSR=[rateANHSR; spikeRate];
    disp(['ANHSR spikeRate= ' num2str(spikeRate)])
    
    % plot AN HSR and LSRrates
    figure(14), subplot(nFrequencies,4, 1)
    plot(signalLevelsSoFar,rateANHSR,'ro'), hold on
    plot(signalLevelsSoFar,rateANLSR,'o'), hold off
    if length(signalLevels)>1,     xlim([min(signalLevels) max(signalLevels)]), end
    theTitle=char(signalType);
    if strcmp(theTitle,'tone')
        theTitle=[theTitle ':  ' num2str(toneFrequency')];
    end
    title([ 'AN:' theTitle])
    xlabel('level dB')
    ylabel('across-channelmean spike rate (s/s)')
    ylim([0 200])
    
    % MacGregor - LSR rates
    [nMacGcells nDataPoints]=size(method.MacGregorData);
    nMacGcells=round(nMacGcells/2); % per fiber type
    spikeRate=sum(sum(method.MacGregorData(1:nMacGcells,:)))...
        / (nMacGcells*overallDuration);
    rateMacGLSR=[rateMacGLSR; spikeRate];
    disp(['MacG LSR spikeRate= ' num2str(spikeRate)])
    % MacGregor - HSR rates
    spikeRate=sum(sum(method.MacGregorData(nMacGcells+1:end,:)))...
        / (nMacGcells*overallDuration);
    rateMacGHSR=[rateMacGHSR; spikeRate];
    disp(['MacG HSR spikeRate= ' num2str(spikeRate)])
    % plot both
    figure(14), subplot(nFrequencies,4, 2)
    plot(signalLevelsSoFar,rateMacGHSR,'ro'), hold on
    plot(signalLevelsSoFar,rateMacGLSR,'o'), hold off
    %         legend({'MacG HSR', 'MacG LSR'}, 'location','northwest')
    %             ylim([0 250])
    if length(signalLevels)>1,     xlim([min(signalLevels) max(signalLevels)]), end
    title([ 'MacG' ])
    xlabel('level dB')
    ylabel('across-channelmean spike rate (s/s)')
    ylim([0 100])

    % MOC ATTdBattenuation
    mocATTdBs=[mocATTdBs; method.mocATTdB(bestCh,:)];
    dt=method.mocATTdBdt;
    time=dt:dt:dt*length(method.mocATTdB);
    figure(14), subplot(nFrequencies,4,4)
    plot(time,-mocATTdBs)
    ylabel('attenuation'),        ylim([0 60])
    xlabel('time')
    xlim([0 overallDuration])
    grid on
    title(['MOC at ' num2str(BFused) ' Hz'])
    legend(num2str(signalLevelsSoFar),'location','northwest')

    figure(15)
    surf(time, method.nonlinCF,  -method.mocATTdB)
    shading interp
    ylabel(' channel BF')
    xlabel(' time (s)')
    title ([signalType  ' MOC efferent attenuation (dB). signal level ' ...
        num2str(leveldB) ' dB'])

    timeOnPTR=round((beginSilenceDuration+.05)/dt);
    timeOffPTR= round((stimulusDuration+beginSilenceDuration)/dt);
    %         To compare with Liberman, the average effect is used
    MOCattenSummary=mean(-mocATTdBs(:,timeOnPTR:timeOffPTR)')';
    ylim([0 BFlist(end)])
    xlim([0 time(end)])
    zlim([0 40])
    
    figure(16), subplot(2,1,2),
    plot(signalLevelsSoFar,MOCattenSummary)
    xlim([signalLevels(1) signalLevels(end)])
    title('MOC')

    % AR attenuation
    dt=method.outerMiddleEardt;
    ARtime=dt:dt:dt*length(method.ARdBAttenuation);
    figure(14), subplot(nFrequencies,4,3)
    %         hold on
    ARdBs=[ARdBs;method.ARdBAttenuation];
    plot(ARtime,-ARdBs)
    ylabel('attenuation'),     ylim([0 60])
    xlabel('time')
    xlim([0 overallDuration])
    grid on
    title('AR:  ' )
    
    ARattenSummary=mean(-ARdBs(:,timeOnPTR:timeOffPTR)')';
    figure(16), subplot(2,1,1),
    plot(signalLevelsSoFar,ARattenSummary)
    title ('AR')
    xlim([signalLevels(1) signalLevels(end)])
    
    % special plot to compare AR with reactance data
    if leveldB==101 && strcmp(signalType, 'noise')
        % at 101 dB SPL compare onset and offset reactance trajectories of
        % model with data from Hung and Dallos (1972) Fig. 1.
        figure(15),
        timeOnPTR=round((beginSilenceDuration+.05)/dt);
        timeOffPTR= round((stimulusDuration+beginSilenceDuration)/dt);
        levelDuringSignal=-method.ARdBAttenuation(timeOnPTR:timeOffPTR);
        plot(ARtime,-method.ARdBAttenuation/mean(levelDuringSignal))
        
        %data
        onset400_110dB=[
            .045	0
            .060	25
            .075	50
            .090	75
            .110	100
            .150	110];
        onset400_110dB(:,1)=onset400_110dB(:,1)+beginSilenceDuration -.025;
        offset=[
            0	1
            .1	0.8
            .2	0.5
            .3	0.25
            .4	0.15
            .5	0.1];
        
        offsetBeginsNow=beginSilenceDuration+stimulusDuration;
        offset(:,1)=offset(:,1)+offsetBeginsNow;
        
        hold on
        plot(onset400_110dB(:,1), ...
            onset400_110dB(:,2)/max(onset400_110dB(:,2)),'ro')
        plot(offset(:,1), offset(:,2), 'ko')
        
        normalseFactor=-method.ARdBAttenuation(timeOffPTR);
        plot(ARtime,-method.ARdBAttenuation/normalseFactor,'r')

        hold off
        title(['AR reactance data at 101 dB SPL. Ramp is ' ...
            num2str(rampDuration, '%6.3f')])
    end
    
UTIL_printTabTable( round([signalLevelsSoFar rateANHSR rateANLSR ...
    rateMacGHSR rateMacGLSR MOCattenSummary ARattenSummary ]),...
    strvcat( 'signalLevelsSoFar', 'rateANHSR', 'rateANLSR', 'rateMacGHSR', ...
    'rateMacGLSR', 'MOCattenSummary', 'ARattenSummary'));
    drawnow
end %signalLevels


% write all parameters to the command window
disp('OME filters')
disp(num2str(outerMiddleEarParams.MElowpass))
nm=UTIL_paramsList(whos);
for i=1:length(nm)
    eval(['UTIL_showStruct(' nm{i} ', ''' nm{i} ''')'])
end

UTIL_printTabTable( round([signalLevelsSoFar rateANHSR rateANLSR ...
    rateMacGHSR rateMacGLSR MOCattenSummary ARattenSummary ]),...
    strvcat( 'signalLevelsSoFar', 'rateANHSR', 'rateANLSR', 'rateMacGHSR', ...
    'rateMacGLSR', 'MOCattenSummary', 'ARattenSummary'));
