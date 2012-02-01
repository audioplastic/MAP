function vectorStrength=testAN(probeFrequency, channleBF, levels, ...
    paramsName,paramChanges)
% testAN generates rate/level functions for AN and brainstem units.
%  also other information like PSTHs, MOC efferent activity levels,
%  phase-locking and coefficient of variation
%
% A full 'spikes' model is used.For probability model use 'testANprob'
%
% Vector strength calculations require the computation of period
%  histograms. For this reason the sample rate must always be an integer
% multiple of the tone frequency. This applies to both the sampleRate and
% the spikes sampling rate.
%
% input arguments (can be omitted from the right):
%  probeFrequency: the test stimulus is a pure tone (default=1000)
%  channleBF: this is a single channel model (like a physiological model
%    with single electrode) default=1000
%  levels: (vector) tone levels to be tested. (default= -10:10:80).
%  paramsName: parameter file name containing model parameters.
%   (default='Normal')
%  paramChanges: cell array contining list of changes to parameters. These
%   are implemented after reading the parameter file (default='')
%
% Output arguments 
% vector strength (phase-locking estimate at the levels requested)
% 
%   Example
%   testAN(1000,1000, -10:10:80,'Normal','');


global IHC_cilia_RPParams IHCpreSynapseParams
global AN_IHCsynapseParams
global ANoutput dtSpikes CNoutput ICoutput ICmembraneOutput ANtauCas
global ARattenuation MOCattenuation

dbstop if error
restorePath=path;
addpath (['..' filesep 'MAP'], ['..' filesep 'utilities'], ...
    ['..' filesep 'parameterStore'],  ['..' filesep 'wavFileStore'])

if nargin<5, paramChanges=[]; end
if nargin<4, paramsName='Normal'; end
if nargin<3, levels=-10:10:80; end
if nargin==0, probeFrequency=1000; channleBF=1000; end
nLevels=length(levels);

% input signal = <silence> <tone>
toneDuration=.2;   
rampDuration=0.002;   
silenceDuration=.02;
localPSTHbinwidth=0.0005;

% guarantee that the sample rate is an interger multiple 
%   of the probe frequency
% We want 5 bins per period for spikes to allow phase-locking to be
%   calculated
spikesSampleRate=5*probeFrequency;
% model sample rate must be in the region of 50000
sampleRate=spikesSampleRate*round(50000/spikesSampleRate);
dt=1/sampleRate;

% avoid very slow *spikes* sampling rate
% NB spikes sample rate is lower than the sample rate used to evaluate the
%   model upt to the level of the IHC.
spikesSampleRate=spikesSampleRate*ceil(10000/spikesSampleRate);

%% pre-allocate storage
AN_HSRonset=zeros(nLevels,1);
AN_HSRsaturated=zeros(nLevels,1);
AN_LSRonset=zeros(nLevels,1);
AN_LSRsaturated=zeros(nLevels,1);
CNLSRsaturated=zeros(nLevels,1);
CNHSRsaturated=zeros(nLevels,1);
ICHSRsaturated=zeros(nLevels,1);
ICLSRsaturated=zeros(nLevels,1);
vectorStrength=zeros(nLevels,1);
AR=zeros(nLevels,1);
MOC=zeros(nLevels,1);

% fix figure locations for the benefit of multiThreshold use
figure(15), clf
set(gcf,'position',[980   356   401   321])
figure(5), clf
set(gcf,'position', [980 34 400 295])
drawnow


%% main computational loop (vary level)
levelNo=0;
for leveldB=levels
    levelNo=levelNo+1;
    fprintf('%4.0f\t', leveldB)

    %% generate tone and silences
    time=dt:dt:toneDuration;
    rampTime=dt:dt:rampDuration;
    ramp=[0.5*(1+cos(2*pi*rampTime/(2*rampDuration)+pi)) ...
        ones(1,length(time)-length(rampTime))];
    ramp=ramp.*fliplr(ramp);
    
    silence=zeros(1,round(silenceDuration/dt));
    
    amp=28e-6*10^(leveldB/20);
    tone=amp*sin(2*pi*probeFrequency'*time);
    tone= ramp.*tone;

    inputSignal=[silence tone];
    
    %% run the model
    AN_spikesOrProbability='spikes';  
    nExistingParamChanges=length(paramChanges);
    paramChanges{nExistingParamChanges+1}=...
        ['AN_IHCsynapseParams.spikesTargetSampleRate=' ...
        num2str(spikesSampleRate) ';'];

    MAP1_14(inputSignal, 1/dt, channleBF, ...
        paramsName, AN_spikesOrProbability, paramChanges);
    
    nFiberTypes=length(ANtauCas);
    if ~isequal(nFiberTypes, 2)
        error('testAN: both LSR and HSR fiber types must be specified')
    end
    
    %% Auditory nerve evaluate and display (Fig. 5)
    %LSR (same as HSR if no LSR fibers present)
    nANFibers=size(ANoutput,1);
    numLSRfibers=nANFibers/nFiberTypes;
    numHSRfibers=numLSRfibers;
    
    % LSR response is at the top of thematrix
    LSRspikes=ANoutput(1:numLSRfibers,:);
    PSTH=UTIL_PSTHmaker(LSRspikes, dtSpikes, localPSTHbinwidth);
    PSTHLSR=mean(PSTH,1)/localPSTHbinwidth;  % across fibers rates
    PSTHtime=localPSTHbinwidth:localPSTHbinwidth:...
        localPSTHbinwidth*length(PSTH);
    AN_LSRonset(levelNo)= max(PSTHLSR); % peak in 5 ms window
    AN_LSRsaturated(levelNo)= mean(PSTHLSR(round(length(PSTH)/2):end));
    
    % AN HSR
    HSRspikes= ANoutput(end- numHSRfibers+1:end, :);
    PSTH=UTIL_PSTHmaker(HSRspikes, dtSpikes, localPSTHbinwidth);
    PSTH=mean(PSTH,1)/localPSTHbinwidth; % sum across fibers (HSR only)
    AN_HSRonset(levelNo)= max(PSTH);
    AN_HSRsaturated(levelNo)= mean(PSTH(round(length(PSTH)/2): end));
    
    figure(5), subplot(2,2,2)
    hold off, bar(PSTHtime,PSTH, 'b')
    hold on,  bar(PSTHtime,PSTHLSR,'r')
    ylim([0 1000])
    xlim([0 length(PSTH)*localPSTHbinwidth])
    grid on
    set(gcf,'name',['PSTH: ' num2str(channleBF), ' Hz: ' num2str(leveldB) ' dB']);
    
    % AN - coefficient of variation
    %  CV is computed 5 times (at different time points)
    %  Use the middle one (3) as most typical
    cvANHSR=  UTIL_CV(HSRspikes, dtSpikes);
    
    % AN - vector strength (phase-locking)
    PSTH=sum(HSRspikes);
    [PH]=UTIL_periodHistogram...
        (PSTH, dtSpikes, probeFrequency);
    VS=UTIL_vectorStrength(PH);
    vectorStrength(levelNo)=VS;
    disp(['sat rate= ' num2str(AN_HSRsaturated(levelNo)) ...
        ';   phase-locking VS = ' num2str(VS)])
    title(['AN HSR: CV=' num2str(cvANHSR(3),'%5.2f') ...
        'VS=' num2str(VS,'%5.2f')])
    
    %% CN - first-order neurons    
    %   CN driven by LSR fibers
    nCNneurons=size(CNoutput,1);
    nLSRneurons=round(nCNneurons/nFiberTypes);
    CNLSRspikes=CNoutput(1:nLSRneurons,:);
    PSTH=UTIL_PSTHmaker(CNLSRspikes, dtSpikes, localPSTHbinwidth);
    PSTH=sum(PSTH)/nLSRneurons;
    CNLSRsaturated(levelNo)=mean(PSTH)/localPSTHbinwidth;
    
    %CN driven by HSR fibers
    MacGregorMultiHSRspikes=...
        CNoutput(end-nLSRneurons+1:end,:);
    PSTH=UTIL_PSTHmaker(MacGregorMultiHSRspikes, dtSpikes, localPSTHbinwidth);
    PSTH=sum(PSTH)/nLSRneurons;
    PSTH=mean(PSTH,1)/localPSTHbinwidth; % sum across fibers (HSR only)
    CNHSRsaturated(levelNo)=mean(PSTH);
    
    figure(5), subplot(2,2,3)
    bar(PSTHtime,PSTH)
    ylim([0 1000])
    xlim([0 length(PSTH)*localPSTHbinwidth])
    cvMMHSR= UTIL_CV(MacGregorMultiHSRspikes, dtSpikes);
    title(['CN    CV= ' num2str(cvMMHSR(3),'%5.2f')])
    
    %% IC - second order brainstem neuron
    %   IC driven by LSR fibers
    nICneurons=size(ICoutput,1);
    nLSRneurons=round(nICneurons/nFiberTypes);
    ICLSRspikes=ICoutput(1:nLSRneurons,:);
    PSTH=UTIL_PSTHmaker(ICLSRspikes, dtSpikes, localPSTHbinwidth);
    ICLSRsaturated(levelNo)=mean(PSTH)/localPSTHbinwidth;
    
    %   IC driven by HSR fibers
    MacGregorMultiHSRspikes=...
        ICoutput(end-nLSRneurons+1:end,:);
    PSTH=UTIL_PSTHmaker(MacGregorMultiHSRspikes, dtSpikes, localPSTHbinwidth);
    ICHSRsaturated(levelNo)= (sum(PSTH)/nLSRneurons)/toneDuration;

    % efferent effects
    AR(levelNo)=min(ARattenuation);
    MOC(levelNo)=min(MOCattenuation(length(MOCattenuation)/2:end));
    
    time=dt:dt:dt*size(ICmembraneOutput,2);
    figure(5), subplot(2,2,4)
    % plot HSR (last of two)
    plot(time,ICmembraneOutput(end-nLSRneurons+1, 1:end),'k')
    ylim([-0.07 0])
    xlim([0 max(time)])
    title(['IC  ' num2str(leveldB,'%4.0f') 'dB'])
    drawnow
    
    figure(5), subplot(2,2,1)
    plot(20*log10(MOC), 'k'), hold on
    plot(20*log10(AR), 'r'), hold off
    title(' MOC'), ylabel('dB attenuation')
    ylim([-30 0])
    
end % level

%% plot with levels on x-axis
figure(5), subplot(2,2,1)
plot(levels,20*log10(MOC), 'k'),hold on
plot(levels,20*log10(AR), 'r'), hold off
title(' MOC'), ylabel('dB attenuation')
ylim([-30 0])
xlim([0 max(levels)])

fprintf('\n')
toneDuration=2;
nRepeats=200;   % no. of AN fibers
fprintf('toneDuration  %6.3f\n', toneDuration)
fprintf(' %6.0f  AN fibers (repeats)\n', nRepeats)
fprintf('levels')
fprintf('%6.2f\t', levels)
fprintf('\n')

%% ---------------------- display summary results (Fig 15)
figure(15), clf
nRows=2; nCols=2;

% AN rate - level ONSET functions
subplot(nRows,nCols,1)
plot(levels,AN_LSRonset,'ro'), hold on
plot(levels,AN_HSRonset,'ko', 'MarkerEdgeColor','k', 'markerFaceColor','k'), hold off
ylim([0 1000])
if length(levels)>1
    xlim([min(levels) max(levels)])
end
myTitle=['tauCa= ' num2str(IHCpreSynapseParams.tauCa)];
title( myTitle)
xlabel('level dB SPL'), ylabel('peak rate (sp/s)'), grid on
text(0, 800, 'AN onset', 'fontsize', 14)

% AN rate - level ADAPTED function
subplot(nRows,nCols,2)
plot(levels,AN_LSRsaturated, 'ro'), hold on
plot(levels,AN_HSRsaturated, 'ko', 'MarkerEdgeColor','k',...
    'markerFaceColor','k'), hold off
ylim([0 400])
set(gca,'ytick',0:50:300)
if length(levels)>1
xlim([min(levels) max(levels)])
end
set(gca,'xtick',levels(1):20:levels(end))
%     grid on
myTitle=[   'spont=' num2str(mean(AN_HSRsaturated(1,:)),'%4.0f')...
    '  sat=' num2str(mean(AN_HSRsaturated(end,1)),'%4.0f')];
title( myTitle)
xlabel('level dB SPL'), ylabel ('adapted rate (sp/s)')
text(0, 340, 'AN adapted', 'fontsize', 14), grid on

% CN rate - level function
subplot(nRows,nCols,3)
plot(levels,CNLSRsaturated, 'ro'), hold on
plot(levels,CNHSRsaturated, 'ko', 'MarkerEdgeColor','k',...
    'markerFaceColor','k'), hold off
ylim([0 400])
set(gca,'ytick',0:50:300)
if length(levels)>1
xlim([min(levels) max(levels)])
end
set(gca,'xtick',levels(1):20:levels(end))
%     grid on
myTitle=[   'spont=' num2str(mean(CNHSRsaturated(1,:)),'%4.0f') ' sat=' ...
    num2str(mean(CNHSRsaturated(end,1)),'%4.0f')];
title( myTitle)
xlabel('level dB SPL'), ylabel ('adapted rate (sp/s)')
text(0, 350, 'CN', 'fontsize', 14), grid on

% IC rate - level function
subplot(nRows,nCols,4)
plot(levels,ICLSRsaturated, 'ro'), hold on
plot(levels,ICHSRsaturated, 'ko', 'MarkerEdgeColor','k',...
    'markerFaceColor','k'), hold off
ylim([0 400])
set(gca,'ytick',0:50:300)
if length(levels)>1
xlim([min(levels) max(levels)])
end
set(gca,'xtick',levels(1):20:levels(end)), grid on
myTitle=['spont=' num2str(mean(ICHSRsaturated(1,:)),'%4.0f') ...
    '  sat=' num2str(mean(ICHSRsaturated(end,1)),'%4.0f')];
title( myTitle)
xlabel('level dB SPL'), ylabel ('adapted rate (sp/s)')
text(0, 350, 'IC', 'fontsize', 14)
set(gcf,'name',' AN CN IC rate/level')

% command window summary
fprintf('\n')
disp('levels vectorStrength')
fprintf('%3.0f \t %6.4f \n', [levels; vectorStrength'])
fprintf('\n')
fprintf('Phase locking, max vector strength=\t %6.4f\n\n',...
    max(vectorStrength))

allData=[ levels'  AN_HSRonset AN_HSRsaturated...
    AN_LSRonset AN_LSRsaturated ...
    CNHSRsaturated CNLSRsaturated...
    ICHSRsaturated ICLSRsaturated];
fprintf('\n levels \tANHSR Onset \tANHSR adapted\tANLSR Onset \tANLSR adapted\tCNHSR\tCNLSR\tICHSR  \tICLSR \n');
UTIL_printTabTable(round(allData))
fprintf('VS (phase locking)= \t%6.4f\n\n',...
    max(vectorStrength))

UTIL_showStruct(IHC_cilia_RPParams, 'IHC_cilia_RPParams')
UTIL_showStruct(IHCpreSynapseParams, 'IHCpreSynapseParams')
UTIL_showStruct(AN_IHCsynapseParams, 'AN_IHCsynapseParams')

path(restorePath)
