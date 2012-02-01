function testANprob (targetFrequency,BFlist, levels, ...
    paramsName, paramChanges)
% testANprob generates rate/level functions for AN.
%  also other information like PSTHs, MOC efferent activity levels.
%
% An AN 'probability' model is used.For spikes model use 'testAN'
%
% Input arguments (can be omitted from the right):
%  probeFrequency: the test stimulus is a pure tone (default=1000)
%  channleBF: this is a single channel model (like a physiological model
%    with single electrode) default=1000.
%    This means that the acoustic reflex is not fully functioning because
%    it depends on a broadband assessment of overall level
%  levels: (vector) tone levels to be tested. (default= -10:10:80).
%  paramsName: parameter file name containing model parameters.
%   (default='Normal')
%    NB the program assumes that two fiber types are nominated, i.e. two
%    values of ANtauCas are specified.
%  paramChanges: cell array contining list of changes to parameters. These
%   are implemented after reading the parameter file (default='')
% 
%   Example
%   testANprob(1000,1000, -10:10:80,'Normal','');

global IHC_cilia_RPParams IHCpreSynapseParams AN_IHCsynapseParams
% these globals pick up the results after running MAP1_14
global ANprobRateOutput dt ANtauCas ARattenuation MOCattenuation

dbstop if error
restorePath=path;
addpath (['..' filesep 'MAP'], ['..' filesep 'utilities'], ...
    ['..' filesep 'parameterStore'],  ['..' filesep 'wavFileStore'],...
    ['..' filesep 'testPrograms'])

% defaults if trailing arguments are omitted.
if nargin<5, paramChanges=[]; end
if nargin<4, paramsName='Normal'; end
if nargin<3, levels=-10:10:80; end
if nargin==0, targetFrequency=1000; BFlist=1000; end

if length(BFlist)>1
    error('testANprob: only one channel is allowed. BFlist is too long')
end

nLevels=length(levels);

% input= <silence> <tone>
toneDuration=.200;          % all in seconds
rampDuration=0.002;
silenceDuration=.02;
localPSTHbinwidth=0.001;    % for AN PSTH

AN_HSRonset=zeros(nLevels,1);
AN_HSRsaturated=zeros(nLevels,1);
AN_LSRonset=zeros(nLevels,1);
AN_LSRsaturated=zeros(nLevels,1);
AR=zeros(nLevels,1);
MOC=zeros(nLevels,1);

figure(15), clf
% Fix position to allow multiThreshold to display all with no overlap
set(gcf,'position',[980   356   401   321])
drawnow

% Guarantee that the sample rate is at least 10 times the frequency
sampleRate=50000;
while sampleRate< 10* targetFrequency
    sampleRate=sampleRate+10000;
end

% Adjust sample rate so that the pure tone stimulus has an integer
%   number of epochs in a period
dt=1/sampleRate;

%% main computational loop (level is varied)
levelNo=0;
for leveldB=levels
    levelNo=levelNo+1;
    fprintf('%4.0f\t', leveldB)

    time=dt:dt:toneDuration;
    rampTime=dt:dt:rampDuration;
    ramp=[0.5*(1+cos(2*pi*rampTime/(2*rampDuration)+pi)) ...
        ones(1,length(time)-length(rampTime))];
    ramp=ramp.*fliplr(ramp);

    % initial silence
    silence=zeros(1,round(silenceDuration/dt));

    % create tone signal (leveldB/ targetFrequency)
    amp=28e-6*10^(leveldB/20);
    tone=amp*sin(2*pi*targetFrequency'*time);
    tone= ramp.*tone;
    inputSignal=[silence tone];

    %% run the model
    AN_spikesOrProbability='probability';

    MAP1_14(inputSignal, 1/dt, BFlist, ...
        paramsName, AN_spikesOrProbability, paramChanges);

    nFiberTypes=length(ANtauCas);

    if ~isequal(nFiberTypes,2)
        error ('testANprob: model should feature both HSR and LSR fibers')
    end

    % by default there is only one fiber per channel of each type
    % LSR fiber is the first row
    LSRspikes=ANprobRateOutput(1,:);
    PSTH=UTIL_PSTHmaker(LSRspikes, dt, localPSTHbinwidth);
    PSTHLSR=PSTH/(localPSTHbinwidth/dt);  % across fibers rates
    PSTHtime=localPSTHbinwidth:localPSTHbinwidth:...
        localPSTHbinwidth*length(PSTH);
    AN_LSRonset(levelNo)= max(PSTHLSR); % peak in 5 ms window
    AN_LSRsaturated(levelNo)= mean(PSTHLSR(round(length(PSTH)/2):end));

    % HSR
    HSRspikes= ANprobRateOutput(2, :);
    PSTH=UTIL_PSTHmaker(HSRspikes, dt, localPSTHbinwidth);
    PSTH=PSTH/(localPSTHbinwidth/dt); % sum across fibers (HSR only)
    AN_HSRonset(levelNo)= max(PSTH);
    AN_HSRsaturated(levelNo)= mean(PSTH(round(length(PSTH)/2): end));

    figure(15), subplot(2,2,4)
    hold off, bar(PSTHtime,PSTH, 'k')
    hold on,  bar(PSTHtime,PSTHLSR,'r')
    ylim([0 1000])
    xlim([0 length(PSTH)*localPSTHbinwidth])
    set(gcf,'name',[num2str(BFlist), ' Hz: ' num2str(leveldB) ' dB']);

    AR(levelNo)=min(ARattenuation);
    MOC(levelNo)=min(MOCattenuation(length(MOCattenuation)/2:end));

    figure(15), subplot(2,2,3)
    plot(20*log10(MOC), 'k'), hold on
    plot(20*log10(AR), 'r'),  hold off
    title(' MOC/AR'), ylabel('dB attenuation')
    ylim([-30 0])

end % level

%% displays
figure(15), subplot(2,2,3)
plot(levels,20*log10(MOC), 'k'), hold on
plot(levels,20*log10(AR), 'r'),  hold off
title(' MOC/AR'), ylabel('dB attenuation')
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

%% ---------------------------------------------------- display parameters
% subplot arrangements
nRows=2; nCols=2;

% AN rate - level ONSET functions
subplot(nRows,nCols,1)
plot(levels,AN_LSRonset,'ro'), hold on
plot(levels,AN_HSRonset,'ko', 'MarkerEdgeColor','k', 'markerFaceColor','k'), hold off
ylim([0 1000]),  xlim([min(levels) max(levels)])
myTitle=['tauCa= ' num2str(IHCpreSynapseParams.tauCa)];
title( myTitle)
xlabel('level dB SPL'), ylabel('peak rate (sp/s)'), grid on
text(0, 800, 'AN onset', 'fontsize', 14)

% AN rate - level ADAPTED function
subplot(nRows,nCols,2)
plot(levels,AN_LSRsaturated, 'ro'), hold on
plot(levels,AN_HSRsaturated, 'ko', 'MarkerEdgeColor','k', 'markerFaceColor','k'), hold off
maxYlim=340;
ylim([0 maxYlim])
set(gca,'ytick',0:50:300)
xlim([min(levels) max(levels)])
set(gca,'xtick',levels(1):20:levels(end))
%     grid on
myTitle=[   'spont=' num2str(mean(AN_HSRsaturated(1,:)),'%4.0f')...
    '  sat=' num2str(mean(AN_HSRsaturated(end,1)),'%4.0f')];
title( myTitle)
xlabel('level dB SPL'), ylabel ('adapted rate (sp/s)')
text(0, maxYlim-50, 'AN adapted', 'fontsize', 14), grid on

% final results table
allData=[ levels'  AN_HSRonset AN_HSRsaturated...
    AN_LSRonset AN_LSRsaturated ];
fprintf('\n levels \tANHSR Onset \tANHSR adapted\tANLSR Onset \tANLSR adapted \n');
UTIL_printTabTable(round(allData))


UTIL_showStruct(IHC_cilia_RPParams, 'IHC_cilia_RPParams')
UTIL_showStruct(IHCpreSynapseParams, 'IHCpreSynapseParams')
UTIL_showStruct(AN_IHCsynapseParams, 'AN_IHCsynapseParams')


allData=[ levels'  AN_HSRonset AN_HSRsaturated...
    AN_LSRonset AN_LSRsaturated ];
fprintf('\n levels \tANHSR Onset \tANHSR adapted\tANLSR Onset \tANLSR adapted\tCNHSR\tCNLSR\tICHSR  \tICLSR \n');
UTIL_printTabTable(round(allData))

path(restorePath)

