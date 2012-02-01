function testEfferent(probeFrequencies,BFlist, levels, paramsName, ...
    paramChanges)
%
% testEfferent generates rate/level functions for AR and MOC
%  A multi-channel model is needed for this
%
%Input arguments:
%  probeFrequencies: vector
%  BFlist: BFs of multichannel model
%  levels: vector of levels to be used in rate/level function
%     all probe frequencies are played at equal levels
%  paramsName: parameter file name containing model parameters.
%   (default='Normal')
%    NB the program assumes that two fiber types are nominated, i.e. two
%    values of ANtauCas are specified.
%  paramChanges: cell array contining list of changes to parameters. These
%   are implemented after reading the parameter file (default='')
%
% Example
%   testEfferent(1000,[250 500 1000 2000 4000], -10:10:80,'Normal',[]);

global  MOCattenuation IHCpreSynapseParams

dbstop if error
restorePath=path;
addpath (['..' filesep 'MAP'], ['..' filesep 'utilities'], ...
    ['..' filesep 'parameterStore'],  ['..' filesep 'wavFileStore'],...
    ['..' filesep 'testPrograms'])

if nargin<5, paramChanges=[]; end
if nargin<4, paramsName='Normal'; end
if nargin<3, levels=-10:10:100; end
if nargin==0,
%     probeFrequencies=1000;
    probeFrequencies=100:100:8000;

    lowestBF=250; 	highestBF= 8000; 	numChannels=21;
    % 21 chs (250-8k)includes BFs at 250 500 1000 2000 4000 8000
    BFlist=round(logspace(log10(lowestBF),log10(highestBF),numChannels));
else
    numChannels=length(BFlist);
end

if numChannels==1
    keyChannel=1;
else
    keyChannel=round(numChannels/2);
end

toneDuration=.2;   
rampDuration=0.002;   
silenceDuration=.02;

sampleRate=64000; dt=1/sampleRate;

%% delare 'showMap' options to control graphical output
showMapOptions.printModelParameters=0;   % prints all parameters
showMapOptions.showModelOutput=0;       % plot of all stages
showMapOptions.printFiringRates=1;      % prints stage activity levels
showMapOptions.showACF=0;               % shows SACF (probability only)
showMapOptions.showEfferent=1;          % tracks of AR and MOC
showMapOptions.surfAN=0;       % 2D plot of HSR response
showMapOptions.surfSpikes=0;            % 2D plot of spikes histogram
showMapOptions.ICrates=0;               % IC rates by CNtauGk

maxMOC=zeros(1,length(levels));

%% main computational loop (vary level)
levelNo=0;
for leveldB=levels
    levelNo=levelNo+1;
    amp=28e-6*10^(leveldB/20);
    fprintf('%4.0f\t', leveldB)

    %% generate tone and silences
    time=dt:dt:toneDuration;
    rampTime=dt:dt:rampDuration;
    ramp=[0.5*(1+cos(2*pi*rampTime/(2*rampDuration)+pi)) ...
        ones(1,length(time)-length(rampTime))];
    ramp=ramp.*fliplr(ramp);

    silence=zeros(1,round(silenceDuration/dt));

    allTones=amp*sin(2*pi*probeFrequencies'*time);
    allTones=sum(allTones,1);
    allTones= ramp.*allTones;
    inputSignal=[silence allTones];

    %% run the model
    AN_spikesOrProbability='spikes';

    MAP1_14(inputSignal, 1/dt, BFlist, ...
        paramsName, AN_spikesOrProbability, paramChanges);
    
    if length(IHCpreSynapseParams.tauCa)<2
        error('testEfferent: efferent test requires that two types of fibers are used (IHCpreSynapseParams.tauCa). Only one found')
    end
    maxMOC(levelNo)= min(MOCattenuation(keyChannel,:));
    UTIL_showMAP(showMapOptions)
    pause(0.1)

end % level

%% MOC atten/ level function
figure(21), subplot(2,1,2)
plot(levels, 20*log10(maxMOC), 'k'), hold off
title(' MOC dB attenuation'), ylabel('dB attenuation')
ylim([-30 0])
figure(21), subplot(2,1,1)
plot(levels, maxMOC, 'k'), hold off
title(' MOC attenuation (scalar)'), ylabel('attenuation (scalar)')
ylim([0 1])

set(gcf,'name','MOC atten/level')

path(restorePath)
