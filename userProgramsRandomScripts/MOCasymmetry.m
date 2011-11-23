close all; clear all; clc

sr = 44.1e3;
fZero = 1000;
pipDuration = 1;
silDuration = 2;
pipLvl = 60; %dB SPL

dtL = 1/sr;
tAxis = dtL:dtL:pipDuration;
stimulus = [zeros(1, ceil(silDuration*sr)),...
            10^(pipLvl/20)*20e-6*sin(2*pi*fZero*tAxis),...
            zeros(1, ceil(silDuration*sr))];

BFlist = fZero;
participant = 'NormalDIFF';
addpath(...
        fullfile('..', 'utilities'),...
        fullfile('..', 'MAP'),...
        fullfile('..', 'parameterStore'));
            
global ANprobRateOutput  ANdt savedBFlist MOCattenuation

paramChanges= { 'DRNLParams.MOCtauR=0.25;',...
                'DRNLParams.MOCtauF=0.25;',...
                'OMEParams.rateToAttenuationFactorProb=0;',...
                'DRNLParams.rateToAttenuationFactorProb=7;',...
                'DRNLParams.MOCrateThresholdProb=110;'};
MAP1_14(stimulus, sr, BFlist, participant, 'probability', paramChanges);

%%
HSRprob = ANprobRateOutput(2,:);
tAxisAN = dtL:dtL:dtL*numel(HSRprob);
subplot(2,1,1); plot(tAxisAN, HSRprob); ylabel('Firing prob (sp/s)')
xlim([tAxisAN(1) tAxisAN(end)])
subplot(2,1,2); plot(tAxisAN, -20*log10( MOCattenuation)); ylabel('MOC attn (dB)')
xlim([tAxisAN(1) tAxisAN(end)])
xlabel('time(s)')
set(gcf,'Position', [680   523   560   575])