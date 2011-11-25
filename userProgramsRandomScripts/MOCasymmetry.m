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
            
global ANprobRateOutput  ANdt savedBFlist MOCattenuation DRNLParams

paramChanges= { 'OMEParams.rateToAttenuationFactorProb=0;',...
                'DRNLParams.MOCtauR=0.250;',...
                'DRNLParams.MOCtauF=0.250;',...                
                'DRNLParams.rateToAttenuationFactorProb=0;',...
                'DRNLParams.MOCrateThresholdProb=100;'};
MAP1_14(stimulus, sr, BFlist, participant, 'probability', paramChanges);



%%
HSRprob = ANprobRateOutput(2,:);
tAxisAN = dtL:dtL:dtL*numel(HSRprob);
subplot(4,1,1); plot(tAxisAN, HSRprob); ylabel('Firing prob (sp/s)')
xlim([tAxisAN(1) tAxisAN(end)])

%%
MOCfilt_aR = dtL/DRNLParams.MOCtauR - 1; % R = Rising edge
MOCfilt_bR = 1.0 + MOCfilt_aR;
MOCfilt_aF = dtL/DRNLParams.MOCtauF - 1; % F = Falling edge
MOCfilt_bF = 1.0 + MOCfilt_aF;
% smoothedRates = filter(MOCfilt_bR, -MOCfilt_aR, HSRprob);
smoothedRates = zeros(size(HSRprob));
rates = HSRprob;
MOCprobBoundary{1}=0;
for idx=1%1chan
    for nn = 1:numel(HSRprob)
        if rates(idx,nn) < MOCprobBoundary{idx} %// - This is line to make smoothing only apply to release
            smoothedRates(nn) = MOCfilt_bF*rates(idx,nn) - MOCfilt_aF*MOCprobBoundary{idx};% // difference eqn for one-pole lpf
        else
            smoothedRates(nn) = MOCfilt_bR*rates(idx,nn) - MOCfilt_aR*MOCprobBoundary{idx};% // difference eqn for one-pole lpf
        end
        MOCprobBoundary{idx} = smoothedRates(nn);
    end
end

%%
demoFactor = 7;
x = -20*log10(  max(smoothedRates/DRNLParams.MOCrateThresholdProb,1)  )*demoFactor; %dB attenuation
%x = 10.^(x/20);
subplot(4,1,2); plot(tAxisAN, smoothedRates); ylabel('Firing prob (sp/s)')
xlim([tAxisAN(1) tAxisAN(end)])


subplot(4,1,3); plot(tAxisAN, max(smoothedRates,DRNLParams.MOCrateThresholdProb)); ylabel('Firing prob (sp/s)')
xlim([tAxisAN(1) tAxisAN(end)])


subplot(4,1,4); plot(tAxisAN, -x); ylabel('Firing prob (sp/s)')
xlim([tAxisAN(1) tAxisAN(end)])

%%
% subplot(3,1,3); plot(tAxisAN, -20*log10( MOCattenuation)); ylabel('Actual MOC attn (dB)')
% xlim([tAxisAN(1) tAxisAN(end)])
% xlabel('time(s)')
% set(gcf,'Position', [680   523   560   575])

