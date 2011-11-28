close all; clear all; clc

% sr = 44.1e3;
% fZero = 1000;
% pipDuration = 1;
% silDuration = 2;
% pipLvl = 60; %dB SPL
% 
% dtL = 1/sr;
% tAxis = dtL:dtL:pipDuration;
% stimulus = [zeros(1, ceil(silDuration*sr)),...
%             10^(pipLvl/20)*20e-6*sin(2*pi*fZero*tAxis),...
%             zeros(1, ceil(silDuration*sr))];

%%
speechLev = 100;
noiseLev = speechLev-10;
[stimulus,sr] = wavread(fullfile('..','userProgramsASR','demo_wavs','microAURORA','MHS_2841A'));
% stimulus = randn(44100,1); sr = 44.1e3;
stimulus = 20e-6*10^(speechLev/20)  *  stimulus./(sqrt(mean(stimulus.^2)));
stimulus = [zeros(44100,1); stimulus];
dtL = 1/sr;

%%
participant = 'NormalDIFF';
addpath(...
        fullfile('..', 'utilities'),...
        fullfile('..', 'MAP'),...
        fullfile('..', 'parameterStore'));
            
global ANprobRateOutput  ANdt savedBFlist MOCattenuation DRNLParams ARattenuation

paramChanges= { 'OMEParams.rateToAttenuationFactorProb=12;',...
                'OMEParams.ARrateThreshold=45;',... %Threshold of 40 makes AR kick off around 65 dB for bb noise
                'DRNLParams.MOCtauR=2;',...
                'DRNLParams.MOCtauF=2;',...                
                'DRNLParams.rateToAttenuationFactorProb=9;',...
                'DRNLParams.MOCrateThresholdProb=85;'};
MAP1_14(stimulus, sr, -1, participant, 'probability', paramChanges);



%%
HSRprob = ANprobRateOutput(numel(savedBFlist)+1:end,:);
LSRprob = ANprobRateOutput(1:numel(savedBFlist),:);
tAxisAN = dtL:dtL:dtL*numel(HSRprob);

YTickIdx = 1:floor(numel(savedBFlist)/6):numel(savedBFlist);
YTickIdxRev = numel(savedBFlist)+1-YTickIdx;
subplot(3,1,1); imagesc(HSRprob); colorbar('peer', gca)
set(gca, 'YTick', YTickIdx);
set(gca, 'YTickLabel', num2str(    savedBFlist(YTickIdxRev)'     ));
ylabel('cf in Hz')

subplot(3,1,2); imagesc(20*log10(MOCattenuation)); colorbar('peer', gca)
set(gca, 'YTick', YTickIdx);
set(gca, 'YTickLabel', num2str(    savedBFlist(YTickIdxRev)'     ));
ylabel('cf in Hz')

subplot(3,1,3); plot(20*log10(ARattenuation));colorbar
