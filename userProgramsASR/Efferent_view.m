% function Efferent_view
%This script is a function now, just so I dont f**k up my path

obj.MAProot = fullfile('..');
addpath(...fullfile(obj.MAProot, 'modules'),...
                fullfile(obj.MAProot, 'utilities'),...
                fullfile(obj.MAProot, 'MAP'),...
                fullfile(obj.MAProot, 'parameterStore'),...
                fullfile('ASR files'));

close all; clear all; clc;

sr = 44100;
dt = 1/sr;
dur = 0.5;
freq = 1000;

nn=0;
for levelSPL = 60;
nn = nn+1;

tAxis = dt:dt:dur;

ipSig = sin(2*pi*freq*tAxis);
% ipSig = randn(size(tAxis));
% ipSig = wavread(fullfile('demo_wavs','noises', 'pink.wav'));
% ipSig = ipSig(1:numel(tAxis));

ipSig = ipSig./sqrt(mean(ipSig.^2));
ipSig = ipSig * 20e-6 * 10 ^ (levelSPL/20);

ipSig = [zeros(1, numel(ipSig)/2) ipSig zeros(1, numel(ipSig)/2)];

paramChanges= { ...
                'OMEParams.rateToAttenuationFactorProb=0;',...
                'OMEParams.ARrateThreshold=30;',... %Threshold of 40 makes AR kick off around 65 dB for bb noise
                'OMEParams.ARtau=0.1;',...
                'DRNLParams.MOCtauR=0.400;',...
                'DRNLParams.MOCtauF=0.100;',...%DRNLParams.MOCtauR;',...                
                'DRNLParams.rateToAttenuationFactorProb=9;',...
                'DRNLParams.MOCrateThresholdProb=85;',...
                };

AN_spikesOrProbability = 'probability';
MAP1_14(ipSig, sr, -1, 'NormalDIFF', AN_spikesOrProbability, paramChanges)

options.showEfferent=1;
UTIL_showMAP(options)
drawnow

%%
global MOCattenuation ANprobRateOutput SAVEsmoothedRates
size(MOCattenuation);

% attFraction = sqrt(mean((MOCattenuation.^2),2));
attdB(nn) = min( mean(20*log10(MOCattenuation(:, ceil(numel(tAxis)/2):end )), 2) )
% attdB(nn) = min( 20*log10(MOCattenuation(:)) )

figure; plot(ANprobRateOutput(31:end)')
figure; plot(SAVEsmoothedRates')
figure; plot(20*log10(MOCattenuation'))

end


