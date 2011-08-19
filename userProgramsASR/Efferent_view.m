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
dur = 0.2;
freq = 1000;

nn=0;
for levelSPL = 10:10:90;
nn = nn+1;

tAxis = dt:dt:dur;

ipSig = sin(2*pi*freq*tAxis);
% ipSig = randn(size(tAxis));

ipSig = ipSig./sqrt(mean(ipSig.^2));
ipSig = ipSig * 20e-6 * 10 ^ (levelSPL/20);

paramChanges = {};
paramChanges{numel(paramChanges)+1} = 'DRNLParams.rateToAttenuationFactorProb = -0.08;';%DEFAULT = 0.005;  % strength of MOC
paramChanges{numel(paramChanges)+1} = 'DRNLParams.MOCrateThresholdProb = 40;';%DEFAULT = 70;
paramChanges{numel(paramChanges)+1} = 'DRNLParams.MOCtau =0.05;' ;%DEFAULT = 20k in new params file
% paramChanges{numel(paramChanges)+1} = 'DRNLParams.a = 1250;';

paramChanges{numel(paramChanges)+1} = 'OMEParams.rateToAttenuationFactorProb = 0.00;';%DEFAULT = 0.01;
paramChanges{numel(paramChanges)+1} = 'OMEParams.ARtau=.06;';
paramChanges{numel(paramChanges)+1} = 'OMEParams.ARrateThreshold=20;';

AN_spikesOrProbability = 'probability';
MAP1_14(ipSig, sr, -1, 'NormalNICK', AN_spikesOrProbability, paramChanges)

options.showEfferent=1;
UTIL_showMAP(options)
drawnow

%%
global MOCattenuation ANprobRateOutput
size(MOCattenuation);

% attFraction = sqrt(mean((MOCattenuation.^2),2));
attdB(nn) = min( mean(20*log10(MOCattenuation(:, ceil(numel(tAxis)/2):end )), 2) )
% attdB(nn) = min( 20*log10(MOCattenuation(:)) )

end


