function Efferent_view
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
dur =1;
freq = 1000;
levelSPL = 65;

tAxis = dt:dt:dur;

ipSig = sin(2*pi*freq*tAxis);
% ipSig = randn(size(tAxis));

ipSig = ipSig./sqrt(mean(ipSig.^2));
ipSig = ipSig * 20e-6 * 10 ^ (levelSPL/20);

paramChanges = {};
paramChanges{numel(paramChanges)+1} = 'DRNLParams.rateToAttenuationFactorProb = 0;';%DEFAULT = 0.005;  % strength of MOC
paramChanges{numel(paramChanges)+1} = 'OMEParams.rateToAttenuationFactorProb = 0.01;';%DEFAULT = 0.01;

AN_spikesOrProbability = 'probability';
MAP1_14(ipSig, sr, -1, 'Normal', AN_spikesOrProbability, paramChanges)

options.showEfferent=1;
UTIL_showMAP(options)
