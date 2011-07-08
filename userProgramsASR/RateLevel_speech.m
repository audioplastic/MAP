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
dur = 0.25;
freq = 1000;

nn=0;
% for levelSPL = 0:10:100;
levelSPL = 60;
SNR = 1000;
preDur = 1;

nn = nn+1;
levelRec(nn) = levelSPL;

tAxis = dt:dt:dur;

% ipSig = sin(2*pi*freq*tAxis);
% ipSig = randn(size(tAxis));
[ipSig, sr] = wavread('C:\corpora\AURORA digits (wav)\TrainingData-Clean\MHS_2841A.wav');

ipSig = ipSig./sqrt(mean(ipSig.^2));
ipSig = ipSig * 20e-6 * 10 ^ (levelSPL/20);

preS = ceil(preDur*sr);
ipSig = [zeros(preS, 1);  ipSig];

[nzSig, ~] = wavread('C:\corpora\noises\pink.wav');
startIdx = randi(numel(nzSig)-numel(ipSig));
nzSig = nzSig(startIdx:startIdx+numel(ipSig)-1);

nzSig = nzSig./sqrt(mean(nzSig.^2));
nzSig = nzSig * 20e-6 * 10 ^ ((levelSPL-SNR)/20);

ipSig = ipSig+nzSig;

paramChanges = {};
paramChanges{numel(paramChanges)+1} = 'DRNLParams.rateToAttenuationFactorProb = 0;';%FIX = -10^(-10/20); %GOOD = 0.012  %DEFAULT = 0.005;  % strength of MOC
paramChanges{numel(paramChanges)+1} = 'DRNLParams.MOCrateThresholdProb = 60;';%GOOD=140 %DEFAULT = 70;

paramChanges{numel(paramChanges)+1} = 'OMEParams.rateToAttenuationFactorProb = 0;';%DEFAULT = 0.01;
paramChanges{numel(paramChanges)+1} = 'DRNLParams.a=1e4;'; %DEFAULT = 5e4;
paramChanges{numel(paramChanges)+1} = 'DRNLParams.MOCtau = 0.45;'; %DEFAULT = 0.1;

AN_spikesOrProbability = 'probability';
BFlist=round(logspace(log10(100),log10(4500),30));
MAP1_14(ipSig, sr, BFlist, 'Normal', AN_spikesOrProbability, paramChanges)

options.showEfferent=1;
UTIL_showMAP(options)
drawnow

%%
global MOCattenuation
size(MOCattenuation);

% attFraction = sqrt(mean((MOCattenuation.^2),2));
attdB(nn) = -min( mean(20*log10(MOCattenuation), 2) )

% end

%%
global ANprobRateOutput
nChans = size(ANprobRateOutput,1)/2;
nSamples = size(ANprobRateOutput,2);
LSR = ANprobRateOutput(1:nChans, preS:end );
HSR = ANprobRateOutput((1+nChans):end,  preS:end );

HSR = jobject.makeANsmooth(HSR, sr);


figure(77);
subplot(2,1,1); imagesc(flipud(LSR)); title('LSR'); colorbar; drawnow; colormap bone
subplot(2,1,2); imagesc(flipud(HSR), [50 120]); title('HSR'); colorbar; drawnow; colormap bone; 

rateLSR(nn) = max(mean(LSR(:,ceil(0.25*nSamples):end),2)) %Last 75% of stimulus for sustained response
rateHSR(nn) = max(mean(HSR(:,ceil(0.25*nSamples):end),2))
% end

%%
% figure(33)
% subplot(3,1,1); plot(levelRec, attdB); title('MOC attenuation')
% subplot(3,1,2); plot(levelRec, rateLSR); title('LSR rate')
% subplot(3,1,3); plot(levelRec, rateHSR); title('HSR rate')

