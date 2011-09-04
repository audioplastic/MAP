% function Efferent_view
% This script is a function now, just so I dont f**k up my path

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
freq = 520;

nn=0;
% for levelSPL = 40:10:70;
for levelSPL = [50:10:100]
% levelSPL = 50;
nn = nn+1;
levelRec(nn) = levelSPL;

tAxis = dt:dt:dur;

ipSig = sin(2*pi*freq*tAxis);

% pink noise code
% ipSig = randn(size(tAxis));
% num_taps = 1024;
% a = zeros(1,num_taps);
% a(1) = 1;
% for ii = 2:num_taps
%     a(ii) = (ii - 2.5) * a(ii-1) / (ii-1);
% end
% ipSig = filter(1,a,ipSig);
% end of pink code

ipSig = wavread(fullfile('demo_wavs','noises', 'pink.wav'));
% ipSig = wavread(fullfile('D:\ASRexperiments\Stimuli\noises', '20TalkerBabble.wav'));
ipSig = ipSig(1:numel(tAxis));

% soundsc(ipSig,sr)

ipSig = ipSig./sqrt(mean(ipSig.^2));
ipSig = ipSig * 20e-6 * 10 ^ (levelSPL/20);

paramChanges = {};
paramChanges{numel(paramChanges)+1} = 'DRNLParams.rateToAttenuationFactorProb =  7;';%GOOD = 0.012  %DEFAULT = 0.005;  % strength of MOC
paramChanges{numel(paramChanges)+1} = 'DRNLParams.MOCrateThresholdProb = 85;';%GOOD=140 %DEFAULT = 70;
paramChanges{numel(paramChanges)+1} = 'DRNLParams.MOCtau = 2;'; %DEFAULT = 0.1;

paramChanges{numel(paramChanges)+1} = 'OMEParams.rateToAttenuationFactorProb = 20;';
paramChanges{numel(paramChanges)+1} = 'OMEParams.ARrateThreshold = 50;';
paramChanges{numel(paramChanges)+1} = 'OMEParams.ARtau=0.1;';



% paramChanges{numel(paramChanges)+1} = 'OMEParams.rateToAttenuationFactorProb = 0.05;';%DEFAULT = 0.01;
% paramChanges{numel(paramChanges)+1} = 'OMEParams.rateToAttenuationFactor = 0;';%DEFAULT = 0.01;
% paramChanges{numel(paramChanges)+1} = 'DRNLParams.a=1e4;'; %DEFAULT = 5e4;



AN_spikesOrProbability = 'probability';
MAP1_14(ipSig, sr, [ freq ], 'NormalDIFF', AN_spikesOrProbability, paramChanges)

options.showEfferent=1;
UTIL_showMAP(options)
drawnow

%%
global MOCattenuation ARattenuation
size(MOCattenuation);

% attFraction = sqrt(mean((MOCattenuation.^2),2));
% attdB(nn) = -min( mean(20*log10(MOCattenuation), 2) )
attdB(nn) = -min( mean(20*log10(MOCattenuation(:, ceil(numel(tAxis)/2):end )), 2) )
ARattdB(nn) = -mean(20*log10(ARattenuation( ceil(numel(tAxis)/2):end ))) 
% end

%%
global ANprobRateOutput ANoutput
% if strcmpi(AN_spikesOrProbability, 'spikes')
% %     ANprobRateOutput = ANoutput;
%     ANprobRateOutput = jobject.makeANsmooth(ANoutput, sr)*sr;
% end

nChans = size(ANprobRateOutput,1)/2;
nSamples = size(ANprobRateOutput,2);
LSR = ANprobRateOutput(1:nChans, : );
HSR = ANprobRateOutput((1+nChans):end,  : );
figure(77);
subplot(2,1,1); imagesc(flipud(LSR)); title('LSR'); colorbar; drawnow
subplot(2,1,2); imagesc(flipud(HSR)); title('HSR'); colorbar; drawnow

rateLSR(nn) = max(mean(LSR(:,ceil(0.75*nSamples):end),2)) %Last 75% of stimulus for sustained response
rateHSR(nn) = max(mean(HSR(:,ceil(0.75*nSamples):end),2))
end

%%
figure(33)
subplot(3,1,1); plot(levelRec, attdB); title('MOC attenuation')
subplot(3,1,2); plot(levelRec, rateLSR); title('LSR rate')
subplot(3,1,3); plot(levelRec, rateHSR); title('HSR rate')

figure(44)
plot(levelRec, attdB, levelRec, ARattdB); title('MOC/AR attenuation'); ylim([0 45])

