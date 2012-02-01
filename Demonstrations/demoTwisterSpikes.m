
% demoTwisterSpikes runs MAP1_14 with a .wav file for stimuus
%  This generates AN, CN and IC spikes for both HSR and LSR streams.
%  Fig. 99 shows all stages of the model.
%   spiking activity is shown in the form of raster plots
%  Sections 1-6 allow the user to change various aspects of the model
%  Remember to make a copy before changing these sections.
%
% More information can be found in the 'MAP1_14 quick reference'
%  document in the documentation folder.


dbstop if error
restorePath=path;
% access to relevant folders
addpath (['..' filesep 'MAP'],    ['..' filesep 'wavFileStore'], ...
    ['..' filesep 'utilities'])


%%  #1 parameter file name
MAPparamsName='Normal';


%% #2 spikes (slow) representation
AN_spikesOrProbability='spikes';


%% #3 speech file input
fileName='twister_44kHz';


%% #4 rms level
% signal details
leveldBSPL=70;                  % dB SPL


%% #5 number of channels in the model
%   21-channel model (log spacing)
numChannels=11;
lowestBF=250; 	highestBF= 8000;
BFlist=round(logspace(log10(lowestBF), log10(highestBF), numChannels));


%% #6 change model parameters
paramChanges={};

%% delare showMap options
showMapOptions=[];  % use defaults

% or (example: show everything including an smoothed SACF output
showMapOptions.printModelParameters=1;
showMapOptions.showModelOutput=1;
showMapOptions.printFiringRates=1;

%% Generate stimuli
[inputSignal sampleRate]=wavread(fileName);
inputSignal(:,1);
targetRMS=20e-6*10^(leveldBSPL/20);
rms=(mean(inputSignal.^2))^0.5;
amp=targetRMS/rms;
inputSignal=inputSignal*amp;


%% run the model
tic

fprintf('\n')
disp(['Signal duration= ' num2str(length(inputSignal)/sampleRate)])
disp([num2str(numChannels) ' channel model'])
disp('Computing ...')

MAP1_14(inputSignal, sampleRate, BFlist, ...
    MAPparamsName, AN_spikesOrProbability, paramChanges);


%% the model run is finished. Now display the results
% delare showMap options
UTIL_showMAP(showMapOptions)

toc
% All done. Now sweep the path!
path(restorePath)

