clear classes

% function Exp_2(isMasterNode)
% This is a template function on which other experiments can be based. If
% you are running this function with just a single node, isMasterNode
% should be set to true. On a distributed system, then set isMasterNode to
% true for the first node initialized and then to false for all subsequent
% nodes that join the party.

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up the basic experiment parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
expName = ['SS_' num2str(randi(1e5))];
if isunix
    expFolderPrefix = '/scratch/nrclark/exps/';
else
    expFolderPrefix = 'D:\Exps';
end

% expFolderPrefix = pwd;
expFolder = fullfile(expFolderPrefix,expName);
hmmFolder = fullfile(expFolder,'hmm');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sort out the training (LEARNING) condition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
learnFolder = fullfile(expFolder,'featL');

xL = jobject('L', learnFolder);

figure(11);
xL.featHaxes = gca;

figure(12);
xL.reconHaxes = gca;

figure(22);
xL.probHaxes = gca;
figure(23);
xL.probHaxesSM = gca;

figure(33);
xL.sacfHaxes = gca;
figure(34);
xL.sacfHaxesSM = gca;

% xL.participant = 'NormalNOEFF';
xL.participant = 'NormalNONE';

xL.noiseLevToUse   =  -300;
xL.speechLevToUse  = 55;

xL.MAPopHSR = 1;
xL.MAPopMSR = 0;
xL.MAPopLSR = 0;

xL.numWavs = 1; %MAx=8440

xL.useAid = 0;

xL.useSACF = 1;
xL.SACFnBins = 128;


mkdir(xL.opFolder);
xL = xL.assignFiles;
xL.wavList  = dir(fullfile(xL.wavFolder, 'MHS_2841A.wav'));

xL.removeEnergyStatic = 1;
xL.useSpectrogram = 0;
xL.numCoeff = 21;

xL.storeSelf;

xL.MAPparamChanges= {' DRNLParams.rateToAttenuationFactorProb = 0;', 'OMEParams.rateToAttenuationFactorProb=0;'};

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ** Generate features **
% This is the time consuming, processing intensive portion of the program.
% Nodes that are not the master node are only interested in the opFolder
% member of the jobjects.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
worker(xL.opFolder);

options.showEfferent=1;
% UTIL_showMAP(options)

global dt ANdt saveAN_spikesOrProbability savedBFlist saveMAPparamsName...
    savedInputSignal TMoutput OMEoutput ARattenuation ...
    DRNLoutput IHC_cilia_output IHCrestingCiliaCond IHCrestingV...
    IHCoutput ANprobRateOutput ANoutput savePavailable tauCas  ...
    CNoutput  ICoutput ICmembraneOutput ICfiberTypeRates MOCattenuation 



%% Look at the SACF

% params.acfTau  = 2;
% nBins = 128;
% 
% params.minLag  = 1 / 4000;
% params.maxLag  = 1 / 50;
% % params.lagStep = (params.maxLag - params.minLag) / 128;
% 
% params.lags = logspace(log10(params.minLag),log10(params.maxLag),nBins);
% 
% params.lambda = 10e-3;
% 
% method.dt = dt;
% method.nonlinCF = savedBFlist(1:30);
% 
% [P, BFlist, sacf, boundaryValue] = filteredSACF(ANprobRateOutput(1:30,:), method, params);
% 
% figure; imagesc(P)

% %% Current IPIH algorithm
% 
% %%%%% IPI analysis %%%%%
% iih=track_formants_from_IPI_guy(ANprobRateOutput, 1/dt);
% 
% %% %%% original VRP pattern
% niih = 50;
% [reduced_iih,ctr_freq] = mapping_IPIs_to_channels(iih,1/dt,logspace(log10(min(savedBFlist)),log10(max(savedBFlist)),niih),niih,0);
% 
% 
% subplot(2,1,1)
% imagesc(flipud(iih))
% colorbar
% subplot(2,1,2)
% imagesc(reduced_iih)
% colorbar


%% Use Tim's reduction on the SACF

% niih = 30;
% [reduced_iih,ctr_freq] = mapping_IPIs_to_channels(P,1/dt,logspace(log10(1/params.minLag),log10(1/params.maxLag),niih),niih,0);
% 
% subplot(2,1,1)
% imagesc(flipud(P))
% colorbar
% subplot(2,1,2)
% imagesc(reduced_iih)
% colorbar
