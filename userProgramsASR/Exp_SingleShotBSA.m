% clear classes

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

% figure(11);
% xL.featHaxes = gca;
% 
% figure(12);
% xL.reconHaxes = gca;

% figure(22);
% xL.probHaxes = gca;
figure(23);
xL.probHaxesSM = gca;

% figure(33);
% xL.sacfHaxes = gca;
% figure(34);
% xL.sacfHaxesSM = gca;

% xL.participant = 'NormalNOEFF';
xL.participant = 'NormalDIFF';

xL.noiseLevToUse   = 0;
xL.speechLevToUse  = 100;

xL.MAPopHSR = 1;
xL.MAPopMSR = 0;
xL.MAPopLSR = 0;

xL.numWavs = 1; %MAx=8440

xL.useAid = 0;

xL.useSACF = 0;
xL.SACFnBins = 128;


mkdir(xL.opFolder);
xL = xL.assignFiles;
xL.wavList  = dir(fullfile(xL.wavFolder, 'MHS_2841A.wav'));
xL.noiseName = '20TalkerBabble';

xL.removeEnergyStatic = 1;
xL.doCMN = 0;

xL.useSpectrogram = 0;
xL.numCoeff = 14;

xL.noisePreDur = 1;
xL.noisePostDur = 0.1;
xL.truncateDur  = xL.noisePreDur-0.1; %Dr. RF used 0.550




xL.MAPparamChanges=         {'DRNLParams.rateToAttenuationFactorProb=7;','DRNLParams.MOCrateThresholdProb=85;', 'DRNLParams.MOCtau=2;',...
                             'OMEParams.rateToAttenuationFactorProb=0.00;',...
                             'DRNLParams.MOCtauR=.4;', 'DRNLParams.MOCtauF=.1;'};% xL.MAPparamChanges= {'OMEParams.rateToAttenuationFactorProb=0;', 'DRNLParams.rateToAttenuationFactorProb = 0.010;', 'DRNLParams.MOCrateThresholdProb =40;','DRNLParams.MOCtau =0.35;'};


xL.storeSelf;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ** Generate features **
% This is the time consuming, processing intensive portion of the program.
% Nodes that are not the master node are only interested in the opFolder
% member of the jobjects.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


options.showEfferent=1;
% UTIL_showMAP(options)

global dt ANdt saveAN_spikesOrProbability savedBFlist saveMAPparamsName...
    savedInputSignal TMoutput OMEoutput ARattenuation ...
    DRNLoutput IHC_cilia_output IHCrestingCiliaCond IHCrestingV...
    IHCoutput ANprobRateOutput ANoutput savePavailable tauCas  ...
    CNoutput  ICoutput ICmembraneOutput ICfiberTypeRates MOCattenuation 

worker(xL.opFolder);

%%
% clims = ;

% figure(12); set(gca, 'CLim', [30 80] ); colorbar
figure(23); set(gca, 'CLim', [60 160] )

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

%%
figure(23)
set(gcf, 'Position', [913-200   544   553   169])

%%
figure(99)
% dt=ANdt
startP = round(xL.truncateDur*(1/dt));
endP = size(MOCattenuation,2);

durS =endP-startP+1;
tAxis = dt:dt:durS*dt;

set(gcf, 'Position', [913   544   553   169])
% imagesc(-20*log10(flipud(MOCattenuation(:, startP:end))), [0 25])




myBFlist = savedBFlist;
YTickIdx = 1:floor(numel(myBFlist)/6):numel(myBFlist);
YTickIdxRev = numel(myBFlist)+1-YTickIdx;

MOClim = MOCattenuation(:, startP:end);
imagesc(tAxis, [],-20*log10(flipud(MOClim)), [0 25])
set(gca, 'YTick', YTickIdx);
set(gca, 'YTickLabel', num2str(    myBFlist(YTickIdxRev)'     ));
ylabel('cf (Hz)')
xlabel('Time (s)')
xlim([0 durS*dt])
colorbar

%%
figure(101)
set(gcf, 'Position', [513   544   553   169])

AR2D = repmat(ARattenuation(startP:end), size(MOClim,1),1);


imagesc(tAxis, [],-20*log10(flipud(AR2D)), [0 25])
set(gca, 'YTick', YTickIdx);
set(gca, 'YTickLabel', num2str(    myBFlist(YTickIdxRev)'     ));
ylabel('cf (Hz)')
xlabel('Time (s)')
xlim([0 durS*dt])
colorbar
