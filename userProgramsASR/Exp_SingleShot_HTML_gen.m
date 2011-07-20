% function Exp_2(isMasterNode)
% This is a template function on which other experiments can be based. If
% you are running this function with just a single node, isMasterNode
% should be set to true. On a distributed system, then set isMasterNode to
% true for the first node initialized and then to false for all subsequent
% nodes that join the party.

aArr = [ 20e3; 10000; 5000; 2500; 1250; 675; 0];
spLevs = [30 60 90];

for kk = 1:numel(spLevs)    
for nn = 1:numel(aArr)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sort out the training (LEARNING) condition
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    learnFolder = fullfile(expFolder,'featL');
    
    xL = jobject('L', learnFolder);
    
    figure(22);
    xL.featHaxes = gca;
    
    % xL.participant = 'NormalNOEFF';
    xL.participant = 'NormalNONE';
    
    xL.noiseLevToUse   =  -300;
    xL.speechLevToUse  =  spLevs(kk);
    
    xL.MAPopHSR = 1;
    xL.MAPopMSR = 0;
    xL.MAPopLSR = 0;
    
    xL.numWavs = 1; %MAx=8440
    
    xL.useAid = 0;
    
    % if isMasterNode
    mkdir(xL.opFolder);
    xL = xL.assignFiles;
    xL.wavList  = dir(fullfile(xL.wavFolder, 'MHS_2841A.wav'));
    
    xL.removeEnergyStatic = 0;
    xL.useSpectrogram = 0;
    xL.numCoeff =9;
    
  
    % end
    
    xL.MAPparamChanges= ...%{'DRNLParams.a=0'};
        {'DRNLParams.rateToAttenuationFactorProb = 0;',...
         'OMEParams.rateToAttenuationFactorProb=0;',...
        ['DRNLParams.a=' num2str(aArr(nn)) ';']};
    


  xL.storeSelf;
    
  %%
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    
   
    disp(['Current a value = ' num2str(aArr(nn))])
    drawnow
    
end
end

