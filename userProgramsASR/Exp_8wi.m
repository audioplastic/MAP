function Exp_8wi(isMasterNode)
% This experiment looks at SRTs.

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up the basic experiment parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
expName = '8wib';
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

xL.participant = 'Normal';%'NormalDIFF';
xL.MAPparamChanges= {' DRNLParams.rateToAttenuationFactorProb = 0;', 'OMEParams.rateToAttenuationFactorProb=0;'};

xL.noiseLevToUse   =  -200;
xL.speechLevToUse  =  60;

xL.MAPopHSR = 1;
xL.MAPopMSR = 0;
xL.MAPopLSR = 0;

xL.removeEnergyStatic = false;
if xL.removeEnergyStatic == true
    xL.numCoeff = 10;
else
    xL.numCoeff = 9;
end

xL.numWavs = 8440; %MAx=8440

xL.noisePreDur = 1;
xL.noisePostDur = 0.1;
xL.truncateDur  = 0.9; %Dr. RF used 0.550

xL.noiseName = 'pink';

if isMasterNode
    mkdir(xL.opFolder);
    xL = xL.assignFiles;
    xL.storeSelf;
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sort out the testing (RECOGNITION) conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spLevel = 40:10:100;
xR=cell(size(spLevel));
recConditions = numel(spLevel);
for nn = 1:recConditions    
    xR{nn} = xL; %simply copy the "Learn" object and change it a bit below
    recFolder = fullfile(expFolder,['NONEfeatR_sp' num2str(spLevel(nn))]);
    xR{nn}.opFolder = recFolder;    
    
    %These are the interesting differences between training and testing
    xR{nn}.numWavs = 358; %MAX = 358
    xR{nn}.speechLevToUse = spLevel(nn);    
    
            
    %Now just to wrap it up ready for processing
    if isMasterNode
        mkdir(xR{nn}.opFolder);
        xR{nn} = xR{nn}.assignWavPaths('R');
        xR{nn} = xR{nn}.assignFiles;
        xR{nn}.storeSelf;
    end
end


tmpIdx=0;
for nn = recConditions+1:2*recConditions   
    tmpIdx=tmpIdx+1;
    xR{nn} = xL; %simply copy the "Learn" object and change it a bit below
    recFolder = fullfile(expFolder,['OHCfeatR_sp' num2str(spLevel(tmpIdx))]);
    xR{nn}.opFolder = recFolder;    
    
    %These are the interesting differences between training and testing
    xR{nn}.numWavs = 358; %MAX = 358
    xR{nn}.speechLevToUse = spLevel(tmpIdx);
    xR{nn}.MAPparamChanges= {'DRNLParams.a=10000;'};
    xR{nn}.useAid = 0; %redundant but here to be explicit
    
    %Now just to wrap it up ready for processing
    if isMasterNode
        mkdir(xR{nn}.opFolder);
        xR{nn} = xR{nn}.assignWavPaths('R');
        xR{nn} = xR{nn}.assignFiles;
        xR{nn}.storeSelf;
    end
end

tmpIdx=0;
for nn = 2*recConditions+1:3*recConditions   
    tmpIdx=tmpIdx+1;
    xR{nn} = xL; %simply copy the "Learn" object and change it a bit below
    recFolder = fullfile(expFolder,['OHCfeatR_sp' num2str(spLevel(tmpIdx))]);
    xR{nn}.opFolder = recFolder;    
    
    %These are the interesting differences between training and testing
    xR{nn}.numWavs = 358; %MAX = 358
    xR{nn}.speechLevToUse = spLevel(tmpIdx);
    xR{nn}.MAPparamChanges= {'DRNLParams.a=5000;'};
    xR{nn}.useAid = 0; %redundant but here to be explicit
    
    %Now just to wrap it up ready for processing
    if isMasterNode
        mkdir(xR{nn}.opFolder);
        xR{nn} = xR{nn}.assignWavPaths('R');
        xR{nn} = xR{nn}.assignFiles;
        xR{nn}.storeSelf;
    end
end

tmpIdx=0;
for nn = 3*recConditions+1:4*recConditions   
    tmpIdx=tmpIdx+1;
    xR{nn} = xL; %simply copy the "Learn" object and change it a bit below
    recFolder = fullfile(expFolder,['OHCfeatR_sp' num2str(spLevel(tmpIdx))]);
    xR{nn}.opFolder = recFolder;    
    
    %These are the interesting differences between training and testing
    xR{nn}.numWavs = 358; %MAX = 358
    xR{nn}.speechLevToUse = spLevel(tmpIdx);
    xR{nn}.MAPparamChanges= {'DRNLParams.a=2500;'};
    xR{nn}.useAid = 0; %redundant but here to be explicit
    
    %Now just to wrap it up ready for processing
    if isMasterNode
        mkdir(xR{nn}.opFolder);
        xR{nn} = xR{nn}.assignWavPaths('R');
        xR{nn} = xR{nn}.assignFiles;
        xR{nn}.storeSelf;
    end
end

tmpIdx=0;
for nn = 4*recConditions+1:5*recConditions   
    tmpIdx=tmpIdx+1;
    xR{nn} = xL; %simply copy the "Learn" object and change it a bit below
    recFolder = fullfile(expFolder,['OHCfeatR_sp' num2str(spLevel(tmpIdx))]);
    xR{nn}.opFolder = recFolder;    
    
    %These are the interesting differences between training and testing
    xR{nn}.numWavs = 358; %MAX = 358
    xR{nn}.speechLevToUse = spLevel(tmpIdx);
    xR{nn}.MAPparamChanges= {'DRNLParams.a=1250;'};
    xR{nn}.useAid = 0; %redundant but here to be explicit
    
    %Now just to wrap it up ready for processing
    if isMasterNode
        mkdir(xR{nn}.opFolder);
        xR{nn} = xR{nn}.assignWavPaths('R');
        xR{nn} = xR{nn}.assignFiles;
        xR{nn}.storeSelf;
    end
end

tmpIdx=0;
for nn = 5*recConditions+1:6*recConditions   
    tmpIdx=tmpIdx+1;
    xR{nn} = xL; %simply copy the "Learn" object and change it a bit below
    recFolder = fullfile(expFolder,['OHCfeatR_sp' num2str(spLevel(tmpIdx))]);
    xR{nn}.opFolder = recFolder;    
    
    %These are the interesting differences between training and testing
    xR{nn}.numWavs = 358; %MAX = 358
    xR{nn}.speechLevToUse = spLevel(tmpIdx);
    xR{nn}.MAPparamChanges= {'DRNLParams.a=750;'};
    xR{nn}.useAid = 0; %redundant but here to be explicit
    
    %Now just to wrap it up ready for processing
    if isMasterNode
        mkdir(xR{nn}.opFolder);
        xR{nn} = xR{nn}.assignWavPaths('R');
        xR{nn} = xR{nn}.assignFiles;
        xR{nn}.storeSelf;
    end
end

% you would normally end here and use a separate script for workers

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ** Generate features **
% This is the time consuming, processing intensive portion of the program.
% Nodes that are not the master node are only interested in the opFolder
% member of the jobjects.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
worker(xL.opFolder);

if ~isMasterNode %dont bother wasting master node effort on generating testing features (for now)
    for nn = 1:6*recConditions
        worker(xR{nn}.opFolder);
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Train and test the recogniser - a job for the master node only
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isMasterNode    
    while(~all(xL.todoStatus==2))        
        disp('Waiting on straggler nodes to complete their jobs before HMM is trained . . .')
        pause(30); %Wait for 30 seconds before looking again
        xL.lockJobList;
        xL = xL.loadSelf; %Reload incase changed
        xL.unlockJobList;
    end
    y = HMMclass(hmmFolder);        
    y.createSCP(xL.opFolder)
    y.createMLF(xL.opFolder)
    y.train(xL.opFolder) %This node can be busy training, even if other jobs are being processed for testing
    
    % ALLOW MASTER NODE TO MUCK IN WITH GENERATING TESTING FEATURES ONCE
    % HMM HAS BEEN TRAINED
    for nn = 1:6*recConditions
        worker(xR{nn}.opFolder);
    end    
    
    xR{end}.lockJobList;
    xR{end} = xR{end}.loadSelf; %Reload changes
    xR{end}.unlockJobList;
    while(~all(xR{end}.todoStatus==2))        
        disp('Waiting on straggler nodes to complete their jobs before HMM is tested . . .')
        pause(30); %Wait for 30 seconds before looking again
        xR{end}.lockJobList;
        xR{end} = xR{end}.loadSelf; %Reload incase changed
        xR{end}.unlockJobList;
    end
      
    for nn = 1:6*recConditions
        y.createSCP(xR{nn}.opFolder);
        y.test(xR{nn}.opFolder);
    end
    
    %Show all of the scores in the command window at the end
    for nn = 1:6*recConditions
        y.score(xR{nn}.opFolder);
    end
end
