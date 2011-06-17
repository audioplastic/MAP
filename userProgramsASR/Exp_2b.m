function Exp_2b(isMasterNode)
% This is a template function on which other experiments can be based. If
% you are running this function with just a single node, isMasterNode
% should be set to true. On a distributed system, then set isMasterNode to
% true for the first node initialized and then to false for all subsequent
% nodes that join the party.

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up the basic experiment parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
expName = 'JunTest_LOUD_0';
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

xL.participant = 'NormalNOEFF';
% xL.participant = 'Normal';

xL.noiseLevToUse   =  50;
xL.speechLevToUse  =  80;

xL.MAPopHSR = 1;
xL.MAPopMSR = 0;
xL.MAPopLSR = 0;
xL.MAPuseEfferent = 0;
xL.numWavs = 1000; %MAx=8440

if isMasterNode
    mkdir(xL.opFolder);
    xL = xL.assignFiles;
    xL.storeSelf;
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sort out the testing (RECOGNITION) conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nzLevel = 50:5:100;
xR=cell(size(nzLevel));
recConditions = numel(nzLevel);
for nn = 1:recConditions    
    xR{nn} = xL; %simply copy the "Learn" object and change it a bit below
    recFolder = fullfile(expFolder,['featR_nz' num2str(nzLevel(nn))]);
    xR{nn}.opFolder = recFolder;    
    
    %These are the interesting differences between training and testing
    xR{nn}.numWavs = 358; %MAX = 358
    xR{nn}.noiseLevToUse = nzLevel(nn);
    
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
for nn = 1:recConditions
    worker(xR{nn}.opFolder);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Train and test the recogniser - a job for the master node only
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isMasterNode    
    y = HMMclass(hmmFolder);    
    y.createSCP(xL.opFolder)
    y.createMLF(xL.opFolder)
    y.train(xL.opFolder) %This node can be busy training, even if other jobs are being processed
    
    xR{end}.lockJobList;
    xR{end} = xR{end}.loadSelf; %Reload changes
    xR{end}.unlockJobList;
    while(~all(xR{end}.todoStatus==2))        
        disp('Waiting on straggler nodes to complete their jobs before HMM is tested . . .')
        pause(2);
        xR{end}.lockJobList;
        xR{end} = xR{end}.loadSelf; %Reload incase changed
        xR{end}.unlockJobList;
    end
      
    for nn = 1:recConditions
        y.createSCP(xR{nn}.opFolder);
        y.test(xR{nn}.opFolder);
    end
    
    %Show all of the scores in the command window at the end
    for nn = 1:recConditions
        y.score(xR{nn}.opFolder);
    end
end
