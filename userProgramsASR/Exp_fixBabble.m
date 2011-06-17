function Exp_fixBabble(isMasterNode, whichGeezer, useEff)
% This is a template function on which other experiments can be based. If
% you are running this function with just a single node, isMasterNode
% should be set to true. On a distributed system, then set isMasterNode to
% true for the first node initialized and then to false for all subsequent
% nodes that join the party.

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up the basic experiment parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
expName = ['fixBabble_nz10R_eff', num2str(useEff), '_', whichGeezer];
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

xL.participant = whichGeezer;

xL.noiseLevToUse   =  30;
xL.speechLevToUse  =  55;
xL.speechDist = 'Uniform';
xL.speechLevStd = 30/sqrt(12); %Range of +/- 15 dB


xL.MAPopHSR = 1;
xL.MAPopMSR = 0;
xL.MAPopLSR = 0;
xL.MAPuseEfferent = useEff;
xL.numWavs = 1000;%8440;

if isMasterNode
    mkdir(xL.opFolder);
    xL = xL.assignFiles;
    xL.storeSelf;
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sort out the testing (RECOGNITION) conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spLevel = 0:5:40;
xR=cell(size(spLevel));
recConditions = numel(spLevel);
for nn = 1:recConditions    
    xR{nn} = xL; %simply copy the "Learn" object and change it a bit below
    recFolder = fullfile(expFolder,['featR_sp' num2str(spLevel(nn))]);
    xR{nn}.opFolder = recFolder;    
    
    %These are the interesting differences between training and testing
    xR{nn}.numWavs = 100;%358;
    xR{nn}.noiseLevToUse = 10;
    xR{nn}.speechLevToUse = spLevel(nn);
    xR{nn}.speechDist = 'none';
    
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
masterReady = 0;

while ~masterReady %this bit of code is added incase a slave begins before its master
    try
        worker(xL.opFolder);
        for nn = 1:recConditions
            worker(xR{nn}.opFolder);
        end
        masterReady = 1;
    catch %#ok<CTCH>
        disp('Master node not ready yet: retrying again in 5 seconds');
        pause(5);
    end
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
    while(~all(xR{end}.todoStatus==2)) %Where 2=done (oh how I yearn for enum)
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
