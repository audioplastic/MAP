function Exp_4_imp(isMasterNode)
% This experiment tests a recogniser on 3 different training sets with
% different efferent conditions

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up the basic experiment parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
expName = '4_imp';
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

xL.participant = 'NormalNONE';%'NormalDIFF';

xL.noiseLevToUse   =  -200;
xL.speechLevToUse  =  60;

xL.MAPopHSR = 1;
xL.MAPopMSR = 0;
xL.MAPopLSR = 0;
xL.MAPuseEfferent = 1;
xL.numWavs = 1000; %MAx=8440

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
nzLevel = [-200 40:5:60];
xR=cell(size(nzLevel));
recConditions = numel(nzLevel);
for nn = 1:recConditions    
    xR{nn} = xL; %simply copy the "Learn" object and change it a bit below
    recFolder = fullfile(expFolder,['NONEfeatR_nz' num2str(nzLevel(nn))]);
    xR{nn}.opFolder = recFolder;    
    
    %These are the interesting differences between training and testing
    xR{nn}.numWavs = 358; %MAX = 358
    xR{nn}.noiseLevToUse = nzLevel(nn);
    xR{nn}.participant =  'NormalNONE'; % This line is redundent but included for documentation purpose
            
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
    recFolder = fullfile(expFolder,['AIDfeatR_nz' num2str(nzLevel(tmpIdx))]);
    xR{nn}.opFolder = recFolder;    
    
    %These are the interesting differences between training and testing
    xR{nn}.numWavs = 358; %MAX = 358
    xR{nn}.noiseLevToUse = nzLevel(tmpIdx);
    xR{nn}.participant = 'OHC';
    
    xR{nn}.useAid = 0;
    
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
    recFolder = fullfile(expFolder,['AUTOfeatR_nz' num2str(nzLevel(tmpIdx))]);
    xR{nn}.opFolder = recFolder;    
    
    %These are the interesting differences between training and testing
    xR{nn}.numWavs = 358; %MAX = 358
    xR{nn}.noiseLevToUse = nzLevel(tmpIdx);
    xR{nn}.participant = 'OHC';
    
    xR{nn}.useAid = 1;
    xR{nn}.mainGain = [10; 10; 10; 10; 10];
    
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

if ~isMasterNode %dont bother wasting master node effort on testing features
    for nn = 1:3*recConditions
        worker(xR{nn}.opFolder);
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
    while(~all(xR{end}.todoStatus==2))        
        disp('Waiting on straggler nodes to complete their jobs before HMM is tested . . .')
        pause(20);
        xR{end}.lockJobList;
        xR{end} = xR{end}.loadSelf; %Reload incase changed
        xR{end}.unlockJobList;
    end
      
    for nn = 1:3*recConditions
        y.createSCP(xR{nn}.opFolder);
        y.test(xR{nn}.opFolder);
    end
    
    %Show all of the scores in the command window at the end
    for nn = 1:3*recConditions
        y.score(xR{nn}.opFolder);
    end
end
