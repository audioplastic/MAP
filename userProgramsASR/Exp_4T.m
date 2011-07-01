function Exp_4T(isMasterNode, num4train)
% This experiment tests a recogniser on 3 different training sets with
% different efferent conditions

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up the basic experiment parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% num4train = 125;

expName = ['4T' num2str(num4train)];
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
xL.numWavs = num4train; %MAx=8440

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
% ** Generate features **
% This is the time consuming, processing intensive portion of the program.
% Nodes that are not the master node are only interested in the opFolder
% member of the jobjects.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
worker(xL.opFolder);



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Train and test the recogniser - a job for the master node only
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isMasterNode    
    y = HMMclass(hmmFolder);    
    y.createSCP(xL.opFolder)
    y.createMLF(xL.opFolder)
    y.train(xL.opFolder) %This node can be busy training, even if other jobs are being processed
    
%     xR{end}.lockJobList;
%     xR{end} = xR{end}.loadSelf; %Reload changes
%     xR{end}.unlockJobList;
%     while(~all(xR{end}.todoStatus==2))        
%         disp('Waiting on straggler nodes to complete their jobs before HMM is tested . . .')
%         pause(20);
%         xR{end}.lockJobList;
%         xR{end} = xR{end}.loadSelf; %Reload incase changed
%         xR{end}.unlockJobList;
%     end
%       
%     for nn = 1:3*recConditions
%         y.createSCP(xR{nn}.opFolder);
%         y.test(xR{nn}.opFolder);
%     end
%     
%     %Show all of the scores in the command window at the end
%     for nn = 1:3*recConditions
%         y.score(xR{nn}.opFolder);
%     end
end
