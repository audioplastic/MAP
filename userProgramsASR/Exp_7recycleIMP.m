function Exp_7recycleIMP(isMasterNode)
% This experiment tests a recogniser on 3 different training sets with
% different efferent conditions.
% This is now using paramChanges and the hearing aid to correct an OHC
% dysfunction.

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up the basic experiment parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
expName = '7hiC_0';
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

xL.participant = 'NormalDIFF';%'NormalDIFF';
xL.MAPparamChanges= {' DRNLParams.rateToAttenuationFactorProb = 0;', 'OMEParams.rateToAttenuationFactorProb=0;'};

xL.noiseLevToUse   =  -200;
xL.speechLevToUse  =  60;

xL.MAPopHSR = 1;
xL.MAPopMSR = 0;
xL.MAPopLSR = 0;


xL.numCoeff = 14;
xL.removeEnergyStatic = 0;

%%%%% Group of params that will influence simulation run time %%%%%%%
xL.numWavs = 8440; %MAx=8440
testWavs = 358; %MAX = 358
nzLevel = 45:5:70;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xL.noisePreDur = 1;
xL.noisePostDur = 0.1;
xL.truncateDur  = 0.9; %Dr. RF used 0.550

xL.noiseName = 'pink';

% if isMasterNode
%     mkdir(xL.opFolder);
%     xL = xL.assignFiles;
%     xL.storeSelf;
% end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sort out the testing (RECOGNITION) conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% xR=cell(size(nzLevel));
recConditions = numel(nzLevel);
% for nn = 1:recConditions    
%     xR{nn} = xL; %simply copy the "Learn" object and change it a bit below
%     recFolder = fullfile(expFolder,['NONEfeatR_nz' num2str(nzLevel(nn))]);
%     xR{nn}.opFolder = recFolder;    
%     
%     %These are the interesting differences between training and testing
%     xR{nn}.numWavs = testWavs; %MAX = 358
%     xR{nn}.noiseLevToUse = nzLevel(nn);    
%     
%             
%     %Now just to wrap it up ready for processing
%     if isMasterNode && ~isdir(xR{nn}.opFolder)
%         mkdir(xR{nn}.opFolder);
%         xR{nn} = xR{nn}.assignWavPaths('R');
%         xR{nn} = xR{nn}.assignFiles;
%         xR{nn}.storeSelf;
%     end
% end



tmpIdx=0;
for nn = 12*recConditions+1:13*recConditions   
    tmpIdx=tmpIdx+1;
    xR{nn} = xL; %simply copy the "Learn" object and change it a bit below
    recFolder = fullfile(expFolder,['IMPfeatR_cond' num2str(nn)]);
    xR{nn}.opFolder = recFolder;    
    
    %These are the interesting differences between training and testing
    xR{nn}.numWavs = testWavs; %MAX = 358
    xR{nn}.noiseLevToUse = nzLevel(tmpIdx);
    xR{nn}.MAPparamChanges= {'DRNLParams.a=400;','DRNLParams.rateToAttenuationFactorProb = 0.0;', 'DRNLParams.MOCrateThresholdProb =57;', 'DRNLParams.MOCtau =0.35;'};
    
    
    %Now just to wrap it up ready for processing
    if isMasterNode && ~isdir(xR{nn}.opFolder)
        mkdir(xR{nn}.opFolder);
        xR{nn} = xR{nn}.assignWavPaths('R');
        xR{nn} = xR{nn}.assignFiles;
        xR{nn}.storeSelf;
    end
end

tmpIdx=0;
for nn = 13*recConditions+1:14*recConditions   
    tmpIdx=tmpIdx+1;
    xR{nn} = xL; %simply copy the "Learn" object and change it a bit below
    recFolder = fullfile(expFolder,['IMPfeatR_cond' num2str(nn)]);
    xR{nn}.opFolder = recFolder;    
    
    %These are the interesting differences between training and testing
    xR{nn}.numWavs = testWavs; %MAX = 358
    xR{nn}.noiseLevToUse = nzLevel(tmpIdx);
    xR{nn}.MAPparamChanges= {'DRNLParams.a=600;','DRNLParams.rateToAttenuationFactorProb = 0.0;', 'DRNLParams.MOCrateThresholdProb =57;', 'DRNLParams.MOCtau =0.35;'};
    
    
    %Now just to wrap it up ready for processing
    if isMasterNode && ~isdir(xR{nn}.opFolder)
        if ~isdir(xR{nn}.opFolder);
            mkdir(xR{nn}.opFolder);
            xR{nn} = xR{nn}.assignWavPaths('R');
            xR{nn} = xR{nn}.assignFiles;
            xR{nn}.storeSelf;
        end
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ** Generate features **
% This is the time consuming, processing intensive portion of the program.
% Nodes that are not the master node are only interested in the opFolder
% member of the jobjects.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% worker(xL.opFolder);

if ~isMasterNode %dont bother wasting master node effort on generating testing features (for now)
    for nn = 12*recConditions+1:14*recConditions
        worker(xR{nn}.opFolder);
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Train and test the recogniser - a job for the master node only
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isMasterNode    
%     while(~all(xL.todoStatus==2))        
%         disp('Waiting on straggler nodes to complete their jobs before HMM is trained . . .')
%         pause(30); %Wait for 30 seconds before looking again
%         xL.lockJobList;
%         xL = xL.loadSelf; %Reload incase changed
%         xL.unlockJobList;
%     end
    y = HMMclass(hmmFolder);    
    y.numCoeff = 14*3;
%     y.createSCP(xL.opFolder)
%     y.createMLF(xL.opFolder)
%     y.train(xL.opFolder) %This node can be busy training, even if other jobs are being processed for testing
    
    % ALLOW MASTER NODE TO MUCK IN WITH GENERATING TESTING FEATURES ONCE
    % HMM HAS BEEN TRAINED
    for nn = 12*recConditions+1:14*recConditions
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
      
    for nn = 12*recConditions+1:14*recConditions
        y.createSCP(xR{nn}.opFolder);
        y.test(xR{nn}.opFolder);
    end
    
    %Show all of the scores in the command window at the end
    for nn = 12*recConditions+1:14*recConditions
        y.score(xR{nn}.opFolder);
    end
end
