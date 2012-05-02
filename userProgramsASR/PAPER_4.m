function PAPER_4(isMasterNode)
% This experiment tests a recogniser on 3 different training sets with
% different efferent conditions.
% This is now using paramChanges and the hearing aid to correct an OHC
% dysfunction.

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up the basic experiment parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
expName = '18ax';
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

xL.participant = 'NormalDIFFstd';%'NormalDIFF';
xL.MAPparamChanges= { 
                'OMEParams.rateToAttenuationFactorProb=0;',...
                'OMEParams.ARrateThreshold=30;',... %Threshold of 40 makes AR kick off around 65 dB for bb noise
                'OMEParams.ARtau=0.1;',...
                'DRNLParams.MOCtauR=2;',...
                'DRNLParams.MOCtauF=DRNLParams.MOCtauR;',...                
                'DRNLParams.rateToAttenuationFactorProb=0;',...
                'DRNLParams.MOCrateThresholdProb=85;'};
            
xL.noiseLevToUse   =  -200;
xL.speechLevToUse  =  60;
xL.speechDist = 'none';
xL.noiseDist  = 'none';
xL.speechLevStd    = 0/sqrt(12);
xL.noiseLevStd    = 0;
%xL.meanSNR = 200;

xL.MAPopHSR = 1;
xL.MAPopMSR = 0;
xL.MAPopLSR = 0;


xL.numCoeff = 14;
xL.removeEnergyStatic = 0;

%%%%% Group of params that will influence simulation run time %%%%%%%
xL.numWavs = 8440; %MAx=8440
testWavs = 358; %MAX = 358
nzLevel = [-200  40:5:70 ];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xL.noisePreDur = 6; %only short lead in needed for SRT type test
xL.noisePostDur = 0.1;
xL.truncateDur  = xL.noisePreDur-0.1;

xL.noiseName = 'BabbleShapedNoise_bp';

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


nn=0;

ATT = 2;
while ATT < 40
    
    tmpIdx=0;
    
    for kk = 1:recConditions
        
        tmpIdx=tmpIdx+1;
        nn=nn+1;
        
        xR{nn} = xL; %simply copy the "Learn" object and change it a bit below
        recFolder = fullfile(expFolder,['FIX' num2str(ATT) '_babn_'  num2str(nn)]);
        xR{nn}.opFolder = recFolder;
        
        %These are the interesting differences between training and testing
        xR{nn}.numWavs = testWavs; %MAX = 358
        xR{nn}.speechLevToUse = 60;%spLevel(tmpIdx);
        xR{nn}.noiseLevToUse = nzLevel(tmpIdx);
        xR{nn}.speechDist = 'None';
        xR{nn}.noiseDist = 'None';
        xR{nn}.MAPparamChanges= {
            'OMEParams.rateToAttenuationFactorProb=0;',...
            'OMEParams.ARrateThreshold=30;',... %Threshold of 40 makes AR kick off around 65 dB for bb noise
            'OMEParams.ARtau=0.1;',...
            'DRNLParams.MOCtauR=2;',...
            'DRNLParams.MOCtauF=DRNLParams.MOCtauR;',...
            ['DRNLParams.rateToAttenuationFactorProb=-10^(-' num2str(ATT) '/20);'],...
            'DRNLParams.MOCrateThresholdProb=85;'};
        
        %Now just to wrap it up ready for processing
        if isMasterNode && ~isdir(xR{nn}.opFolder)
            mkdir(xR{nn}.opFolder);
            xR{nn} = xR{nn}.assignWavPaths('R');
            xR{nn} = xR{nn}.assignFiles;
            xR{nn}.storeSelf;
        end                
        
    end
    
    ATT = ATT + 2;
    if ATT == 10 % Don't bother as we already done that one!
        ATT = ATT + 2;
    end
    
end





%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ** Generate features **
% This is the time consuming, processing intensive portion of the program.
% Nodes that are not the master node are only interested in the opFolder
% member of the jobjects.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% worker(xL.opFolder);
maxConds = nn;
if ~isMasterNode %dont bother wasting master node effort on generating testing features (for now)
    for nn = 1:maxConds
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
    for nn = 1:maxConds
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
      
    for nn = 1:maxConds
        y.createSCP(xR{nn}.opFolder);
        y.test(xR{nn}.opFolder);
    end
    
    %Show all of the scores in the command window at the end
    for nn = 1:maxConds
        y.score(xR{nn}.opFolder);
    end
end