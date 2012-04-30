function audio = callEssexAid(audio,MAPparamsName,sampleRate,paramChanges)
addpath(['..' filesep 'userProgramsASRforDummies']); %path to the Aid
addpath(['..' filesep 'Aidparameters']); %path to the Aidparameters
%look, if there is a paramChange enabling the hearing aid
nChanges=length(paramChanges);
for idx=1:nChanges
    tmp = paramChanges{idx};
    if strfind(tmp,'inputStimulusParams.useAid = 1;')
        useAid = 1;
    end
end
%only call the aid if it is set in the paramChanges
if exist('useAid','var')
    if useAid
%set the aid's parameters
        aidparameterfile = ['Aidparameters' MAPparamsName];
        %xT = cJob;
        aidInstance = cEssexAid;
        eval(['aidInstance = ' aidparameterfile ';']);
        aidInstance.sr = sampleRate;
        
        %call the aid
        aidInstance.stimulusUSER = audio;
        tic
        aidInstance = aidInstance.processStim;
        toc
        audio = aidInstance.aidOPnice;
    end
end