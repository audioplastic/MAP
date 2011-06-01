% doSCPorMLF.m
%
% Rob's sneaky program to fix his latest cock-up
% 13th August 2009
%

clc

experimentName = 'MAP1_10impairedInSilence';
ASRpath = fullfile('C:\Users\rferry\Uni\Post-Doc\Projects\SpeechRecognition\ASRexperiments', filesep);
experimentPath = fullfile(ASRpath, 'MAPv110', filesep, experimentName);
ListsAndLabels = fullfile(experimentPath, filesep, ['ListsAndLabels_', experimentName], filesep);
FeaturePath = fullfile(experimentPath, ['Features_' experimentName]);
testingFeaturePath = fullfile(FeaturePath, filesep, 'TestFeatures', filesep);

if ~isdir(ListsAndLabels)
    mkdir(ListsAndLabels)
end

%currentParticipant = 'ECrARMOC';
currentParticipant = 'ECr';
currentLevel = 65;
identifier = ['Impaired-Test_', currentParticipant];
noiseLevelsToUse = {'noiseLevel-200', 'noiseLevel40', 'noiseLevel45', 'noiseLevel50', 'noiseLevel55', 'noiseLevel60', 'noiseLevel65', 'noiseLevel70', 'noiseLevel75', 'noiseLevel80', 'noiseLevel85', 'noiseLevel90', 'noiseLevel95'};
speechLevelsToUse = {'speechLevel40', 'speechLevel45', 'speechLevel50', 'speechLevel55', 'speechLevel60', 'speechLevel65', 'speechLevel70','speechLevel75', 'speechLevel80', 'speechLevel85', 'speechLevel90',  'speechLevel95'};

for n = 1:length(noiseLevelsToUse)
    noise = noiseLevelsToUse{n};
    
    for s = 1:length(speechLevelsToUse)
        speech= speechLevelsToUse{s};
        
        outputPath = fullfile(testingFeaturePath, identifier);
        targetPath = fullfile(outputPath, [noise, filesep, speech], filesep);
        
        createSCP([ListsAndLabels [identifier '_' noise '_' speech '_testFiles']], targetPath); % filename list (testing)
        createMLFFiles([ListsAndLabels [identifier '_' noise '_' speech '_testLabels']], targetPath, 0); % labels for testing data
        
    end
end

fprintf('\n\n All done!... \n')
