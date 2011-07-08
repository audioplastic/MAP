function testPhysiologyProb(BF,paramsName, paramChanges)

restorePath=path;
addpath (['..' filesep 'MAP'])

if nargin<3
    paramChanges=[];
end

testOME(paramsName, paramChanges)
relativeFrequencies=[0.25    .5   .75  1  1.25 1.5    2];
testBM (BF, paramsName,relativeFrequencies,'probability', paramChanges)
testRP(BF,paramsName, paramChanges)
testSynapse(BF,paramsName, paramChanges)
testANprob(BF,BF, -10:10:80,paramsName, paramChanges);

figure(4)

path(restorePath)
