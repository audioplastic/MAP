function testPhysiologyProb(BF,paramsName, paramChanges)
% testPhysiologyProb is a portmanteau programs that
%   invokes a number of physiological assessments to test the
%   'probability' model. Tests invoked are 
%   testOME, testBM testRP2 testSynapse testFM testANprob
%Input arguments:
%  BF: single channel model best frequency
%  paramsName: parameter file name containing model parameters.
%   (default='Normal')
%    NB the program assumes that two fiber types are nominated, i.e. two
%    values of ANtauCas are specified.
%  paramChanges: cell array contining list of changes to parameters. These
%   are implemented after reading the parameter file (default='')
%
% e.g.
%  testPhysiologyProb(1000,'Normal', '');

restorePath=path;
addpath (['..' filesep 'MAP'])

if nargin<3
    paramChanges=[];
end

if nargin==0
    error('testPhysiologyProb must be called from the command line')
end


disp('testPhysiologyProb...........computing')

disp('testOME...........computing')
testOME(paramsName, paramChanges)

disp('testBM...........computing')
relativeFrequencies=[0.25    .5   .75  1  1.25 1.5    2];
testBM (BF, paramsName,relativeFrequencies,'probability', paramChanges)

disp('testRP...........computing')
testRP2(paramsName,paramChanges)

disp('testSynapse...........computing')
testSynapse(BF,paramsName, 'probability', paramChanges)

disp('testFM...........computing')
testFM(BF,paramsName,'probability', paramChanges)

disp('testANprob...........computing')
testANprob(BF,BF, -10:10:80,paramsName, paramChanges);

figure(4)

path(restorePath)
