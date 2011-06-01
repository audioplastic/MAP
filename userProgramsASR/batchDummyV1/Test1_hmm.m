close all; clear all; clc; clear classes

%%
featL     = fullfile('d:', 'exps', '_testExp','featL');
featR{1}  = fullfile('d:', 'exps', '_testExp','featR');
hmmFolder = fullfile('d:', 'exps', '_testExp','hmm');

y = HMMclass(hmmFolder);

y.createSCP(featL)
y.createMLF(featL)
y.train(featL)

%%
nn=1;
y.createSCP(featR{nn})
y.test(featR{nn})
%%
y.score(featR{nn});