close all; clear all; clc; clear classes;

x = EssexAid_WrapClass;

% Disable feedback loops
x.MOCfactor = 0;
x.ARthreshold_dB = 200;


dt = 1/x.sr;
tAxis = dt:dt:0.2;
f = [x.channelBFs(3)*2^0.2 x.channelBFs(4)*2^-0.2];
dBlev = 60;
s1 = sin(2*pi*f(1)*tAxis);


s2 = sin(2*pi*f(2)*tAxis);


figure; plot(s1+s2)


s1 = s1/sqrt(mean(s1.^2));
s2 = s2/sqrt(mean(s2.^2));

s1 = s1*20e-6*10^(dBlev/20);
s2 = s2*20e-6*10^(dBlev/20);

x.stimulusUSER = s1+s2;
x = x.processStim;
figure; plot(x)

% Beating peaks @50 dB SPL