close all; clear all; clc; clear classes;

x = EssexAid_WrapClass;

% Disable feedback loops
x.MOCfactor = 0;
x.ARthreshold_dB = 200;


dt = 1/x.sr;
tAxis = dt:dt:0.1;
f = [x.channelBFs(3)*2^0.2 x.channelBFs(4)*2^-0.2];
dBlev = 60;
s1 = sin(2*pi*f(1)*tAxis);
s1 = s1/max(s1);%/sqrt(mean(s1.^2));

s2 = sin(2*pi*f(2)*tAxis);
s2 = s2/max(s2);%/sqrt(mean(s2.^2));

x = x.processStim