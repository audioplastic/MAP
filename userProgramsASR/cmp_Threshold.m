close all; clear all; clc

a = 1000;
b = 8e-6;
c = 0.2;

cmpA = 10^((1/(1-c))* log10(b/a)) %THIS IS WRONG!!!!!
cmpB = ((b^5)/a)^0.25 %This is valid only for DRNLc = 0.2
cmpC = a * 10^((1/(1-c))* log10(b/a)) %THIS IS CORRECT!!!!!


% OK, so the previous version of MAP was missing the factor a in the
% calculation of the compression threshold.

att = 1;
xOrig = 0:1e-14:100e-12;

x = xOrig*att;

aPath = a*x;
bPath = b*(x).^c;

plot(xOrig*1e12,aPath*1e12,xOrig*1e12,bPath*1e12); ylim([0 10e4])
% set(gca,'XScale', 'log', 'YScale', 'log')

y= x.* a;  % linear section.
% compress parts of the signal above the compression threshold
abs_x = abs(x);
idx=find(abs_x>cmpA);
if ~isempty(idx)>0
    y(idx)=sign(y(idx)).* (b*abs_x(idx).^c);
end
nonlinOutput=y;

figure; plot(xOrig*1e12,nonlinOutput*1e12); ylim([0 10e4])
set(gca,'XScale', 'log', 'YScale', 'log')