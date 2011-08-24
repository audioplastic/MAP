close all; clear all; clc

a = 10000;
b = 8e-6;
c = 0.2;

cmpA = 10^((1/(1-c))* log10(b/a)) %THIS IS CORRECT ON THE INPUT AXIS (Stapes vel)
cmpC = a * 10^((1/(1-c))* log10(b/a)) %THIS IS CORRECT on the output axis (BM displacement)


% OK, so the previous version of MAP was missing the factor a in the
% calculation of the compression threshold.

att = 1/10;
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




CtBM = 4.2546e-008; %Compression threshold in units of basilar membrane disp
CtS  = CtBM/a;      %Compression threshold in units of stapes disp

y = zeros(size(x));
abs_x = abs(x);
y(abs_x<CtS)   = a * x(abs_x<CtS);
y(abs_x>=CtS)  = sign(x(abs_x>=CtS)) * a * CtS .* exp(   c * log(  abs_x(abs_x>=CtS)/CtS  )   );



figure;  plot(xOrig*1e12, y*1e12); ylim([0 10e4])
hold on; plot(xOrig([1 end])*1e12, [CtBM CtBM]*1e12, ':r')
hold on; plot(xOrig([1 end])*1e12, [CtBM CtBM]*1e12/4, ':k')
xlabel('Stapes Displacement (nm)'); ylabel('BM Displacement (nm)')


