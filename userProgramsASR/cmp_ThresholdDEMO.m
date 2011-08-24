function cmp_ThresholdDEMO

clc; close all

a = 10000;
b = 8e-6; %here for legacy reasons to calculate a reasonable BM cmp threshold
c = 0.2;
MOC = 1/1; %Fraction between 0 and 1


% These are not used - only here to display thresholds based on old
% parameters
cmpS  =     10^((1/(1-c))* log10(b/a)) % (Stapes vel)
cmpBM = a * 10^((1/(1-c))* log10(b/a)) % (BM displacement)

% First test is to make an I/O function
xOrig = 0:1e-14:100e-12;
x = xOrig*MOC;
CtBM = 4.2546e-008; %Compression threshold in units of basilar membrane disp
y = NewNonLinFunc(x,a,c,CtBM);

figure;  plot(xOrig*1e12, y*1e12); ylim([0 10e4])
hold on; plot(xOrig([1 end])*1e12, [CtBM CtBM]*1e12, ':r')
% hold on; plot(xOrig([1 end])*1e12, [CtBM CtBM]*1e12/4, ':k')
xlabel('Stapes Displacement (nm)'); ylabel('BM Displacement (nm)')


%NOw quickly test with a sine wave that peaks above BM compression thresh
dt = 1/25e3;
tAxis = dt:dt:0.005;
x = 3 * CtBM/a * sin(2*pi*1000*tAxis);
y = NewNonLinFunc(x,a,c,CtBM);
figure;  plot(1000*tAxis, y*1e12);
hold on; plot(1000*tAxis([1 end]), [CtBM CtBM]*1e12, ':r')
hold on; plot(1000*tAxis([1 end]), -[CtBM CtBM]*1e12, ':r')
xlabel('Time (ms)'); ylabel('BM Displacement (nm)')


function y = NewNonLinFunc(x,a,c,CtBM)
CtS  = CtBM/a;      %Compression threshold in units of stapes disp
y = zeros(size(x));
abs_x = abs(x);
y(abs_x<CtS)   = a * x(abs_x<CtS);
y(abs_x>=CtS)  = sign(x(abs_x>=CtS)) * a * CtS .* exp(   c * log(  abs_x(abs_x>=CtS)/CtS  )   );
