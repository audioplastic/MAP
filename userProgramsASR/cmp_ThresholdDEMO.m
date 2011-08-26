function cmp_ThresholdDEMO

clc; close all

a = 10000;
b = 8e-6; %here for legacy reasons to calculate a reasonable BM cmp threshold
c = 0.2;
MOC = 0.08; %Fraction between 0 and 1


% These are not used - only here to display thresholds based on old
% parameters
cmpS  =     10^((1/(1-c))* log10(b/a)) % (Stapes vel)
cmpBM = a * 10^((1/(1-c))* log10(b/a)) % (BM displacement)

% First test is to make an I/O function
xOrig = 0:1e-14:400e-12;
x = xOrig*MOC;

%using old method
y= x.* a;  % linear section.
% compress parts of the signal above the compression threshold
abs_x = abs(x);
idx=find(abs_x>cmpS);
if ~isempty(idx)>0
    y(idx)=sign(y(idx)).* (b*abs_x(idx).^c);
end

figure;  plot(xOrig*1e12, y*1e12); ylim([1000 10e4]); xlim([1e-12 max(xOrig)]*1e12)
hold on; plot(xOrig([1 end])*1e12, [cmpBM cmpBM]*1e12, ':r')
xlabel('Stapes Displacement (nm)'); ylabel('BM Displacement (nm)')

% Using new method
CtBM = 4.2546e-008; %Compression threshold in units of basilar membrane disp
% CtBM = 8.0000e-008;
% CtBM = cmpBM;
y = NewNonLinFunc(x,a,c,CtBM);

figure;  plot(xOrig*1e12, y*1e12); ylim([1000 10e4]); xlim([1e-12 max(xOrig)]*1e12)
hold on; plot(xOrig([1 end])*1e12, [CtBM CtBM]*1e12, ':r')
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
