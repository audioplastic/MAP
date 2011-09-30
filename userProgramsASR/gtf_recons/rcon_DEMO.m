function rcon_DEMO


close all; clear all; clc;

bwOct = 1/6;

fNmax = 5/bwOct;
cf = 250 * 2.^((0:fNmax)'*bwOct);
lowerCutOff=cf*2^(-bwOct/2);
upperCutOff=cf*2^( bwOct/2);
bwHz = upperCutOff - lowerCutOff;

[x, sr] = wavread('MHS_2841A');

% soundsc(x,sr)

%% Do the old method
fOrd = 6;
nChans = numel(cf)
y = zeros(size(x));
for nn = 1:nChans
    [b,a] = gammatone(bwHz(nn), cf(nn), 1/sr);
    yTmp = x;
    for kk = 1:fOrd
        yTmp = filter(b,a,yTmp);
    end
    
    %Do compression
    
    for kk = 1:fOrd
        yTmp = filter(b,a,yTmp);
    end
    y=y+yTmp;
end
% pause   

%% Do the new zero-phase method
iLen = 1500;
% figure
z = zeros(size(x));
for nn = 1:nChans
    imp = [1 zeros(1,iLen-1)];
    impTmp = imp;
    [b,a] = gammatone(bwHz(nn), cf(nn), 1/sr);
    for kk = 1:fOrd
        impTmp = filter(b,a,impTmp);
    end

%     [b,a] = butter(2, [lowerCutOff(nn) upperCutOff(nn)] ./ (sr/2));
%     impTmp = filter(b,a,impTmp);
    
if nn ==1
     figure; plot(impTmp); drawnow
end

    impTmp = conv(impTmp, fliplr(impTmp));
%     impTmp = impTmp.*hanning(numel(impTmp))';

if nn ==1
     figure; plot(impTmp); drawnow
end
    
%     pause
% %     impTmp = impTmp./sum(impTmp.^2);
%     freqz(impTmp, 1) 
    disp(sum(abs(impTmp)))
%     pause  
    
    
    
    zTmp = x;
    zTmp = filter(impTmp,1,zTmp);
    zTmp = filter(impTmp,1,zTmp);
    z=z+zTmp;
end
figure
[b,a] = butter(4, [lowerCutOff(1) upperCutOff(end)] ./ (sr/2));
x = filter(b,a,x);
plot(z, 'k')
hold on; plot(x, 'r')
% hold on; plot(y,'g')



soundsc([x; y; z],sr)
end

function [b,a] = gammatone(bw, cf, dt)
phi = 2 * pi * bw * dt;
theta = 2 * pi * cf * dt;
cos_theta = cos(theta);
sin_theta = sin(theta);
alpha = -exp(-phi) * cos_theta;
b0 = 1.0;
b1 = 2 * alpha;
b2 = exp(-2 * phi);
z1 = (1 + alpha * cos_theta) - (alpha * sin_theta) * 1i;
z2 = (1 + b1 * cos_theta) - (b1 * sin_theta) * 1i;
z3 = (b2 * cos(2 * theta)) - (b2 * sin(2 * theta)) * 1i;
tf = (z2 + z3) / z1;
a0 = abs(tf);
a1 = alpha * a0;

a = [b0, b1, b2];
b = [a0, a1];
end% ------ OF gammatone