close all; clear all; clc

sr = 48e3;

onTC = 500e-3;
offTC = 50e-3;

durTA = 5;
durTP = 3;

dt = 1/sr;
ponT = (durTA-durTP)/2;
poffT = ponT+durTP;
durS = ceil(sr*durTA);
ponSidx = round(sr*ponT):round(sr*poffT)-1;

sig = zeros(durS,1);
sig(ponSidx) = 1;

%Rise coefficients
aR = dt/onTC - 1;
bR = 1.0 + aR;

%Fall coefficients
aF = dt/offTC - 1;
bF = 1.0 + aR;

tic
sigOUT = zeros(size(sig));
sigOLD=0;
for nn = 1:durS
    if sig(nn) < sigOLD %// - This is line to make smoothing only apply to release
        sigOUT(nn) = bF*sig(nn) - aF*sigOLD;% // difference eqn for one-pole lpf
    else
        sigOUT(nn) = bR*sig(nn) - aR*sigOLD;% // difference eqn for one-pole lpf
    end
    sigOLD = sigOUT(nn);
end
toc                

plot(sigOUT)                
                
                
