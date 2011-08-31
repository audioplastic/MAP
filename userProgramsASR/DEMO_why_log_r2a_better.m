close all; clear all; clc

rate = 50:250;
MOCrateThresholdProb = 85;
rateToAttenuationFactorProb = 6;


x = -20*log10(  max(rate/MOCrateThresholdProb,1)  )*rateToAttenuationFactorProb; %dB attenuation
x = max(x,-35);


rateToAttenuationFactorProb = 0.010;
rateOverThreshold = rate - MOCrateThresholdProb;
rateOverThreshold(rateOverThreshold<0)=0;
y =  (1- rateOverThreshold* rateToAttenuationFactorProb);
y(y<0)=0.001;
y = 20*log10(y);
y = max(y,-35);

plot(rate, -x, rate, -y)
xlabel('Firing rate')
ylabel('Attenuation (dB)')
ylim([-5 40])

legend('New', 'Old',2)