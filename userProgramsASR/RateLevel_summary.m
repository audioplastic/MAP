close all

% paramChanges{numel(paramChanges)+1} = 'DRNLParams.rateToAttenuationFactorProb = 0.00;';%GOOD = 0.012  %DEFAULT = 0.005;  % strength of MOC
% paramChanges{numel(paramChanges)+1} = 'DRNLParams.MOCrateThresholdProb = 140;';%GOOD=140 %DEFAULT = 70;
% paramChanges{numel(paramChanges)+1} = 'OMEParams.rateToAttenuationFactorProb = 0;';%DEFAULT = 0.01;
% paramChanges{numel(paramChanges)+1} = 'DRNLParams.a=1e4;'; %DEFAULT = 5e4;

data(1,:) = [57.7536   63.4700  107.1336  179.2254  187.4560  191.0853  194.3707  197.7813  201.2181  204.3290  206.2312];
data(2,:) = [57.1659   57.7662   63.5498  107.1202  178.7683  187.5097  191.7263  195.8818  200.1848  203.9475  206.1200];
data(3,:) = [57.7536   63.4700   94.9067  128.0028  157.8979  174.8033  179.8676  181.6341  197.4968  204.3665  205.7559];

levAx = 0:10:100;

figure(50)
plot(levAx, data')
legend(...
    'Normal',...
    'Fixed 10 dB att''n',...
    'Auto (T=) (F=)',...
    4);
title('1 kHz')

