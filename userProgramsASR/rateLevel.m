close all; clear all; clc; clear classes; clear path

addpath ('..\modules', '..\utilities',  '..\parameterStore',  '..\wavFileStore' , '..\testPrograms')


durT = 5; %seconds
f = 1000; % Hz
% cLev = 110; % dB SPL
participant        = 'DEMO2_multiSpont';
MAPuseEfferent = 1;

MAPopLSR = 0;
MAPopMSR = 0;
MAPopHSR = 1;

sr = 25e3;
dt = 1/sr;
tAxis = dt:dt:durT;

% wavFile = fullfile('batchDummyV1','demo_wavs','microAURORA','MFF_3B.wav');
wavFile = fullfile('batchDummyV1','demo_wavs','microAURORA','MHS_2841A.wav');

stream = RandStream('mt19937ar','Seed',0);
% x = randn(stream,size(tAxis));
x = sin(2*pi*f*tAxis); %Make sine wave
% [x,sr] = wavread(wavFile);
x=x';

dt = 1/sr;
tAxis = dt:dt:durT;

%**********************************************************
% Map related bits
%**********************************************************
global AN_IHCsynapseParams inputStimulusParams IHCpreSynapseParams DRNLParams; %#ok<NUSED>

lowestBF=250; 	highestBF= 8000; 	numChannels=21;
    % 21 chs (250-8k)includes BFs at 250 500 1000 2000 4000 8000
    BFlist=round(logspace(log10(lowestBF),log10(highestBF),numChannels));

paramsToLoad = ['method=MAPparams', participant, '([' num2str(BFlist) '],' num2str(sr) ');'];
eval(paramsToLoad);
%             DRNLParams.fixedMOCAttenuation= -20;

% now that the correct parameters have been loaded, change them!

method.plotGraphs =	0;

%**********************************************************
% If using efferent change some params pre processing
%**********************************************************
if MAPuseEfferent
    method.saveANprobability = 1; %new
    method.useEfferent = 1;
    MAPmoduleSequence       = [1 2 3 4 5 6 7 8 9];
    AN_IHCsynapseParams.mode= 'spikes'; %now define here and hide from class interface (the global)
else
    method.saveANprobability = 0; %new
    method.useEfferent = 0;
    MAPmoduleSequence       = [1 2 3 4 5 6 7];
    AN_IHCsynapseParams.mode= 'probability';
    
    %Below is some stuff related to HSR / LSR
    %This could be done lower down but it makes little
    %difference. Probably more efficient to just not calculate
    %it for probability mode
    
    if numel(IHCpreSynapseParams.tauCa)==2
        %                     LSRtauCa=IHCpreSynapseParams.tauCa(1);
        HSRtauCa=IHCpreSynapseParams.tauCa(end);
    elseif numel(IHCpreSynapseParams.tauCa)==3
        LSRtauCa=IHCpreSynapseParams.tauCa(1);
        MSRtauCa=IHCpreSynapseParams.tauCa(2);
        HSRtauCa=IHCpreSynapseParams.tauCa(end);
    else
        LSRtauCa=[];
        HSRtauCa=IHCpreSynapseParams.tauCa; %incase the model used only has one fibre type
    end
    
    IHCpreSynapseParams.tauCa= []; %wipe it
    if MAPopLSR
        assert(~isempty(LSRtauCa), 'No LSR data available for this parameter file')
        IHCpreSynapseParams.tauCa = LSRtauCa;
    end
    if MAPopMSR
        assert(~isempty(MSRtauCa), 'No MSR data available for this parameter file')
        IHCpreSynapseParams.tauCa = [IHCpreSynapseParams.tauCa MSRtauCa];
    end
    if MAPopHSR
        IHCpreSynapseParams.tauCa = [IHCpreSynapseParams.tauCa HSRtauCa];
    end
end

% dBaxis = (0:8:120)';
dBaxis = 80;
rateLev = zeros(size(dBaxis));
for nn = 1:numel(dBaxis);
    
    x = x./sqrt(mean(x.^2)); %Normalize RMS to 1
    x = x * 20e-6*10^(dBaxis(nn)/20); %Convert RMS to pascals at desired level
    
    tic
    [ANprobabilityResponse, method, ANrate] = MAPsequenceSeg(x, method, MAPmoduleSequence);
    toc
    
    %**********************************************************
    % If using efferent do some different post processing
    %**********************************************************
    if MAPuseEfferent
        ANprobabilityResponse = method.ANprobabilities;%(logical([MAPopLSR MAPopMSR MAPopHSR]), :);
        fRate = ANprobabilityResponse/method.AN_IHCsynapsedt;
        if MAPopLSR
            fRate = fRate(1:numChannels,:);
        elseif MAPopMSR
            fRate = fRate(numChannels+1:2*numChannels,:);
        elseif MAPopHSR
            fRate = fRate(2*numChannels+1:end,:);
        end
    else
        fRate = ANrate;
    end
    ANtAxis = method.AN_IHCsynapsedt:method.AN_IHCsynapsedt:method.AN_IHCsynapsedt*length(ANprobabilityResponse);
    plot(ANtAxis,fRate(9,:));  drawnow;
    
    %
    % fRate = ANprobabilityResponse*1000;
    % fRate = ANprobabilityResponse;
    
    rateLev(nn) = mean(fRate(9,100:end));
end

figure; plot(dBaxis, rateLev); xlim([0 max(dBaxis)]); ylim([0 400])
xlabel('Input level in dB SPL')
ylabel('Spikes per second')

op4excel = [dBaxis rateLev];


