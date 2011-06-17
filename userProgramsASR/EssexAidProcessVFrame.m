%#eml
function [ frameBuffer, filterStates, ARamp, MOCcontrol ] =...
    EssexAidProcessVFrame( ...
    frameBuffer,...
    filterStates,...
    filterCoeffs,...
    numChannels,...
    numSamples,...
    ARamp,...
    MOCcontrol,...
    ARthresholdPa,...
    filterOrder,...
    DRNLaBaseline,...
    DRNLb,...
    DRNLc,...
    MOCthreshold,...
    MOCfactor,...
    stapesScalar)

%ESSEXAIDPROCESSFRAME Essex aid algorithm in frame processing mode
%   This code will look a bit odd to most Matlab programmers. This is
%   because the intended target is a C function that will be called on a
%   sub-millisecond basis. The bizzare enumerations assist the
%   pass-by-reference functionality that allows this function to fly in
%   real-time. This function works on a need to know basis, eliminating any
%   unnecessary data copying or parameter calculation.

% eml.varsize('frameBuffer', 6192);

%Fake enumeration - must be kept up to date with juce enum
enumC_ARb = 0;
enumC_ARa = 2;
enumC_MOCb = 4;
enumC_MOCa = 6;

% enumC_BPb1 = 8;
% enumC_BPa1 = 13;
% enumC_BPb2 = 18;
% enumC_BPa2 = 23;
% enumC_BPb3 = 28;
% enumC_BPa3 = 33;
% enumC_BPb4 = 38;
% enumC_BPa4 = 43;

enumS_AR   = 0;

% enumS_MOC1  = 1;
% enumS_BPin_1_1 = 2;
% enumS_BPin_2_1 = 6;
% enumS_BPout_1_1 = 10;
% enumS_BPout_2_1 = 14;
% 
% enumS_MOC2 = 18;
% enumS_BPin_1_2 = 19;
% enumS_BPin_2_2 = 23;
% enumS_BPout_1_2 = 27;
% enumS_BPout_2_2 = 31;
% ...

% rmsLev[0] = iunput RMS from AR smoothed response

%Initial gain
% frameBuffer(1:numSamples) = frameBuffer(1:numSamples)*ipScalar;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ACOUSTIC REFLEX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find  rms of smoothed input signal
%  this will be used to trigger the AR reflex
[y,  filterStates(enumS_AR+1)] = filter(filterCoeffs(enumC_ARb+1:enumC_ARb+2), filterCoeffs(enumC_ARa+1:enumC_ARa+2) , frameBuffer(1:numSamples).^2, filterStates(enumS_AR+1));
% restore Pa scale
smoothedARrms = sqrt(y);  %confusing name for parameter - it is a short term RMS.

% attenuate input (NB cross product used)
frameBuffer(1:numSamples) = frameBuffer(1:numSamples)./ARamp(1:numSamples);

%CALC ARamp FOR NEXT FRAME
% compare levels in the previous segment with AR threshold
ARamp(1:numSamples) = smoothedARrms/ARthresholdPa;
% all sub-treshold values are set to 1
ARamp(ARamp(1:numSamples)<1)=1;

% tympanic membrane response in meters
stapesAR = frameBuffer(1:numSamples)*stapesScalar;
frameBuffer(1:numSamples) = zeros(size(frameBuffer(1:numSamples)));

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DRNL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for filterCount=1:numChannels
    y=stapesAR;    
    for nn = 1:filterOrder/2;
        [y, filterStates(17*(filterCount-1)+3+(nn-1)*4   :  17*(filterCount-1)+6+(nn-1)*4)] =  filter(filterCoeffs(10*(filterCount-1)+9:10*(filterCount-1)+13), filterCoeffs(10*(filterCount-1)+14:10*(filterCount-1)+18), y, filterStates(17*(filterCount-1)+3+(nn-1)*4   :  17*(filterCount-1)+6+(nn-1)*4));
    end
    
    
    %                            NonLin compression function
    % first identify MOC attenuation to be applied (see below)
    %  this will be different for each channel
    MOC=MOCcontrol(filterCount,1:numSamples);
    DRNLa = DRNLaBaseline(filterCount)./MOC;
    y = DRNL_brokenstick_nl (y, DRNLa, DRNLb(filterCount), DRNLc(filterCount));
    
    assert(filterOrder == 2 || filterOrder == 4, 'filterOrder must be 2 or 4')
    for nn = 1:filterOrder/2;
        [y, filterStates(17*(filterCount-1)+11+(nn-1)*4   :  17*(filterCount-1)+14+(nn-1)*4)] =  filter(filterCoeffs(10*(filterCount-1)+9:10*(filterCount-1)+13), filterCoeffs(10*(filterCount-1)+14:10*(filterCount-1)+18), y, filterStates(17*(filterCount-1)+11+(nn-1)*4   :  17*(filterCount-1)+14+(nn-1)*4));
    end
    
    frameBuffer(1:numSamples) = frameBuffer(1:numSamples) + y;

    % compute MOC control (smoothed squared DRNLoutput)     
    [MOCcontrol(filterCount,1:numSamples),  filterStates(17*(filterCount-1)+2)]=  filter(filterCoeffs(enumC_MOCb+1:enumC_MOCb+2), filterCoeffs(enumC_MOCa+1:enumC_MOCa+2) , y.^2, filterStates(17*(filterCount-1)+2));
    
    MOCcontrol(filterCount,1:numSamples) = sqrt(MOCcontrol(filterCount,1:numSamples));    % restore to meaningful scale (meters)
    
    % compare MOC control signal to specified threshold
    MOCcontrol(filterCount,1:numSamples)= MOCcontrol(filterCount,1:numSamples) - MOCthreshold(filterCount);    
    
end  % BF channel

% and zero all values below threshold (HWR)
MOCcontrol = max(MOCcontrol,0);

% default attenuation is 1 ( i.e. no MOC gives unit gain)
MOCcontrol= 1 + MOCcontrol*MOCfactor ;

% OUTPUT=> combine across all channels and restore to meaningful Pa
DRNLa1to1 = 1e4;
frameBuffer(1:numSamples) = frameBuffer(1:numSamples) / (DRNLa1to1*stapesScalar);
end

%nick modified broken stick function
function [x] = DRNL_brokenstick_nl (x, a, b, c)
% y = sign(x).* min(a*abs_x,  b*abs_x .^ c);
% This function could be replaced by a lookup table

abs_x = abs(x);
% linear (low amplitude) response
x=a.*x;

% compressed high amplitude
compressionThreshold=10.^((1/(1-c)).*log10(b./a));
% only values outside the compression threshold
%  need be subject to compression
idx=find(abs_x>compressionThreshold);
x(idx) = sign(x(idx)).* ( b*abs_x(idx) .^ c);

end %of DRNL_brokenstick_nl