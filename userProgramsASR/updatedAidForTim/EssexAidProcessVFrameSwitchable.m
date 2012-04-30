%#eml
function [ frameBufferL, frameBufferR, filterStatesL, filterStatesR,  ARampL, ARampR, MOCend, peakIPL, peakOPL, rmsIPL, rmsOPL, peakIPR, peakOPR, rmsIPR, rmsOPR, MOCcontrol ] =...
    EssexAidProcessVFrameSwitchable( ...
    frameBufferL,...
    frameBufferR,...
    filterStatesL,...
    filterStatesR,...
    filterCoeffs,...
    numChannels,...
    numSamples,...
    ARampL,...
    ARampR,...
    ARthresholdPa,...
    filterOrder,...
    DRNLb,...
    DRNLc,...
    MOCthreshold_dBOP,...
    MOCfactor,...
    peakIPL,...
    peakOPL,...
    rmsIPL,...
    rmsOPL,...
    peakIPR,...
    peakOPR,...
    rmsIPR,...
    rmsOPR,...
    MOCend,...
    MOCcontrol,...
    mainGain,...
    useGTF)

%ESSEXAIDPROCESSFRAME Essex aid algorithm in frame processing mode
%   This code will look a bit odd to most Matlab programmers. This is
%   because the intended target is a C function that will be called on a
%   sub-millisecond basis. The bizzare enumerations assist the
%   pass-by-reference functionality that allows this function to fly in
%   real-time. This function works on a need to know basis, eliminating any
%   unnecessary data copying or parameter calculation.

% eml.varsize('frameBuffer', 6192);

%Fake enumeration - must be kept up to date with c++ enum
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
% Place where conversion to velocity once lived
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% tympanic membrane response in meters
stapesVelL = frameBufferL(1:numSamples);
stapesVelR = frameBufferR(1:numSamples);
frameBufferL(1:numSamples) = zeros(size(frameBufferL(1:numSamples)));
frameBufferR(1:numSamples) = zeros(size(frameBufferL(1:numSamples)));

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ACOUSTIC REFLEX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find  rms of smoothed ip signal
%  this will be used to trigger the AR reflex
% % % [yL,  filterStatesL(enumS_AR+1)] = filter(filterCoeffs(enumC_ARb+1:enumC_ARb+2), filterCoeffs(enumC_ARa+1:enumC_ARa+2) , stapesVelL.^2, filterStatesL(enumS_AR+1));
% % % [yR,  filterStatesR(enumS_AR+1)] = filter(filterCoeffs(enumC_ARb+1:enumC_ARb+2), filterCoeffs(enumC_ARa+1:enumC_ARa+2) , stapesVelR.^2, filterStatesR(enumS_AR+1));
% % % 
% % % % restore Pa scale
% % % smoothedARrmsL = sqrt(yL);  %confusing name for parameter - it is a short term RMS.
% % % smoothedARrmsR = sqrt(yR);  %confusing name for parameter - it is a short term RMS.

% attenuate input (NB cross product used)
stapesVelL = stapesVelL./ARampL(1:numSamples);
stapesVelR = stapesVelR./ARampR(1:numSamples);

% % % %CALC ARamp FOR NEXT FRAME
% % % % compare levels in the previous segment with AR threshold
% % % ARampL(1:numSamples) = smoothedARrmsL/ARthresholdPa;
% % % ARampR(1:numSamples) = smoothedARrmsR/ARthresholdPa;
% % % % all sub-treshold values are set to 1
% % % ARampL(ARampL(1:numSamples)<1)=1;
% % % ARampR(ARampR(1:numSamples)<1)=1;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DRNL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for filterCount=1:numChannels
    yL=stapesVelL;
    yR=stapesVelR;
    MOCnow = MOCcontrol(filterCount, 1:numSamples); %syntactic shorthand
    
    if useGTF
        assert(filterOrder < 5, 'filterOrder must be less than 5')
        for nn = 1:filterOrder;
            [yL, filterStatesL(17*(filterCount-1)+3+(nn-1)*2   :  17*(filterCount-1)+4+(nn-1)*2)] =  filter(filterCoeffs(10*(filterCount-1)+9:10*(filterCount-1)+10), filterCoeffs(10*(filterCount-1)+14:10*(filterCount-1)+16), yL, filterStatesL(17*(filterCount-1)+3+(nn-1)*2   :  17*(filterCount-1)+4+(nn-1)*2));
            [yR, filterStatesR(17*(filterCount-1)+3+(nn-1)*2   :  17*(filterCount-1)+4+(nn-1)*2)] =  filter(filterCoeffs(10*(filterCount-1)+9:10*(filterCount-1)+10), filterCoeffs(10*(filterCount-1)+14:10*(filterCount-1)+16), yR, filterStatesR(17*(filterCount-1)+3+(nn-1)*2   :  17*(filterCount-1)+4+(nn-1)*2));
        end
    else
        assert(filterOrder == 2 || filterOrder == 4, 'filterOrder must be 2 or 4')
        for nn = 1:filterOrder/2;
            [yL, filterStatesL(17*(filterCount-1)+3+(nn-1)*4   :  17*(filterCount-1)+6+(nn-1)*4)] =  filter(filterCoeffs(10*(filterCount-1)+9:10*(filterCount-1)+13), filterCoeffs(10*(filterCount-1)+14:10*(filterCount-1)+18), yL, filterStatesL(17*(filterCount-1)+3+(nn-1)*4   :  17*(filterCount-1)+6+(nn-1)*4));
            [yR, filterStatesR(17*(filterCount-1)+3+(nn-1)*4   :  17*(filterCount-1)+6+(nn-1)*4)] =  filter(filterCoeffs(10*(filterCount-1)+9:10*(filterCount-1)+13), filterCoeffs(10*(filterCount-1)+14:10*(filterCount-1)+18), yR, filterStatesR(17*(filterCount-1)+3+(nn-1)*4   :  17*(filterCount-1)+6+(nn-1)*4));
        end
    end
    
    peakIPL(filterCount) = 20*log10(   max(abs(yL)) /2e-5   ); %peak in in dB SPL
    peakIPR(filterCount) = 20*log10(   max(abs(yR)) /2e-5   ); %peak in in dB SPL
    
%     rmsIPL(filterCount) = 20*log10(  sqrt(mean(yL.^2)) /2e-5 ); %rms in in dB SPL
%     rmsIPR(filterCount) = 20*log10(  sqrt(mean(yR.^2)) /2e-5 ); %rms in in dB SPL

    rmsIPL(filterCount) = 20*log10(  abs(yL(end)) /2e-5 ); %rms in in dB SPL for GUI (bit of a hack, but it is smoothed in GUI)
    rmsIPR(filterCount) = 20*log10(  abs(yR(end)) /2e-5 ); %rms in in dB SPL  
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compression time
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % MOC compression first
    yL = yL.*MOCnow;
    yR = yR.*MOCnow;
    
    % Instantaneous after
    DRNLa = 1;
    yL = DRNL_brokenstick_nl (yL, DRNLa, DRNLb(filterCount), DRNLc(filterCount));
    yR = DRNL_brokenstick_nl (yR, DRNLa, DRNLb(filterCount), DRNLc(filterCount));
    
    if useGTF
        for nn = 1:filterOrder;
            [yL, filterStatesL(17*(filterCount-1)+11+(nn-1)*2   :  17*(filterCount-1)+12+(nn-1)*2)] =  filter(filterCoeffs(10*(filterCount-1)+9:10*(filterCount-1)+10), filterCoeffs(10*(filterCount-1)+14:10*(filterCount-1)+16), yL, filterStatesL(17*(filterCount-1)+11+(nn-1)*2   :  17*(filterCount-1)+12+(nn-1)*2));
            [yR, filterStatesR(17*(filterCount-1)+11+(nn-1)*2   :  17*(filterCount-1)+12+(nn-1)*2)] =  filter(filterCoeffs(10*(filterCount-1)+9:10*(filterCount-1)+10), filterCoeffs(10*(filterCount-1)+14:10*(filterCount-1)+16), yR, filterStatesR(17*(filterCount-1)+11+(nn-1)*2   :  17*(filterCount-1)+12+(nn-1)*2));
        end
    else
        for nn = 1:filterOrder/2;
            [yL, filterStatesL(17*(filterCount-1)+11+(nn-1)*4   :  17*(filterCount-1)+14+(nn-1)*4)] =  filter(filterCoeffs(10*(filterCount-1)+9:10*(filterCount-1)+13), filterCoeffs(10*(filterCount-1)+14:10*(filterCount-1)+18), yL, filterStatesL(17*(filterCount-1)+11+(nn-1)*4   :  17*(filterCount-1)+14+(nn-1)*4));
            [yR, filterStatesR(17*(filterCount-1)+11+(nn-1)*4   :  17*(filterCount-1)+14+(nn-1)*4)] =  filter(filterCoeffs(10*(filterCount-1)+9:10*(filterCount-1)+13), filterCoeffs(10*(filterCount-1)+14:10*(filterCount-1)+18), yR, filterStatesR(17*(filterCount-1)+11+(nn-1)*4   :  17*(filterCount-1)+14+(nn-1)*4));
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Deal with MOC control signal for next frame
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    yLdB = 20*log10(  ((2^-16) + abs(yL)) / 20e-6  ); %dB peak level with 16 bit offset
    yRdB = 20*log10(  ((2^-16) + abs(yR)) / 20e-6  ); %dB peak level with 16 bit offset
    
    [MOCnow,  filterStatesL(17*(filterCount-1)+2)]=  filter(filterCoeffs(enumC_MOCb+1:enumC_MOCb+2), filterCoeffs(enumC_MOCa+1:enumC_MOCa+2) , mean([yLdB+6; yRdB+6]), filterStatesL(17*(filterCount-1)+2));    
%     [MOCnow,  filterStatesL(17*(filterCount-1)+2)]=  filter(filterCoeffs(enumC_MOCb+1:enumC_MOCb+2), filterCoeffs(enumC_MOCa+1:enumC_MOCa+2) , mean([yL.^2; yR.^2]), filterStatesL(17*(filterCount-1)+2));
%     MOCnow = sqrt(MOCnow);    % restore to meaningful scale (meters - not anymore now in velocity mode)

    MOCnow = MOCnow - MOCthreshold_dBOP(filterCount);
%     MOCnow = 20*log10( MOCnow /  20e-6 ) - MOCthreshold_dBOP(filterCount);

    MOCnow = max(MOCnow, 0) * MOCfactor + 0.000000001; %the tiny offset is due to a crummy Embedded Matlab bug (This caused much misery)!
    MOCnow = 10.^(-MOCnow/20);
    
    MOCcontrol(filterCount, 1:numSamples) = MOCnow; %Replace the shortened version    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SOme calculations for the GUI monitoring
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    MOCend(filterCount) = MOCcontrol(filterCount, numSamples);
    
    peakOPL(filterCount) = 20*log10( max(abs(yL)) /2e-5 ); 
    peakOPR(filterCount) = 20*log10( max(abs(yR)) /2e-5 ); 
    
%     rmsOPL(filterCount) = 20*log10( sqrt(mean(yL.^2)) /2e-5  ); 
%     rmsOPR(filterCount) = 20*log10( sqrt(mean(yR.^2)) /2e-5 ); 

    rmsOPL(filterCount) = yLdB(end); 
    rmsOPR(filterCount) = yRdB(end);
%     rmsOPL(filterCount) = 20*log10( abs(yL(end)) /2e-5 ); 
%     rmsOPR(filterCount) = 20*log10( abs(yR(end)) /2e-5 ); 
        
    frameBufferL(1:numSamples) = frameBufferL(1:numSamples) + yL*mainGain(filterCount);
    frameBufferR(1:numSamples) = frameBufferR(1:numSamples) + yR*mainGain(filterCount);
end  % BF channel

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ACOUSTIC REFLEX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find  rms of smoothed ip signal
%  this will be used to trigger the AR reflex
[yL,  filterStatesL(enumS_AR+1)] = filter(filterCoeffs(enumC_ARb+1:enumC_ARb+2), filterCoeffs(enumC_ARa+1:enumC_ARa+2) , frameBufferL(1:numSamples).^2, filterStatesL(enumS_AR+1));
[yR,  filterStatesR(enumS_AR+1)] = filter(filterCoeffs(enumC_ARb+1:enumC_ARb+2), filterCoeffs(enumC_ARa+1:enumC_ARa+2) , frameBufferR(1:numSamples).^2, filterStatesR(enumS_AR+1));

% restore Pa scale
smoothedARrmsL = sqrt(yL);  %confusing name for parameter - it is a short term RMS.
smoothedARrmsR = sqrt(yR);  %confusing name for parameter - it is a short term RMS.

% % % % attenuate input (NB cross product used)
% % % stapesVelL = stapesVelL./ARampL(1:numSamples);
% % % stapesVelR = stapesVelR./ARampR(1:numSamples);

%CALC ARamp FOR NEXT FRAME
% compare levels in the previous segment with AR threshold
ARampL(1:numSamples) = smoothedARrmsL/ARthresholdPa;
ARampR(1:numSamples) = smoothedARrmsR/ARthresholdPa;
% all sub-treshold values are set to 1
ARampL(ARampL(1:numSamples)<1)=1;
ARampR(ARampR(1:numSamples)<1)=1;

end

%nick modified broken stick function
function [x] = DRNL_brokenstick_nl (x, a, b, c)
% y = sign(x).* min(a*abs_x,  b*abs_x .^ c);
% This function could be replaced by a lookup table

abs_x = abs(x);
% linear (low amplitude) response
x=a*x;

% compressed high amplitude
compressionThreshold=10.^((1/(1-c)).*log10(b/a));
% only values outside the compression threshold
%  need be subject to compression
idx=find(abs_x>compressionThreshold);
x(idx) = sign(x(idx)).* ( b*abs_x(idx) .^ c);

end %of DRNL_brokenstick_nl

