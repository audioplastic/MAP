clear all; clc; clf
close all;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HA parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TAspl = [30; 30; 35; 50]; %Thresholds actual (in dB SPL)
TDspl = [0; 0; 0; 0]; %Thresholds desired (in dB SPL)
TCspl = TDspl + 50;       %Compression thresholds (in dB SPL)
TMspl = TCspl + 10;       %MOC thresholds (in dB SPL)
x.ARthresholddB = 80;     %Acoustic Reflex threshold (just 1 global value - now in dB SPL)


x.useAid = 1;
x.ARtau=0.03;             % 0.03

x.MOCfactor = 1e7;       % now that the conversions between velocity and displacement have been removed, internal values are 100 times greater (10kHz / 100 Hz) and so defaults need to be at least 100 times less than the old value of 2e9
x.MOCtau = 0.06;         % 0.06;

x.DRNLc=0.2;

x.channelBFs= [400 800 1600 3200];
x.bwOct     = [1 1 1 1]; %Octaves

x.filterOrder  = 4; %This sounds better than 2nd order - prevents anomolies seen in report also - I will work to make this happen in the hardware

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HA params - DO NOT EDIT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stapesScalar =6e-8;
DRNLa1to1 = 1e4;
x.DRNLaBaseline = DRNLa1to1 * 10.^((TAspl-TDspl)/20);  
x.DRNLb    =   x.DRNLaBaseline  .*  (stapesScalar * 2e-5 .* 10.^(TCspl/20)) .^ (1-x.DRNLc)  ;
x.MOCthreshold =  min([...
    x.DRNLaBaseline *  stapesScalar *  2e-5 .* 10.^(TMspl/20)...
    x.DRNLb        .* (stapesScalar .* 2e-5 .* 10.^(TMspl/20)  ) .^ (x.DRNLc)...
    ], [], 2); 

x.opScaling_dB = 20*log10(  1  /  (DRNLa1to1*stapesScalar)  );

%%
Levs = 0:25:100;

for jj = 1:numel(Levs)
    for kk = 1:numel(x.channelBFs)
        INdb  = Levs(jj);
        
        f = x.channelBFs(kk);
        
        dur = .4;
        fs = 50e3;
        dt = 1/fs;
        tAxis = dt:dt:dur;
        
        stimulusIN = sin(2*pi*f*tAxis);
        
        stimulusIN = stimulusIN./sqrt(mean(stimulusIN.^2)); %Normalize RMS to 1
        stimulusIN = stimulusIN * 2e-5*10^(INdb/20); %Convert RMS to pascals at desired level
        
        stimulusOUT = EssexAid_Fcn2(stimulusIN, fs, x);
        
        INdbx(jj,kk) = 20*log10(sqrt(mean(stimulusIN(1e3:end).^2)) / 2e-5)
        OUTdb(jj,kk) = 20*log10(sqrt(mean(stimulusOUT(1e3:end).^2)) / 2e-5)
        if kk ==1; figure; plot(stimulusOUT); title(num2str(Levs(jj))); ylim([-0.1 0.1]); drawnow; end
    end
end

%%
hearingAidGain = x.opScaling_dB;



inLo = 2e-5; %20e-6*10^(0/20);
inHi = 2.00; %20e-6*10^(100/20);
    
figure(33)
for kk = 1:numel(x.channelBFs)
    DRNLa = x.DRNLaBaseline(kk);%*(10^(x.DRNLaBaselineRelative_dB(kk)/20));
    DRNLb = x.DRNLb(kk);%*(10^(x.DRNLbRelative_dB(kk)/20));
    DRNLc = x.DRNLc;
    
%     if x.channelBFs(kk)<100
%         dispScalar = 1;
%     else        
%         dispScalar = 0.5^log2(x.channelBFs(kk)/lfRef); %re the 100 Hz cutoff
%     end
%     
%     if (x.doHpFilt)
%         hpFilterScalar= 0.5^log2(x.channelBFs(kk)/hfRef);%6 dB per octave
%     else
%         hpFilterScalar = 1;
%     end       
    
    
%     aLo(kk) = 20*log10(DRNLa * inLo * dispScalar * stapesScalar * (1/hpFilterScalar) /2e-5 ) + hearingAidGain;
%     aHi(kk) = 20*log10(DRNLa * inHi * dispScalar * stapesScalar * (1/hpFilterScalar) /2e-5 ) + hearingAidGain;
%     bLo(kk) = 20*log10(DRNLb * ((inLo * dispScalar * stapesScalar)^DRNLc) * (1/hpFilterScalar) /2e-5 ) + hearingAidGain;
%     bHi(kk) = 20*log10(DRNLb * ((inHi * dispScalar * stapesScalar)^DRNLc) * (1/hpFilterScalar) /2e-5 ) + hearingAidGain;
    
    aLo(kk) = 20*log10(DRNLa * inLo * stapesScalar /2e-5 ) + hearingAidGain;
    aHi(kk) = 20*log10(DRNLa * inHi * stapesScalar /2e-5 ) + hearingAidGain;
    bLo(kk) = 20*log10(DRNLb * ((inLo * stapesScalar)^DRNLc)  /2e-5 ) + hearingAidGain;
    bHi(kk) = 20*log10(DRNLb * ((inHi * stapesScalar)^DRNLc)  /2e-5 ) + hearingAidGain;
    
    
    subplot(1,4,kk)
    plot(Levs,OUTdb(:,kk),'x')
    hold on; plot([Levs(1) Levs(end)], [aLo(kk) aHi(kk)], 'k')
    hold on; plot([Levs(1) Levs(end)], [bLo(kk) bHi(kk)], 'k')
    xlabel('Input level in dB SPL')
    ylabel('Output level in dB SPL assuming op gain is set correctly')
    title(['cf=' num2str(x.channelBFs(kk)/1000) 'kHz'])
    ylim([0 100]); 
    xlim([0 100])
    hold on; plot([Levs(1) Levs(end)], [TAspl(kk) TAspl(kk)], ':g')
    cmpLev = 0;
    cmpThreshold(kk) = 20*log10(  (  (DRNLa/DRNLb)^(1/(DRNLc-1))  ) * 1/(2e-5*stapesScalar) );
    hold on; plot([cmpThreshold(kk) cmpThreshold(kk)], [0 100], ':r')
    ipLev = ( 2e-5  *  10^( (TAspl(kk)-hearingAidGain)/20 )  )  /  (DRNLa * stapesScalar  );
    ipThreshold(kk) = 20*log10(ipLev/2e-5);
    hold on; plot([ipThreshold(kk) ipThreshold(kk)], [0 100], ':g')
    
    MOCa(kk) = 20*log10((x.MOCthreshold(kk) / (DRNLa*stapesScalar) ) / 2e-5);
    MOCb(kk) = 20*log10( ( (  (x.MOCthreshold(kk) / DRNLb)^(1/DRNLc) ) / (stapesScalar) ) / 2e-5);
    ipMOCthr(kk) = max(MOCa(kk),MOCb(kk));
    hold on; plot([ipMOCthr(kk) ipMOCthr(kk)], [0 100], ':r')
    hold on; plot([x.ARthresholddB x.ARthresholddB], [0 100], ':y')
    
    hold on; plot([0 100], [0 100], ':k', 'LineWidth', 2)
    
%     hold on; plot([0 100], [TA-TDspl-MOCattLo(kk) 100], ':k');
    
end

