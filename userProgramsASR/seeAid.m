clear all; clc; clf
% close all;


thresholds = [40 40 40 40];

TA = thresholds;        %Thresholds actual
TD = thresholds /2;    %Thresholds desired
TC = TD + 20;           %Compression threshold
stapesScalar =6e-8;
hfRef = 1e4;
lfRef = 1e2; 


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HA parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x.useAid = 0;



x.ARthresholddB=170;         % 60
x.ARtau=0.2;               % 0.03

x.DRNLaBaseline=1e4;        % 1e4
x.DRNLb=1;            % 0.5e-5 ?? --> YES but misleading. This = 5e-6 as in the readme (watch out fot the 0.5 vs 5)
x.DRNLc=0.2;

x.MOCfactor = 2e7;%2e9;            % 2e9
x.MOCtau = 0.06;              % 0.06;
x.MOCthreshold = 5e-10 * 10^(60/20);

x.channelBFs= [400 800 1600 3200];%[250 630 1587 4000];

x.bwOct                       = [1 1 1 1];%[4/3 4/3 4/3 4/3]; %Octaves
x.DRNLaBaselineRelative_dB    = TA - TD;

d = 0.5.^log2(x.channelBFs/lfRef);
x.DRNLbRelative_dB            =  20*log10( 10.^((TA-TD)/20)  .*  x.DRNLaBaseline  .*  (d.*stapesScalar*2e-5 .* 10.^(TC/20)) .^ (1-x.DRNLc)  );
x.filterOrder                 = 2; %now changed --> See experiment 6 in log book

x.doHpFilt                    = 1;    %TRUE


x.opScaling_dB = 20*log10(  1  /  (x.DRNLaBaseline*stapesScalar*(lfRef/hfRef))  );

%%
Levs = 0:10:100;




for jj = 1:numel(Levs)
    for kk = 1:numel(x.channelBFs)
        INdb  = Levs(jj);
        
        f = x.channelBFs(kk);
        
        dur = 0.4;
        fs = 50e3;
        dt = 1/fs;
        tAxis = dt:dt:dur;
        
        stimulusIN = sin(2*pi*f*tAxis);
        
        stimulusIN = stimulusIN./sqrt(mean(stimulusIN.^2)); %Normalize RMS to 1
        stimulusIN = stimulusIN * 2e-5*10^(INdb/20); %Convert RMS to pascals at desired level
        
        stimulusOUT = EssexAid_Fcn(stimulusIN, fs, x);
        
        
        OUTdb(jj,kk) = 20*log10(sqrt(mean(stimulusOUT(1e3:end).^2)) / 2e-5)
    end
end

%%
clf


hearingAidGain = x.opScaling_dB;



inLo = 2e-5; %20e-6*10^(0/20);
inHi = 2.00; %20e-6*10^(100/20);
    

for kk = 1:numel(x.channelBFs)
    DRNLa = x.DRNLaBaseline*(10^(x.DRNLaBaselineRelative_dB(kk)/20));
    DRNLb = x.DRNLb*(10^(x.DRNLbRelative_dB(kk)/20));
    DRNLc = x.DRNLc;
    
    if x.channelBFs(kk)<100
        dispScalar = 1;
    else        
        dispScalar = 0.5^log2(x.channelBFs(kk)/lfRef); %re the 100 Hz cutoff
    end
    
    if (x.doHpFilt)
        hpFilterScalar= 0.5^log2(x.channelBFs(kk)/hfRef);%6 dB per octave
    else
        hpFilterScalar = 1;
    end       
    
    
    aLo(kk) = 20*log10(DRNLa * inLo * dispScalar * stapesScalar * (1/hpFilterScalar) /2e-5 ) + hearingAidGain;
    aHi(kk) = 20*log10(DRNLa * inHi * dispScalar * stapesScalar * (1/hpFilterScalar) /2e-5 ) + hearingAidGain;
    bLo(kk) = 20*log10(DRNLb * ((inLo * dispScalar * stapesScalar)^DRNLc) * (1/hpFilterScalar) /2e-5 ) + hearingAidGain;
    bHi(kk) = 20*log10(DRNLb * ((inHi * dispScalar * stapesScalar)^DRNLc) * (1/hpFilterScalar) /2e-5 ) + hearingAidGain;
    
    
    subplot(1,4,kk)
    plot(Levs,OUTdb(:,kk),'x')
    hold on; plot([Levs(1) Levs(end)], [aLo(kk) aHi(kk)], ':k')
    hold on; plot([Levs(1) Levs(end)], [bLo(kk) bHi(kk)], ':k')
    xlabel('Input level in dB SPL')
    ylabel('Output level in dB SPL assuming op gain is set correctly')
    title(['cf=' num2str(x.channelBFs(kk)/1000) 'kHz'])
    ylim([0 100]); xlim([0 100])
    hold on; plot([Levs(1) Levs(end)], [thresholds(kk) thresholds(kk)], ':g')
    cmpLev = 0;
    cmpThreshold(kk) = 20*log10(  (  (DRNLa/DRNLb)^(1/(DRNLc-1))  ) * 1/(2e-5*stapesScalar*dispScalar) );
    hold on; plot([cmpThreshold(kk) cmpThreshold(kk)], [0 100], ':r')
    ipLev = ( 2e-5  *  10^( (thresholds(kk)-hearingAidGain)/20 )  )  /  (DRNLa * stapesScalar * (lfRef/hfRef) );
    ipThreshold(kk) = 20*log10(ipLev/2e-5);
    hold on; plot([ipThreshold(kk) ipThreshold(kk)], [0 100], ':r')
    
    MOCa(kk) = 20*log10((x.MOCthreshold / (DRNLa*d(kk)*stapesScalar) ) / 2e-5);
    MOCb(kk) = 20*log10( ( (  (x.MOCthreshold / DRNLb)^(1/DRNLc) ) / (stapesScalar*d(kk)) ) / 2e-5);
    ipMOCthr(kk) = max(MOCa(kk),MOCb(kk));
    hold on; plot([ipMOCthr(kk) ipMOCthr(kk)], [0 100], '-b')
    
    
    hold on; plot([0 100], [0 100], 'b', 'LineWidth', 4)
    
end

