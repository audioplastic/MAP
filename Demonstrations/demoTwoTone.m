% function demoTwoTone
%
% Demonstration of two-tone suppression in auditory nerve
%  A steady probe tone is presented and a second tone (suppressor)
%  is added half way through the presentation.
%  The suppressor level is gradually increased between trials
%  Between runs the suppressor frequency is incremented.
%
% The display (Fig. 6) shows (after each trial)
%  1. the original signal
%  2. the AN response in the channel corresponding to the probe frequency.
%  3. the channel x time matrix of the auditory nerve response
%  4. a contour plot of the peak AN response during the tone + suppressor
%      as a function of suppressor level and frequency.
%      A level of 1 corresponds to the probeAlone rate.
%      The white dot shows the frequency and level of the probe tone
%
% Fig. 8 shows
%  1. a 3D plot of the firing rate in response to both tones
%       relative to probeAlone rate
%  2. a contour plot of 2-tone inhibition as defined by
%       Abbas and Sachs (1976). This is the
%         2tone driven rate/ probe tone driven rate
%         (The driven rate is the rate -spontaneous rate).
%
% By altering the code below, the user can set the channel BFs in the
%  model, the probe frequency, probe level, suppressor frequencies and
%  levels,stimulus duration and the duration of the initial silence.
%   For best results probe and suppressor frequencies should have corresponding
%   channel BFs. More channels than stimulus frequencies is OK, however.
%
% Model parameter changes can also be made using the paramChanges cell
%  array of strings (see manual for more information on how to do this)
%

global  ANprobRateOutput IHCpreSynapseParams
restorePath=path;
addpath (['..' filesep 'MAP'],    ['..' filesep 'wavFileStore'], ...
    ['..' filesep 'utilities'])
dbstop if error

%% #1 number of channels in the model
%   21-channel model (log spacing)
lowestBF=250; 	highestBF= 8000; 	numChannels=21;
BFlist=round(logspace(log10(lowestBF), log10(highestBF), numChannels));
% This list can be set manually, e.g. BFlist=[250 500 1000...] etc.

%% #2 probe frequency
probeFrequency=2000;

%% #3 probe levels
probedB=15;

%% #4 suppressor levels
suppressorLevels=00:10:80;

%% #5 suppressor frequencies
% Make sure that suppressor frequencies correspond to available channels
suppressorFrequencies=BFlist;
BFchannel=find(BFlist==probeFrequency);
if isempty(BFchannel)
    error(['probeFrequency must be set to an existing channel BF. BFs = '...
        num2str(BFlist)])
end

%% #6 Other stimulus characteristics
duration=.050;		      % seconds
tone2Duration=duration/2; % s
startSilenceDuration=0.010;

%% #7 parameter changes

paramChanges={};
% example parameter change
% paramChanges={'DRNLParams.g=10;'};

%% #8 contours used in suppression summary plot
suppressionContours=[.1 :.1: .9, .95, 1.05,]; 

%%
sampleRate= 44100; % Hz (higher sample rate needed for BF>8000 Hz)
dt=1/sampleRate; % seconds

%% hand over to demonstration
% sit back and enjoy the movie
figure(6), clf
startSilence=...
    zeros(1,startSilenceDuration*sampleRate);
suppressorStartsPTR=round((duration/2)*sampleRate);
drivenRatePTR=      round((duration/4)*sampleRate);


twoToneRelativeResponse= zeros(length(suppressorLevels),length(suppressorFrequencies));
twoToneMeanResponse=zeros(size(twoToneRelativeResponse));
singleToneMeanResponse=zeros(size(twoToneRelativeResponse));
twoToneGain=zeros(size(twoToneRelativeResponse));

suppFreqCount=0;
for suppressorFrequency=suppressorFrequencies
    suppFreqCount=suppFreqCount+1;
    
    suppressorLevelCount=0;
    for suppressorDB=suppressorLevels
        pause(.5)
        suppressorLevelCount=suppressorLevelCount+1;
        
        % primary + suppressor
        suppressorLevelDB=suppressorDB;
        
        % primary BF tone
        time1=dt: dt: duration;
        amp=10^(probedB/20)*28e-6;
        inputSignal=amp*sin(2*pi*probeFrequency*time1);
        rampDuration=.005; rampTime=dt:dt:rampDuration;
        ramp=[0.5*(1+cos(2*pi*rampTime/(2*rampDuration)+pi)) ones(1,length(time1)-length(rampTime))];
        inputSignal=inputSignal.*ramp;
        inputSignal=inputSignal.*fliplr(ramp);
        
        % suppressor
        time2= dt: dt: tone2Duration;
        % B: tone on
        amp=10^(suppressorLevelDB/20)*28e-6;
        inputSignal2=amp*sin(2*pi*suppressorFrequency*time2);
        rampDuration=.005; rampTime=dt:dt:rampDuration;
        ramp=[0.5*(1+cos(2*pi*rampTime/(2*rampDuration)+pi)) ones(1,length(time2)-length(rampTime))];
        inputSignal2=inputSignal2.*ramp;
        inputSignal2=inputSignal2.*fliplr(ramp);
        
        % add tone and suppressor components
        inputSignal(suppressorStartsPTR: ...
            suppressorStartsPTR+length(inputSignal2)-1)= ...
            inputSignal(suppressorStartsPTR: ...
            suppressorStartsPTR+length(inputSignal2)-1)+ inputSignal2 ;
        
        % run MAP
        MAPparamsName='Normal';
        AN_spikesOrProbability='probability';
        MAP1_14(inputSignal, sampleRate, BFlist, ...
            MAPparamsName, AN_spikesOrProbability, paramChanges);
        
        % fetch AN response
        nFiberTypess=length(IHCpreSynapseParams.tauCa);
        if nFiberTypess==1
            % only HSR fibers are used
            response=ANprobRateOutput;
        elseif nFiberTypess==2
            [nRow nCols]=size(ANprobRateOutput);
            HSRvaluesStartAt=floor(nRow/2)+1;
            response=ANprobRateOutput(HSRvaluesStartAt:end,:);
        end
        
        
        twoToneMeanResponse(suppressorLevelCount,suppFreqCount)=...
            mean(response(BFchannel,suppressorStartsPTR:end));
        singleToneMeanResponse(suppressorLevelCount,suppFreqCount)=...
            mean(response(BFchannel,...
            drivenRatePTR: suppressorStartsPTR-1));
        % express this relative to probeAlone activity.
        if suppressorLevelCount==1 && suppFreqCount==1
            probeAlone=twoToneMeanResponse(1,1);
            spontaneousRate=mean(response(1, startSilenceDuration*sampleRate)) ;
        end
        
        % Abbas and Sachs (using driven rates)
        twoToneGain(suppressorLevelCount,suppFreqCount)=...
            (twoToneMeanResponse(suppressorLevelCount,suppFreqCount)-...
            spontaneousRate)/...
            (singleToneMeanResponse(suppressorLevelCount,suppFreqCount)-...
            spontaneousRate);
        
        twoToneRelativeResponse(suppressorLevelCount,suppFreqCount)=...
            twoToneMeanResponse(suppressorLevelCount,suppFreqCount)/...
            probeAlone;
        
        disp(['F2 level= ', num2str([suppressorFrequency suppressorDB ...
            twoToneRelativeResponse(suppressorLevelCount,suppFreqCount)])])
        
        % Fig. 6 (results trial by trial)
        % signal
        figure(6), subplot(6,1,1)
        plot(dt:dt:dt*length(inputSignal), inputSignal, 'k')
        title('signal:probe       with added           suppressor')
        ylim([-.1 .1]), ylabel('Pascals')
        xlim([0 duration])
        
        % AN response
        figure(6), subplot(6,1,2)
        probChannelResponse=response(BFchannel,:);
        plot(dt:dt:dt*length(probChannelResponse), ...
            probChannelResponse, 'k')
        title(['AN response at BF=probe frequency (' ...
            num2str(probeFrequency) ' Hz)'])
        ylabel('sp/s'),             ylim([0 500])
        xlim([0 duration])
        set(gca, 'xtick',[])
        
        figure(6), subplot(3,1,2)
        PSTHbinWidth=0.010;
        PSTH= UTIL_shrinkBins(response, dt, PSTHbinWidth);
        [nY nX]=size(PSTH);
        time=PSTHbinWidth*(0:nX-1);
        surf(time, BFlist, PSTH)
        zlim([0 500]),
        caxis([0 500])
        shading interp
        colormap(jet)
        set(gca, 'yScale','log')
        xlim([0 max(time)+dt])
        set(gca, 'ytick',[500 1000 2000 4000])
        ylim([0 max(BFlist)]), ylabel('Channel BF')
        view([0 90]) % view([-8 68])
        title(['AN firing rate,  probe: '...
            num2str(probedB) 'dB   Suppressor: '...
            num2str(suppressorFrequency) 'Hz,  '...
            num2str(suppressorLevelDB) 'dB'])
        pause(0.05)
    end
    
    % start contour plot
    figure(6),   subplot(3,1,3), cla
    contourf(suppressorFrequencies,suppressorLevels, twoToneRelativeResponse,...
        suppressionContours);
    hold on
    plot(probeFrequency, probedB,'ok','markerfacecolor','w')
    set(gca,'Xscale','log')
    set(gca,'Xtick', [500 1000 2000 4000],'xticklabel',{'500','1000',...
        '2000', '4000'})
    set(gcf, 'name',['probe tone frequency= ' ...
        num2str(probeFrequency)])
    title('AN suppression contour: 1=probeAlone rate')
    xlabel ('suppressor tone frequency'), ylabel ('suppressor tone dB')
    colorbar
    drawnow
    
end

fprintf('\n')
disp(['AN response during both tones relative to probeAlone rate (=' ...
    num2str(probeAlone,'%4.1f') ' spikes/s)'])
UTIL_printTabTable([[0 suppressorLevels]' [suppressorFrequencies; ...
    twoToneRelativeResponse]], [], '%6.2f')
% sweep the path
path=restorePath;

return

%% plot summaries (not included in main demo)

% 3d image of relative rates 2-tone rate/probe rate.
figure(5), subplot(2,1,1)
surf(suppressorFrequencies,suppressorLevels,twoToneRelativeResponse)
shading interp
view([-9 40])
set(gca,'Xscale','log')
xlim([min(suppressorFrequencies) max(suppressorFrequencies)])
set(gca,'Xtick', [1000  4000],'xticklabel',{'1000', '4000'})
xlabel('suppressor frequency')
ylim([min(suppressorLevels) max(suppressorLevels)])
ylabel(' level')
zlim([min(min(twoToneRelativeResponse)) max(max(twoToneRelativeResponse))])
title('AN two tone rate/ probeAlone rate')

% suppression measure after abbas and Sachs (relative *driven* rates)
figure(5), subplot(2,1,2)
contourf(suppressorFrequencies,suppressorLevels,twoToneGain)
set(gca,'Xscale','log')
set(gca,'Xtick', [1000  4000],'xticklabel',{'1000', '4000'})
set(gcf, 'name',['probeFrequency= ' num2str(probeFrequency)])
xlabel ('suppressor  frequency'), ylabel ('suppressor level dB SPL')
title([' 2tone driven/1tone driven.   Probe' num2str(probeFrequency) 'Hz  '...
    num2str(probedB) ' dB'])
colorbar

