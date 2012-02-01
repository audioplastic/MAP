%% test2toneSuppression
% Systematic demonstration of two-tone suppressionuppression
%  at BM and AN levels.
%  Use 'demoTwoTone' instead for a simpler, AN only, probe probe level demo
%
% A probe tone is played at a probe level and a second tone is introduced
%  half way through the presentation. The response to the combined signal
%  is recorded and analysed.
%
% The probe frequency is set at the BF of the BM location (AN fiber) being
%  studied (important). All measurements are made at this location.
%
% The second tone called the seeep tone is presented at a reange of
%  frequencies and levels. In a linear system we might expect the addition
%  of a second tone to increase the magnitude of the output.
%  In fact, it often results in a *reduction* in the response.
%  This is called two-tone suppression.
%
% The effect of the sweep tone is represented in a countour plot showing
%  both reductions and increases in response
%  relative to the response in quiet.
%
% The background colour in this plot is the response to the
%  probe (at BF) tone alone
%
% The test is based on the paradigm of Ruggero et al 1992 (Fig 2)
%  where they used a probe suppressor level and
%  varied the probe (at BF) level. this not the same as many AN studies
%  where the probe tone is probe and the suppressor (sweep) tone is varied.

global  ANprobRateOutput DRNLoutput
dbstop if error
restorePath=path;
addpath (['..' filesep 'MAP'],    ['..' filesep 'wavFileStore'], ...
    ['..' filesep 'utilities'])


% Probe tone
probeToneFrequency=8000;
probeToneLevelsdB=[45 :5: 90];
probeToneLevelsdB=[0 :10: 90];
% probeToneLevelsdB=[30];
nFixedToneLevels=length(probeToneLevelsdB);
probeToneDuration=.050;		      % seconds

% Suppressor tone (called 'sweep tone')
sweeptoneLevelsdB=[ 0 : 10:90];
% nSweepToneFrequencies=21;
% lowestSweepFrequency=probeToneFrequency/6;
% highestSweepFrequency=probeToneFrequency*3;
% sweepToneFrequencies=round(logspace(log10(lowestSweepFrequency), ...
%     log10(highestSweepFrequency), nSweepToneFrequencies));
% sampleSweepFrequency=10600;

% or
sweepToneFrequencies=10000;
nSweepToneFrequencies=1;
sampleSweepFrequency=sweepToneFrequencies;
highestSweepFrequency=sampleSweepFrequency;

startSilenceDuration=0.010;

figure(5), clf
figure(87), clf

% All set. Run test
% key channels for snapshots
[a BFchannel]=min((sweepToneFrequencies-probeToneFrequency).^2);
[a sampleChannel]=min((sweepToneFrequencies-sampleSweepFrequency).^2);

sampleRate= max(44100, 10*highestSweepFrequency);
dt=1/sampleRate; % seconds

sweepStartPTR=...
    round((startSilenceDuration + probeToneDuration/2)*sampleRate);

BF_BMresponse=zeros(length(sweepToneFrequencies), ...
    length(probeToneLevelsdB), length(sweeptoneLevelsdB));
startSilence= zeros(1,startSilenceDuration*sampleRate);

probeTonedBCount=0;
for probeTonedB=probeToneLevelsdB
    probeTonedBCount=probeTonedBCount+1;

    BMpeakResponse= zeros(length(sweeptoneLevelsdB),length(sweepToneFrequencies));
    ANpeakResponse= zeros(length(sweeptoneLevelsdB),length(sweepToneFrequencies));
    sweepTonedBCount=0;
    for sweepToneDB=sweeptoneLevelsdB
        sweepTonedBCount=sweepTonedBCount+1;
        suppFreqCount=0;
        for sweepToneFrequency=sweepToneFrequencies
            suppFreqCount=suppFreqCount+1;

            % probeTone tone
            time1=dt: dt: probeToneDuration;
            amp=10^(probeTonedB/20)*28e-6;
            inputSignal=amp*sin(2*pi*probeToneFrequency*time1);
            rampDuration=.005; rampTime=dt:dt:rampDuration;
            ramp=[0.5*(1+cos(2*pi*rampTime/(2*rampDuration)+pi)) ones(1,length(time1)-length(rampTime))];
            inputSignal=inputSignal.*ramp;
            inputSignal=inputSignal.*fliplr(ramp);
            nsignalPoints=length(inputSignal);
            sweepStart=round(nsignalPoints/2);
            nSweepPoints=nsignalPoints-sweepStart;

            % sweepTone
            time2= dt: dt: dt*nSweepPoints;
            % B: tone on
            amp=10^(sweepToneDB/20)*28e-6;
            inputSignal2=amp*sin(2*pi*sweepToneFrequency*time2);
            rampDuration=.005; rampTime=dt:dt:rampDuration;
            ramp=[0.5*(1+cos(2*pi*rampTime/(2*rampDuration)+pi)) ones(1,length(time2)-length(rampTime))];
            inputSignal2=inputSignal2.*ramp;
            inputSignal2=inputSignal2.*fliplr(ramp);
            % add tone and sweepTone components
            inputSignal(sweepStart+1:end)= inputSignal(sweepStart+1:end)+inputSignal2;

            inputSignal=...
                [startSilence inputSignal ];

            %% run MAP
                MAPparamsName='Normal';
                AN_spikesOrProbability='probability';
                % only use HSR fibers (NB no acoustic reflex)
                paramChanges={'IHCpreSynapseParams.tauCa=80e-6;'};
                %                 paramChanges={'IHCpreSynapseParams.tauCa=80e-6;',...
                %                                 'DRNLParams.g=0;'};
                MAP1_14(inputSignal, sampleRate, probeToneFrequency, ...
                    MAPparamsName, AN_spikesOrProbability, paramChanges);

                % find toneAlone response
                if sweepTonedBCount==1 && suppFreqCount==1
                    BMprobeToneAloneRate=...
                        mean(abs(DRNLoutput(sweepStartPTR:end)));
                    ANprobeToneAloneRate=...
                        mean(abs(ANprobRateOutput(sweepStartPTR:end)));
                end
                BF_BMresponse(suppFreqCount,probeTonedBCount, ...
                    sweepTonedBCount)=...
                    mean(abs(DRNLoutput(sweepStartPTR:end)));

                BMpeakResponse(sweepTonedBCount,suppFreqCount)=...
                    mean(abs(DRNLoutput(sweepStartPTR:end)))...
                    /BMprobeToneAloneRate;
                ANpeakResponse(sweepTonedBCount,suppFreqCount)=...
                    mean(abs(ANprobRateOutput(sweepStartPTR:end)))...
                    /ANprobeToneAloneRate;
                disp(['F2: ', num2str([sweepToneFrequency sweepToneDB ...
                    BMpeakResponse(sweepTonedBCount,suppFreqCount)...
                    ANpeakResponse(sweepTonedBCount,suppFreqCount)])...
                    ' dB'])

                figure(5)
                time=dt:dt:dt*length(inputSignal);
                subplot(3,1,1), plot(time, inputSignal)
                title(['stimulus: Suppressor=' ...
                    num2str([sweepToneFrequency, sweepToneDB]) ' Hz/ dB'])

                time=dt:dt:dt*length(DRNLoutput);
                subplot(3,1,2), plot(time, DRNLoutput)
                title('BM displacement')
                xlim([0 probeToneDuration])

                time=dt:dt:dt*length(ANprobRateOutput);
                subplot(3,1,3), plot(time, ANprobRateOutput)
                xlim([0 probeToneDuration])
                ylim([0 500])
                title(['ANresponse: probeTone' num2str([probeToneFrequency, probeTonedB]) ' Hz/ dB'])

                if nSweepToneFrequencies>1
                    subplot(3,2,5)
                    contourf(sweepToneFrequencies,sweeptoneLevelsdB,BMpeakResponse, ...
                        [.1:.1:.9 1.05] )
                    set(gca,'Xscale','log')
                    set(gca,'Xtick', [1000  4000],'xticklabel',{'1000', '4000'})
                    set(gcf, 'name',['probeToneFrequency= ' num2str(probeToneFrequency)])
                    title(['BM' num2str(probeTonedB) ' dB'])

                    subplot(3,2,6)
                    contourf(sweepToneFrequencies,sweeptoneLevelsdB,ANpeakResponse, ...
                        [.1:.1:.9 1.05] )
                    set(gca,'Xscale','log')
                    set(gca,'Xtick', [1000  4000],'xticklabel',{'1000', '4000'})
                    set(gcf, 'name',['probeToneFrequency= ' num2str(probeToneFrequency)])
                    title(['AN:' num2str(probeTonedB) ' dB'])
                    drawnow
                end
        end
    end

    %% plot matrix

    if nSweepToneFrequencies>1

        figure (87)
        subplot(3, nFixedToneLevels, probeTonedBCount)
        surf(sweepToneFrequencies,sweeptoneLevelsdB,BMpeakResponse)
        set(gca,'Xscale','log')
        zlabel('gain')
        xlim([lowestSweepFrequency highestSweepFrequency])
        ylim([min(sweeptoneLevelsdB) max(sweeptoneLevelsdB)])
        title('BM response')
        view([-11 52])

        subplot(3, nFixedToneLevels, nFixedToneLevels+probeTonedBCount)
        contourf(sweepToneFrequencies,sweeptoneLevelsdB,BMpeakResponse, ...
            [.1:.5:.9 0.99 1.05] )
        hold on
        plot(probeToneFrequency, probeTonedB, 'or', 'markerfacecolor','w')
        set(gca,'Xscale','log')
        set(gca,'Xtick', [1000  5000],'xticklabel',{'1000', '5000'})
        ylabel('(BM) sweep level')
        xlabel('(BM) sweep freq')
        title(['probe tone level=' num2str(probeTonedB) ' dB'])
        %     colorbar

        subplot(3, nFixedToneLevels, 2*nFixedToneLevels+probeTonedBCount)
        contourf(sweepToneFrequencies,sweeptoneLevelsdB,ANpeakResponse, ...
            [.1:.5:.9 0.99 1.05] )
        hold on
        plot(probeToneFrequency, probeTonedB, 'or', 'markerfacecolor','w')
        set(gca,'Xscale','log')
        set(gca,'Xtick', [1000  5000],'xticklabel',{'1000', '5000'})
        ylabel('(AN) sweep level')
        xlabel('(AN) sweep freq')
        title(['probe tone level=' num2str(probeTonedB) ' dB'])
        %     colorbar
    end
end

%% Ruggero fig 2 (probe tone level is x-axis, sweep tone level is y-axis
%          BF_BMresponse=zeros(length(sweepToneFrequencies), ...
%     length(probeToneLevelsdB), length(sweeptoneLevelsdB));
%

figure(1),semilogy(probeToneLevelsdB,...
    squeeze(BF_BMresponse(sampleChannel,:,:)))
ylim([10e-9 10e-8])
legend(num2str(sweeptoneLevelsdB'),'location','southeast')
xlabel('probe SPL')
ylabel ('displacement (m)')
title(['Probe ' num2str(probeToneFrequency) ' Hz. Sweep ' ...
    num2str(sweepToneFrequencies(sampleChannel)) ' Hz.'])
path=restorePath;
