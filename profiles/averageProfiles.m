% average profiles
clear all, dbstop if error, figure(10),clf

profileType='N';
averageProfileFileName='averageNormal';

load participantCompendium

restorePath=path;
addpath(['..' filesep 'multiThreshold 1.46'])

noParticipants=length(participant);
allLongtoneThresholds=NaN(1,2*noParticipants);
allShortToneThresholds=NaN(1,2*noParticipants);

probeFrequencies=[250 500 1000 2000 4000 6000];
probeNo=0;
leftEar=participant.leftEarData;
TMCgaps=leftEar.Gaps;
longToneMean=NaN(1,length(probeFrequencies));
shortToneMean=NaN(1,length(probeFrequencies));
longTonestdev=NaN(1,length(probeFrequencies));
shortTonestdev=NaN(1,length(probeFrequencies));
meanDepth=NaN(1,length(probeFrequencies));
meanIFMC=NaN(length(probeFrequencies), 7);
meanTMC=NaN(length(probeFrequencies),length(TMCgaps));
stdevIFMC=NaN(length(probeFrequencies), 7);
stdevTMC=NaN(length(probeFrequencies),length(TMCgaps));
TMCfittedSlope=NaN(1,length(probeFrequencies));
sampleSize=NaN(1,length(probeFrequencies));
for  probeFrequency=probeFrequencies
    probeNo=probeNo+1;

    allIFMCreTip=NaN(2*noParticipants,7);
    allTMCs=NaN(2*noParticipants,9);
    tipThresholds=NaN(2*noParticipants,7);
    depths=NaN(2*noParticipants,7);

    earCount=0;
    for participantNo=1:noParticipants
        code=participant(participantNo).code;
        if strcmp(code(1),profileType)
            leftEar=participant(participantNo).leftEarData;
            IFMCs=leftEar.IFMCs;
            IFMCprobeFrequencies=leftEar.IFMCFreq;
            TMCs=leftEar.TMC;

            probePTR=find(IFMCprobeFrequencies==probeFrequency);
            IFMCatProbe=IFMCs(probePTR,:);
            TMCatProbe=TMCs(probePTR,:);
            longToneThreshold=leftEar.LongTone(probePTR);
            allLongtoneThresholds(earCount)=longToneThreshold;
            shorttoneThreshold=leftEar.ShortTone(probePTR);
            allShortToneThresholds(earCount)=shorttoneThreshold;
            allTMCs(earCount,:)=TMCatProbe;

            tipThresholds(earCount)=IFMCatProbe(4);
            % relative to tip frequency
            lIFMCreTip=IFMCatProbe-tipThresholds(earCount);
            allIFMCreTip(earCount,:)=lIFMCreTip;
            depths(earCount)=mean([lIFMCreTip(2) lIFMCreTip(6)]);
            %
            earCountR=earCount+1;
            rightEar=participant(participantNo).rightEarData;
            IFMCs=rightEar.IFMCs;
            IFMCprobeFrequencies=rightEar.IFMCFreq;
            IFMCmaskerFrequencies=probeFrequency*rightEar.MaskerRatio;
            TMCs=leftEar.TMC;
            TMCgaps=leftEar.Gaps;

            probePTR=find(IFMCprobeFrequencies==probeFrequency);
            IFMCatProbe=IFMCs(probePTR,:);
            TMCatProbe=TMCs(probePTR,:);
            longToneThreshold=rightEar.LongTone(probePTR);
            allLongtoneThresholds(earCountR)=longToneThreshold;
            shorttoneThreshold=rightEar.ShortTone(probePTR);
            allShortToneThresholds(earCountR)=shorttoneThreshold;
            allTMCs(earCountR,:)=TMCatProbe;

            tipThresholds(earCountR)=IFMCatProbe(4);
            % relative to tip frequency
            lIFMCreTip=IFMCatProbe-tipThresholds(earCountR);
            allIFMCreTip(earCountR,:)=lIFMCreTip;
            depths(earCountR)=mean([lIFMCreTip(2) lIFMCreTip(6)]);

        end
        earCount=earCount+2;
    end % participant

    % compute means ignoring missing data
    for i=1:7
        x=allIFMCreTip(:,i); x=x(~isnan(x));
     meanIFMC(probeNo,i)=mean(x);
     stdevIFMC(probeNo,i)=std(x);
    end
    for i=1:9
        x=allTMCs(:,i); x=x(~isnan(x));
        meanTMC(probeNo,i)=mean(x);
        stdevTMC(probeNo,i)=std(x);
    end

    x=allIFMCreTip(:,4); x=x(~isnan(x));
    sampleSize(probeNo)=length(x);

    meanAllLongToneAtProbe=mean(allLongtoneThresholds(~isnan(allLongtoneThresholds)));
    stdAllLongToneAtProbe=std(allLongtoneThresholds(~isnan(allLongtoneThresholds)));
    meanAllShortToneAtProbe=mean(allShortToneThresholds(~isnan(allShortToneThresholds)));
    stdAllShortToneAtProbe=std(allShortToneThresholds(~isnan(allShortToneThresholds)));
    meanTipThresholds=mean(tipThresholds(~isnan(tipThresholds)));
    meanIFMC(probeNo,:)=meanIFMC(probeNo,:) + meanTipThresholds;
    stdevIFMC(probeNo,:)=std(meanIFMC(probeNo,:) );
    
    figure(10),subplot(2,1,2)
    plot(IFMCmaskerFrequencies,meanIFMC(probeNo,:) )
    hold on

    figure(10), subplot(2, length(probeFrequencies)+1,probeNo+1)
    errorbar(TMCgaps,meanTMC(probeNo,:),stdevTMC(probeNo,:),'o'), ylim([0 100])
    P=polyfit(TMCgaps,meanTMC(probeNo,:),1);
    TMCfittedSlope(probeNo)=P(1)/10;
    text(0.01,10,num2str(TMCfittedSlope(probeNo),'%4.0f'))

    longToneMean(probeNo)=meanAllLongToneAtProbe;
    shortToneMean(probeNo)=meanAllShortToneAtProbe;
    longTonestdev(probeNo)=stdAllLongToneAtProbe;
    shortTonestdev(probeNo)=stdAllShortToneAtProbe;
    meanDepth(probeNo)=mean(depths(~isnan(depths)));
end % probe frequency

meanIFMC_STD=mean(mean(stdevIFMC))

figure(10),subplot(2,1,2),hold on
errorbar(probeFrequencies, longToneMean,longTonestdev)
figure(10),subplot(2,1,2),hold on
errorbar(probeFrequencies*1.1, shortToneMean, shortTonestdev)
ylim([0 100]), xlim([100 10000])
set(gca,'Xscale','log')
title([' average standard deviation for IFMC is ' num2str(meanIFMC_STD) ' dB'])

text(150,80,'N=')
for i=1:length(probeFrequencies)
    text(probeFrequencies(i),80,num2str(sampleSize(i)))
end

maskerRatios=leftEar.MaskerRatio;
profile2mFile(longToneMean', shortToneMean', TMCgaps', probeFrequencies, ...
    meanTMC', maskerRatios', ...
    meanIFMC', averageProfileFileName, 'allParticipants')
path(restorePath)
