function [foregroundResults meanSummary]= ...
    plotProfile(location, fgName, bgName, figureNumber)
% plotProfile plots an auditory profile from a profile.m file
% 
% Input arguments:
%   location, is the relative file path to the folder containing the profile
%   fgName, is the full name of the file to be used
%   bgName, (optional) is the name of a comparison file to be plotted 
%          in the background, defaujlt is ''
%   figureNumber, MATLAB figure number default is #90
%
% Example:
%   plotProfile('allParticipants','profile_I10_L')
%   plotProfile('allParticipants','profile_I10_L', 'profile_N06_L', 1)
%
% Output arguments:
%   foregroundResults=[BFs', LongTone', ShortTone',  IFMCFreq', ...
%                      TMCfittedSlope', TMCFreq', fgDepth];
%                     NB TMCfittedSlope and fgDepth are the only new information
%    meanSummary=[meanAbsThreshold meanDepth meanFittedSlope];

restorePath=path;
addpath(['..\profiles\' location])
dbstop if error

%% plot profile
if nargin<4, figureNumber=90;   end
if nargin<3, bgName = '';       end
if nargin<1, error('plotProfile is not a script'), end

% fetch data from foreground file
cmd=['foreground = ' fgName ';'];
eval(cmd)

% fetch data from background file, if required
if nargin>=2 && ~isempty(bgName)
    cmd=['background = ' bgName ';'];
    eval(cmd)
end

%   EITHER
% figure(figureNumber), hold on  % multiple superimposed profiles
%   OR
figure(figureNumber), clf
set(gcf, 'name', fgName)

%% absolute thresholds absolute thresholds absolute thresholds
subplot(2,1,2)

% longTone
% remove NaNs for joined up plot
[x y]=stripNaNsfromPairedVariables(foreground.BFs,foreground.LongTone);
semilogx(x, y,'ko-','lineWidth',2,'markerSize', 2); hold on

% shortTone
% remove NaNs for joined up plot
[x y]=stripNaNsfromPairedVariables(foreground.BFs,foreground.ShortTone);
semilogx(x,y,'bo-','lineWidth',2,'markerSize', 2); hold on

meanAbsThreshold=mean(y);
    text(14000, 40, 'mean')
    text(14000, 30, 'absThr')
text(15000, 18, int2str(round(meanAbsThreshold)))

absThresholds=[x; y]; % for the record

% plot background thresholds
if ~isempty(bgName)
    % longTone
    % remove NaNs for joined up plot
    [x y]=stripNaNsfromPairedVariables(background.BFs,background.LongTone);
    semilogx(x, y,'k:'); hold on

    % shortTone
    [x y]=stripNaNsfromPairedVariables(background.BFs,background.ShortTone);
    semilogx(x, y,'k:'); hold on
end
ylim([0 100])


%% TMC   TMC   TMC   TMC   TMC   TMC   TMC   TMC   TMC   TMC   TMC
TMCfittedSlope=NaN(size(foreground.TMCFreq));

for BFno=1:length(foreground.TMCFreq) 
    subplot(2, length(foreground.TMCFreq)+1, BFno+1)

    x=foreground.TMCFreq(BFno);
    y=foreground.TMC(BFno,:);
    % remove NaNs for joined up plot
    [x y]=stripNaNsfromPairedVariables(foreground.Gaps,foreground.TMC(BFno,:));

    plot(1000*x,y,'o-b','lineWidth',2,'markerSize', 2), hold on
    ylim([-10 110]),     xlim([0 100])

    if BFno==1  % label only the first figure
        ylabel('masker dB SPL')
        xlabel('gap (ms)')
        title([num2str(foreground.TMCFreq(BFno)) ' Hz'])
    else
        set(gca,'YTickLabel',[])
        set(gca,'xTickLabel',[])
        title([num2str(foreground.TMCFreq(BFno))]) % NB no 'Hz'
    end

    % publish slope at the bottom of the chart
    if ~isempty(x)
        P=polyfit(x,y,1);
        TMCfittedSlope(BFno)=P(1)/10;
        text(40,10,num2str(TMCfittedSlope(BFno),'%4.0f'))
    end
end
meanFittedSlope=mean(TMCfittedSlope(~isnan(TMCfittedSlope)));
text(150, 48, num2str(meanFittedSlope,'%4.0f'))
    text(130, 70, 'mean')
    text(130, 60, 'slope')

% background (plot only)
if ~isempty(bgName)
    for BFno=1:length(background.TMCFreq)-1
        BF = background.TMCFreq(BFno);
        idx = find(BF == foreground.TMCFreq);
        if ~isempty(idx);
            subplot(2, length(foreground.TMCFreq)+1, BFno+1)
            [x y]=stripNaNsfromPairedVariables(background.Gaps,background.TMC(BFno,:));
            plot(1000*x,y,'k:')
        end

        % publish slope at the bottom of the chart
        if ~isempty(x)
            P=polyfit(x,y,1);
            TMCfittedSlope(BFno)=P(1)/10;
            text(40,0,num2str(TMCfittedSlope(BFno),'%4.0f'))
        end
    end
end



%% IFMCs   IFMCs   IFMCs   IFMCs   IFMCs   IFMCs   IFMCs   IFMCs
% fgDepth is height of wings above tip
tipThreshold=foreground.IFMCs(:,foreground.MaskerRatio==1);
fgDepth=    mean([foreground.IFMCs(:,2) foreground.IFMCs(:,6)], 2)...
    -tipThreshold;

for BFno=1:length(foreground.IFMCFreq)
    subplot(2,1,2)
    % convert from frequency ratio to frequency
    probeFreq=foreground.IFMCFreq(BFno);
    maskerFrequencies=foreground.MaskerRatio'* probeFreq;
    %     % multiple superimposed profiles
    %     foreground.IFMCs(BFno,:)=foreground.IFMCs(BFno,:)-tipThreshold(BFno)+30;
    % remove NaNs for joined up plots
    [x y]=stripNaNsfromPairedVariables(maskerFrequencies,foreground.IFMCs(BFno,:));
    semilogx(x,y,'o-r','lineWidth',2,'markerSize', 2), hold on

    % white circles for probe frequency
    if ~isempty(tipThreshold) && ~isnan(tipThreshold(BFno))
        plot(probeFreq, tipThreshold(BFno),'ro','markerFaceColor','w')
        if meanAbsThreshold<50
            % print at the top of the plot
            text(probeFreq,90,num2str(fgDepth(BFno),'%5.0f'))
        else
            text(probeFreq,10,num2str(fgDepth(BFno),'%5.0f'))
        end
    end

    ylim([-10 100]),  xlim([100 12000])
    set(gca,'XScale','log') % force log even when plot is empty
    set(gca,'xTick',[250; 500; 1000; 2000; 4000; 8000]);
    set(gca,'xTickLabel',{'250', '500', '1000', '2000', '4000', '8000'});
    meanDepth=mean(fgDepth(~isnan(fgDepth)));
    text(14000, 80, 'mean')
    text(14000, 70, 'depth')
    text(15000, 58, num2str(meanDepth, '%4.0f'))
    %     grid on
end
xlabel('frequency (Hz)')
ylabel(' dB SPL')

% background plot only
if ~isempty(bgName)
    tipThreshold=background.IFMCs(:,background.MaskerRatio==1);
    bgDepth= mean([background.IFMCs(:,2) background.IFMCs(:,6)], 2)...
        -tipThreshold;
    for BFno=1:length(background.IFMCFreq)
        % convert from frequency ratio to frequency
        probeFreq=background.IFMCFreq(BFno);
        maskerFrequencies=background.MaskerRatio'*background.IFMCFreq(BFno);
        subplot(2,1,2)
        [x y]=stripNaNsfromPairedVariables(maskerFrequencies,background.IFMCs(BFno,:));
        semilogx(x,y,'k:')
        if meanAbsThreshold<50
            % print at the top of the plot
            text(probeFreq,80,num2str(bgDepth(BFno),'%5.0f'))
        else
            text(probeFreq,0,num2str(bgDepth(BFno),'%5.0f'))
        end
    end
end

set(get(gca,'title'),'interpreter','None') % no funny characters
ttl=fgName(9:end);
% add second file name
if ~isempty(bgName),  ttl=[ttl ' / ' bgName(9:end)]; end
title(ttl,'HorizontalAlignment','Left')

% collect together results in a matrix ready for return
foregroundResults=[ foreground.BFs', foreground.LongTone', ...
    foreground.ShortTone',  foreground.TMCFreq', TMCfittedSlope', ...
     foreground.IFMCFreq', fgDepth];
meanSummary=[meanAbsThreshold  meanFittedSlope meanDepth];
pause(.1)

path(restorePath);


function [a b]=stripNaNsfromPairedVariables(a,b)
idx=find(~isnan(b)); a=a(idx); b=b(idx);

