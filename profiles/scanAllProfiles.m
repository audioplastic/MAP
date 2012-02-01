function meanSDsummary=scanAllProfiles(subfolder)
% scanAllProfiles reads every file in a folder
%  and plots it using 'plotProfile'.
% The folder should be a subfolder of 'profiles' in MAP1_14
%
% Each plot can be accompanied by a second reference plot based on
%  backgroundPlotName which can be '' to ignore
%
% It also summarises the results in 'allForegroundResults'
% A profile is generated for each ear examined in 'profile'
%  allForegroundResults & profile are both saved in a .mat file with the
%  same name as the original folder.
%
% This program normally resides in 'profiles' but can sit in any
%  subfolder of MAP1_14.
%
% scanAllProfiles('allParticipants')
% scanAllProfiles('allParticipants', 'profile_CMA_L')

% locate all files
folder='profiles';
if nargin<1
    % presume PUBLISH
    subfolder='normalHearing'
    subfolder='impairedHearing'
%    subfolder='allParticipants'
end
fileFolderName= ['..' filesep folder filesep subfolder];

% backgroundPlotName='profile_CMA_L';
% backgroundPlotName= ['..\' backgroundPlotName];
backgroundPlotName='';

rowNo=90; figure(rowNo), clf

dbstop if error
restorePath=path;
addpath (['..' filesep 'MAP'],    ['..' filesep 'multiThreshold 1.46'])
addpath    (['..' filesep 'utilities'])

% input the directory
fileList=dir(fileFolderName);
if isempty(fileList)
    error(['scanAllFiles: No files found in ' fileFolderName])
end

allForegroundResults=[];
allMeanSummaries=[];
filesProcessed=0;
for fileNo=1:length(fileList)
    if fileList(fileNo).isdir
        continue
    end
    filesProcessed=filesProcessed+1;

    foreGroundName=fileList(fileNo).name;
    foreGroundName=foreGroundName(1:end-2); % remove '.m'
    %%
    [foregroundResults meanSummary]= ...
        plotProfile(subfolder, foreGroundName, backgroundPlotName, ...
        rowNo);

    % add fileNumber to first column of foregroundResults 
    [r c]=size(foregroundResults);
    foregroundResults=[repmat(fileNo,r,1) foregroundResults ];
    
    % combine results across files
    allForegroundResults=[allForegroundResults; foregroundResults];
    allMeanSummaries=[allMeanSummaries; meanSummary];
    
end


%% display results
save (folder, 'allForegroundResults')

% allForegroundResults matrix has the following columns
%        1       2     3    4   5       6       7       8
% participantNo BFs short long TMCfr  TMCslope IFMfr IFMCdepth

maxDepth= max(allForegroundResults(:,8));
maxSlope= max(allForegroundResults(:,6))/100;
meanAbsThresholds=[]; meanDepths=[]; meanSlopes=[];
stdevAbsThresholds=[]; stdevSlopes=[]; stdevDepths=[];
rowNo=0;
figure(5),clf
figure(6), clf
for BF=[250 500 1000 2000 4000 6000];
    rowNo=rowNo+1;
    %     figure(rowNo), clf, set(gcf,'name', [int2str(BF) ' Hz'])
    figure(5)
    % plot abs thresholds vs slope
    idx=allForegroundResults(:,2)==BF;
    selectedForegroundResults=allForegroundResults(idx,:);
    subplot(6,3,(rowNo-1)*3+1)
    plot(selectedForegroundResults(:,4),selectedForegroundResults(:,6),'o')
    x=[selectedForegroundResults(:,4) selectedForegroundResults(:,6)];
    [r sampleSize] =UTIL_correlatePairWithNaN (x);
    xlim([0 100]), ylim([-10 100])
    xlabel(['abs threshold(r= ' ...
        num2str(r,'%4.2f') ' N=' int2str(sampleSize) ')'])
    ylabel('TMC slope')
    grid on

    % plot abs threshold vs IFMC depth
    subplot(6,3,(rowNo-1)*3+ 2)
    plot(selectedForegroundResults(:,4),selectedForegroundResults(:,8),'o')
    x=[selectedForegroundResults(:,4) selectedForegroundResults(:,8)];
    [r sampleSize] =UTIL_correlatePairWithNaN (x);
    xlim([0 100])
    xlabel(['abs threshold(r= ' ...
        num2str(r,'%4.2f') ' N=' int2str(sampleSize) ')'])
    ylim([-10 100]), ylabel('IFMC depth')
    grid on
    if BF==250
        title(subfolder)
    end

    % plot slope vs depth
    subplot(6,3,(rowNo-1)*3+3)
    plot(selectedForegroundResults(:,6),selectedForegroundResults(:,8),'o')
    x=[selectedForegroundResults(:,6) selectedForegroundResults(:,8)];
    xlim([-10 100])
    [r sampleSize] =UTIL_correlatePairWithNaN (x);
    if sampleSize>10
        xlabel(['TMC slope (r= ' ...
            num2str(r,'%4.2f') ' N=' int2str(sampleSize) ')'])
    else
        xlabel('TMC slope ')
    end
    ylim([-10 100]),ylabel('IFMC depth')
    grid on

    % Histograms Histograms Histograms Histograms 
    figure (6)
    edges=-10:10:90;
    % hist abs thresholds
    subplot(6,3,(rowNo-1)*3+1)
    N=histc(selectedForegroundResults(:,4),edges);
    bar(edges,N)
    idx=isnan(selectedForegroundResults(:,4));
    N=sum(~idx);
    if BF==6000, xlabel('abs threshold '), end
    if BF==6000, set(gca,'xtick', 0:20:100), else set(gca,'xtick', []),end
    xlim([-10 100])
    y=ylim;
    x=xlim;
    ave=mean(selectedForegroundResults(~idx,4));
    stdev=std(selectedForegroundResults(~idx,4));
    sampleSize=length(selectedForegroundResults(~idx,4));
    if BF==250, text(-70, 1.2*y(2), 'BF (Hz)'), end
    text(0.2*x(2),1*y(2), ...
        [num2str(ave,'%4.0f') ' sd' num2str(stdev,'%3.0f') ' n='...
        num2str(sampleSize,'%3.0f')],...
        'backgroundcolor','w')
    ylabel([num2str(BF) '     '],'rotation',0 )
    meanAbsThresholds=[meanAbsThresholds ave];
    stdevAbsThresholds=[stdevAbsThresholds stdev];

    % hist TMC slope
    subplot(6,3,(rowNo-1)*3+2)
    N=histc(selectedForegroundResults(:,6),edges);
    bar(edges,N,'histc')
    idx=isnan(selectedForegroundResults(:,6));
%     N=sum(~idx);
    xlim([-10 100])
    y=ylim;

    ave=mean(selectedForegroundResults(~idx,6));
    stdev=std(selectedForegroundResults(~idx,6));
    sampleSize=length(selectedForegroundResults(~idx,6));
    text(0.2*x(2),1*y(2), ...
        [num2str(ave,'%4.0f') ' sd' num2str(stdev,'%3.0f') ' n='...
        num2str(sampleSize,'%3.0f')],...
        'backgroundcolor','w')
    meanSlopes=[meanSlopes ave];
    stdevSlopes=[stdevSlopes stdev];

    if BF==6000
        set(gca,'xtick', 0:20:100)
        xlabel('TMC slopes')
    else
        set(gca,'xtick', [])
    end
    if BF==250
        title(subfolder)
    end

    % hist IFMC depth
    subplot(6,3,(rowNo-1)*3+3)
    N=histc(selectedForegroundResults(:,8),edges);
    bar(edges,N,'histc')
    if BF==6000, set(gca,'xtick', 0:20:100), else set(gca,'xtick', []),end
    idx=isnan(selectedForegroundResults(:,8));
    N=sum(~idx);
    if BF==6000, xlabel('IFMC depth'), end
    xlim([-10 100])
    y=ylim;
    x=xlim;

    ave=mean(selectedForegroundResults(~idx,8));
    stdev=std(selectedForegroundResults(~idx,8));
    sampleSize=length(selectedForegroundResults(~idx,8));
    text(0.2*x(2),1*y(2), ...
        [num2str(ave,'%4.0f') ' sd=' num2str(stdev,'%3.0f') ' n='...
        num2str(sampleSize,'%3.0f')],...
        'backgroundcolor','w')
    meanDepths=[meanDepths ave];
    stdevDepths=[stdevDepths stdev];

end

% means cols 1-3, stdevs cols 4-6
meanSDsummary=[meanAbsThresholds;  meanSlopes; meanDepths; ...
    stdevAbsThresholds; stdevSlopes; stdevDepths]';
path(restorePath)
