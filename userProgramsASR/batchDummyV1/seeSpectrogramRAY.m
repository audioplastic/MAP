function seeSpectrogramRAY(normal, useAid)
% close all;clear all;clc
% clear classes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not edit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dbstop if error

if nargin<2
    normal=0;
    useAid=1;
end

% status = fclose('all');
% expName = ['v7_demo' num2str(rand)];

x = jobject('L', fullfile('_testExp','featL'));

x.MAPopHSR = 1;
x.MAPopMSR = 0;
x.MAPopLSR = 0;


% Add some noise to let hearing aid adjust like the efferent would
x.noisePreDur = 0.500;
x.truncateDur = 0.500;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HA parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x.useAid = useAid;

x.opScaling_dB = 105;   % 106.58 --> See experiment 8 in log book - 152?

x.ARthresholddB=70;         % 60
x.ARtau=0.2;               % 0.03

x.DRNLaBaseline=1e3;        % 1e4
x.DRNLb=0.5e-5 ;            % 0.5e-5 ?? --> 
% YES but misleading. This = 5e-6 as in the readme % 
% (watch out fot the 0.5 vs 5)
x.DRNLc=0.2;

x.MOCfactor=0;            % 2e9
x.MOCtau=0.03;              % 0.06;
x.MOCthreshold = 5e-10 * 40^(0/20);

x.channelBFs= [400 800 1600 3200];%[250 630 1587 4000];

x.bwOct                       =  [1 1 1 1];%[4/3 4/3 4/3 4/3]; %Octaves
x.DRNLaBaselineRelative_dB    = 45* [1 1 1 1];
x.DRNLbRelative_dB            = 45* [1 1 1 1];
x.filterOrder                 =  2; %now changed %
% --> See experiment 6 in log book

x.doHpFilt                    = 1;    %TRUE


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAP parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% normal=0;
if normal
    x.useAid=0;
end

x.MAPuseEfferent = 0;

x.currentSpeechLevel = 60;
x.currentNoiseLevel  = 50;

if normal
    x.participant = 'DEMO2_multiSpont'; 
    x.participant = 'MPa';
else
    % choose impaired participant
    x.participant = 'DEMO2_multiSpontIMP'; 
    x.participant = 'KFr'; 
    x.participant = 'JJo';
    x.participant = 'JSan';
end

x.wavFolder = fullfile(pwd, 'demo_wavs', 'microAURORA');
x.noiseFolder = fullfile(pwd, 'demo_wavs', 'noises');
currentWav = 'MHS_2841A.wav';

x.noiseName = '8TalkerBabble';


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do the time consuming bit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
[~, ANresponse, method] =  x.processWavs(currentWav);
toc

% Winsize and hopsize control the amount of auditory spectrogram smoothing
% and spectrogram smoothing
winSize = 25;
hopSize = 10;

% Plain old spectrogram
x.noisePreDur = 0;
x.truncateDur = 0;
[stimulus, fs] = x.getStimulus(currentWav);

winSizeSamples = round(winSize*fs/1000);
hopSizeSamples = round(hopSize*fs/1000);
window = hanning(winSizeSamples);

F = method.nonlinCF;
[S,F,T,powerSpectrum] = spectrogram(stimulus,window,hopSizeSamples,F,fs);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
labSize = 14; % font size
ANrate = ANresponse/method.AN_IHCsynapsedt; % convert rate into spikes/sec

smoothSpect = ...
    x.makeANsmooth(ANrate, 1/method.AN_IHCsynapsedt, winSize, hopSize);

if normal
        specialLabel=[x.participant ];
    hFig(1) = figure(88);
    nSS=smoothSpect;
    save 'normalAuditorySpectrum' nSS
else
    if ~x.useAid
        specialLabel=[x.participant ' - unaided '];
        % unaided impaired
        hFig(1) = figure(89);
    else
        load 'normalAuditorySpectrum'
            hFig(1) = figure(90);
        specialLabel=[x.participant ' - with aid '];
        if size(smoothSpect)==size(nSS)
            % aided impaired
            iSS=smoothSpect;         
            save 'aidedAuditorySpectrum'iSS
            diff=iSS-nSS;
            hFig(3)=figure(1);
            hAxis(3) = gca;
            timeAxis = (hopSize:hopSize:hopSize*size(smoothSpect,2) ); 
            surf(timeAxis, 1:size(smoothSpect,1), diff); shading interp
            %         figure(1), surf(diff),shading interp
            nAs=mean(mean(diff.^2))^0.5
            title(['impaired -normal:   rms= ' num2str(nAs) ...
                ' sum= ' num2str(mean(mean(iSS)))],...
                'Color', [1 1 1], 'FontSize', labSize)
            zlabel('spikes/s)'), ylabel('channel'), xlabel('time in ms')
        end
    end
end

% index of hAxis: 1--AN spikes, 2-powerSpecrum, 3-diff matrix
% plot AN response
figure(hFig(1))
hAxis(1) = gca;
timeAxis = (hopSize:hopSize:hopSize*size(smoothSpect,2) ); % time axis
hSurf(1) = surf(timeAxis, 1:size(smoothSpect,1), smoothSpect); 
shading interp

maxRate = 300;
zlim(hAxis(1), [0 maxRate])
set(hAxis(1), 'CLim', [0 maxRate])
ylabel(hAxis(1), 'CF in kHz')

set(hFig(1), 'Name', 'Auditory spectrogram');
zlabel(hAxis(1), 'Spikes/s')

title(hAxis(1),specialLabel,'Color', [1 1 1], 'FontSize', labSize );;
% title(hAxis(1),[specialLabel ' speech = ' num2str(x.currentSpeechLevel) ...
%         ', noise = ' num2str(x.currentNoiseLevel) ...
%         ', aid = ' num2str(x.useAid) ...
%         ', efferent = ' num2str(x.MAPuseEfferent)], ...
%         'Color', [1 1 1], 'FontSize', labSize );

% plot power spectreum
hFig(2) = figure(66);
hAxis(2) = gca;
hSurf(2) = surf(1000*T, 1:size(smoothSpect,1), ...
    94 + 10*log10(powerSpectrum) ); shading interp

maxLevel = 80;
zlim(hAxis(2), [0 maxLevel])
set(hAxis(2), 'CLim', [0 maxLevel])

set(hFig(2), 'Name', 'Regular spectrogram');
zlabel(hAxis(2), 'Level in dB SPL')

title(hAxis(2),['Spectrogram:  speech = ' num2str(x.currentSpeechLevel) ...
    ', noise = ' num2str(x.currentNoiseLevel)], ...
    'Color', [1 1 1], 'FontSize', labSize );

ylabel(hAxis(2), 'Freq in kHz')

% Here are modifications that we want to apply to all figures
for nn = 1:length(hAxis)
    view (hAxis(nn), [26 58])
    
    xlabel(hAxis(nn), 'Time in ms')
    ylim(hAxis(nn),[1 size(smoothSpect,1)])
    
    %The following lines put meaningful numbers on the frequency axis
    numYTicks = 6;
    tickLocs = round(linspace( 1, size(smoothSpect,1), numYTicks));
    tickLocs = (linspace( 1, size(smoothSpect,1), numYTicks));
    set(hAxis(nn), 'YTick', tickLocs);
    set(hAxis(nn), 'YTIckLabel', ...
        num2str((method.nonlinCF(round(tickLocs))/1000  )' , '%0.2f'  ));
    
    set(hAxis(nn), 'XColor', [1 1 1], 'YColor', [1 1 1], ...
        'ZColor', [1 1 1], 'Color', [0 0 0])
    set(hFig(nn), 'Color', [0 0 0])
    
    set(get(hAxis(nn), 'Xlabel'), 'FontSize', labSize)
    set(get(hAxis(nn), 'Ylabel'), 'FontSize', labSize)
    set(get(hAxis(nn), 'Zlabel'), 'FontSize', labSize)
    
    %----------------------------------------------------------------
    % Arty lighting effects (comment this section out if it runs slow)
    % This is probably a bit excessive and can be safely commented out
    % from here down.
    % -----------------------------------------------------------------
    %     colormap(hAxis(nn), jet)
    %     set(hSurf(nn), ...
    %         'EdgeColor','none', ...
    %         'FaceLighting','phong', ...
    %         'AmbientStrength',0.9, ...
    %         'DiffuseStrength',0.4, ...
    %         'Clipping','off',...
    %         'BackFaceLighting','lit', ...
    %         'SpecularStrength',1.1, ...
    %         'SpecularColorReflectance',1, ...
    %         'SpecularExponent',7);
    %
    %     l1 = light('Position',[40 100 20], ...
    %         'Style','local', ...
    %         'Color',[0.8 0.8 0.6], ...
    %         'parent',hAxis(nn));
    %
    %     l2 = light('Position',[.5 -1 .4], ...
    %         'Color',[0.6 0.8 0.8], ...
    %         'parent',hAxis(nn));
end


