close all;clear all;clc
clear classes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not edit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

status = fclose('all');
expName = ['v7_demo' num2str(rand)];

x = jobject('L', fullfile('_testExp','featL'));

x.MAPopHSR = 0;
x.MAPopMSR = 0;
x.MAPopLSR = 1;


% Add some noise to let hearing aid adjust like the efferent would
x.noisePreDur = 0.000;
x.truncateDur = x.noisePreDur;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HA parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x.channelBFs= [250; 500; 1000; 2000; 4000;]; %MUST BE A COLUMN!!!!
        
x.TAspl = 50 * ones(5, 1); %Thresholds actual (in dB SPL)
x.TDspl = 5 * ones(5, 1); %Thresholds desired (in dB SPL)
x.TDspl = [10; 10; 5; 5; 15];

x.TCspl = 35 * ones(5, 1);      %Compression thresholds (in dB SPL)
x.TMspl = 30 * ones(5, 1);     %MOC thresholds (in dB SPL)

x.ARtau  = 0.03;            % decay time constant
x.ARthresholddB = 80;      % dB SPL (input signal level) =>200 to disable

x.MOCtau = 0.3;  %0.06
x.MOCfactor = 1e6;       % now that the conversions between velocity and displacement have been removed, internal values are 100 times greater (10kHz / 100 Hz) and so defaults need to be at least 100 times less than the old value of 2e9

x.DRNLc = 0.2 * ones(5, 1);
x.bwOct =    1 * ones(5, 1); %Octaves

x.filterOrder  = 4; %This sounds better than 2nd order

x.useAid = 0;


% for Lev = -10:10:100
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MAP parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    normal=1;
    
    x.MAPuseEfferent = 1;
    x.MAPDRNLSave   = 1;
    
    x.currentSpeechLevel = 65;
    x.currentNoiseLevel  = -200;
    
    if normal
        x.participant = 'Normal'; 
    else
        x.participant = 'OHC'; 
    end
    
    x.wavFolder = fullfile(pwd, 'demo_wavs', 'microAURORA');
    x.noiseFolder = fullfile(pwd, 'demo_wavs', 'noises');
%     currentWav = 'MHS_2841A.wav';
    currentWav = 'MFF_3B';
    
    x.noiseName = 'pink';%'pink';%
    x.freezeNoise = true;
    
    
   
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Do the time consuming bit
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic
    [~, ANresponse, method] =  x.processWavs(currentWav);
    toc
    
    % Winsize and hopsize control the amount of auditory spectrogram smoothing
    % and spectrogram smoothing
    winSize = 25; %25
    hopSize = 10; %10
    
    % Plain old spectrogram
    x.noisePreDur = 0;
    x.truncateDur = 0;
    [stimulus, fs] = x.getStimulus(currentWav);
    %disp(20*log10(sqrt(mean(stimulus.^2))/20e-6))
    
    winSizeSamples = round(winSize*fs/1000);
    hopSizeSamples = round(hopSize*fs/1000);
    window = hanning(winSizeSamples);
    
    F = method.nonlinCF;
    [S,F,T,P] = spectrogram(stimulus,window,hopSizeSamples,F,fs);
    
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Make figures
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    global ANdt
    ANrate = ANresponse;%/ANdt; % convert rate into spikes/sec
    
    smoothSpect = x.makeANsmooth(ANrate, 1/ANdt, winSize, hopSize);
    
    if normal
        hFig(1) = figure(88);
    else
        hFig(1) = figure(89);
    end
    
    labSize = 14;
    hAxis(1) = gca;
    tAx = (hopSize:hopSize:hopSize*size(smoothSpect,2) ); % time axis
    hSurf(1) = surf(tAx, 1:size(smoothSpect,1), smoothSpect); shading interp
    
    hFig(2) = figure(66);
    hAxis(2) = gca;
    hSurf(2) = surf(1000*T, 1:size(P,1), 94 + 10*log10(P) ); shading interp
    
    maxRate = 300;
    zlim(hAxis(1), [0 maxRate])
    set(hAxis(1), 'CLim', [0 maxRate])
    
    maxLevel = 80;
    zlim(hAxis(2), [0 maxLevel])
    set(hAxis(2), 'CLim', [0 maxLevel])
    
    set(hFig(1), 'Name', 'Auditory spectrogram');
    set(hFig(2), 'Name', 'Regular spectrogram');
    zlabel(hAxis(1), 'Spikes/s')
    zlabel(hAxis(2), 'Level in dB SPL')
    
    title(hAxis(1),['speech = ' num2str(x.currentSpeechLevel) ...
        ', noise = ' num2str(x.currentNoiseLevel) ...
        ', aid = ' num2str(x.useAid) ...
        ', efferent = ' num2str(x.MAPuseEfferent)], ...
        'Color', [1 1 1], 'FontSize', labSize );
    
    title(hAxis(2),['speech = ' num2str(x.currentSpeechLevel) ...
        ', noise = ' num2str(x.currentNoiseLevel)], ...
        'Color', [1 1 1], 'FontSize', labSize );
    
    ylabel(hAxis(1), 'CF in kHz')
    ylabel(hAxis(2), 'Freq in kHz')
    
    % Here are modifications that we want to apply to both figures
    for nn = [1 2]
        view (hAxis(nn), [26 58])
        
        xlabel(hAxis(nn), 'Time in ms')
        ylim(hAxis(nn),[1 size(smoothSpect,1)])
        
        %The following lines put meaningful numbers on the frequency axis
        numYTicks = 6;
        tickLocs = round(linspace( 1, numel(method.nonlinCF), numYTicks));
        set(hAxis(nn), 'YTick', tickLocs);
        set(hAxis(nn), 'YTIckLabel', num2str(  (  method.nonlinCF(tickLocs)/1000  )' , '%0.2f'  ));
        
        set(hAxis(nn), 'XColor', [1 1 1], 'YColor', [1 1 1], 'ZColor', [1 1 1], 'Color', [0 0 0])
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
        
        drawnow
    end
    
%     %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Do some rate-level type stuff
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     chan1k = find(method.nonlinCF>=1000,1);
%     rateNow = (mean(smoothSpect(  chan1k , 1:end-5  )));
%     dbrmsDRNLNow = 20*log10(sqrt((mean(method.DRNLData(  chan1k , 2000:end-2000  ).^2))));
%     if exist('rateLev')
%         rateLev = [rateLev; rateNow]
%     else
%         rateLev =  rateNow
%     end
%     
%     if exist('dispLev')
%         dispLev = [dispLev; dbrmsDRNLNow]
%     else
%         dispLev =  dbrmsDRNLNow
%     end
%     
%     if exist('Levs')
%         Levs = [Levs; Lev];
%     else
%         Levs =  Lev;
%     end
% end
% 
% %%
% figure(11);
% subplot(2,1,1); plot(Levs,rateLev); xlabel('Input level in dB SPL'); ylabel('Spikes per second')
% subplot(2,1,2); plot(Levs,dispLev); xlabel('Input level in dB SPL'); ylabel('RMS displacement in dB')
% 

%disp(20*log10(sqrt(mean(stimulus.^2))/20e-6))
% showMAP
