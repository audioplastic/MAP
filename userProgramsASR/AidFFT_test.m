clear variables; clear classes; clc
% 
% x = EssexAid_WrapClass;
% 
% % [sig,srIN] = wavread(fullfile('demo_wavs','microAURORA','MHS_2841A'));
% [sig,srIN] = wavread('x1');
% 
% sig = sig(2000:numel(sig)/2 + 1000);
% 
% [p,q] = rat(x.sr/srIN,0.0001);
% sig = resample(sig, p, q);
% 
% 
% % % soundsc(sig,x.sr)
% 
% x.stimulusUSER = sig;
% x.TM_dBHL = ones(size(x.TM_dBHL))*1000;
% x.ARthreshold_dB = 1000;
% x.TC_dBHL = ones(size(x.TM_dBHL))*20;
% x.bwOct = 1/5;
% x = x.processStim;
% 
% soundsc(x.aidOPnice,x.sr)
% 





x = EssexAid_WrapClassFFT;

% [sig,srIN] = wavread(fullfile('demo_wavs','microAURORA','MHS_2841A'));
% sig = [1; zeros(2^12-1,1)]; srIN = 48e3;
% srIN = 48e3; sig = sin(2*pi*1000*(0:1/srIN:0.4))';
[sig,srIN] = wavread('x1');

[p,q] = rat(x.sr/srIN,0.0001);
sig = resample(sig, p, q);


x.stimulusUSER = sig;
x.TM_dBHL = ones(size(x.TM_dBHL))*1000;
x.TC_dBHL = ones(size(x.TM_dBHL))*20;
x.ARthreshold_dB = 1000;
x.bwOct = 1/5;

x.nfft = 2048;
x.numSamples = x.nfft/2;

x = x.processStim;
soundsc(x)
plot(x.aidOPnice)

%
xNorm =  x.aidOPnice./sqrt(mean(x.aidOPnice.^2));
xNorm = xNorm/10;
wavwrite(xNorm, x.sr, 'S26ChanFFT_20')


%%


% %% C Implementablbe OLA rectangle (fast - NO Overlap)
% Nsig = numel(sig);
% nfft = x.nfft;
% M = x.numSamples;
% nChans = x.numChannels;
% 
% yI = zeros(Nsig,1); % allocate output
% FILO = zeros(nfft,1);
% 
% Nframes = 1+floor((Nsig-M)/M);  % no. complete frame
% 
% load('imps')
% 
% impsT = imps(:,1:nfft/2);
% H = fft([impsT zeros(size(impsT))]');
% 
% 
% tic
% for m = 0:(Nframes-1)
%     index = m*M+1:min(m*M+M,Nsig); % indices for the mth frame
%     xm = sig(index);  % windowed mth frame (rectangular window)
%     xmzp = [xm; zeros(nfft-length(xm),1)]; % zero pad the signal
%     Xm = fft(xmzp);    
%     Ym = repmat(Xm,1,nChans) .* H;   % freq domain multiplication
%     ym = ifft(Ym);         % inverse transform    
%     FILO = [FILO(M+1:end); zeros(M,1)] + sum(ym,2); %first in, last out buffer
%     yI(index) = FILO(1:M); % overlap add
% end
% 
% 
% clf; plot(yI,'k');hold on; plot(x.aidOPnice,'r');


