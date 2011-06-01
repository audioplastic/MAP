% smooth
hann = GJB_hanning(winSizeSamples);

rate=[];
for chan=1:numChannels,
    f=enframe(earObject(chan,:),hann,hopSizeSamples);
    rate(chan,:)=mean(f,2)';
end

% make feature vectors
features = GJB_dct(rate);

features = features(1:14,:);
frameshift = 10;                    % shift between frames (ms)
sampPeriodFromMsFactor = 10000;     % appropriate for 10ms frame rate
paramKind = 9;                      % HTK USER format for user-defined features
byteOrder = 'be';                   % byte order is big endian
writeHTK(targetpath,features,size(features,2),frameshift*sampPeriodFromMsFactor,size(features,1)*4,paramKind,byteOrder);