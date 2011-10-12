clear all;
clc

wName = 'pink'

[x,fs] = wavread([wName '.wav']);

x = [x; x; x];

[b,a] = butter(4, [100 3.999e3]/(fs/2), 'bandpass');

y = filter(b,a,x);

soundsc(y(1:fs),fs)

wavwrite(y,fs,[wName '_bp'])