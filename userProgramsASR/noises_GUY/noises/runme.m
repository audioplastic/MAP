% aurora files are 16 bit ints, big-endian

close all

files = dir('*.raw');

for i=1:length(files),
    
    filename = files(i).name;
    
fs=8000;
ifp = fopen(filename,'r','b'); 
x = fread(ifp,inf,'int16');
fclose(ifp);



win = hanning(1024);

frames = enframe(x,win);

h = abs(fft(frames,[],2));

m = mean(h);

figure; plot(10*log10(m(1:512)));
axis tight; box on; 
title(filename);
xlabel('frequency');
ylabel('dB re 1');

soundsc(x(1:fs*4),fs);
pause
end

