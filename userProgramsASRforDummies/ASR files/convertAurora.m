addpath(pwd)
cd  D'igit Training Data'\Clean\

d = dir('*.08');
mkdir('wavs')
for I = 1:length(d)
    x = readAurora(d(I).name);
    x = 0.99*(x./max(abs(x)));
    fileName = strrep(d(I).name,'.08','.wav')
    wavwrite(x,44100,['.' filesep 'wavs' filesep fileName]);
end;

cd ..
cd ..

cd Di'git Triplets Test Data'\;


d = dir('*.08');
mkdir('wavs')
for I = 1:length(d)
    x = readAurora(d(I).name);
    x = 0.99*(x./max(abs(x)));
    fileName = strrep(d(I).name,'.08','.wav')
    wavwrite(x,44100,['.' filesep 'wavs' filesep fileName]);
end;

cd ..
rmpath(pwd)
