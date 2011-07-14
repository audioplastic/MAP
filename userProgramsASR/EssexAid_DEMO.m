clc; close all; clear classes

% First off, load up the essex aid
x = EssexAid_WrapClass

% Ommiting the semicolon shows a number of properties that are accesible to
% the user of the class.

% The default stimulus for the aid is a sequence of tone pips with
% increasing level that are separated by silent intervals. This is a useful
% stimulus for debugging and inspection of the aid performance. However, it
% is easy to override this default to use any user defined mono or stereo
% stimulus. The wrapper also overrides Matlab's default plot function,
% giving a quick way to inspect what the aid is doing. Plotting the wrapper
% object before running the aid algorithm manually will just show the
% envelope of the input stimulus.

plot(x)

% The wrapper class contains a handy helper function for generating
% different types of tone sequences that can be used as follows.

silDur = 0.1; %seconds
pulseDur = 0.1; %seconds
dBlevs = [40 50 60 40 50 60];
freq = 1000;
sampleRate = x.sr;

x.stimulusUSER = EssexAid_WrapClass.pipSequence(sampleRate, freq, dBlevs, pulseDur, silDur);

% This new stimulus envelope can be viewed by using the plot method once again

plot(x)

% It is also possible to input any stimulus type by setting the properties
% manually. The wrapper expects the input stimulus to be in units of
% pascals. Therefore, to input a 70 dB Gaussian noise with a 25 kHz sample
% rate we would do the following . . 

newSr = 25e3;
nz = randn(newSr/2, 1);
nz = nz./sqrt(mean(nz.^2)); %normalize RMS to 1
nz = nz * 20e-6 * 10^(70/20); %scale to 70 dB SPL

x.stimulusUSER = nz;
x.sr = newSr;

plot(x)

% It is also possible to specify the stimulus and sample rate when
% instantiating a new instance of the class as shown by the following
% example.

y = EssexAid_WrapClass(newSr, nz)
plot(y)

% But let us return to the default test stimulus for the first example of
% aid processing in this demo

x = EssexAid_WrapClass;

% First set up a presciption for a typical 40 dB flat loss

x.audiometry_dB = 40*ones(6,1);

% The audiometry_dB property is the pure tone thresholds of the listener in
% dB SPL at frequencies octave spaced between 0.25 and 8 kHz.
%
% Next we can set the gain of the aid to 50% of the loss.

x.mainGain_dB = x.audiometry_dB * 0.5;

% We shall assume that this imaginary person has no residual compression
% and so we shall set the instantaneous compression threshold to 20 dB
% above the listener's detection threshold.

x.TC_dBHL = 20*ones(6,1);

% To keep things simple to understand, for the first example, we will
% disable the MOC response by setting the MOC threshold to a high level.

x.TM_dBHL = 200*ones(6,1);

% Now process the stimulus using the hearing aid algorithm and view the
% output

x = x.processStim;
plot(x)

% The aid output is shown in red and the original stimulus is shown in
% black. The listener has a 40 dB loss, but with 20 dB gain applied.
% therefore, 0-dB HL = 20-dB SPL

% The first pulse was input at 20 dB spl and comes out at 40 dB SPL
% after the gain was applied. The compression threshold in this
% prescription was set to 20 dB HL. The dummy listener here has pure tone
% thresholds of 40 dB SPL and a 20 dB gain has been applied, therefore, the
% compression threshold
