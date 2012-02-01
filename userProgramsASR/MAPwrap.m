% This function wraps up whatever version of MAP I want to call. It is
% implemented partly because I want to avoid messing with jobject too much
% and partly because I dont want to declare globals in my class.
function [myANprobRateOutput, mydt, myBF] = MAPwrap(stimulus, sampleRate, BFlist, participant, AN_spikesOrProbability, paramChanges)


global ANprobRateOutput ANoutput dt dtSpikes savedBFlist
% disp(20*log10(sqrt(mean(stimulus.^2))/20e-6))
MAP1_14(stimulus, sampleRate, BFlist, participant, AN_spikesOrProbability, paramChanges);
% disp(20*log10(sqrt(mean(stimulus.^2))/20e-6))

if strcmpi(AN_spikesOrProbability, 'spikes')
    myANprobRateOutput   = ANoutput; 
    mydt = dtSpikes; 
else
    myANprobRateOutput   = ANprobRateOutput;
    mydt = dt;
end


myBF = savedBFlist;
