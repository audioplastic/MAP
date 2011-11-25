close all; clear all; clc

sr = 44100;
dt = 1/sr;
tc = 150e-3; %seconds
fc = 1/(2*pi*tc);
[b,a] = butter(1, fc/(sr/2), 'low');

% The time constant is the time required to charge the capacitor, through
% the resistor, to 63.2 (? 63) percent of full charge; or to discharge it
% to 36.8 (? 37) percent of its initial voltage. These values are derived
% from the mathematical constant e, specifically 1?exp(-1) and exp(-1)
% respectively. It is completely symmetric when applied to a linear system.
for nn = 1:3
    x = [zeros(sr,1); ones(2*sr,1); zeros(sr,1);]*nn;
    y = filter(b,a,x);
    figure(1); hold on; plot(y); title('Limitless RC integrator')     
    %set(gca,'YScale', 'log')
end

% However, if we introduce a cap, then the rise time depends on the input
% level, while the decay is a pure function of the time constant.
%
% BUT THERE SHOULD BE A WAY OF NOT INTRODUCING A CAP TO GET CONSTANT OUTPUT
% LIKE THE BUCKET EXAMPLE
for nn = 1:3
    limit = 1;
    x = [zeros(sr,1); ones(2*sr,1); zeros(sr,1);]*nn;
    y(1) = 0;
    for t = 2:numel(x)
        y(t) = b(1)*x(t) + b(2)*x(t-1) - a(2)*y(t-1);%diff eqn for matlab filter
        y(t) = min (y(t), limit);
    end
    figure(2); hold on; plot(y); title('Limited RC integrator')     
end


% In mathematics, a leaky integrator equation is a specific differential
% equation, used to describe a component or system that takes the integral
% of an input, but gradually leaks a small amount of input over time. It
% appears commonly in hydraulics, electronics, and neuroscience in
% single-neuron models
%
% The general solution is x(t)=k*exp(-A*t)+C/A
% Where C is the input, A is the rate of leak, and k is a constant
% THE RESULT IS IDENTICAL TO THE LIMITLESS INTEGRATOR

A=0.0002;
x = [zeros(sr,1); ones(2*sr,1); zeros(sr,1);];
for nn = 1:3
    y=0;
    yHistory = zeros(size(x));
    for t = 2:numel(x)
        y = y + nn*x(t);  %calculate the volume in bucket + new fluid volume
        y = y - A*y;      %calculate the leaked amount
        yHistory(t) = y;  %just for plotting
    end
    figure(3); hold on; plot(yHistory); title('"Leaky" integrator')     
end



