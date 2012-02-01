function UTIL_cascadePlot(toPlot, colValues, waveHeight)
% UTIL_cascadePlot is cunning code to represent rows as parallel lines
%  in the current axes. No titles or labels are applied here
% input arguments:
%   toPlot is matrix to be displayed
%   colValues is x-axis values (e.g. lags in ACF)

if nargin<1,
    error('UTIL_cascadePlot: no data found')
end
[nRows nCols]=size(toPlot);
if nRows<2
    error('UTIL_cascadePlot: only one row found')
end

if nargin<2, colValues=1:nCols; end

if nargin<3, waveHeight=0.1; end
peakGain=nRows*waveHeight;

% a is the height to be added to each channel
a=max(max(toPlot))*(0:nRows-1)';

% peakGain emphasises the peak height
% lower channels obscure higher channels
x=peakGain*toPlot+repmat(a,1,nCols);
x=nRows*x/max(max(x));

for row=2:nRows
    x(row,:)=max(x(row-1,:), x(row,:));
end

plot(colValues,   x','k')
ylim([0 nRows])

