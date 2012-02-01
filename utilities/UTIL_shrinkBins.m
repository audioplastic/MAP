function [reducedBins adjustedBinWidth]=UTIL_shrinkBins(inputData, dt, binWidth)
% UTIL_shrinkBins accumulates values across adjacent matrix columns.
%  and returns the *average* values per bin
% The bin width is adjusted to make it an integer multiple of dt
% usage:
%	[reducedBins adjustedBinWidth]=UTIL_shrinkBins(inputData, method)
%
% arguments
%	inputData is a channel x time matrix
%	reducedBins is the reduced matrix (channel x time/binWidth)
%

% the binwidth must be an integer multiple of dt and >0
adjustedBinWidth=dt*ceil(binWidth/dt);

[numChannels numDataPoints]= size(inputData);

% Multiple fibers is the same as repeat trials
% Consolidate data into a histogram 
dataPointsPerBin=round(adjustedBinWidth/dt);

if dataPointsPerBin<=1;
% 	Too few data points
	reducedBins=inputData;
	return
end

numBins=floor(numDataPoints/dataPointsPerBin);
reducedBins=zeros(numChannels,numBins);

% take care that signal length is an integer multiple of bin size
%  by dropping the last unuseable values
useableDataLength=numBins*dataPointsPerBin;
inputData=inputData(:,1:useableDataLength);

for ch=1:numChannels
	% Convert each channel into a matrix where each column represents 
	% the content of a single reducedBins bin
	reducedBins2D=reshape (inputData(ch,:), dataPointsPerBin, numBins );
	% and sum within each bin (across the rows
	reducedBins(ch,:)=mean (reducedBins2D,1);
end
