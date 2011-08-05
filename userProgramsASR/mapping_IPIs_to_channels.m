function [y,ctr_freq] = mapping_IPIs_to_channels(ipihist,sfreq,BFlist,numChannels,useMax)%(x,factor,useMax)

%function that maps the IPI histogram output of an IPI analysis to
%frequency channels with a certain spacing defined by BFlist
%first, channel borders are calculated and mapped to the IPI values
%returning indices
%second, the mean or max between these indices is taken as a value for the
%ouput ipi histogram with fewer channels
%
% Tim Jürgens, March 2011, extension of code from Matthew Robertson
%
% input: ipihist: IPIhistogram with dimensions IPIs (1) and time step (2)
%        sfreq: sampling frequency
%        BFlist: list of best frequencies
%        numChannels: number of channels desired for the output
%        usemax: 1 -> use maximum in each IPI bin
%                0 -> use mean in each IPI bin
% output: y: adjusted IPI pattern with dimensions channels (1), which are a
%            total of numChannels and time step (2)
%            ctr_freq: center_frequencies of the channels

ipihistlength = size(ipihist,1);
channelstep = length(BFlist)/numChannels;

%preallocation due to speed
channelborders = zeros(1,numChannels+1);
channelborders_IPI_index = zeros(1,numChannels+1);
y = zeros(numChannels,size(ipihist,2));
ctr_freq = zeros(1,numChannels);

%the first values (lower border)
channelborders(1) = BFlist(1);
[tmp, channelborders_IPI_index(1)] = min(abs(channelborders(1)-1./([1:ipihistlength]./sfreq)));

%loop for every channel
for iCounter = 1:numChannels
    channelborders(iCounter+1) = BFlist(round(iCounter*channelstep));
    ctr_freq(iCounter) = mean([channelborders(iCounter) channelborders(iCounter+1)]);
    [tmp,channelborders_IPI_index(iCounter+1)] = min(abs(channelborders(iCounter+1)-1./([1:ipihistlength]./sfreq)));
    if ~useMax
        y(iCounter,:) = mean(ipihist(channelborders_IPI_index(iCounter+1):channelborders_IPI_index(iCounter),:));
    else
        y(iCounter,:) = max(ipihist(channelborders_IPI_index(iCounter+1):channelborders_IPI_index(iCounter),:));
    end
end





% [chans,frames] = size(x);
% 
% numfeatures = chans/factor;
% y = zeros(numfeatures,frames);
% 
% for i = 1:numfeatures
%     st=1+(i-1)*factor;
%     fn=st+factor-1;
%     %fprintf('%d %d\n',st,fn);
%     if ~useMax
%         y(i,:)=mean(x(st:fn,:)); %previously max
%     else 
%         y(i,:)=max(x(st:fn,:)); %previously max
%     end;
%     
% end