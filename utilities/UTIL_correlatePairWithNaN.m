function [r sampleSize] =UTIL_correlatePairWithNaN (x)
% x is a two column matrix
% r is the correlation between pairs with no NaNs
% sampleSize is the number of pairs that were found
sampleSize=0;
r=NaN;

    idx=find(~isnan(x(:,1)));     x=x(idx,:);
    idx=find(~isnan(x(:,2)));     x=x(idx,:);
    
    [sampleSize c]=size(x);
    if sampleSize>2
    r=corrcoef(x); 
    r=r(2,1);
    end
    