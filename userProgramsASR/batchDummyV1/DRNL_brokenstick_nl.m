function [y] = DRNL_brokenstick_nl (x, a, b, c)
% y = sign(x).* min(a*abs_x,  b*abs_x .^ c);
% This function could be replaced by a lookup table

% linear (low amplitude) response
y=a.*x;

% compressed high amplitude
compressionThreshold=10.^((1/(1-c)).*log10(b./a));
% only values outside the compression threshold
%  need be subject to compression
abs_x = abs(x);
idx=find(abs_x>compressionThreshold);
if ~isempty(idx)>0
    y(idx) = sign(x(idx)).* ( b*abs_x(idx) .^ c);
end
