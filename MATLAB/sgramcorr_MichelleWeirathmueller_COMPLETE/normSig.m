function xHat = normSig(x, alpha, r, p)
%normSig        Exopnential-decay normalization for time series
%
% xHat = normSig(x, alpha, r, p)
%    x       input signal
%    alpha   time constant, close to 1.0
%    r       scaling factor so output signal is around +-r; default 1.0
%    p       use L_p norm for long-term averaging (default p=1; p=2 gives power)

if (nargin < 3 || isnan(r)), r = 1.0; end        % contant multiplier
if (nargin < 4 || isnan(p)), p = 1; end        % for L_p norm; p=2 gives power

xBarP = mean(abs(x(1 : min(100, length(x))))) ^ p;
%xBarP= mean(   (x(1 : min(100, length(x))))) ^ p;  % for making figure
xHat = zeros(1, length(x));
for i = 1 : length(x)
  xHat(i) = r * x(i) / xBarP^(1/p);
  xBarP = alpha * xBarP + (1 - alpha) * abs(x(i))^p;
end
