function spect = normPct(spect, binBW, perFreq, fMin, fMax, pct0, pct1, mul0, mul1)
% normPct	Percentile normalization (aka low order statistics).
%
% spect = normPct(spect, binBW, perFreq, fMin, fMax, pct0, pct1, mul0, mul1)
%   Any of the pct or mul arguments may be NaN, in which case it is ignored.
%   If perFreq is 0, fMin and fMax are ignored.
% 
% See the 'p' option to normGram.

if (nargin < 6), pct0 = NaN; end
if (nargin < 7), pct1 = NaN; end
if (nargin < 8), mul0 = NaN; end
if (nargin < 9), mul1 = NaN; end

pcts = [iff(isnan(pct0), 0, pct0) iff(isnan(pct1), 0, pct1)];

% Compute v = normalization levels.
if (perFreq)
  v = percentile(spect.', pcts);
else
  % Non-per-frequency: Do it on the sub-gram specified by fMin and fMax.
  loBin = max(1,            round(fMin / binBW));
  hiBin = min(nRows(spect), round(fMax / binBW));
  x = spect(loBin:hiBin, :);
  v = percentile(x(:), pcts);
end

% Multipliers, for doing low-order statistics.
if (~isnan(mul0)), v(1,:) = v(1,:) * mul0; end;
if (~isnan(mul1)), v(2,:) = v(2,:) * mul1; end;

if (~isnan(pct0))
  spSub = repmat(v(1,:)', size(spect) ./ size(v(1,:)'));
  %if (nCols(v) == 1), printf('normPct: subtracting %g', v(1,1)); end
  spect = spect - spSub;
else
  spSub = 0;
end
if (~isnan(pct1))
  spDiv = repmat(v(2,:)', size(spect) ./ size(v(2,:)')) - spSub;
  %printf('normPct: dividing by %g (= %g - %g)', v(2,:) - v(1,:), v(2,:),v(1,:))
  spect = spect ./ spDiv;
end

