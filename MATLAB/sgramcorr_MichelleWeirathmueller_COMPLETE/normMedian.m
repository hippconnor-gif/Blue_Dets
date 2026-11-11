function sp1 = normMedian(spect, fRate, medwidth, pctile)
%normMedian	Apply a median filter to a spectrogram.
%
% spect = normMedian(spect, binBW, fRate, medwidth [,pctile])
%   Apply a median filter horizontally across each frequency bin of a 
%   spectrogram.  medwidth specifies the distance, in seconds, over which
%   the median is taken.  pctile specifies what percentile value to
%   use; it defaults to 0.50, which is the median.

if (nargin < 4),    pctile = NaN; end
if (isnan(pctile)), pctile = 0.50; end

w = round(fRate * medwidth);
sp1 = morphomedian(spect, [(1 : max(w,1))', zeros(w, 1)], pctile);

% Make sp1 the same shape as spect by duplicating the last column.
sp1 = [sp1 (spect(:, nCols(spect)) * ones(1, nCols(spect) - nCols(sp1)))];
