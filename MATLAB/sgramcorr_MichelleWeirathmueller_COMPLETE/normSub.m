function spect = normSub(spect, binBW, f0, f1, pctile, mult)
% normSub	Subtraction normalization.  See normGram.
%
% spect = normSub(spect, binBW, fRate, prev, f0, f1, mult)
%
% See normGram for an explanation.

b = round([f0 f1] / binBW) + 1;	% frequency band to use as average power

p = exp(spect);			% power spectrogram; is >0 since spect is real
if (b(1) == b(2))
  g = spect;			% handle case where spect(b0:b1,:) is a vector
else
  g = percentile(spect(min(b):max(b),:), pctile);
end

spect = spect - mult * repmat(g, nRows(spect), 1);

% Deal with zero and negative values.
x = find(spect <= 0);
spect(x) = repmat(1e-50, length(x), 1);
