function [spect, prev] = normGram(spect, sRate, fRate, np, prev)
%normGram		Normalize a spectrogram, i.e., apply gain control.
%
% Normalize a spectrogram, i.e., scale it so that its absolute levels
% fall into some desired range.  This is usually performed on
% logarithmic, not quadratic, spectrograms, though either kind will work.
%
% Different types of normalization are performed depending on the 
% value of the normParams argument. All types use at least these arguments
% (though some use more; see 'prev'):
%
%     spectOut = normGram(spectIn, sRate, fRate, normParams)
% 
% where spectIn is the input spectrogram, spectOut is the result, 
% sRate is the signal sampling rate, fRate is the spectrogram frame rate,
% and normParams specifies the normalization to perform.
%
% Some parts of this routine can also be used on time series.  See the
% descriptions of the individual types of normalization for details.
% For time series, use the sampling rate for fRate and use NaN for sRate.
% (NB: Most types of normalization have not yet been tested with time series.)
%
% normParams is usually a row vector.  If it has more than one row, then
% it specifies successive passes of normalization, one after another.
% In each row, the first column is a letter specifying what general type of
% normalization algorithm to use, and the remaining values are usually 
% parameters used by the algorithm.  Here are the types:
%
% normParams = ['p' perFreq freqMin freqMax pct0 pct1]
%   Percentile normalization, also called low order statistics.  This is
%   useful for short-duration spectrograms to remove a background level
%   of white noise.  It doesn't work so well for long-duration spectrograms
%   because the sort involved in the percentile computation takes too long.
%   
%   Operation with perFreq=0: Compute the spectrogram, then extract the 
%   frequency range from freq0 to freq1.  Over all cells in this block,
%   compute the values at the pct0'th and pct1'th percentiles; call these
%   v0 and v1.  pct0 and pct1 should be between 0 and 1 inclusive.
%   Subtract v0 from the spectrogram, then scale it linearly by 1/(v1-v0).
%   Thus values in the original spectrogram at the pct0'th and pct1'th
%   percentiles get transformed to 0 and 1, respectively.  (Often after 
%   this it's useful to change spectrogram negative values to zero.)
% 
%   With perFreq=1, the operation is similar, except that the processing
%   happens on each frequency bin of the spectrogram (each horizontal line)
%   independently.
% 
%   If pct0 is a NaN, no subtraction happens (i.e., v0 is set to 0).  
%   If pct1 is a NaN, no scaling (no division by v1) happens.
%
%   Works on spectrograms, and on time series if perFreq is 1 (!).
% 
% normParams = ['p' perFreq freqMin freqMax pct0 pct1 mul0 mul1]
%   Percentile normalization can be done with two additional
%   multipliers.  The input spectrogram value that gets translated 
%   to an output value of 0 is 
%             v0 * mul0
%   (instead of just v0) and the input value that is scaled to 1 is
%             v1 * mul1
%   (instead of just v1).
%
%   Again, mul0 or mul1 may be NaN, in which case it is ignored.
%
%   Works on spectrograms, and on time series iff perFreq is true (!).
% 
% normParams = ['m' medianWidth pctile]
%   Apply a median filter with a horizontal-strip kernel.  Useful for
%   removing snapping shrimp noise.  medianWidth is in seconds.  It
%   doesn't actually have to be a median value (50th percentile); it can
%   be any percentile, specified by pctile.  pctile defaults to 0.50.
%
%   Works on time series.
%
% normParams = ['s' freqMin freqMax pctile mult]
%   For each frame (vertical slice) in the spectrogram:  Extract the
%   spectrogram cells between freqMin and freqMax.  Compute the value V
%   at the pctile'th percentile.  Subtract V*mult from each cell
%   in the frame.  Note that the power is subtracted in linear space,
%   not log space.  If mult is missing, 0, or NaN, it defaults to 1.
%
%   Works only with spectrograms, not with time series.
%
% normParams = ['e' withPrev scaleByMax bandMax decay1 decay2 freqMin freqMax noIncrease]
%   IIR exponential-decay normalization.  This is useful for long-duration
%   spectrograms to remove constant-frequency noise like 60-Hz hum or boat 
%   noise, and to whiten the background noise (Wiener filter).
%
%   The parameters up to and including decay1 are always required.  The
%   ones after that are used only by some combinations of flags (see 
%   below); parameters are positional, so if you need, say, the last
%   parameter but not some of the ones before it, put in NaN's for the 
%   unused parameters.
%
%   At each time frame in the gram, the value in each frequency bin
%   is added to an average that decays exponentially with a 1/e time 
%   constant given by decay1 -- that is, an impulse decays to 1/e of
%   its original value in decay1 seconds.  The running total is 
%   subtracted from each corresponding frame of spect.
%
%   Works with time series; use bandMax=0 and freqMin=freqMax=NaN.
%
%   Four flags affect the algorithm:
%
%   withPrev: If this flag is non-zero, initialize the IIR decay
%       process with the 'prev' argument.  In this case the call to
%       normGram is like this:
%
%       [spectOut,prev] = normGram(spectIn, sRate, fRate, normParams, prev)
%
%       If bandMax is 0, then prev should be a column vector with as
%       many rows as spect.  If bandMax is 1, then prev should be a
%       scalar.
%
%       Zeros are used if prev is empty, so you can initialize prev to
%       [] and always use the same call to normGram. prev is returned,
%       so that it may be passed in on the next call; this allows long
%       spectrograms to be processed in separate, successive segments.
%
%   scaleByMax: If this flag is non-zero, then a second exponentially-
%       decaying process is run simultaneously with the first.  This
%       second process has a time constant of decay2.  This process
%       tracks the MAXIMUM value, and the usual AVERAGE value is
%       divided by this maximum.  Thus not only is the noise floor
%       removed, but the value of the loudest sounds are normalized
%       to 1.  If bandMax=0, the maximum of each frequency bin is
%       tracked (and prev, if present, has three columns instead 
%       of just one), and if bandMax=1, the maximum over all the
%       freq bins is tracked (and prev is a 3-element row vector).
%
%       decay2 is usally a much larger time constant (by a factor 
%       of 10, say) than decay1.  If scaleByMax is 0, then decay2
%       is not used.
%
%   bandMax: If this flag is non-zero, the exponentially-decaying value
%       used in the computation is the maximum over all the frequencies
%       in the band from freqMin to freqMax.  If this flag is zero (the
%       usual case), then each frequency bin is treated separately, as
%       its own exponentially decaying process; in this case freqMin 
%       and freqMax are not used (and need not be specified).  For time
%       series, set bandMax to 0.
%
%   noIncrease
%       If this flag is non-zero, then the exponentially-decaying
%       process is never allowed to increase, only to decrease.
%       This is a hack, added to help remove snapping shrimp noise.
%
%   NB!  The combination of scaleByMax=1 and bandMax=0 is not
%   implemented yet (see normExp).
%
% normParams = ['v' ... ]
%   Vertical processing.  Gain control is applied along single time 
%   slices (vertical strips) instead of horizontal frequency bins.
%   Following the 'v' should be another letter indicating what type
%   of processing to do, like this:
%           normParams=['ve' withPrev ... ]
%   All frequency parameters should now become time parameters, and
%   vice versa.  For instance, decay1 should specify the number of 
%   hertz over which the decay process falls to 1/e of its value.
%   Normalizations that are usually dependent on time direction, like 
%   normExp, operate vertically from the top down (high frequency to low).
%
%   Does not work with time series.
%
% normParams = ['z' minval maxval subval]
%   Zeroing out.  Change all values in the gram less than minval to minval.
%   Change all values greater than maxval to maxval.  Then subtract subval.
%   If minval, maxval, or subval is NaN, it is ignored.
%
%   Works with time series.
%
% See also expdecay, normPct, normExp, normMedian, normSub.

if (nargin < 5), prev = []; end

% Handle multi-row np's recursively.
if (nRows(np) > 1)
  for i = 1 : nRows(np)
    if (nargout == 1)
      spect = normGram(spect, sRate, fRate, np(i,:));
    else
      [spect,prev] = normGram(spect, sRate, fRate, np(i,:), prev);
    end
  end
  return
end

typ = char(np(1));
np = double(np);
binBW = sRate/2 / nRows(spect);

if (typ == 'p')			% percentile normalization
  if (~ismatrix(spect)), error('Not ready for multi-track spect yet.'); end

  % Order of values in np:      'p' perFreq fMin fMax pct0 pct1 mul0 mul1
  np(length(np)+1 : 8) = NaN;	% use NaN for missing values
  spect = normPct(spect, binBW, np(2), np(3), np(4), np(5), np(6),np(7),np(8));

elseif (typ == 'm')		% horizontal median
  if (~ismatrix(spect)), error('Not ready for multi-track spect yet.'); end
  
  % Order of values in np:	'm' timeSpan pctile
  np(length(np)+1 : 3) = NaN;	% NaN makes normMedian use median
  spect = normMedian(spect, fRate, np(2), np(3));

elseif (typ == 's')		% subtract average power in a time slice
  if (~ismatrix(spect)), error('Not ready for multi-track spect yet.'); end
  
  % Order of values in np:	's' freqMin freqMax pctile [mult]
  if (length(np) < 5), mult = 1;
  elseif (np(5) == 0 || isnan(np(5))), mult = 1;
  else mult = np(5);
  end
  spect = normSub(spect, binBW, np(2), np(3), np(4), mult);

elseif (typ == 'e')		% IIR exponential decay normalization.
  % Order of values in np:
  %  ['e' withPrev scaleByMax bandMax decay1 decay2 fMin fMax noIncr scaleBySD]
  %
  np(length(np)+1 : 10) = NaN;        % use NaN for missing values
  [spect,prev] = normExp(spect, binBW, fRate, prev, ...
      np(2), np(3), np(4), np(5), np(6), np(7), np(8), np(9), np(10));

elseif (typ == 'v')		% vertical processing

  % WARNING: THIS CODE HASN'T BEEN CHECKED SINCE CONVERSION TO 3-D
  % SPECTROGRAMS.
  % All we do here is adjust the parameters appropriately, transpose spect, 
  % and make a recursive call.
  newSRate = (1/fRate) * nCols(spect) * 2;
  newFRate = 1/binBW;
  % Flip order of first 2 dims; then do fliplr so exp decaying happens top-down.
  % Before multi-channel spects, this was "spect = fliplr(spect.')".
  % Do NOT use prev, which is used only for horizontal processing.
  spect = flipdim(permute(spect,[2 1 3]), 2);
  spect = normGram(spect, newSRate, newFRate, np(2:length(np)));
  spect = ipermute(flipdim(spect, 2), [2 1 3]); % was "spect = fliplr(spect).'"
 
elseif (typ == 'z')		% zero out gram (apply min and/or max)
  
  % Order of values in np:	'z' minval maxval subval
  if (~isnan(np(2)))
    ix = (spect < np(2));
    spect(ix) = np(2);
  end
  if (length(np) >= 3 && ~isnan(np(3)))
    ix = (spect > np(3));
    %spect(ix) = repmat(np(3), sum(sum(ix)), 1);
    spect(ix) = np(3);
  end
  if (length(np) >= 4 && ~isnan(np(4)))
    spect = spect - np(4);
  end


else
  error('Unknown normalization type: %d == ''%c''', typ, typ)
end
