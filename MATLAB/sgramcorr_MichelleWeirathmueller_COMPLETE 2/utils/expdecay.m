function [x,prev,sd] = expdecay(x, sRate, decayTime, prev, warmTime, ...
  noIncrease, scaleBySD)
%EXPDECAY	Compute an exponentially-decaying average.
%
% y = expdecay(x, sRate, decayTime)
%    Given a vector x, compute an exponentially decaying average of x.
%    Note that the average itself is returned, not x minus the average.
%    This is the running average   y(i) = alpha * y(i-1) + (1-alpha) * x(i)
%    where alpha is a decay constant.  alpha is set such that a unit impulse
%    will decay to 1/e in time decayTime.  Before the averaging process is
%    started, the filter is "warmed up" for decayTime*3 seconds.
%    On output, y is a column (!) vector the same length as x.
%
%    If x is an array, do the averaging separately on each column (!) of x.
%
% [y,prev] = expdecay(x, sRate, decayTime, prev, warmTime)
%    If prev is used as an input argument and is non-empty, use it as the
%    initial value of the average y.  When prev is an output argument, 
%    it is the value to use next time around.  If it's empty, use warmTime
%    seconds of the signal to initialize the averaging process.  warmTime
%    defaults to decayTime*3.
%
%    If x is an array, prev is a row vector.  warmTime is still scalar.
%
%    The usual way to run this function on several successive chunks of
%    samples from one long signal is to use prev=[] for the first call, 
%    then use whatever prev value was returned on succeeding calls.
%
% [y,prev] = expdecay(x, sRate, decayTime, prev, warmTime, noIncrease)
%    If noIncrease is non-zero, the decay process is never allowed to
%    increase -- it can only decrease.  noIncrease defaults to 0.  (This
%    is a total hack invented for snapping shrimp noise removal.)
%
% [y,prev] = expdecay(x, sRate, decayTime, prev, warmTime,noIncrease,scaleBySD)
%    If scaleBySD is present and non-zero, also divide the output by the
%    running standard deviation.  In this case prev has 2 rows, the usual first
%    one with the running mean and a second one with the running variance.
%
% See also normExp, normGram.

% Since x may be quite large, the calculation is done in place.  x's rows
% get replaced by output values as processing proceeds.

if (nargin < 4), prev = []; end
if (nargin < 5), warmTime = decayTime * 3; end
if (nargin < 6), noIncrease = 0; end
if (nargin < 7), scaleBySD = 0; end

if (nRows(x) == 1 && ~(nargin >= 4 && nCols(prev) == nCols(x)))
  x = permute(x,[2 1 3]);            % make into col vector
end

nsamp = nRows(x);

reshapedX = (ndims(x) == 2);
if (reshapedX)
  x    = reshape(x,    [size(x)    1]);
  prev = reshape(prev, [size(prev) 1]);
end

if (isempty(prev) && warmTime > 0)
  % Warm up for n samples to initialize prev.  May need multiple iterations.
  nwarm = round(warmTime * sRate);			% # of warmup frames
  chunk = x(1 : max(1, min(nsamp, nwarm)), :, :);	% what to warm up on
  prev = mean(chunk, 1);				% set up mean
  if (scaleBySD)
    % Also set up variance in second row of prev.
    prev = [prev; mean((chunk - repmat(prev,[nRows(chunk) 1 1])) .^ 2, 1)];
  end
  for i = 1 : ceil(nwarm / nsamp)
    [~,prev] = expdecay(chunk, sRate, decayTime, prev, 0, ...
	noIncrease, scaleBySD);
  end
end

alpha = 1 - exp(-1 / (decayTime * sRate));	% decay per sample

runMn = prev(1,:,:);		% running mean
if (scaleBySD && nRows(prev) > 1)
  runVar = prev(2,:,:);		% running variance
  sd = zeros(size(x));		% output value -- std. dev.
else
  sd = [];			% still have an output
end

if (~noIncrease && ~scaleBySD && exist('filter', 'builtin'))
  % This way is faster, but requires the DSP toolbox.
  runMnSize = size(runMn);
  [x,runMn] = filter([alpha 0], [1 -1+alpha], x, runMn, 1);
  runMn = reshape(runMn, runMnSize);	% filter leaves it the wrong shape!
else
  % This way works in the absence of filter().
  ix = ones(1, size(x,2), size(x,3));	% gets overridden if noIncrease is true
  nIx = numel(ix);			% ditto
  for i = 1 : nRows(x)
    if (noIncrease)
      ix = (x(i,:,:) < runMn);		% override!
      nIx = sum(sum(ix));
    end
    
    % Calculate running mean and (optionally) running variance.
    if (nIx > 0)
      runMn(ix) = (1-alpha) * runMn(ix) + alpha * x(i,ix,:);
      if (scaleBySD)
        runVar(ix) = (1-alpha)*runVar(ix) + alpha * (x(i,ix,:) - runMn(ix)).^2;
      end
    end
    if (scaleBySD), sd(i,:,:) = sqrt(runVar); end
    x(i,:,:) = runMn;
  end
end

if (~scaleBySD), prev = runMn; 
else prev = [runMn; runVar];
end

if (reshapedX)
  if (nargout >= 1), sz = size(x   ); x    = reshape(x,    sz(1:2)); end
  if (nargout >= 2), sz = size(prev); prev = reshape(prev, sz(1:2)); end
  if (nargout >= 3), sz = size(sd  ); sd   = reshape(sd,   sz(1:2)); end
end
