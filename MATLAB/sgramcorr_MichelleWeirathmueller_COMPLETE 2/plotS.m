function plotS(S, fRate, filename, noRect, varargin)
% plotS(S [,fRate [,filename [,noRect [,plotargs]]]])
%    Plot the results of the transient detector.
%    The filename is used for generating an x-axis label
%    fRate (the frame rate) defaults to 1.
%    noRect, if present and non-zero, says don't do halfwave rectification.
%    args after noRect are passed to plot (q.v.).
%
% plotS(filename)
%    Load and plot the detector output from the file.

if (nargin < 4), noRect = 0; end

if (ischar(S))
  filename = S;
  compatname = filename; loadCompat;
  S = manyS;
  if (~exist('fRate', 'var'))
    fRate = 1;
  end
else
  if (nargin < 2), fRate = 1; end
  if (nargin < 3), filename = ''; end
end

if (nRows(S) > nCols(S))
  S = S.';
end

% Use only the last row of S.
S = S(nRows(S), :, :);
S = shiftdim(S, 1);               % eliminate first dimension
time = repmat((0:nRows(S)-1).' / fRate, 1, size(S,2));
if (noRect)
  plot(time, S, varargin{:});
else
  plot(time, max(0, S), varargin{:});	% halfwave-rectified
end
set(gca, 'xlim', [0 nRows(S)/fRate], 'tickdir', 'out');
set(xlabel(filename), 'Interpreter', 'none');	% filename may be empty
xlabel('time, s')
ylabel('Detection score');
orient landscape

ylims('fit');				% may get overridden
if (~isempty(filename))
  ext = pathExt(filename);
  if     (strmtch1(ext, 'comma', 'blue-A', 'blueA')) axis([sub(axis,1:3) 100]);
  elseif (strmtch1(ext, 'fin-C', 'finC'))            axis([sub(axis,1:3)  30]);
  elseif (strmtch1(ext, 'minke-B', 'minkeB'))        axis([sub(axis,1:3) 150]);
  end
end

drawnow
