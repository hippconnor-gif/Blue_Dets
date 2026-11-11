% loadCompat      Load compatfile, correcting for changes in storage format.
%
% Provide backward compatibility for detector-output files (see crunch.m).
% This is a script because it loads a lot of variables into the current 
% context.
%
% Input:  compatname
% Output: whatever is in the file loaded.
%
% Changes, at various times, in storage format for detector output files.
%    fftSize	 replaced by dataSize and zeroPad
%    mean	 renamed to avg because of existing function named mean
%    hoffset	 added
%    processType added but not dealt with here yet.
%    A,B         changed to A1, B1 because of conflict with filter coefficients
%
% This module has to be a script because so many variables get loaded.

% First check for files created by the ~/sc/sc program.
scSuffixes = str2mat('finC','finQ','finR','finU','blueA','blueC','blueD',...
    'minkeB','minkeC','minkeD', 'blueT');
if (any(strmtch1(pathExt(compatname), scSuffixes)))
  [sRate, dataSize, zeroPad, hopSize, winType, fRate, hoffset, manyS]= ...
      loadSC(compatname);
  A1 = [];
  B1 = [];
  return
end

clear fftSize
clear mean
clear hoffset
if (exist('A')), Asave = A; clear A; 		% save filter coeffs
else clear Asave
end
if (exist('B')), Bsave = B; clear B; 
else clear Bsave
end

load(compatname, '-mat');

if (exist('fftSize'))
  dataSize = fftSize;
  zeroPad = 0;
  clear fftSize
end
if (exist('mean') == 1)
  avg = mean;
  clear mean;
end
x = exist('hoffset');		% MATLAB bug, do not combine these two lines
if (x ~= 1)
  suf = pathExt(compatname);
  hoffset = 0;
  if     (strcmp(suf, 'fin-C')), hoffset = 3.36;	% start of pulse
  elseif (strcmp(suf, 'fin-G')), hoffset = 3.32;	% start of pulse
  elseif (strcmp(suf, 'fin-I')), hoffset = 3.28;	% start of pulse
  end
end
if (exist('A')), A1 = A; 
else A1 = []; 
end
if (exist('B')), B1 = B;
else B1 = []; 
end
