function [sRate,dataSize,zeroPad,hopSize,winType,fRate,hoffset,manyS] = ...
    loadSC(filename)
%[sRate,dataSize,zeroPad,hopSize,winType,fRate,hoffset,manyS]= loadSC(filename)
%
% Load a file output by the C program sc, which is in ../sc/sc.c.
%
% Values returned:
%   sRate		sampling rate
%   dataSize		number of samples in each data frame
%   zeroPad		number of zeroes added to data frame before FFT
%   hopSize		distance in samples between data frames
%   winType		code for Hamming/Hanning/Blackman/etc.
%   fRate		frames per second of S
%   hoffset		filter "warm-up" time, in seconds (offset
%   			between the start of S and the start of the
%   			sound signal)
%   manyS		the detector output (row vector of detection scores)

[fd,msg] = fopen(filename, 'r');
if (fd < 0)
  error([filename, ': Can''t open file; error message was ''', msg, '''']); 
end

% global header
scVers = fread(fd, 1, 'float');
if (scVers ~= 1)
  fclose(fd);
  error([filename, ': Unknown sc.c version number ', num2str(scVers)]); 
end

% gram params
x = fread(fd, 9, 'float');
gcVers	 = x(1);
sRate	 = x(2);
dataSize = x(3);
zeroPad	 = x(4);
hopSiz	 = x(5);
winType	 = x(6);
fRate	 = x(7);
% x(8) and x(9) are ignored here
if (gcVers ~= 1)
  fclose(fd);
  error([filename, ': Unknown gramCorr.c version number ', num2str(gcVers)]); 
end

% kernel params; many elements of x are ignored here
x = fread(fd, 11, 'float');
ksVers = x(1);
if (ksVers > 3)
  fclose(fd);
  error([filename, ': Unknown kernel.c version number ', num2str(ksVers)]); 
end

% Mistake: Files labeled version 1 are really version 2 format!
typ = x(3);
if (ksVers < 3 & typ == 3)
  timeSpread = x(9);
  hoffset = timeSpread * 4.5;		% see kernel.c to explain the 4.5
else
  hoffset = x(2);
end
  
% the result
manyS = fread(fd, Inf, 'float').';
x = isinf(manyS);			% Infs in sound data cause problems
manyS(find(x)) = NaN * ones(1,sum(x));

fclose(fd);
