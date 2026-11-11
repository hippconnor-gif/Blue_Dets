function res = morphomedian(im, ker, pctile)
%MORPHOMEDIAN	Find median value in each point's nbd in image (morphology)
%
% result = morphomedian(image, kernel)
%    Process an image with a median filter.  ker is a 2-column array
%    [deltaX deltaY] which specifies the shape of a neighborhood.
%    For each point at (x,y) in the image, replace the point by a
%    median taken over the points in the neighborhood:
%
%       result(x,y) = MEDIAN {image(x+deltaX(i), y+deltaY(i))}
%                      i=1:n
%
%    where n is the number of rows in ker.
%    Usually ker contains a row with [0 0], but it doesn't have to.
%    
%    If ker has more than one point in it, 'result' will be smaller 
%    than 'image'.
%
% result = morphomedian(image, kernel, pctile)
%    You can use a percentile value instead of the median (which is
%    the 50th percentile).  pctile should be between 0 and 1.
%
%    If there is no value exactly at the specified percentile, a
%    linear interpolation of the two adjacent values is used.

if (nargin < 3), pctile = 0.5; end

% Don't care where the origin is.  Make all coords zero-based.
kx = ker(:,1);  kx = kx - min(kx);  
ky = ker(:,2);  ky = ky - min(ky);  
nker = length(kx);

if (nker <= 1)
  res = im; 
  return;
end

transp = nRows(im) > nCols(im);
if (transp)
  im = im.';
  [kx,ky] = swap(kx,ky);
end

% This gets complicated only because the processing is broken into
% separate chunks in order to avoid out-of-memory problems.

chunksize = 500000;			  % max # of doubles to process at once
tileN = round(chunksize/nRows(im)/nker);  % vert. slices to process per chunk

% Create output array.
res = zeros(nRows(im) - max(ky), nCols(im) - max(kx));

for tileX = 1 : tileN - max(kx) : nCols(res)
  tile = im(:, tileX : min(tileX + tileN - 1, nCols(im)));
  nx = nCols(tile) - max(kx); xi = 1 : nx;
  ny = nRows(tile) - max(ky); yi = 1 : ny;
  
  % Construct res with one layer for each translation.
  res0 = zeros(nker, ny, nx);
  for i = 1 : nker
    res0(i,yi,xi) = tile(yi + ky(i), xi + kx(i));
  end
  
  res0 = percentile(res0, pctile);

  if (nRows(res) > 0)				% remove 1st dim, singleton
    res(:, tileX : tileX+nx-1) = reshape(res0,ny,nx);
  end
end

if (transp)
  res = res.';
end
