function A1 = multiCorr(spect,kerspec,sRate,fRate,dataSize,zeroPad,nOverlap,...
    showPlot)
%MULTICORR	Run spectrogram-correlation filter(s).
%
% A1 = multiCorr(spectrogram, kerspec, srate, fRate, dataSize, zeroPad, ...
%                                                          nOverlap, showPlot)
%
% Run spectrogram-convolution filter(s) on spect, a spectrogram.
%
% Arguments:
%  spectrogram	(MxN array) a spectrogram, with time on the horizontal axis
%  kerspec	(array) see below
%  srate	(scalar) sample rate, samples/s
%  fRate	(scalar) frame rate of spectrogram (rate of successive slices)
%  dataSize	(scalar) number of samples that went into making each slice
%  zeroPad	(scalar) number of zeros added to the dataSize samples before 
%                        the FFT is calculated
%  nOverlap	(scalar) number of samples of overlap in time series between 
%                        the start of successive slices
%  showPlot	(scalar) whether to plot the result; if >1, leave extra blank plot
%  A1 [output]	(row vector or array) output of spectrogram correlation, with
%		         one sample per spectrogram frame
%
%
% The input arg kerspec contains a series of rows specifying kernel parameters,
% each row having
%
%	type  start-freq  end-freq  start-time end-time  bandwidth
%
% that define each frequency sweep kernel.  Each kernel is cross-correlated
% with the spectrogram, yielding a time series with one sample per spectrogram
% frame.  The 'type' value should be '1' for each frequency sweep.
%
% Multiple kernels may be used, in which case multiple cross-correlations are
% done.  A matrix of several row vectors is returned, one per cross-correlation.
% The rows may be different lengths (if the kernels are different lengths), 
% in which case the short rows are zero-padded.
%
% To define multiple kernels, put each kernel definition into row(s)
% of kerspec; the kernels are separated by a row of kerspec that does
% NOT have a '1' in the first column (the type column).
%
% See also model, backtrack, multiKernel, munch, the old cFilter.

addExtraBlankPlot = (showPlot > 1);  % extra plot is for the combined value

n = nRows(kerspec);
A1 = [];

t = (kerspec(:,1) == 1);		% rows with a '1' type value
nsubplots = sum(diff([0; t]) == 1) + addExtraBlankPlot;
subplotnum = 1;
k = 1;
while (k <= n)
  if (t(k) == 1)
    % Find consecutive rows with a '1' type value.
    zeropos = find([t(k:n); 0] == 0);
    i = k : k+zeropos(1)-2;		% set of consecutive type '1' rows

    %disp(['    Kernel row(s) ', num2str(i), '...']);
    minT = min(kerspec(i,4));		      % for normalizing start-time to 0
    [ker,voff,hoff] = multiKernel(kerspec(i,2), kerspec(i,3), ...   % f0, f1
	kerspec(i,4)-minT, kerspec(i,5)-minT, kerspec(i,6), ...     % t0,t1,BW
	sRate, fRate, dataSize, zeroPad, nOverlap, 0);

    %%%%%%Remove negative lobes from kernel%%%%%%
    kerinds=find(ker < 0);
    %ker(kerinds)=0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    x = max(dumbConv(spect, ker, voff), 0);
    %x = [zeros(1, round(hoff * fRate))  x];                % before 10/04
    x = [zeros(1, round(hoff * fRate + nCols(ker)/2))  x];  % after 10/04
    
    % Plot it.
    if (showPlot)
      subplot(nsubplots, 1, subplotnum)
      subplotnum = subplotnum + 1;
      plotS(x, fRate)
      %if (k ~= 5), ylabel(' '); else disp('multiCorr: Hacking ylabel!'); end
      %title(['Kernel-filter ' int2str(i) ' output']);
      drawnow
      if (addExtraBlankPlot)
        subplot(nsubplots, 1, subplotnum)	% makes crunch's plot come out right
      end
    end
    
    % Expand result matrix A1 to be wide enough.
    if (~isempty(A1) & nCols(A1) < nCols(x))
      A1(1, nCols(x)) = 0;
    end
    
    A1 = [A1; x, zeros(1, nCols(A1) - nCols(x))];
    k = k + length(i);
  else
    % Not a kernel row -- ignore.
    k = k + 1;
  end
end
