function [spect,prev] = normExp(spect, binBW, fRate, prev, ...
    withPrev, scaleByMax, bandMax, decay1, decay2, fMin, fMax, noIncrease, ...
    scaleBySD)
% normExp		Exponential-decay (IIR) normalization.  See normGram.
%
% [spect,prev] = normExp(spect, binBW, fRate, prev, withPrev, scaleByMax, ...
%                  bandMax, decay1, decay2, fMin, fMax, noIncrease, scaleBySD)
%
% See normGram for an explanation.

%disp('normExp: Normalizing (IIR)...');

if (isnan(noIncrease)), noIncrease = 0; end
if (isnan(scaleBySD )), scaleBySD  = 0; end

loBin = max(1,            round(fMin / binBW));	% handles NaN okay!
hiBin = min(nRows(spect), round(fMax / binBW));	% handles NaN okay!

if (~scaleByMax)	% most common option: just subtract running average

  % Get the matrix or vector to run the exponential decay process on.
  s = spect;
  if (bandMax)
    if (loBin > hiBin)
      error('Bad (maybe backwards?) frequency bounds.');
    end
    s = max(s(loBin:hiBin, :, :), 1);
  end
  
  % Run the exponential decay process.
  s = permute(max(-10, s), [2 1 3]);
  if (withPrev)
    [s1,prev,sd] = expdecay(s, fRate, decay1, prev.', decay1*3, ...
	noIncrease, scaleBySD);
  else
    disp('Using non-warmed-up AGC.');
    if (noIncrease), error('Does not work with noIncrease ~= 0.'); end
    %[s1,prev,sd] = expdecay(s, fRate, decay1, 0); I think this is wrong
    [s1,prev,sd] = expdecay(s, fRate, decay1, []);
  end
  s1   = ipermute(s1,   [2 1 3]);
  sd   = ipermute(sd,   [2 1 3]);
  prev = ipermute(prev, [2 1 3]);
  
  % If we used the max, convert back to a matrix.
  if (bandMax)
    s1 = repmat(s1, [nRows(spect) 1 1]);
    sd = repmat(sd, [nRows(spect) 1 1]);
  end

  % Plot stuff for debugging.
  if (0 && nRows(spect) == 1)
    disp('normExp: Smoothing for the debug plot!')
    figure(2)
    t = (0 : nCols(spect)-1) / fRate;
    smoothSams = fRate * 2.0;          % time (in sams) to smooth over
    s1A = filter(ones(1, round(smoothSams)) / smoothSams, 1, s1);
    sdA = filter(ones(1, round(smoothSams)) / smoothSams, 1, sd);
    plot(t, spect, 'g-', t, s1A, 'b-', t, s1A+sdA, 'b:', t, s1A-sdA, 'b:');
    ylims fit+
    figure(1)
  end

  % At last!  Apply the mean-subtraction and maybe the sd-division.
  spect = spect - s1;
  if (scaleBySD)
    spect = spect ./ sd;
  end

else		% scaleByMax: subtract running avg, divide by running max

  disp('WARNING: THIS CODE HAS NOT BEEN CHECKED SINCE CONVERSION TO 3-D')
  disp('SPECTROGRAMS.')
  
  if (~bandMax)
    error(['The combination of scaling by a maximum process (scaleByMax=1) '...
	'and no frequency-band maximizing (bandMax=0) is not implemented.']);
  end

  if (~withPrev) 
    prev = []; 
  else
    disp('normExp: Warning: no warm-up processing.');
  end

  % Subtract exponential-decay value, divide by running maximum.
  if (length(prev) == 0)
    prev = [0 0 0]; 
  end
  alpha  = 1 - exp(-1 / (decay1 * fRate));	% decay per spect frame
  alpha2 = 1 - exp(-1 / (decay2 * fRate));	% decay per spect frame
  
  prev1 = prev(1);	% tracks the average value
  prev2 = prev(2);	% tracks the maximum value
  prev3 = prev(3);	% tracks the maximum value; slower transient response
  
  m = max(spect(loBin:hiBin, :));	% max value in band at each time step
  
  x = zeros(nCols(spect), 4, size(s,3));
  for i = 1:nCols(spect)
    spect(:,i,:) = max(0, (spect(:,i,:) - prev1) ./ (prev3 - prev1));
    prev1 = max(0,    (1-alpha) * prev1 + alpha * m(i), 1);
    prev2 = max(m(i), (1-alpha2) * prev2 + alpha2 * m(i), 1);
    prev3 = max([prev1, ((1-alpha2) * prev3 + alpha2 * prev2), m(i)], 1);
    x(i,1:3,:) = [prev1 prev2 prev3];
    x(i,4,:) = max(spect(loBin:hiBin,i), 1);
  end
  prev = [prev1 prev2 prev3];
  
  figure(4)
  for i = 1:nCols(x)
    subplot(nCols(x)*100 + 10 + i); plot((1:nRows(x))/fRate, x(:,i,:));
  end
  subplot(412); y = get(gca,'ylim');
  subplot(413); set(gca, 'ylim', y);
end
