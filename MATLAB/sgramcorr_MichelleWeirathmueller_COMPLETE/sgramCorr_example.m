
%%%%%%%%%%%%%%%%%%%%%%%%%%% start of configuration %%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input file:
f = 'G03BBHZ.20121110-20121115';	% root of filename for Cascadia files
filename = ['C:\Dave\sounds\Cascadia\' f '\' f '.wav'];
sRate_override = 125;

% Spectrogram parameters for a frame size of ~4.0 s and overlap fraction of 1/4:
gramParams = struct(...
    'frameSizeS',  4.0,...   % amount of time-series data per gram time 
    ...                      % slice, in seconds
    'overlapFrac', 0.25, ... % amount of overlap between succesive frames
    'zeroPadFrac', 0.0, ...  % amount of zero-padding added to the frameSize 
    ...                      % data points before the FFT is calculated
    'winType', 'hanning');   % FFT window type

% Normalization parameters
normTime = 20.0;	% exponential decay time constant, seconds

% This defines a s'gram corrlation kernel with a single downsweep from 15.7 Hz
% to 14.4 Hz, with 10 s duration and bandwidth 0.5 Hz:
kerspec = [1  15.7  14.4  0  10  0.5];

%%%%%%%%%%%%%%%%%%%%%%%%%%%% end of configuration %%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Read in the samples and sample rate.
[x,sRate] = soundIn(filename);
sRate = sRate_override;

% Calculate gram. This also adds a bunch of fields to gramParams with values
% in samples instead of seconds; see davespect.m for details.
[gram,fRate,gramParams] = davespect(x, gramParams, sRate);

% Normalize using a decaying-exponential subtraction. This makes all background
% noise values close to 0 and makes transient animal calls stand out.
normParams = ['e' 0 0 0 normTime];	% see normGram.m for syntax
normG = normGram(gram, sRate, fRate, normParams);

% Calculate spectrogram cross-correlation. 
sgc = multiCorr(normG, kerspec, sRate, fRate, gramParams.frameSize, ...
    gramParams.zeroPad, gramParams.nOverlap, 1);
