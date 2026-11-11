function hIndex = readHarufileDateIndex(indexfilename)
%readHarufileDateIndex    Extract start-time of Haru files from file_dates file
%
% hIndex = readHarufileDateIndex(indexfilename)
%   Given an indexfilename like 'file_dates_Goa99_N50W135.txt' or 'Bilan*.txt'
%   or 'Bilan_969F_20071022200715.txt' that contains the starting time of
%   a list of Haru-format files or AURAL or .wav files, read the lines in the
%   file one by one and get the start time of each file. The return value
%   hIndex is a structure:
%       time      vector containing the start time of each file encoded as
%                 in datenum (q.v.). The time is the time extracted directly
%                 from the Haru file, which should be GMT.  If a file is
%                 missing in the index (sometimes sequence numbers are
%                 missing), its time is NaN.
%       hname     cell vector with the name of each Harufile, as given in the
%                 index file; RELATIVE TO DIRECTORY OF INDEX FILE in fname
%	nSamples  vector with the number of samples in each file
%       dur       (scalar) duration (in s) of each file; same for all files
%       fname     the index file name (the input file); useful for error 
%                 reporting
%
%   This routine is a lot like readBilanIndex (q.v.), and uses that one
%   when called with a 'Bilan*' file.
%
% [This is the new version of readHarufileDateIndex, the one that stores
% information about multiple disks sequentially in the return variables
% and includes the 'diskN/' string in hname. There is no disk number anymore.]
% 
% See also makeHarufileDateIndex, datenum, harufileFromDate, haruSoundAtTime,
%   goaspermdates, readDetLog, readBilanIndex.

if (strindex(indexfilename, 'Bilan') > 0)
  hIndex = readBilanIndex(indexfilename);
  return
end

fd = fopen(indexfilename);
if (fd < 0)
  error('Can''t open file_dates file %s.', indexfilename); 
end

filetimes = zeros(0,1);
nSamples  = zeros(0,1);
hname     = cell(0,1);
foundNSam = false;
fnum = 0; % numbering used to be taken from the file name, but now all files 
          % are just put into the index sequentially.
while (1)
  % Read next line and check that it's formatted right, not blank, etc.
  ln = fgetl(fd);
  if (length(ln) == 1 && ln == -1), break; end      % EOF?
  if (length(ln) <= 20 || ~any(strcmp(ln(1:2), {'19' '20'}))) % 20XX 19XX
    continue
  end
  
  % Store values in filestimes, nSamples, hname. First ensure there's space.
  fnum = fnum + 1;
  while (fnum > nRows(filetimes))   % ensure filetimes has enough rows
    filetimes(nRows(filetimes) + (1:1000), :) = NaN;
    nSamples (nRows(nSamples)  + (1:1000), :) = NaN;
    hname(    nRows(hname)     + (1:1000), :) = {'none'};  %...hname too
  end
  [t,~,~,nextI] = sscanf(ln, '%4d-%3d (%*6c) %d:%d:%g', 5); % extract date/time
  filetimes(fnum) = datenum(t(1), 1, t(2), t(3), t(4), t(5));
  [hname{fnum},~,~,nextI2] = sscanf(ln(nextI:end), '%s', 1);
  
  % If number of samples is provided, get it too.
  nSam = sscanf(ln(nextI+nextI2:end), '%d', 1);
  if (~isempty(nSam))
    nSamples(fnum) = nSam;
    foundNSam = true;
  end
end		% end while (1)
fclose(fd);

% Trim vectors to the actual number of files found.
filetimes = filetimes(1 : fnum, :);
nSamples  =  nSamples(1 : fnum, :);
hname     =     hname(1 : fnum, :);

% File duration: mean of all durations that are within 2 s of median duration.
d = diff(filetimes(:)) * 24*60*60;		% durations, in seconds
d = d(~isnan(d));
filedur = mean(d(abs(d - median(d)) < 2));

% Construct the return value.
hIndex = struct( ...
  'time',	filetimes, ...
  'dur',	filedur, ...
  'hname',	{hname}, ...
  'fname',	indexfilename);
if (foundNSam)
  hIndex.nSamples = nSamples;
end
