function hIndex = readBilanIndex(indexfilename)
%readBilanIndex    Extract start-time of AURAL files from 'Bilan*' file
%
% hIndex = readBilanIndex(indexfilename)
%   Given an indexfilename like 'Bilan_969F_20071022200715.txt' that contains
%   the starting time of a list of AURAL files, read the lines in the
%   file one by one and get the start time of each file.  The return
%   value hIndex is a structure:
%       time      Fx1 array containing the start time of each file encoded as
%                 in datenum (q.v.); the array is indexed by
%                 F=(file number + 1), since file numbers start at 0.
%                 The time is the time extracted directly from the Bilan file,
%                 which should be GMT.  If a file is missing in the index 
%                 (sometimes sequence numbers are missing), its time is NaN.
%       hname     Fx1 cell array: the name of each AURAL soundfile, as given in
%                 the index file; RELATIVE TO DIRECTORY OF INDEX FILE
%       dur       (scalar) duration (in s) of each file; same for all files
%       fname     the input file name (the index file); useful for error 
%                 reporting
%	nSamples  set to NaN
%
%   This routine is a lot like readHarufileDateIndex (q.v.).
%
% See also datenum, harufileDateIndex, goaspermdates, atlminkedates,
%   readDetLog.

% Differences from harufileDateIndex:
%   - format of the string on each line
%   - duration is calculated from the index file (which has end-times in it),
%     not by calculating the difference between start-times of successive files

if (strindex(indexfilename, 'file_dates') > 0)
  hIndex = readHarufileIndex(indexfilename);
  return
end

realindexname = indexfilename;
fd = fopen(realindexname);
if (fd < 0)
  % Bilan file names are allowed to have wildcards.  Try this.
  if (strindex(indexfilename, '*') > 0)		% is there a wildcard?
    filelist = dir(indexfilename);		% dir() allows wildcards
    % Use only if there's exactly 1 match, and it's a file.
    fn = [pathDir(indexfilename) filesep filelist(1).name];
    if (length(filelist) == 1 && exist(fn, 'file'))
      realindexname = fn;
      fd = fopen(realindexname);
    elseif (length(filelist) > 1)
      error('%s: The index name %s resolves to more than one file:\n%s', ...
	  mfilename, indexfilename, char(filelist.name));
    end
  end
  if (fd < 0)
    error('Can''t open "Bilan" file %s.', realindexname); 
  end
end

filetimes = zeros(0,1);
hname     = cell(0,1);
durs      = zeros(0,1);		% duration of each file
nsams     = zeros(0,1);		% number of samples in each file
maxFnum1 = 0;
count = 0;
while (1)
  % Read next line and check that it's formatted right, not blank, etc.
  ln = fgetl(fd);
  if (length(ln) == 1 && ln == -1), break; end      % check for EOF
  if (length(ln) > 20)
    % First find separators -- semicolons or tab characters.
    semis = find(ln == ';');		% semicolons
    isSemi = true;
    if (length(semis) < 12 || semis(1) < 12)
      semis = find(ln == 9);		% tab characters
      isSemi = false;
      if (length(semis) < 12 || semis(1) < 14)
        error('Index file %s: Can''t find sound file name at this line:%c%s', ...
          realindexname, 10, ln)
      end
    end
    s1 = semis(1);
    if (~(strcmpi(ln(s1-4 : s1-1), '.WAV') || strcmpi(ln(s1-4 : s1-1), '.WVE'))...
        || any(~isdigit(ln(s1-8 : s1-5))))
      error('File %s: No file number at this line%c%s', realindexname, 10, ln)
    end

    % Get start of filename.
    fn0 = find(~isletter(ln(1:s1-5)) & ~isdigit(ln(1:s1-5)), 1, 'last') + 1;
    
    fnum = sscanf(ln(fn0+4 : s1-5), '%d');
    % Handle the case in which file indexes wrap around at 10000.
    while (fnum+1 <= nRows(filetimes) && ~isnan(filetimes(fnum+1)))
      fnum = fnum + 10000;
    end

    % Get number of samples.
    [x,nscan] = sscanf(ln(semis(8)+1 : semis(9)-1), '%d');
    nSam = iff(nscan == 1, x, NaN);
    
    if (isSemi)
      date0 = find(ln(      s1:end) == '/', 1) +       s1-1 - 4;
      date1 = find(ln(date0+10:end) == '/', 1) + date0+10-1 - 4;
    else
      date0 = find(ln(s1:end) == '/', 1) + s1-1 - 2;
      date1 = find(ln(date0+3:end) == ':', 2) + date0 + 4;
      date1 = date1(2);
    end
    if (date0 < s1 || date1 < date0)
      error('Index file %s: Can''t find date string at this line:%c%s', ...
	  realindexname, 10, ln)
    end

    % For valid lines, extract date and time, store in filetimes.
    if (~isnan(fnum))
      count = count + 1;
      if (isSemi)
        [t,cnt] = sscanf(ln(date0:date1), '%d/%d/%d;%f:%f:%f;%f');
      else
        [t,cnt] = sscanf(ln(date0:date1), '%d/%d/%d %f:%f:%f %f');
        t = t([3 1 2 4 5 6]);       % rearrange into yr-mon-day-hr-min-sec
      end
      if (cnt < 7)
	error('Index file %s: Bad format date string at this line:%c%s', ...
	    realindexname, 10, ln)
      end
      dn = datenum(t(1), t(2), t(3), t(4), t(5), t(6) + t(7)/1000);

      if (isSemi)
        [t1,cnt] = sscanf(ln(date1:end), '%d/%d/%d;%f:%f:%f');
      else
        [t1,cnt] = sscanf(ln(date1+3:end), '%d/%d/%d %f:%f:%f');
        t1 = t1([3 1 2 4 5 6]);
      end
      if (cnt < 6)
	error('Index file %s: Bad format date string at this line:%c%s', ...
	    realindexname, 10, ln)
      end
      dn1 = datenum(t1(1), t1(2), t1(3), t1(4), t1(5), t1(6));

      % Store results.  First ensure the return vectors have enough room.
      while (fnum + 1 > nRows(filetimes))   % ensure filetimes has enough rows
	filetimes(nRows(filetimes) + (1:1000), :) = NaN;
	hname(    nRows(hname)     + (1:1000), :) = {'none'};  %...hname too
	durs(     nRows(durs)      + (1:1000), :) = NaN;
	nsams(    nRows(nsams)     + (1:1000), :) = NaN;
      end
      filetimes(fnum+1) = dn;      % +1 since datafile numbers start at 0000
      hname{fnum+1} = ln(1:s1-1);
      durs(fnum+1) = dn1 - dn;
      maxFnum1 = max(maxFnum1, fnum + 1);
      nsams(fnum+1) = nSam;
    end
  end		% end if (~isnan(fnum))
end		% end while (1)
filetimes = filetimes(1 : maxFnum1, :);
hname     = hname(    1 : maxFnum1, :);
durs      = durs(     1 : maxFnum1, :);
nsams     = nsams(    1 : maxFnum1, :);

fclose(fd);

% File duration: Floor of mean of durations that are within 1 s of median.
d = durs(~isnan(durs)) * 24*60*60;		% non-NaN durations, in seconds
filedur = floor(mean(d(abs(d - median(d)) < 1)));

hIndex = struct('time', filetimes, 'dur', filedur, 'hname', {hname}, ...
  'fname', realindexname, 'nSamples', nsams);

