% SSDB_DATA_IMPORT   Reads the scattering data of a folder
%
%   The reading can be restrictred in frequency and temperature by the two
%   first optional arguments.
%
%   The default behaviour is to perform a plain reading and then returning
%   *SSD* as an unsorted struct vector. *freq* and *temp* are in this case
%   also vectors, having the same length as *SSD* (and *nodata* contains
%   just false).
%
%   By setting *grid_sort* to true, the scattering data are sorted in terms
%   of frequency (row dimension) and temperature (column dimension), and
%   *SSD* is returned as structure matrix. For example SSD(2,3) is the data
%   for frequency 2 and temperature 3. *freq* and *temp* are in this case
%   returned as the frequency and temperature grid of the data. If any
%   element of SSD is unfilled, the corresponding element in *nodata* is
%   set to true.
%
% FORMAT [SSD,freq,temp,nodata] = ssdb_data_import(data_folder,
%                                                    [freq_range,temp_range,grid_sort])
%
% OUT SSD           Read single scattering data.
%     freq          The frequency of each SSD element.
%     temp          The temperature of each SSD element.
%     nodata        Boolean vector/matrix, flagging empty SSD elements.
%  IN data_folder   Full path to folder holding scattering netcdf data files
% OPT freq_range    Frequency limits. Ignore data outside these limits.
%                   Default is [0,Inf]. If set to [], the default is used.
%     temp range    Temperature limits. Ignore data outside these limits.
%                   Default is [0,Inf]. If set to [], the default is used.
%     grid_sort     Boolean to to trigger that data are pre-sorted according
%                   to frequency and temperature. See further above. Default is false.

% 2016-11-02 Patrick Eriksson
% 2016-11-03 Robin Ekelund: l97: changed to ssdb_read_ncfile.m


function [SSD,freq,temp,nodata] = ...
                     ssdb_data_import( data_folder, freq_range, temp_range, grid_sort )
%
if nargin < 2  |  isempty(freq_range)
  freq_range = [ 0 Inf ];
end
if nargin < 3  |  isempty(temp_range)
  temp_range = [ 0 Inf ];
end
if nargin < 4  |  isempty(grid_sort)
  grid_sort = false;
end
%
if ~isnumeric(freq_range)  |  length(freq_range) ~= 2
  error( '*freq_range* must be a numeric vector of length 2.*' );
end
if ~isnumeric(temp_range)  |  length(temp_range) ~= 2
  error( '*temp_range* must be a numeric vector of length 2.*' );
end


% Locate existing files and find index of ones inside given freq and temp
% ranges
%
[data_files,freqs,temps] = ssdb_data_files( data_folder );
%
ind = find( freqs >= freq_range(1)  &  freqs <= freq_range(2)  & ... 
            temps >= temp_range(1)  &  temps <= temp_range(2) );


% Return [] of no data at all
%
if isempty(ind)
  [SSD,freq,temp,nodata] = deal( [] );
  return
end


% Without grid sorting
%
if ~grid_sort
  %
  freq   = freqs(ind);
  temp   = temps(ind);
  nodata = repmat( false, size(freq) );
  %
  for i = 1 : length(ind)  
    SSD(i) = ssdb_read_ncfile( data_files{ind(i)} );
  end
  
% With grid sorting
%
else
  %
  freq = unique( freqs(ind) );
  temp = unique( temps(ind) );
  %
  nodata = repmat( false, [length(freq),length(temp)] );
  
  for f = 1 : length(freq)
    for t = 1 : length(temp)
      
      i = find( freq(f) == freqs  &  temp(t) == temps );

      if length(i) == 0
        nodata(f,t) = true;
      elseif length(i) == 1
        SSD(f,t) = ssdb_read_ncfile( data_files{i} );
      else
        error( 'Duplicate of one frequency-temperature found!' );
      end
    end
  end
end

