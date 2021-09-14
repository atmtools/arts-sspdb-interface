% SSDB_IMPORT   Reads and compiles data from SSDB database for one habit.
%
%   Default is to read all data for specified habit and orientation
%   combination, but the reading can be restricted in terms of size, frequency
%   and temperature. For frequency and temperature, data are selected as: 
%      limit1 <= data <= limit2
%   while for size, data at the lower limit is excluded:
%      limit1 < data <= limit2.
%
%   Default is to issue an error as soon as data are missing for a frequency
%   and temperature combination. Allow missing data by setting the optional
%   argument *allow_nodata* to true. Note that set of frequencies and
%   temperatures are allowed to differ between sizes, independently of how
%   *allow_nodata* is set. For example, for one size there can be a single
%   temperature, while other sizes have a temperature grid with several
%   elements.
%
% FORMAT  [ data ] = ssdb_import_habit( habit_id, orientation, ...
%                    [allow_nodata, size_range, size_type, freq_range, temp_range] )
%
% OUT data          Struct array that contains SSD data.
%  IN habit_id      Habit id nr.
%     orientation   String describing orientation.
% OPT allow_nodata  See above. Default is false.
%     size_range    Size limits. Ignore data outside these limits.
%                   Default is [0,Inf]. If set to [], the default is used.
%     size_type     Quantity used for size cropping. Allowed options are
%                   'dveq', 'dmax' and 'mass'. Default is 'dveq'.
%     freq_range    Frequency limits. Ignore data outside these limits.
%                   Default is [0,Inf]. If set to [], the default is used.
%     temp range    Temperature limits. Ignore data outside these limits.
%                   Default is [0,Inf]. If set to [], the default is used.
%
% 2017-09-05 Robin Ekelund: Created
% 2017-09-08 Robin Ekelund: Bug fix. Couldn't handle missing data.
% 2017-10-20 Robin Ekelund: Changed to account for changes in other
%                               functions
% 2017-10-20 Robin Ekelund: Minor fix. Did not account for empty data in
%                               ssdb_read_ncfile.
% 2017-10-30 Robin Ekelund: Removed various superfluous code.
% 2018-02-13 Robin Ekelund: Fixed bug when allow_nodata is allowed.
%
function [ data ] = ssdb_import_habit( habit_id, orientation, allow_nodata, ...
                                   size_range, size_type, freq_range, temp_range )
%
if nargin < 3  ||  isempty(allow_nodata)
  allow_nodata = false;
end
if nargin < 4  ||  isempty(size_range)
  size_range = [ 0 Inf ];
end
if nargin < 5  ||  isempty(size_type)
  size_type = 'dveq';
end
if nargin < 6  ||  isempty(freq_range)
  freq_range = [ ];   % To trigger default set in *ssdb_data_import*
end
if nargin < 7  ||  isempty(temp_range)
  temp_range = [ ];   % To trigger default set in *ssdb_data_import*
end
%
if length(habit_id) ~= 1
  error( 'Length of *habit* must be one.' );
end
%
if ~ischar( orientation )
  error( '*orientation* must be a string.' );
end
%
size_type = lower( size_type );
%
if ~isnumeric(size_range)  ||  length(size_range) ~= 2
  error( '*size_range* must be a numeric vector of length 2.' );
end
if ~( strcmp(size_type,'dveq') || strcmp(size_type,'dmax') || ...
      strcmp(size_type,'mass') )
  error( ['Valid options for *size_type* are ''dveq'', ''dmax'', and ''mass''. ' ...
          'You selected %s.'], size_type );
end
%
% freq_range and temp_range are checked in ssdb_data_import


% Make an inventory of existing data
%
[orientation_folders, orientations, tilt_angles] = ...
                           ssdb_particle_folders( ssdb_habits(habit_id) );
%
if isempty(orientation_folders)
  error( 'No orientation folders could be located for specified habit (%d).', habit_id );
end

% Filter out selected orientation
%
ind = find(strcmp(orientations, orientation));
%
assert( (length(ind)<2), ...
        ['More than one matching orientation found. Should not happen.\n' ... 
         'Complain to interface developers!'] )
if isempty(ind)
  fprintf( 'The following orientations are at hand:\n' );
  os = cell(0,0);
  for k = 1 : length(orientations)
    os{end+1} = orientations{k};
  end
  os = unique( os );
  for i = 1 : length(os)
    fprintf( '   %s\n', os{i} );
  end
  error( 'No data could be located for specified orientation (%s).', orientation );
end
%
% Crop data according to orientation
orientation_folder = orientation_folders{ind};
tilt_angle = tilt_angles( ind );

%   
[data_files,dmax,dveq,mass] = ssdb_data_files( orientation_folder );
% Crop data according to "size"
%
if strcmp(size_type,'dveq')
  size = dveq;
elseif strcmp(size_type,'dmax')
  size = dmax;
elseif strcmp(size_type,'mass')
  size = mass;
end
ind =  size > size_range(1)  &  size <= size_range(2);
data_files = data_files(ind);
dmax = dmax(ind);
for n = length(data_files):-1:1
  [SSD_temp, freq, temp, nodata] =  ssdb_read_ncfile( data_files{n}, freq_range, temp_range);
  if isempty(SSD_temp)
      continue
  end
  data(n).SSD = SSD_temp;
  ii = find(~nodata(:),1);
  data(n).dveq = SSD_temp(ii).ShapeData.diameter_vol_eq;
  data(n).dmax = SSD_temp(ii).ShapeData.diameter_max;
  data(n).mass = SSD_temp(ii).ShapeData.mass;
  data(n).freq = freq;
  data(n).temp = temp;
  data(n).nodata = nodata;
  data(n).orientation = orientation;
  data(n).tilt_angle = tilt_angle;
  if ~allow_nodata  &&  any(nodata(:))
    fprintf( ['At least one frequency and temperature combination lacks data ' ...
            'for id/orientation %d/%s and Dmax = %.2f um.\n'], habit_id, ...
           orientation, dmax(n)*1e6 );
    error( ['Lacking data for a frequency-temperature combination, see above ' ...
          'for details. To allow lacking data, set allow_nodata=true.'] );
  end
end


