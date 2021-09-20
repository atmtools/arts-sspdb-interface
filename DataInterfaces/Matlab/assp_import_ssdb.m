% ASSP_IMPORT_SSDB   Reads and compiles data from SSDB database for one habit
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
% FORMAT  [S,M] = assp_import_ssdb( habit_id, orientation, ...
%                 [allow_nodata, size_range, size_type, freq_range, temp_range] )
% 
% OUT S             Vector of SingleScatteringData.
%     M             Vector of ScatteringMetaData.
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

% 2016-11-02 Patrick Eriksson
% 2016-11-03 Robin Ekelund: Added support for dveq (l77).
% 2017-08-29 Robin Ekelund: Added SSD_out output.
% 2017-08-31 Robin Ekelund: Bug fix.
% 2017-09-05 Robin Ekelund: Major change. Now utilizes ssdb_import_habit.

function [S,M] = assp_import_ssdb( habit_id, orientation, allow_nodata, ...
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

% import data
 data = ssdb_import_habit( habit_id, orientation, allow_nodata, ...
     size_range, size_type, freq_range, temp_range );

% Loop particle sizes and convert data
n = 0;
for i = 1 : length(data) 
  if ~isempty(data(i).SSD)
    n = n + 1;
    [S(n),M(n)] = ssdb2assp( data(i).SSD, data(i).freq, data(i).temp, data(i).nodata );
  end
end
if n == 0
  [S,M] = deal( [] );
end    
