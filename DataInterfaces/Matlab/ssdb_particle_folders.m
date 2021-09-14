% SSDB_PARTICLE_FOLDERS   Finds folders for individual particles
%
%   Explores the particle data of a given habit folder.
%
% FORMAT [orientation_folder, orientations, tilt_angles] = ssdb_particle_folders( habit_folder )
%
% OUT    orientation_folder  Full path to each orientation folder found. An array of strings
%        orientations  Orientation type cell
%        tilt_angles   tilt angle vector, according to folder name
%  IN    habit_folder  Full path to a habit folder.
%
% 2016-10-30 Patrick Eriksson
% 2016-11-03 Robin Ekelund: Added dveq as output. Textscan now used to
%                               extract file-name data in a more flexible
%                               manner.
% 2017-10-20 Robin Ekelund: Changed. Now looks for orientation folders.
% 2017-10-30 Robin Ekelund: Fix. Syntax problem for ARO relevant code.
% 2017-10-31 Robin Ekelund: Cleaned up a bit.
%
function [orientation_folder, orientations, tilt_angles] = ssdb_particle_folders( habit_folder )


content = dir( habit_folder );

ihit = find([content.isdir] & ~ismember({content.name}, {'.' '..'}) );

if isempty( ihit )
  error( ['Not a single orientation folder was found. Have you really provided a ' ...
          'habit folder?'] );
end

n    = length(ihit);
content = content(ihit);
orientation_folder = fullfile( habit_folder, { content.name } );
orientations = cell(n, 1);
tilt_angles = cell(n, 1);

for i = 1 : n
  folder_name = content(i).name;
  if strcmp(folder_name, 'TotallyRandom')
    orientations{i} = 'totally_random';
  elseif strncmp(folder_name, 'AzimuthallyRandom', 17)
    tmp = textscan(folder_name, 'AzimuthallyRandom_beta%fdeg');
    tilt_angles{i} = tmp{1};
    orientations{i} = strrep(folder_name, 'AzimuthallyRandom', 'azimuthally_random');
  else
    error('Folder %s not recognized as an orientation folder.', folder_name);
  end
end
end
