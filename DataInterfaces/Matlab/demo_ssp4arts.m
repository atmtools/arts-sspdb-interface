% NOTE: 1) Requires atmlab.
%       2) Currently produces one scat and meta data file with all scattering
%          elements (packed as one scattering species) included. In ARTS, you can
%          use ReadXML (and Append in case you want more than one scattering
%          species) to import these data.


%%%%%%%%%%%%%%%%%%%%
% Input. Modify according to your needs.
%%%%%%%%%%%%%%%%%%%%

% location of SSDB in your system
SSDBpath = '/your/database/location/SSD/';

% location (and name) of ARTS scat_data and scat_meta output
%   NOTE: as is, this extracts data for a single scattering species (with
%         multiple scattering elements) and stores them in *scat_data* and
%         *scat_meta* equivalent containers (ie if you need further scattering
%         species, you have to prepared (append) them inside ARTS.
outfolder = pwd;
sdname = 'scat_data.xml';
smname = 'scat_meta.xml';

% properties of particles to extract SSP for
%   use ssdb_display to see what data (habits, orientations, sizes,
%   frequencies, temperatures) are available.
%   see assp_import_ssdb for further details on the parameters

habID = 1; % habit ID
orient = 'totally_random'; % orientation

minD = 0.; % minimum size to extract
maxD = Inf; % maximum size to extract
sizeparam = 'dmax'; % size parameter to apply extraction limits on.
                    % available (param [unit]): 'dmax' [m], 'dveq' [m], 'mass' [kg]

fmin = 0.; % minimum frequency to extract
fmax = Inf; % maximum freuency to extract

tmin = 0.; % minimum temperature to extract
tmax = Inf; % maximum temperature to extract



%%%%%%%%%%%%%%%%%%%%
% Data extraction. (Likely) no need to modify.
%%%%%%%%%%%%%%%%%%%%

% Init database
ssdb_init( SSDBpath );

% Import data, with some cropping in size and freq
[S,M] = assp_import_ssdb( habID, orient, false, ...
                          [minD maxD], sizeparam, ...
                          [fmin fmax], ...
                          [tmin tmax] );

% Some processing could be necessary. For example, to ensure that data are
% ordered in Dmax: 
%
[dmax,ind] = unique( [ M.diameter_max ] );
%
S = S(ind);
M = M(ind);

% Convert S and M to ARTS internal format (assuming this is scat species #1)
for i = 1 : length(S)
  scat_data{1}{i} = S(i);
  scat_meta{1}{i} = M(i);
end

% Create files
xmlStore( fullfile( outfolder, sdname ), scat_data, ...
            'ArrayOfArrayOfSingleScatteringData', 'binary' );
xmlStore( fullfile( outfolder, smname ), scat_meta, ...
            'ArrayOfArrayOfScatteringMetaData', 'binary' );
