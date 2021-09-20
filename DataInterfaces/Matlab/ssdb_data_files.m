% SSDB_DATA_FILES   Finds data files in given folder
%
%   Explores the actual data of a given data folder.
%
% FORMAT  [data_file,dmax,dveq,mass] = ssdb_data_files( data_folder )
%
% OUT    data_file     Full file name to each data file found. An array of strings.
%        dmax          Dmax, according to file name.
%        dveq          Dveq, according to file name.
%        mass          Mass, according to file name.
%  IN    data_folder   Full path to a data folder.
%
% 2016-10-30 Patrick Eriksson
% 2016-11-03 Robin Ekelund: Textscan now used to extract file-name data in 
%                               a more flexible manner.
% 2017-10-20 Robin Ekelund: Changed to account for changes in folder
%                               stucture
%
function [data_file,dmax,dveq,mass] = ssdb_data_files( data_folder )


content = dir( data_folder );


ihit = find( strncmp( {content.name}, 'Dmax', 4 )  &  ~[content.isdir] );

if isempty( ihit )
  error( ['Not a single data file was found. Have you really provided a ' ...
          'data folder?.'] );
end

n    = length(ihit);

data_file = fullfile( data_folder, { content(ihit).name } );
dmax      = zeros( n, 1 );
dveq      = zeros( n, 1 );
mass      = zeros( n, 1 );  


for i = 1 : n
  tmp = textscan(content(ihit(i)).name, 'Dmax%fum_Dveq%fum_Mass%fkg.nc');
  dmax(i) = round(tmp{1}*1e-6, 10, 'significant');
  dveq(i) = round(tmp{2}*1e-6, 10, 'significant');
  mass(i) = round(tmp{3}, 10, 'significant');
end
end
