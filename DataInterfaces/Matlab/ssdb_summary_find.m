% SSDB_SUMMRY_FIND   Locates the summary file inside a habit folder
%
%    The summary file must be placed in *habitfolder*. Sub-folders are not
%    searched for any summary file.
%
% FORMAT   summary_file = ssdb_summary_find( habit_folder )
%
% OUT   summary_file  Full path to the summary file.
%  IN   habit_folder  Full path to a habit folder.

% 2016-10-29 Patrick Eriksson
% 2017-04-07 Robin Ekelund: l.21. Changed dataSummary name.


function summary_file = ssdb_summary_find( habit_folder )


content = dir( habit_folder );


i = find( strncmp( 'DataSummary.', {content.name}, 12 ) );

if length(i) == 0
  error( 'No data summary file could be found in folder %s', habit_folder );
end
if length(i) > 1
  error( 'Multiple data summary files were found in folder %s', habit_folder );
end


summary_file = fullfile( habit_folder, content(i).name );


