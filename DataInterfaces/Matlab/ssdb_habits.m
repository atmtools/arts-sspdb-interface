% SSDB_HABITS   Keeps track on basic habit information.
%
%   The function keeps track of defined habit id-s, and the folder for each
%   habit. This information is stored by two internal (persistent)
%   information variables. To create these information variables, use the
%   first format version below, that is call the function with the path to
%   the database top folder. 
%
%   To obtain the folder for a habit, use the second format version.
%   If the specified habit_id is not defined, there is an error. If
%   the function is called without an output argument, information on
%   *folderpath* is printed to the screen.
%
%   Setting habit id to -1 has a special meaning. In this case, the
%   the function returns all defined habit id-s as a vector.
%   
% FORMAT   ssdb_habits( topfolder )
%
% IN   topfolder   Full path to top folder of the database
%
%   or
%
% FORMAT   folder = ssdb_habits( habit_id )
%
% OUT   folder     Full path to folder for specified habit
% IN    habit_id   Habit id

% 2016-10-28 Patrick Eriksson
% 2016-11-17 Robin Ekelund: Bug fix l.126 (error msg)
% 2017-09-15 Robin Ekelund: Changed so that only one of the Ice, Mixed and
%                           Liquid folders are required in the SSD folder.
% 2017-09-18 Robin Ekelund: Fixed error message.
% 2017-10-23 Robin Ekelund: Changed. Changed according to the new folder structure.
function folder = ssdb_habits( habit_id );
%
if nargin ~= 1
  error( 'The function requires/handles only a single input argument.' );
end


persistent habits folders


%- Fill variables linking habits with a folder
%
if ischar( habit_id )

  if nargout
    error( 'No output argument is provided when *topfolder* is input.' );
  end
    
  % Scan the database folder structures (local sub-function)
  %
  [hs,folders] = scan_topfolder( habit_id );
  %
  if ~length(hs)
    error( 'Not a single folder was found. Something must be wrong.' );
  end

  % Repack habit id-s
  %
  n      = length( hs );
  habits = zeros( n, 1, 'int16' );
  %
  for i = 1 : n
    habits(i) = hs{i};
  end
    
  % Check that no id numbers is duplicated 
  %
  hs = unique ( habits );
  %
  if length(hs) < n
    for i = 1 : length(hs)
      ii = find( hs(i) == habits );
      if length( ii ) > 1
        fprintf( 'Habit id %d is used %d times.\n', hs(i), length(ii) );
      end
    end
    error( 'One or several habit numbers are duplicated. See above.' );
  end

  
%- Extract folder for specified habit
%
elseif isnumeric( habit_id )  &  length(habit_id) == 1

  if isempty( habits )
    error( ['The function has not yet been initialised. Do this by calling ' ...
            'the function with the top folder path as input.'] );
  end
    
  if habit_id == -1

    folder = habits;

  else
    
    i = find( habit_id == habits );
   
    if isempty(i)
      error( 'There is no habit with id = %d.', habit_id );
    end
     
    folder = folders{i};
   
    if nargout == 0
      fprintf( 'Habit %d is found in folder:\n%s\n', habit_id, folder );
    end
  end
  
  
else
  error( 'Unvalid input argument.' );
end

return


%--------------------------------------------------------------------------------

function [hs,folders] = scan_topfolder( topfolder )

  main_folders = { 'TotallyRandom', 'AzimuthallyRandom'};

  % Check that topfolder contains any of the "main folders"
  main_folders_path = fullfile( topfolder, main_folders );
  bool_test = cellfun(@(x)exist(x, 'dir' ) == 7, main_folders_path);
  if all(~bool_test) 
    error( [ 'There should be a folder in *topfolder* named either ''TotallyRandom'' or ''AzimuthallyRandom''.'] );
  end


  % Scan the main folders
  %
  hs      = [];
  folders = [];
  %
  for i = 1 : length( main_folders )
    [hs,folders] = scan_folder( fullfile(topfolder,main_folders{i}), hs, folders );
  end

return


function [hs,folders] = scan_folder( datafolder, hs, folders )

  content = dir( datafolder );
  
  % Use report.pdf as indicator if we are in a habit folder, or higher up in
  % the folder tree
  
  if( any( strcmp( {content.name}, 'report.pdf' ) ) )

    % We have reached a habit folder
    %
    summary_file = ssdb_summary_find( datafolder );
    %
    S         = ssdb_summary_read( summary_file );
    hs{end+1} = S.HABIT_IDENT;
    %
    folders{end+1} = datafolder;
          
  else
    
    % We have to continue downwards
    for i = 1 : length(content)
      % Any possible file and . and .. should be ignored
      if content(i).isdir  &&  content(i).name(1) ~= '.'
        [hs,folders] = scan_folder( fullfile(datafolder,content(i).name), hs, folders );
      end
    end  
  end 

return
