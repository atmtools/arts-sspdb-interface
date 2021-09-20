% SSDB_SUMMARY   Habit summary data
%
%     The functions obtains the summary data for given habit id. 
%
%     If the function is called without output argument, the summary data
%     are displayed on screen.
%   
% FORMAT   S = ssdb_summary( habit_id )
%
% OUT   S          Structure with summary data
% IN    habid_id   Habit id number

% 2016-10-29 Patrick Eriksson


function S = ssdb_summary( habit_id )

habit_folder = ssdb_habits( habit_id );

summary_file = ssdb_summary_find( habit_folder );

S = ssdb_summary_read( summary_file );


if ~nargout
   
  fields = fieldnames( S );
  
  l2col = 15;
  
  fprintf( '\n' );
  
  for i = 1 : length(fields)
    
    for j = 1 : (l2col - length(fields{i}))
      fprintf( ' ' );
    end
    if strcmp( fields{i}, 'HABIT_IDENT' )
      fprintf( '%s: %d\n', fields{i}, S.(fields{i}) ); 
    else
      fprintf( '%s: %s\n', fields{i}, S.(fields{i}) ); 
    end
  end
  fprintf( '\n' );
end
