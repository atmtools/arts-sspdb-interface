% SSDB_SUMMARY_READ   Reads a habit summary file.
%
%     The functions performs the actual reading of an identified summary file.
%
%     The data are returned as a structure. The habit id is converted to a
%     number, while all other fields are kept as the strings extracted from
%     the file.
%   
% FORMAT   S = ssdb_summary_read( summary_file )
%
% OUT   S             Structure with summary data
% IN    summary_file  Full path to a summary file

% 2016-10-29 Patrick Eriksson
% 2017-04-07 Robin Ekelund: Fixed file not being closed


function S = ssdb_summary_read( summary_file );


if ~strcmp( summary_file(end+[-3:0]), '.txt' )
  error( 'Data summary files are expected to have extension .txt.' );
end


fid = fopen( summary_file, 'r' );
cleanupObj = onCleanup(@()fclose(fid));
%
if fid < 0
  error( 'Could not open %s for reading.', summary_file );
end


line = fgetl( fid );
%
while ischar(line)

  i = find( line == '=' );
  
  if ~length(i)
    error( ['Each line in a summary must contain at least one =. This is not ',...
            'the case for this line\n%s\nfound in file\n%s'], line, summary_file ); 
  end
    
  var  = deblank( line(1:i(1)-1) );
  info = line(i(1)+1:end);
  
  S.(var) = info;
  
  line = fgetl( fid );
end    


st = fopen( fid );
%
if st < 0
  error( 'Could not close %s after reading', summary_file );
end


if ~isfield( S, 'HABIT_IDENT' )
  error( 'No setting of HABIT_IDENT was found in %s', summary_file );    
end


S.HABIT_IDENT = str2num( S.HABIT_IDENT );


if isempty( S.HABIT_IDENT )
  error( 'Invalid setting of HABIT_IDENT in %s', summary_file );    
end
