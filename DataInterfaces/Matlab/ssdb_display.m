% SSDB_DISPLAY   Explore the database content
%
%     The function summarises the database content on the screen. Different
%     levels of the database are displayed depending on the number of input
%     arguments. Three format versions ar at hand:
%   
% To get an overview of habits and orientations at hand:  
%
% FORMAT   ssdb_display
%   
%
% To list the particle sizes at hand for a habit and orientation combination:
%
% FORMAT   ssdb_display( habit_id, orientation )
%   
%
% To list the frequencies and temperature combinations at hand for a habit,
% orientation and size combination:
%
% FORMAT   ssdb_display( habit_id, orientation, dveq )

% 2016-10-29 Patrick Eriksson
% 2016-11-21 Robin Ekelund : l.109
% 2017-05-03 Robin Ekelund : l.67, header added. l.117, added status
%                               markers to display output.
% 2017-10-20 Robin Ekelund : Changed in order to account for changes in
%                                other functions.
% 2017-10-23 Robin Ekelund : Fix. Did not work for 2 variable input.
% 2017-10-30 Robin Ekelund : Change in orientation output.
% 2017-10-30 Robin Ekelund : Fixed problem with displaying frequencies and
%                               temperatures.

function ssdb_display( habit_id, orientation, dveq )


if nargin == 0

  main_folders = { 'Ice', 'Mixed', 'Liquid' };

                
  % Get all habit id-s
  habits = ssdb_habits( -1 );

  % Obtain folders names and extract the four levels for each habit
  for i = 1 : length(habits)

    fullpath = ssdb_habits( habits(i) );

    if i == 1
        % Find length of path that shall be removed
      l = [ strfind( fullpath, main_folders{1} ), ...
            strfind( fullpath, main_folders{2} ), ...
            strfind( fullpath, main_folders{3} ) ] - 1;
    end
    
    [partialpath,level4{i}] = fileparts( fullpath(l+1:end) );
    [partialpath,level3{i}] = fileparts( partialpath );
    [partialpath,level2{i}] = fileparts( partialpath );
    [partialpath,level1{i}] = fileparts( partialpath );
    
  end


  % Print top lines
  fprintf(['\n---------------------- ARTS single scattering data ------'...
            '----------------\n\n'] );
  fprintf(['  The table below gives an overview of habits at hand. The'...
            'first three\n'] );
  fprintf(['  levels give a rough classification of each habit. The'...
            'information found\n'] );
  fprintf( '  in the two last levels is: \n\n' );
  fprintf( '        habit id: habit name                           a         / b \n' );
  fprintf( '          orientations\n' );
  fprintf( '  where mass = a * Dmax^b.\n\n' );
  fprintf( '  Habits marked with "(x)" are still work in progress. Use with caution.\n\n' );


  % Traverse levels and print info
  for i1 = 1 : length( main_folders )
    
      
    hits1 = find( strcmp( level1, main_folders{i1} ) );
    
    names2 = unique( { level2{hits1} } );
    
    if length(names2)
      fprintf( '%s\n---\n', main_folders{i1} );
    end
    
    for i2 = 1 : length( names2 )

        
      hits2 = find( strcmp( level1, main_folders{i1} )  &  ...
                    strcmp( level2, names2{i2} ) );
      names3 = unique( { level3{hits2} } );

      if length(names3)
        fprintf( ' %s\n', names2{i2} );
      end
      
      for i3 = 1 : length( names3 )
          
        hits3 = find( strcmp( level1, main_folders{i1} )  &  ...
                      strcmp( level2, names2{i2} )        &  ...
                      strcmp( level3, names3{i3} ) );
        names4 = unique( { level4{hits3} } );
        
        if length(names4)
          fprintf( '  %s\n', names3{i3} );
        end
        
        for i4 = 1 : length( names4 )

          hit4 = find( strcmp( level1, main_folders{i1} )  &  ...
                       strcmp( level2, names2{i2} )        &  ...
                       strcmp( level3, names3{i3} )        &  ...
                       strcmp( level4, names4{i4} ) );

          habit_id = habits(hit4);
          fullpath = ssdb_habits(habit_id);
          S = ssdb_summary( habit_id );
          
          % status
          if strcmp(S.STATUS, 'WORKING')
                status_symbol = '(x)';
          elseif strcmp(S.STATUS, 'COMPLETED')
                status_symbol = '';
          end
          fprintf( '   %3s %3d: %s ', status_symbol, habit_id, S.HABIT_NAME );
          for i = 1:42-length(S.HABIT_NAME)
            fprintf( ' ' );
          end
          fprintf( 'a=%.1e / b=%.2f\n', str2num(S.ALPHA), str2num(S.BETA) );

          % Orientations
          [~, orientations, tilt_angles] = ssdb_particle_folders( fullpath );
          %
          for o = 1 : length(orientations)
            fprintf( '          %s\n', orientations{o} );
          end
        end
      end
    end  
  end
  fprintf( ...
      '\n-------------------------------------------------------------------------\n\n' );
  
  

elseif nargin == 2  |  nargin == 3
    
  % We get the basic information by combining existing functions  
  %
  S = ssdb_summary( habit_id );
  %
  [orientation_folders, orientations, tilt_angles] = ssdb_particle_folders( ...
                                                            ssdb_habits(habit_id) );  
    
  % Filter out selected orientation
  ind = strcmp( orientation, orientations );

  %
  assert( (sum(ind)<2), ...
        ['More than one matching orientation found. Should not happen.\n' ... 
         'Complain to interface developers!'] )
  if ind == 0
    error( 'No data could be located for specified orientation (%s).', orientation );
  end

  % Crop data according to orientation
  orientation_folder = orientation_folders{ind};
  
  [data_files,dmax_vec,dveq_vec,mass_vec] = ssdb_data_files( orientation_folder );
  %
  clear orientations

  % Habit + orientation
  if nargin == 2
    fprintf( '\n--- %d: %s - %s ---\n\n', habit_id, S.HABIT_NAME, orientation );
    fprintf( '    a = %.1e / b = %.2f\n', str2num(S.ALPHA), str2num(S.BETA) );
    fprintf( '    N = %d\n\n', length(dmax_vec) );
    
    fprintf( '    Dveq [um]   Dmax [um]   Mass [g]\n' );
    fprintf( '    ---------   ---------   --------\n' );
  
    for i = 1 : length(dmax_vec)
      fprintf( '%10.0f%12.0f%14.2e\n', dveq_vec(i)*1e6, dmax_vec(i)*1e6, mass_vec(i)*1e3 )
    end
    fprintf( '\n' );
  
  % Habit + orientation + size
  else
 
    [dd,i] = min( abs( dveq - dveq_vec ) );

    if dd > 0.5e-6
      error( ['No particle found close to the given Dveq value. Allowed ' ...
              'difference is +-0.5 um.'] );
    end
   
    nc_folders = ssdb_read_ncfile( data_files{i}, [], [], 1);
    freq_vec = [nc_folders(:).freq];
    temp_vec = [nc_folders(:).temp];
    
    if isempty(freq_vec)
      error( ['No data files could be found for given habit and orientation ' ...
              'combination!'] );
    end

    freq = unique(freq_vec);

    fprintf( '\n--- %d: %s - %s - %.0f um ---\n\n', habit_id, S.HABIT_NAME, ...
                                                orientation, dveq_vec(i)*1e6 );
    fprintf( [' Frequency [GHz]   x [-]   Temperatures [K] \n ---------------   ',...
              '-----   ----------------\n'] );
    
    for f = 1 : length( freq )
      ii    = find( freq(f) == freq_vec );
      x     = pi * dveq_vec(i) * freq(f) / 3e8;
      temp = temp_vec(ii);
        
      fprintf( '%16.3f  %6.2f   %.2f', freq(f)/1e9, x, temp(1) );
      
      for t = 2 : length(temp)
        fprintf( '  %.2f', temp(t) );
      end
      fprintf( '\n' );
    end
    fprintf( '\n' );
  end
  
  
else
  error( 'The function shall be used with either 0, 2 or 3 input arguments.' );
end
