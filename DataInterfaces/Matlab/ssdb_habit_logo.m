% SSDB_HABIT_LOGO   Displays habit logo image
%
%     The function displays the image (but can be deactivated) and returns
%     the path to the logo file.
%
% FORMAT   logo_file = ssdb_summary( habit_id [,show_image])
%
% OUT   logo_file   Full path to habit logo file
% IN    habit_id    Habit id number
% OPT   show_image  Flag to control if logo image to be displayed or not.
%                   Default is true.

% 2016-10-29 Patrick Eriksson


function logo_file = ssdb_habit_logo( habit_id, show_image )

habit_folder = ssdb_habits( habit_id );


logo_file = fullfile( habit_folder, 'shape_img.png' );

if nargin == 1  ||  show_image

  I = imread( logo_file );

  imshow( I );

end

