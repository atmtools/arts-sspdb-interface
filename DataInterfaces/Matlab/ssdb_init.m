% SSDB_INIT   Initialises the SSDB interface
%
%   The function adds the interface folder to Matlab's search path and
%   sets up the habit-folder table used by *ssdb_habits*.
%
% FORMAT   ssdb_init( topfolder )
%
% IN   topfolder   Full path to top folder of the database

% 2016-10-28 Patrick Eriksson


function ssdb_init( topfolder )

%- Add folder to Matlab's search path
%
addpath( pwd );


%- Create habit id table
%
ssdb_habits( topfolder );
