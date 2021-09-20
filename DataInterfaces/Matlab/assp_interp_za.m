% ASSP_INTERP_ZA   Zenith angle interpolation of ARTS single scattering properties
%
%   Input data must cover zenith angles of 0 and 180 deg. The new grid can
%   either be set by the user, or interpolation to a common grid is made.
%
% FORMAT   S = assp_interp_za(S,[new_za_grid,interpm])
%
% OUT   S              As input but data interpolated to *new_za_grid*.
%  IN   S              Input ASSP structure(s). I.e. *S* can hold multiple structures.
% OPT   new_za_grid    New zenith angle grid. Default is to create *new_za_grid* 
%                      as the union of all S.za_grid. An interpolation is always
%                      applied and time is wasted if the data are already on a
%                      common zenith angle grid.
%       interpm        Interpolation method. Default is 'linear'.

% 2017-03-14 Patrick Eriksson


function S = assp_interp_za( S, new_za_grid, interpm )
%
if nargin < 2 
  new_za_grid  = [];   
end
if nargin < 3  |  isempty(interpm)
  interpm = 'linear';   
end


if isempty( new_za_grid )
  for i = 1 : size(S,1)
    for j = 1 : size(S,2)
      assert( S(i,j).za_grid(1) == 0 );
      assert( S(i,j).za_grid(end) == 180 );
      new_za_grid = union( new_za_grid, S(i,j).za_grid );
    end
  end
end

nza = length( new_za_grid );


for i = 1 : size(S,1)
  for j = 1 : size(S,2)

    % Copy existing data 
    %
    old_za_grid = S(i,j).za_grid;
    %
    p = S(i,j).pha_mat_data;

    % Reallocate 
    %
    S(i,j).za_grid = new_za_grid;
    %
    s    = size( p );
    s(3) = nza;
    S(i,j).pha_mat_data = zeros( s );
    
    for f = 1 : length( S(i,j).f_grid )

      for t = 1 : length( S(i,j).T_grid )
          
        D = squeeze( p(f,t,:,:,:,:,:) );

        S(i,j).pha_mat_data(f,t,:,:,:,:,:) = interp1( old_za_grid, D, ...
                                                      new_za_grid, interpm );
      end
    end
  end
end
