% ASSP_INTERP_F   Frequency interpolation of ARTS single scattering properties
%
%   The function allows you to select interpolation method and if extrapolation
%   is allowed or not. The input data are NOT allowed to have have missing data.
%
% FORMAT   S = assp_interp_f(S,[new_f_grid,interpm,allow_extrap])
%
% OUT   S              As input but data interpolated to *new_t_grid*.
%  IN   S              Input ASSP structure(s). I.e. *S* can hold multiple structures.
% OPT   new_f_grid     New frequency grid. Default is to create *new_f_grid*
%                      as the union of all S.f_grid. An interpolation is always
%                      applied and time is wasted if the data are already on a
%                      common frequency grid.
%       interpm        Interpolation method. Default is 'pchip'.
%       allow_extrap   Flag to allow extrapolation. Default is false.

% 2016-12-28 Patrick Eriksson
% 2017-10-08 Robin Ekelund: Cleaned up header.


function S = assp_interp_f( S, new_f_grid, interpm, allow_extrap )
%
if nargin < 2
  new_f_grid  = [];
end
if nargin < 3  |  isempty(interpm)
  interpm = 'pchip';
end
if nargin < 4  |  isempty(allow_extrap)
  allow_extrap = false;
end


if isempty( new_f_grid )
  for i = 1 : size(S,1)
    for j = 1 : size(S,2)
      new_f_grid = union( new_f_grid, S(i,j).f_grid );
    end
  end
end


for i = 1 : size(S,1)
  for j = 1 : size(S,2)

    if ~allow_extrap
      if any( new_f_grid < S(i,j).f_grid(1) )
        error( ['Lowest frequency in S(%d,%d) is %.1f GHz and an extrapolation ' ...
                'is needed, while *allow_extrap* is set to be false.'], ...
               i, j, S(i,j).f_grid(1)/1e9 );
      end
      if any( new_f_grid > S(i,j).f_grid(end) )
        error( ['Highest frequency in S(%d,%d) is %.1f GHz and an extrapolation ' ...
                'is needed, while *allow_extrap* is set to be false.'], ...
               i, j, S(i,j).f_grid(end)/1e9 );
      end
    end

    % Replace f_grid
    %
    old_f_grid    = S(i,j).f_grid;
    S(i,j).f_grid = new_f_grid;
    
    % Interpolate
    %
    S(i,j).abs_vec_data = interp1( old_f_grid, S(i,j).abs_vec_data, ...
                                               new_f_grid, interpm, 'extrap' );
    S(i,j).ext_mat_data = interp1( old_f_grid, S(i,j).ext_mat_data, ...
                                               new_f_grid, interpm, 'extrap' ); 
    S(i,j).pha_mat_data = interp1( old_f_grid, S(i,j).pha_mat_data, ...
                                               new_f_grid, interpm, 'extrap' );
  end
end
