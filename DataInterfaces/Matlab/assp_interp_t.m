% ASSP_INTERP_T   Temperature interpolation of ARTS single scattering properties
%
%   The function allows you to select interpolation method and if extrapolation
%   is allowed or not. The input data are allowed to have have missing data, as
%   long as missing data are set to Nan in the ext_mat_data field. It is
%   allowed that the temperature corresponding to missing data differs
%   between the frequencies. The number of missing data points is
%   considered when checking if *min_tpoints* is met.
%
% FORMAT   S = assp_interp_t(S,[new_t_grid,interpm,allow_extrap,min_tpoints])
%
% OUT   S              As input but data interpolated to *new_t_grid*.
%  IN   S              Input ASSP structure(s). I.e. *S* can hold multiple structures.
% OPT   new_t_grid     New temperature grid. Default is to create *new_t_grid*
%                      as the union of all S.T_grid. An interpolation is always
%                      applied and time is wasted if the data are already on a
%                      common temperature grid.
%       interpm        Interpolation method. Default is 'pchip'.
%       allow_extrap   Flag to allow extrapolation. Default is false.
%       min_tpoints    Minimum number of temperatures with data. Default is 3.

% 2016-11-03 Patrick Eriksson


function S = assp_interp_t( S, new_t_grid, interpm, allow_extrap, min_tpoints )
%
if nargin < 2 
  new_t_grid  = [];   
end
if nargin < 3  |  isempty(interpm)
  interpm = 'pchip';   
end
if nargin < 4  |  isempty(allow_extrap)
  allow_extrap = false;   
end
if nargin < 5  |  isempty(min_tpoints)
  min_tpoints = 3;   
end


if isempty( new_t_grid )
  for i = 1 : size(S,1)
    for j = 1 : size(S,2)
      new_t_grid = union( new_t_grid, S(i,j).T_grid );
    end
  end
end


nt = length( new_t_grid );


for i = 1 : size(S,1)
  for j = 1 : size(S,2)

    if ~allow_extrap
      if any( new_t_grid < S(i,j).T_grid(1) )
        error( ['Lowest temperature in S(%d,%d) is %.1f and an extraplotation is ' ...
                'needed, while *allow_extrap* is set to be false.'], ...
               i, j, S(i,j).T_grid(1) );
      end
      if any( new_t_grid > S(i,j).T_grid(end) )
        error( ['Highest temperature in S(%d,%d) is %.1f and an extraplotation is ' ...
                'needed, while *allow_extrap* is set to be false.'], ...
               i, j, S(i,j).T_grid(end) );
      end
    end
    if length( S(i,j).T_grid ) < min_tpoints
      error( ['You want at least %d temperature points to allow an interpolation ' ...
              'but length of S(%d,%d).T_grid is just %d'], min_tpoints, ...
             i, j, length(S(i,j).T_grid) );
    end
    
    % Copy existing data 
    %
    old_t_grid = S(i,j).T_grid;
    %
    e = S(i,j).ext_mat_data;
    a = S(i,j).abs_vec_data;
    p = S(i,j).pha_mat_data;

    % Reallocate 
    %
    S(i,j).T_grid = new_t_grid;
    %
    s    = size( e );
    s(2) = nt;
    S(i,j).ext_mat_data = zeros( s );
    s    = size( a);
    s(2) = nt;
    S(i,j).abs_vec_data = zeros( s );
    s    = size( p );
    s(2) = nt;
    S(i,j).pha_mat_data = zeros( s );
    
    for f = 1 : length( S(i,j).f_grid )

      if strcmp( S(i,j).ptype, 'totally_random' ) 
        D = e(f,:)';      
      elseif strcmp( S(i,j).ptype, 'azimuthally_random' ) 
        D = squeeze( e(f,:,:,:,:) );
      else
        error( 'Unknow *ptype* found in S(%d,%d): %s', i, j, S(i,j).ptype );
      end  

      % Non-NaN temperature index     
      %
      iok = find( ~isnan( D(:,1,1,1) ) );
      %
      if length(iok) < min_tpoints
      keyboard
        error( ['You have demanded at least %d temperature points, but for S(%d,%d) ' ...
                'and frequency %.1f GHz there are only %d temperatures.'], ...
               min_tpoints, i, j, S(i,j).f_grid(f)/1e9, length(iok) );
      end
      
      if allow_extrap
        S(i,j).ext_mat_data(f,:,:,:,:) = interp1( old_t_grid(iok), D(iok,:,:,:), ...
                                                new_t_grid, interpm, 'extrap' );
      else
        S(i,j).ext_mat_data(f,:,:,:,:) = interp1( old_t_grid(iok), D(iok,:,:,:), ...
                                                new_t_grid, interpm );
      end

      if strcmp( S(i,j).ptype, 'totally_random' ) 
        D = a(f,:)';      
      elseif strcmp( S(i,j).ptype, 'azimuthally_random' ) 
        D = squeeze( a(f,:,:,:,:) );
      else
        error( 'Unknow *ptype* found in S(%d,%d): %s', i, j, S(i,j).ptype );
      end  
      if allow_extrap
        S(i,j).abs_vec_data(f,:,:,:,:) = interp1( old_t_grid(iok), D(iok,:,:,:), ...
                                                new_t_grid, interpm, 'extrap' );
      else
        S(i,j).abs_vec_data(f,:,:,:,:) = interp1( old_t_grid(iok), D(iok,:,:,:), ...
                                                new_t_grid, interpm );
      end
      D = squeeze( p(f,:,:,:,:,:,:) );
      if allow_extrap
        S(i,j).pha_mat_data(f,:,:,:,:,:,:) = interp1( old_t_grid(iok), ...
                                       D(iok,:,:,:,:,:), new_t_grid, interpm, 'extrap' );
      else
        S(i,j).pha_mat_data(f,:,:,:,:,:,:) = interp1( old_t_grid(iok), ...
                                       D(iok,:,:,:,:,:), new_t_grid, interpm );
      end
    end
  end
end
