% ASSP_INTERP_SIZE   Size interpolation of ARTS single scattering properties
%
%   The function allows to perform an interpolation in size, with some
%   constraints on the input data. All elements in *S* must have common
%   frequency, temperature and angular grids, and have the same "ptype".
%
%   The interpolation can be done for three size variables
%     dveq : M.diameter_volume_equ is used as size parameter
%     dmax : M.diameter_max is used as size parameter
%     mass : M.mass is used as size parameter
%
%   Input data can be unsorted in size.
%
%   Data matching zero size are added automatically in the interpolation.
%   All optical properties, as well as diameters and mass, are set to 0 for
%   size zero. This zero size is not considered when checking sizes with
%   *allow_extrap* set to false. 
%
% FORMAT   [S,M] = assp_interp_size( S, M, [ new_size_grid, 
%                                    size_type, interpm, allow_extrap] )
%
% OUT   S              As input but data interpolated to *new_size_grid*.
%       M              As input but data interpolated to *new_size_grid*.
%  IN   S              Input ASSP structures. Must be of vector shape.
%       M              Meta data corresponding to *S*.
% OPT   new_size_grid  New size grid. Values and unit shall match selection
%                      of *size_type*. Default is to create *new_size_grid*
%                      as the union of all incoming sizes. An interpolation is
%                      always applied and time is wasted if the data are already
%                      on a common size grid.
%       size_type      Size variable to use, see above. Default is 'dveq'.
%       interpm        Interpolation method. Default is 'pchip'.
%       allow_extrap   Flag to allow extrapolation. Default is false.


% 2016-12-30 Patrick Eriksson


function [S,M] = assp_interp_size( S, M, new_size_grid, size_type, interpm, allow_extrap )
%
if nargin < 3
  new_size_grid  = [];   
end
if nargin < 4  |  isempty(size_type)
  size_type = 'dveq';   
end
if nargin < 5  |  isempty(interpm)
  interpm = 'pchip';   
end
if nargin < 6  |  isempty(allow_extrap)
  allow_extrap = false;   
end


% Check and compile data
%
if min(size(S)) > 1
  error( '*S* must be a column or row vector of structure (not a matrix).' );
end
%
if length(S) ~= length(M)
  error( '*S* and *M* must have the same length.' );
end
%
if ~any( strcmp( size_type, {'dveq','dmax','mass'} ) )
  error( 'Allowed options for *size_type* are: ''dveq'', ''dmax'' and ''mass''.' );
end


if isempty( new_size_grid )
  for i = 1 : length(S)
    if strcmp( size_type, 'dveq' )
      new_size_grid = union( new_size_grid, M(i).diameter_volume_equ );
    elseif strcmp( size_type, 'dmax' )
      new_size_grid = union( new_size_grid, M(i).diameter_max );
    else
      new_size_grid = union( new_size_grid, M(i).mass );
    end
  end
end



% Sort in size
%
if strcmp( size_type, 'dveq' )
  s = [M.diameter_volume_equ];
elseif strcmp( size_type, 'dmax' )
  s = [M.diameter_max];
else
  s = [M.mass];
end
%
[~,ind] = unique( s );
%
S = S(ind);
M = M(ind);


% Perform interpolation
%
old_size_grid = zeros( size(S) );
%
n         = length(S);
ext_data  = zeros( [n+1 size(S(1).ext_mat_data)] );
abs_data  = zeros( [n+1 size(S(1).abs_vec_data)] );
pha_data  = zeros( [n+1 size(S(1).pha_mat_data)] );
size_data = zeros( n+1, 3 );
area_data = zeros( n+1, 1 );
%
for i = 1 : n

  if strcmp( size_type, 'dveq' )
    old_size_grid(i) = M(i).diameter_volume_equ;
  elseif strcmp( size_type, 'dmax' )
    old_size_grid(i) = M(i).diameter_max;
  else
    old_size_grid(i) = M(i).mass;
  end
  
  if i > 1
    if ~strcmp( S(1).ptype, S(i).ptype )
      error( 'The particle type of all elements of *S* must be the same.' );
    end
    if ~isequal( S(1).f_grid, S(i).f_grid )
      error( 'The frequency grid of all elements of *S* must be the same.' );
    end
    if ~isequal( S(1).T_grid, S(i).T_grid )
      error( 'The frequency grid of all elements of *S* must be the same.' );
    end
    if ~isequal( S(1).za_grid, S(i).za_grid )
      error( 'The zenith angle grid of all elements of *S* must be the same.' );
    end
    if ~isequal( S(1).aa_grid, S(i).aa_grid )
      error( 'The azimuth angle grid of all elements of *S* must be the same.' );
    end
  end
  %
  ext_data(i+1,:,:,:,:,:)     = S(i).ext_mat_data;
  abs_data(i+1,:,:,:,:,:)     = S(i).abs_vec_data;
  pha_data(i+1,:,:,:,:,:,:,:) = S(i).pha_mat_data;
  size_data(i+1,:)            = [ M(i).mass, 
                                  M(i).diameter_max, 
                                  M(i).diameter_volume_equ ]';
  area_data(i+1,:)            = [ M(i).diameter_area_equ_aerodynamical ];
end
%
if ~allow_extrap
  if any( new_size_grid < old_size_grid(1) )
    error( ['Smallest size in S is %.2e and an extraplotation ' ...
            'is needed, while *allow_extrap* is set to be false.'], ...
           old_size_grid(1) );
  end
  if any( new_size_grid > old_size_grid(end) )
    error( ['Largest size in S is %.2e and an extraplotation ' ...
            'is needed, while *allow_extrap* is set to be false.'], ...
           old_size_grid(end) );
  end
end


% Interpolate
%
old_with0 = [ 0; vec2col(old_size_grid) ];
%
ext_data  = interp1( old_with0, ext_data, new_size_grid, interpm, 'extrap' );
abs_data  = interp1( old_with0, abs_data, new_size_grid, interpm, 'extrap' );
pha_data  = interp1( old_with0, pha_data, new_size_grid, interpm, 'extrap' );
size_data = interp1( old_with0, size_data, new_size_grid, interpm, 'extrap' );
%
% Special treatment here as M(i).diameter_area_equ_aerodynamical often is
% set to be NaN
if any(isnan( area_data ) )
  area_data = repmat( NaN, length(new_size_grid), 1 );
else
  area_data = interp1( old_with0, area_data, new_size_grid, interpm, 'extrap' );
end


% Create new S and M
%
S1 = S(1);
M1 = M(1);
%
clear S M
%
for i = 1 : length(new_size_grid)

  S(i).version       = S1(1).version;
  S(i).ptype         = S1(1).ptype;
  S(i).description   = S1(1).description;
  S(i).f_grid        = S1(1).f_grid;
  S(i).T_grid        = S1(1).T_grid;
  S(i).za_grid       = S1(1).za_grid;
  S(i).aa_grid       = S1(1).aa_grid;
  S(i).abs_vec_data  = shiftdim( abs_data(i,:,:,:,:,:), 1 );   
  S(i).ext_mat_data  = shiftdim( ext_data(i,:,:,:,:,:), 1 );   
  S(i).pha_mat_data  = shiftdim( pha_data(i,:,:,:,:,:,:,:), 1 );   


  M(i).version       = M1(1).version;
  M(i).description   = M1(1).description;
  M(i).source        = M1(1).source;
  M(i).refr_index    = M1(1).refr_index;
  M(i).mass          = size_data(i,1);
  M(i).diameter_max                    = size_data(i,2);
  M(i).diameter_volume_equ             = size_data(i,3);
  M(i).diameter_area_equ_aerodynamical = area_data(i);

  if strcmp( size_type, 'dveq' )
    M(i).diameter_volume_equ = new_size_grid(i);
  elseif strcmp( size_type, 'dmax' )
    M(i).diameter_max = new_size_grid(i);
  else
    M(i).mass = new_size_grid(i);
  end
  
end


