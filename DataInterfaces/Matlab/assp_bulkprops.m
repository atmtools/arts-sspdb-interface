% ASSP_BULKPROPS   Calculates bulk properties for given PSD
%
%   Particle properties are specified by vectors of scattering data and
%   matching meta data. These data can be provided unsorted in size. 
%
%   The quantity to consider as size is selected by *size_unit*. The options
%   are: 
%     'dveq' : Size is taken from M.diameter_volume_equ
%     'dmax' : Size is taken from M.diameter_max
%     'area' : Size is taken from M.diameter_area_equ_aerodynamical
%     'mass' : Size is taken from M.mass
%
%   The particle size distribution (PSD) is specified by a function handle. The
%   function selected shall take size as first input, and can have an arbitrary
%   number of additional input arguments. The PSD arguments are given after
%   *psd*.
%   For example, if an exponential PSD exists having the format
%      n = exp_psd(x,n0,la)
%   this function could be called as
%      B = assp_bulkprops(S,M,'dveq',@exp_psd,n0,la);
%   Another example, with PSD defined as an anonymous function 
%      f = @(x,n0,la)n0*exp(-la*x);
%      B = assp_bulkprops(S,M,'dmax',@f,1e6,2000);
%
%   Of course, the PSD function and arguments selected must be consistent
%   with *size_unit*. This consistency is up to the user, it can not be
%   checked by the function.
%
%   The derived bulk properties are returned using the SingleScatteringData
%   format. The only difference to monodisperse single scattering data is that
%   the extinction, absorption and scattering data now are not in terms of
%   cross-sections anymore, but are integrated bulk properties in terms of
%   ext/abs/sca coefficient (i.e. extinction per meter for B.ext_mat_data).
%
%   No interpolations are performed, and all data for each individual size must
%   have common grids and the same *ptype*.
%   (Note: currently, only *ptype* 'totally_random' allowed).
%
% FORMAT B = assp_bulkprops( S, M, size_unit, psd, psd_arg1, psd_arg2, ... )
%
% OUT   B          Bulk properties, in the format of SingleScatteringData.
% IN    S          Vector of SingleScatteringData.
%       M          Vector of ScatteringMetaData
%       size_unit  See above.
%       psd        See above.
%       psd_arg1   PSD input argument 1
%       ...

% 2017-10-21  Patrick Eriksson
% 2017-10-27  Robin Ekelund: Bug fixes.


function B = assp_bulkprops( S, M, size_unit, psd, varargin )
%
if ~strcmp( S(1).ptype, 'totally_random' )
  error( 'This function handles so far only totally random orientation data.' );
end
if min(size(S)) > 1
  error( 'S must be a vector of scattering data.' );
end
if min(size(M)) > 1
  error( 'M must be a vector of scattering meta data.' );
end
if length(S) ~= length(M)
  error( 'S and M must have the same length.' );
end
if length(S) <= 1
  error( 'Length of S must be > 1.' );
end


%- Size grid
%
switch lower(size_unit)
  case 'dveq'
    x = [M.diameter_volume_equ];
  case 'dmax'
    x = [M.diameter_max];
  case 'area'
    x = [M.diameter_area_equ_aerodynamical];
  case 'mass'
    x = [M.mass];
  otherwise
    error( ['Valid options for *size_unit* are: ''dveq'', ''dmax'', ''area'' ' ...
            'and ''mass''.'] );
end
%
if any( isnan( x ) )
  error( 'There is at least one Nan in data matching size_unit = %s.', size_unit );
end
%
% We want x to be a column vector
if size(x,1) == 1
  x = x';
end


%- Sort in size (otherwise trapz will fail)
%
[x,i] = sort( x ); 
S     = S(i);


%- Loop sizes and compile data
%
ns  = length( S );
nf  = length( S(1).f_grid );
nt  = length( S(1).T_grid );
nza = length( S(1).za_grid );
%
abs_vec = zeros( ns, nf, nt );
ext_mat = zeros( ns, nf, nt );
pha_mat = zeros( ns, nf, nt, nza, 1, 1, 1, 6 );
%
for i = 1 : ns
  abs_vec(i,:,:) = S(i).abs_vec_data;
  ext_mat(i,:,:) = S(i).ext_mat_data;
  pha_mat(i,:,:,:,:,:,:) = S(i).pha_mat_data;
end
%
if any( isnan( abs_vec(:) ) )
  error( 'There is at least one Nan in S.abs_vec_data.' );
end
if any( isnan( ext_mat(:) ) )
  error( 'There is at least one Nan in S.ext_mat_data.' );
end
if any( isnan( pha_mat(:) ) )
  error( 'There is at least one Nan in S.pha_mat_data.' );
end


%- Get PSD
%
n = psd( x, varargin{:} );


%- Init output argument
%
B = S(1);


%- Calculate bulk properties
%
D              = trapz( x, abs_vec .* repmat( n, [1 nf nt] ) );
B.abs_vec_data = reshape( D, [nf nt] );
D              = trapz( x, ext_mat .* repmat( n, [1 nf nt] ) );
B.ext_mat_data = reshape( D, [nf nt] );
D              = trapz( x, pha_mat .* repmat( n, [1 nf nt nza 1 1 1 6] ) );
B.pha_mat_data = reshape( D, [nf nt nza 1 1 1 6] );

