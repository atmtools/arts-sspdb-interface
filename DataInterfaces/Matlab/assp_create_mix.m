% ASSP_CREATE_MIX   Sums up scattering data 
%
%   This function can be used to create scattering data for a habit mix. The
%   data are summed up according to provided weights. 
%
%   The input arguments *S0*, *M0* and *W* must have the same size. Each row is
%   assumed to correspond to a habit, and each column to a size. That is, the
%   scattering data for habit i and size j are S0(i,j) and all *S0* in a column
%   should have common size.
%   
%   Each S0(i,j) is given weight W(i,j). Accordingly, each column of *W* should
%   sum up to 1. Diameters and masses in *M0* are given the same weights, to
%   create the matching data in final *M*.
%
%   No interpolations are performed, and all data for each size must have
%   common grids and have the same *ptype*. But note that grids and *ptype*
%   are allowed to differ between sizes.
%
% FORMAT   [S,M] = assp_create_mix( S0, M0, W )
%
% OUT   S    Scattering data for habit mix
%       M    Meta data for habit
%  IN   S0   Scattering data for original particles 
%       M0   Meta data for original particles 
%       W    Weight matrix


% 2016-12-30 Patrick Eriksson


function [S,M] = assp_create_mix( S0, M0, W )


% Basic input checks
%
nhabits = size(S0,1);
%
if nhabits == 1
  error( 'Your data just contain one habit. Why use this function?' );
end
if nhabits ~= size(M0,1)  |  size(S0,2) ~= size(M0,2) 
  error( '*S* and *M* must have the same size.' );
end
if nhabits ~= size(W,1)  |  size(S0,2) ~= size(W,2) 
  error( '*S* and *W* must have the same size.' );
end
if any( abs(sum(W)-1) > 0.001 )
  error( 'The sum over each column in *W* shall be 1.' );
end


% Check and harmonise content
%
for j = 1 : size(S0,2)

  ptype = zeros( nhabits, 1 );

  for i = 1 : nhabits
    
    if i > 1
      if ~isequal( S0(1,j).f_grid, S0(i,j).f_grid )
        error( 'The frequency grid of all elements of *S0* must be the same.' );
      end
      if ~isequal( S0(1,j).T_grid, S0(i,j).T_grid )
        error( 'The temperature grid of all elements of *S0* must be the same.' );
      end
      if ~isequal( S0(1,j).za_grid, S0(i,j).za_grid )
        error( 'The zenith angle grid of all elements of *S0* must be the same.' );
      end
      if ~isequal( S0(1,j).aa_grid, S0(i,j).aa_grid )
        error( 'The azimuth angle grid of all elements of *S0* must be the same.' );
      end 
    end

    if strcmp( S0(i,j).ptype, 'totally_random' )
      ptype(i) = 1;
    elseif strcmp( S0(i,j).ptype, 'azimuthally_random' )
      ptype(i) = 2;
    else
      error( 'Unknow *ptype* found in S0(%d,%d): %s', i, j, S0(i,j).ptype );
    end
  end 
  
  % Do we have a mix of ptypes?
  if length( unique( ptype ) ) > 1
    error( 'Mixed ptypes are not yet handled.' );
  end
end



% Create data for weighted mx
%
% Init output argument
S = S0(1,:);
M = M0(1,:);
%
for j = 1 : size(S0,2)

  S(j).description     = 'A habit mix';
  S(j).abs_vec_data(:) = 0;
  S(j).ext_mat_data(:) = 0;
  S(j).pha_mat_data(:) = 0;
  %
  M(j).description     = 'A habit mix';
  M(j).source          = '';
  M(j).refr_index      = '';
  M(j).mass                            = 0;
  M(j).diameter_max                    = 0;
  M(j).diameter_volume_equ             = 0;
  M(j).diameter_area_equ_aerodynamical = 0;
  
  for i = 1 : nhabits

    S(j).abs_vec_data = S(j).abs_vec_data + W(i,j) * S0(i,j).abs_vec_data;
    S(j).ext_mat_data = S(j).ext_mat_data + W(i,j) * S0(i,j).ext_mat_data;
    S(j).pha_mat_data = S(j).pha_mat_data + W(i,j) * S0(i,j).pha_mat_data;
    %
    M(j).mass                            = M(j).mass                            + ...
                               W(i,j) * M0(i,j).mass;
    M(j).diameter_max                    = M(j).diameter_max                    + ...
                               W(i,j) * M0(i,j).diameter_max;
    M(j).diameter_volume_equ             = M(j).diameter_volume_equ             + ...
                               W(i,j) * M0(i,j).diameter_volume_equ;
    M(j).diameter_area_equ_aerodynamical = M(j).diameter_area_equ_aerodynamical + ...
                               W(i,j) * M0(i,j).diameter_area_equ_aerodynamical;
  end
end
