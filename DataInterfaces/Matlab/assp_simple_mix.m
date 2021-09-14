% ASSP_SIMPLE_MIX   Combines two habits using a linear transition region
%
%    Creates habit mix scattering data. The main input is data for two
%    habit, where the first one represents small particles (the small habit)
%    and the second represents large particles.  
%
%    Dveq is used as size descriptor. The data of the mix contain all Dveq of
%    the large habit, as well as all Dveq of the small habit that are 20 um
%    smaller than the minimum Dveq of the large habit. 
%
%    A transition region is defined by setting *n_trans*. The n_trans first
%    Dveq of the large habit constitutes the transition region. The large habit
%    gets weight *w_edge* at the first of these sizes, 1-w_edge at the last
%    size, and a linear variation is applied between these two sizes (linear as
%    a function of Dveq). The weight of the small habit is 1 minus the weight
%    of the large habit. Above/below the transition region, the large/small
%    habit has weight 1.
%
% FORMAT [S,M] = assp_simple_mix( S_s, M_s, S_l, M_l, n_trans, w_edge )
%
% OUT   S        S of habit mix
%       M        M of habit mix
% IN    S_s      S of "small habit"
%       M_s      M of "small habit"
%       S_l      S of "large habit"
%       M_l      M of "large habit"
%       n_trans  Width of transition region, in terms of sizes in M_l.
%       w_edge   Weight of large/small habit at first and last transition size

% 2017-11-30   Patrick Eriksson

function [S,M] = assp_simple_mix( S_s, M_s, S_l, M_l, n_trans, w_edge )


% Dveqs of each habit
%
dveq_small = [ M_s.diameter_volume_equ ];
dveq_large = [ M_l.diameter_volume_equ ];


% Sort in Dveq
%
[dveq_small,ind] = sort( dveq_small );
S_s              = S_s( ind );
M_s              = M_s( ind );
[dveq_large,ind] = sort( dveq_large );
S_l              = S_l( ind );
M_l              = M_l( ind );


% Checks
%
if dveq_small(end) < dveq_large(n_trans)
  error('The small habit does not cover the transition region in terms of Dveq.'); 
end
    

% Determine what small particles to put in with weight 1
%
i_small = find( dveq_small < dveq_large(1) - 20e-6 );
n_small = length( i_small );
n_large = length( S_l );


% Fill output S and M, so far with weight 1 for large habit
%
S(n_small+1:n_small+n_large) = S_l;
M(n_small+1:n_small+n_large) = M_l;
%
S(1:n_small) = S_s(1:n_small);
M(1:n_small) = M_s(1:n_small);




% Loop tranistion region and create weighted data
%
for i = 1 : n_trans 
    
  % Extract small habit data just covering the Dveq of concern
  %
  ind       = max( find( dveq_small<=dveq_large(i) ) ) + [0:1];
  S_i       = S_s(ind);
  M_i       = M_s(ind);

  % Interpolate small habit to Dvewq and angles of large habit
  %
  S_i       = assp_interp_za( S_i, S_l(i).za_grid );
  [S_i,M_i] = assp_interp_size( S_i, M_i, dveq_large(i), 'dveq', 'pchip' ); 
  
  % Check that t and f grids agree
  %
  if length(S_i.f_grid) ~= length(S_s(i).f_grid) 
    error( 'Different lengths of f_grid found.' );
  end
  if length(S_i.T_grid) ~= length(S_s(i).T_grid) 
    error( 'Different lengths of T_grid found.' );
  end
  if max(abs( S_i.f_grid - S_s(i).f_grid )) > 10e3
    error( 'The two f_grids contain different frequencies.' );
  end
  if max(abs( S_i.T_grid - S_s(i).T_grid )) > 0.2
    error( 'The two T_grids contain different temperatures.' );
  end
  
  w = w_edge + (1-2*w_edge) * (dveq_large(i)-dveq_large(1)) /  ...    
                              (dveq_large(n_trans)-dveq_large(1));
  
  S(n_small+i).pha_mat_data = w * S(n_small+i).pha_mat_data + ... 
                                   (1-w) * S_i.pha_mat_data;
  S(n_small+i).ext_mat_data = w * S(n_small+i).ext_mat_data + ... 
                                   (1-w) * S_i.ext_mat_data;
  S(n_small+i).abs_vec_data = w * S(n_small+i).abs_vec_data + ... 
                                   (1-w) * S_i.abs_vec_data;
  M(n_small+i).diameter_max = w * M(n_small+i).diameter_max + ... 
                                   (1-w) * M_i.diameter_max;
end