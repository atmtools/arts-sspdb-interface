% ASSP_EXTEND_F  Extension of frequency coverage of ASSP data
%
%   The function provides a fast way to ensure that the frequency coverage of
%   the ASSP data reaches the given limit. If there is a gap between highest
%   frequency and *f_limit*, the data are extended in one of two simple
%   manners.
%
%   Default is to copy the scattering properties at highest existing
%   frequency, and include with *f_limit* as the assigned frequency.
%
%   If *zerofill* is set to true, the data are instead extended with zeros.
%   Here two frequencies are added. The first point is the highest existing
%   frequency + 1 GHz, and the second point is *f_limit*. The value 1 GHz is
%   adjusted when necessary to always get an increasing frequency grid.
%
% FORMAT   S = assp_extend_f( S, f_limit [, zerofill ] )
%
% OUT   S          As input but possibly with extended frequency coverage.
%  IN   S          Input ASSP structure(s). I.e. *S* can hold multiple structures.
%       f_limit    The frequency limit to meet.
% OPT   zerofill   Flag to instead fill with zeros, see above. Default is false.

% 2016-11-03 Patrick Eriksson


function S = assp_extend_f( S, f_limit, zerofill )
%
if nargin < 3  |  isempty(zerofill)
  zerofill = false;   
end


for i = 1 : size(S,1)
  for j = 1 : size(S,2)
    
    if S(i,j).f_grid < f_limit 
        
      % Copy existing data 
      %
      e = S(i,j).ext_mat_data;
      a = S(i,j).abs_vec_data;
      p = S(i,j).pha_mat_data;
      %
      nf = length( S(i,j).f_grid );

      if zerofill
          
        S(i,j).f_grid(end+1) = min( [ S(i,j).f_grid(end) + 1e9, ...
                                    mean([S(i,j).f_grid(end),f_limit]) ] );
        S(i,j).f_grid(end+1) = f_limit;
        
        s    = size( e );
        s(1) = nf + 2;
        S(i,j).ext_mat_data = zeros( s );
        S(i,j).ext_mat_data(1:nf,:,:,:,:) = e;

        s    = size( a );
        s(1) = nf + 2;
        S(i,j).abs_vec_data = zeros( s );
        S(i,j).abs_vec_data(1:nf,:,:,:,:) = a;

        s    = size( p );
        s(1) = nf + 2;
        S(i,j).pha_mat_data = zeros( s );
        S(i,j).pha_mat_data(1:nf,:,:,:,:,:,:) = p;
        
      else
          
        S(i,j).f_grid(end+1) = f_limit;
        
        s    = size( e );
        s(1) = nf + 1;
        S(i,j).ext_mat_data = zeros( s );
        S(i,j).ext_mat_data(1:nf,:,:,:,:) = e;
        S(i,j).ext_mat_data(end,:,:,:,:)  = e(end,:,:,:,:);

        s    = size( a );
        s(1) = nf + 1;
        S(i,j).abs_vec_data = zeros( s );
        S(i,j).abs_vec_data(1:nf,:,:,:,:) = a;
        S(i,j).abs_vec_data(end,:,:,:,:)  = a(end,:,:,:,:);

        s    = size( p );
        s(1) = nf + 1;
        S(i,j).pha_mat_data = zeros( s );
        S(i,j).pha_mat_data(1:nf,:,:,:,:,:,:) = p;
        S(i,j).pha_mat_data(end,:,:,:,:,:,:)  = p(end,:,:,:,:,:,:);
      
      end
        
    end
  
  end
end
