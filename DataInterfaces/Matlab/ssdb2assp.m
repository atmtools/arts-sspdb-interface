% SSDB2ASSP   Converts internal SSD to ARTS single scattering properties
%
%    The input SSD data should be imported using the grid_sort option.
%
% FORMAT  [S,M] = ssdb2assp( SSD, freq, temp, nodata [, interpm ])
%
% OUT   S        The data as ARTS-XML SingleScatteringData
%       M        The data as ARTS-XML ScatteringMetaData
%  IN   SSD      Input single scattering data.
%       freq     The frequency of each SSD element.
%       temp     The temperature of each SSD element.
%       nodata   Boolean vector/matrix, flagging empty SSD elements.
% OPT   interpm  Method for zenith angle interpolation. Default is 'linear'.

% 2016-11-02 Patrick Eriksson
% 2016-11-03 Robin Ekelund: Paramater name change.
% 2016-11-17 Robin Ekelund: Changed due to attribute name change
%                               (orient_type).
% 2017-02-13 Robin Ekelund: Fixed bug when nodata is allowed.

function [S,M] = ssdb2assp( SSD, freq, temp, nodata, interpm )
%
if nargin < 5  |  isempty(interpm)
  interpm = 'linear';   
end


% Basic sanity check of input
%
if ~isfield( SSD, 'SingleScatteringData' )
  error( 'Input *SSD* does not seem to hold expected data formnat.' ); 
end   
if size(SSD,1) ~= length(freq)
  error( 'Mismatch in size between *SSD* and *freq*.' ); 
end
if size(SSD,2) ~= length(temp)
  error( 'Mismatch in size between *SSD* and *temp*.' ); 
end


% M
%
ii = find(~nodata(:),1);
M.version             = 3;
M.description         = [ 'Meta data for ', SSD(ii).ShapeData.description ];
M.source              = SSD(ii).ShapeData.source;
M.refr_index          = SSD(ii).ShapeData.refrIndex_model;
M.mass                = SSD(ii).ShapeData.mass;
M.diameter_max        = SSD(ii).ShapeData.diameter_max;
M.diameter_volume_equ = SSD(ii).ShapeData.diameter_vol_eq;
M.diameter_area_equ_aerodynamical = SSD(ii).ShapeData.diameter_area_eq_aerodynamical;


% Basic data of S
%
S.version     = 3;
S.ptype       = SSD(ii).SingleScatteringData.orient_type;
S.description = SSD(ii).ShapeData.description;

% Fill S according to ptype
%
nf = length( freq );
nt = length( temp );
%
if ismember( S.ptype, {'totally_random' 'azimuthally_random'} )
  %
  S.f_grid  = freq;
  S.T_grid  = temp;
  %
  % Create a common za_grid
  S.za_grid = [];
  S.aa_grid = [];
  for f = 1 : nf
    for t = 1 : nt
      if ~nodata(f,t)
        assert( SSD(f,t).SingleScatteringData.za_scat(1) == 0 );
        assert( SSD(f,t).SingleScatteringData.za_scat(end) == 180 );
        S.za_grid =  union( S.za_grid, SSD(f,t).SingleScatteringData.za_scat );
        if strcmp( S.ptype, 'totally_random' )
          S.aa_grid = [];
        elseif strcmp( S.ptype, 'azimuthally_random' )
          S.aa_grid =  union( S.aa_grid, SSD(f,t).SingleScatteringData.aa_scat );
        end
      end
    end
  end
  %
  nza            = length( S.za_grid );
  if strcmp( S.ptype, 'totally_random' )
    npha = 6;
    next = 1;
    nabs = 1;
    nza_in = 1;
    naa = 1;
  elseif strcmp( S.ptype, 'azimuthally_random' )
    npha = 16;
    next = 3;
    nabs = 2;
    nza_in = nza;
    naa = length( S.aa_grid );
  end
  S.pha_mat_data = nan( nf,nt,nza,naa,nza_in,1,npha );  
  S.ext_mat_data = nan( nf,nt,nza_in,next );
  S.abs_vec_data = nan( nf,nt,nza_in,nabs );
  %
  for f = 1 : nf
    for t = 1 : nt
      if ~nodata(f,t)

        % A number of assert to double-check the input
        %
        assert( abs(M.mass-SSD(f,t).ShapeData.mass) < 1e-12  ); 
        %
        assert( strcmp( S.ptype, SSD(f,t).SingleScatteringData.orient_type ) );
        assert( abs(freq(f)-SSD(f,t).SingleScatteringData.frequency) < 1e3 );
        assert( abs(temp(t)-SSD(f,t).SingleScatteringData.temperature) < 0.001 );
        if strcmp( S.ptype, 'totally_random' )
          assert( all( int8([1 2 0 0 2 3 0 0 0 0 4 -5 0 0 5 6])' - ...
                     SSD(f,t).SingleScatteringData.phaMat_index(:) == 0 ) );
        elseif strcmp( S.ptype, 'azimuthally_random' )
          assert( all( int8([1 5 9 13 2 6 10 14 3 7 11 15 4 8 12 16])' - ...
                     SSD(f,t).SingleScatteringData.phaMat_index(:) == 0 ) );
        end
        
        % Fill scattering fields
        %
        if strcmp( S.ptype, 'totally_random' )
          S.pha_mat_data(f,t,:,:,:,:,:) = interp1( ...
              SSD(f,t).SingleScatteringData.za_scat, ...
              SSD(f,t).SingleScatteringData.phaMat_data, S.za_grid, interpm );
          S.ext_mat_data(f,t,:,:) = SSD(f,t).SingleScatteringData.extMat_data;
          S.abs_vec_data(f,t,:,:) = SSD(f,t).SingleScatteringData.absVec_data;

        elseif strcmp( S.ptype, 'azimuthally_random' )
          [za, aa, zain] = ndgrid(SSD(f,t).SingleScatteringData.za_scat,...
            SSD(f,t).SingleScatteringData.aa_scat,...
            SSD(f,t).SingleScatteringData.za_inc);
          [zaq, aaq, zainq] = ndgrid(S.za_grid, S.aa_grid,S.za_grid);
          for ind = 1:npha
            S.pha_mat_data(f,t,:,:,:,:,ind) = interpn( za,aa,zain, ...
                SSD(f,t).SingleScatteringData.phaMat_data(:,:,:,:,ind),...
                zaq,aaq,zainq,interpm );
          end
          S.ext_mat_data(f,t,:,:) = interp1(SSD(f,t).SingleScatteringData.za_scat,...
            SSD(f,t).SingleScatteringData.extMat_data,S.za_grid,interpm);
          S.abs_vec_data(f,t,:,:) = interp1(SSD(f,t).SingleScatteringData.za_scat,...
            SSD(f,t).SingleScatteringData.absVec_data,S.za_grid,interpm);
        end
      end
    end
  end
  
else
  error( 'Found a ptype that is not handled: %s', S.ptype );
end


