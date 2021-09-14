%SSPNC_DEFINITION Returns a .nc format schema.
%
%   This function defines the format of the netCDF database.
%
%   Can be used either for checking extracted data or writing.
%
%   ncwriteschema(path, schema_nc) can be used to create a .nc-file with this
%   format. However, empty dimension and size-fields need to be filled
%   in advance (located at group level and in variables). They need to 
%   consistent everywhere as well. 
%
% FORMAT [ nc_schema ] = assd_ncFormat( orient_type )
%
% OUT   nc_schema       -nc-format definition. Struct format is equivalent
%                           to output of ncinfo.
%
% OPT   orient_type     Orient type. By giving this, dimensions length of
%                           relevant scattering matrices will be set in
%                           advance.
%

% 2016-11-03 Robin Ekelund: Created. Format is alpha v.1
% 2017-06-04 Robin Ekelund: Modyfied azimuthally_random matrix size.
%   Added 'shape_file' entry. Changed name of 'ADDA_version' to
%   'software_version'. Format is in alphaV.2
% 2017-08-31 Robin Ekelund: Is now in beta, version 0.9.
% 2017-11-20 Robin Ekelund: Added deflate level, version 0.9.2.
% 2018-02-15 Robin Ekelund: Added aspect ratio, version 1.0.0.
% 2020-10-11 Robin Ekelund: Added db_version as input. Replaced version
% with db_version and format_version, format version 1.1.0. Fixed wrong
% values for azimuthally_random
%
%%
function [ nc_schema ] = ssdb_ncFormat( orient_type, db_version )
% Change format_version only if format is changed. Keep as is if a new 
% database version (with new habits for instance) is released, but the 
% format is unchanged. If format is changed, set to the same version number
% as the new database version (if new database version is 1.3.4, format_version
% is 1.3.4 as well.)
format_version = '1.1.0'; 

if exist('db_version', 'var') ~= 1 || isempty(db_version)
  db_version = '';
else
  assert(ischar(db_version), 'db_version must be a string.');
end

if exist('orient_type', 'var') == 1
    otype_allowed = {'totally_random' 'azimuthally_random' 'general'};
    if ~any(strcmp(orient_type, otype_allowed))
        error('sspNc_definition:BadInput', 'Input orient_type not supported')
    end
else
    orient_type = '';
end

switch orient_type
    case ''
        n.phaMatElem = [];
        n.extMatElem = [];
        n.absVecElem = [];
    case 'totally_random'
        n.phaMatElem = 6;
        n.extMatElem = 1;
        n.absVecElem = 1;
    case 'azimuthally_random'
        n.phaMatElem = 16;
        n.extMatElem = 3;
        n.absVecElem = 2;
    case 'general'
        n.phaMatElem = 16;
        n.extMatElem = 7;
        n.absVecElem = 4;
    otherwise
        error('orient_type not supported')
end
    
%% Top level
nc_schema = createGroup('/');
nc_schema.Filename = [];
nc_schema.Format = 'netcdf4';
% Atributes
nc_schema = putAtt(nc_schema, 'date');  % Date of file creation
nc_schema = putAtt(nc_schema, 'format_version', format_version);   % Format version of the file
nc_schema = putAtt(nc_schema, 'db_version', db_version);   % Database version
% Variables

%% SSP level
struct_ssp = createGroup('SingleScatteringData');

% Atributes
struct_ssp = putAtt(struct_ssp, 'orient_type', orient_type);    % Orientation type

% Dimensions
dim_aaScat = createDim('aa_scat', []);
dim_zaScat = createDim('za_scat', []);
dim_aaInc = createDim('aa_inc', []);
dim_zaInc = createDim('za_inc', []);
dim_scatMatRow = createDim('scatMat_row', 4);
dim_scatMatCol = createDim('scatMat_Col', 4);


dim_phaMatElem = createDim('phaMatElem', n.phaMatElem);
dim_extMatElem = createDim('extMatElem', n.extMatElem);
dim_absVecElem = createDim('absVecElem', n.absVecElem);

struct_ssp = putDim(struct_ssp, [dim_aaScat dim_zaScat dim_aaInc...
    dim_zaInc dim_scatMatRow dim_scatMatCol dim_phaMatElem...
    dim_extMatElem dim_absVecElem]);

% Variables
struct_ssp = putVar(struct_ssp,  'frequency', [], 'double', 'Hz', '');
struct_ssp = putVar(struct_ssp,  'temperature', [], 'double', 'K', '');

struct_ssp = putVar(struct_ssp,  'aa_scat', dim_aaScat, 'double',...
    'degree', 'Azimuth scattering angle array');
struct_ssp = putVar(struct_ssp,  'za_scat', dim_zaScat, 'double',...
    'degree', 'Zenith scattering angle array');
struct_ssp = putVar(struct_ssp,  'aa_inc', dim_aaInc, 'double',...
    'degree', 'Azimuth incident angle array');
struct_ssp = putVar(struct_ssp,  'za_inc', dim_zaInc, 'double',...
    'degree', 'Zenith incident angle array');

struct_ssp = putVar(struct_ssp,  'phaMat_index', [dim_scatMatRow...
    dim_scatMatCol], 'int8', '', 'Phase matrix index', '', -1);
struct_ssp = putVar(struct_ssp,  'extMat_index', [dim_scatMatRow...
    dim_scatMatCol], 'int8', '', 'Extinction matrix index', '', -1);
struct_ssp = putVar(struct_ssp,  'absVec_index', [dim_scatMatRow],...
    'int8', '', 'Absorption vector index', '', -1);

struct_ssp = putVar(struct_ssp,  'phaMat_data', [dim_zaScat, dim_aaScat,...
    dim_zaInc, dim_aaInc, dim_phaMatElem], 'double', 'm^2',...
    'Phase matrix data', [], [], true, 4);
struct_ssp = putVar(struct_ssp,  'extMat_data', [dim_zaInc, dim_aaInc,...
    dim_extMatElem], 'double', 'm^2', 'Extinction matrix data');
struct_ssp = putVar(struct_ssp,  'absVec_data', [dim_zaInc, dim_aaInc,...
    dim_absVecElem], 'double', 'm^2', 'Absorption vector data');

% put group in output
nc_schema = putGrp(nc_schema, struct_ssp);

%% shape data level
struct_shape = createGroup('ShapeData');

% Atributes
struct_shape = putAtt(struct_shape, 'description'); % Description of shape
struct_shape = putAtt(struct_shape, 'source');  % Description of shape file source.
struct_shape = putAtt(struct_shape, 'shape_file');  % Name of shape_file
struct_shape = putAtt(struct_shape, 'refrIndex_model'); % Name of refractive index model
struct_shape = putAtt(struct_shape, 'habit_file_id');   % Habit string identifieŕ
struct_shape = putAtt(struct_shape, 'habit_id');    % Habit integer identifier
struct_shape = putAtt(struct_shape, 'phase');   % Phase of particle (ice, liquid or melting)
struct_shape = putAtt(struct_shape, 'refrIndex_homogenous_bool');   % Bool, true if refractive index is homogenous
struct_shape = putAtt(struct_shape, 'density_homogenous_bool'); % Bool, true if density is homogenous

% Dimensions

% Variables
struct_shape = putVar(struct_shape,  'diameter_max', [], 'double', 'm', 'Maximimum diameter, calculated as the minimum circumsphere diameter');
struct_shape = putVar(struct_shape,  'diameter_vol_eq', [], 'double', 'm', 'Volume equivalent diameter');
struct_shape = putVar(struct_shape,  'diameter_area_eq_aerodynamical', [], 'double', 'm', '');
struct_shape = putVar(struct_shape,  'mass', [], 'double', 'kg', '');
struct_shape = putVar(struct_shape,  'aspect_ratio', [], 'double', '', 'Aspect ratio of particle');
struct_shape = putVar(struct_shape,  'dpl', [], 'int32', '', 'Dipoles per wavelength', '', -1);
struct_shape = putVar(struct_shape,  'N_dipoles', [], 'int32', '', 'Number of dipoles', '', -1);
struct_shape = putVar(struct_shape,  'refrIndex_real', [], 'double', '', 'Real part of refractive index', 'Value only given for homogenous particles');
struct_shape = putVar(struct_shape,  'refrIndex_imag', [], 'double', '', 'Imaginary part of refractive index', 'Value only given for homogenous particles');
struct_shape = putVar(struct_shape,  'alpha', [], 'double', 'degree', 'Initial rotation angle alpha (zyz-notation)');
struct_shape = putVar(struct_shape,  'beta', [], 'double', 'degree', 'Initial rotation angle beta (zyz-notation)');
struct_shape = putVar(struct_shape,  'gamma', [], 'double', 'degree', 'Initial rotation angle gamma (zyz-notation)');

% put group in output
nc_schema = putGrp(nc_schema, struct_shape);
%% Calculation data level
struct_calc = createGroup('CalculationData');

% Atributes
struct_calc = putAtt(struct_calc, 'method');    % What type of calculation was used to derive the SSP
struct_calc = putAtt(struct_calc, 'software');  % What software were used for the calculations
struct_calc = putAtt(struct_calc, 'software_version');  % What was the software version
struct_calc = putAtt(struct_calc, 'system');    % On what computer system was the calculations performed
struct_calc = putAtt(struct_calc, 'n_nodes');   % How many computer cluster nodes were used
struct_calc = putAtt(struct_calc, 'n_cores');   % How many computer cores were used
struct_calc = putAtt(struct_calc, 'date_completion');   % Date of completion for ther calculations
struct_calc = putAtt(struct_calc, 'ADDA_eps');  % EPS if ADDA was used
struct_calc = putAtt(struct_calc, 'ADDA_avgParam_file');    % avgParams file if ADDA was used
struct_calc = putAtt(struct_calc, 'ADDA_scatParam_file');   % scatParam file if ADDA was used

% Dimensions

% Variables

% put group in output
nc_schema = putGrp(nc_schema, struct_calc);

end

function [struct_group] = createGroup(grp_name)
    struct_group = struct('Name', grp_name, 'Dimensions', [],...
        'Variables', [], 'Attributes', [], 'Groups', []);
end

function [struct_out] = putAtt(struct_out, name, value)
    if exist('value', 'var') ~= 1 || strcmp(value, '') || isempty(value)
        value = [];  
    end
    attribute = struct('Name', name, 'Value', value); 
    struct_out.Attributes = [struct_out.Attributes attribute];
end

function [dimension] = createDim(name, length)
    if exist('length', 'var') ~= 1 || isempty(length) || isnan(length) 
        length = [];
    end
    unlim_bool = false;
    dimension = struct('Name', name, 'Length', length, 'Unlimited',...
        unlim_bool); 
end

function [struct_out] = putVar(struct_out, name, dim_array, dataType,...
    unit, description, note, fillVal, shuffle, deflate_level)
    % checking
    dataType_allowed = {'double' 'int32' 'int8' 'char'};
    if ~any(strcmp(dataType, dataType_allowed))
        error('sspNc_definition:BadInput', 'Given datatype not supported')
    end
    if exist('fillVal', 'var') ~= 1 || isempty(fillVal)
        fillVal = NaN;
    end
    if exist('shuffle', 'var') ~= 1 || isempty(shuffle)
        shuffle = false;
    end
    if isempty(dim_array)
        size_vec = [];
    else
        size_vec = [dim_array.Length];
    end
    if exist('deflate_level', 'var') ~= 1 || isempty(deflate_level)
        deflate_level = [];
        chunkSize = [];
    else
        chunkSize = size_vec;
    end


    % Create variable
    variable = struct('Name', name, 'Dimensions', dim_array,...
        'Size', size_vec, 'Datatype', dataType,...
        'Attributes', [], 'ChunkSize', chunkSize, 'FillValue', fillVal,...
        'DeflateLevel', deflate_level, 'Shuffle', shuffle); 
    
    % Put attributes
    if exist('unit', 'var') == 1 && ~strcmp(unit, '') && ~isempty(unit)
        variable = putAtt(variable, 'unit', unit);  
    end
    if exist('description', 'var') == 1 && ~strcmp(description, '') &&...
            ~isempty(description)
        variable = putAtt(variable, 'description', description);  
    end
    if exist('note', 'var') == 1 && ~strcmp(note, '') && ~isempty(note)
        variable = putAtt(variable, 'note', note);  
    end
    % Create output
    struct_out.Variables = [struct_out.Variables variable];
end

function [struct_out] = putGrp(struct_out, struct_group)
    struct_out.Groups = [struct_out.Groups struct_group];  
end

function [struct_out]= putDim(struct_out, dim_array)
    struct_out.Dimensions = [struct_out.Dimensions dim_array];  
end
