% ASSD_READ_NCFILE   Extracts SSD data from a nc-file
%
%   This functions returns single scattering data from input netCDF file
%   with supported format.
%
%   NOTE: No format check implemented yet!
%
% FORMAT  [ SSD, freq, T, nodata ] = ssdb_read_ncfile( nc_file,...
%    freq_range, temp_range, bool_onlyFolders, nc_schema )
%
% OUT   SSD     Struct array containing SSD data. Fields names should be self explanatory.
%       freq    Frequency vector, corresponding to rows in SSD.
%       T       Temperature vector, corresponding to columns SSD.
%       nodata  Matrix corresponding to SSD. True if data is missing.
%
% IN    nc_file     Full path to nc-file.
%
% OPT   nc_schema       Format definition. By default, the funtion
%                       assd_ncFormat is used to generate the format.
%       freq_range      Frequency limits. Ignore data outside these limits.
%                       Default is [0,Inf]. If set to [], the default is used.
%       temp range      Temperature limits. Ignore data outside these limits.
%                       Default is [0,Inf]. If set to [], the default is used.
%       bool_onlyFolders    Option to only extract frequency and
%                           temperatures of folders

% 2016-10-31 Robin Ekelund: Created
% 2016-11-03 Robin Ekelund: Subfunction cleaned up. Reading is now based
%                               on the format definition in assd_ncFormat. 
% 2016-11-03 Robin Ekelund: Changed from using ncread to netcdf.getVar. 
% 2017-10-15 Robin Ekelund: Changed to return the internal netCDF internal
%                               structure. Included only functionality to
%                               filter out freuquencies and temperature,
%                               and option to only extract folder names.
% 2017-10-20 Robin Ekelund: Changed output format. Added output.
% 2017-10-23 Robin Ekelund: Fixed. Can now return empty data.
% 2017-10-23 Robin Ekelund: Fixed bug with bool_onlyFolders option.
% 2017-10-30 Robin Ekelund: Changed output for "folders only"-option.
% 2017-10-31 Robin Ekelund: Update. Can deal with files containing no internal folders.
%
function [ SSD, freq, T, nodata ] = ssdb_read_ncfile( nc_file,...
    freq_range, temp_range, bool_onlyFolders, nc_schema )
% setup
persistent nc_schema_p
if exist('nc_schema', 'var') == 1
    nc_schema_p = nc_schema;   
elseif isempty(nc_schema_p)
    nc_schema_p = ssdb_ncFormat;     
end
if exist('freq_range', 'var') == 0 || isempty(freq_range)
  freq_range = [ 0 Inf ];
end
if exist('temp_range', 'var') == 0 || isempty(temp_range)
  temp_range = [ 0 Inf ];
end
if exist('bool_onlyFolders', 'var') == 0
  bool_onlyFolders = 0;
end
% extract data
nc_id = netcdf.open(nc_file);
cleanup_obj = onCleanup(@()netcdf.close(nc_id));

folders = get_internalStruct(nc_id, nc_schema_p, freq_range, temp_range, bool_onlyFolders);
if isempty(folders)
    SSD = [];
    freq = [];
    T = [];
    nodata = [];
elseif isfield(folders, 'date') && length(folders) == 1
    SSD = folders;
else
    [SSD, freq, T, nodata] = order_data(folders);
end
end

function [SSD, freq_ref, T_ref, nodata] = order_data(folders)
    freq = [folders.freq];
    T = [folders.temp];
    freq_ref = unique(freq);
    T_ref = unique(T);
    N_f = length(freq_ref);
    N_T = length(T_ref);
    fi = fieldnames(folders(1).data);
    tmp = cell(1, 2*length(fi));
    tmp(1:2:end) = fi;
    SSD(N_f, N_T) = struct(tmp{:});
    nodata = true(N_f, N_T);
    for i = 1:length(freq_ref)
        for n = 1:length(T_ref)
            I = freq_ref(i) == freq & T_ref(n) == T;
            if sum(I) == 1
                nodata(i, n) = false;
                SSD(i, n) = folders(I).data;
            elseif sum(I) > 1
                error('More than 1 frequency and temperature combination found.')
            end
        end
    end
end

function [folders] = get_internalStruct(nc_id, nc_schema, freq_range, temp_range, bool_onlyFolders)
% Look for internal folders:
grp_id_veq = netcdf.inqGrps(nc_id);
i = 0;
folders(length(grp_id_veq)) = struct('name', [], 'data', [], 'freq', [], 'T', []);
for grp_id = grp_id_veq
    i = i + 1;
    grp_name = netcdf.inqGrpName(grp_id);

    if ismember(grp_name, {nc_schema.Groups.Name}) == 1    % If group is matching nc_schema, stop looking for internal folders.
        folders = get_ncData(nc_id, nc_schema);
        return
    end
    tmp = textscan(grp_name, 'Freq%fGHz_T%fK');
    if isempty(tmp{1}) == 0
        freq = tmp{1} * 1e9;
        temp = tmp{2};
        folders(i).temp = temp;
        folders(i).freq = freq;
        if freq < freq_range(1) || freq > freq_range(2) ||...
            temp < temp_range(1) || temp > temp_range(2)
            continue 
        end
    end
    folders(i).name = grp_name;
    if bool_onlyFolders == 0
        folders(i).data = get_internalStruct(grp_id, nc_schema, freq_range, temp_range);
    else
        folders(i).data = struct('freq', freq, 'temp', temp, 'name', grp_name);
    end
end
I = cellfun(@(x)0==isempty(x), {folders.name});
folders = folders(I);
end

function [nc_Data] = get_ncData(nc_id, nc_schema)
nc_global = netcdf.getConstant('NC_GLOBAL');
% Read attributes
for att_ref = nc_schema.Attributes
    nc_Data.(att_ref.Name) = netcdf.getAtt(nc_id, nc_global, att_ref.Name);
end
% Read variables
for var_ref = nc_schema.Variables
    var_id = netcdf.inqVarID(nc_id, var_ref.Name);
    nc_Data.(var_ref.Name) = netcdf.getVar(nc_id,var_id);
end
% Recursive call group
for grp_ref = nc_schema.Groups
    grp_id = netcdf.inqNcid(nc_id, grp_ref.Name);    
    nc_Data.(grp_ref.Name) = get_ncData(grp_id, grp_ref);
end
end

