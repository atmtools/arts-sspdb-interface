% ASSP_LOAD   Loads a combination of S and M from a file
%
% FORMAT   [S,M] = assp_load(filename)
%
% OUT  S          Vector of SingleScatteringData.
%      M          Vector of ScatteringMetaData
%  IN  filename   Full path to file to load

% 2016-12-25 Patrick Eriksson


function [S,M] = assp_load(filename)


load( filename );
