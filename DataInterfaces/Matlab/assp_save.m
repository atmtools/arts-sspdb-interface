% ASSP_SAVE   Stores a combination of S and M to a file
%
% FORMAT   assp_save(S,M,filename)
%
% IN   S          Vector of SingleScatteringData.
%      M          Vector of ScatteringMetaData
%      filename   Full path of file to create.

% 2016-12-25 Patrick Eriksson


function assp_save(S,M,filename)


save( filename, 'S', 'M' );
