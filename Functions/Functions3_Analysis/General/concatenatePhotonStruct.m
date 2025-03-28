%%  concatenatePhotonStruct.m
%   David Thompson 03-11-2022, last update 13-3-2024
%
%   Concatenates two structs of detectedPhotons field by field
%
%   Inputs:     [struct]        detectedPhotonsSlice
%               [struct]        slicePhotons
%
%   Outputs:    [struct]        detectedPhotonsSlice, with slicePhotons added to each field

%%
function concatenatedStruct = concatenatePhotonStruct(concatenatedStruct,structToBeAdded)
arguments
    concatenatedStruct    struct
    structToBeAdded       struct
end
fieldNames = fieldnames(structToBeAdded);
for n = 1:length(fieldNames)
    concatenatedStruct.(fieldNames{n}) = horzcat(concatenatedStruct.(fieldNames{n}),structToBeAdded.(fieldNames{n}));
end
end