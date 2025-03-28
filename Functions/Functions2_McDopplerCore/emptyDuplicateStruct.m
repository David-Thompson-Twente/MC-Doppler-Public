%%  emptyDuplicateStruct.m
%   David Thompson 23-02-2024, last update 13-3-2024
%
%   Generates empty struct with same fields as input struct
%
%   Inputs:     [struct]        inputStruct
%
%   Optional:   [cell]          omitFields, cell array of field names not to copy
%
%   Outputs:    [struct]        outputStruct 

%%
function outputStruct = emptyDuplicateStruct(inputStruct,options)
arguments
    inputStruct         struct
    options.omitFields  cell    =   {};
end
fieldNames = fieldnames(inputStruct)';
fieldNames(ismember(fieldNames,options.omitFields)) = [];
fieldNames{2,1} = [];
outputStruct = struct(fieldNames{:});
end