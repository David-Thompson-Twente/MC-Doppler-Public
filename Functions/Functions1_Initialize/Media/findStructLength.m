%%  findStructLength.m
%   Wietske Verveld 27-03-2024, last update 27-03-2024
%
%   Finds maximum length of any property of the struct
%
%   Inputs:     [struct]        inputStruct
%
%   Optional:   [cell]          omitFields, cell array of field names not to copy
%
%   Outputs:    [1x1 double]    maxLength 

%%
function maxLength = findStructLength(inputStruct,options)
arguments
    inputStruct         struct
    options.omitFields  cell    =   {};
end
fieldNames = fieldnames(inputStruct)';

maxLength = 0;
for n = 1:length(fieldNames)
    if ~ismember(fieldNames(n),options.omitFields)
        maxLength = max(maxLength,size(inputStruct.(fieldNames{n}),2));
    end
end