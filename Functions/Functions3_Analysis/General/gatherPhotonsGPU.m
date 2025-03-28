%%  gatherPhotonsGPU.m
%   David Thompson 05-06-2023, last update 13-3-2024
%
%   Gathers the individual entries of a McDoppler photons struct from the GPU back into normal RAM
%
%   Inputs:     [struct]        detectedPhotons
%   Outputs:    [struct]        detectedPhotons

%%
function detectedPhotons = gatherPhotonsGPU(detectedPhotons)
arguments
    detectedPhotons struct
end
fieldNames = fieldnames(detectedPhotons);
for n = 1:length(fieldNames)
    detectedPhotons.(fieldNames{n}) = gather(detectedPhotons.(fieldNames{n}));
end
end