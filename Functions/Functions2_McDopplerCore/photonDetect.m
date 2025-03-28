%%  photonDetect.m
%   David Thompson 29-04-2022, last update 31-10-2024
%   Removes any dead photons from photons struct, places them in
%   detectedPhotons struct;
%
%   Inputs:     [struct]    photons
%               [struct]    detectedPhotons
%
%   Outputs:    [struct]    photons
%               [struct]    detectedPhotons

%%
function [photons, detectedPhotons] = photonDetect(photons,detectedPhotons)
arguments
    photons struct
    detectedPhotons struct
end
criterion = photons.alive == 0;
fieldNames = fieldnames(detectedPhotons);
ndetect = sum(criterion);
for n = 1:length(fieldNames)
    detectedPhotons.(fieldNames{n})(:,(end+1):(end+ndetect)) = photons.(fieldNames{n})(:,criterion);
    photons.(fieldNames{n}) = photons.(fieldNames{n})(:,~criterion);
end
photons.alive = photons.alive(~criterion);
end