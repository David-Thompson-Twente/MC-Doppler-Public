%%  updateDirections.m
%   David Thompson 19-04-2022, last update 19-03-2025
%   Updates directions of reflected, refracted and scattered photons
%
%   Inputs:     [struct]        photons
%               [struct]        media
%               [3xn double]    normalVectors
%               [struct]        anglePickFunction
%               [string]        polarisation
%               [string]        classArgument
%
%   Outputs:    [struct]        photons, with updated   photons.direction
%                                                       photons.DopplerShift
%                                                       photons.timesScattered
%                                                       photons.timesReflected
%                                                       photons.timesRefracted

%%
function photons = updateDirections(photons,media,normalVectors,anglePickFunction,polarisation,classArgument)
arguments
    photons             struct
    media               struct
    normalVectors       (3,:)   {mustBeNumeric(normalVectors)}
    anglePickFunction   struct
    polarisation        string  {mustBeMember(polarisation,{'unpolarised','polarised'})}
    classArgument       string
end
% section for handling reflection and refraction at a medium interface
crossingPhotons = photons.nextMedium ~= photons.mediumNum & photons.alive;
if ~isempty(crossingPhotons(crossingPhotons == 1))
    [reflectionCoefficients,indexRatio] = reflectionCoeffs(photons,media,normalVectors,crossingPhotons,classArgument);
    photons = updateDirectionsReflectRefract(photons,normalVectors,crossingPhotons,reflectionCoefficients,indexRatio,polarisation,classArgument);
end
% section for handling scattering
scatteringPhotons = ~crossingPhotons & [media.scatterCoeff(photons.mediumNum)] > 0 & photons.alive;
if sum(scatteringPhotons)>0
    switch polarisation
        case 'unpolarised'
            scatteringAngles = scatteringAngleQuantile(anglePickFunction,photons,scatteringPhotons,media,classArgument);
        case 'polarised'
            scatteringAngles = scatteringAnglePolarised(anglePickFunction,photons,scatteringPhotons,media,classArgument);
    end
    if ~isempty(scatteringAngles(scatteringAngles < 0 | isnan(scatteringAngles)))
        disp('negative or NaN angles found and corrected')
        photons.alive(~(scatteringAngles(1,:) >= 0)) = 0;
        photons.alive(~(scatteringAngles(2,:) >= 0)) = 0;
        scatteringPhotons(~(scatteringAngles(1,:) >= 0)) = 0;
        scatteringPhotons(~(scatteringAngles(2,:) >= 0)) = 0;
        scatteringAngles(1,~(scatteringAngles(1,:) >= 0)) = 0;
        scatteringAngles(2,~(scatteringAngles(2,:) >= 0)) = 0;
    end
    photons = updateDirectionsScatter(photons,scatteringPhotons,scatteringAngles,media,classArgument);
end
photons.pathLengths(photons.pathLengths ~= 0) = 0;
end