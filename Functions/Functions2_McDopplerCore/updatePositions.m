%%  updatePositions.m
%   David Thompson 13-04-2022, last update 17-5-2024
%   Updates photon positions and path lengths, to full propagation path
%   length or up to an interface
%
%   Inputs:     [struct]        photons
%               [struct]        media
%               [struct]        sim
%               [string]        classArgument
%
%   Outputs:    [struct]        photons with updated:   photons.position
%                                                       photons.pathLengths
%                                                       photons.alive
%                                                       photons.amplitude
%                                                       photons.distanceTravelled
%                                                       photons.opticalPathLength
%               [3xn double]    normalVectors

%%
function [photons,normalVectors] = updatePositions(photons,media,sim,classArgument)
arguments
    photons             struct
    media               struct
    sim                 struct
    classArgument       string
end

photons = pathLengths(media,photons,classArgument);

if length(media.name)>1   % Only calculate intersection if there is a 2nd medium to intersect with
    photons = interfaceIntersect(photons,media,sim,classArgument);
end
normalVectors = interfaceNormalVectors(photons,media,classArgument);
photons = photonEscapeCheck(photons,sim);
photons = photonEscapeTruncate(photons,sim,media,classArgument);

photons.position(:,photons.pathLengths ~= 0) = photons.position(:,photons.pathLengths ~= 0) ...
    + photons.direction(:,photons.pathLengths ~= 0).*photons.pathLengths(:,photons.pathLengths ~= 0);

photons.amplitude(photons.pathLengths ~= 0) = photons.amplitude(photons.pathLengths ~= 0)...
    .*exp(horzcat(media.absorptionCoeff(photons.mediumNum(photons.pathLengths ~= 0))).*-photons.pathLengths(photons.pathLengths ~= 0));

photons.distanceTravelled(photons.pathLengths ~= 0) = photons.distanceTravelled(photons.pathLengths ~= 0) + photons.pathLengths(photons.pathLengths ~= 0);
photons.opticalPathLength(photons.pathLengths ~= 0) = photons.opticalPathLength(photons.pathLengths ~= 0) + photons.pathLengths(photons.pathLengths ~= 0).*media.refractiveIndex(photons.mediumNum(photons.pathLengths ~= 0));
% Ensure that all vectors that should be unit are corrected for rounding errors to prevent later crashes
photons.direction = photons.direction./vecnorm(photons.direction);
photons.polarisationVector = photons.polarisationVector./vecnorm(photons.polarisationVector);
end