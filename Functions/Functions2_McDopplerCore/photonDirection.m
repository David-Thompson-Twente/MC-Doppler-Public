%%  photonDirection.m
%   David Thompson 03-10-2022, last update 28-3-2024
%   Generates direction unit vectors for sim.numPhotons photons
%
%   Inputs:     [struct]        source with relevant properties source.opticalAxis
%                                                               source.maxAngle
%                                                               source.sourceType
%                                                               source.focusPosition
%                                                               source.focusRatio
%                                                               source.position
%               [struct]        sim with relevant property:     sim.numPhotons
%               [3xn double]    photonPosition
%               [string]        classArgument
%
%   Outputs:    [3xn double]    directions

%%
function directions = photonDirection(source,sim,photonPosition,classArgument)
arguments
    source          struct
    sim             struct
    photonPosition  (3,:)   {mustBeNumeric}
    classArgument   string  {mustBeMember(classArgument,{'single','double','gpuArray'})}
end
switch source.sourceType
    case 'point'
        directionsNonUnit = 2*rand(3,sim.numPhotons) - [1;1;1];
        directions = directionsNonUnit./vecnorm(directionsNonUnit);
    case 'focussed'
        directionNonUnit = sign(dot(photonPosition - source.focusPosition,source.opticalAxis*ones(1,sim.numPhotons,classArgument),1)).*(photonPosition - source.focusPosition);
        directions = directionNonUnit./vecnorm(directionNonUnit);
    case 'ratio'
        directionNonUnit = (source.focusPosition + source.focusRatio.*(photonPosition - source.position)) - photonPosition;
        directions = directionNonUnit./vecnorm(directionNonUnit);
    case 'truncated ratio'
        directionNonUnit = (source.focusPosition + source.focusRatio.*(photonPosition - source.position)) - photonPosition;
        directions = directionNonUnit./vecnorm(directionNonUnit);
    case 'pencil'
        directions = source.opticalAxis.*ones(1,sim.numPhotons,classArgument);
end
end