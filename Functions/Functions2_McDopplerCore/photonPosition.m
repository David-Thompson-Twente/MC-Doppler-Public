%%  photonPosition.m
%   David Thompson 03-10-2022, last update 28-3-2024
%   Inputs:     [struct]        source, with relevant properties:   source.position
%                                                                   source.opticalAxis
%                                                                   source.posPDF
%               [struct]        sim, with relevant property:        sim.numPhotons
%               [string]        classArgument
%
%   Outputs:    [3xn double]    position

%%
function position = photonPosition(source,sim,classArgument)
arguments
    source          struct
    sim             struct
    classArgument   string  {mustBeMember(classArgument,{'single','double','gpuArray'})}
end
switch source.sourceType
    case {'focussed', 'ratio'}
        photonCoord1 = source.posPDF(sim.numPhotons);
        photonCoord2 = source.posPDF(sim.numPhotons);
        sourcePlaneVec1 = rotateRodrigues(source.opticalAxis,pi/2,0,[],classArgument);
        sourcePlaneVec2 = rotateRodrigues(source.opticalAxis,pi/2,pi/2,sourcePlaneVec1,classArgument);
        position = source.position + sourcePlaneVec1.*photonCoord1 + sourcePlaneVec2.*photonCoord2;
    case 'truncated ratio'
        photonCoord = ppval(source.posPDF,rand(1,sim.numPhotons));
        sourcePlaneVec1 = rotateRodrigues(source.opticalAxis,pi/2,0,[],classArgument);
        sourcePlaneVec2 = rotateRodrigues(source.opticalAxis,pi/2,pi/2,sourcePlaneVec1,classArgument);
        position = source.position + photonCoord.*rotateRodrigues(sourcePlaneVec1.*ones(1,sim.numPhotons),2*pi*rand(1,sim.numPhotons),0,sourcePlaneVec2,classArgument);
    case {'point', 'pencil'}
        position = source.position.*ones(1,sim.numPhotons,classArgument);
end
end