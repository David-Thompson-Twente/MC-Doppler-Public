%%  interfaceNormalVectors.m
%   David Thompson 11-04-2022, last update 17-5-2024
%   
%   Gives incidence angles for any photons located on the medium interface.
%
%   Inputs:     [struct]        photons, with relevant properties:  photons.position
%                                                                   photons.direction
%                                                                   photons.mediumNum
%                                                                   photons.nextMedium
%
%               [struct]        media, with relevant properties:    media.axis
%                                                                   media.axisorigin
%               [string]        classArgument
%
%   Outputs:    [nx3 double]    normalVectors, with the normals to the interface at the n positions of all intersecting photons
%%
function normalVectors = interfaceNormalVectors(photons,media,classArgument)
arguments
    photons         struct
    media           struct
    classArgument   string
end

normalVectors = zeros(size(photons.position),classArgument);
interfacePhotons = ~(photons.mediumNum == photons.nextMedium);

if any(interfacePhotons,2)
    outsideNormals = photons.mediumNum < photons.nextMedium;
    insideNormals = photons.mediumNum > photons.nextMedium;

    secondMediumNum(outsideNormals) = photons.nextMedium(outsideNormals);
    secondMediumNum(insideNormals) = photons.mediumNum(insideNormals);
    
    posShifted = photons.position(:,interfacePhotons) - media.axisorigin(:,secondMediumNum(interfacePhotons));
    posShifted = posShifted - dot(posShifted,media.axis(:,secondMediumNum(interfacePhotons)),1).*media.axis(:,secondMediumNum(interfacePhotons));
    normalVectors(:,interfacePhotons) = posShifted./vecnorm(posShifted).*(outsideNormals(interfacePhotons) - insideNormals(interfacePhotons));
    planarInterfaces = (interfacePhotons & strcmp(media.shape(photons.nextMedium),'plane')) | (interfacePhotons & strcmp(media.shape(photons.mediumNum),'plane'));
    
    normalVectors(:,planarInterfaces) = media.axis(:,strcmp(media.shape,'plane')).*ones(1,sum(planarInterfaces),classArgument); %now very particular and for only one planar medium
    flippedNormalsPlane = dot(photons.direction,normalVectors) > 0;
    
    normalVectors(:,flippedNormalsPlane & planarInterfaces) = -normalVectors(:,flippedNormalsPlane & planarInterfaces);
end

end