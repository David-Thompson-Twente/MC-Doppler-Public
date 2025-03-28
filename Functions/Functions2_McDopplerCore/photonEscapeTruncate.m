%%  photonEscapeTruncate.m
%   David Thompson 05-04-2022, Last update 17-5-2024
%   Sets position of escaped photons to simulation boundary crossing point.
%
%   Inputs:     [struct]    photons with relevant properties:   photons.position
%                                                               photons.direction
%                                                               photons.alive
%                                                               photons.amplitude
%                                                               photons.pathLength
%                                                               photons.mediumNum
%
%               [struct]    sim with relevant property:         sim.radius
%                                                               sim.shape
%                                                               sim.halfLength
%                                                               sim.halfWidth
%                                                               sim.halfDepth
%               [struct]    media with relevant property:       media.absorptionCoeff
%               [string]    classArgument
%
%   Outputs:    [struct]    photons, with updated positions to truncate escaped 
%                                    photon paths at the simulation boundary

%%
function photons = photonEscapeTruncate(photons,sim,media,classArgument)
arguments
    photons             struct
    sim                 struct
    media               struct
    classArgument       string
end
criterion = photons.alive == 0 & photons.amplitude > 0 & vecnorm(photons.position) <= sim.radius;
switch sim.shape
    case 'sphere'
        photons.pathLengths(criterion) = -dot(photons.direction(:,criterion),photons.position(:,criterion),1)...
            +sqrt(dot(photons.direction(:,criterion),photons.position(:,criterion),1).^2 ...
            -(vecnorm(photons.position(:,criterion)).^2 - sim.radius^2));
    case 'cylinder'
        vecXZ = photons.direction;
        vecXZ(2,:) = 0;
        unitVecXZ = vecXZ./vecnorm(vecXZ);
        projectionRatio = dot(photons.direction,unitVecXZ,1);
        posXZ = photons.position;
        posXZ(2,:) = 0;
        
        dotProd = dot(unitVecXZ,posXZ,1);
        discriminant = dotProd.^2 - (vecnorm(posXZ).^2 - sim.radius.^2);
        discriminant(discriminant < 0) = NaN;
        radialDistances =  (-dotProd + sqrt(discriminant))./projectionRatio;

        negPlaneCentres = ones(size(photons.position),classArgument).*[0;-sim.halfLength;0];
        negNormalVectors = ones(size(photons.position),classArgument).*[0;1;0];
        posPlaneCentres = ones(size(photons.position),classArgument).*[0;sim.halfLength;0];
        posNormalVectors = ones(size(photons.position),classArgument).*[0;-1;0];
        negAxialDistances = dot((negPlaneCentres - photons.position),negNormalVectors)./dot(photons.direction,negNormalVectors,1);
        posAxialDistances = dot((posPlaneCentres - photons.position),posNormalVectors)./dot(photons.direction,posNormalVectors,1);
        
        distanceMat = [radialDistances ; negAxialDistances ; posAxialDistances];
        distanceMat(distanceMat < 0) = NaN;        
        photons.pathLengths(criterion) = min(distanceMat(:,criterion));        
    case 'cube'       
        xNegPlaneCentres = ones(size(photons.position)).*[-sim.halfWidth;0;0];
        xPosPlaneCentres = ones(size(photons.position)).*[sim.halfWidth;0;0];
        yNegPlaneCentres = ones(size(photons.position)).*[0;-sim.halfLength;0];
        yPosPlaneCentres = ones(size(photons.position)).*[0;sim.halfLength;0];
        zNegPlaneCentres = ones(size(photons.position)).*[0;0;-sim.halfDepth];
        zPosPlaneCentres = ones(size(photons.position)).*[0;0;sim.halfDepth];
        
        xNegNormalVectors = ones(size(photons.position)).*[1;0;0];
        xPosNormalVectors = ones(size(photons.position)).*[-1;0;0];
        yNegNormalVectors = ones(size(photons.position)).*[0;1;0];
        yPosNormalVectors = ones(size(photons.position)).*[0;-1;0];
        zNegNormalVectors = ones(size(photons.position)).*[0;0;1];
        zPosNormalVectors = ones(size(photons.position)).*[0;0;-1];

        xNegAxialDistances = dot((xNegPlaneCentres - photons.position),xNegNormalVectors,1)./dot(photons.direction,xNegNormalVectors,1);
        xPosAxialDistances = dot((xPosPlaneCentres - photons.position),xPosNormalVectors,1)./dot(photons.direction,xPosNormalVectors,1);
        yNegAxialDistances = dot((yNegPlaneCentres - photons.position),yNegNormalVectors,1)./dot(photons.direction,yNegNormalVectors,1);
        yPosAxialDistances = dot((yPosPlaneCentres - photons.position),yPosNormalVectors,1)./dot(photons.direction,yPosNormalVectors,1);
        zNegAxialDistances = dot((zNegPlaneCentres - photons.position),zNegNormalVectors,1)./dot(photons.direction,zNegNormalVectors,1);
        zPosAxialDistances = dot((zPosPlaneCentres - photons.position),zPosNormalVectors,1)./dot(photons.direction,zPosNormalVectors,1);

        distanceMat = [xNegAxialDistances ; xPosAxialDistances ; yNegAxialDistances; yPosAxialDistances ; zNegAxialDistances ; zPosAxialDistances];
        distanceMat(distanceMat < 0) = NaN;        
        photons.pathLengths(criterion) = min(distanceMat(:,criterion));        
end
photons.amplitude(criterion) = photons.amplitude(criterion)...
    .*exp(horzcat(media.absorptionCoeff(photons.mediumNum(criterion))).*-photons.pathLengths(criterion));
end
