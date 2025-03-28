%%  interfaceIntersect.m
%   David Thompson 11-04-2022, last update 19-03-2025
%
%   Places any intersecting photons on medium interface with updated value for pathLengths
%
%   Inputs:     [struct]        photons, with relevant properties:  photons.position
%                                                                   photons.direction
%                                                                   photons.alive
%                                                                   photons.amplitude
%                                                                   photons.mediumNum
%                                                                   photons.nextMedium
%                                                                   photons.pathLengths
%                                                                   photons.distanceTravelled
%                                                                   photons.opticalPathLength
%               [struct]        media, with relevant properties:    media.name
%                                                                   media.relevantDimension
%                                                                   media.axis
%                                                                   media.axisorigin
%                                                                   media.absorptionCoeff
%                                                                   media.refractiveIndex
%                                                                   media.surroundingMedium
%               [struct]        sim with relevant properties:       sim.GPU
%                                                                   sim.shape
%                                                                   sim.radius      (for sphere)
%                                                                   sim.halfLength  (for cylinder and cube shape)
%                                                                   sim.halfWidth   (for cube shape)
%                                                                   sim.halfDepth   (for cube shape)
%               [string]        classArgument
%
%   Optional:   [1x1 double]    toleranceLevel
%
%   Outputs:    [struct]        photons, with updated               photons.position
%                                                                   photons.alive
%                                                                   photons.nextMedium
%                                                                   photons.pathLengths set to zero for any interfacing photons
%                                                                   photons.distanceTravelled
%                                                                   photons.opticalPathLength

%%
function photons = interfaceIntersect(photons,media,sim,classArgument,options)
arguments
    photons                 struct
    media                   struct
    sim                     struct
    classArgument           string
    options.toleranceLevel  (1,1) {mustBeNumeric} = 1000*eps(1)
end

photons.nextMedium = photons.mediumNum;
medNums = repelem((2:length(media.name)),2);
mediumDims = vertcat(media.relevantDimension(:,medNums));
mediumAxes = vertcat(media.axis(:,medNums));
mediumAxOrigins = vertcat(media.axisorigin(:,medNums));
medInd = medNums;
medInd(2:2:end) = media.surroundingMedium(medNums(2:2:end));

if strcmp(classArgument,'gpuArray')
    medInd = gpuArray(medInd);
end

distanceMat = zeros(length(medNums),length(photons.mediumNum),classArgument);

for n_mediumRow = 1:length(medNums)
    switch media.shape{medNums(n_mediumRow)}
        case 'cylinder'
            posShifted = photons.position-mediumAxOrigins(:,n_mediumRow);
            f = dot(photons.direction,mediumAxes(:,n_mediumRow).*ones(1,length(photons.alive),classArgument),1); % extra help variable to reduce double dot products in a and b
            g = dot(posShifted,mediumAxes(:,n_mediumRow).*ones(1,length(photons.alive),classArgument),1); % extra help variable to reduce double dot products in b and c
            a = 1 - f.^2;
            b = 2.*(dot(photons.direction,posShifted,1) - f.*g);
            c = dot(posShifted,posShifted,1) - g.^2 - (mediumDims(1,n_mediumRow)).^2;
            discriminant = b.^2 - 4.*a.*c;
            discriminant(discriminant <= 0) = NaN;
            distanceMat(n_mediumRow,:) = (-b + (-1).^n_mediumRow.*sqrt(discriminant))./(2.*a);
        case 'plane'
            planeNormal = media.axis(:,medNums(n_mediumRow));
            signs = sign(dot(photons.direction,planeNormal.*ones(1,length(photons.alive))));
            planePosition = media.axisorigin(:,medNums(n_mediumRow)) + (-1).^n_mediumRow.*signs.*media.axis(:,medNums(n_mediumRow)).*media.relevantDimension(1,medNums(n_mediumRow))/2;
            distanceMat(n_mediumRow,:) = dot(planePosition - photons.position,planeNormal.*ones(1,length(photons.alive),classArgument),1)./dot(photons.direction,planeNormal.*ones(1,length(photons.alive),classArgument),1);
    end
end
distanceMat(distanceMat <= options.toleranceLevel) = NaN;

[distances, indices] = min(distanceMat,[],"omitnan");
newPos = photons.position + distances.*photons.direction;
switch sim.shape
    case 'cylinder'
        crossingPhotons = photons.pathLengths >= distances & abs(newPos(2,:)) < sim.halfLength;
    case 'sphere'
        crossingPhotons = photons.pathLengths >= distances & vecnorm(newPos) < sim.radius;
    case 'cube'
        crossingPhotons = photons.pathLengths >= distances & abs(newPos(1,:)) < sim.halfWidth & abs(newPos(2,:)) < sim.halfLength & abs(newPos(3,:)) < sim.halfDepth;
end

if ~isempty(crossingPhotons(crossingPhotons == 1))
    photons.nextMedium(crossingPhotons) = medInd(indices(crossingPhotons));
    if sum(abs(photons.nextMedium(crossingPhotons)-photons.mediumNum(crossingPhotons)))~=sum(crossingPhotons)
        n0 = abs(photons.nextMedium(crossingPhotons)-photons.mediumNum(crossingPhotons))==0;
        n2 = abs(photons.nextMedium(crossingPhotons)-photons.mediumNum(crossingPhotons))==2;
        photons.alive(n0|n2) = 0;
        crossingPhotons(n0|n2) = 0;       
    end
    if length(distances(crossingPhotons)) ~= size(photons.direction(:,crossingPhotons),2)
        disp(['mismatch d:' num2str(length(distances(crossingPhotons))) ', k:' num2str(size(photons.direction(:,crossingPhotons),2))])
    end
    if size(photons.position(:,crossingPhotons),2) ~= size(photons.direction(:,crossingPhotons),2)
        disp(['mismatch r:' num2str(size(photons.position(:,crossingPhotons),2)) ', k:' num2str(size(photons.direction(:,crossingPhotons),2))])
    end
    photons.position(:,crossingPhotons) = newPos(:,crossingPhotons);
    photons.distanceTravelled(crossingPhotons) = photons.distanceTravelled(crossingPhotons) + distances(crossingPhotons);
    photons.opticalPathLength(crossingPhotons) = photons.opticalPathLength(crossingPhotons) + distances(crossingPhotons).*media.refractiveIndex(photons.mediumNum(crossingPhotons));
    photons.pathLengths(crossingPhotons) = 0;
    photons.amplitude(crossingPhotons) = photons.amplitude(crossingPhotons).*exp(horzcat(media.absorptionCoeff(photons.mediumNum(crossingPhotons))).*-distances(crossingPhotons));
end

end