%%  updateDirectionsReflectRefract.m
%   David Thompson 12-4-2022, last update 19-03-2025
% 
%   Updates photons direction and medium at an interface, either reflecting or refracting based on Tkaczyk 2012:
%   https://opg.optica.org/ol/fulltext.cfm?uri=ol-37-5-972&id=229633
%
%   Inputs:     [struct]        photons, with relevant properties:      photons.direction
%                                                                       photons.alive
%                                                                       photons.polarisationVector
%                                                                       photons.timesScattered
%                                                                       photons.timesReflected
%                                                                       photons.timesRefracted
%                                                                       photons.mediumNum
%                                                                       photons.nextMedium
%               [3xn double]    normalVectors, giving normal vectors to surface for any interfacing photons
%               [1xn double]    crossingPhotons, indicating which photons actually need to reflect or refract
%               [2xn double]    reflectionCoefficients
%               [1xn double]    indexRatio
%               [string]        polarisation
%               [string]        classArgument
%
%   Optional:   [1x1 double]    tolerance
%
%   Outputs:    [struct]        photons, with updated properties        photons.direction
%                                                                       photons.alive
%                                                                       photons.polarisationVector
%                                                                       photons.timesScattered
%                                                                       photons.timesReflected
%                                                                       photons.timesRefracted
%                                                                       photons.mediumNum
%                                                                       photons.nextMedium

%%
function photons = updateDirectionsReflectRefract(photons,normalVectors,crossingPhotons,reflectionCoefficients,indexRatio,polarisation,classArgument,options)
arguments
    photons                 struct
    normalVectors           (3,:)   {mustBeNumeric(normalVectors)}
    crossingPhotons         (1,:)   {mustBeNumericOrLogical(crossingPhotons),mustBeInteger(crossingPhotons),mustBeLessThanOrEqual(crossingPhotons,1),mustBeNonnegative(crossingPhotons)}
    reflectionCoefficients  (2,:)   {mustBeNumeric(reflectionCoefficients),mustBeLessThanOrEqual(reflectionCoefficients,1.00001),mustBeNonnegative(reflectionCoefficients)}
    indexRatio              (1,:)   {mustBeNumeric(indexRatio),mustBePositive(indexRatio)}
    polarisation            string  {mustBeMember(polarisation,{'unpolarised','polarised'})}
    classArgument           string
    options.tolerance       (1,1)   {mustBeNumeric} = 20*eps(1)
end
normalCheck = abs(abs(dot(photons.direction,normalVectors,1)) - 1) < options.tolerance;
if ~isempty(normalCheck(normalCheck == 1))  % If photons.direction is perfectly normal to the interface, a minor shift is applied
    parallelVec(:,normalCheck) = cross(photons.direction(:,normalCheck) + options.tolerance.*rand(3,1,classArgument),normalVectors(:,normalCheck),1);
    parallelVec(:,~normalCheck) = cross(photons.direction(:,~normalCheck),normalVectors(:,~normalCheck),1);
else
    parallelVec = cross(photons.direction,normalVectors,1);
end
numCrossing = length(crossingPhotons(crossingPhotons == 1));
parallelVec(:,crossingPhotons) = parallelVec(:,crossingPhotons)./vecnorm(parallelVec(:,crossingPhotons));
perpendicularIn = zeros(size(parallelVec),classArgument);
perpendicularIn(:,crossingPhotons) = rotateRodrigues(photons.direction(:,crossingPhotons),pi/2.*ones(1,numCrossing,classArgument),pi/2.*ones(1,numCrossing,classArgument),parallelVec(:,crossingPhotons),classArgument);%
polarisationComponents = [dot(photons.polarisationVector,parallelVec,1);dot(photons.polarisationVector,perpendicularIn,1)];

switch polarisation
    case 'unpolarised'
        R = (reflectionCoefficients(1,:)+reflectionCoefficients(2,:))/2;
    case 'polarised'
        R = (abs(polarisationComponents(1,:)).*reflectionCoefficients(1,:)+abs(polarisationComponents(2,:)).*reflectionCoefficients(2,:))./sum(abs(polarisationComponents));
end

randomNumbers = rand(1,length(R),classArgument);

%reflection
reflectedPhotons = randomNumbers <= R & crossingPhotons & indexRatio~=1;
if sum(reflectedPhotons) > 0
    photons.direction(:,reflectedPhotons) = photons.direction(:,reflectedPhotons) - 2*dot(photons.direction(:,reflectedPhotons),normalVectors(:,reflectedPhotons),1).*normalVectors(:,reflectedPhotons);
    perpendicularReflected(:,reflectedPhotons) = rotateRodrigues(photons.direction(:,reflectedPhotons),pi/2.*ones(1,sum(reflectedPhotons),classArgument),pi/2.*ones(1,sum(reflectedPhotons),classArgument),parallelVec(:,reflectedPhotons),classArgument);
    newPolarisationComponents(:,reflectedPhotons) = [-reflectionCoefficients(1,reflectedPhotons).*polarisationComponents(1,reflectedPhotons);-reflectionCoefficients(2,reflectedPhotons).*polarisationComponents(2,reflectedPhotons)];
    photons.polarisationVector(:,reflectedPhotons) = newPolarisationComponents(1,reflectedPhotons).*parallelVec(:,reflectedPhotons) + newPolarisationComponents(2,reflectedPhotons).*perpendicularReflected(:,reflectedPhotons);
end
%refraction
refractedPhotons = randomNumbers > R & crossingPhotons & indexRatio~=1;
if sum(refractedPhotons) > 0
    photons.direction(:,refractedPhotons) = photons.direction(:,refractedPhotons)./indexRatio(refractedPhotons) + (dot(-normalVectors(:,refractedPhotons),photons.direction(:,refractedPhotons),1)./indexRatio(refractedPhotons) - sqrt(1-(1./indexRatio(refractedPhotons)).^2.*(1-dot(-normalVectors(:,refractedPhotons),photons.direction(:,refractedPhotons),1).^2))).*normalVectors(:,refractedPhotons);
    perpendicularRefracted(:,refractedPhotons) = rotateRodrigues(photons.direction(:,refractedPhotons),pi/2.*ones(1,sum(refractedPhotons),classArgument),pi/2.*ones(1,sum(refractedPhotons),classArgument),parallelVec(:,refractedPhotons),classArgument);
    newPolarisationComponents(:,refractedPhotons) = [(1-reflectionCoefficients(1,refractedPhotons)).*polarisationComponents(1,refractedPhotons);(1-reflectionCoefficients(2,refractedPhotons)).*polarisationComponents(2,refractedPhotons)];
    photons.polarisationVector(:,refractedPhotons) = newPolarisationComponents(1,refractedPhotons).*parallelVec(:,refractedPhotons) + newPolarisationComponents(2,refractedPhotons).*perpendicularRefracted(:,refractedPhotons);
end
photons.polarisationVector(:,crossingPhotons) = photons.polarisationVector(:,crossingPhotons)./vecnorm(photons.polarisationVector(:,crossingPhotons));

photons.mediumNum(randomNumbers > R & crossingPhotons) = photons.nextMedium(randomNumbers > R & crossingPhotons);
photons.nextMedium(randomNumbers <= R & crossingPhotons) = photons.mediumNum(randomNumbers <= R & crossingPhotons);
%Update reflected/refracted counters
photons.timesReflected(randomNumbers <= R & crossingPhotons) = photons.timesReflected(randomNumbers <= R & crossingPhotons) + 1;
photons.timesRefracted(randomNumbers > R & crossingPhotons) = photons.timesRefracted(randomNumbers > R & crossingPhotons) + 1;
end