%%  stripDetector.m
%   David Thompson 05-09-2024
%
%   Define a strip detector within the simulation medium 
%   create a struct containing any photons that passed through the strip on the way out during their final propagation step
%
%   Inputs:     [struct]        photons
%               [struct]        sim
%               [struct]        media
%               [3x1 double]    detectorPosition
%               [2x1 double]    detectorDimensions
%
%   Optional:   [3x1 double]    detectorNormal, unit normal vector to detector surface, defaults to negative z-axis
%               [1x1 double]    detectorNA, numerical aperture of detector, defining the acceptance angle, defaults to 1
%
%   Outputs:    [struct]        photonsOnDetector

%%
function photonsOnDetector = stripDetector(photons,sim,media,detectorPosition,detectorDimension,options)
arguments
    photons                 struct
    sim                     struct
    media                   struct
    detectorPosition        (3,1)   {mustBeNumeric}
    detectorDimension       (1,2)   {mustBeNumeric, mustBePositive}
    options.detectorNormal  (3,1)   {mustBeNumeric}                     = [0;0;-1]
    options.detectorNA      (1,1)   {mustBeNumeric, mustBeNonnegative}  = 1
    options.omitFields      cell                                        = {};
end
if abs(detectorPosition(2)) > sim.halfLength || vecnorm(detectorPosition([1,3])) > sim.radius
    error('Detector placed outside simulation boundaries')
end
detectorStruct.position = detectorPosition;
detectorMedium = currentMedium(detectorStruct,media);
acceptanceAngle = asin(options.detectorNA);

distancesToDetectorPlane = dot((detectorPosition.*ones(1,length(photons.position)) - photons.position),options.detectorNormal.*ones(1,length(photons.position)))./dot(-photons.direction,options.detectorNormal.*ones(1,length(photons.position)));
projectedPositions = photons.position - photons.direction.*distancesToDetectorPlane;
relativePositions = projectedPositions - detectorPosition;
onDetector = relativePositions(1,:) >= -detectorDimension(1) & relativePositions(1,:) <= detectorDimension(1) & relativePositions(2,:) >= -detectorDimension(2) & relativePositions(2,:) <= detectorDimension(2);

incidenceAngles = acos(dot(-photons.direction,options.detectorNormal.*ones(1,length(photons.position))));
withinNA = dot(-photons.direction,options.detectorNormal.*ones(1,length(photons.position))) > 0 & incidenceAngles <= acceptanceAngle;

detected = onDetector & withinNA;
photons.position(:,detected) = projectedPositions(:,detected);
photons.mediumNum(detected) = horzcat(detectorMedium.mediumNum);
photons.distanceTravelled = photons.distanceTravelled - distancesToDetectorPlane;
photons.opticalPathLength = photons.opticalPathLength - media.refractiveIndex(photons.mediumNum).*distancesToDetectorPlane;
photonsOnDetector = emptyDuplicateStruct(photons,'omitFields',options.omitFields);
fieldNames = fieldnames(photonsOnDetector);
for n = 1:length(fieldNames)
    photonsOnDetector.(fieldNames{n}) = photons.(fieldNames{n})(:,detected);
end
end