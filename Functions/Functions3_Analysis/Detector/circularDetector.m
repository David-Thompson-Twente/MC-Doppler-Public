%%  circularDetector.m
%   David Thompson 17-04-2023, last update 31-10-2024
%
%   Define an arbitrarily oriented disc-shaped detector within the simulation medium, 
%   create a struct containing any photons that passed through the disc on the way out during their final propagation step
%
%   Inputs:     [struct]        photons
%               [struct]        sim
%               [struct]        media
%               [3x1 double]    detectorPosition
%               [1x1 double]    detectorRadius
%
%   Optional:   [3x1 double]    detectorNormal, unit normal vector to detector surface, defaults to negative z-axis
%               [1x1 double]    detectorNA, numerical aperture of detector, defining the acceptance angle, defaults to 1
%
%   Outputs:    [struct]        photonsOnDetector

%%
function photonsOnDetector = circularDetector(photons,sim,media,detectorPosition,detectorRadius,options)
arguments
    photons                 struct
    sim                     struct
    media                   struct
    detectorPosition        (3,1)   {mustBeNumeric}
    detectorRadius          (1,1)   {mustBeNumeric, mustBePositive}
    options.detectorNormal  (3,1)   {mustBeNumeric}                     = [0;0;-1]
    options.detectorNA      (1,1)   {mustBeNumeric, mustBeNonnegative}  = 1
    options.omitFields      cell                                        = {};
end
if abs(detectorPosition(2)) > sim.halfLength || vecnorm(detectorPosition([1,3])) > sim.radius
    error('Detector placed outside simulation boundaries')
end

distancesToDetectorPlane = dot((detectorPosition.*ones(1,length(photons.position)) - photons.position),options.detectorNormal.*ones(1,length(photons.position)))./dot(-photons.direction,options.detectorNormal.*ones(1,length(photons.position)));
projectedPositions = photons.position - photons.direction.*distancesToDetectorPlane;
onDetector = vecnorm(projectedPositions - detectorPosition) <= detectorRadius & photons.amplitude > 0;

detectorStruct.position = detectorPosition;
detectorMedium = currentMedium(detectorStruct,media);
acceptanceAngle = asin(options.detectorNA);
incidenceAngles = acos(dot(-photons.direction,options.detectorNormal.*ones(1,length(photons.position))));
withinNA = dot(-photons.direction,options.detectorNormal.*ones(1,length(photons.position))) > 0 & incidenceAngles <= acceptanceAngle;

detected = onDetector & withinNA;
photons.position(:,detected) = projectedPositions(:,detected);
photons.mediumNum(detected) = horzcat(detectorMedium.mediumNum);

if isfield(photons,'distanceTravelled')
    photons.distanceTravelled = photons.distanceTravelled - distancesToDetectorPlane;
end
if isfield(photons,'opticalPathLength')
    photons.opticalPathLength = photons.opticalPathLength - gather(media.refractiveIndex(photons.mediumNum)).*distancesToDetectorPlane;
end

photonsOnDetector = emptyDuplicateStruct(photons,'omitFields',options.omitFields);
fieldNames = fieldnames(photonsOnDetector);
for n = 1:length(fieldNames)
    photonsOnDetector.(fieldNames{n}) = photons.(fieldNames{n})(:,detected);
end
end