%%  cylinderSliceDetector.m
%   David Thompson 03-11-2022, last update 17-09-2024
%
%   After running a McDoppler sim only extracts those photons which landed in the target slice, a cylinder segment
%   of arbitrary location, orientation and radius, removing the rest
%
%   Inputs:     [struct]        detectedPhotons
%               [struct]        sim
%
%   Optional:   [3x1 double]    detectionPlaneNormal 
%               [1x1 double]    sliceThickness
%               [1x1 double]    sliceRadius
%               [3x1 double]    sliceOrigin
%               [1x1 double]    detectorNA
%
%   Outputs:    [struct]        detectedPhotonsSlice
%%
function detectedPhotonsSlice = cylinderSliceDetector(detectedPhotons, sim, options)
arguments
    detectedPhotons                 struct
    sim                             struct
    options.detectionPlaneNormal    (3,1) {mustBeNumeric}                       = [0;1;0]
    options.sliceThickness          (1,1) {mustBeNumeric, mustBePositive}       = 400e-6
    options.sliceRadius             (1,1) {mustBeNumeric, mustBeNonnegative}    = 0
    options.sliceOrigin             (3,1) {mustBeNumeric}                       = [0;0;0]
    options.detectorNA              (1,1) {mustBeNumeric, mustBePositive}       = 1
    options.toleranceLevel          (1,1) {mustBeNumeric, mustBePositive}       = 10*eps(1)
end
if options.sliceRadius == 0
    options.sliceRadius = sim.radius;
end

options.detectionPlaneNormal = options.detectionPlaneNormal./vecnorm(options.detectionPlaneNormal); % should be normalized

rVec = detectedPhotons.position - options.sliceOrigin;
detVec = options.detectionPlaneNormal.*ones(1,size(detectedPhotons.position,2));
a = 1 - dot(detectedPhotons.direction,detVec,1).^2;
b = 2.*(dot(detectedPhotons.direction,rVec,1) - dot(detectedPhotons.direction,detVec,1).*dot(rVec,detVec,1));
c = dot(rVec,rVec,1) - dot(rVec,detVec,1).^2 - options.sliceRadius.^2;
D = b.^2-4.*a.*c;
distanceMat = (-b + [1;-1].*sqrt(D))./(2*a);
[~,minDistanceIndex] = min(abs(distanceMat));
detectedPhotons.position = detectedPhotons.position + distanceMat(sub2ind(size(distanceMat),minDistanceIndex,1:length(minDistanceIndex))).*detectedPhotons.direction;

normalVectors = ((detectedPhotons.position - options.sliceOrigin) - dot(detectedPhotons.position - options.sliceOrigin,options.detectionPlaneNormal.*ones(1,size(detectedPhotons.position,2))).*options.detectionPlaneNormal);
normalVectors = -normalVectors./vecnorm(normalVectors);

photonsToKeep  = pi/2-acos(abs(dot(options.detectionPlaneNormal.*ones(1,length(detectedPhotons.amplitude)),(detectedPhotons.position-options.sliceOrigin)./vecnorm((detectedPhotons.position-options.sliceOrigin))))) <= atan(options.sliceThickness./(2*options.sliceRadius))...
    & D > 0 & acos(abs(dot(normalVectors,detectedPhotons.direction))) <= asin(options.detectorNA);

detectedPhotonsSlice = emptyDuplicateStruct(detectedPhotons);
fieldNames = fieldnames(detectedPhotons);
for n = 1:length(fieldNames)
    detectedPhotonsSlice.(fieldNames{n}) = detectedPhotons.(fieldNames{n})(:,photonsToKeep);
end
end