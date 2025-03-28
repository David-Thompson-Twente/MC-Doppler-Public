%%  testPhotonLaunch.m
%   David Thompson, 19-03-2024
%   Tests accuracy of position, initial propagation direction and polarisation for photonLaunch with different sources:
%   makePointSource, makePencilSource, makeFocussedSource, makeRatioSource, makeTruncatedRatioSource
%
%   Inputs:     struct fail_pass carrying any information about previously tested functions
%               classArgument string to determine class to be used, usually double, can also be string or gpuArray
%
%   Outputs:    FAIL or PASS per function
%%
function fail_pass = testPhotonLaunch(fail_pass,classArgument)
arguments
    fail_pass       struct
    classArgument   string  {mustBeMember(classArgument,{'single','double','gpuArray'})}
end
testedFunctions = {'photonLaunchPointSource','photonLaunchPencilSource','photonLaunchFocussedSource','photonLaunchRatioSource','photonLaunchTruncatedRatioSource','photonLaunchPolarisation'};
fail_pass.photonLaunchPointSource = testPhotonLaunchPointSource(classArgument);
fail_pass.photonLaunchPencilSource = testPhotonLaunchPencilSource(classArgument);
fail_pass.photonLaunchFocussedSource = testPhotonLaunchFocussedSource(classArgument);
fail_pass.photonLaunchRatioSource = testPhotonLaunchRatioSource(classArgument);
fail_pass.photonLaunchTruncatedRatioSource = testPhotonTruncatedLaunchRatioSource(classArgument);
fail_pass.photonLaunchPolarisation = testPhotonLaunchPolarisation(classArgument);
passCounter = 0;
fail_pass.photonLaunch = '';
for n = 1:length(testedFunctions)
    if strcmp(fail_pass.(testedFunctions{n}),'PASS')
        fail_pass = rmfield(fail_pass,testedFunctions{n});
        passCounter = passCounter + 1;
    else
        fail_pass.photonLaunch = [fail_pass.photonLaunch 9 testedFunctions{n} '_' fail_pass.(testedFunctions{n})];
        fail_pass = rmfield(fail_pass,testedFunctions{n});
    end
end
if passCounter == length(testedFunctions)
    fail_pass.photonLaunch = 'PASS_ALL';
end
end
%% Unit test for position and direction, photonLaunch with pointSource
function fail_pass = testPhotonLaunchPointSource(classArgument,options)
arguments
    classArgument   string  {mustBeMember(classArgument,{'single','double','gpuArray'})}
    options.sourcePosition      (3,1) {mustBeNumeric}                   = rand(3,1)
    options.wavelength          (1,1) {mustBeNumeric, mustBePositive}   = 666e-9
    options.numPhotons          (1,1) {mustBeInteger, mustBePositive}   = 1e6
    options.killThreshold       (1,1) {mustBeNumeric, mustBePositive}   = 0.1
    options.positionTolerance   (1,1) {mustBeNumeric, mustBePositive}   = 1e3*eps(1)
    options.directionTolerance  (1,1) {mustBeNumeric, mustBePositive}   = .005
end
source = makePointSource(options.sourcePosition,options.wavelength);
sim.numPhotons = options.numPhotons;
photons = photonLaunch(source,sim,options.killThreshold,classArgument);
%Check if all photons start in the same position
positionCheck = photons.position == options.sourcePosition.*ones(1,sim.numPhotons);
%Check if all sides of a unit cube with the point source at the centre receive the same number of photons within 0.5%
zDot = dot(photons.direction,[0;0;1].*ones(1,sim.numPhotons));
xDot = dot(photons.direction,[1;0;0].*ones(1,sim.numPhotons));
yDot = dot(photons.direction,[0;1;0].*ones(1,sim.numPhotons));
positiveZ = zDot(zDot > 0 & abs(zDot) < sqrt(2/3));
negativeZ = zDot(zDot < 0 & abs(zDot) < sqrt(2/3));
positiveX = xDot(xDot > 0 & abs(xDot) < sqrt(2/3));
negativeX = xDot(xDot < 0 & abs(xDot) < sqrt(2/3));
positiveY = yDot(yDot > 0 & abs(yDot) < sqrt(2/3));
negativeY = yDot(yDot < 0 & abs(yDot) < sqrt(2/3));
directionCheck = [length(positiveZ),length(negativeZ),length(positiveX),length(negativeX),length(positiveY),length(negativeY)];
directionCheck = (max(directionCheck) - min(directionCheck))/options.numPhotons <= options.directionTolerance;

fail_pass = 'FAIL';
if min(positionCheck) == 0
    fail_pass = [fail_pass '_POSITION'];
end
if directionCheck == 0
    fail_pass = [fail_pass '_DIRECTION'];
end
if length(fail_pass) > 4
    return
end
fail_pass = 'PASS';
end
%% Unit test for position and direction, photonLaunch with makePencilSource
function fail_pass = testPhotonLaunchPencilSource(classArgument,options)
arguments
    classArgument               string  {mustBeMember(classArgument,{'single','double','gpuArray'})}
    options.sourcePosition      (3,1) {mustBeNumeric}                   = rand(3,1)
    options.focusPosition       (3,1) {mustBeNumeric}                   = rand(3,1)
    options.wavelength          (1,1) {mustBeNumeric, mustBePositive}   = 666e-9
    options.numPhotons          (1,1) {mustBeInteger, mustBePositive}   = 1e6
    options.killThreshold       (1,1) {mustBeNumeric, mustBePositive}   = 0.1
    options.tolerance           (1,1) {mustBeNumeric, mustBePositive}   = 1e3.*eps(1)
end
source = makePencilSource(options.sourcePosition,options.focusPosition,options.wavelength);
sim.numPhotons = options.numPhotons;
photons = photonLaunch(source,sim,options.killThreshold,classArgument);
%Check if all photons start in the same position
positionCheck = photons.position - options.sourcePosition.*ones(1,sim.numPhotons) < options.tolerance;
%Check if all photons propagate in the same direction
directionCheck = dot(photons.direction,source.opticalAxis.*ones(1,options.numPhotons)) - 1 < options.tolerance;
fail_pass = 'FAIL';
if min(positionCheck) == 0
    fail_pass = [fail_pass '_POSITION'];
end
if min(directionCheck) == 0
    fail_pass = [fail_pass '_DIRECTION'];
end
if length(fail_pass) > 4
    return
end
fail_pass = 'PASS';
end
%% Unit test for position and direction, photonLaunch with makeFocussedSource
function fail_pass = testPhotonLaunchFocussedSource(classArgument,options)
arguments
    classArgument               string  {mustBeMember(classArgument,{'single','double','gpuArray'})}
    options.sourcePosition      (3,1) {mustBeNumeric}                   = rand(3,1)
    options.focusPosition       (3,1) {mustBeNumeric}                   = rand(3,1)
    options.focalLength         (1,1) {mustBeNumeric}                   = 50e-3*rand - 25e-3
    options.pupilDiameter       (1,1) {mustBeNumeric, mustBePositive}   = 19.9e-3*rand + 0.1e-3
    options.wavelength          (1,1) {mustBeNumeric, mustBePositive}   = 666e-9
    options.numPhotons          (1,1) {mustBeInteger, mustBePositive}   = 1e6
    options.killThreshold       (1,1) {mustBeNumeric, mustBePositive}   = 0.1
    options.refractiveIndex     (1,1) {mustBeNumeric, mustBePositive}   = 1.42
    options.tolerance           (1,1) {mustBeNumeric, mustBePositive}   = 0.005
end
source = makeFocussedSource(options.sourcePosition,options.refractiveIndex,options.focalLength,options.focusPosition,options.pupilDiameter,options.wavelength);
sim.numPhotons = options.numPhotons;
photons = photonLaunch(source,sim,options.killThreshold,classArgument);
NA = options.refractiveIndex*options.pupilDiameter/(2*sqrt(options.focalLength^2 + (options.pupilDiameter/2)^2));
beamWaist = options.wavelength./(pi*NA);
RayleighRange = pi*beamWaist^2*options.refractiveIndex/options.wavelength;
distanceFromFocus = vecnorm(options.sourcePosition - options.focusPosition);
spotRadius = beamWaist*sqrt(1 + (distanceFromFocus/RayleighRange)^2);
%Check position precision as fraction of spot radius from mean of photon positions
positionCheck = vecnorm(options.sourcePosition - mean(photons.position,2))/spotRadius < options.tolerance;
sourcePlaneVec1 = rotateRodrigues(source.opticalAxis,pi/2,0,[],classArgument);
sourcePlaneVec2 = rotateRodrigues(source.opticalAxis,pi/2,pi/2,sourcePlaneVec1,classArgument);
radiusVector = options.sourcePosition - photons.position;
%Check distribution of photons by looking at std of position - spotRadius/2 as fraction of spotRadius along 2 axes
distributionCheck = (std(dot(radiusVector,sourcePlaneVec1*ones(1,sim.numPhotons))) - spotRadius/2)/spotRadius < options.tolerance...
    && (std(dot(radiusVector,sourcePlaneVec2*ones(1,sim.numPhotons))) - spotRadius/2)/spotRadius < options.tolerance;
%Check directionality of photons against expected directionality from geometric calculation
testAngle = acos(dot(photons.direction,source.opticalAxis.*ones(1,sim.numPhotons)));
referenceAngle = atan(vecnorm(source.position - photons.position)./vecnorm(options.focusPosition - options.sourcePosition));

fail_pass = 'FAIL';
if positionCheck == 0
    fail_pass = [fail_pass '_POSITION'];
end
if distributionCheck == 0
    fail_pass = [fail_pass '_DISTRIBUTION'];
end
if min(abs(testAngle - referenceAngle) < options.tolerance) == 0
    fail_pass = [fail_pass '_DIRECTION'];
end
if length(fail_pass) > 4
    return
end
fail_pass = 'PASS';
end
%% Unit test for position and direction, photonLaunch with makeRatioSource
function fail_pass = testPhotonLaunchRatioSource(classArgument,options)
arguments
    classArgument               string  {mustBeMember(classArgument,{'single','double','gpuArray'})}
    options.sourcePosition      (3,1) {mustBeNumeric}                   = rand(3,1)
    options.focusPosition       (3,1) {mustBeNumeric}                   = rand(3,1)
    options.sourceRadius        (1,1) {mustBeNumeric, mustBePositive}   = rand
    options.focusRadius         (1,1) {mustBeNumeric, mustBePositive}   = rand
    options.wavelength          (1,1) {mustBeNumeric, mustBePositive}   = 666e-9
    options.numPhotons          (1,1) {mustBeInteger, mustBePositive}   = 1e6
    options.killThreshold       (1,1) {mustBeNumeric, mustBePositive}   = 0.1
    options.tolerance           (1,1) {mustBeNumeric, mustBePositive}   = 0.005
end
source = makeRatioSource(options.sourcePosition, options.focusPosition, options.sourceRadius, options.focusRadius, options.wavelength);
sim.numPhotons = options.numPhotons;
photons = photonLaunch(source,sim,options.killThreshold,classArgument);
positionCheck = vecnorm(options.sourcePosition - mean(photons.position,2))/options.sourceRadius < options.tolerance;
sourcePlaneVec1 = rotateRodrigues(source.opticalAxis,pi/2,0,[],classArgument);
sourcePlaneVec2 = rotateRodrigues(source.opticalAxis,pi/2,pi/2,sourcePlaneVec1,classArgument);
radiusVector = options.sourcePosition - photons.position;
distributionCheck = (std(dot(radiusVector,sourcePlaneVec1*ones(1,sim.numPhotons))) - options.sourceRadius)/options.sourceRadius < options.tolerance...
    && (std(dot(radiusVector,sourcePlaneVec2*ones(1,sim.numPhotons))) - options.sourceRadius)/options.sourceRadius < options.tolerance;

testAngle = acos(dot(photons.direction,source.opticalAxis.*ones(1,sim.numPhotons)));
r = vecnorm(source.position - photons.position);
radialFactor = r./options.sourceRadius;
referenceAngle = atan(abs(r - radialFactor.*options.focusRadius)./vecnorm(options.focusPosition - options.sourcePosition));

fail_pass = 'FAIL';
if positionCheck == 0
    fail_pass = [fail_pass '_POSITION'];
end
if distributionCheck == 0
    fail_pass = [fail_pass '_DISTRIBUTION'];
end
if min(abs(testAngle - referenceAngle) < options.tolerance) == 0
    fail_pass = [fail_pass 'DIRECTION'];
end
if length(fail_pass) > 4
    return
end
fail_pass = 'PASS';
end
%% Unit test for position and direction, photonLaunch with makeTruncatedRatioSource
function fail_pass = testPhotonTruncatedLaunchRatioSource(classArgument,options)
arguments
    classArgument               string  {mustBeMember(classArgument,{'single','double','gpuArray'})}
    options.sourcePosition      (3,1) {mustBeNumeric}                   = rand(3,1)
    options.focusPosition       (3,1) {mustBeNumeric}                   = rand(3,1)
    options.sourceRadius        (1,1) {mustBeNumeric, mustBePositive}   = rand
    options.focusRadius         (1,1) {mustBeNumeric, mustBePositive}   = rand
    options.wavelength          (1,1) {mustBeNumeric, mustBePositive}   = 666e-9
    options.numPhotons          (1,1) {mustBeInteger, mustBePositive}   = 1e6
    options.killThreshold       (1,1) {mustBeNumeric, mustBePositive}   = 0.1
    options.tolerance           (1,1) {mustBeNumeric, mustBePositive}   = 0.005
end
source = makeTruncatedRatioSource(options.sourcePosition, options.focusPosition, options.sourceRadius, options.focusRadius, options.wavelength);
sim.numPhotons = options.numPhotons;
photons = photonLaunch(source,sim,options.killThreshold,classArgument);
positionCheck = vecnorm(options.sourcePosition - mean(photons.position,2))/options.sourceRadius < options.tolerance;
sourcePlaneVec1 = rotateRodrigues(source.opticalAxis,pi/2,0,[],classArgument);
sourcePlaneVec2 = rotateRodrigues(source.opticalAxis,pi/2,pi/2,sourcePlaneVec1,classArgument);
radiusVector = options.sourcePosition - photons.position;
distributionCheck = (std(dot(radiusVector,sourcePlaneVec1*ones(1,sim.numPhotons))) - options.sourceRadius)/options.sourceRadius < options.tolerance...
    && (std(dot(radiusVector,sourcePlaneVec2*ones(1,sim.numPhotons))) - options.sourceRadius)/options.sourceRadius < options.tolerance;

testAngle = acos(dot(photons.direction,source.opticalAxis.*ones(1,sim.numPhotons)));
r = vecnorm(source.position - photons.position);
radialFactor = r./options.sourceRadius;
referenceAngle = atan(abs(r - radialFactor.*options.focusRadius)./vecnorm(options.focusPosition - options.sourcePosition));

fail_pass = 'FAIL';
if positionCheck == 0
    fail_pass = [fail_pass '_POSITION'];
end
if distributionCheck == 0
    fail_pass = [fail_pass '_DISTRIBUTION'];
end
if min(abs(testAngle - referenceAngle) < options.tolerance) == 0
    fail_pass = [fail_pass 'DIRECTION'];
end
if length(fail_pass) > 4
    return
end
fail_pass = 'PASS';
end
%% Unit test for polarisation vector, photonLaunch with all source types
function fail_pass = testPhotonLaunchPolarisation(classArgument,options)
arguments
    classArgument               string  {mustBeMember(classArgument,{'single','double','gpuArray'})}
    options.sourcePosition      (3,1) {mustBeNumeric}                   = rand(3,1)
    options.focusPosition       (3,1) {mustBeNumeric}                   = rand(3,1)
    options.sourceRadius        (1,1) {mustBeNumeric, mustBePositive}   = rand
    options.focusRadius         (1,1) {mustBeNumeric, mustBePositive}   = rand
    options.wavelength          (1,1) {mustBeNumeric, mustBePositive}   = 666e-9
    options.numPhotons          (1,1) {mustBeInteger, mustBePositive}   = 1e6
    options.killThreshold       (1,1) {mustBeNumeric, mustBePositive}   = 0.1
    options.focalLength         (1,1) {mustBeNumeric}                   = 50e-3*rand - 25e-3
    options.pupilDiameter       (1,1) {mustBeNumeric, mustBePositive}   = 19.9e-3*rand + 0.1e-3
    options.refractiveIndex     (1,1) {mustBeNumeric, mustBePositive}   = 1.42
    options.polarisationAngle   (1,1) {mustBeNumeric}                   = 2*pi*rand
    options.tolerance           (1,1) {mustBeNumeric, mustBePositive}   = eps(10)
end
sim.numPhotons = options.numPhotons;
fail_pass = 'FAIL';
for n = 1:4
    if n == 1
        source = makePointSource(options.sourcePosition,options.wavelength);
        failMessage = '_POINT_SOURCE';
    elseif n == 2
        source = makePencilSource(options.sourcePosition,options.focusPosition,options.wavelength);
        failMessage = '_PENCIL_SOURCE';
    elseif n == 3
        source = makeFocussedSource(options.sourcePosition,options.refractiveIndex,options.focalLength,options.focusPosition,options.pupilDiameter,options.wavelength);
        failMessage = '_FOCUSSED_SOURCE';
    elseif n == 4
        source = makeRatioSource(options.sourcePosition, options.focusPosition, options.sourceRadius, options.focusRadius, options.wavelength);
        failMessage  = '_RATIO_SOURCE';
    end
    photons = photonLaunch(source,sim,options.killThreshold,classArgument,'polarisationAngle',options.polarisationAngle);
    polarisationReference = cross(photons.direction,[1;0;0]*ones(1,sim.numPhotons));
    polarisationReference = polarisationReference./vecnorm(polarisationReference);
    orthogonalityCheck = abs(dot(photons.polarisationVector,photons.direction)) < options.tolerance;
    polAngleCheck = abs(dot(photons.polarisationVector,polarisationReference) - cos(options.polarisationAngle)) < options.tolerance;
    if min(orthogonalityCheck) == 0
        fail_pass = [fail_pass failMessage '_ORTHOGONALITY'];
    end
    if min(polAngleCheck) == 0
        fail_pass = [fail_pass failMessage '_POLARISATION_ANGLE'];
    end
end
if length(fail_pass) > 20
    return
end
fail_pass = 'PASS';
end