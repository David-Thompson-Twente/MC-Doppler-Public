%%  testUpdateDirectionsScatter.m
%   David Thompson 26-03-2024, last updat 19-03-2025
%   Tests the following outputs of updateDirectionsScatter:
%   photons.direction, photons.timesScattered, photons.DopplerShift
%
%   Inputs:     struct fail_pass carrying any information about previously tested functions
%               media struct with all required media properties
%               classArgument string to determine class to be used, usually double, can also be string or gpuArray
%
%   Outputs:    PASS_ALL if all subtests succesfull
%               FAIL message indicating which part failed
%%
function fail_pass = testUpdateDirectionsScatter(fail_pass,media,classArgument,options)
arguments
    fail_pass                   struct
    media                       struct
    classArgument               string  {mustBeMember(classArgument,{'single','double','gpuArray'})}
    options.numPhotons          (1,1)   {mustBeNumeric, mustBePositive}     = 1e6
    options.photonsDirection    (3,1)   {mustBeNumeric}                     = rand(3,1)
    options.wavelength          (1,1)   {mustBeNumeric, mustBePositive}     = 666e-9
end
photons.direction = (options.photonsDirection./vecnorm(options.photonsDirection)).*ones(1,options.numPhotons);
photons.polarisationVector = rotateRodrigues(photons.direction,pi/2,0,[],classArgument);
photons.position = zeros(3,options.numPhotons);
photons.wavelength = options.wavelength;
photons.DopplerShift = zeros(1,options.numPhotons);
photons.timesScattered = zeros(1,options.numPhotons);
photons.likelihood = zeros(1,options.numPhotons);
photons.alive = ones(1,options.numPhotons);
scatteringMediumNumber = find(media.scatterCoeff ~= 0);
media.S1(scatteringMediumNumber) = csaps(linspace(0,pi,length(media.phaseFunction{scatteringMediumNumber}(1,:))),media.phaseFunction{scatteringMediumNumber}(1,:));
media.S2(scatteringMediumNumber) = csaps(linspace(0,pi,length(media.phaseFunction{scatteringMediumNumber}(2,:))),media.phaseFunction{scatteringMediumNumber}(2,:));

testedFunctions = {'scatterDirectionsCounter','scatterPolarisation','scatterDoppler'};
fail_pass.scatterDirectionsCounter = testUpdateDirectionsScatterDirection(photons,media,classArgument);
fail_pass.scatterPolarisation = testUpdateDirectionsScatterPolarisation(photons,media,classArgument);
fail_pass.scatterDoppler = testUpdateDirectionsScatterDoppler(photons,media,classArgument);
passCounter = 0;
fail_pass.updateDirectionsScatter = '';
for n = 1:length(testedFunctions)
    if strcmp(fail_pass.(testedFunctions{n}),'PASS')
        fail_pass = rmfield(fail_pass,testedFunctions{n});
        passCounter = passCounter + 1;
    else
        fail_pass.updateDirectionsScatter = [fail_pass.updateDirectionsScatter 9 testedFunctions{n} '_' fail_pass.(testedFunctions{n})];
        fail_pass = rmfield(fail_pass,testedFunctions{n});
    end
end
if passCounter == length(testedFunctions)
    fail_pass.updateDirectionsScatter = 'PASS_ALL';
end
end
%% Unit test for directionality and scatter counter of updateDirectionsScatter
function fail_pass = testUpdateDirectionsScatterDirection(photons,media,classArgument,options)
arguments
    photons                     struct
    media                       struct
    classArgument               string  {mustBeMember(classArgument,{'single','double','gpuArray'})}    
    options.showPlots           (1,1)   {mustBeNumericOrLogical}            = 0
    options.tolerance           (1,1)   {mustBeNumeric, mustBePositive}     = 1e3*eps(1)
end
numPhotons = length(photons.DopplerShift);
scatteringPhotons = logical(zeros(1,numPhotons));
scatteringPhotons(rand(1,numPhotons) < 0.5) = logical(1);
scatteringAngles = zeros(2,numPhotons);
numScatter = length(scatteringPhotons(scatteringPhotons == 1));
scatteringAngles(:,logical(scatteringPhotons)) = [linspace(0,4*pi,numScatter);linspace(0,16*pi,numScatter)];
scatteringMediumNumber = find(media.scatterCoeff ~= 0);
photons.mediumNum = find(media.scatterCoeff == 0,1,'first').*ones(1,numPhotons);
photons.mediumNum(scatteringPhotons == 1) = scatteringMediumNumber.*ones(1,length(scatteringPhotons(scatteringPhotons == 1)));
photons.alive = ones(1,numPhotons);
originalDirections = photons.direction;


reference_direction = cos(scatteringAngles(2,:)).*photons.direction(:,1);
reference_polarisation = sin(scatteringAngles(2,:)).*(cos(scatteringAngles(1,:)).*photons.polarisationVector(:,1) + cos(scatteringAngles(1,:) - pi/2).*cross(photons.direction(:,1),photons.polarisationVector(:,1)));
referenceVectors = reference_direction + reference_polarisation;
photons = updateDirectionsScatter(photons,scatteringPhotons,scatteringAngles,media,classArgument);

testAngles = acos(dot(photons.direction,originalDirections));
testCriterion = abs(photons.direction - referenceVectors) < options.tolerance;
counterTest = photons.timesScattered(scatteringPhotons) == ones(1,numScatter);
if options.showPlots
    figure;
    quiver3(0,0,0,originalDirections(1),originalDirections(2),originalDirections(3),'AutoScale','off')
    hold on
    quiver3(originalDirections(1).*ones(1,numScatter),originalDirections(2).*ones(1,numScatter),originalDirections(3).*ones(1,numScatter),photons.direction(1,scatteringPhotons),photons.direction(2,scatteringPhotons),photons.direction(3,scatteringPhotons),'AutoScale','off')
    daspect([1 1 1])
    plot3(referenceVectors(1,scatteringPhotons)+originalDirections(1).*ones(1,numScatter),referenceVectors(2,scatteringPhotons)+originalDirections(2).*ones(1,numScatter),referenceVectors(3,scatteringPhotons)+originalDirections(3).*ones(1,numScatter),'LineWidth',2)
    figure;
    plot(scatteringAngles(2,scatteringPhotons))
    hold on
    plot(testAngles(scatteringPhotons),'--')
end
fail_pass = 'FAIL';
if min(testCriterion(:)) == 0
    fail_pass = [fail_pass '_DIRECTION'];
end
if min(counterTest) == 0
    fail_pass = [fail_pass '_COUNTER'];
end
if length(fail_pass) > 4
    return
end
fail_pass = 'PASS';
end
%% Unit test for polarisation of updateDirectionsScatter
function fail_pass = testUpdateDirectionsScatterPolarisation(photons,media,classArgument,options)
arguments
    photons                     struct
    media                       struct
    classArgument               string  {mustBeMember(classArgument,{'single','double','gpuArray'})}    
    options.showPlots           (1,1)   {mustBeNumericOrLogical}            = 0
    options.tolerance           (1,1)   {mustBeNumeric, mustBePositive}     = 1e3*eps(1)
end
numPhotons = length(photons.DopplerShift);
originalDirections = photons.direction;
scatteringPhotons = ones(1,numPhotons);
photons.mediumNum = find(media.scatterCoeff ~= 0).*ones(1,numPhotons);
scatteringAngles = [linspace(0,2*pi,numPhotons-2), 0, pi/2;linspace(0,pi,numPhotons-2), pi*rand(1,2)];
oldPolarisation = photons.polarisationVector;
rotatedOldInPlane = rotateRodrigues(oldPolarisation(:,end-1),-scatteringAngles(2,end-1),scatteringAngles(1,end-1),originalDirections(:,1),classArgument);
photons = updateDirectionsScatter(photons,scatteringPhotons,scatteringAngles,media,classArgument);
perpendicularityTest = dot(photons.polarisationVector,photons.direction) < options.tolerance;
unitVectorTest = abs(vecnorm(photons.polarisationVector) - 1) < options.tolerance;
inPlaneTest = abs(dot(cross(oldPolarisation(:,end-1),originalDirections(:,1)),cross(photons.polarisationVector(:,end-1),photons.direction(:,end-1))) - 1) < options.tolerance;
outOfPlaneTest = abs(oldPolarisation(:,end) - photons.polarisationVector(:,end)) < options.tolerance;

if options.showPlots
    figure;
    quiver3(zeros(1,numPhotons),zeros(1,numPhotons),zeros(1,numPhotons),photons.direction(1,:),photons.direction(2,:),photons.direction(3,:),'AutoScale','off')
    hold on
    quiver3(zeros(1,numPhotons),zeros(1,numPhotons),zeros(1,numPhotons),photons.polarisationVector(1,:),photons.polarisationVector(2,:),photons.polarisationVector(3,:),'AutoScale','off')
    daspect([1 1 1])
    figure;
    quiver3(0,0,0,originalDirections(1),originalDirections(2),originalDirections(3),'AutoScale','off')
    hold on
    quiver3(originalDirections(1),originalDirections(2),originalDirections(3),photons.direction(1,end-1),photons.direction(2,end-1),photons.direction(3,end-1),'AutoScale','off')
    quiver3(0,0,0,oldPolarisation(1,end-1),oldPolarisation(2,end-1),oldPolarisation(3,end-1),'AutoScale','off')
    quiver3(originalDirections(1),originalDirections(2),originalDirections(3),rotatedOldInPlane(1),rotatedOldInPlane(2),rotatedOldInPlane(3),'AutoScale','off')
    quiver3(originalDirections(1),originalDirections(2),originalDirections(3),photons.polarisationVector(1,end-1),photons.polarisationVector(2,end-1),photons.polarisationVector(3,end-1),'--','AutoScale','off')
    daspect([1 1 1])
    figure;
    quiver3(0,0,0,originalDirections(1),originalDirections(2),originalDirections(3),'AutoScale','off')
    hold on
    quiver3(originalDirections(1),originalDirections(2),originalDirections(3),photons.direction(1,end),photons.direction(2,end),photons.direction(3,end),'AutoScale','off')
    quiver3(0,0,0,oldPolarisation(1,end),oldPolarisation(2,end),oldPolarisation(3,end),'AutoScale','off')
    quiver3(originalDirections(1),originalDirections(2),originalDirections(3),photons.polarisationVector(1,end),photons.polarisationVector(2,end),photons.polarisationVector(3,end),'--','AutoScale','off')
    daspect([1 1 1])
end
fail_pass = 'FAIL';
if min(perpendicularityTest) == 0
    fail_pass = [fail_pass '_NOT_PERPENDICULAR_TO_PROPAGATION'];
end
if min(unitVectorTest) == 0
    fail_pass = [fail_pass '_VECTORS_NOT_UNIT'];
end
if min(inPlaneTest) == 0
    fail_pass = [fail_pass '_IN_PLANE_WRONG'];
end
if min(outOfPlaneTest) == 0
    fail_pass = [fail_pass '_OUT_OF_PLANE_WRONG'];
end
if length(fail_pass) > 4
    return
end
fail_pass = 'PASS';
end
%% Unit test for Doppler Shift of updateDirectionsScatter
function fail_pass = testUpdateDirectionsScatterDoppler(photons,media,classArgument,options)
arguments
    photons                     struct
    media                       struct
    classArgument               string  {mustBeMember(classArgument,{'single','double','gpuArray'})}    
    options.showPlots           (1,1)   {mustBeNumericOrLogical}            = 0
    options.tolerance           (1,1)   {mustBeNumeric, mustBePositive}     = 5e-10
end
numPhotons = length(photons.DopplerShift);
photons.direction = [0;0;1].*ones(1,numPhotons);
photons.polarisationVector = [0;1;0].*ones(1,numPhotons);
originalDirections = photons.direction;
scatteringPhotons = zeros(1,numPhotons);
photons.position = media.relevantDimension(1,1).*rand(1,numPhotons).*[1;0;0] + media.relevantDimension(1,2).*rand(1,numPhotons).*[0;1;0];
elevationAngle = acos(dot(photons.position./vecnorm(photons.position),[0;1;0].*ones(1,numPhotons)));
photons.position = vecnorm(photons.position).*rotateRodrigues([0;1;0].*ones(1,numPhotons),elevationAngle,2*pi*rand(1,numPhotons),[1;0;0],classArgument);
photons = currentMedium(photons,media);
numScatter = sum(photons.mediumNum == find(media.scatterCoeff ~= 0));
scatteringAngles = zeros(2,numPhotons);
scatteringPhotons(photons.mediumNum == find(media.scatterCoeff ~= 0)) = 1;
scatteringAngles(:,scatteringPhotons == 1) = [linspace(0,4*pi,numScatter-2), 0, pi/2;linspace(0,2*pi,numScatter-2), pi*rand(1,2)];
newDirections = rotateRodrigues(photons.direction,scatteringAngles(2,:),scatteringAngles(1,:),photons.polarisationVector,classArgument);
qUnit = (originalDirections - newDirections)./vecnorm(originalDirections - newDirections);
flowVelocity = zeros(3,numPhotons);
for n = 1:length(media.name)
    flowVelocity(:,photons.mediumNum == n) = media.flowVelocity{n}(photons.position(1,photons.mediumNum == n),photons.position(2,photons.mediumNum == n),photons.position(3,photons.mediumNum == n));
end

flowDirectionUnit = flowVelocity./vecnorm(flowVelocity);
flowSpeeds = vecnorm(flowVelocity);
expectedDopplerShift = -2*abs(flowSpeeds).*media.refractiveIndex(media.scatterCoeff ~= 0)./photons.wavelength.*sin(scatteringAngles(2,:)./2).*dot(flowDirectionUnit,qUnit);

expectedDopplerShift(photons.mediumNum ~= find(media.scatterCoeff ~= 0)) = 0;

photons = updateDirectionsScatter(photons,logical(scatteringPhotons),scatteringAngles,media,classArgument);

if options.showPlots
    figure;plot(photons.DopplerShift,'.')
    hold on;plot(expectedDopplerShift,'.')
end
fail_pass = 'FAIL';
if max(abs((expectedDopplerShift - photons.DopplerShift)./photons.DopplerShift)) > options.tolerance
    fail_pass = [fail_pass '_DOPPLER_SHIFT'];
end
if length(fail_pass) > 4
    return
end
fail_pass = 'PASS';
end