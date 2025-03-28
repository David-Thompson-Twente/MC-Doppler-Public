%%  testDetectors.m
%   David Thompson, 28-03-2024
%   Tests detector functions circularDetector and cylinderSliceDetector
%
%   Inputs:     struct fail_pass carrying any information about previously tested functions
%               media struct with all required media properties
%               classArgument string to determine class to be used, usually double, can also be string or gpuArray
%
%   Outputs:    PASS_ALL if all subtests succesfull
%               FAIL message indicating which part failed
%%
function fail_pass = testDetectors(fail_pass,media,classArgument,options)
arguments
    fail_pass           struct
    media               struct
    classArgument       string  {mustBeMember(classArgument,{'single','double','gpuArray'})}
    options.numPhotons  (1,1)   {mustBeNumeric, mustBePositive}     = 1e6
end
sim.halfLength = media.relevantDimension(2,1);
sim.radius = media.relevantDimension(1,1);
sim.numPhotons = options.numPhotons;
photonDirection = rand(3,options.numPhotons)-0.5;
photons.direction = photonDirection./vecnorm(photonDirection);
photons.position = sim.radius.*rand(1,options.numPhotons).*photons.direction;% + 2.*rand.*photons.direction.*media.relevantDimension(1,1);
photons.amplitude = ones(1,options.numPhotons);

testedFunctions = {'circularDetector','cylinderSliceDetector'};
fail_pass.circularDetector = testCircularDetector(photons,sim,media,classArgument);
fail_pass.cylinderSliceDetector = testCylinderSliceDetector(photons,sim);
passCounter = 0;
fail_pass.testDetectors = '';
for n = 1:length(testedFunctions)
    if strcmp(fail_pass.(testedFunctions{n}),'PASS')
        fail_pass = rmfield(fail_pass,testedFunctions{n});
        passCounter = passCounter + 1;
    else
        fail_pass.testDetectors = [fail_pass.testDetectors 9 testedFunctions{n} '_' fail_pass.(testedFunctions{n})];
        fail_pass = rmfield(fail_pass,testedFunctions{n});
    end
end
if passCounter == length(testedFunctions)
    fail_pass.testDetectors = 'PASS_ALL';
end
end
%% Unit test for circularDetector
function fail_pass = testCircularDetector(photons,sim,media,classArgument,options)
arguments
    photons                 struct
    sim                     struct
    media                   struct  
    classArgument           string  {mustBeMember(classArgument,{'single','double','gpuArray'})}
    options.detectorNA      (1,1)   {mustBeNumeric,mustBePositive}  = 0.1
    options.showPlots       (1,1)   {mustBeNumericOrLogical}        = 0
    options.tolerance       (1,1)   {mustBeNumeric,mustBePositive}  = 1e3*eps(1)
end
detectorPosition = [0;0;mean(media.relevantDimension(1,[2,3]))];
detectorNormal = rotateRodrigues([0;0;-1],pi/8.*rand,2*pi*rand,[],classArgument);
detectorRadius = 0.1*media.relevantDimension(1,1);

detectedPhotons = circularDetector(photons,sim,media,detectorPosition,detectorRadius,'detectorNormal',detectorNormal);

onDetectorTest = vecnorm(detectedPhotons.position - detectorPosition) <= detectorRadius;
distancesToDetectorPlane = dot((detectorPosition.*ones(1,sim.numPhotons) - photons.position),detectorNormal.*ones(1,sim.numPhotons))./dot(-photons.direction,detectorNormal.*ones(1,sim.numPhotons));
photons.position = photons.position - photons.direction.*distancesToDetectorPlane;
onDetectorNumbers = vecnorm(photons.position - detectorPosition) <= detectorRadius & dot(-photons.direction,detectorNormal.*ones(1,sim.numPhotons)) >= 0;

acceptanceAngle = asin(options.detectorNA);
incidenceAngles = acos(dot(-detectedPhotons.direction,detectorNormal.*ones(1,size(detectedPhotons.position,2))));
withinNA = dot(-detectedPhotons.direction,detectorNormal.*ones(1,length(detectedPhotons.position))) > 0 & incidenceAngles <= acceptanceAngle;
detectedPhotonsNA = circularDetector(photons,sim,media,detectorPosition,detectorRadius,'detectorNormal',detectorNormal,'detectorNA',options.detectorNA);
NAtest = sum(withinNA & vecnorm(detectedPhotons.position - detectorPosition) <= detectorRadius) == size(detectedPhotonsNA.position,2);

fail_pass = 'FAIL';
if min(onDetectorTest) == 0
    fail_pass = [fail_pass '_UNDETECTED_PHOTONS_PASSED'];
end
if sum(onDetectorNumbers) ~= sum(onDetectorTest)
    fail_pass = [fail_pass '_PHOTONS_MISSED'];
end
if ~NAtest
    fail_pass = [fail_pass '_NA'];
end
if options.showPlots
    figure;
    plot3(detectedPhotons.position(1,:),detectedPhotons.position(2,:),detectedPhotons.position(3,:),'.')
    hold on
    scale = 0.025;
    quiver3(detectedPhotons.position(1,:),detectedPhotons.position(2,:),detectedPhotons.position(3,:),scale*detectedPhotons.direction(1,:),scale*detectedPhotons.direction(2,:),scale*detectedPhotons.direction(3,:),'AutoScale','off')
    quiver3(detectedPhotonsNA.position(1,:),detectedPhotonsNA.position(2,:),detectedPhotonsNA.position(3,:),scale*detectedPhotonsNA.direction(1,:),scale*detectedPhotonsNA.direction(2,:),scale*detectedPhotonsNA.direction(3,:),'AutoScale','off')
    plot3(detectedPhotonsNA.position(1,:),detectedPhotonsNA.position(2,:),detectedPhotonsNA.position(3,:),'x')
    plot3(detectorPosition(1),detectorPosition(2),detectorPosition(3),'*');
    detectorEdge = detectorPosition + detectorRadius.*rotateRodrigues(detectorNormal.*ones(1,1e3),pi/2.*ones(1,1e3),linspace(0,2*pi,1e3),[],classArgument);
    plot3(detectorEdge(1,:),detectorEdge(2,:),detectorEdge(3,:),'LineWidth',2)
    detectorNormal = detectorNormal.*detectorRadius;
    quiver3(detectorPosition(1),detectorPosition(2),detectorPosition(3),detectorNormal(1),detectorNormal(2),detectorNormal(3),'AutoScale','off')
    xlabel('x (m)');ylabel('y (m)');zlabel('z(m)')
    daspect([1 1 1])
end
if length(fail_pass) > 4
    return
end
fail_pass = 'PASS';
end
%% Unit test for cylinderSliceDetector
function fail_pass = testCylinderSliceDetector(photons,sim,options)
arguments
    photons                         struct
    sim                             struct
    options.showPlots               (1,1)   {mustBeNumericOrLogical}        = 0
    options.tolerance               (1,1)   {mustBeNumeric,mustBePositive}  = 1e3*eps(1)
    options.detectionPlaneNormal    (3,1)   {mustBeNumeric}                 = rand(3,1)
    options.detectorNA              (1,1)   {mustBeNumeric,mustBePositive}  = 0.2 + 0.8.*rand
end
detectionPlaneNormal = options.detectionPlaneNormal./vecnorm(options.detectionPlaneNormal);
sliceRadius = sim.radius.*rand;
sliceThickness = rand*sliceRadius;
sliceOrigin = sim.radius.*rand(3,1);

detectedPhotonsSlice = cylinderSliceDetector(photons,sim,'sliceThickness',sliceThickness,'sliceRadius',sliceRadius,'detectionPlaneNormal',detectionPlaneNormal,'sliceOrigin',sliceOrigin);
detectedPhotonsNA = cylinderSliceDetector(photons,sim,'sliceThickness',sliceThickness,'sliceRadius',sliceRadius,'detectionPlaneNormal',detectionPlaneNormal,'sliceOrigin',sliceOrigin,'detectorNA',options.detectorNA);

normalVectors = ((detectedPhotonsSlice.position - sliceOrigin) - dot(detectedPhotonsSlice.position - sliceOrigin,detectionPlaneNormal.*ones(1,size(detectedPhotonsSlice.position,2))).*detectionPlaneNormal);
normalVectors = -normalVectors./vecnorm(normalVectors);

positionTest = vecnorm(detectedPhotonsSlice.position - sliceOrigin) <= sqrt(sliceRadius.^2 + (sliceThickness/2).^2);
NAtest = sum(acos(abs(dot(normalVectors,detectedPhotonsSlice.direction))) <= asin(options.detectorNA)) == size(detectedPhotonsNA.position,2);

fail_pass = 'FAIL';
if min(positionTest) == 0
    fail_pass = [fail_pass '_POSITION_OUTSIDE_SLICE'];
end
if ~NAtest
    fail_pass = [fail_pass '_NA'];
end

if options.showPlots
    figure;
    plot3(photons.position(1,:),photons.position(2,:),photons.position(3,:),'.')
    daspect([1 1 1])
    title('Original photon positions')
    figure;
    plot3(detectedPhotonsSlice.position(1,:),detectedPhotonsSlice.position(2,:),detectedPhotonsSlice.position(3,:),'.')
    hold on
    plot3(sliceOrigin(1),sliceOrigin(2),sliceOrigin(3),'*')
    plot3(detectedPhotonsNA.position(1,:),detectedPhotonsNA.position(2,:),detectedPhotonsNA.position(3,:),'.')
    circumference = sliceOrigin + sliceRadius.*rotateRodrigues(detectionPlaneNormal,pi/2,linspace(0,2*pi,1e3),[],'double');
    upperEdge = circumference + 0.5*sliceThickness.*detectionPlaneNormal;
    lowerEdge = circumference - 0.5*sliceThickness.*detectionPlaneNormal;
    plot3(circumference(1,:),circumference(2,:),circumference(3,:),'LineWidth',3)
    plot3(lowerEdge(1,:),lowerEdge(2,:),lowerEdge(3,:),'LineWidth',3)
    plot3(upperEdge(1,:),upperEdge(2,:),upperEdge(3,:),'LineWidth',3)
    xlabel('x (m)');ylabel('y (m)');zlabel('z(m)')
    daspect([1 1 1])
    title('Detected photon positions and detector dimensions')
end
if length(fail_pass) > 4
    return
end
fail_pass = 'PASS';
end
