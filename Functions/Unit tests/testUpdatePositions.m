%%  testUpdatePositions.m
%   David Thompson, last update 19-03-2024
%   Tests updatePositions and the following sub-functions:
%   pathLengths, interfaceIntersect, interfaceNormalVectors, photonEscapeCheck, photonEscapeTruncate
%
%   Inputs:     struct fail_pass carrying any information about previously tested functions
%               media struct with all required media properties
%               classArgument string to determine class to be used, usually double, can also be string or gpuArray
%
%   Outputs:    PASS_ALL if all subtests succesfull
%               FAIL message indicating which part failed
%%
function fail_pass = testUpdatePositions(fail_pass,media,classArgument)
arguments
    fail_pass       struct
    media           struct
    classArgument   string  {mustBeMember(classArgument,{'single','double','gpuArray'})}
end
testedFunctions = {'pathLengths','interfaceIntersect','interfaceNormalVectors','photonEscapeCheck','photonEscapeTruncate','unitTestUpdatePositions'};
fail_pass.pathLengths = testPathLengths(media,classArgument);
fail_pass.interfaceIntersect = testInterfaceIntersect(media,classArgument);
fail_pass.interfaceNormalVectors = testInterfaceNormalVectors(media,classArgument);
fail_pass.photonEscapeCheck = testPhotonEscapeCheck(classArgument);
fail_pass.photonEscapeTruncate = testPhotonEscapeTruncate(classArgument);
fail_pass.unitTestUpdatePositions = unitTestUpdatePositions(media,classArgument);
passCounter = 0;
fail_pass.updatePositions = '';
for n = 1:length(testedFunctions)
    if strcmp(fail_pass.(testedFunctions{n}),'PASS')
        fail_pass = rmfield(fail_pass,testedFunctions{n});
        passCounter = passCounter + 1;
    else
        fail_pass.updatePositions = [fail_pass.updatePositions 9 testedFunctions{n} '_' fail_pass.(testedFunctions{n})];
        fail_pass = rmfield(fail_pass,testedFunctions{n});
    end
end
if passCounter == length(testedFunctions)
    fail_pass.updatePositions = 'PASS_ALL';
end
end
%% Unit test for pathLengths
function fail_pass = testPathLengths(media,classArgument,options)
arguments
    media                   struct
    classArgument           string  {mustBeMember(classArgument,{'single','double','gpuArray'})}
    options.numPhotons      (1,1) {mustBeNumeric, mustBePositive} = 1e6
    options.numBins         (1,1) {mustBeNumeric, mustBePositive} = 1e3;
end
scatterCoeff = media.scatterCoeff(media.scatterCoeff ~= 0);
photons.mediumNum = ones(1,options.numPhotons);
photons.pathLengths = zeros(1,options.numPhotons);
photons.alive = ones(1,options.numPhotons);
photons = pathLengths(media,photons,classArgument);
maxPathLength = max(photons.pathLengths);
LambertBeer = exp(-scatterCoeff*linspace(0,maxPathLength,options.numBins));
pathLengthHist = histcounts(photons.pathLengths,options.numBins);
CC = corrcoef(LambertBeer,pathLengthHist);
if min(CC(:)) > 0.99
    fail_pass = 'PASS';
else
    fail_pass = 'FAIL';
end
end
%% Unit test for interfaceIntersect
function fail_pass = testInterfaceIntersect(media,classArgument,options)
arguments
    media                   struct
    classArgument           string  {mustBeMember(classArgument,{'single','double','gpuArray'})}
    options.numPhotons          (1,1) {mustBeNumeric, mustBePositive} = 1e6
    options.tolerance           (1,1) {mustBeNumeric, mustBePositive} = 1e3*eps(1);
end
sim.GPU = 0;
sim.halfLength = media.relevantDimension(2,1);
sim.numPhotons = options.numPhotons;
sim.shape = 'cylinder';
source = makePointSource([0;0;0],666e-9);
photonStartPositions = {[0;0;-media.relevantDimension(1,2)-(media.relevantDimension(1,1) - media.relevantDimension(1,2))/2]*ones(1,options.numPhotons),...
    [0;0;-media.relevantDimension(1,3)-(media.relevantDimension(1,2) - media.relevantDimension(1,3))/2]*ones(1,options.numPhotons), [0;0;0]*ones(1,options.numPhotons)};
startMedium = [1;2;3];
photons.direction = photonDirection(source,sim,[],classArgument);
photons.amplitude = ones(1,options.numPhotons);
fail_pass = 'FAIL';
for n = 1:length(photonStartPositions)
    photons.position = photonStartPositions{n};
    photons.pathLengths = 20*ones(1,options.numPhotons);
    photons.mediumNum = startMedium(n)*ones(1,options.numPhotons);
    photons.alive = ones(1,options.numPhotons);
    photons.distanceTravelled = zeros(1,options.numPhotons);
    photons.opticalPathLength = zeros(1,options.numPhotons);
    photons = interfaceIntersect(photons,media,sim,classArgument);
    switch startMedium(n)
        case 1
            criterion = abs(vecnorm(photons.position([1,3],photons.pathLengths == 0)) - media.relevantDimension(1,2)) < options.tolerance;
            if min(criterion) == 0
                fail_pass = [fail_pass '_MEDIUM_1'];
            end
        case 2
            criterion = abs(vecnorm(photons.position([1,3],photons.pathLengths == 0)) - media.relevantDimension(1,2)) < options.tolerance | abs(vecnorm(photons.position([1,3],photons.pathLengths == 0)) - media.relevantDimension(1,3)) < options.tolerance;
            if min(criterion) == 0
                fail_pass = [fail_pass '_MEDIUM_2'];
            end
        case 3
            criterion = abs(vecnorm(photons.position([1,3],photons.pathLengths == 0)) - media.relevantDimension(1,3)) < options.tolerance;
            if min(criterion) == 0
                fail_pass = [fail_pass '_MEDIUM_3'];
            end
    end
end
if length(fail_pass) > 4
    return
end
fail_pass = 'PASS';
end
%% Unit test for interfaceNormalVectors
function fail_pass = testInterfaceNormalVectors(media,classArgument,options)
arguments
    media               struct
    classArgument       string  {mustBeMember(classArgument,{'single','double','gpuArray'})}
    options.numPhotons  (1,1) {mustBeNumeric, mustBePositive} = 1e6
    options.simRadius   (1,1) {mustBeNumeric, mustBePositive} = rand
end
sim.radius = options.simRadius;
sim.halfLength = options.simRadius;
sim.numPhotons = options.numPhotons;
source = makePointSource([0;0;0],666e-9);
photons.position = [rand(1,options.numPhotons);rand(1,options.numPhotons);rand(1,options.numPhotons)];
photons.direction = photonDirection(source,sim,photons.position,classArgument);
photons.direction(2,1:10) = 0;
photons.direction(:,1:10) = photons.direction(:,1:10)./vecnorm(photons.direction(:,1:10));
photons.direction(:,11) = [1;0;0];
photons.direction(:,12) = [-1;0;0];
photons.direction(:,13) = [0;1;0];
photons.direction(:,14) = [0;-1;0];
photons.direction(:,15) = [0;0;1];
photons.direction(:,16) = [0;0;-1];
photons.mediumNum = 2*ones(1,options.numPhotons);
photons.nextMedium = round(rand(1,options.numPhotons) + 1);
photons.nextMedium(photons.nextMedium == 2) = 3;
photons.alive = ones(1,options.numPhotons);
normalVectors = interfaceNormalVectors(photons,media,classArgument);
fail_pass = 'FAIL';
if min(abs(vecnorm(normalVectors) - 1) < eps(5)) == 0
    fail_pass = [fail_pass '_VECTORS_NOT_UNIT'];
end
normalCriterion = (abs(photons.position([1,3],:) - vecnorm(photons.position([1,3],:)).*normalVectors([1,3],:) .* sign(normalVectors([3,1],:))) < eps(1));
if min(normalCriterion(:)) == 0
    fail_pass = [fail_pass '_VECTORS_NOT_NORMAL'];
end
inwardCriterion = sign(normalVectors(1,photons.nextMedium < photons.mediumNum)) == -1;
outwardCriterion = sign(normalVectors(1,photons.nextMedium > photons.mediumNum)) == 1;
if min(outwardCriterion) == 0 || min(inwardCriterion) == 0
    fail_pass = [fail_pass '_INWARD_OUTWARD'];
end
if length(fail_pass) > 4
    return
end
fail_pass = 'PASS';
end
%% Unit test for photonEscapeCheck
function fail_pass = testPhotonEscapeCheck(classArgument,options)
arguments
    classArgument       string  {mustBeMember(classArgument,{'single','double','gpuArray'})}
    options.numPhotons  (1,1) {mustBeNumeric, mustBePositive} = 1e6
    options.simRadius   (1,1) {mustBeNumeric, mustBePositive} = rand
end
simShape = {'sphere','cylinder'};
sim.radius = options.simRadius;
sim.halfLength = options.simRadius;
sim.numPhotons = options.numPhotons;
testPathLengths = [sqrt(sim.radius^2 + sim.halfLength^2) sim.radius 0.8*sim.radius];
photons.position = ones(3,1)*zeros(1,options.numPhotons);
source = makePointSource([0;0;0],666e-9);
photons.direction = photonDirection(source,sim,photons.position,classArgument);
photons.direction(2,1:10) = 0;
photons.direction(:,1:10) = photons.direction(:,1:10)./vecnorm(photons.direction(:,1:10));
photons.direction(:,11) = [1;0;0];
photons.direction(:,12) = [-1;0;0];
photons.direction(:,13) = [0;1;0];
photons.direction(:,14) = [0;-1;0];
photons.direction(:,15) = [0;0;1];
photons.direction(:,16) = [0;0;-1];
fail_pass = 'FAIL';
for n = 1:2
    sim.shape = simShape{n};
    for m = 1:3
        photons.alive = ones(1,options.numPhotons);
        photons.pathLengths = testPathLengths(m).*ones(1,options.numPhotons)+eps(1);
        photons = photonEscapeCheck(photons,sim);
        if n == 1 && m == 1 && max(photons.alive) == 1
            fail_pass = [fail_pass '_SPHERE_ALL_OUTSIDE'];
        end
        if n == 1 && m == 2 && max(photons.alive) == 1
            fail_pass = [fail_pass '_SPHERE_ALL_ON_SURFACE'];
        end
        if n == 1 && m == 3 && min(photons.alive) == 0
            fail_pass = [fail_pass '_SPHERE_ALL_INSIDE'];
        end
        if n == 2 && m == 1 && max(photons.alive) == 1
            fail_pass = [fail_pass '_CYLINDER_ALL_OUTSIDE'];
        end
        if n == 2 && m == 2 && (max(photons.alive(1:16)) == 1 || min(photons.alive(17:end)) == 0)
            fail_pass = [fail_pass '_CYLINDER_PATHLENGTH_RADIUS'];
        end
        if n == 2 && m == 3 && min(photons.alive) == 0
            fail_pass = [fail_pass '_CYLINDER_ALL_INSIDE'];
        end
    end
end
if length(fail_pass) > 4
    return
end
fail_pass = 'PASS';
end
%% Unit test for photonEscape
function fail_pass = testPhotonEscapeTruncate(classArgument,options)
arguments
    classArgument       string  {mustBeMember(classArgument,{'single','double','gpuArray'})}
    options.numPhotons  (1,1) {mustBeNumeric, mustBePositive} = 1e5
    options.simRadius   (1,1) {mustBeNumeric, mustBePositive} = rand
    options.tolerance   (1,1) {mustBeNumeric, mustBePositive} = 1e3*eps(1)
end
simShape = {'sphere','cylinder'};
sim.numPhotons = options.numPhotons;
sim.radius = options.simRadius;
sim.halfLength = options.simRadius/2;
media.absorptionCoeff(1) = 0.1;
photons.position = [0;0;0]*ones(1,options.numPhotons);
photons.mediumNum = ones(1,options.numPhotons);
source = makePointSource([0;0;0],666e-9);
photons.direction = photonDirection(source,sim,[],classArgument);
fail_pass = 'FAIL';
for n = 1:2
    sim.shape = simShape{n};
    photons.amplitude = ones(1,options.numPhotons);
    photons.pathLengths = zeros(1,options.numPhotons);
    photons.alive = [zeros(1,round(options.numPhotons/2)) ones(1,options.numPhotons - round(options.numPhotons/2))];
    photons = photonEscapeTruncate(photons,sim,media,classArgument);
    if n == 1 && min(abs(photons.pathLengths(photons.alive == 0) - sim.radius) < options.tolerance) == 0
        fail_pass = [fail_pass '_SPHERE'];
    end
    criterion = abs(vecnorm(photons.pathLengths.*photons.direction([1,3],:)) - options.simRadius) < options.tolerance | abs(abs(photons.pathLengths.*photons.direction(2,:)) - sim.halfLength) < options.tolerance;
    if n == 2 && min(criterion(photons.alive == 0)) == 0
        fail_pass = [fail_pass '_CYLINDER'];
    end
    absorptionCheck = abs(photons.amplitude(photons.alive == 0) - exp(horzcat(media(photons.mediumNum(photons.alive == 0)).absorptionCoeff).*-photons.pathLengths(photons.alive == 0))) < options.tolerance;
    if min(absorptionCheck) == 0
        fail_pass = [fail_pass '_LAMBERT_BEER'];
    end
end
if length(fail_pass) > 4
    return
end
fail_pass = 'PASS';
end
%% Unit test for updatePositions
function fail_pass = unitTestUpdatePositions(media,classArgument,options)
arguments
    media                   struct
    classArgument           string  {mustBeMember(classArgument,{'single','double','gpuArray'})}
    options.numPhotons      (1,1) {mustBeNumeric, mustBePositive} = 1e6
    options.tolerance       (1,1) {mustBeNumeric, mustBePositive} = 1e3*eps(1)
end
simShape = {'sphere','cylinder'};
sim.radius = media.relevantDimension(1,1);
sim.halfLength = media.relevantDimension(1,2);
sim.numPhotons = options.numPhotons;
sim.GPU = 0;

source = makePointSource([0;0;0],666e-9);
photons.direction = photonDirection(source,sim,[],classArgument);
photons.alive = ones(1,options.numPhotons);
photons.mediumNum = 3*ones(1,options.numPhotons);
photons.pathLengths = zeros(1,options.numPhotons);
photons.polarisationVector = ones(3,options.numPhotons).*[1,0,0]';
media.scatterCoeff(3) = 1e6;
for n = 1:2
    sim.shape = simShape{n};
    photons.position = [0;0;0].*ones(1,options.numPhotons);    
    photons.amplitude = ones(1,options.numPhotons);
    photons.distanceTravelled = zeros(1,options.numPhotons);
    photons.opticalPathLength = zeros(1,options.numPhotons);
    photons = updatePositions(photons,media,sim,classArgument);
    propagationCheck = abs(photons.position - photons.pathLengths.*photons.direction) < eps(1);
    absorptionCheck = abs(photons.amplitude - exp(horzcat(media.absorptionCoeff(photons.mediumNum)).*-photons.pathLengths)) < eps(1);
    fail_pass = 'FAIL';
    if min(propagationCheck(:)) == 0
        fail_pass = [fail_pass '_PROPAGATION_VECTOR'];
    end
    if min(absorptionCheck) == 0
        fail_pass = [fail_pass '_LAMBERT_BEER'];
    end
end
if length(fail_pass) > 4
    return
end
fail_pass = 'PASS';
end