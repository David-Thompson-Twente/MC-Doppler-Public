%%  testUtilityFunctions.m
%   David Thompson, 18-03-2024
%   Tests basic functions underpinning McDoppler: rotateRodrigues, cylinderSpatial, emptyDuplicateStruct
%
%   Inputs:     struct fail_pass carrying any information about previously tested functions
%               media struct to test cylinderSpatial and emptyDuplicateStruct
%               classArgument string to determine class to be used, usually double, can also be string or gpuArray
%
%   Outputs:    FAIL or PASS per function
%%
function fail_pass = testUtilityFunctions(fail_pass,media,classArgument)
arguments
    fail_pass       struct
    media           struct
    classArgument   string  {mustBeMember(classArgument,{'single','double','gpuArray'})}
end
testedFunctions = {'rotateRodrigues','cylinderSpatial','emptyDuplicateStruct'};
fail_pass.rotateRodrigues = testRotateRodrigues(classArgument);
fail_pass.cylinderSpatial = testCylinderSpatial(media,classArgument);
fail_pass.emptyDuplicateStruct = testEmptyDuplicateStruct(media);
passCounter = 0;
fail_pass.utilityFunctions = '';
for n = 1:length(testedFunctions)
    if strcmp(fail_pass.(testedFunctions{n}),'PASS')
        fail_pass = rmfield(fail_pass,testedFunctions{n});
        passCounter = passCounter + 1;
    else
        fail_pass.utilityFunctions = [fail_pass.utilityFunctions 9 testedFunctions{n} '_' fail_pass.(testedFunctions{n})];
        fail_pass = rmfield(fail_pass,testedFunctions{n});
    end
end
if passCounter == length(testedFunctions)
    fail_pass.utilityFunctions = 'PASS_ALL';
end
end
%% Unit test for rotateRodrigues
function fail_pass = testRotateRodrigues(classArgument,options)
arguments
    classArgument           string  {mustBeMember(classArgument,{'single','double','gpuArray'})}
    options.numVectors      (1,1)   {mustBeNumeric, mustBePositive} = 1e6
    options.tolerance       (1,1)   {mustBeNumeric, mustBePositive} = 1e3.*eps(1)
end
axisVector = rand(3,options.numVectors,classArgument);
elevationAngle = linspace(eps(1),2*pi,options.numVectors);
azimuthAngle = linspace(0,2*pi,options.numVectors);
axisVector = axisVector./vecnorm(axisVector);
vec_parallel = cos(elevationAngle).*axisVector;
assistVector = axisVector+rand(3,1);
vec_perpendicular = cross(axisVector,assistVector);
vec_perpendicular = vec_perpendicular./vecnorm(vec_perpendicular);
parallelnessCheck = dot(axisVector,assistVector);
if ~isempty(parallelnessCheck(parallelnessCheck == 1))
    warning([num2str(length(parallelnessCheck(parallelnessCheck == 1))) ' unit vector pairs are parallel'])
end
inputVectors = vec_parallel + sin(elevationAngle).*vec_perpendicular;
outputVectors = rotateRodrigues(axisVector,elevationAngle,azimuthAngle,vec_perpendicular,classArgument);
normalVecIn = inputVectors - vec_parallel;
normalVecOut = outputVectors - vec_parallel;

elevationCheck = dot(inputVectors,axisVector) - dot(outputVectors,axisVector) < options.tolerance;%eps(10) a bit arbitrary
azimuthCheck = dot(normalVecOut./vecnorm(normalVecOut),normalVecIn./vecnorm(normalVecIn)) - cos(azimuthAngle) < options.tolerance;

fail_pass = 'FAIL';
if min(elevationCheck) == 0
    fail_pass = [fail_pass '_ELEVATION'];
end
if min(azimuthCheck) == 0
    fail_pass = [fail_pass '_AZIMUTH'];
end
if length(fail_pass) > 4
    return
end
fail_pass = 'PASS';
end
%% Unit test for cylinderSpatial
function fail_pass = testCylinderSpatial(media,classArgument,options)
arguments
    media               struct
    classArgument       string  {mustBeMember(classArgument,{'single','double','gpuArray'})}
    options.numPhotons  (1,1)   {mustBeNumeric, mustBePositive} = 10
end
assistVector = media.axis + rand(3,1);
vec_perpendicular = cross(media.axis,assistVector);
vec_perpendicular = vec_perpendicular./vecnorm(vec_perpendicular);
minRadii = [media.relevantDimension(1,2), media.relevantDimension(1,3), 0];
photonPositions = cell(1,length(media.name));
for n = 1:length(media.name)
photonPositions{n} = (minRadii(n) + 0.9*(media.relevantDimension(1,n)-minRadii(n)).*rand(1,options.numPhotons)).*rotateRodrigues(media.axis(:,n).*ones(1,options.numPhotons),pi/2*ones(1,options.numPhotons),2*pi*rand(1,options.numPhotons),vec_perpendicular(:,n),classArgument)...
    + (2.*media.relevantDimension(2,n).*rand(1,options.numPhotons) - media.relevantDimension(2,n)).*media.axis(:,n);
end
media = cylinderSpatial(media,classArgument);
positionCheck = zeros(length(media.name),options.numPhotons);
testMatrix = ones(3,options.numPhotons);
fail_pass = 'FAIL';
for n = 1:length(media.name)
    for m = 1:length(media.name)
        posVec = photonPositions{m};
        positionCheck(m,:) = media.spatial{n}(posVec(1,:),posVec(2,:),posVec(3,:));
    end
    testCrit = min(positionCheck == testMatrix);
    if min(testCrit(:)) == 0
        fail_pass = [fail_pass '_MEDIUM ' num2str(n)];
    end
    testMatrix(n,:) = 0;
end
if length(fail_pass) > 4
    return
end
fail_pass = 'PASS';
end
%% Unit test for emptyDuplicateStruct - check that an empty struct with the same fields as the input struct is generated
function fail_pass = testEmptyDuplicateStruct(media,options)
arguments
    media               struct
    options.omitFields  cell    = {'name'}
end
emptyCopy = emptyDuplicateStruct(media);
emptyCopy_omitName = emptyDuplicateStruct(media,'omitFields',options.omitFields);
fail_pass = 'FAIL';
copyTest = rmfield(emptyCopy,intersect(fieldnames(emptyCopy),fieldnames(media)));
omitTest = rmfield(media,intersect(fieldnames(media),fieldnames(emptyCopy_omitName)));
if ~isempty(fieldnames(copyTest))
    fail_pass = [fail_pass '_COPY'];
end
if ~sum(strcmp(fieldnames(omitTest),fieldnames(media))) == 1 || ~strcmp(fieldnames(omitTest),options.omitFields)
    fail_pass = [fail_pass '_OMIT'];
end
if length(fail_pass) > 4
    return
end
fail_pass = 'PASS';
end