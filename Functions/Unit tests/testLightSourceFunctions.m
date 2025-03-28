%%  testLightSourceFunctions.m
%   David Thompson, 19-03-2024
%   Tests functions that generate light source properties for McDoppler:
%   makePointSource, makePencilSource,makeFocussedSource,makeRatioSource,makeTruncatedRatioSource
%
%   Inputs:     struct fail_pass carrying any information about previously tested functions
%
%   Outputs:    FAIL or PASS per function
%%
function fail_pass = testLightSourceFunctions(fail_pass)
arguments
    fail_pass   struct
end
testedFunctions = {'makePointSource','makePencilSource','makeFocussedSource','makeRatioSource','makeTruncatedRatioSource'};
fail_pass.makePointSource = testMakePointSource;
fail_pass.makePencilSource = testMakePencilSource;
fail_pass.makeFocussedSource = testMakeFocussedSource;
fail_pass.makeRatioSource = testMakeRatioSource;
fail_pass.makeTruncatedRatioSource = testMakeTruncatedRatioSource;
passCounter = 0;
fail_pass.lightSourceFunctions = '';
for n = 1:length(testedFunctions)
    if strcmp(fail_pass.(testedFunctions{n}),'PASS')
        fail_pass = rmfield(fail_pass,testedFunctions{n});
        passCounter = passCounter + 1;
    else
        fail_pass.lightSourceFunctions = [fail_pass.anglePickFunctions 9 testedFunctions{n} '_' fail_pass.(testedFunctions{n})];
        fail_pass = rmfield(fail_pass,testedFunctions{n});
    end
end
if passCounter == length(testedFunctions)
    fail_pass.lightSourceFunctions = 'PASS_ALL';
end
end
%% Unit test for makePointSource
function fail_pass = testMakePointSource(options)
arguments
    options.numRuns         (1,1) {mustBeInteger, mustBePositive}   = 100
    options.wavelength      (1,1) {mustBeNumeric, mustBePositive}   = 666e-9
end
for n = 1:options.numRuns
    sourcePosition = rand(3,1);
    source = makePointSource(sourcePosition,options.wavelength);
    if max(source.position ~= sourcePosition) == 1 && max([0;0;1] ~= source.opticalAxis) == 0
        fail_pass = ['FAIL_POSITION_N = ' num2str(n)];
        return
    elseif max([0;0;1] ~= source.opticalAxis) == 1 && max(source.position ~= sourcePosition) == 0
        fail_pass = ['FAIL_OPTICAL_AXIS_N = ' num2str(n)];
        return
    elseif max([0;0;1] ~= source.opticalAxis) == 1 && max(source.position ~= sourcePosition) == 1
        fail_pass = ['FAIL_ALL_N = ' num2str(n)];
        return
    end
end
fail_pass = 'PASS';
end
%% Unit test for makePencilSource
function fail_pass = testMakePencilSource(options)
arguments
    options.numRuns         (1,1) {mustBeInteger, mustBePositive}   = 100
    options.wavelength      (1,1) {mustBeNumeric, mustBePositive}   = 666e-9
end
for n = 1:options.numRuns
    sourcePosition = rand(3,1);
    focusPosition = rand(3,1);
    testAxis = (focusPosition - sourcePosition)./vecnorm(focusPosition - sourcePosition);
    source = makePencilSource(sourcePosition,focusPosition,options.wavelength);
    if max(source.position ~= sourcePosition) == 1 && max(testAxis ~= source.opticalAxis) == 0
        fail_pass = ['FAIL_POSITION_N = ' num2str(n)];
        return
    elseif max(testAxis ~= source.opticalAxis) == 1 && max(source.position ~= sourcePosition) == 0
        fail_pass = ['FAIL_OPTICAL_AXIS_N = ' num2str(n)];
        return
    elseif max(testAxis ~= source.opticalAxis) == 1 && max(source.position ~= sourcePosition) == 1
        fail_pass = ['FAIL_ALL_N = ' num2str(n)];
        return
    end
end
fail_pass = 'PASS';
end
%% Unit test for makeFocussedSource
function fail_pass = testMakeFocussedSource(options)
arguments
    options.numRuns             (1,1) {mustBeInteger, mustBeGreaterThan(options.numRuns,5)} = 100;
    options.refractiveIndex     (1,1) {mustBeNumeric, mustBePositive}                       = 1.4;
    options.pupilDiameter       (1,1) {mustBeNumeric, mustBePositive}                       = 5e-3;
    options.wavelength          (1,1) {mustBeNumeric, mustBePositive}                       = 666e-9;
end
sourcePosMatrix = [[1,1,0,0,0,0;0,0,1,1,0,0;0,0,0,0,1,1] rand(3,options.numRuns-6)];
focusPosMatrix = [zeros(3,6) rand(3,options.numRuns - 6)];
compareAxisMatrix = [[-1,1,0,0,0,0;0,0,-1,1,0,0;0,0,0,0,-1,1] (-1).^(6:options.numRuns-1).*(focusPosMatrix(:,7:options.numRuns) - sourcePosMatrix(:,7:options.numRuns))];
for n = 1:options.numRuns
    focalLength = (-1)^(n-1).*rand;
    sourcePosition = sourcePosMatrix(:,n);
    focusPosition = focusPosMatrix(:,n);
    compareAxis = compareAxisMatrix(:,n)./vecnorm(compareAxisMatrix(:,n));
    source = makeFocussedSource(sourcePosition,options.refractiveIndex,focalLength,focusPosition,options.pupilDiameter,options.wavelength);
    if min(compareAxis == source.opticalAxis) == 0
        fail_pass = ['FAIL_N = ' num2str(n)];
        return
    end
end
fail_pass = 'PASS';
end
%% Unit test for makeRatioSource
function fail_pass = testMakeRatioSource(options)
arguments
    options.numRuns     (1,1) {mustBeInteger, mustBeGreaterThan(options.numRuns,5)} = 100;
    options.wavelength  (1,1) {mustBeNumeric, mustBePositive}                       = 666e-9;
end
sourcePosMatrix = [[1,-1,0,0,0,0;0,0,1,-1,0,0;0,0,0,0,1,-1] rand(3,options.numRuns-6)];
focusPosMatrix = [zeros(3,6) rand(3,options.numRuns - 6)];
compareAxisMatrix = [[-1,1,0,0,0,0;0,0,-1,1,0,0;0,0,0,0,-1,1] (focusPosMatrix(:,7:options.numRuns) - sourcePosMatrix(:,7:options.numRuns))];
for n = 1:options.numRuns
    sourcePosition = sourcePosMatrix(:,n);
    focusPosition = focusPosMatrix(:,n);
    compareAxis = compareAxisMatrix(:,n)./vecnorm(compareAxisMatrix(:,n));
    sourceRadius = rand;
    focusRadius = rand;
    compareRatio = focusRadius/sourceRadius;
    source = makeRatioSource(sourcePosition, focusPosition, sourceRadius, focusRadius, options.wavelength);
    if min(compareAxis == source.opticalAxis) == 0 && compareRatio == source.focusRatio
        fail_pass = ['FAIL_AXIS_N = ' num2str(n)];
        return
    elseif min(compareAxis == source.opticalAxis) == 1 && compareRatio ~= source.focusRatio
        fail_pass = ['FAIL_RATIO_N = ' num2str(n)];
        return
    elseif min(compareAxis == source.opticalAxis) == 0 && compareRatio ~= source.focusRatio
        fail_pass = ['FAIL_ALL_N = ' num2str(n)];
        return
    end
end
fail_pass = 'PASS';
end
%% Unit test for makeTruncatedRatioSource
function fail_pass = testMakeTruncatedRatioSource(options)
arguments
    options.numRuns     (1,1) {mustBeInteger, mustBeGreaterThan(options.numRuns,5)} = 100;
    options.wavelength  (1,1) {mustBeNumeric, mustBePositive}                       = 666e-9;
end
sourcePosMatrix = [[1,-1,0,0,0,0;0,0,1,-1,0,0;0,0,0,0,1,-1] rand(3,options.numRuns-6)];
focusPosMatrix = [zeros(3,6) rand(3,options.numRuns - 6)];
compareAxisMatrix = [[-1,1,0,0,0,0;0,0,-1,1,0,0;0,0,0,0,-1,1] (focusPosMatrix(:,7:options.numRuns) - sourcePosMatrix(:,7:options.numRuns))];
for n = 1:options.numRuns
    sourcePosition = sourcePosMatrix(:,n);
    focusPosition = focusPosMatrix(:,n);
    compareAxis = compareAxisMatrix(:,n)./vecnorm(compareAxisMatrix(:,n));
    sourceRadius = rand;
    focusRadius = rand;
    compareRatio = focusRadius/sourceRadius;
    source = makeTruncatedRatioSource(sourcePosition, focusPosition, sourceRadius, focusRadius, options.wavelength);
    if min(compareAxis == source.opticalAxis) == 0 && compareRatio == source.focusRatio
        fail_pass = ['FAIL_AXIS_N = ' num2str(n)];
        return
    elseif min(compareAxis == source.opticalAxis) == 1 && compareRatio ~= source.focusRatio
        fail_pass = ['FAIL_RATIO_N = ' num2str(n)];
        return
    elseif min(compareAxis == source.opticalAxis) == 0 && compareRatio ~= source.focusRatio
        fail_pass = ['FAIL_ALL_N = ' num2str(n)];
        return
    end
end
fail_pass = 'PASS';
end