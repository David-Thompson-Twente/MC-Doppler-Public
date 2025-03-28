%%  testAnglePickFunctions.m
%   David Thompson, last update 19-03-2025
%   Tests functions that generate inverse CDFs from input phase functions from MatScat:
%   scatterInverseCDF, scatterInverseCDFpolarised, mediaAnglePickFunction
%
%   Inputs:     struct fail_pass carrying any information about previously tested functions
%               struct media with medium information, specifically the phase function
%
%   Outputs:    FAIL or PASS per function
%%
function fail_pass = testAnglePickFunctions(fail_pass,media)
arguments
    fail_pass   struct
    media       struct
end
testedFunctions = {'scatterInverseCDF','scatterInverseCDFpolarised','mediaAnglePickFunction'};
fail_pass.scatterInverseCDF = testScatterInverseCDF(media);
fail_pass.scatterInverseCDFpolarised = testScatterInverseCDFpolarised(media);
fail_pass.mediaAnglePickFunction = testMediaAnglePickFunction(media);
passCounter = 0;
fail_pass.anglePickFunctions = '';
for n = 1:length(testedFunctions)
    if strcmp(fail_pass.(testedFunctions{n}),'PASS')
        fail_pass = rmfield(fail_pass,testedFunctions{n});
        passCounter = passCounter + 1;
    else
        fail_pass.anglePickFunctions = [fail_pass.anglePickFunctions 9 testedFunctions{n} '_' fail_pass.(testedFunctions{n})];
        fail_pass = rmfield(fail_pass,testedFunctions{n});
    end
end
if passCounter == length(testedFunctions)
    fail_pass.anglePickFunctions = 'PASS_ALL';
end
end
%% Unit test for scatterInverseCDF
function fail_pass = testScatterInverseCDF(media,options)
arguments
    media               struct
    options.numRandom   (1,1) {mustBeInteger, mustBePositive} = 1e6;
end
Smatrix = media.phaseFunction{media.scatterCoeff ~= 0};
phaseFunction = (Smatrix(1,:) + Smatrix(2,:))./2;
anglePickFunction = scatterInverseCDF(phaseFunction);
testPhaseFunction = histcounts(ppval(anglePickFunction.scatterPickFunction,rand(1,options.numRandom)),length(phaseFunction));
correlationCoefficients = corrcoef(phaseFunction,testPhaseFunction);
CC = min(correlationCoefficients);
if min(CC) < 0.95
    fail_pass = 'FAIL';
else
    fail_pass = 'PASS';
end
end
%% Unit test for scatterInverseCDFpolarised
function fail_pass = testScatterInverseCDFpolarised(media,options)
arguments
    media                           struct
    options.numAzimuth              (1,1) {mustBeInteger, mustBePositive}       = 1e3
    options.numRandom               (1,1) {mustBeInteger, mustBePositive}       = 1e7
    options.smoothKernelSize        (1,1) {mustBeInteger, mustBePositive}       = 10
    options.showPlots               (1,1) {mustBeNumericOrLogical}              = 0
    options.polarisationAngle       (1,1) {mustBeNumeric, mustBeNonnegative}    = 0
end
phi = linspace(0,2*pi,options.numAzimuth);
Smatrix = media.phaseFunction{media.scatterCoeff ~= 0};
phaseFunction = cos(phi).^2.*Smatrix(1,:)' + sin(phi).^2.*Smatrix(2,:)';
anglePickFunction = scatterInverseCDFpolarised(phaseFunction);
scatteringAngles(1,:) = interp1(anglePickFunction.azimuthPickFunction,phi,rand(1,options.numRandom));
scatteringAngles(2,:) = interp2(linspace(0,1,size(anglePickFunction.scatterPickMatrix,2)),linspace(0,2*pi,size(anglePickFunction.scatterPickMatrix,1)),anglePickFunction.scatterPickMatrix,rand(1,options.numRandom),scatteringAngles(1,:),'spline');
testPhaseFunction = conv2(hist3(scatteringAngles',[options.numAzimuth length(Smatrix)]),(1./options.smoothKernelSize^2)*ones(options.smoothKernelSize),'same');
correlationCoefficients = corrcoef(sum(phaseFunction./max(phaseFunction(:)),2),sum(testPhaseFunction'./max(testPhaseFunction(:)),2));
if options.showPlots
    figure;
    plot(log10(sum(phaseFunction./max(phaseFunction(:)),2)))
    hold on
    plot(log10(sum(testPhaseFunction'./max(testPhaseFunction(:)),2)))
    title(['Correlation ' num2str(min(correlationCoefficients(:)))])
end
CC = min(correlationCoefficients(:));
if min(CC) < 0.95
    fail_pass = 'FAIL';
else
    fail_pass = 'PASS';
end
end
%% Unit test for mediaAnglePickFunction
function fail_pass = testMediaAnglePickFunction(media,options)
arguments
    media               struct
    options.numPoints   (1,1)   {mustBeNumeric,mustBePositive}  = 1e3
end
[~,mediaUnpolarised] = mediaAnglePickFunction(media);

S1 = ppval(mediaUnpolarised.S1(media.scatterCoeff ~= 0),linspace(0,pi,options.numPoints));
S2 = ppval(mediaUnpolarised.S2(media.scatterCoeff ~= 0),linspace(0,pi,options.numPoints));
corrCoeff_S1 = corrcoef(S1,media.phaseFunction{media.scatterCoeff ~= 0}(1,:));
corrCoeff_S2 = corrcoef(S2,media.phaseFunction{media.scatterCoeff ~= 0}(2,:));
fail_pass = 'FAIL';
if min(corrCoeff_S1(:)) >= 0.98 || min(corrCoeff_S2(:)) >= 0.98
    fail_pass = 'PASS';
end
end