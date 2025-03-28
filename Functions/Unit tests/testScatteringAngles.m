%%  testScatteringAngles.m
%   David Thompson 22-03-2024, last update 19-03-2025
%   Tests functions that pick scattering angles for any photons that are to be scattered:
%   scatteringAngleQuantile, scatterinAngleQuantileBiased, scatterinAnglePolarised, scatterinAnglePolarisedBiased
%
%   Inputs:     struct fail_pass carrying any information about previously tested functions
%               media struct with all required media properties
%               classArgument string to determine class to be used, usually double, can also be string or gpuArray
%
%   Outputs:    PASS_ALL if all subtests succesfull
%               FAIL message indicating which part failed
%%
function fail_pass = testScatteringAngles(fail_pass,media,classArgument)
arguments
    fail_pass       struct
    media           struct
    classArgument   string  {mustBeMember(classArgument,{'single','double','gpuArray'})}
end
testedFunctions = {'scatteringAngleQuantile','scatteringAnglePolarised'};
fail_pass.scatteringAngleQuantile = testScatteringAngleQuantile(media,classArgument);
fail_pass.scatteringAnglePolarised = testScatteringAnglePolarised(media,classArgument);
passCounter = 0;
fail_pass.scatteringAngles = '';
for n = 1:length(testedFunctions)
    if strcmp(fail_pass.(testedFunctions{n}),'PASS')
        fail_pass = rmfield(fail_pass,testedFunctions{n});
        passCounter = passCounter + 1;
    else
        fail_pass.scatteringAngles= [fail_pass.scatteringAngles 9 testedFunctions{n} '_' fail_pass.(testedFunctions{n})];
        fail_pass = rmfield(fail_pass,testedFunctions{n});
    end
end
if passCounter == length(testedFunctions)
    fail_pass.scatteringAngles = 'PASS_ALL';
end
end
%% Unit test for scatteringAngleQuantile
function fail_pass = testScatteringAngleQuantile(media,classArgument,options)
arguments
    media                   struct
    classArgument           string  {mustBeMember(classArgument,{'single','double','gpuArray'})}
    options.numPhotons      (1,1)   {mustBeNumeric, mustBePositive}     = 1e6
    options.wavelength      (1,1)   {mustBeNumeric, mustBePositive}     = 666e-9
    options.showPlots       (1,1)   {mustBeNumericOrLogical}            = 0
    options.fitTolerance    (1,1)   {mustBeNumeric, mustBePositive}     = 1e-4;
end
theta = linspace(0,pi,length(media.phaseFunction{media.scatterCoeff ~= 0}));
photons.alive = ones(1,options.numPhotons);
photons.mediumNum = find(media.scatterCoeff ~= 0).*ones(1,options.numPhotons);
scatteringPhotons = ones(1,options.numPhotons);
Smatrix = media.phaseFunction{media.scatterCoeff ~= 0};
phaseFunction = (Smatrix(1,:) + Smatrix(2,:))./2;
[anglePickFunction,media] = mediaAnglePickFunction(media);
scatteringAngles = scatteringAngleQuantile(anglePickFunction,photons,scatteringPhotons,media,classArgument);
testPhaseFunction = histcounts(scatteringAngles(2,:),length(phaseFunction),'BinLimits',[0,pi]);
elevationCorrelation = corrcoef(phaseFunction.*sin(theta),testPhaseFunction);
azimuthTest = histcounts(scatteringAngles(1,:),length(phaseFunction),'BinLimits',[0,2*pi]);
azimuthFit = fit((1:length(phaseFunction))',azimuthTest','poly1');
if options.showPlots
    figure;
    subplot(211)
    plot(log10((sin(theta).*phaseFunction)/max(sin(theta).*phaseFunction)))
    hold on
    plot(log10(testPhaseFunction/max(testPhaseFunction)))
    subplot(212)
    plot(azimuthTest)
    hold on
    plot(azimuthFit(1:length(phaseFunction)))
end
fail_pass = 'FAIL';
if min(elevationCorrelation) < 0.98
    fail_pass = [fail_pass '_ELEVATION'];
end
if abs(azimuthFit.p1./azimuthFit.p2) > options.fitTolerance
    fail_pass = [fail_pass '_AZIMUTH'];
end
if length(fail_pass) > 4
    return
end
fail_pass = 'PASS';
end
%% Unit test for scatteringAnglePolarised
function fail_pass = testScatteringAnglePolarised(media,classArgument,options)
arguments
    media                       struct
    classArgument               string  {mustBeMember(classArgument,{'single','double','gpuArray'})}
    options.numPhotons          (1,1)   {mustBeNumeric, mustBePositive}     = 1e6
    options.wavelength          (1,1)   {mustBeNumeric, mustBePositive}     = 666e-9
    options.showPlots           (1,1)   {mustBeNumericOrLogical}            = 0
    options.smoothKernelSize    (1,1)   {mustBeInteger, mustBePositive}     = 10;
end
numBins = length(media.phaseFunction{media.scatterCoeff ~= 0});
theta = linspace(0,pi,numBins);
photons.alive = ones(1,options.numPhotons);
photons.mediumNum = find(media.scatterCoeff ~= 0).*ones(1,options.numPhotons);
scatteringPhotons = ones(1,options.numPhotons);
phi = linspace(0,2*pi,numBins);
polarisationAngle = 2*pi*rand(1,options.numPhotons);
Smatrix = media.phaseFunction{media.scatterCoeff ~= 0};
phaseFunction = cos(phi).^2.*Smatrix(1,:)' + sin(phi).^2.*Smatrix(2,:)';
[anglePickFunction,media] = mediaAnglePickFunction(media,'polarisation','polarised','polarisationAngle',polarisationAngle);
scatteringAngles = scatteringAnglePolarised(anglePickFunction,photons,scatteringPhotons,media,classArgument);
phaseFunction = phaseFunction.*sin(theta)';
elevationPF = sum(phaseFunction,2);
testElevation = histcounts(scatteringAngles(2,:),length(theta),'BinLimits',[0,pi]);
azimuthPF = sum(phaseFunction,1);
testAzimuth = histcounts(scatteringAngles(1,:),length(phi),'BinLimits',[0,2*pi]);
azimuthFitFun = fittype(@(a,b,x) a.*cos(x).^2 + b.*sin(x).^2);
testAzimuthFit = fit(phi',(max(azimuthPF).*testAzimuth./max(testAzimuth))',azimuthFitFun,'StartPoint',[sum(Smatrix(1,:).*sin(theta)),sum(Smatrix(2,:).*sin(theta))]);
phaseFunctionCorrelationTheta = corrcoef(elevationPF./max(elevationPF),testElevation./max(testElevation));
normAmplitudeTestAzimuth = max(testAzimuthFit(phi)./max(testAzimuthFit(phi))) - min(testAzimuthFit(phi)./max(testAzimuthFit(phi)));
normAmplitudeAzimuthPF = max(azimuthPF./max(azimuthPF)) - min(azimuthPF./max(azimuthPF));
deviationPercentage = (abs(normAmplitudeTestAzimuth-normAmplitudeAzimuthPF))*1e2;
if options.showPlots
    figure;
    subplot(121)
    plot(elevationPF./max(elevationPF))
    hold on
    plot(testElevation./max(testElevation))
    title(['Correlation ' num2str(min(phaseFunctionCorrelationTheta(:)))])
    subplot(122)
    plot(azimuthPF./max(azimuthPF))
    hold on
    plot(testAzimuth./max(testAzimuth),'LineWidth',1.5)
    plot(testAzimuthFit(phi)./max(testAzimuthFit(phi)),'LineWidth',1.5,'LineStyle','--')
    title(['Deviation percentage ' num2str(deviationPercentage) '%'])
end
fail_pass = 'FAIL';
if min(phaseFunctionCorrelationTheta(:)) < 0.95
    fail_pass = [fail_pass '_ELEVATION'];
end
if deviationPercentage > 0.5
    fail_pass = [fail_pass '_AZIMUTH'];
end
if ~strcmp(fail_pass,'FAIL')
    return
end
fail_pass = 'PASS';
end