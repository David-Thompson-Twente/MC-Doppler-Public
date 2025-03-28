%%  McDopplerUnitTests.m
%   David Thompson, last update 19-03-2025
%   Unit test suite for McDoppler
%   
%   Inputs:     (optional)  classArgument, in which class are all variables to be
%               (optional)  randomRadii, radii of the three media (cylindrical), standard values random [0,1]
%               (optional)  mediaLength, length of all media and simulation itself
%               (optional)  scatterLottery, random number that picks which is the scattering medium (only one)
%               (optional)  volumeFlowRate
%               (optional)  particleSize, size of scattering particles used to generate phase function
%               (optional)  particleIndex, particle refractive index
%               (optional)  wavelength of the light
%
%   Outputs:    PASS/FAIL per function / subfunction
%%
function fail_pass = McDopplerUnitTests(options)
arguments
    options.classArgument   string  {mustBeMember(options.classArgument,{'single','double','gpuArray'})}    = 'double'
    options.randomRadii     (3,1)   {mustBeNumeric,mustBePositive}                                          = rand(3,1)
    options.mediaLength     (3,1)   {mustBeNumeric, mustBePositive}                                         = rand.*ones(3,1)
    options.scatterLottery  (1,1)   {mustBeNumericOrLogical}                                                = floor(3*rand + 1)
    options.volumeFlowRate  (1,1)   {mustBeNumeric,mustBePositive}                                          = 1e-6*rand
    options.particleSize    (1,1)   {mustBeNumeric,mustBePositive}                                          = 5e-6*rand
    options.particleIndex   (1,1)   {mustBeNumeric,mustBePositive}                                          = rand + 1
    options.wavelength      (1,1)   {mustBeNumeric,mustBePositive}                                          = 666e-9
end
%Define media to be used for all tests
media.name                  = [{'medium1'} {'medium2'} {'medium3'}];
media.shape                 = [{'cylinder'} {'cylinder'} {'cylinder'}];
media.relevantDimension     = [flip(sort(options.randomRadii)) options.mediaLength]';
media.axis                  = [0 0 0;1 1 1;0 0 0];%ambiguity introduced if media cross one another's boundaries within sim
media.axisorigin            = zeros(3,3);
media.scatterCoeff          = zeros(1,3);
media.absorptionCoeff       = rand(1,3);
media.refractiveIndex       = rand(1,3) + 1;
media.shiftX                = zeros(1,3);
media.flowVelocity          = [{@(x,y,z) [zeros(size(x));zeros(size(y));zeros(size(z))]} {@(x,y,z) [zeros(size(x));zeros(size(y));zeros(size(z))]} {@(x,y,z) [zeros(size(x));zeros(size(y));zeros(size(z))]}];
media.phaseFunction         = cell(1,3);
media.surroundingMedium     = [0, 1, 2];

tubeInnerRadius = media.relevantDimension(1,options.scatterLottery);
media.scatterCoeff(options.scatterLottery) = 1 + rand;%quite low, not a big deal for most tests though, could consider changing it
media.flowVelocity(options.scatterLottery) = {@(x,y,z) [zeros(size(x));(2*options.volumeFlowRate)/(pi*tubeInnerRadius^2).*(1 - (sqrt(x.^2+z.^2)/tubeInnerRadius).^2);zeros(size(z))]};
media = cylinderSpatial(media,options.classArgument);
media.phaseFunction{options.scatterLottery} = simplePhaseFunction(options.particleSize,options.particleIndex,media.refractiveIndex(options.scatterLottery),options.wavelength,1e3);
fail_pass = struct();

%% Tests for simulation prep functions
%Test utility functions underlying simulation set up
disp([9 'testing utility functions'])
fail_pass = testUtilityFunctions(fail_pass,media,options.classArgument);
%Test generation of light sources
disp([9 'testing light sources'])
fail_pass = testLightSourceFunctions(fail_pass);
%Test generation of scattering inverse CDFs
disp([9 'testing scattering inverse CDFs'])
fail_pass = testAnglePickFunctions(fail_pass,media);
%% Tests for simulation startup
%Test photonLaunch for all source types
disp([9 'testing photonLaunch'])
fail_pass = testPhotonLaunch(fail_pass,options.classArgument);
%Test for starting medium
disp([9 'testing starting medium'])
fail_pass = testCurrentMedium(fail_pass,media);
%% Tests for updatePositions and related functions
disp([9 'testing updatePositions'])
fail_pass = testUpdatePositions(fail_pass,media,options.classArgument);
%% Tests for updateDirections and related functions
disp([9 'testing reflection and refraction'])
fail_pass = testReflectRefract(fail_pass,options.classArgument);
%Test scattering angle selection
disp([9 'testing scattering angle calculations'])
fail_pass = testScatteringAngles(fail_pass,media,options.classArgument);
%Test scattering
disp([9 'testing scattering behaviour'])
fail_pass = testUpdateDirectionsScatter(fail_pass,media,options.classArgument);
%Test general behaviour of updateDirections
disp([9 'testing updateDirections'])
fail_pass = testUpdateDirections(fail_pass,media,options.classArgument);
%% Simulation wrap-up tests
%Test removal of photons that are either absorbed or that ecape the simulation, subsequent storage in separate struct 
disp([9 'testing photon removal and gathering of photons from GPU'])
fail_pass = testPhotonRemoval(fail_pass,options.classArgument);
%Test gathering of photons struct from GPU and concatenation of structs in multi-run simulations
fail_pass = testGatherConcatenate(fail_pass);
%Test selection of subset of photons based on finite-sized detectors
disp([9 'testing detector functions'])
fail_pass = testDetectors(fail_pass,media,options.classArgument);
end