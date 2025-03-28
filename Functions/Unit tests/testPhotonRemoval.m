%%  testPhotonRemoval.m
%   David Thompson, 27-03-2024
%   Tests the functions photonKill and photonDetect, which remove photons from a McDoppler simulation
%
%   Inputs:     struct fail_pass carrying any information about previously tested functions
%               classArgument string to determine class to be used, usually double, can also be string or gpuArray
%
%   Outputs:    FAIL or PASS per sub-function
%%
function fail_pass = testPhotonRemoval(fail_pass,classArgument,options)
arguments
    fail_pass           struct
    classArgument       string  {mustBeMember(classArgument,{'single','double','gpuArray'})}
    options.numPhotons  (1,1)   {mustBeNumeric,mustBePositive}  = 1e6;
    options.omitFields  cell                                    = {'killThreshold','alive'}
end
fail_pass.photonKill = 'FAIL';
fail_pass.photonDetect = 'FAIL';

photons.alive = ones(1,options.numPhotons);
photons.killThreshold = 0.1;
photons.amplitude = [ones(1,floor(options.numPhotons/2)) (photons.killThreshold/2).*ones(1,options.numPhotons - floor(options.numPhotons/2))];
photons.testNumbers = 1:options.numPhotons;

killedPhotons = photonKill(photons,classArgument);

killTestAmplitude = sum(killedPhotons.alive(1:floor(options.numPhotons/2))) == floor(options.numPhotons/2) && sum(killedPhotons.alive(options.numPhotons - floor(options.numPhotons/2):end)) < options.numPhotons - floor(options.numPhotons/2);
killTestRoulette = sum(killedPhotons.alive(options.numPhotons - floor(options.numPhotons/2):end))./(options.numPhotons - floor(options.numPhotons/2)) - photons.killThreshold < .01.*photons.killThreshold;
if ~killTestAmplitude
    fail_pass.photonKill = [fail_pass.photonKill '_PHOTONS_ABOVE_THRESHOLD_KILLED'];
end
if ~killTestRoulette
    fail_pass.photonKill = [fail_pass.photonKill '_ROULETTE_KILLTHRESHOLD'];
end
if strcmp(fail_pass.photonKill,'FAIL')
    fail_pass.photonKill = 'PASS';
end

detectedPhotons = emptyDuplicateStruct(photons,'omitFields',options.omitFields);
[photons,detectedPhotons] = photonDetect(killedPhotons,detectedPhotons);

correctRemovalTest = detectedPhotons.testNumbers == killedPhotons.testNumbers(killedPhotons.alive == 0);
correctRemainingTest = photons.testNumbers == killedPhotons.testNumbers(killedPhotons.alive == 1);
detectedNames = fieldnames(detectedPhotons);
inputNames = fieldnames(photons);
inputNames(ismember(inputNames,options.omitFields)) = [];
fieldsTest = horzcat(detectedNames{:}) == horzcat(inputNames{:});

if min(correctRemovalTest) == 0
    fail_pass.photonDetect = [fail_pass.photonDetect '_WRONG_PHOTON_REMOVED'];
end
if min(correctRemainingTest) == 0
    fail_pass.photonDetect = [fail_pass.photonDetect '_WRONG_PHOTON_REMAINING'];
end
if min(fieldsTest) == 0
    fail_pass.photonDetect = [fail_pass.photonDetect '_FIELDS_OMITTED'];
end
if strcmp(fail_pass.photonDetect,'FAIL')
    fail_pass.photonDetect = 'PASS';
end

testedFunctions = {'photonKill','photonDetect'};
passCounter = 0;
fail_pass.photonRemoval = '';
for n = 1:length(testedFunctions)
    if strcmp(fail_pass.(testedFunctions{n}),'PASS')
        fail_pass = rmfield(fail_pass,testedFunctions{n});
        passCounter = passCounter + 1;
    else
        fail_pass.photonRemoval = [fail_pass.photonRemoval 9 testedFunctions{n} '_' fail_pass.(testedFunctions{n})];
        fail_pass = rmfield(fail_pass,testedFunctions{n});
    end
end
if passCounter == length(testedFunctions)
    fail_pass.photonRemoval = 'PASS_ALL';
end
end