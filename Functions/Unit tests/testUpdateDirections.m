%%  testUpdateDirections.m
%   David Thompson 27-03-2024
%   Tests if the function updateDirections does in fact only change the directions of those photons that are
%   supposed to reflect, refract or scatter, leaving the others unchanged
%
%   Inputs:     struct fail_pass carrying any information about previously tested functions
%               media struct with all required media properties
%               classArgument string to determine class to be used, usually double, can also be string or gpuArray
%
%   Outputs:    PASS or FAIL per component
%%
function fail_pass = testUpdateDirections(fail_pass, media, classArgument,options)
arguments
    fail_pass                   struct
    media                       struct
    classArgument               string  {mustBeMember(classArgument,{'single','double','gpuArray'})}
    options.numPhotons          (1,1)   {mustBeNumeric, mustBePositive}     = 1e6
    options.wavelength          (1,1)   {mustBeNumeric, mustBePositive}     = 666e-9
    options.tolerance           (1,1)   {mustBeNumeric, mustBePositive}     = 1e3*eps(1)
end
fail_pass.updateDirections = 'FAIL';
polarisation = {'unpolarised','polarised'};
for n = 1:2
    scatteringMediumNum = find(media.scatterCoeff ~= 0);
    photons.position = zeros(3,options.numPhotons);
    photons.pathLengths = rand(1,options.numPhotons);
    photons.DopplerShift = zeros(1,options.numPhotons);
    photons.mediumNum = find(media.scatterCoeff == 0,1,'first').*ones(1,options.numPhotons);
    photons.nextMedium = photons.mediumNum;
    photons.nextMedium(1:floor(2*end/3)) = scatteringMediumNum;
    photons.mediumNum(ceil(end/3):floor(2*end/3)) = scatteringMediumNum;
    randomDirections = rand(3,options.numPhotons);
    photons.direction = randomDirections./vecnorm(randomDirections);
    photons.polarisationVector = rotateRodrigues(photons.direction,pi/2,0,[],classArgument);
    photons.timesScattered = zeros(1,options.numPhotons);
    photons.timesReflected = zeros(1,options.numPhotons);
    photons.timesRefracted = zeros(1,options.numPhotons);
    photons.alive = ones(1,options.numPhotons);
    photons.likelihood = zeros(1,options.numPhotons);
    photons.wavelength = options.wavelength;
    normalVectors = zeros(3,options.numPhotons);
    normalVectors(:,1:floor(end/3)) = rotateRodrigues(-photons.direction(:,1:floor(end/3)),pi/2.*rand(1,floor(options.numPhotons/3)),2*pi.*rand(1,floor(options.numPhotons/3)),[],classArgument);
    [anglePickFunction,media] = mediaAnglePickFunction(media,'polarisation',polarisation{n});

    processedPhotons = updateDirections(photons,media,normalVectors,anglePickFunction,polarisation{n},classArgument);
    
    scatterNumberCheck = sum(processedPhotons.timesScattered) == length(ceil(options.numPhotons/3):2*floor(options.numPhotons/3));
    reflectRefractNumberCheck = sum(processedPhotons.timesReflected) + sum(processedPhotons.timesRefracted) == floor(options.numPhotons./3);
    nonUpdateCheck = dot(photons.direction(:,ceil(2*end/3):end),processedPhotons.direction(:,ceil(2*end/3):end)) - 1 < options.tolerance;
    scatterCheck = dot(photons.direction(:,ceil(end/3):2*floor(end/3)),processedPhotons.direction(:,ceil(end/3):2*floor(end/3))) < 1;
    reflectRefractCheck = dot(photons.direction(:,1:floor(end/3)),processedPhotons.direction(:,1:floor(end/3))) < 1;
    if min(scatterNumberCheck) == 0
        fail_pass.updateDirections = [fail_pass.updateDirections '_' polarisation{n} '_TIMES_SCATTERED'];
    end
    if min(reflectRefractNumberCheck) == 0
        fail_pass.updateDirections = [fail_pass.updateDirections '_' polarisation{n} '_TIMES_REFLECTED_REFRACTED'];
    end
    if min(nonUpdateCheck) == 0
        fail_pass.updateDirections = [fail_pass.updateDirections '_' polarisation{n} '_WRONG_PHOTONS_UPDATED'];
    end
    if min(scatterCheck) == 0
        fail_pass.updateDirections = [fail_pass.updateDirections '_' polarisation{n} '_PHOTONS_NOT_SCATTERED'];
    end
    if min(reflectRefractCheck) == 0
        fail_pass.updateDirections = [fail_pass.updateDirections '_' polarisation{n} '_PHOTONS_NOT_REFLECTED_REFRACTED'];
    end
end
if length(fail_pass.updateDirections) > 4
    return
end
fail_pass.updateDirections = 'PASS';
end