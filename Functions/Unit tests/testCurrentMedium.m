%%  testCurrentMedium.m
%   David Thompson, 20-03-2024
%   Tests if photons placed in a particular medium are indeed classified as such
%   by the function currentMedium
%
%   Inputs:     struct fail_pass carrying any information about previously tested functions
%               media struct to test cylinderSpatial and emptyDuplicateStruct
%
%   Outputs:    FAIL of PASS per function
%%
function fail_pass = testCurrentMedium(fail_pass,media,options)
arguments
    fail_pass                   struct
    media                       struct
    options.numPhotons          (1,1) {mustBeNumeric} = 1e5;
end
mediumRadii = media.relevantDimension(1,:);
mediumHalfLengths = media.relevantDimension(2,:);
radiusVec = [0;0;1]*ones(1,options.numPhotons);
photons.position = mediumRadii(1).*rand(1,options.numPhotons).*rotateRodrigues(radiusVec,2*pi*rand(1,options.numPhotons),0,[1;0;0],'double') + (2*mediumHalfLengths(3).*rand(1,options.numPhotons) - mediumHalfLengths(3)).*[0;1;0];
photons.mediumNum = zeros(1,options.numPhotons);
photons = currentMedium(photons,media);

medium1 = sqrt(photons.position(1,:).^2 + photons.position(3,:).^2) <= mediumRadii(1) & abs(photons.position(2,:)) <= mediumHalfLengths(1);
medium2 = sqrt(photons.position(1,:).^2 + photons.position(3,:).^2) <= mediumRadii(2) & abs(photons.position(2,:)) <= mediumHalfLengths(2);
medium3 = sqrt(photons.position(1,:).^2 + photons.position(3,:).^2) <= mediumRadii(3) & abs(photons.position(2,:)) <= mediumHalfLengths(3);
medium1crit = medium1 & ~medium2 & ~medium3;
medium2crit = medium2 & ~medium3;
medium3crit = medium3;

if ~min(photons.mediumNum(medium1crit) == 1) == 0 && ~min(photons.mediumNum(medium2crit) == 2) == 0 && ~min(photons.mediumNum(medium3crit) == 3) == 0
    fail_pass.currentMedium = 'PASS';
else
    fail_pass.currentMedium = 'FAIL';
end
end