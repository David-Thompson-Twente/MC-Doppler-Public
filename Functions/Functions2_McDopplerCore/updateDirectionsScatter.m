%%  updateDirectionsScatter.m
%   David Thompson 28-04-2022, last update 19-03-2025
%   Updates photons.direction for scattered photons
%
%   Inputs:     [struct]        photons
%               [1xn logical]   crossingPhotons
%               [2xn double]    scatteringAngles
%               [struct]        media
%               [string]        classArgument
%
%   Outputs:    [struct]        photons, with updated   photons.direction
%                                                       photons.polarisationVector
%                                                       photons.DopplerShift
%                                                       photons.timesScattered
%                                                       photons.lastPositionScattered

%%
function photons = updateDirectionsScatter(photons,scatteringPhotons,scatteringAngles,media,classArgument)
arguments
    photons             struct
    scatteringPhotons   (1,:) {mustBeNumericOrLogical(scatteringPhotons),mustBeInteger(scatteringPhotons),mustBeLessThanOrEqual(scatteringPhotons,1),mustBeNonnegative(scatteringPhotons)}
    scatteringAngles    (2,:) {mustBeNumeric(scatteringAngles),mustBeNonnegative(scatteringAngles)}
    media               struct
    classArgument       string
end
rotationAngles = scatteringAngles(1,:);
scatteringAngles = scatteringAngles(2,:);
oldDirection = photons.direction;
S1 = zeros(1,length(photons.alive),classArgument);
S2 = zeros(1,length(photons.alive),classArgument);
DopplerShift = zeros(1,length(photons.alive),classArgument);
scatterVector = zeros(3,length(photons.alive),classArgument);
if ~isempty(scatteringPhotons(scatteringPhotons == 1))
    %direction change from scattering
    photons.direction(:,scatteringPhotons) = rotateRodrigues(oldDirection(:,scatteringPhotons),scatteringAngles(scatteringPhotons),rotationAngles(scatteringPhotons),photons.polarisationVector(:,scatteringPhotons),classArgument);

    %update polarisation 1: define unit vectors relative to scattering plane
    inPlaneIncoming = rotateRodrigues(oldDirection,pi/2.*ones(1,length(scatteringAngles),classArgument),rotationAngles,photons.polarisationVector,classArgument);
    outOfPlaneScattered = cross(inPlaneIncoming,oldDirection,1);
    inPlaneScattered = cross(photons.direction,outOfPlaneScattered,1);

    %update polarisation 2: calculate polarisation components, see eqs. 7a and 7b in Negus & Drain, 1981 (DOI 10.1088/0022-3727/15/3/003)
    %calculate flow velocity and Doppler shift per medium
    for n = 1:length(media.name)
        if media.scatterCoeff(n) ~= 0 && sum(photons.mediumNum == n) > 0
            S1(photons.mediumNum == n) = interp1(linspace(0,pi,length(media.phaseFunction{n}(1,:))),media.phaseFunction{n}(1,:),scatteringAngles(photons.mediumNum == n),'cubic');
            S2(photons.mediumNum == n) = interp1(linspace(0,pi,length(media.phaseFunction{n}(2,:))),media.phaseFunction{n}(2,:),scatteringAngles(photons.mediumNum == n),'cubic');
            scatterVector(:,photons.mediumNum == n) = media.refractiveIndex(n).*oldDirection(:,photons.mediumNum == n)./photons.wavelength - media.refractiveIndex(n).*photons.direction(:,photons.mediumNum == n)./photons.wavelength;
            DopplerShift(:,photons.mediumNum == n) = dot(-media.flowVelocity{n}(photons.position(1,photons.mediumNum == n),photons.position(2,photons.mediumNum == n),photons.position(3,photons.mediumNum == n)),scatterVector(:,photons.mediumNum == n));
        end
    end
    polarisationInPlane = cos(rotationAngles).*S1.*inPlaneScattered;
    polarisationOutOfPlane = sin(rotationAngles).*S2.*outOfPlaneScattered;

    %update polarisation 3: add components together and normalise
    photons.polarisationVector(:,scatteringPhotons) = polarisationInPlane(:,scatteringPhotons) + polarisationOutOfPlane(:,scatteringPhotons);
    photons.polarisationVector(:,scatteringPhotons) = photons.polarisationVector(:,scatteringPhotons)./vecnorm(photons.polarisationVector(:,scatteringPhotons));

    %update Doppler shift
    photons.DopplerShift(scatteringPhotons) = photons.DopplerShift(scatteringPhotons) + DopplerShift(scatteringPhotons);

    %Update scatter counter
    photons.timesScattered(scatteringPhotons) = photons.timesScattered(scatteringPhotons) + 1;
    photons.lastPositionScattered(:,scatteringPhotons) = photons.position(:,scatteringPhotons);
end
end