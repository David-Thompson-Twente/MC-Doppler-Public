%%  scatteringAnglePolarised.m
%   David Thompson 28-06-2022, Last updated 19-3-2024
%   Obtains scattering azimuth and elevation angles based on
%   polarisation-dependent phase function
%
%   Inputs:     [struct]        anglePickFunction, with relevant properties:
%                               anglePickFunction.azimuthPickFunction - quantile function for picking scattering azimuth
%                               anglePickFunction.scatterPickMatrix - matrix of 2D quantile function as function of azimuth
%               [struct]        photons, with relevant properties:
%                               photons.alive
%                               photons.mediumNum
%               [1xn logical]   scatteringPhotons
%               [struct]        media, with relevant properties:
%                               media.name
%                               media.scatterCoeff
%               [string]        classArgument
%
%   Outputs:    [2xn double]    scatteringAngles, first row corresponding
%                                    to azimuth, second row to elevation

%%
function scatteringAngles = scatteringAnglePolarised(anglePickFunction,photons,scatteringPhotons,media,classArgument)
arguments
    anglePickFunction   struct
    photons             struct
    scatteringPhotons   (1,:)   {mustBeNumericOrLogical(scatteringPhotons),mustBeInteger(scatteringPhotons),mustBeLessThanOrEqual(scatteringPhotons,1),mustBeNonnegative(scatteringPhotons)}
    media               struct
    classArgument       string
end
scatteringAngles = zeros(2,length(photons.alive),classArgument);
for n = 1:length(media.name)
    scatteringPhotonsMedium = scatteringPhotons & photons.alive == 1 & photons.mediumNum == n;
    numScatter = length(scatteringPhotonsMedium(scatteringPhotonsMedium == 1));
    if media.scatterCoeff(n) ~= 0        
        scatteringAngles(1,scatteringPhotonsMedium) = interp1(anglePickFunction(n).azimuthPickFunction,linspace(0,2*pi,size(anglePickFunction(n).azimuthPickFunction,2)),rand(1,numScatter,classArgument),'spline');
        scatteringAngles(1,scatteringAngles(1,:) < 0) = 0;
        scatteringAngles(2,scatteringPhotonsMedium) = interp2(linspace(0,1,size(anglePickFunction(n).scatterPickMatrix,2)),linspace(0,2*pi,size(anglePickFunction(n).scatterPickMatrix,1)),anglePickFunction(n).scatterPickMatrix,rand(1,numScatter,classArgument),scatteringAngles(1,scatteringPhotonsMedium),'cubic');
        scatteringAngles(2,scatteringAngles(2,:) > pi) = pi;
        scatteringAngles(2,scatteringAngles(2,:) < 0) = 0;
    end
end
end