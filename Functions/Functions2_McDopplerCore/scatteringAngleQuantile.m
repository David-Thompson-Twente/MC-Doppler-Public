%%  scatteringAngleQuantile.m
%   David Thompson 26-04-2022, last update 19-3-2024
%   Obtains scattering angle based on a given phase function
%
%   Inputs:     [struct]        anglePickFunction, with relevant properties:
%                               anglePickFunction.scatterPickFunction
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
function scatteringAngles = scatteringAngleQuantile(anglePickFunction,photons,scatteringPhotons,media,classArgument)
arguments
    anglePickFunction   struct
    photons             struct
    scatteringPhotons   (1,:)   {mustBeNumericOrLogical(scatteringPhotons),mustBeInteger(scatteringPhotons),mustBeLessThanOrEqual(scatteringPhotons,1),mustBeNonnegative(scatteringPhotons)}
    media               struct
    classArgument       string
end
scatteringAngles = zeros(1,length(photons.alive),classArgument);
for n = 1:length(media.name)
    scatteringPhotonsMedium = scatteringPhotons & photons.alive == 1 & photons.mediumNum == n;
    numScatter = length(scatteringPhotonsMedium(scatteringPhotonsMedium == 1));
    if media.scatterCoeff(n) ~= 0
        scatteringAngles(1,scatteringPhotonsMedium) = rand(1,numScatter,classArgument)*2*pi;
        scatteringAngles(2,scatteringPhotonsMedium) = ppval(anglePickFunction(n).scatterPickFunction,rand(1,numScatter,classArgument));
        scatteringAngles(2,scatteringAngles(2,:) > pi) = pi;
        scatteringAngles(2,scatteringAngles(2,:) < 0) = 0;
    end
end
end