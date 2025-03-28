%%  reflectionCoeffs.m
%   David Thompson 12-04-2022, last update 19-3-2024
%   Determines reflection coefficients for all photons on interface, from
%   normal vector and direction unit vector
%
%   Inputs:     [struct]        photons, with relevant properties:  photons.direction
%                                                                   photons.mediumNum
%                                                                   photons.nextMedium
%                                                                   photons.alive
%               [struct]        media, with relevant properties:    media.refractiveIndex
%               [3xn double]    normalVectors - containing normal vectors to the interface
%               [1xn logical]   crossingPhotons
%               [string]        classArgument
%
%   Optional    [1x1 double]    tolerance
%
%   Outputs:    [2xn double]    reflectionCoefficients - giving a reflection 
%                                       coefficient per photon on the interface
%               [1xn double]    indexRatio - n2/n1 to later calculate refraction

%%
function [reflectionCoefficients,indexRatio] = reflectionCoeffs(photons,media,normalVectors,crossingPhotons,classArgument,options)
arguments
    photons             struct
    media               struct
    normalVectors       (3,:) {mustBeNumeric(normalVectors)}
    crossingPhotons     (1,:) {mustBeNumericOrLogical(crossingPhotons),mustBeInteger(crossingPhotons),mustBeLessThanOrEqual(crossingPhotons,1),mustBeNonnegative(crossingPhotons)}
    classArgument       string
    options.tolerance   (1,1) {mustBeNumeric,mustBePositive}    = 200*eps(1)
end
reflectionCoefficients = zeros(2,length(photons.alive),classArgument);
incidenceAngles = zeros(1,length(photons.alive),classArgument);

n1 = horzcat(media.refractiveIndex(photons.mediumNum));
n2 = horzcat(media.refractiveIndex(photons.nextMedium));
% clip dot product value to 1, as any larger values should not happen with unit vectors, and rounding errors lead to trouble using acos on the GPU
dotProduct = dot(-photons.direction(:,crossingPhotons),normalVectors(:,crossingPhotons),1);
dotProduct(abs(dotProduct) > 1) = 1;
incidenceAngles(crossingPhotons) = acos(dotProduct);
if ~isempty(incidenceAngles(incidenceAngles > pi/2))
    if sum(incidenceAngles(incidenceAngles > pi/2)) - sum(incidenceAngles > pi/2).*pi/2 < sum(incidenceAngles > pi/2)*options.tolerance
        incidenceAngles(incidenceAngles > pi/2) = pi/2;
    else
        error([num2str(sum(incidenceAngles > pi/2)) ' angles > pi/2 detected'])                
    end
end

R_s = abs((n1.*cos(incidenceAngles) - n2.*sqrt(complex(1-(n1./n2.*sin(incidenceAngles)).^2)))./...
    (n1.*cos(incidenceAngles) + n2.*sqrt(complex(1-(n1./n2.*sin(incidenceAngles)).^2)))).^2;

R_p = abs((n1.*sqrt(complex(1-(n1./n2.*sin(incidenceAngles)).^2)) - n2.*cos(incidenceAngles))./...
    (n1.*sqrt(complex(1-(n1./n2.*sin(incidenceAngles)).^2)) + n2.*cos(incidenceAngles))).^2;

R = [R_s ; R_p];
reflectionCoefficients(:,crossingPhotons) = R(:,crossingPhotons);
indexRatio = n2./n1;

if ~isempty(R(R - 1 > options.tolerance))
    numOffenders = length(R(R>1));
    plot(R_s)
    hold on
    plot(R_p)
    error(['R > 1, ' num2str(numOffenders) ' instances'])
end

end