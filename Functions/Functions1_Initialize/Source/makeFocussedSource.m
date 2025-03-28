%%  makeFocussedSource.m
%   David Thompson 26-10-2022, last update 28-3-2024
%   Creates a focussed Gaussian beam source object as input for launchPhoton.m
%
%   From focal length, focus position, input pupil size and initial photon
%   position a diverging or converging beam can be generated
%
%   Inputs:     [3x1 double]    sourcePosition
%               [1x1 double]    refractiveIndex
%               [1x1 double]    focalLength [m]
%               [3x1 double]    focusPosition [m] distance from
%               [1x1 double]    pupilDiameter
%               [1x1 double]    wavelength
%
%   Optional:   [1x1 double]    lineWidth
%
%   Outputs:    [struct]        source, with properties: 
%                               source.opticalAxis      (3x1 double) 
%                               source.position         (3x1 double)    
%                               source.focusPosition    (3x1 double)    
%                               source.posPDF           (function)      
%                               source.wavelength       (1x1 double)    
%                               source.lineWidth        (1x1 double)    
%                               source.sourceType       (char)   

%%
function source = makeFocussedSource(sourcePosition,refractiveIndex,focalLength,focusPosition,pupilDiameter,wavelength,options)
arguments
    sourcePosition      (3,1)   {mustBeNumeric(sourcePosition)}
    refractiveIndex     (1,1)   {mustBeNumeric, mustBePositive}
    focalLength         (1,1)   {mustBeNumeric}
    focusPosition       (3,1)   {mustBeNumeric}
    pupilDiameter       (1,1)   {mustBeNumeric, mustBePositive}
    wavelength          (1,1)   {mustBeNumeric(wavelength),mustBePositive(wavelength),mustBeLessThanOrEqual(wavelength,2e-6)}    
    options.lineWidth   (1,1)   {mustBeNumeric, mustBeNonnegative} = 0
end
NA = refractiveIndex*pupilDiameter/(2*sqrt(focalLength^2 + (pupilDiameter/2)^2));
beamWaist = wavelength/(pi*NA);
RayleighRange = pi*beamWaist^2*refractiveIndex/wavelength;
r = vecnorm(sourcePosition - focusPosition);
spotRadius = beamWaist*sqrt(1+(r/RayleighRange)^2);

source.opticalAxis = sign(focalLength)*(focusPosition - sourcePosition)/vecnorm(focusPosition - sourcePosition);
source.position = sourcePosition;
source.focusPosition = focusPosition;
source.posPDF = @(numPhotons) spotRadius/2.*randn(1,numPhotons);
source.wavelength = wavelength;
source.lineWidth = options.lineWidth;
source.sourceType = 'focussed';
end