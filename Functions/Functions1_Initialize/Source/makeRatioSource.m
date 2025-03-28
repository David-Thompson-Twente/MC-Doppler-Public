%%  makeRatioSource.m
%   David Thompson 01-11-2022, last update 28-3-2024
%   Creates a source with an initial radius and a focal radius, where a
%   photon in the Gaussian source spot is projected towards the eventual 
%   position in the Gaussian focal spot
%
%   Inputs:     [3x1 double]    sourcePosition
%               [3x1 double]    focusPosition
%               [1x1 double]    sourceRadius [m]
%               [1x1 double]    focusRadius [m]
%               [1x1 double]    wavelength
%
%   Optional:   [1x1 double]    lineWidth
%
%   Outputs:    [struct]        source, with properties: 
%                               source.opticalAxis      (3x1 double) 
%                               source.position         (3x1 double)    
%                               source.focusPosition    (3x1 double)    
%                               source.radius           (1x1 double)
%                               source.posPDF           (function)      
%                               source.focusRatio       (1x1 double)      
%                               source.wavelength       (1x1 double)    
%                               source.lineWidth        (1x1 double)    
%                               source.sourceType       (char)   

%%
function source = makeRatioSource(sourcePosition, focusPosition, sourceRadius, focusRadius, wavelength, options)
arguments
    sourcePosition      (3,1) {mustBeNumeric}
    focusPosition       (3,1) {mustBeNumeric}
    sourceRadius        (1,1) {mustBeNumeric, mustBePositive}
    focusRadius         (1,1) {mustBeNumeric, mustBePositive}
    wavelength          (1,1) {mustBeNumeric, mustBePositive}
    options.lineWidth   (1,1) {mustBeNumeric, mustBeNonnegative} = 0
end
source.opticalAxis = (focusPosition - sourcePosition)/vecnorm(focusPosition - sourcePosition);
source.position = sourcePosition;
source.focusPosition = focusPosition;
source.radius = sourceRadius;
source.posPDF = @(numPhotons) sourceRadius.*randn(1,numPhotons);
source.focusRatio = focusRadius/sourceRadius;
source.wavelength = wavelength;
source.lineWidth = options.lineWidth;
source.sourceType = 'ratio';
end