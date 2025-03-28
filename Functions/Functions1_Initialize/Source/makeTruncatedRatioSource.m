%%  makeTruncatedRatioSource.m
%   David Thompson 20-10-2023, last update 03-05-2024
%   Creates a source with an initial radius and a focal radius, where a
%   photon in the Gaussian source spot is projected towards the eventual 
%   position in the Gaussian focal spot, truncated at sourceRadius.
%   Focal spot standard deviation has a default value of sourceRadius but
%   can be tuned separately as well.
%
%   Inputs:     [3x1 double]    sourcePosition
%               [3x1 double]    focusPosition
%               [1x1 double]    sourceRadius [m]
%               [1x1 double]    focusRadius [m]
%               [1x1 double]    wavelength
%
%   Optional:   [1x1 double]    lineWidth
%               [1x1 double]    sourceStandardDeviation
%               [1x1 double]    coordNumPoints
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
function source = makeTruncatedRatioSource(sourcePosition, focusPosition, sourceRadius, focusRadius, wavelength, options)
arguments
    sourcePosition                  (3,1) {mustBeNumeric}
    focusPosition                   (3,1) {mustBeNumeric}
    sourceRadius                    (1,1) {mustBeNumeric, mustBePositive}
    focusRadius                     (1,1) {mustBeNumeric, mustBePositive}
    wavelength                      (1,1) {mustBeNumeric, mustBePositive}
    options.lineWidth               (1,1) {mustBeNumeric, mustBeNonnegative} = 0
    options.sourceStandardDeviation  (1,1) {mustBeNumeric, mustBeNonnegative} = sourceRadius
    options.coordNumPoints          (1,1) {mustBeNumeric, mustBeNonnegative} = 1e3
end
source.opticalAxis = (focusPosition - sourcePosition)/vecnorm(focusPosition - sourcePosition);
source.position = sourcePosition;
source.focusPosition = focusPosition;
source.radius = sourceRadius;
radialCoordinate = linspace(0,sourceRadius,options.coordNumPoints);
truncatedDistribution = radialCoordinate.*exp(-(radialCoordinate.^2./(2*options.sourceStandardDeviation.^2)));
truncatedCDF = cumtrapz(truncatedDistribution);
source.posPDF = csaps(truncatedCDF./max(truncatedCDF),radialCoordinate);
source.focusRatio = focusRadius/sourceRadius;
source.wavelength = wavelength;
source.lineWidth = options.lineWidth;
source.sourceType = 'truncated ratio';
end