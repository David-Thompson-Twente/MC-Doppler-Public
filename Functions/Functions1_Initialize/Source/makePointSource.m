%%  makePointSource.m
%   David Thompson 14-10-2022, last update 28-3-2024
%   Creates a simple point source at sourcePosition
%
%   Inputs:     [3x1 double]    sourcePosition
%               [1x1 double]    wavelength
%
%   Optional:   [1x1 double]    lineWidth
%
%   Outputs:    [struct]        source, with properties: 
%                               source.opticalAxis      (3x1 double) 
%                               source.position         (3x1 double)    
%                               source.wavelength       (1x1 double)    
%                               source.lineWidth        (1x1 double)    
%                               source.sourceType       (char)   

%%
function source = makePointSource(sourcePosition,wavelength,options)
arguments
    sourcePosition      (3,1) {mustBeNumeric}
    wavelength          (1,1) {mustBeNumeric, mustBePositive} 
    options.lineWidth   (1,1) {mustBeNumeric, mustBeNonnegative} = 0
end
source.opticalAxis = [0;0;1];
source.position = sourcePosition;
source.wavelength = wavelength;
source.lineWidth = options.lineWidth;
source.sourceType = 'point';
end