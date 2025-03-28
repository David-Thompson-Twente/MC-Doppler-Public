%%  makePencilSource.m
%   Wietske Verveld 13-1-2023, last update 19-3-2024
%   Creates a pencil source at the sourcePosition (no radius) pointing to 
%   the focusPosition
%
%   Inputs:     [3x1 double]    sourcePosition
%               [3x1 double]    focusPosition
%               [1x1 double]    wavelength
%
%   Optional:   [1x1 double]    lineWidth
%
%   Outputs:    [struct]        source, with properties: 
%                               source.opticalAxis      (3x1 double) 
%                               source.position         (3x1 double)    
%                               source.focusPosition    (3x1 double)    
%                               source.wavelength       (1x1 double)    
%                               source.lineWidth        (1x1 double)    
%                               source.sourceType       (char)   

%%
function source = makePencilSource(sourcePosition, focusPosition, wavelength,options)
arguments
    sourcePosition      (3,1) {mustBeNumeric}
    focusPosition       (3,1) {mustBeNumeric}
    wavelength          (1,1) {mustBeNumeric, mustBePositive}
    options.lineWidth   (1,1) {mustBeNumeric, mustBeNonnegative} = 0
end
source.opticalAxis = (focusPosition - sourcePosition)/vecnorm(focusPosition - sourcePosition);
source.position = sourcePosition;
source.focusPosition = focusPosition;
source.wavelength = wavelength;
source.lineWidth = options.lineWidth;
source.sourceType = 'pencil';
end