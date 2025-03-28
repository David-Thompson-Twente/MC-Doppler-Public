%%  MediaAnglePickFunction.m
%   Wietske Verveld 20-07-2023, last update 19-03-2025
%   Overarching function to obtain anglePickFunction for any polarisation case
%
%   Inputs:     [struct]        media, with relevant properties:    media.name
%                                                                   media.scatterCoeff
%                                                                   media.phaseFunction
%   Optional:   [string]        polarisation
%               [1x1 double]    polarisationAngle
%               [1x1 double]    numAzimuth
%
%   Outputs:    [struct]        anglePickFunction (fields depend on polarization)
%               [struct]        media, with updated properties  media.S1
%                                                               media.S2

%%
function [anglePickFunction,media] = mediaAnglePickFunction(media,options)
arguments
    media                           struct
    options.polarisation            string  {mustBeMember(options.polarisation,{'unpolarised','polarised'})} = 'unpolarised'
    options.polarisationAngle       double = pi/2
    options.numAzimuth              (1,1)   {mustBeInteger}             = 1e3
end

for n = 1:length(media.name)
    if media.scatterCoeff(n) ~= 0
        theta = linspace(0,pi,size(media.phaseFunction{n},2));
        switch options.polarisation
            case 'unpolarised'
                phaseFunction = ((media.phaseFunction{n}(1,:) + media.phaseFunction{n}(2,:))./2).*sin(theta);
                anglePickFunction(n) = scatterInverseCDF(phaseFunction);
            case 'polarised'
                phi = linspace(0,2*pi,options.numAzimuth);
                phaseFunction2D = (cos(phi).^2.*media.phaseFunction{n}(1,:)' + sin(phi).^2.*media.phaseFunction{n}(2,:)').*sin(theta)';
                anglePickFunction(n) = scatterInverseCDFpolarised(phaseFunction2D);
        end
        media.S1(n) = csaps(linspace(0,pi,length(media.phaseFunction{n}(1,:))),media.phaseFunction{n}(1,:));
        media.S2(n) = csaps(linspace(0,pi,length(media.phaseFunction{n}(2,:))),media.phaseFunction{n}(2,:));
    else
        switch options.polarisation
            case 'unpolarised'
                anglePickFunction(n) = struct('scatterPickFunction',[],'phaseFunction',[]);
            case 'polarised'
                anglePickFunction(n) = struct('azimuthPickFunction',[],'scatterPickMatrix',[]);
        end
    end
end
end