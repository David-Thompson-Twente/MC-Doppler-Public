%%  scatterInverseCDFpolarised.m
%   David Thompson 17-06-2022, last update 19-03-2025
%   Generates bivariate quantile function for scatering angle of polarised light
%
%   Inputs:     [2xn double]    phaseFunction2D - matrix of P(\phi,\theta)
%
%   Outputs:    [struct]        anglePickFunction, with the properties:
%                               anglePickFunction.azimuthPickFunction - struct containing fit to quantile function for picking scattering azimuth
%                               anglePickFunction.scatterPickMatrix - undersampled matrix of azimuth-dependent quantile functions for elevation

%%
function anglePickFunction = scatterInverseCDFpolarised(phaseFunction2D)
arguments
   phaseFunction2D {mustBeNumeric,mustBeNonnegative}
end
phi = linspace(0,2*pi,size(phaseFunction2D,2));
theta = linspace(0,pi,size(phaseFunction2D,1));
azimuthPDF = sum(phaseFunction2D);
azimuthCDF = cumtrapz(azimuthPDF);
azimuthPickFunction = azimuthCDF/max(azimuthCDF);

scatterPickMatrix = zeros(length(phi),length(theta));

if any(azimuthPDF==0)
    error('0 in azimuthPDF will give errors')
end

for n = 1:length(phi)
    CDF = cumtrapz(phaseFunction2D(:,n)./azimuthPDF(n));
    scatterPickMatrix(n,:) = csaps(CDF/max(CDF),theta,[],linspace(0,1,length(theta)));
end
anglePickFunction.azimuthPickFunction = azimuthPickFunction;
anglePickFunction.scatterPickMatrix = scatterPickMatrix;
end