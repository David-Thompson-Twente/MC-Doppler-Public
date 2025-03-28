%%  compositePhaseFunction.m
%   David Thompson 19-01-2024, last update 13-03-2024
% 
%   Calculates generic composite phase function from fit to a measured or simulated PSD in terms of particle number
%
%   Inputs:     [1x1 double]    wavelength  [m]
%               [1x1 double]    mediumIndex
%               [2x1 double]    sizeLimits  (minimum and maximum particle radius in m)
%               [cfit]          PSDfit      (fit to measured PSD in terms of radius)
%               [string]        fittp       (indication if the PSD is by volume or number of particles)
%
%   Optional:   [1x1 double]    numAngles   (number of angles in phase function)
%               [1x1 double]    numFitPoints
%               [1x1 double]    particleDensity
%               [1x1 double]    particleIndex
%               [1x1 double]    particleVolumeFraction
%
%   Outputs:    [2xn double]    phaseFunctions
%               [1x1 double]    mu_s

%%
function [phaseFunctions, mu_s] = compositePhaseFunction(wavelength,mediumIndex,sizeLimits,PSDfit,options)
arguments
    wavelength                          (1,1) double %m
    mediumIndex                         (1,1) double
    sizeLimits                          (2,1) double
    PSDfit                              cfit
    options.numAngles                   (1,1) double = 1e3;
    options.numFitPoints                (1,1) double = 1e3;
    options.particleDensity             (1,1) double = 1; %kg/L
    options.particleIndex               (1,1) double = 1;
    options.particleVolumeFraction      (1,1) double = 0.01;
end
if ~exist('MatScat','dir')
    error('MatScat package required to use this function, add it to the path or download it at: https://nl.mathworks.com/matlabcentral/fileexchange/36831-matscat')
else
    addpath('MatScat','MatScat\bessel','MatScat\expcoeff','MatScat\test','MatScat\util')
end

% Light characteristics
k = 2*pi*mediumIndex/wavelength;

% PSD characteristics
particleRadiusRange = linspace(sizeLimits(1),sizeLimits(2),options.numFitPoints);
particleVolumePerSize = 4/3.*pi.*(particleRadiusRange).^3*1000;
PSD = PSDfit(particleRadiusRange)'.*particleVolumePerSize;

% Amount of particles per L
N_PSD = options.particleVolumeFraction./sum(PSD);       % number of PSD's to get total particle volume per L fluid
numberOfParticles = N_PSD.*PSD./particleVolumePerSize;  % array of number of particles per L fluid per diameter of the particles
N = length(particleRadiusRange);

% Calculate phase functions and scattering coefficients using MatScat
scatteringCrossSection = zeros(1,N);
I_perpendicular = zeros(N,options.numAngles);
I_parallel = zeros(N,options.numAngles);
for n = 1:N
    [S, C, ~] = calcmie(particleRadiusRange(n), options.particleIndex, mediumIndex, wavelength, options.numAngles);
    scatteringCrossSection(n) = C.sca;
    Int_parallel = abs(S(1,1,:)).^2;
    Int_perpendicular = abs(S(2,2,:)).^2;
    I_perpendicular(n,:) = Int_perpendicular./(scatteringCrossSection(n).*k^2);
    I_parallel(n,:) = Int_parallel./(scatteringCrossSection(n).*k^2);
end
pPolarised = sum(I_parallel.*scatteringCrossSection'.*numberOfParticles')./sum(scatteringCrossSection.*numberOfParticles);
sPolarised = sum(I_perpendicular.*scatteringCrossSection'.*numberOfParticles')./sum(scatteringCrossSection.*numberOfParticles);
phaseFunctions = [sPolarised;pPolarised];
mu_s = sum(numberOfParticles*1000.*scatteringCrossSection);
end