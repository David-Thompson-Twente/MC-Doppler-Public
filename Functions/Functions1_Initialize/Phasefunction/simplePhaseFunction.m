%%  simplePhaseFunction.m
%   David Thompson 26-04-2022, last update 19-03-2025
%   Calculates normalized phase function using MatScat package for a simple monodisperse particle solution.
%
%   Inputs:     [1x1 double]        particleSize [m]
%               [1x1 double]        particleRefractiveIndex
%               [1x1 double]        mediumRefractiveIndex
%               [1x1 double]        wavelength [m]
%               [1x1 double]        numAngles, number of angles to sample phase function at
%
%   Outputs:    [2xn double]        phaseFunction array with n=numAngles
%               [struct]            scatterCrossSection

%%
function [Smatrix,scatterCrossSection] = simplePhaseFunction(particleSize,particleRefractiveIndex,mediumRefractiveIndex,wavelength,numAngles)
arguments
    particleSize (1,1) {mustBeNumeric(particleSize),mustBePositive(particleSize),mustBeLessThan(particleSize,100e-6)}
    particleRefractiveIndex (1,1) {mustBeNumeric(particleRefractiveIndex),mustBePositive(particleRefractiveIndex)}
    mediumRefractiveIndex (1,1) {mustBeNumeric(mediumRefractiveIndex),mustBePositive(mediumRefractiveIndex)}
    wavelength (1,1) {mustBeNumeric(wavelength),mustBePositive(wavelength),mustBeLessThanOrEqual(wavelength,2e-6)}    
    numAngles (1,1) {mustBeNumeric(numAngles),mustBePositive(numAngles)}
end
if ~exist('MatScat','dir')
    error('MatScat package required to use this function, add it to the path or download it at: https://nl.mathworks.com/matlabcentral/fileexchange/36831-matscat')
else
    addpath('MatScat','MatScat\bessel','MatScat\expcoeff','MatScat\test','MatScat\util')
end
k = 2*pi*mediumRefractiveIndex/wavelength;
[scatterMat,scatterCrossSection] = calcmie(particleSize,particleRefractiveIndex,mediumRefractiveIndex,wavelength,numAngles);
Smatrix = squeeze([abs(scatterMat(1,1,:).^2);abs(scatterMat(2,2,:)).^2])./(scatterCrossSection.sca.*k^2);
end