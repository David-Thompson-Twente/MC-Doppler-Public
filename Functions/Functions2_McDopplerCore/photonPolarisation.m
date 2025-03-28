%%  photonPolarisation.m
%   David Thompson 03-10-2022, last update 19-3-2024
%   Generates polarisation unit vectors
%
%   Inputs:     [struct]        sim, with relevant property         sim.numPhotons
%               [struct]        photons, with relevant property     photons.direction
%               [1x1 double]    polarisationAngle
%               [string]        classArgument
%
%   Outputs:    [3xn double]    polarisation 

%%
function polarisation = photonPolarisation(sim,photons,polarisationAngle,classArgument)
arguments
    sim                 struct
    photons             struct
    polarisationAngle   (1,1) {mustBeNumeric}
    classArgument       string  {mustBeMember(classArgument,{'single','double','gpuArray'})}
end
polarisationReference = cross(photons.direction,[1;0;0].*ones(1,sim.numPhotons,classArgument),1);
polarisationReference = polarisationReference./vecnorm(polarisationReference);
polarisationReference(:,sum(polarisationReference) == 0) = [0;1;0]*ones(1,length(polarisationReference(:,sum(polarisationReference) == 0)),classArgument);
polarisation = rotateRodrigues(photons.direction,pi/2,polarisationAngle,polarisationReference,classArgument);
end