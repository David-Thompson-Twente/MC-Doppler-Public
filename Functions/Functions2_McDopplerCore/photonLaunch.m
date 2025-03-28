%%  photonLaunch.m
%   David Thompson 30-03-2022, last updated 19-3-2024
%   Launches photon struct for specified source and sim geometry
%
%   Inputs:     [struct]        source, with relevant properties 
%                               source.position - [x,y,z] origin of source in [m]
%                               source.PDF -  probability distribution function for photon initial location 
%                               source.wavelength - source wavelength in [m]
%               [struct]        sim, with relevant properties:  
%                               sim.numPhotons - the number of photons to generate
%               [1x1 double]    killThreshold
%               [string]        classArgument
%
%   Optional:   [1x1 double]    polarisationAngle
%
%   Output:     [struct]        photons

%%
function photons = photonLaunch(source,sim,killThreshold,classArgument,options)
arguments
    source                      struct
    sim                         struct
    killThreshold               (1,1)   {mustBePositive(killThreshold),mustBeLessThan(killThreshold,1)}
    classArgument               string  {mustBeMember(classArgument,{'single','double','gpuArray'})}
    options.polarisationAngle   (1,1)   {mustBeNumeric} = 0 %in z-y plane
end
photons.amplitude = ones(1,sim.numPhotons,classArgument);
photons.position = photonPosition(source,sim,classArgument);
photons.direction = photonDirection(source,sim,photons.position,classArgument);
photons.wavelength = source.wavelength;
photons.alive = ones(1,sim.numPhotons,classArgument);
photons.mediumNum = ones(1,sim.numPhotons,classArgument);
photons.nextMedium = ones(1,sim.numPhotons,classArgument);
photons.pathLengths = zeros(1,sim.numPhotons,classArgument);
photons.DopplerShift = source.lineWidth.*randn(1,sim.numPhotons,classArgument);
photons.killThreshold = killThreshold;
photons.polarisationVector = photonPolarisation(sim,photons,options.polarisationAngle,classArgument);
photons.timesScattered = zeros(1,sim.numPhotons,classArgument);
photons.timesReflected = zeros(1,sim.numPhotons,classArgument);
photons.timesRefracted = zeros(1,sim.numPhotons,classArgument);
photons.likelihood = ones(1,sim.numPhotons,classArgument);
photons.distanceTravelled = zeros(1,sim.numPhotons,classArgument);
photons.opticalPathLength = zeros(1,sim.numPhotons,classArgument);
photons.lastPositionScattered = photons.position;
end