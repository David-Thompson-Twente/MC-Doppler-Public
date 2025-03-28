%%  photonKill.m
%   David Thompson 29-04-2022, last update 05-03-2025
%   Selects at random which photons with below-threshold amplitudes are
%   killed and which survive another iteration
%
%   Inputs:   [struct]      photons, with relevant properties   photons.amplitude
%                                                               photons.killThreshold
%                                                               photons.alive
%             [1x1 double]  rouletteAmplitude, amplitude non-absorbed photons are given, default 0.5              
%
%   Outputs:  [struct]      photons, with updated properties    photons.amplitude
%                                                               photons.alive

%%
function photons = photonKill(photons,classArgument,options)
arguments
    photons                     struct
    classArgument               string
    options.rouletteAmplitude   (1,1)   {mustBeNumeric, mustBePositive} = 0.5
end
lowAmplitudePhotons = photons.amplitude <= photons.killThreshold & photons.amplitude ~=0 & photons.alive == 1;
roulette = rand(1,length(lowAmplitudePhotons),classArgument);

killCriterion = lowAmplitudePhotons & roulette > photons.killThreshold;

photons.amplitude(killCriterion) = 0;
photons.alive(killCriterion) = 0;

photons.amplitude(~killCriterion & lowAmplitudePhotons) = options.rouletteAmplitude;
end