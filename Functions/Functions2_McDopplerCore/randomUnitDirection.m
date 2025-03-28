%%  randomUnitDirection.m
%   David Thompson 26-08-2024
%   Creates uniformly directed unit vectors in 3D space
%
%   Inputs: [1x1 double] numVectors
%
%   Outputs [3 x numVectors double] randomDirections
%%
function randomDirections = randomUnitDirection(numVectors)
arguments
    numVectors  (1,1) {mustBeNumeric,mustBePositive}
end
randomDirections = randn(3,numVectors);
randomDirections = randomDirections./vecnorm(randomDirections);
end