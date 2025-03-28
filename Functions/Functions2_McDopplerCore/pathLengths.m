%%  pathLengths.m
%   David Thompson 30-03-2022, last update 19-3-2024
%   Evaluates quantile function of photon scattering path length probability distribution -ln(1-x)/µ_s
%
%   Input:      [struct]        media, with relevant property   media.scatterCoeff
%               [struct]        photons, with relevant property photons.alive
%                                                               photons.pathLengths
%                                                               photons.mediumNum
%               [string]        classArgument
%
%   Output:     [struct]        photons with updated property pathLengths
%%
function photons = pathLengths(media,photons,classArgument)
arguments
    media           struct
    photons         struct
    classArgument   string  
end
randomNumber = rand(1,length(photons.alive),classArgument);
scatterCoeff = horzcat(media.scatterCoeff(photons.mediumNum));
photons.pathLengths(photons.pathLengths == 0) = -log(1-randomNumber(photons.pathLengths == 0))./scatterCoeff(photons.pathLengths == 0);
photons.pathLengths(photons.pathLengths == Inf) = realmax;
end