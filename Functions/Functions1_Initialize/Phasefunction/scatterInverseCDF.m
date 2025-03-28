%%  scatterInverseCDF.m
%   David Thompson 26-04-2022, last update 19-03-2025
%   Creates inverse cumulative distribution function from phase function
% 
%   Inputs:   [1xn double]      phaseFunction
%
%   Outputs:  [struct]          anglePickFunction, with relevant fields
%                               anglePickFunction.scatterPickFunction (cubic spline object)
%                               anglePickFunction.phaseFunction (cubic spline object)

%%
function anglePickFunction = scatterInverseCDF(phaseFunction)
arguments
    phaseFunction (:,1) {mustBeNumeric(phaseFunction),mustBeNonnegative(phaseFunction)}
end
angles = linspace(0,pi,length(phaseFunction));
CDF = cumtrapz(phaseFunction);
anglePickFunction.scatterPickFunction = csaps(CDF./max(CDF),angles);

anglePickFunction.phaseFunction = csaps(angles,phaseFunction./mean(phaseFunction));
end