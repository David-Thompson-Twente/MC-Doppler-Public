%%  rotateRodrigues.m
%   David Thompson 30-06-2022, last update 05-03-2025
%   Rotates a vector or array of vectors around a given axis
%
%   Inputs:     [3xn double]    axisVector - axes to be rotated about (unit vector)
%               [1xn double]    elevationAngle - elevation angles of subject vector to axisvector
%               [1xn double]    azimuth - angle of rotation about axisVector
%               [3xn double]    vec_perpendicular - vector perpendicular to axisVector (reference for azimuth)
%               [string]        classArgument
%
%   Outputs:    [3xn double]    rotatedVector - subject vector rotated by azimuthAngle radians about axisVector

%% 
function rotatedVector = rotateRodrigues(axisVector, elevationAngle, azimuth, vec_perpendicular, classArgument)
arguments
    axisVector          (3,:)   {mustBeNumeric,mustBeNonNan} %should be unit vector
    elevationAngle      (1,:)   {mustBeNumeric}
    azimuth             (1,:)   {mustBeNumeric}
    vec_perpendicular   (3,:)   {mustBeNumeric,mustBeNonNan}
    classArgument       string  {mustBeMember(classArgument,{'single','double','gpuArray'})} %only here for the one use of rand, for the rest, ensure inputs are gpuArrays if you want to run on GPU
end
vec_parallel = cos(elevationAngle).*axisVector;
if isempty(vec_perpendicular)
    assistVector = axisVector+rand(3,1,classArgument);
    vec_perpendicular = cross(axisVector,assistVector,1);
    vec_perpendicular = sin(elevationAngle).*vec_perpendicular./vecnorm(vec_perpendicular);
    
    parallelnessCheck = dot(axisVector,assistVector);    
    if ~isempty(parallelnessCheck(parallelnessCheck == 1))
        warning([num2str(length(parallelnessCheck(parallelnessCheck == 1))) ' unit vector pairs are parallel'])
    end
else
    vec_perpendicular = sin(elevationAngle).*vec_perpendicular./vecnorm(vec_perpendicular);
end
subjectVector = vec_parallel + vec_perpendicular;

rotatedVector = vec_parallel + cos(azimuth).*vec_perpendicular + sin(azimuth).*cross(axisVector,subjectVector,1);
end