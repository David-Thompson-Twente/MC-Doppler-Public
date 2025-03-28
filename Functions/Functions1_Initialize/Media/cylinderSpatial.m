%%  cylinderSpatial.m
%   Wietske Verveld 23-06-2023, last update 27-03-2024
%
%   Formulate the spatial function for a cylinder with an axis through an axis origin
%   Based on https://math.stackexchange.com/questions/3518495/check-if-a-general-point-is-inside-a-given-cylinder
%
%   Inputs:     [struct]        media, with relevant properties:    media.name
%                                                                   media.shape
%                                                                   media.spatial
%                                                                   media.relevantDimension
%                                                                   media.axis
%                                                                   media.axisorigin
%               [string]        classArgument
%
%   Outputs:    [struct]        media, with property media.spatial for cylindrical shape

%%
function media = cylinderSpatial(media,classArgument)
arguments
    media           struct
    classArgument   string
end

for n_medium = 1:length(media.name)%try to get rid of loop
    if strcmp(media.shape{n_medium},'cylinder')
        media.spatial(n_medium) = {@(x,y,z) ...
            vecnorm(cross(media.axis(:,n_medium).*ones(1,length(x),classArgument),([x;y;z]-media.axisorigin(:,n_medium))))./vecnorm(media.axis(:,n_medium)) <= media.relevantDimension(1,n_medium) & ...  % check if distance to axis is less than cylinder radius
            abs(y) <= media.relevantDimension(2,n_medium)}; % check y distance
    end
end

