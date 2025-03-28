%%  currentMedium.m
%   David Thompson 06-04-2022, last update 13-3-2024
%
%   Checks current medium for all photons, updates optical properties as needed for propagation.
%
%   Inputs:   [struct]      photons, with relevant property photons.position
%             [struct]      media,   with relevant property media.spatial                                    
%
%   Outputs:  [struct]      photons, with updated property photons.mediumNum

%%
function photons = currentMedium(photons,media)
arguments
    photons struct
    media   struct
end
x = photons.position(1,:);
y = photons.position(2,:);
z = photons.position(3,:);
for n = 1:length(media.spatial)
    photons.mediumNum(media.spatial{n}(x,y,z)) = n;
end
end