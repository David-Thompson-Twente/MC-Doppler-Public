%%  checkPhotonEscape.m
%   David Thompson 04-04-2022, last update 13-3-2024
%
%   Checks if photons are still within the simulation ROI after propagation.
%
%   Inputs:     [struct]        photons with relevant properties:   photons.position [m]
%                                                                   photons.pathLengths [m]
%                                                                   photons.direction
%                                                                   photons.alive (binary array)
%               [struct]        sim with relevant properties:       sim.shape
%                                                                   sim.radius      (for sphere and cylinder shape)
%                                                                   sim.halfLength  (for cylinder and cube shape)
%                                                                   sim.halfWidth   (for cube shape)
%                                                                   sim.halfDepth   (for cube shape)
%
%   Outputs:    [struct]        photons, with photons.alive updated to give zeros for photons which will escape the simulation in the next propagation step.

%%
function photons = photonEscapeCheck(photons,sim)
arguments
    photons struct
    sim struct
end
switch sim.shape
    case 'sphere'
        escapees = vecnorm(photons.position + photons.pathLengths.*photons.direction) >= sim.radius;
    case 'cylinder'
        propagatedPos = photons.position + photons.pathLengths.*photons.direction;
        escapees = vecnorm(propagatedPos([1,3],:)) >= sim.radius | abs(propagatedPos(2,:)) >= sim.halfLength;
    case 'cube'
        propagatedPos = photons.position + photons.pathLengths.*photons.direction;
        escapees_x = abs(propagatedPos(1,:)) >= sim.halfWidth;
        escapees_y = abs(propagatedPos(2,:)) >= sim.halfLength;
        escapees_z = abs(propagatedPos(3,:)) >= sim.halfDepth;
        escapees = max([escapees_x;escapees_y;escapees_z]);
end
photons.alive(escapees) = 0;
end