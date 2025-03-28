%%  plotPositions.m
%   Wietske Verveld 19-3-2024, last update 19-3-2024
%
%   Plots positions during McDopplerSim
%
%   Inputs:   [struct]      photons, with relevant property photons.position
%                                                           photons.alive
%             [struct]      sim,     with relevant property sim.radius
%                                                           sim.halfLength
%             [handle]      pointPlot
%   
%   Optional: [1x1 double]  photonPlotLimit

%%
function plotPositions(photons,sim,pointPlot,options)
arguments
    photons         struct
    sim             struct
    pointPlot       handle
    options.photonPlotLimit         (1,1)   {mustBeNumeric}             = 5000
end
    
if length(photons.alive) > options.photonPlotLimit
    photonSubSample = 1:options.photonPlotLimit;
else
    photonSubSample = 1:length(photons.alive);
end

figure(pointPlot)
hold on
plot3(photons.position(1,photonSubSample),photons.position(2,photonSubSample),photons.position(3,photonSubSample),'.')
view(45,45)
daspect([1 1 1]);

xlabel('x (m)');ylabel('y (m)');zlabel('z (m)');
xlim(2*[-sim.radius,sim.radius]);
ylim(2*[-sim.halfLength,sim.halfLength]);
zlim([-sim.radius,sim.radius]);
shg
drawnow

end