%%  plotVectors.m
%   Wietske Verveld 19-3-2024, last update 19-3-2024
%
%   Plots vectors during McDopplerSim
%
%   Inputs:     [struct]      photons, with relevant property photons.alive
%               [struct]      sim,     with relevant property sim.radius
%                                                           sim.halfLength
%               [handle]      vectorPlot
%               [3xn double]  oldPositions
%               [3xn double]  propVec2
%               [3xn double]  polVec
%
%   Optional:   [1x1 double]  photonPlotLimit
%               [logical]     plotPolarisation
%               [logical]     plotGIF
%               [logical]     GIFname

%%
function plotVectors(photons,sim,vectorPlot,oldPositions,propVec2,polVec,options)
arguments
    photons         struct
    sim             struct
    vectorPlot      handle
    oldPositions    (3,:)       {mustBeNumeric}
    propVec2        (3,:)       {mustBeNumeric}
    polVec          (3,:)       {mustBeNumeric}
    options.photonPlotLimit         (1,1)   {mustBeNumeric}             = 5000
    options.plotPolarisation        (1,1)   {mustBeNumericOrLogical}    = false
    options.plotGIF                 (1,1)   {mustBeNumericOrLogical}    = false
    options.GIFname                 char                                = ['GIF_' datestr(now,'dd_mm_HH_MM_SS') '.gif']
end

if length(photons.alive) > options.photonPlotLimit
    photonSubSample = 1:options.photonPlotLimit;
else
    photonSubSample = 1:length(photons.alive);
end

figure(vectorPlot)
hold on
quiver3(oldPositions(1,photonSubSample),oldPositions(2,photonSubSample),oldPositions(3,photonSubSample)...
    ,propVec2(1,photonSubSample),propVec2(2,photonSubSample),propVec2(3,photonSubSample),'AutoScale','off','LineWidth',1.5);
if options.plotPolarisation
    quiver3(oldPositions(1,photonSubSample),oldPositions(2,photonSubSample),oldPositions(3,photonSubSample)...
        ,polVec(1,photonSubSample),polVec(2,photonSubSample),polVec(3,photonSubSample),'AutoScale','off','Color','green','LineWidth',1.5);
end
view(45,45)
daspect([1 1 1]);
xlabel('x (m)');ylabel('y (m)');zlabel('z (m)');
xlim(2*[-sim.radius,sim.radius]);
ylim(2*[-sim.radius,sim.radius]);
ylim(2*[-sim.halfLength,sim.halfLength]); % why is this line double?
zlim([-sim.radius,sim.radius]);
grid on
shg
drawnow

if options.plotGIF
    if ~exist('gif.m','file')
        error('You require the file gif.m on the path: https://nl.mathworks.com/matlabcentral/fileexchange/63239-gif')
    end
    figure(vectorPlot)
    if ~exist(options.GIFname,'file')
        gif(options.GIFname,'DelayTime',0.25,'LoopCount',1);
    else
        gif
    end
end

end