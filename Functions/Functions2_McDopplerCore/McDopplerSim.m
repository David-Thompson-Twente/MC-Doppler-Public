%%  McDopplerSim.m
%   David Thompson 17-05-2022, last update 19-03-2025
%   Function to run a full McDoppler simulation
%
%   Inputs:     [struct]            media, containing n media defined from the outside in,
%                                   media have properties:  
%                                   media(n).name [char]
%                                   media(n).shape [char] cylinder/sphere
%                                   media(n).relevantDimension [double]
%                                   media(n).spatial function handle
%                                   defining medium geometry
%                                   media(n).scatterCoeff [double]
%                                   media(n).absorptionCoeff [double]
%                                   media(n).refractiveIndex [double]
%                                   media(n).flowVelocity function handle defining flow profile
%                                   media(n).surroundingMedium [double]
%
%               [struct]            sim, defining basic properties of the simulation:
%                                   sim.shape [char] sphere/cylinder
%                                   sim.numPhotons [double]
%                                   sim.dimensions [double] radius for sphere, 
%                                                  [1x2 double] [radius, length] for cylinder

%               [struct]            source, defining the source properties:
%                                   source.position [3x1 double]
%                                   source.wavelength [double]
%                                   source.maxAngle [double] max divergence angle of source
%                                   source.opticalAxis [3x1 double]
%                                   source.radiusPDF [function handle]
%                                   source.anglePDF [prob.UniformDistribution]
%
% Optional:     [1x1 double]        killThreshold
%               [1x1 double]        photonPlotLimit
%               [logical]           plotVectors
%               [logical]           plotPolarisation
%               [logical]           plotInteractPositions
%               [logical]           plotGIF
%               [string]            GIFname
%               [logical]           showProgressBar
%               [logical]           indicateMemoryUse
%               [logical]           saveTrajectories
%               [string]            polarisation
%               [1x1 double]        polarisationAngle
%               [string]            classArgument
%
%   Outputs:    [struct]            detectedPhotons, containing positional, directional
%                                   properties of photons escaped from the medium or absorbed, 
%                                   as well as Doppler shift properties and amplitude.

%%
function detectedPhotons = McDopplerSim(media,sim,source,anglePickFunction,options)
arguments
    media                           struct
    sim                             struct
    source                          struct
    anglePickFunction               struct
    options.killThreshold           (1,1)   {mustBeNumeric}             = 0.1
    options.photonPlotLimit         (1,1)   {mustBeNumeric}             = 5000
    options.plotVectors             (1,1)   {mustBeNumericOrLogical}    = false
    options.plotPolarisation        (1,1)   {mustBeNumericOrLogical}    = false
    options.plotInteractPositions   (1,1)   {mustBeNumericOrLogical}    = false
    options.plotGIF                 (1,1)   {mustBeNumericOrLogical}    = false
    options.GIFname                 char                                = ['GIF_' datestr(now,'dd_mm_HH_MM_SS') '.gif']
    options.showProgressBar         (1,1)   {mustBeNumericOrLogical}    = 1
    options.indicateMemoryUse       (1,1)   {mustBeNumericOrLogical}    = 0
    options.saveTrajectories        (1,1)   {mustBeNumericOrLogical}    = 0
    options.polarisation            string  {mustBeMember(options.polarisation,{'unpolarised','polarised'})} = 'unpolarised'
    options.polarisationAngle       (1,1)   {mustBeNumeric}             = 0
    options.classArgument           string  {mustBeMember(options.classArgument,{'single','double','gpuArray'})} = 'double'
end

if options.plotVectors == 1
    vectorPlot = figure;
end
if options.plotInteractPositions == 1
    pointPlot = figure;
end
if ~strcmp(class(media.scatterCoeff),options.classArgument)
    warning(['One or more fields in the media struct are of class: ' class(media.scatterCoeff) ' this does not match the indicated classArgument: ' options.classArgument ', this could negatively impact performance.'])
end
if options.showProgressBar
    progressBar = waitbar(0, 'Running McDoppler, fraction of photons escaped:');
end

photons = photonLaunch(source,sim,options.killThreshold,options.classArgument,'polarisation',options.polarisation,'polarisationAngle',options.polarisationAngle);
photons = currentMedium(photons,media);
if options.saveTrajectories
    photons.positionHistory = reshape(gather(photons.position),3,1,[]);
    photons.ID = linspace(1,sim.numPhotons,sim.numPhotons);
end
detectedPhotons = emptyDuplicateStruct(photons,'omitFields',{'wavelength','killThreshold','alive','positionHistory','directionHistory'});

while max(photons.alive) == 1    
    oldPositions = photons.position; 
    polVec = photons.polarisationVector;
    [photons, normalVectors] = updatePositions(photons,media,sim,options.classArgument);
    photons = updateDirections(photons,media,normalVectors,anglePickFunction,options.polarisation,options.classArgument);
    photons = photonKill(photons,options.classArgument);
    propVec2 = photons.position - oldPositions;
    polVec = polVec./(1./vecnorm(propVec2));

    if options.saveTrajectories
        photons.positionHistory = [photons.positionHistory photons.positionHistory(:,end,:)];        
        photons.positionHistory(:,end,photons.ID) = photons.position(:,photons.ID == photons.ID);
    end

    if options.plotVectors
        plotVectors(photons,sim,vectorPlot,oldPositions,propVec2,polVec,'photonPlotLimit',options.photonPlotLimit,'plotPolarisation',options.plotPolarisation,'plotGIF',options.plotGIF,'GIFname',options.GIFname);
    end
    if options.plotInteractPositions
        plotPositions(photons,sim,pointPlot,'photonPlotLimit',options.photonPlotLimit);
    end

    [photons, detectedPhotons] = photonDetect(photons,detectedPhotons);

    if options.showProgressBar
        photonFraction = (sim.numPhotons - length(photons.alive))./sim.numPhotons;
        waitbar(photonFraction,progressBar,['Running McDoppler, ' num2str(sim.numPhotons - length(photons.alive)) '/' num2str(sim.numPhotons) ' photons escaped'])
    end
    if options.indicateMemoryUse
        disp(memory)
    end
    if max(isnan(photons.position(:))) == 1
        error('NaN-valued positions found')
    end
    if max(abs(imag(photons.position(:)))) > 0
        numComplex = length(photons.position(max(abs(imag(photons.position(:)))) > 0));
        error([num2str(numComplex) ' complex-valued positions detected']);
    end
end

if options.saveTrajectories
    detectedPhotons.positionHistory = photons.positionHistory(:,:,detectedPhotons.ID);
end
if options.showProgressBar
    close(progressBar)
end

end
