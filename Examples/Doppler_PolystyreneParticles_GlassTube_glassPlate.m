%% McDoppler_PolystyreneDryForm_GlassTube_glassPlate.m
% David Thompson 22-08-2024
% Script for simulating a polystyrene microsphere solution flowing through a glass tube, with an excitation fibre (SM600, Thorlabs) on one side, detector (DET10A2,Thorlabs) on the other
% Includes a brownian motion component in the flow velocity
%%
clear
reset(gpuDevice)

%% 1. Initialization

numAngles = 1e3;                                                %number of angles in phase function
desiredNumPhotons = 5e6;                                        %stop running once this many photons are detected
[~,commitNumber] = system('git -C MC-Doppler rev-parse HEAD');  %save git hash so you always know with which version you ran the simulation

% Define geometry
tubeOuterRadius = 0.675e-3;     %outer radius of flow tube [m]
tubeInnerRadius = 0.475e-3;     %inner radius of flow tube [m]
tubeRefractiveIndex = 1.5216;   %soda lime glass, refractiveindex.info (633 nm)
windowRefractiveIndex = 1.5151; %Fused Silica, refractiveindex.info (633 nm)
plateThickness = 100e-6;        %thickness of glass window in front of detector [m]
glassPlatePosZ = 1.8e-3;        %z-position of glass plate [m]

% Define source
wavelength = 1e-9.*633;                                                                             %[m]
sourceElevation = 0;                                                                                %[m]
sourceLateral = 0;                                                                                  %[m]
sourceFibreRotation = 0;                                                                            %[radians]
sourceAxialOffset = 5e-3;                                                                           %[m]
sourcePosition = [sourceElevation;sourceLateral;-tubeOuterRadius - sourceAxialOffset];              %[m]
focusPosition = [sourceElevation;sourceLateral;0];                                                  %[m]
focusPosition(1) = abs(sourcePosition(3) - focusPosition(3)).*tan(sourceFibreRotation);             %[m]
focusPosition(2) = 0;                                                                               %[m]
inputFibreNA = 2.5e-3;                                                                              %[dimensionless]
sourceOpeningHalfAngle = asin(inputFibreNA);                                                        %[radians]
sourceRadius = 170e-6;                                                                              %[m]
focusRadius = sourceRadius + vecnorm(focusPosition - sourcePosition)*tan(sourceOpeningHalfAngle);   %[m]
sourceLineWidth = 0;                                                                                %[m]

% Define detector
detectorRotation = 0;                                                                               %[radians]
detectorAxialOffset = 2.5e-3;                                                                       %[m]
detectorLateral = 0;                                                                                %[m]
detectorPosition = [0;detectorLateral;tubeOuterRadius + detectorAxialOffset];                       %[m]
detectorNormal = rotateRodrigues([0;0;-1],detectorRotation,0,[1;0;0],'gpuArray');                   %[unit vector, m]
detectorNormal = rotateRodrigues(detectorNormal,-detectorRotation,0,[0;1;0],'gpuArray');            %[unit vector, m]
detectorDimension = 0.5e-3;                                                                         %[m] 
detectorNA = 1;                                                                                     %[dimensionless]

%% 2. Build necessary structs
% Build source struct
source = makeTruncatedRatioSource(sourcePosition,focusPosition,sourceRadius,focusRadius,wavelength,'lineWidth',sourceLineWidth); %Truncated Gaussian spot, focussed

% Build sim struct
sim.shape = 'cylinder';                         %shape of simulation boundary
sim.numPhotons = 5e6;                           %number of photons per iteration
sim.radius = vecnorm(detectorPosition)+10e-3;   %radius of simulation space
sim.halfLength = 2.5e-3;                        %half of simulation space length
sim.GPU = 1;                                    %GPU-enabled, if not desired set to 0, remove gpuArray in what follows

% Build media struct
media.name                  = {'air','PVC','Polystyrene sphere suspension','glass plate'};
media.shape                 = {'cylinder','cylinder','cylinder','plane'};
media.surroundingMedium     = [0 1 2 1];                                                                                            %which medium is each medium contained within?
media.relevantDimension     = gpuArray([sim.radius, tubeOuterRadius, tubeInnerRadius plateThickness;sim.halfLength.*ones(1,4)]);    %medium dimension, radius and length for cylinders, thickness and lateral extent for planes
media.axis                  = gpuArray([0 0 0 0;1 1 1 0;0 0 0 -1]);                                                                 %for cylinders, define cylinder axis orientation
media.axisorigin            = gpuArray([zeros(3,3) [0;0;tubeOuterRadius + glassPlatePosZ + media.relevantDimension(1,4)/2]]);       %origin of cylinder axis for cylindrical media, center of plate for planes
media.scatterCoeff          = zeros(1,4,'gpuArray');                                                                                %scattering coefficient [m^-1]
media.absorptionCoeff       = gpuArray([0 1.4 0.3 1.4]);                                                                            %absorption coefficient [m^-1]
media.refractiveIndex       = gpuArray([1 tubeRefractiveIndex 1.3313 windowRefractiveIndex]);                                       %refractive index [dimensionless]
media.shiftX                = zeros(1,4,'gpuArray');                                                                                %lateral offset of media
media                       = cylinderSpatial(media,'gpuArray');                                                                    %generate functions describing media boundaries

%% 3. Specify samples
%For Brownian motion
particleDensity = 1050;                                             %[kg/m^3]
sampleTemperature = 300;                                            %[K]
fluidDensity = 1000;                                                %[kg/m^3]
kB = 1.3806e-23;                                                    %[J/K]
%Flow rate
flowrate_mL_per_min = 15;                                           %[mL/min]
flowRate = 1e-6.*flowrate_mL_per_min./60;                           %[m^3s^-1]
%Particle properties
particleRefractiveIndex = 1.58654;                                  %[dimensionless]
mu_s = 1e3.*[0.82 5.67 8.99 1.18 5.53 9.43 1.11 5.18 9.11];         %scattering coefficients of 9 different samples [m^-1]   
beadDiameters = 1e-6.*[1 1 1 1.5 1.5 1.5 3 3 3];                    %diameters of polystyrene beads in 9 different samples [m]

%% 4. Iterate over the samples

for sampleNumber = 1:length(beadDiameters) %iterate over all samples

    % Sample specifics
    particleDiameter = beadDiameters(sampleNumber);                     %[m]
    particleMass = particleDensity.*(4/3*pi*(particleDiameter/2)^3);    %[kg]
    fluidMass = fluidDensity.*(4/3*pi*(particleDiameter/2)^3);          %[kg]
    effectiveMass = particleMass + 0.5*fluidMass;                       %[kg], from Mo et al, 2015: https://doi.org/10.1364/OE.23.001888
    speedStandardDeviation = sqrt(kB*sampleTemperature/effectiveMass);  %For Brownian motion component
    media.flowVelocity(1:3)     = {@(x,y,z) [zeros(size(x),'gpuArray');zeros(size(y),'gpuArray');zeros(size(z),'gpuArray')]}; %for non-flowing media
    media.flowVelocity(3)       = {@(x,y,z) gpuArray([zeros(size(x));(2*flowRate)/(pi*tubeInnerRadius^2).*(1 - (sqrt(x.^2+z.^2)/tubeInnerRadius).^2);zeros(size(z))] + (speedStandardDeviation.*randn(size(x))).*randomUnitDirection(length(x)))}; %Poiseuille flow + Brownian motion [m/s]

    % Add phasefunctions to media struct
    [phaseFunction, scatterCrossSection] = simplePhaseFunction(particleDiameter/2,particleRefractiveIndex,gather(media.refractiveIndex(3)),source.wavelength,numAngles); %needs MatScat: https://nl.mathworks.com/matlabcentral/fileexchange/36831-matscat
    media.phaseFunction(3) = {phaseFunction};
    media.scatterCoeff(3) = mu_s(sampleNumber);
    [anglePickFunction,media] = mediaAnglePickFunction(media,'polarisation','polarised');

    % Run simulation
    detectedPhotons.timesScattered = []; %used to measure detected number of photons, needs to be empty at start of each run
    progressBar = waitbar(0, ['Detected 0/' num2str(desiredNumPhotons)],'OuterPosition',[820 430 280 80]);
    simText = ['Run ' num2str(sampleNumber) ' of ' num2str(length(beadDiameters)) 9 char(datetime)];
    disp(simText)
    n = 1;
    while length(detectedPhotons.timesScattered) < desiredNumPhotons %keep running simulations of sim.numPhotons photons, until desiredNumPhotons are detected
        waitbar(length(detectedPhotons.timesScattered)/desiredNumPhotons,progressBar,['Detected ' num2str(length(detectedPhotons.timesScattered)) '/' num2str(desiredNumPhotons)])
        simPhotons = gatherPhotonsGPU(McDopplerSim(media,sim,source,anglePickFunction,'plotVectors',0,'polarisation','polarised','polarisationAngle',pi/2,'saveTrajectories',0,'classArgument','gpuArray'));    %one run of MC-Doppler
        %optional decide which fields to save, omitFieldNames is given to circularDetector() as a list of field *not* to be saved, so any fields removed from omitFieldNames as below *is* saved
        omitFieldNames = fieldnames(simPhotons);
        omitFieldNames(ismember(omitFieldNames,'DopplerShift')) = [];
        omitFieldNames(ismember(omitFieldNames,'position')) = [];
        omitFieldNames(ismember(omitFieldNames,'direction')) = [];
        omitFieldNames(ismember(omitFieldNames,'opticalPathLength')) = [];
        omitFieldNames(ismember(omitFieldNames,'amplitude')) = [];
        omitFieldNames(ismember(omitFieldNames,'distanceTravelled')) = [];
        omitFieldNames(ismember(omitFieldNames,'timesScattered')) = [];
        omitFieldNames(ismember(omitFieldNames,'polarisationVector')) = [];
        photonsOnDetector = circularDetector(simPhotons,sim,media,detectorPosition,detectorDimension,'detectorNA',detectorNA,'detectorNormal',detectorNormal,'omitFields',omitFieldNames); %Which photons are detected
        if n == 1
            detectedPhotons = emptyDuplicateStruct(photonsOnDetector);
        end
        detectedPhotons = concatenatePhotonStruct(detectedPhotons,photonsOnDetector);
        n = n+1;
    end
    close(progressBar)
    disp(['Run ' num2str(sampleNumber) ' complete in ' num2str(n-1) ' iterations'])

    % Save results
    saveFolder = '';
    saveFileName = ['simResults_glassTube_' num2str(particleDiameter*1e6) 'micronParticles_mu_s_' num2str(media.scatterCoeff(3)*1e-3) 'perMM_uniformFlowIndexMatched.mat'];
    save([saveFolder saveFileName],'commitNumber','sim','source','detectedPhotons','flowrate_mL_per_min','detectorAxialOffset','-v7.3');
    clear detectedPhotons
end