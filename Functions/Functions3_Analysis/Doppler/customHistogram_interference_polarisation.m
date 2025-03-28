%%  customHistogram_interference_polarisation.m
%   Brute force calculation of Doppler power spectra based on the interference between the detected photons from a McDoppler simulation.
%%
function bruteForceOutcomesPol = customHistogram_interference_polarisation(detectedPhotons, frequencyStep, maxFreq, numberOfLateralBands, numberOfElevationalBands, detectorBoundsLateral, detectorBoundsElevational, options)
arguments
    detectedPhotons             struct
    frequencyStep               (1,1)   {mustBeNumeric, mustBePositive}                 %frequency step size
    maxFreq                     (1,1)   {mustBeNumeric, mustBePositive}                 %Maximum Doppler shift frequency taken into account
    numberOfLateralBands        (1,1)   {mustBeNumeric, mustBePositive}
    numberOfElevationalBands    (1,1)   {mustBeNumeric, mustBePositive}
    detectorBoundsLateral       (2,1)   {mustBeNumeric}
    detectorBoundsElevational   (2,1)   {mustBeNumeric}
    options.maxRadiusDetector   (1,1)   {mustBeNumeric, mustBePositive}     = 0.4e-3    %Maximum radial position on a circular detector
    options.minimumPhotonCount  (1,1)   {mustBeNumeric, mustBeNonnegative}  = 0
    options.saveFolder          string                                      = ''
end
bandCoordinatesLateral = linspace(detectorBoundsLateral(1),detectorBoundsLateral(2),numberOfLateralBands+1);
bandCoordinatesElevational = linspace(detectorBoundsElevations(1),detectorBoundsElevational(2),numberOfElevationalBands+1);

frequencyMax = round(maxFreq/frequencyStep).*frequencyStep;
freqAx = -frequencyMax/2:frequencyStep:frequencyMax/2;
numBinsDoppler = length(freqAx);

for l = 1:length(detectedPhotons)
    powerSpectra = zeros(numBinsDoppler,numberOfLateralBands,numberOfElevationalBands);
    batchPhotons = detectedPhotons(l);
    phase = mod(2*pi*batchPhotons.opticalPathLength./633e-9,2*pi);
    incidenceAngles = acosd(dot(batchPhotons.direction,[0;0;1].*ones(1,length(batchPhotons.DopplerShift))));
    transmittanceValues = glassPlateTransmittance(incidenceAngles);
    newAmplitudes = batchPhotons.amplitude.*transmittanceValues;
    batchFields = gpuArray(newAmplitudes.*batchPhotons.polarisationVector.*real(exp(-1i*phase)));
    batchPhotons.DopplerShift = gpuArray(batchPhotons.DopplerShift);
    batchDoppler = histcounts(batchPhotons.DopplerShift(sqrt(batchPhotons.position(1,:).^2 + batchPhotons.position(2,:).^2) < options.maxRadiusDetector),numBinsDoppler,'BinLimits',[-frequencyMax/2 - frequencyStep/2,frequencyMax/2 + frequencyStep/2]);
    lowerLimitFrequency = freqAx(find(batchDoppler >= options.minimumPhotonCount,1,'first'));
    upperLimitFrequency = freqAx(find(batchDoppler >= options.minimumPhotonCount,1,'last'));
    lowerLimitBin = find(freqAx == lowerLimitFrequency);
    upperLimitBin = find(freqAx == upperLimitFrequency);
    numBinsPowerSpec = upperLimitBin - lowerLimitBin;

    k = 1;
    for latBand = 1:numberOfLateralBands
        for elBand = 1:numberOfElevationalBands
            progressBar = waitbar(0, ['Run ' num2str(l) ', band ' num2str(latBand) '/' num2str(numberOfLateralBands.*numberOfElevationalBands)],'OuterPosition',[820 430 280 80]);            
            bandPhotons = gpuArray(batchPhotons.position(2,:) >= bandCoordinatesLateral(latBand) & batchPhotons.position(2,:) <= bandCoordinatesLateral(latBand+1) & batchPhotons.position(1,:) >= bandCoordinatesElevational(elBand) & batchPhotons.position(1,:) <= bandCoordinatesElevational(elBand+1) & isfinite(phase));
            for q = 0:numBinsPowerSpec - 1
                waitbar(q/numBinsPowerSpec,progressBar,['Run ' num2str(l) ', band ' num2str(k) '/' num2str(numberOfLateralBands.*numberOfElevationalBands)]);
                I = gpuArray(0);
                for n = lowerLimitBin:upperLimitBin - q
                    shiftCrit1 = gpuArray(batchPhotons.DopplerShift >= (-frequencyMax/2 - frequencyStep/2) + (n-1)*frequencyStep & batchPhotons.DopplerShift < (-frequencyMax/2 - frequencyStep/2) + (n)*frequencyStep & bandPhotons);
                    shiftCrit2 = gpuArray(batchPhotons.DopplerShift >= (-frequencyMax/2 - frequencyStep/2) + (n-1+q)*frequencyStep & batchPhotons.DopplerShift < (-frequencyMax/2 - frequencyStep/2) + (n+q)*frequencyStep & bandPhotons);
                    I = I + dot(sum([batchFields(shiftCrit1) batchFields(shiftCrit2)],2),sum([batchFields(shiftCrit1) batchFields(shiftCrit2)],2)) - dot(sum(batchFields(shiftCrit1),2),sum(batchFields(shiftCrit1),2)) - dot(sum(batchFields(shiftCrit2),2),sum(batchFields(shiftCrit2),2));
                end
                powerSpectra(q+1,latBand,elBand) = I;
            end
            close(progressBar)
            k = k+1;
            disp(['Band ' num2str(k-1) ', # of photons: ' num2str(sum(bandPhotons)) ' ' char(datetime)])
        end
    end
    powerSpectra = gather(reshape(powerSpectra,[numBinsDoppler,numberOfLateralBands*numberOfElevationalBands]));
    correctedPowerSpectra(l).freqAx = freqAx;
    correctedPowerSpectra(l).spectrum = mean(powerSpectra,2);
end
bruteForceOutcomesPol(simNumber).powerSpectra = correctedPowerSpectra;
save([options.saveFolder 'bruteForceOutcomesPol_' num2str(numberOfLateralBands) '_lat_' num2str(numberOfElevationalBands) '_el.mat'],'bruteForceOutcomesPol')
end