%%  bandedPowerSpectrum.m
%   calculates Doppler power spectra from McDoppler simulations using autocorrelation of Doppler spectrum. Divide up the detector in bands or pixels by changing the values of numberOfLateralBands and numberOfElevationalBands.
%% 
function crossCorrOutcomes = bandedPowerSpectrum(detectedPhotons, frequencyStep, maxFreq, numberOfLateralBands, numberOfElevationalBands, detectorBoundsLateral, detectorBoundsElevational, options)
arguments
    detectedPhotons             struct
    frequencyStep               (1,1)   {mustBeNumeric, mustBePositive}             %frequency step size
    maxFreq                     (1,1)   {mustBeNumeric, mustBePositive}             %Maximum Doppler shift frequency taken into account
    numberOfLateralBands        (1,1)   {mustBeNumeric, mustBePositive}
    numberOfElevationalBands    (1,1)   {mustBeNumeric, mustBePositive}
    detectorBoundsLateral       (2,1)   {mustBeNumeric}
    detectorBoundsElevational   (2,1)   {mustBeNumeric}
    options.maxRadiusDetector   (1,1)   {mustBeNumeric, mustBePositive} = 0.4e-3    %Maximum radial position on a circular detector
    options.saveFolder          string                                  = ''
end
    bandCoordinatesLateral = linspace(detectorBoundsLateral(1),detectorBoundsLateral(2),numberOfLateralBands+1);
    bandCoordinatesElevational = linspace(detectorBoundsElevations(1),detectorBoundsElevational(2),numberOfElevationalBands+1);
    frequencyMax = round(maxFreq/frequencyStep).*frequencyStep;
    freqAx = -frequencyMax/2:frequencyStep:frequencyMax/2;
    for l = 1:length(detectedPhotons)
        batchPhotons = detectedPhotons(l);
        photonsToInclude = batchPhotons.timesScattered <= 3 & sqrt(batchPhotons.position(1,:).^2 + batchPhotons.position(2,:).^2) <= options.maxRadiusDetector;
        batchPhotons.position = batchPhotons.position(:,photonsToInclude);
        batchPhotons.DopplerShift = batchPhotons.DopplerShift(photonsToInclude);
        powerSpectra = zeros(length(freqAx),numberOfLateralBands);
        k = 1;
        for n = 1:numberOfLateralBands
            for m = 1:numberOfElevationalBands
                bandPhotons = batchPhotons.position(2,:) >= bandCoordinatesLateral(n) & batchPhotons.position(2,:) <= bandCoordinatesLateral(n+1) & batchPhotons.position(1,:) >= bandCoordinatesElevational(m) & batchPhotons.position(1,:) <= bandCoordinatesElevational(m+1);
                DopplerDist = histcounts(batchPhotons.DopplerShift(bandPhotons),length(freqAx),'BinLimits',[-frequencyMax/2 - frequencyStep/2,frequencyMax/2 + frequencyStep/2]);
                powerSpectra(:,k) = fftshift(ifft(abs(fft(DopplerDist)).^2));
                k = k+1;
                disp(['Band ' num2str(k-1) ', # of photons: ' num2str(sum(bandPhotons))])

            end
        end
        correctedPowerSpectraOld(l).freqAx = freqAx;
        correctedPowerSpectraOld(l).spectrum = mean(powerSpectra,2);
        correctedPowerSpectraOld(l).DopplerSpectrum = DopplerDist;
    end
    crossCorrOutcomes(simNumber).powerSpectra = correctedPowerSpectraOld;
save([options.saveFolder 'crossCorrOutcomes_' num2str(particleSize) 'micron_' num2str(numberOfLateralBands) '_lat_' num2str(numberOfElevationalBands) '_el.mat'],'crossCorrOutcomes')
end