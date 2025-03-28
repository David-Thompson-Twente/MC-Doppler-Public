%%  testReflectRefract.m
%   David Thompson, last update 19-03-2025
%   Tests the functions related to reflection and refraction: reflectionCoeffs and updateDirectionsReflectRefract
%
%   Inputs:     struct fail_pass carrying any information about previously tested functions
%               classArgument string to determine class to be used, usually double, can also be string or gpuArray
%
%   Outputs:    PASS_ALL if all subtests succesfull
%               FAIL message indicating which part failed otherwise
%%
function fail_pass = testReflectRefract(fail_pass,classArgument)
arguments
    fail_pass       struct
    classArgument   string  {mustBeMember(classArgument,{'single','double','gpuArray'})}
end
testedFunctions = {'reflectionCoeffs','updateDirectionsReflectRefract'};
fail_pass.reflectionCoeffs = testReflectionCoeffs(classArgument);
fail_pass.updateDirectionsReflectRefract = testUpdateDirectionsReflectRefract(classArgument);
passCounter = 0;
fail_pass.reflectRefract = '';
for n = 1:length(testedFunctions)
    if strcmp(fail_pass.(testedFunctions{n}),'PASS')
        fail_pass = rmfield(fail_pass,testedFunctions{n});
        passCounter = passCounter + 1;
    else
        fail_pass.reflectRefract = [fail_pass.reflectRefract 9 testedFunctions{n} '_' fail_pass.(testedFunctions{n})];
        fail_pass = rmfield(fail_pass,testedFunctions{n});
    end
end
if passCounter == length(testedFunctions)
    fail_pass.reflectRefract = 'PASS_ALL';
end
end
%% Unit test for reflectionCoeffs (for planar interface under a few angles, check R_s and R_p for incidence 0-pi on both sides of interface. Check critical angle. Check Brewster angle)
function fail_pass = testReflectionCoeffs(classArgument,options)
arguments
    classArgument               string  {mustBeMember(classArgument,{'single','double','gpuArray'})}
    options.refractiveIndex     (2,1) {mustBeNumeric, mustBePositive}   = 1 + rand(2,1)
    options.numPhotons          (1,1) {mustBeInteger, mustBePositive}   = 1e6
    options.showPlots           (1,1) {mustBeNumericOrLogical}          = 0
    options.tolerance           (1,1) {mustBeNumeric, mustBePositive}   = eps(1)
    options.normalElevation     (1,1) {mustBeNumeric}                   = 2*pi*rand
end
[~,maxIndex] = max(options.refractiveIndex);
normalVectors = rotateRodrigues([0;0;-1],options.normalElevation,0,[],classArgument).*ones(1,options.numPhotons);
photons.direction = zeros(3,options.numPhotons);
media.refractiveIndex(1) = options.refractiveIndex(1);
media.refractiveIndex(2) = options.refractiveIndex(2);
tCrit = asin(min(options.refractiveIndex)/max(options.refractiveIndex));
tBrewster1 = atan(media.refractiveIndex(2)/media.refractiveIndex(1));
tBrewster2 = atan(media.refractiveIndex(1)/media.refractiveIndex(2));
photons.alive = ones(1,options.numPhotons);
crossingPhotons = logical(ones(1,options.numPhotons));
photons.mediumNum = 1 + round(rand(1,options.numPhotons));
photons.nextMedium(photons.mediumNum == 1) = 2;
photons.nextMedium(photons.mediumNum == 2) = 1;
numPhotonsMedium1 = length(normalVectors(photons.mediumNum == 1));
numPhotonsMedium2 = options.numPhotons - numPhotonsMedium1;
normalVectors(:,photons.mediumNum == 2) = -normalVectors(:,photons.mediumNum == 2);
anglesMedium(1).angles = [linspace(0,tBrewster1,floor(numPhotonsMedium1/2)) linspace(tBrewster1+eps(1),tCrit,floor(ceil(numPhotonsMedium1/2)/2)) linspace(tCrit+eps(1),pi/2,ceil(ceil(numPhotonsMedium1/2)/2))];
anglesMedium(2).angles = [linspace(0,tBrewster2,floor(numPhotonsMedium2/2)) linspace(tBrewster2+eps(1),tCrit,floor(ceil(numPhotonsMedium2/2)/2)) linspace(tCrit+eps(1),pi/2,ceil(ceil(numPhotonsMedium2/2)/2))];
photons.direction(:,photons.mediumNum == 1) = rotateRodrigues(-normalVectors(:,photons.mediumNum == 1),anglesMedium(1).angles,0,[],classArgument);
photons.direction(:,photons.mediumNum == 2) = rotateRodrigues(-normalVectors(:,photons.mediumNum == 2),anglesMedium(2).angles,0,[],classArgument);
reflectionCoefficients = reflectionCoeffs(photons,media,normalVectors,crossingPhotons,classArgument);
media(1).reflectionCoefficients = reflectionCoefficients(:,photons.mediumNum == 1);
media(2).reflectionCoefficients = reflectionCoefficients(:,photons.mediumNum == 2);
brewsterCheck = media(1).reflectionCoefficients(2,anglesMedium(1).angles == tBrewster1) < options.tolerance & media(2).reflectionCoefficients(2,anglesMedium(2).angles == tBrewster2) < options.tolerance;
critAngleCheck = media(maxIndex).reflectionCoefficients(:,anglesMedium(maxIndex).angles > tCrit) -1 < options.tolerance*100;
if options.showPlots
    figure;
    vecNum = randi(options.numPhotons);
    quiver3(0,0,0,normalVectors(1,vecNum),normalVectors(2,vecNum),normalVectors(3,vecNum));
    hold on
    quiver3(0,0,0,photons.direction(1,vecNum),photons.direction(2,vecNum),photons.direction(3,vecNum));
    daspect([1 1 1])
    figure;
    plot(anglesMedium(1).angles*180/pi,media(1).reflectionCoefficients(1,:),'.')
    hold on
    plot(anglesMedium(1).angles*180/pi,media(1).reflectionCoefficients(2,:),'.')
    plot(anglesMedium(2).angles*180/pi,media(2).reflectionCoefficients(1,:),'.')
    plot(anglesMedium(2).angles*180/pi,media(2).reflectionCoefficients(2,:),'.')
    title(['n_1 = ' num2str(media(1).refractiveIndex) ', n_2 = ' num2str(media(2).refractiveIndex), ', \theta_{crit} = ' num2str(tCrit*180/pi) '°, \theta_{B} = ' num2str(tBrewster1*180/pi) '/' num2str(tBrewster2*180/pi)])
end
fail_pass = 'FAIL';
if ~brewsterCheck
    fail_pass = [fail_pass '_BREWSTER_ANGLE'];
end
if min(critAngleCheck(:)) == 0
    fail_pass = [fail_pass '_CRITICAL_ANGLE'];
end
if length(fail_pass) > 4
    return
end
fail_pass = 'PASS';
end
%% Unit test for updateDirectionsReflectRefract
function fail_pass = testUpdateDirectionsReflectRefract(classArgument,options)
arguments
    classArgument               string  {mustBeMember(classArgument,{'single','double','gpuArray'})}
    options.reflectivity        (2,1) {mustBeNumeric, mustBePositive, mustBeLessThanOrEqual(options.reflectivity,1), mustBeNonnegative} = rand(2,1)
    options.refractiveIndex     (2,1) {mustBeNumeric, mustBePositive}                                                                   = 1 + rand(2,1)
    options.numPhotons          (1,1) {mustBeInteger, mustBePositive}                                                                   = 1e6
    options.tolerance           (1,1) {mustBeNumeric, mustBePositive}                                                                   = 1e-5
    options.normalElevation     (1,1) {mustBeNumeric}                                                                                   = 2*pi*rand
    options.showPlots           (1,1) {mustBeNumericOrLogical}                                                                          = 0
    options.biasReflect         (1,1) {mustBeNumeric}                                                                                   = 0
end
normalVectors = rotateRodrigues([0;0;-1],options.normalElevation,0,[],classArgument).*ones(1,options.numPhotons);
tCrit = asin(min(options.refractiveIndex)/max(options.refractiveIndex));
media(1).refractiveIndex = options.refractiveIndex(1);
media(2).refractiveIndex = options.refractiveIndex(2);
photons.timesReflected = zeros(1,options.numPhotons);
photons.timesRefracted = zeros(1,options.numPhotons);
photons.mediumNum = 1 + round(rand(1,options.numPhotons));
photons.nextMedium(photons.mediumNum == 1) = 2;
photons.nextMedium(photons.mediumNum == 2) = 1;
photons.likelihood = ones(1,options.numPhotons);
photons.alive = ones(1,options.numPhotons);
n1 = horzcat(media(photons.mediumNum).refractiveIndex);
n2 = horzcat(media(photons.nextMedium).refractiveIndex);
indexRatio = n2./n1;
crossingPhotons = logical(ones(1,options.numPhotons));
reflectionCoefficients = options.reflectivity.*ones(1,options.numPhotons);
numPhotonsMedium1 = length(normalVectors(photons.mediumNum == 1));
numPhotonsMedium2 = options.numPhotons - numPhotonsMedium1;
normalVectors(:,photons.mediumNum == 2) = -normalVectors(:,photons.mediumNum == 2);
anglesMedium(1).angles = linspace(0,tCrit,numPhotonsMedium1);
anglesMedium(2).angles = linspace(0,tCrit,numPhotonsMedium2);
photons.direction(:,photons.mediumNum == 1) = rotateRodrigues(-normalVectors(:,photons.mediumNum == 1),anglesMedium(1).angles,0,[],classArgument);
photons.direction(:,photons.mediumNum == 2) = rotateRodrigues(-normalVectors(:,photons.mediumNum == 2),anglesMedium(2).angles,0,[],classArgument);

polCross = cross(photons.direction, normalVectors);
polCross(:,sum(polCross) == 0) = rotateRodrigues(photons.direction(:,sum(polCross) == 0),pi/2,0,[],classArgument);
polarisationStates(1).polarisationVecs = polCross./vecnorm(polCross);
polarisationStates(2).polarisationVecs = rotateRodrigues(photons.direction,pi/2,pi/2,polarisationStates(1).polarisationVecs,classArgument);
polarisationStates(3).polarisationVecs = rotateRodrigues(photons.direction,pi/2,2*pi*rand,polarisationStates(1).polarisationVecs,classArgument);

for n = 1:3
    photons.polarisationVector = polarisationStates(n).polarisationVecs;
    processedPhotons = updateDirectionsReflectRefract(photons,normalVectors,crossingPhotons,reflectionCoefficients,indexRatio,'polarised',classArgument);
    reflectedPhotons = processedPhotons.nextMedium == photons.mediumNum;
    refractedPhotons = ~reflectedPhotons;
    reflectionCheck = dot(photons.direction(:,reflectedPhotons),normalVectors(:,reflectedPhotons)) + dot(processedPhotons.direction(:,reflectedPhotons),normalVectors(:,reflectedPhotons)) < options.tolerance;
    refractionCheck = n1(refractedPhotons).*sin(acos(dot(photons.direction(:,refractedPhotons),normalVectors(:,refractedPhotons)))) - n2(refractedPhotons).*sin(acos(dot(processedPhotons.direction(:,refractedPhotons),-normalVectors(:,refractedPhotons)))) < options.tolerance;

    if n == 1
        polarisationCheck = abs(dot(photons.polarisationVector,processedPhotons.polarisationVector)) - 1 < options.tolerance;
        reflectivityCheck = length(reflectedPhotons(reflectedPhotons == 1))/options.numPhotons - options.reflectivity(1) < 1e-2;
    elseif n == 2
        polarisationCheck = abs(dot(photons.polarisationVector,processedPhotons.polarisationVector)) - abs(dot(photons.direction,processedPhotons.direction)) < options.tolerance;
        reflectivityCheck = length(reflectedPhotons(reflectedPhotons == 1))/options.numPhotons - options.reflectivity(2) < 1e-2;
    elseif n == 3
        polarisationCheck = abs(dot(photons.polarisationVector,photons.direction)) < options.tolerance;
        expectedR = (options.reflectivity(1).*abs(dot(polarisationStates(1).polarisationVecs,photons.polarisationVector)) + options.reflectivity(2)*abs(dot(polarisationStates(2).polarisationVecs,photons.polarisationVector)))./(abs(dot(polarisationStates(1).polarisationVecs,photons.polarisationVector)) + abs(dot(polarisationStates(2).polarisationVecs,photons.polarisationVector)));
        reflectivityCheck = length(reflectedPhotons(reflectedPhotons == 1))/options.numPhotons - expectedR < 1e-2;
    end
    fail_pass = 'FAIL';
    if min(reflectivityCheck) == 0
        fail_pass = [fail_pass '_REFLECTIVITY'];
    end
    if min(reflectionCheck) == 0
        fail_pass = [fail_pass '_REFLECTION'];
    end
    if min(refractionCheck) == 0
        fail_pass = [fail_pass '_REFRACTION'];
    end
    if min(polarisationCheck) == 0
        fail_pass = [fail_pass '_POLARISATION'];
    end

    if options.showPlots
        figure
        plot(dot(photons.direction(:,reflectedPhotons),normalVectors(:,reflectedPhotons)))
        hold on
        plot(dot(processedPhotons.direction(:,reflectedPhotons),normalVectors(:,reflectedPhotons)),'--')
        plot(dot(photons.direction(:,reflectedPhotons),normalVectors(:,reflectedPhotons)) + dot(processedPhotons.direction(:,reflectedPhotons),normalVectors(:,reflectedPhotons)))
        figure
        plot(n1(refractedPhotons).*sin(acos(dot(photons.direction(:,refractedPhotons),normalVectors(:,refractedPhotons)))))
        hold on
        plot(n2(refractedPhotons).*sin(acos(dot(processedPhotons.direction(:,refractedPhotons),-normalVectors(:,refractedPhotons)))),'--')
        plot(n1(refractedPhotons).*sin(acos(dot(photons.direction(:,refractedPhotons),normalVectors(:,refractedPhotons)))) - n2(refractedPhotons).*sin(acos(dot(processedPhotons.direction(:,refractedPhotons),-normalVectors(:,refractedPhotons)))))
        title(['n_1 = ' num2str(options.refractiveIndex(1)), ', n_2 = ' num2str(options.refractiveIndex(2))])
    end
end
if length(fail_pass) > 4
    return
end
fail_pass = 'PASS';
end