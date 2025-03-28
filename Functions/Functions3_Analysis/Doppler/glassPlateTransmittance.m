%%  glassPlateTransmittance.m
%   Calculates angle-dependent transmittance through a finite thickness glass plate based on http://www.jedsoft.org/physics/notes/multilayer.pdf
%   Outputs an array of transmittance values for each photon.
%%
function transmittance = glassPlateTransmittance(incidenceAngles,options)
arguments
    incidenceAngles             (1,:)   {mustBeNumeric}                 %incidence angles in degrees relative to normal
    options.mu_vac              (1,1)   {mustBeNumeric,mustBePositive}  =   4*pi*1e-7
    options.e_vac               (1,1)   {mustBeNumeric,mustBePositive}  =   8.85e-12
    options.mu_air              (1,1)   {mustBeNumeric,mustBePositive}  =   1.26e-6
    options.n_air               (1,1)   {mustBeNumeric,mustBePositive}  =   1
    options.n_glass             (1,1)   {mustBeNumeric,mustBePositive}  =   1.457
    options.e_glass             (1,1)   {mustBeNumeric,mustBePositive}  =   3.75
    options.plateThickness      (1,1)   {mustBeNumeric,mustBePositive}  =   500e-6
    options.wavelength          (1,1)   {mustBeNumeric,mustBePositive}  =   633e-9
    options.plotTransmittance   (1,1)   {mustBeNumericOrLogical}        =   0   
end
mu_glass = options.mu_vac*(options.n_glass^2/options.e_glass);
layerThickness = [0 options.plateThickness 0];

kz_1 = 2*pi/options.wavelength.*sqrt(options.n_glass^2 - cosd(90-incidenceAngles).^2);
kz_0 = 2*pi/options.wavelength.*sqrt(options.n_air^2 - cosd(90-incidenceAngles).^2);
k = [kz_0;kz_1;kz_0];
mu = [options.mu_air mu_glass options.mu_air];

R = zeros(3,length(incidenceAngles));
f = k./mu';
a = exp(2.*1i.*k.*layerThickness');
F = (f - circshift(f,[-1,0]))./(f + circshift(f,[-1,0]));

for n = 1:2
    R(end-n,:) = a(end-n,:).*(F(end-n,:) + R(end-n,:))./(1 + F(end-n,:).*R(end-(n-1),:));
end
T = exp(1i.*circshift(k,[-1,0]).*circshift(layerThickness',-1)).*(1+F)./(1+F.*circshift(R,[-1,0]));
transmittance = abs(prod(T(1:end-1,:))).^2;
if options.plotTransmittance
    figure; plot(incidenceAngles,transmittance);
end
end