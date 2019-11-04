function [ output_args ] = Proposed_path_loss_model(distance2, moisture, sand, clay)
%we consider three distances: 
%distance1 for the first aboveground area
%distance2 for underground place
%distance3 for  the last aboveground area

% the fixed operating frequency is 433MHz
frequency = 433e+6;
frequency_dielectric = frequency;
temperature = 40;
distance1 = 0.13;
distance3 = 0.13;

%distance2 = distance/cosd(30) ;
%--------------------------------------------------------------------------
%                  Peplinski's complex dielectric prediction
%--------------------------------------------------------------------------
bulk_density = 1.5 ; % bulk density take from the modified friis model in g/cm3
epsilon_infinite = 4.9;
epsilon_vaccum = 8.854187e-12;
epsilon_water = 80.1;
relaxation_time = 0.58e-10; %2*pi*relation_time
particle_density = 2.66; % specific density of the solid soil particles in g/cm3 
beta_prime = 1.2748 - 0.519 * sand - 0.152 * clay ;%the value of beta prime which depend of the sand and clay portions
beta_prime_prime = 1.33797 - 0.603*sand - 0.166*clay ;%the value of beta prime prime which depend of the sand and clay portions
relative_CDC = ((1.01 + 0.44*particle_density)^2) - 0.062; % The relative complex dielectric constant
sigma = 0.65 ;% A predefined constant

effective_conductivity = 0.0467 + (0.2204 * bulk_density) - (0.4111 * sand) + (0.6614 * clay); %Value of the effective conductive in terms of the frequency and the sand, clay portion

% Real and imaginary parts of the CDC of free water
epsilon_real_water = epsilon_infinite + (epsilon_water - epsilon_infinite)/(1 + (relaxation_time*frequency_dielectric)^2) ;
%epsilon_imaginary_water = ((relaxation_time*frequency_dielectric*(epsilon_vaccuum - epsilon_infinite))/(1 + (relaxation_time*frequency_dielectric)^2)) + ((effective_conductivity)*(particle_density - bulk_density))/(2*pi*epsilon_vaccuum*particle_density*moisture);
epsilon_imaginary_water = ((relaxation_time * frequency_dielectric * (epsilon_water - epsilon_infinite))/(1 + (relaxation_time * frequency_dielectric)^2)) + ((effective_conductivity) * (particle_density - bulk_density)) / (2 * pi * frequency_dielectric * epsilon_vaccum * particle_density * moisture);
    
%Real and imaginary parts of the complex mixing dielectric constant of soil-water
epsilon_real_CDC_peplinski = 1.15*(1 + ((bulk_density*((relative_CDC)^(sigma) - 1))/particle_density) + moisture^(beta_prime)*epsilon_real_water^(sigma) - moisture)^(1/sigma) - 0.68
epsilon_imaginary_CDC_peplinski = (((moisture)^(beta_prime_prime)) * ((epsilon_imaginary_water)^(sigma)))^(1/sigma)

%--------------------------------------------------------------------------
%                   End of the pleplinski's DC prediction
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%               Start of the MBSDM dielectric prediction
%--------------------------------------------------------------------------
RI_dry = 1.634 - (0.539e-2)*clay + (0.2748e-4)*clay^2; %RI_dry is the refractive index of dry soil
NAC_dry = 0.03952 - (0.04038e-2)*clay; %NAC_dry is the Normalized Attenuation Coefficient known as k for dry soil 
MBWF = 0.02863 + (0.30673e-2)*clay; %MBWF is the Maximum Bound Water Fraction

epsilon_starting_bound = 79.8 - (85.4e-2)*clay + (32.7e-4)*clay^2 ;
relaxation_bound = 1.062e-11 + (3.450e-12)*10^(-2)*clay; % relaxation time of bound water
sigma_bound = 0.3112 + (0.467e-2)*clay; %Conductivity of bound water

epsilon_starting_free = 100;
relaxation_free = 8.5e-12 ; % relaxation time of free water
sigma_free = 0.3631 + (1.217e-2)*clay; %Conductivity of free water

CDC_real_bound = epsilon_infinite + (epsilon_starting_bound - epsilon_infinite)/(1+(2*pi*frequency_dielectric*relaxation_bound)^2);%real part of the CDC of bound water 
CDC_imaginary_bound = (((epsilon_starting_bound - epsilon_infinite)*(2*pi*frequency_dielectric*relaxation_bound))/(1+(2*pi*frequency_dielectric*relaxation_bound)^2)) + ((sigma_bound)/(2*pi*frequency_dielectric*epsilon_vaccum));%imaginary part of the CDC of bound water 

CDC_real_free = epsilon_infinite + (epsilon_starting_free - epsilon_infinite)/(1+(2*pi*frequency_dielectric*relaxation_free)^2);%real part of the CDC of free water 
CDC_imaginary_free = (((epsilon_starting_free - epsilon_infinite)*(2*pi*frequency_dielectric*relaxation_free))/(1+(2*pi*frequency_dielectric*relaxation_free)^2)) + ((sigma_free)/(2*pi*frequency_dielectric*epsilon_vaccum));%imaginary part of the CDC of free water 

RI_bound = sqrt(sqrt((CDC_real_bound)^2 + (CDC_imaginary_bound)^2) + CDC_real_bound)/sqrt(2); %RI_bound is the refractive index of bound water
NAC_bound = sqrt(sqrt((CDC_real_bound)^2 + (CDC_imaginary_bound)^2) - CDC_real_bound)/sqrt(2); %NAC_bound is the Normalized Attenuation Coefficient known as k for bound water 

RI_free = sqrt(sqrt((CDC_real_free)^2 + (CDC_imaginary_free)^2) + CDC_real_free)/sqrt(2); %RI_free is the refractive index of free water
NAC_free = sqrt(sqrt((CDC_real_free)^2 + (CDC_imaginary_free)^2) - CDC_real_free)/sqrt(2); %NAC_free is the Normalized Attenuztion Coefficient known as k for free water 

if moisture<MBWF
    RI_mbsdm = RI_dry +(RI_bound - 1)*moisture ; %resulting Refractive index
    NAC_mbsdm = NAC_dry + NAC_bound*moisture ; %resulting Normalized attenuation 
else
    RI_mbsdm = RI_dry +(RI_bound - 1)*MBWF + (RI_free - 1)*(moisture - MBWF); %resulting Refractive index
    NAC_mbsdm = NAC_dry + NAC_bound*MBWF + NAC_free*(moisture - MBWF); %resulting Normalized attenuation
end

epsilon_real_CDC_mbsdm = RI_mbsdm^2 - NAC_mbsdm^2  %Real part of the mixing CDC
epsilon_imaginary_CDC_mbsdm = 2*RI_mbsdm*NAC_mbsdm %Imaginary part of the mixing CDC

%--------------------------------------------------------------------------
%                   End of the MBSDM dielectric prediction
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%                   Start of the TMDM dielectric prediction
%--------------------------------------------------------------------------
temperature_starting = 20; %Starting temperature is set to 10¤ Celsius
beta_starting_bound = 8.67e-19 - 0.00126e-2*clay + 0.00184e-4*clay^2 - 9.77e-10*clay^3 - 1.39e-15*clay^4; % Volumetric expansion coefficient of bound water
beta_starting_free = 1.11e-4 - 1.603e-7*clay + 1.239e-9*clay^2 + 8.33e-13*clay^3 -1.007e-14*clay^4; % Volumetric expansion coefficient of free water

activation_energy_bound = 1467 + 2697e-2*clay - 980e-4*clay^2 + 1.368e-10*clay^3 - 8.61e-13*clay^4; %Activation energy of bound water over the universal gas constant
entropy_activation_bound = 0.888 + 9.7-2*clay - 4.262e-4*clay^2 + 6.79e-21*clay^3 +4.263e-22*clay^4; %Zntropy of activation of bound water over the universal gas constant

activation_energy_free = 2231 - 143.1*10^(-2)*clay + 223.2*10^(-4)*clay^2 - 142.1*10^(-6)*clay^3 + 27.14*10^(-8)*clay^4; %Activation energy of free water over the universal gas constant
entropy_activation_free = 3.649 - 0.4894*10^(-2)*clay + 0.763*10^(-4)*clay^2 - 0.4859*10^(-6)*clay^3 + 0.0928*10^(-8)*clay^4; %Zntropy of activation of free water over the universal gas constant

sigma_bound_starting = 0.3112 + 0.467*10^(-2)*clay; %Conductivity of bound water
beta_temperature_bound = 0.0028 + 0.02094*10^(-2)*clay - 0.01229*10^(-4)*clay^2 - 5.03*10^(-22)*clay^3 + 4.163*10^(-24)*clay^4; %Temperature incrementation coefficient for conductivity of bound water

sigma_free_starting = 0.05 + 14*(1-(1-clay*10^(-2))^(4.664)); %Conductivity of free water
beta_temperature_free = 0.00108 + 0.1413e-2*clay - 0.2555e-4*clay^2 + 0.2147e-6*clay^3 - 0.0711e-8*clay^4; %Temperature incrementation coefficient for conductivity of free water

sigma_bound_tmdm = sigma_bound_starting + beta_temperature_bound * (temperature - temperature_starting); %Conductivity of bound water at a temperature t
sigma_free_tmdm = sigma_free_starting + beta_temperature_free * (temperature - temperature_starting); %Conductivity of free water at a temperature t 
x1 = temperature + 273.15;

relaxation_time_bound = (48*10^(-12))/(x1)*exp(activation_energy_bound/x1 - entropy_activation_bound)*10^(-12); %relaxation time of bound water
relaxation_time_free = (48*10^(-12))/(x1)*exp(activation_energy_free/x1 - entropy_activation_free)*10^(-12); %relaxation time of free water

f_starting_bound = log((epsilon_starting_bound - 1)/(epsilon_starting_bound +2)); %Computation of the function F in terms of the starting temperature of bound water
f_starting_free = log((epsilon_starting_free - 1)/(epsilon_starting_free +2)); %Computation of the function F in terms of the starting temperature of free water

%Calculation of the epsilon zero of bound water at the temperature t
epsilon_zero_bound = (1 + 2*exp(f_starting_bound - beta_starting_bound*(temperature - temperature_starting)))/(1 - exp(f_starting_bound -  beta_starting_bound*(temperature - temperature_starting))); 

%Calculation of the epsilon zero of free water at the temperature t
epsilon_zero_free = (1 + 2*exp(f_starting_free - beta_starting_free*(temperature - temperature_starting)))/(1 - exp(f_starting_free -  beta_starting_free*(temperature - temperature_starting))); 

angular_velocity = 2*pi*frequency_dielectric*epsilon_vaccum;

%Computation of the CDC of bound water
CDC_real_bound_tmdm = epsilon_infinite + ((epsilon_zero_bound - epsilon_infinite)/(1 + (2*pi*frequency_dielectric*relaxation_time_bound)^2));
CDC_imaginary_bound_tmdm = ((2*pi*frequency_dielectric*relaxation_time_bound*(epsilon_zero_bound - epsilon_infinite))/(1 + (2*pi*frequency_dielectric*relaxation_time_bound)^2))+((2*pi*frequency_dielectric*sigma_bound_tmdm*relaxation_time_bound*relaxation_time_bound)/(epsilon_vaccum*(1 + (2*pi*frequency_dielectric*relaxation_time_bound)^2)))+((sigma_bound_tmdm)/(angular_velocity*(1 + (2*pi*frequency_dielectric*relaxation_time_bound)^2))) ;

%Computation of the CDC of free water
CDC_real_free_tmdm = epsilon_infinite + ((epsilon_zero_free - epsilon_infinite)/(1 + (2*pi*frequency_dielectric*relaxation_time_free)^2));
CDC_imaginary_free_tmdm = ((2*pi*frequency_dielectric*relaxation_time_free*(epsilon_zero_free - epsilon_infinite))/(1 + (2*pi*frequency_dielectric*relaxation_time_free)^2))+((2*pi*frequency_dielectric*sigma_free_tmdm*relaxation_time_free*relaxation_time_free)/(epsilon_vaccum*(1 + (2*pi*frequency_dielectric*relaxation_time_free)^2)))+((sigma_free_tmdm)/(angular_velocity*(1 + (2*pi*frequency_dielectric*relaxation_time_free)^2))) ;

%RI and NAC of bound and free water 
RI_bound_tmdm = sqrt(sqrt((CDC_real_bound_tmdm)^2 + (CDC_imaginary_bound_tmdm)^2) + CDC_real_bound_tmdm)/sqrt(2); %RI_bound is the refractive index of bound water
NAC_bound_tmdm = sqrt(sqrt((CDC_real_bound_tmdm)^2 + (CDC_imaginary_bound_tmdm)^2) - CDC_real_bound_tmdm)/sqrt(2); %NAC_bound is the Normalized Attenuation Coefficient known as k for bound water 

RI_free = sqrt(sqrt((CDC_real_free_tmdm)^2 + (CDC_imaginary_free_tmdm)^2) + CDC_real_free_tmdm)/sqrt(2); %RI_free is the refractive index of free water
NAC_free = sqrt(sqrt((CDC_real_free_tmdm)^2 + (CDC_imaginary_free_tmdm)^2) - CDC_real_free_tmdm)/sqrt(2); %NAC_free is the Normalized Attenuztion Coefficient known as k for free water 

if moisture<MBWF
    RI_tmdm = RI_dry + (RI_bound_tmdm - 1)*moisture;
    NAC_tmdm = NAC_dry + NAC_bound_tmdm*moisture
else
    RI_tmdm = RI_dry + (RI_bound_tmdm - 1)*MBWF + (RI_free - 1)*(moisture - MBWF);
    NAC_tmdm = NAC_dry + NAC_bound_tmdm*MBWF + NAC_free*(moisture - MBWF);
end

epsilon_real_CDC_tmdm = RI_tmdm^2 - NAC_tmdm^2 
epsilon_imaginary_CDC_tmdm = 2*RI_tmdm*NAC_tmdm 
    
%--------------------------------------------------------------------------
%                   End of the TMDM dielectric prediction
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%                Path Loss prediction with the Path losses
%--------------------------------------------------------------------------
mu_factor = pi * 4 * 10^(-7) * 8.854 * 10^(-12);
omega = 2*pi*frequency;    
%--------------------------------------------------------------------------
%                        Peplinski permittivity prediction
%--------------------------------------------------------------------------
permea_peplinski = (mu_factor*epsilon_real_CDC_peplinski)/2 ;
reflect_peplinski = ((1-sqrt(epsilon_real_CDC_peplinski))/(1+sqrt(epsilon_real_CDC_peplinski)))^2;
reflect_attenuation_peplinski = 10*log((2*reflect_peplinski)/(1+reflect_peplinski));

loss_expo_peplinski = sqrt(1+((epsilon_imaginary_CDC_peplinski)/(epsilon_real_CDC_peplinski))^2);

alph_peplinski = omega * (sqrt(permea_peplinski * (loss_expo_peplinski - 1)));
beta_peplinski = omega * (sqrt(permea_peplinski * (loss_expo_peplinski + 1)));

Conventional_modified_friis = 6.4 + (20*log10(distance2*beta_peplinski)) + (8.68 * alph_peplinski * distance2) 
NC_modified_friis = Conventional_modified_friis + reflect_attenuation_peplinski
%--------------------------------------------------------------------------
%                        MBSDM permittivity prediction
%--------------------------------------------------------------------------
permea_mbsdm = (mu_factor*epsilon_real_CDC_mbsdm)/2 ;
reflect_mbsdm = ((1-sqrt(epsilon_real_CDC_mbsdm))/(1+sqrt(epsilon_real_CDC_mbsdm)))^2;
reflect_attenuation_mbsdm = 10*log((2*reflect_mbsdm)/(1+reflect_mbsdm));

loss_expo_mbsdm = sqrt(1+((epsilon_imaginary_CDC_mbsdm)/(epsilon_real_CDC_mbsdm))^2);

alph_mbsdm = omega * (sqrt(permea_mbsdm * (loss_expo_mbsdm - 1)));
beta_mbsdm = omega * (sqrt(permea_mbsdm * (loss_expo_mbsdm + 1)));
refractive = sqrt((sqrt(epsilon_real_CDC_mbsdm^2 + epsilon_imaginary_CDC_mbsdm^2) + epsilon_real_CDC_mbsdm) / 2);
refractive_attenuation = 20 * log10((refractive + 1)/4);

MBSDM_modified_friis = 6.4 + (20*log10(distance2*beta_mbsdm)) + (8.68 * alph_mbsdm * distance2)
Proposed_LPL1 = reflect_attenuation_mbsdm + 20*log10(distance1 * distance3 * distance2 * beta_mbsdm) - 288.8 + 40*log10(frequency_dielectric) + (8.68 * alph_mbsdm * distance2)
Proposed_LPL2 = 20*log10(distance1 * distance3 * distance2 * beta_mbsdm) - 288.8 + 40*log10(frequency_dielectric) + (8.68 * alph_mbsdm * distance2)

%--------------------------------------------------------------------------
%                        TMDM permittivity prediction
%--------------------------------------------------------------------------
permea_tmdm = (mu_factor*epsilon_real_CDC_tmdm)/2 ;
reflect_tmdm = ((1-sqrt(epsilon_real_CDC_tmdm))/(1+sqrt(epsilon_real_CDC_tmdm)))^2;
reflect_attenuation_tmdm = 10*log((2*reflect_tmdm)/(1+reflect_tmdm));
refractive_tmdm = sqrt((sqrt((epsilon_real_CDC_tmdm)^2 + (epsilon_imaginary_CDC_tmdm)^2) + epsilon_real_CDC_tmdm) / 2);
refractive_attenuation_tmdm = 20 * log10((refractive_tmdm + 1)/4);

loss_expo_tmdm = sqrt(1+((epsilon_imaginary_CDC_tmdm)/(epsilon_real_CDC_tmdm))^2);

alph_tmdm = omega * (sqrt(permea_tmdm * (loss_expo_tmdm - 1)));
beta_tmdm = omega * (sqrt(permea_tmdm * (loss_expo_tmdm + 1)));

Proposed_LPL1_tmdm = refractive_attenuation_tmdm + reflect_attenuation_tmdm + 20*log10(distance1 * distance3 * distance2 * beta_tmdm) - 288.8 + 40*log10(frequency_dielectric) + (8.68 * alph_tmdm * distance2);
Proposed_LPL2_tmdm = refractive_attenuation_tmdm + 20*log10(distance1 * distance3 * distance2 * beta_tmdm) - 288.8 + 40*log10(frequency_dielectric) + (8.68 * alph_tmdm * distance2);

end

