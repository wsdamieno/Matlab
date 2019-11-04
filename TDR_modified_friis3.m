function [ ] = TDR_modified_friis3( moisture, sand, clay )

% the fixed operating frequency is 433MHz
frequency_dielectric = 433e+6;
temperature = 20 ;


%Distance varies from 0 to 20m
distance = 0.010:0.0001:6;

%--------------------------------------------------------------------------
%                  Peplinski's complex dielectric prediction
%--------------------------------------------------------------------------
bulk_density = 1.5 ; % bulk density take from the modified friis model in g/cm3
epsilon_infinite = 4.9;
epsilon_vaccum = 8.854187e-12;
epsilon_water = 80.1;
relaxation_time = 0.58e-10; %2*pi*relation_time
particle_density = 2.66; % specific density of the solid soil particles in g/cm3 
beta_prime = 1.2748 - 0.519*sand - 0.152*clay ;%the value of beta prime which depend of the sand and clay portions
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
NAC_dry = 0.03952 - (0.04038e-2)*clay; %NAC_dry is the Normalized Attenuztion Coefficient known as k for dry soil 
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
omega = 2*pi*frequency_dielectric;    
distance_boite1 = 0.12;
distance_boite2 = 0.12;
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%                            MODIFIED FRIIS COMPUTATIONS
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
%                               TDR prediction
%--------------------------------------------------------------------------
% permea_tdr = (mu_factor*epsilon_real_CDC_tdr)/2 ;
% reflect_tdr = ((1-sqrt(epsilon_real_CDC_tdr))/(1+sqrt(epsilon_real_CDC_tdr)))^2;
% reflect_attenuation_tdr = 10*log((2*reflect_tdr)/(1+reflect_tdr));
% 
% loss_expo_tdr = sqrt(1+((epsilon_imaginary_CDC_tdr)/(epsilon_real_CDC_tdr))^2);
% 
% alph_tdr = omega * (sqrt(permea_tdr * (loss_expo_tdr - 1)));
% beta_tdr = omega * (sqrt(permea_tdr * (loss_expo_tdr + 1)));
% 
% Loss_Path_modified_friis_tdr = 6.4 + (20*log10(distance*beta_tdr)) + (8.68 * alph_tdr * distance) ;
% Loss_Path_NC_tdr = Loss_Path_modified_friis_tdr + reflect_attenuation_tdr;
%--------------------------------------------------------------------------
%                        Peplinski permittivity prediction
%--------------------------------------------------------------------------
permea_peplinski = (mu_factor*epsilon_real_CDC_peplinski)/2 ;
reflect_peplinski = ((1-sqrt(epsilon_real_CDC_peplinski))/(1+sqrt(epsilon_real_CDC_peplinski)))^2;
reflect_attenuation_peplinski = 10*log((2*reflect_peplinski)/(1+reflect_peplinski));

loss_expo_peplinski = sqrt(1+((epsilon_imaginary_CDC_peplinski)/(epsilon_real_CDC_peplinski))^2);

alph_peplinski = omega * (sqrt(permea_peplinski * (loss_expo_peplinski - 1)))
beta_peplinski = omega * (sqrt(permea_peplinski * (loss_expo_peplinski + 1)))

distance_15 = 0.40;
Loss_Path_modified_friis_peplinski = 6.4 + (20*log10(distance*beta_peplinski)) + (8.68 * alph_peplinski * distance) ;

Loss_Path_modified_friis_peplinski_15 = 6.4 + (20*log10(distance_15*beta_peplinski)) + (8.68 * alph_peplinski * distance_15) ;

Loss_Path_NC_peplinski = Loss_Path_modified_friis_peplinski + reflect_attenuation_peplinski;
% CRIM-Fresnel Path loss
electrical_conduct = 2.32;
epsilon_vacuum = 8.854 * 10^(-12);

alph_num = 8.68 * 60 * pi * ((2 * pi * frequency_dielectric * epsilon_vacuum * epsilon_real_CDC_peplinski) + (electrical_conduct) );
alph_denum = sqrt(( epsilon_real_CDC_peplinski / 2 )*( 1 + (sqrt( 1 + (( epsilon_imaginary_CDC_peplinski + (electrical_conduct) / (2 * pi * frequency_dielectric * epsilon_vacuum))/(epsilon_real_CDC_peplinski))^2))));
alph_CRIM = alph_num / alph_denum;

Loss_Path_CRIM_peplinski = alph_CRIM * distance + reflect_attenuation_peplinski;
%--------------------------------------------------------------------------
Loss_Path_proposed_peplinski = reflect_attenuation_peplinski + 20*log10(distance_boite1 * distance_boite2) - 288.8 + 40*log10(frequency_dielectric) + 20*log10(distance) + 20*log10(beta_peplinski) + (8.68 * alph_peplinski * distance);
%--------------------------------------------------------------------------
%                        MBSDM permittivity prediction
%--------------------------------------------------------------------------
permea_mbsdm = (mu_factor*epsilon_real_CDC_tmdm)/2 ;
reflect_mbsdm = ((1-sqrt(epsilon_real_CDC_tmdm))/(1+sqrt(epsilon_real_CDC_tmdm)))^2;
reflect_attenuation_mbsdm = 10*log((2*reflect_mbsdm)/(1+reflect_mbsdm));

loss_expo_mbsdm = sqrt(1+((epsilon_imaginary_CDC_tmdm)/(epsilon_real_CDC_tmdm))^2);

refractive = sqrt((sqrt(epsilon_real_CDC_tmdm^2 + epsilon_imaginary_CDC_tmdm^2) + epsilon_real_CDC_tmdm) / 2)
refractive_attenuation = 20 * log10((refractive + 1)/4)

alph_mbsdm = omega * (sqrt(permea_mbsdm * (loss_expo_mbsdm - 1)))
beta_mbsdm = omega * (sqrt(permea_mbsdm * (loss_expo_mbsdm + 1)))

Loss_Path_modified_friis_mbsdm = 6.4 + (20*log10(distance*beta_mbsdm)) + (8.68 * alph_mbsdm * distance) ;
Loss_Path_NC_mbsdm = reflect_attenuation_mbsdm + 6.4 + (20*log10(distance*beta_mbsdm)) + (8.68 * alph_mbsdm * distance) ;
Loss_Path_proposed_mbsdm = 20*log10(distance_boite1 * distance_boite2 * distance * beta_mbsdm) - 288.8 + 40*log10(frequency_dielectric) + (8.68 * alph_mbsdm * distance)
Loss_Path_proposed_mbsdm_reflect = reflect_attenuation_mbsdm + 20*log10(distance_boite1 * distance_boite2 * distance * beta_mbsdm) - 288.8 + 40*log10(frequency_dielectric) + (8.68 * alph_mbsdm * distance)
% Loss_Path_proposed_mbsdm_reflect_inv = 20*log10(distance_boite1 * distance_boite2 * distance * beta_mbsdm) - refractive_attenuation - 288.8 + 40*log10(frequency_dielectric) + (8.68 * alph_mbsdm * distance)
% Loss_Path_proposed_mbsdm_refract =  refractive_attenuation + 20*log10(distance_boite1 * distance_boite2 * distance * beta_mbsdm) - 288.8 + 40*log10(frequency_dielectric) + (8.68 * alph_mbsdm * distance)
% Loss_Path_proposed_mbsdm_reflect_refract =  refractive_attenuation + reflect_attenuation_mbsdm + 20*log10(distance_boite1 * distance_boite2 * distance * beta_mbsdm) - 288.8 + 40*log10(frequency_dielectric) + (8.68 * alph_mbsdm * distance)

%--------------------------------------------------------------------------
%                        TMDM permittivity prediction
%--------------------------------------------------------------------------
permea_tmdm = (mu_factor*epsilon_real_CDC_tmdm)/2 ;
reflect_tmdm = ((1-sqrt(epsilon_real_CDC_tmdm))/(1+sqrt(epsilon_real_CDC_tmdm)))^2;
reflect_attenuation_tmdm = 10*log((2*reflect_tmdm)/(1+reflect_tmdm));

loss_expo_tmdm = sqrt(1+((epsilon_imaginary_CDC_tmdm)/(epsilon_real_CDC_tmdm))^2);

alph_tmdm = omega * (sqrt(permea_tmdm * (loss_expo_tmdm - 1)))
beta_tmdm = omega * (sqrt(permea_tmdm * (loss_expo_tmdm + 1)))

Loss_Path_modified_friis_tmdm = 6.4 + (20*log10(distance*beta_tmdm)) + (8.68 * alph_tmdm * distance) ;
Loss_Path_NC_tmdm = Loss_Path_modified_friis_tmdm + reflect_attenuation_tmdm
Loss_Path_proposed_tmdm = reflect_attenuation_tmdm + 20*log10(distance_boite1 * distance_boite2) - 288.8 + 40*log10(frequency_dielectric) + 20*log10(distance) + 20*log10(beta_tmdm) + (8.68 * alph_tmdm * distance);

%--------------------------------------------------------------------------
%                            Free Space Path Loss
%--------------------------------------------------------------------------
Loss_Path_Free_Space = -147.6 + 20*log10(frequency_dielectric) + 20*log10(distance);
Loss_Path_Free_Space2 = -147.6 + 20*log10(frequency_dielectric) + 20*log10(distance) + abs(reflect_attenuation_mbsdm + refractive_attenuation) ;
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%                            Graphical representation
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

figure(1)
    p1 = plot(distance,Loss_Path_modified_friis_peplinski,'b--', 'LineWidth',1);
    hold on
    p2 = plot(distance,Loss_Path_modified_friis_mbsdm,'g--', 'LineWidth',1);
    p3 = plot(distance,Loss_Path_modified_friis_tmdm,'r--', 'LineWidth',1);
    %p4 = plot(distance,Loss_Path_modified_friis_tdr,'k', 'LineWidth',2);
    p5 = plot(distance,Loss_Path_proposed_peplinski,'y--', 'LineWidth',1);
    p6 = plot(distance,Loss_Path_proposed_mbsdm,'LineWidth',1);
    p7 = plot(distance,Loss_Path_proposed_tmdm,'y.--', 'LineWidth',1);
    grid on
    grid minor
    legend([p1 p2 p3 p5 p6 p7],{' Peplinski ', ' MBSDM ', ' TMDM ', 'Proposed Peplinski', 'Proposed MBSDM', 'Proposed TMDM'},'Location','northwest')
    ylabel('Path Loss (dB)','FontWeight','bold')
    xlabel('Distance (m)','FontWeight','bold')
    title('Modified Friis (Gravelly sand)','FontWeight','bold', 'FontSize', 13)
    
figure(2);
    %p1 = plot(distance,Loss_Path_modified_friis_tdr,'k--', 'LineWidth',1);
    p2 = plot(distance,Loss_Path_modified_friis_mbsdm,'g--', 'LineWidth',1);
    hold on
    %p3 = plot(distance,Loss_Path_modified_friis_tdr - Loss_Path_modified_friis_mbsdm,'b', 'LineWidth',2);
    grid on
    grid minor
    legend(p2,{' MBSDM '},'Location','northwest')
    ylabel('Path Loss (dB)','FontWeight','bold')
    xlabel('Distance (m)','FontWeight','bold')
    title('Modified Friis (Clayey Silt)','FontWeight','bold', 'FontSize', 13)
    
figure(3)
    p1 = plot(distance,Loss_Path_NC_peplinski,'b--', 'LineWidth',1);
    hold on
    p2 = plot(distance,Loss_Path_NC_mbsdm,'g--', 'LineWidth',1);
    p3 = plot(distance,Loss_Path_NC_tmdm,'r--', 'LineWidth',1);
    %p4 = plot(distance,Loss_Path_NC_tdr,'k', 'LineWidth',2);
    grid on
    grid minor
    legend([p1 p2 p3],{' TDR ', ' Peplinski ', ' MBSDM ', ' TMDM '},'Location','northwest')
    ylabel('Path Loss (dB)','FontWeight','bold')
    xlabel('Distance (m)','FontWeight','bold')
    title('NC - Modified Friis (Clayey Silt)','FontWeight','bold', 'FontSize', 13)
    
figure(4)
    %p1 = plot(distance,Loss_Path_modified_friis_tdr, 'k', 'LineWidth',2);
    p2 = plot(distance,Loss_Path_CRIM_peplinski,'--', 'LineWidth',1);
    hold on
    p3 = plot(distance,Loss_Path_NC_peplinski,'--', 'LineWidth',1);
    p4 = plot(distance,Loss_Path_proposed_peplinski,'--', 'LineWidth',1);
    p5 = plot(distance,Loss_Path_modified_friis_mbsdm,'g', 'LineWidth',1);
    grid on
    grid minor
    legend([p2 p3 p4 p5],{' CRIM-Fresnel ', ' Modified Friis ', ' NC Modified Friis ', ' Proposed approach '},'Location','northwest')
    ylabel('Path Loss (dB)','FontWeight','bold')
    xlabel('Distance (m)','FontWeight','bold')
    title('Modified Friis (Clayey Silt)','FontWeight','bold', 'FontSize', 13)
    
figure(5)
    %p1 = plot(distance,Loss_Path_modified_friis_tdr, 'k', 'LineWidth',2);
    p4 = plot(distance,Loss_Path_NC_peplinski,'--', 'LineWidth',1);
    hold on
    p3 = plot(distance,Loss_Path_proposed_peplinski,'--', 'LineWidth',1);
    p5 = plot(distance,Loss_Path_modified_friis_mbsdm,'g', 'LineWidth',1);
    p6 = plot(distance,Loss_Path_NC_mbsdm, 'g--.', 'LineWidth',1);
    grid on
    grid minor
    legend([p3 p4 p5 p6],{' Modified Friis ', ' NC Modified Friis ', ' Proposed approach ', ' Proposed with reflection '},'Location','northwest')
    ylabel('Path Loss (dB)','FontWeight','bold')
    xlabel('Distance (m)','FontWeight','bold')
    moisutre_percent = moisture * 100;
    title_main = sprintf('Path loss comparison (Gravelly sand %d %)', moisutre_percent)
    title(title_main)
    %title('Path loss comparison (Gravelly sand 17.02%)','FontWeight','bold', 'FontSize', 13)
    
figure(6)
    p1 = plot(distance,Loss_Path_modified_friis_peplinski,'r', 'LineWidth',1);
        hold on
    p2 = plot(distance,Loss_Path_NC_peplinski,'b', 'LineWidth',1);
    p3 = plot(distance,Loss_Path_modified_friis_mbsdm,'g', 'LineWidth',1);
    
    plot(3.1 ,112, 'ko', 'LineWidth',1);
    plot(3 ,110, 'ko', 'LineWidth',1);
    
    plot(4.1 ,115, 'ko', 'LineWidth',1);
    plot(4 ,117, 'ko', 'LineWidth',1);
    plot(4 ,113, 'ko', 'LineWidth',1);
    
    plot(5.1 ,140, 'ko', 'LineWidth',1);
    p4 = plot(5 ,137, 'ko', 'LineWidth',1);
    plot(5 ,131, 'ko', 'LineWidth',1);
    
    grid on
    grid minor
    legend([p1 p2 p3 p4],{' Modified Friis ', ' NC Modified Friis ', ' Proposed path loss approach ', ' Measured values'},'Location','northwest')
    ylabel('Path Loss (dB)','FontWeight','bold')
    xlabel('Distance (m)','FontWeight','bold')
    title('Path loss comparison (Sandy clay)','FontWeight','bold', 'FontSize', 13)
hold off
antenna_gain = 2.5;
transmitted_power = 12;
rssi_peplinski = transmitted_power + 2*antenna_gain - Loss_Path_modified_friis_peplinski;
rssi_NC_peplinski = transmitted_power + 2*antenna_gain - Loss_Path_NC_peplinski;
%rssi_tmdm = transmitted_power + 2*antenna_gain - Loss_Path_modified_friis_tmdm;
rssi_mbsdm = transmitted_power + 2*antenna_gain - Loss_Path_modified_friis_mbsdm;

figure(7)
    p1 = plot(distance,rssi_peplinski,'r', 'LineWidth',1);
        hold on
    p2 = plot(distance,rssi_NC_peplinski,'b', 'LineWidth',1);
    %p3 = plot(distance,rssi_tmdm,'k', 'LineWidth',1);
    p3 = plot(distance,rssi_mbsdm,'g', 'LineWidth',1);
    
    %plot(1 ,-79, 'ko', 'LineWidth',1);
    
    %plot(2 ,-90, 'ko', 'LineWidth',1);
     %plot(2.1 ,-96, 'ko', 'LineWidth',1);
    
%    plot(3 ,-94, 'ko', 'LineWidth',1);
    plot(3.1 ,-90, 'ko', 'LineWidth',1);
    plot(3 ,-88, 'ko', 'LineWidth',1);
%     plot(3 ,-94, 'ko', 'LineWidth',1);
%     plot(3 ,-96, 'ko', 'LineWidth',1);
    
%     plot(4 ,-88, 'ko', 'LineWidth',1);
    plot(4.1 ,-93, 'ko', 'LineWidth',1);
%     plot(4 ,-94, 'ko', 'LineWidth',1);
    plot(4 ,-95, 'ko', 'LineWidth',1);
%     plot(4 ,-97, 'ko', 'LineWidth',1);
    plot(4 ,-91, 'ko', 'LineWidth',1);
    
     plot(5.1 ,-118, 'ko', 'LineWidth',1);
%     plot(5 ,-126, 'ko', 'LineWidth',1);
%     plot(5 ,-129, 'ko', 'LineWidth',1);
%     plot(5 ,-130, 'ko', 'LineWidth',1);
    p4 = plot(5 ,-115, 'ko', 'LineWidth',1);
    plot(5 ,-109, 'ko', 'LineWidth',1);
%     plot(5 ,-108, 'ko', 'LineWidth',1);
    %regress = -8.5 * distance - 75.5 ;
    %p5 = plot(distance,regress,'k', 'LineWidth',1);
    %plot(0.1297 ,1.88, 'kx', 'LineWidth',1);
    
    grid on
    grid minor
    legend([p1 p2 p3 p4],{' Modified Friis ', ' NC Modified Friis ', ' Proposed path loss ', ' Measured RSSI'},'Location','northwest')
    ylabel('RSSI (dBm)','FontWeight','bold')
    xlabel('Distance (m)','FontWeight','bold')
    title('Power received (Sandy clay)','FontWeight','bold', 'FontSize', 13)
    hold off

figure(8)
    p1 = plot(distance,Loss_Path_modified_friis_peplinski,'-', 'LineWidth',1);
        hold on
    p2 = plot(distance,Loss_Path_NC_peplinski,'-', 'LineWidth',1);
    p3 = plot(distance,Loss_Path_modified_friis_mbsdm,'-', 'LineWidth',1);
    p4 = plot(distance,Loss_Path_proposed_mbsdm_reflect, 'k-','LineWidth',1);
    p5 = plot(distance,Loss_Path_proposed_mbsdm,'-', 'LineWidth',1);
       
%    plot([0.1297 ,1.88], [5.1297 ,1.96], 'kx', 'LineWidth',1);
%     p5 = plot(distance,Loss_Path_proposed_mbsdm_reflect,'--','color', [0.9290 0.6940 0.1250], 'LineWidth',1);
%     p6 = plot(distance,Loss_Path_proposed_mbsdm_refract,'--','color', [0.4940 0.1840 0.5560], 'LineWidth',1);
%     p7 = plot(distance,Loss_Path_Free_Space,'m','LineWidth',1);
%     p8 = plot(distance,Loss_Path_Free_Space2,'k','LineWidth',1);
    
    %p7 = plot(distance,Loss_Path_proposed_mbsdm_reflect_refract,'--','color', [0.3010 0.7450 0.9330], 'LineWidth',1);
    %p8 = plot(distance,Loss_Path_proposed_mbsdm_reflect_inv,'k--', 'LineWidth',1);
    
%     plot(5 ,[44.45 21.84 59.35 72.78 50.64], 'kx');
%     plot(10 ,[51.24 28.63 80.94 94.37 72.23], 'kx');
%     plot(15 ,[55.53 32.91 99.99 113.4 91.3], 'kx');
%     plot(20 ,[58.78 36.17 118.1 131.5], 'kx');
%    p6 = plot(20 ,109.4, 'kx'); 
    %plot(distance(2),Loss_Path_NC_peplinski,'r*', 'LineWidth',1);
    
    grid on
    grid minor
    legend([p1 p2 p3 p4 p5],{' Modified Friis ', ' NC Modified Friis ', ' MBSDM Modified Friis ', ' Proposed LPL#1 ', ' Proposed LPL#2 '},'Location','northwest')
    ylabel('Path Loss (dB)','FontWeight','bold')
    xlabel('Distance (m)','FontWeight','bold')
    title('UG2UG Path Loss at dry soil (Sandy clay #2)','FontWeight','bold', 'FontSize', 13)
    hold off
end

