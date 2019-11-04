function [ ] = TDR_modified_friis2( moisture, sand, clay )

% the fixed operating frequency is 433MHz
frequency_dielectric = 433e+6;
temperature = 20 ;


%Distance varies from 0 to 20m
distance2 = 0:0.000001:0.50;
distance = distance2/cosd(15) ;
distance_boite1 = 0.13;
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
%--------------------------------------------------------------------------
%                Path Loss prediction with the Path losses
%--------------------------------------------------------------------------
mu_factor = pi * 4 * 10^(-7) * 8.854 * 10^(-12);
omega = 2*pi*frequency_dielectric;    
% distance_boite1 = 0.12;
% distance_boite2 = 0.12;
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

distance_boite2 = 5;
distance_boite2_10 = 10;
distance_boite2_15 = 15;
distance_boite2_20 = 20;

%--------------------------------------------------------------------------
%                        MBSDM permittivity prediction
%--------------------------------------------------------------------------
permea_mbsdm = (mu_factor*epsilon_real_CDC_mbsdm)/2 ;
reflect_mbsdm = ((1-sqrt(epsilon_real_CDC_mbsdm))/(1+sqrt(epsilon_real_CDC_mbsdm)))^2;
reflect_attenuation_mbsdm = 10*log((2*reflect_mbsdm)/(1+reflect_mbsdm));

refractive = sqrt((sqrt(epsilon_real_CDC_mbsdm^2 + epsilon_imaginary_CDC_mbsdm^2) + epsilon_real_CDC_mbsdm) / 2)
refractive_attenuation = 20 * log10((refractive + 1)/4)

loss_expo_mbsdm = sqrt(1+((epsilon_imaginary_CDC_mbsdm)/(epsilon_real_CDC_mbsdm))^2);

alph_mbsdm = omega * (sqrt(permea_mbsdm * (loss_expo_mbsdm - 1)))
beta_mbsdm = omega * (sqrt(permea_mbsdm * (loss_expo_mbsdm + 1)))


% Loss_Path_modified_friis_mbsdm = 6.4 + (20*log10(distance*beta_mbsdm)) + (8.68 * alph_mbsdm * distance) ;
% Loss_Path_NC_mbsdm = reflect_attenuation_mbsdm + 6.4 + (20*log10(distance*beta_mbsdm)) + (8.68 * alph_mbsdm * distance) ;
Loss_Path_proposed_mbsdm = 20*log10(distance_boite1 * distance_boite2 * distance * beta_mbsdm) - 288.8 + 40*log10(frequency_dielectric) + (8.68 * alph_mbsdm * distance)
Loss_Path_proposed_mbsdm_refrac = refractive_attenuation + 20*log10(distance_boite1 * distance_boite2 * distance * beta_mbsdm) - 288.8 + 40*log10(frequency_dielectric) + (8.68 * alph_mbsdm * distance)

Loss_Path_proposed_mbsdm_reflect = reflect_attenuation_mbsdm + 20*log10(distance_boite1 * distance_boite2 * distance * beta_mbsdm) - 288.8 + 40*log10(frequency_dielectric) + (8.68 * alph_mbsdm * distance)
Loss_Path_proposed_mbsdm_reflect_refrac = refractive_attenuation + reflect_attenuation_mbsdm + 20*log10(distance_boite1 * distance_boite2 * distance * beta_mbsdm) - 288.8 + 40*log10(frequency_dielectric) + (8.68 * alph_mbsdm * distance)


Loss_Path_proposed_mbsdm_10 = 20*log10(distance_boite1 * distance_boite2_10 * distance * beta_mbsdm) - 288.8 + 40*log10(frequency_dielectric) + (8.68 * alph_mbsdm * distance)
Loss_Path_proposed_mbsdm_10_refrac = refractive_attenuation + 20*log10(distance_boite1 * distance_boite2_10 * distance * beta_mbsdm) - 288.8 + 40*log10(frequency_dielectric) + (8.68 * alph_mbsdm * distance)

Loss_Path_proposed_mbsdm_reflect_10 = reflect_attenuation_mbsdm + 20*log10(distance_boite1 * distance_boite2_10 * distance * beta_mbsdm) - 288.8 + 40*log10(frequency_dielectric) + (8.68 * alph_mbsdm * distance)
Loss_Path_proposed_mbsdm_reflect_10_refrac = refractive_attenuation + reflect_attenuation_mbsdm + 20*log10(distance_boite1 * distance_boite2_10 * distance * beta_mbsdm) - 288.8 + 40*log10(frequency_dielectric) + (8.68 * alph_mbsdm * distance)

Loss_Path_proposed_mbsdm_15 = 20*log10(distance_boite1 * distance_boite2_15 * distance * beta_mbsdm) - 288.8 + 40*log10(frequency_dielectric) + (8.68 * alph_mbsdm * distance)
Loss_Path_proposed_mbsdm_15_refrac = refractive_attenuation + 20*log10(distance_boite1 * distance_boite2_15 * distance * beta_mbsdm) - 288.8 + 40*log10(frequency_dielectric) + (8.68 * alph_mbsdm * distance)

Loss_Path_proposed_mbsdm_reflect_15 = reflect_attenuation_mbsdm + 20*log10(distance_boite1 * distance_boite2_15 * distance * beta_mbsdm) - 288.8 + 40*log10(frequency_dielectric) + (8.68 * alph_mbsdm * distance)
Loss_Path_proposed_mbsdm_reflect_15_refrac = refractive_attenuation + reflect_attenuation_mbsdm + 20*log10(distance_boite1 * distance_boite2_15 * distance * beta_mbsdm) - 288.8 + 40*log10(frequency_dielectric) + (8.68 * alph_mbsdm * distance)

Loss_Path_proposed_mbsdm_20 = 20*log10(distance_boite1 * distance_boite2_20 * distance * beta_mbsdm) - 288.8 + 40*log10(frequency_dielectric) + (8.68 * alph_mbsdm * distance)
Loss_Path_proposed_mbsdm_20_refrac = refractive_attenuation + 20*log10(distance_boite1 * distance_boite2_20 * distance * beta_mbsdm) - 288.8 + 40*log10(frequency_dielectric) + (8.68 * alph_mbsdm * distance)

Loss_Path_proposed_mbsdm_reflect_20 = reflect_attenuation_mbsdm + 20*log10(distance_boite1 * distance_boite2_20 * distance * beta_mbsdm) - 288.8 + 40*log10(frequency_dielectric) + (8.68 * alph_mbsdm * distance)
Loss_Path_proposed_mbsdm_reflect_20_refrac = refractive_attenuation + reflect_attenuation_mbsdm + 20*log10(distance_boite1 * distance_boite2_20 * distance * beta_mbsdm) - 288.8 + 40*log10(frequency_dielectric) + (8.68 * alph_mbsdm * distance)

% Loss_Path_proposed_mbsdm_reflect_inv = 20*log10(distance_boite1 * distance_boite2 * distance * beta_mbsdm) - refractive_attenuation - 288.8 + 40*log10(frequency_dielectric) + (8.68 * alph_mbsdm * distance)
% Loss_Path_proposed_mbsdm_refract =  refractive_attenuation + 20*log10(distance_boite1 * distance_boite2 * distance * beta_mbsdm) - 288.8 + 40*log10(frequency_dielectric) + (8.68 * alph_mbsdm * distance)
% Loss_Path_proposed_mbsdm_reflect_refract =  refractive_attenuation + reflect_attenuation_mbsdm + 20*log10(distance_boite1 * distance_boite2 * distance * beta_mbsdm) - 288.8 + 40*log10(frequency_dielectric) + (8.68 * alph_mbsdm * distance)

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%                            Graphical representation
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%% UG2AG COMMUNICATION
%figure(1)
%     subplot (2,2,1);
%     p1 = plot(distance2,Loss_Path_proposed_mbsdm_reflect, '-','LineWidth',1);
%         hold on
%     p2 = plot(distance2,Loss_Path_proposed_mbsdm, '-','LineWidth',1);
%      plot(0.15 ,[38.51 60.65], 'kx');
%      plot(0.20 ,[41.17 63.31], 'kx');
%      plot(0.30 ,[45.01 67.15], 'kx');
%      plot(0.40 ,[47.83 69.97], 'kx');
%     grid on
%     grid minor
%     legend([p1 p2],{' Proposed WUSN-PLM #1 ', ' Proposed WUSN-PLM #2 '},'Location','northwest')
%     ylabel('Path Loss (dB)','FontWeight','bold')
%     %xlabel('Underground distance (m)','FontWeight','bold')
%     title('5m UG2AG','FontWeight','bold', 'FontSize', 13)
%     hold off
%     
%     subplot (2,2,2);
%     p1 = plot(distance2,Loss_Path_proposed_mbsdm_reflect_10, '-','LineWidth',1);
%         hold on
%     p2 = plot(distance2,Loss_Path_proposed_mbsdm_10, '-','LineWidth',1);
%     
%      plot(0.15 ,[44.53 66.66], 'kx');
%      plot(0.20 ,[47.19 69.33], 'kx');
%      plot(0.30 ,[51.03 73.17], 'kx');
%      plot(0.40 ,[53.85 75.99], 'kx');
%     grid on
%     grid minor
%     legend([p1 p2],{' Proposed WUSN-PLM #1 ', ' Proposed WUSN-PLM #2 '},'Location','northwest')
%     %ylabel('Path Loss (dB)','FontWeight','bold')
%     %xlabel('Underground distance (m)','FontWeight','bold')
%     title('10m UG2AG','FontWeight','bold', 'FontSize', 13)
%     hold off
%     
%     subplot (2,2,3);
%     p1 = plot(distance2,Loss_Path_proposed_mbsdm_reflect_15, '-','LineWidth',1);
%         hold on
%     p2 = plot(distance2,Loss_Path_proposed_mbsdm_15, '-','LineWidth',1);
%      plot(0.15 ,[48.05 70.18], 'kx');
%      plot(0.20 ,[50.71 72.84], 'kx');
%      plot(0.30 ,[54.55 76.69], 'kx');
%      plot(0.40 ,[57.38 79.51], 'kx');
%     grid on
%     grid minor
%     legend([p1 p2],{' Proposed WUSN-PLM #1 ', ' Proposed WUSN-PLM #2 '},'Location','northwest')
%     ylabel('Path Loss (dB)','FontWeight','bold')
%     xlabel('Burial depth (m)','FontWeight','bold')
%     title('15m UG2AG','FontWeight','bold', 'FontSize', 13)
%     hold off
%     
%     subplot (2,2,4);
%     p1 = plot(distance2,Loss_Path_proposed_mbsdm_reflect_20, '-','LineWidth',1);
%         hold on
%     p2 = plot(distance2,Loss_Path_proposed_mbsdm_20, '-','LineWidth',1);
%      plot(0.15 ,[50.55 72.68], 'kx');
%      plot(0.20 ,[53.21 75.34], 'kx');
%      plot(0.30 ,[57.05 79.19], 'kx');
%      plot(0.40 ,[59.87 82.01], 'kx');
%     grid on
%     grid minor
%     legend([p1 p2],{' Proposed WUSN-PLM #1 ', ' Proposed WUSN-PLM #2 '},'Location','northwest')
%     %ylabel('Path Loss (dB)','FontWeight','bold')
%     xlabel('Burial depth (m)','FontWeight','bold')
%     title('20m UG2AG','FontWeight','bold', 'FontSize', 13)
%     hold off
%    

%% PLOT OF AG2UG COMMUNICATION TYPE %%
% figure(1)
%     subplot (2,2,1);
%     p1 = plot(distance2,Loss_Path_proposed_mbsdm_reflect_refrac, '-','LineWidth',1);
%         hold on
%     p2 = plot(distance2,Loss_Path_proposed_mbsdm_refrac, '-','LineWidth',1);
%      plot(0.15 ,[34.88 57.01], 'kx');
%      plot(0.20 ,[37.54 59.67], 'kx');
%      plot(0.30 ,[41.38 63.52], 'kx');
%      plot(0.40 ,[44.2 66.34], 'kx');
%      
%     grid on
%     grid minor
%     legend([p1 p2],{' Proposed WUSN-PLM #1 ', ' Proposed WUSN-PLM #2 '},'Location','northwest')
%     ylabel('Path Loss (dB)','FontWeight','bold')
%     %xlabel('Underground distance (m)','FontWeight','bold')
%     title('5m AG2UG','FontWeight','bold', 'FontSize', 13)
%     hold off
%     
%     subplot (2,2,2);
%     p1 = plot(distance2,Loss_Path_proposed_mbsdm_reflect_10_refrac, '-','LineWidth',1);
%         hold on
%     p2 = plot(distance2,Loss_Path_proposed_mbsdm_10_refrac, '-','LineWidth',1);
%     
%      plot(0.15 ,[40.9 63.03], 'kx');
%      plot(0.20 ,[43.56 65.69], 'kx');
%      plot(0.30 ,[47.4 69.54], 'kx');
%      plot(0.40 ,[50.22 72.36], 'kx');
%     grid on
%     grid minor
%     legend([p1 p2],{' Proposed WUSN-PLM #1 ', ' Proposed WUSN-PLM #2 '},'Location','northwest')
%     %ylabel('Path Loss (dB)','FontWeight','bold')
%     %xlabel('Underground distance (m)','FontWeight','bold')
%     title('10m AG2UG','FontWeight','bold', 'FontSize', 13)
%     hold off
%     
%     subplot (2,2,3);
%     p1 = plot(distance2,Loss_Path_proposed_mbsdm_reflect_15_refrac, '-','LineWidth',1);
%         hold on
%     p2 = plot(distance2,Loss_Path_proposed_mbsdm_15_refrac, '-','LineWidth',1);
%      plot(0.15 ,[44.42 66.56], 'kx');
%      plot(0.20 ,[47.08 69.22], 'kx');
%      plot(0.30 ,[50.92 73.06], 'kx');
%      plot(0.40 ,[53.74 75.88], 'kx');
%     grid on
%     grid minor
%     legend([p1 p2],{' Proposed WUSN-PLM #1 ', ' Proposed WUSN-PLM #2 '},'Location','northwest')
%     ylabel('Path Loss (dB)','FontWeight','bold')
%     xlabel('Burial depth (m)','FontWeight','bold')
%     title('15m AG2UG','FontWeight','bold', 'FontSize', 13)
%     hold off
%     
%     subplot (2,2,4);
%     p1 = plot(distance2,Loss_Path_proposed_mbsdm_reflect_20_refrac, '-','LineWidth',1);
%         hold on
%     p2 = plot(distance2,Loss_Path_proposed_mbsdm_20_refrac, '-','LineWidth',1);
%      plot(0.15 ,[46.92 69.06], 'kx');
%      plot(0.20 ,[49.58 71.72], 'kx');
%      plot(0.30 ,[53.42 75.56], 'kx');
%      plot(0.40 ,[56.24 78.38], 'kx');
%     grid on
%     grid minor
%     legend([p1 p2],{' Proposed WUSN-PLM #1 ', ' Proposed WUSN-PLM #2 '},'Location','northwest')
%     %ylabel('Path Loss (dB)','FontWeight','bold')
%     xlabel('Burial depth (m)','FontWeight','bold')
%     title('20m AG2UG','FontWeight','bold', 'FontSize', 13)
%     hold off
    
%% figure(2)
%     p1 = plot(distance2,Loss_Path_proposed_mbsdm_reflect, '-','LineWidth',1);
%         hold on
%     p2 = plot(distance2,Loss_Path_proposed_mbsdm, '-','LineWidth',1);
%     grid on
%     grid minor
%     legend([p1 p2],{' Proposed LPL#1 ', ' Proposed LPL#2 '},'Location','northwest')
%     ylabel('Path Loss (dB)','FontWeight','bold')
%     xlabel('Underground distance (m)','FontWeight','bold')
%     title('5m UG2AG','FontWeight','bold', 'FontSize', 13)
%     hold off
%     
% figure(3)
%     p1 = plot(distance2,Loss_Path_proposed_mbsdm_reflect_10, '-','LineWidth',1);
%         hold on
%     p2 = plot(distance2,Loss_Path_proposed_mbsdm_10, '-','LineWidth',1);
%     grid on
%     grid minor
%     legend([p1 p2],{' Proposed LPL#1 ', ' Proposed LPL#2 '},'Location','northwest')
%     ylabel('Path Loss (dB)','FontWeight','bold')
%     xlabel('Underground distance (m)','FontWeight','bold')
%     title('10m UG2AG','FontWeight','bold', 'FontSize', 13)
%     hold off
%     
% figure(4)
%     p1 = plot(distance2,Loss_Path_proposed_mbsdm_reflect_15, '-','LineWidth',1);
%         hold on
%     p2 = plot(distance2,Loss_Path_proposed_mbsdm_15, '-','LineWidth',1);
%     grid on
%     grid minor
%     legend([p1 p2],{' Proposed LPL#1 ', ' Proposed LPL#2 '},'Location','northwest')
%     ylabel('Path Loss (dB)','FontWeight','bold')
%     xlabel('Burial depth (m)','FontWeight','bold')
%     title('15m UG2AG','FontWeight','bold', 'FontSize', 13)
%     hold off
%     
% figure(5)
%     p1 = plot(distance2,Loss_Path_proposed_mbsdm_reflect_20, '-','LineWidth',1);
%         hold on
%     p2 = plot(distance2,Loss_Path_proposed_mbsdm_20, '-','LineWidth',1);
%     grid on
%     grid minor
%     legend([p1 p2],{' Proposed LPL''#1 ', ' Proposed LPL''#2 '},'Location','northwest')
%     ylabel('Path Loss (dB)','FontWeight','bold')
%     xlabel('Burial depth (m)','FontWeight','bold')
%     title('20m AG2UG','FontWeight','bold', 'FontSize', 13)
%     hold off
    
% figure(6)
%     subplot (1,2,1);
%     delivery_rate = [100 81.87; 100 100; 100 100; 100 100];
%     diag_bar = bar(delivery_rate, 'EdgeColor','k');
%     grid on
%     grid minor
%     diag_bar(1).FaceColor = [0.9290 0.6940 0.1250];    
%     diag_bar(2).FaceColor = [0.4660 0.6740 0.1880];
%     set(gca,'xticklabel',{'15cm-ground','20cm-ground','30cm-ground','40cm-ground'},'FontWeight','bold');
%     ylabel('Date delivery rate (%)','FontWeight','bold')
%     legend([diag_bar(1) diag_bar(2)],{' UG -> AG ', ' AG -> UG '},'Location','northwest')
%     title('Data delivery rate UG2AG - AG2UG (5m)','FontWeight','bold', 'FontSize', 13)
%     hold off
%     
%     subplot (1,2,2);
%     delivery_rate = [97.08 76.1; 95.32 97.64; 97.66 92.39];
%     diag_bar = bar(delivery_rate, 'EdgeColor','k');
%     grid on
%     grid minor
%     diag_bar(1).FaceColor = [0.9290 0.6940 0.1250];    
%     diag_bar(2).FaceColor = [0.4660 0.6740 0.1880];
%     set(gca,'xticklabel',{'15cm-ground','20cm-ground','40cm-ground'},'FontWeight','bold');
%     ylabel('Date delivery rate (%)','FontWeight','bold')
%     legend([diag_bar(1) diag_bar(2)],{' UG -> AG ', ' AG -> UG '},'Location','northwest')
%     title('Data delivery rate UG2AG - AG2UG (15m)','FontWeight','bold', 'FontSize', 13)
%     hold off
%    

%% PLOT OF UG2AG AND AG2UG COMMUNICATION ON A SAME PLOT 

figure(1)
    subplot (2,2,1);
    p1 = plot(distance2,Loss_Path_proposed_mbsdm_reflect_refrac, '--','LineWidth',1);
        hold on
    p2 = plot(distance2,Loss_Path_proposed_mbsdm_refrac, '--','LineWidth',1);
     plot(0.15 ,[34.88 57.01], 'kx');
     plot(0.20 ,[37.54 59.67], 'kx');
     plot(0.30 ,[41.38 63.52], 'kx');
     plot(0.40 ,[44.2 66.34], 'kx');
     
     p1_1 = plot(distance2,Loss_Path_proposed_mbsdm_reflect, '-','LineWidth',1);
        hold on
     p2_1 = plot(distance2,Loss_Path_proposed_mbsdm, '-','LineWidth',1);
     plot(0.15 ,[38.51 60.65], 'k+');
     plot(0.20 ,[41.17 63.31], 'k+');
     plot(0.30 ,[45.01 67.15], 'k+');
     plot(0.40 ,[47.83 69.97], 'k+');
    grid on
    grid minor
    legend([p1_1 p1 p2_1 p2],{' WUSN-PLM #1 (UG2AG) ',' WUSN-PLM #1 (AG2UG) ', ' WUSN-PLM #2 (UG2AG) ', ' WUSN-PLM #2 (AG2UG) '},'Location','northwest')
    ylabel('Path Loss (dB)','FontWeight','bold')
    %xlabel('Underground distance (m)','FontWeight','bold')
    title(' 5m ','FontWeight','bold', 'FontSize', 13)
    hold off
    
    subplot (2,2,2);
    p1 = plot(distance2,Loss_Path_proposed_mbsdm_reflect_10_refrac, '--','LineWidth',1);
        hold on
    p2 = plot(distance2,Loss_Path_proposed_mbsdm_10_refrac, '--','LineWidth',1);
    
     plot(0.15 ,[40.9 63.03], 'kx');
     plot(0.20 ,[43.56 65.69], 'kx');
     plot(0.30 ,[47.4 69.54], 'kx');
     plot(0.40 ,[50.22 72.36], 'kx');
     
     p1_1 = plot(distance2,Loss_Path_proposed_mbsdm_reflect_10, '-','LineWidth',1);
        hold on
     p2_1 = plot(distance2,Loss_Path_proposed_mbsdm_10, '-','LineWidth',1);
    
     plot(0.15 ,[44.53 66.66], 'k+');
     plot(0.20 ,[47.19 69.33], 'k+');
     plot(0.30 ,[51.03 73.17], 'k+');
     plot(0.40 ,[53.85 75.99], 'k+');

    grid on
    grid minor
    legend([p1_1 p1 p2_1 p2],{' WUSN-PLM #1 (UG2AG) ',' WUSN-PLM #1 (AG2UG) ', ' WUSN-PLM #2 (UG2AG) ', ' WUSN-PLM #2 (AG2UG) '},'Location','northwest')
    %ylabel('Path Loss (dB)','FontWeight','bold')
    %xlabel('Underground distance (m)','FontWeight','bold')
    title(' 10m ','FontWeight','bold', 'FontSize', 13)
    hold off
    
    subplot (2,2,3);
    p1 = plot(distance2,Loss_Path_proposed_mbsdm_reflect_15_refrac, '--','LineWidth',1);
        hold on
    p2 = plot(distance2,Loss_Path_proposed_mbsdm_15_refrac, '--','LineWidth',1);
     plot(0.15 ,[44.42 66.56], 'kx');
     plot(0.20 ,[47.08 69.22], 'kx');
     plot(0.30 ,[50.92 73.06], 'kx');
     plot(0.40 ,[53.74 75.88], 'kx');
     
     p1_1 = plot(distance2,Loss_Path_proposed_mbsdm_reflect_15, '-','LineWidth',1);
        hold on
     p2_1 = plot(distance2,Loss_Path_proposed_mbsdm_15, '-','LineWidth',1);
     plot(0.15 ,[48.05 70.18], 'k+');
     plot(0.20 ,[50.71 72.84], 'k+');
     plot(0.30 ,[54.55 76.69], 'k+');
     plot(0.40 ,[57.38 79.51], 'k+');
    
    grid on
    grid minor
    legend([p1_1 p1 p2_1 p2],{' WUSN-PLM #1 (UG2AG) ',' WUSN-PLM #1 (AG2UG) ', ' WUSN-PLM #2 (UG2AG) ', ' WUSN-PLM #2 (AG2UG) '},'Location','northwest')
    ylabel('Path Loss (dB)','FontWeight','bold')
    xlabel('Burial depth (m)','FontWeight','bold')
    title(' 15m ','FontWeight','bold', 'FontSize', 13)
    hold off
    
    subplot (2,2,4);
    p1 = plot(distance2,Loss_Path_proposed_mbsdm_reflect_20_refrac, '--','LineWidth',1);
        hold on
    p2 = plot(distance2,Loss_Path_proposed_mbsdm_20_refrac, '--','LineWidth',1);
     plot(0.15 ,[46.92 69.06], 'kx');
     plot(0.20 ,[49.58 71.72], 'kx');
     plot(0.30 ,[53.42 75.56], 'kx');
     plot(0.40 ,[56.24 78.38], 'kx');
     
     p1_1 = plot(distance2,Loss_Path_proposed_mbsdm_reflect_20, '-','LineWidth',1);
        hold on
     p2_1 = plot(distance2,Loss_Path_proposed_mbsdm_20, '-','LineWidth',1);
     plot(0.15 ,[50.55 72.68], 'k+');
     plot(0.20 ,[53.21 75.34], 'k+');
     plot(0.30 ,[57.05 79.19], 'k+');
     plot(0.40 ,[59.87 82.01], 'k+');
    grid on
    grid minor
    legend([p1_1 p1 p2_1 p2],{' WUSN-PLM #1 (UG2AG) ',' WUSN-PLM #1 (AG2UG) ', ' WUSN-PLM #2 (UG2AG) ', ' WUSN-PLM #2 (AG2UG) '},'Location','northwest')
    %ylabel('Path Loss (dB)','FontWeight','bold')
    xlabel('Burial depth (m)','FontWeight','bold')
    title(' 20m ','FontWeight','bold', 'FontSize', 13)
    hold off
end

