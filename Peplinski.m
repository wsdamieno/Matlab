%we try to simulate the Peplinski moodel for esstimating the CDC according
%to the Sand and Clay partion, bulk density and the moisture (VWC)

function [epsilon_real_CDC, epsilon_imaginary_CDC] = Peplinski(moisture, sand, clay)

close all
i = 1;
for frequency = 300000000:1000000:900000000

    bulk_density = 1.5 ; % bulk density take from the modified friis model in g/cm3
    epsilon_infinite = 4.9;
    epsilon_vaccuum = 80.1;
    dielectric_vaccum = 8.854187e-12
    relaxation_time = 0.58e-10; %2*pi*relation_time
    particle_density = 2.66; % specific density of the solid soil particles in g/cm3 
    beta_prime = 1.2748 - 0.519*sand - 0.152*clay ;%the value of beta prime which depend of the sand and clay portions
    beta_prime_prime = 1.33797 - 0.603*sand - 0.166*clay ;%the value of beta prime prime which depend of the sand and clay portions
    relative_CDC = ((1.01 + (0.44 * particle_density))^2) - 0.062; % The relative complex dielectric constant
    sigma = 0.62 ;% A predefined constant

    if (300000000 <= frequency) && (frequency < 1400000000)
        effective_conductivity = 0.0467 + (0.2204 * bulk_density) - (0.4111 * sand) + (0.6614 * clay); %Value of the effective conductive in terms of the frequency and the sand, clay portion
    elseif (1400000000 <= frequency) && (frequency < 18000000000)
        effective_conductivity = -1.645 + (1.939 * bulk_density) - (2.25622 * sand) + (1.594 * clay);
    end

    % Real and imaginary parts of the CDC of free water
    epsilon_real_water = epsilon_infinite + ((epsilon_vaccuum - epsilon_infinite)/(1 + (relaxation_time * frequency)^2));
    epsilon_imaginary_water = ((relaxation_time * frequency * (epsilon_vaccuum - epsilon_infinite))/(1 + (relaxation_time * frequency)^2)) + ((effective_conductivity) * (particle_density - bulk_density)) / (2 * pi * frequency * dielectric_vaccum * particle_density * moisture);
    
    %Ratio between the bulk density and tge particle density
    density_particle = bulk_density/particle_density;
    
    %Real and imaginary parts of the complex mixing dielectric constant of soil-water
    if (300000000 <= frequency) && (frequency < 1390000000)
        epsilon_real_CDC = 1.15 * (1 + (density_particle * ((relative_CDC)^(sigma) - 1)) + moisture ^(beta_prime) * epsilon_real_water ^(sigma) - moisture)^(1/sigma) - 0.68
    elseif (1400000000 <= frequency) && (frequency < 18000000000)
        epsilon_real_CDC = (1 + (density_particle*((relative_CDC)^(sigma) - 1)) + moisture^(beta_prime)*epsilon_real_water^(sigma) - moisture)^(1/sigma)
    end

    epsilon_imaginary_CDC = (((moisture) ^(beta_prime_prime)) * (epsilon_imaginary_water) ^(sigma))^(1/sigma)
    
    % graphical representation of the peplinski model
    subplot (1,2,1);
    plot(frequency ,epsilon_real_CDC, '-o', 'MarkerIndices', 1:5:length(epsilon_real_CDC));
    ylabel('Dielectric Constant (DC)','FontWeight','bold')
    xlabel('Frequency (Hz)','FontWeight','bold')
    hold on
    
    subplot (1,2,2);
    plot(frequency ,epsilon_imaginary_CDC, '-o', 'MarkerIndices', 1:5:length(epsilon_imaginary_CDC));
    xlabel('Frequency (Hz)','FontWeight','bold')
    ylabel('Loss Factor (FL)','FontWeight','bold')
    hold on
    
    X(i) = epsilon_real_CDC;
    Y(i) = epsilon_imaginary_CDC;
    
    i = i+1;

end

frequency_mean = 900e+6;

epsilon_real_water_mean = epsilon_infinite + ((epsilon_vaccuum - epsilon_infinite)/(1 + (relaxation_time * frequency_mean)^2));
epsilon_imaginary_water_mean = ((relaxation_time * frequency_mean * (epsilon_vaccuum - epsilon_infinite))/(1 + (relaxation_time * frequency_mean)^2)) + ((effective_conductivity) * (particle_density - bulk_density)) / (2 * pi * frequency_mean * dielectric_vaccum * particle_density * moisture);
  
epsilon_real_CDC_mean = 1.15 * (1 + (density_particle * ((relative_CDC)^(sigma) - 1)) + moisture ^(beta_prime) * epsilon_real_water_mean ^(sigma) - moisture)^(1/sigma) - 0.68
epsilon_imaginary_CDC_mean = (((moisture) ^(beta_prime_prime)) * (epsilon_imaginary_water_mean) ^(sigma))^(1/sigma)
          
DC = max(X)
LF = max(Y)

end

