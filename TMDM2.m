%we try to simulate the Temperature and Mineralogy based model for
%predicting the dielectric constant according to the moisture, frequency,
%temperature and clay content

function [epsilon_real_CDC, epsilon_imaginary_CDC] = TMDM2(moisture,temperature,clay)

close all
i = 1;
for frequency = 300000000:1000000:900000000

    temperature_starting = 20; %Starting temperature is set to 10¤ Celsius

    RI_dry = 1.634 - 0.539e-2*clay + 0.2748e-4*clay^2; %RI_dry is the refractive index of dry soil
    NAC_dry = 0.03952 - 0.04038e-2*clay; %NAC_dry is the Normalized Attenuztion Coefficient known as k for dry soil 

    MBWF = 0.02863 + 0.30673e-2*clay; %MBWF is the Maximum Bound Water Fraction
    epsilon_starting_bound = 79.8 - 85.4e-2*clay + 32.7e-4*clay*clay ;
    beta_starting_bound = 8.67e-19 - 0.00126e-2*clay + 0.00184e-4*clay^2 - 9.77e-10*clay^3 - 1.39e-15*clay^4; % Volumetric expansion coefficient of bound water

    epsilon_starting_free = 100;
    beta_starting_free = 1.11e-4 - 1.603e-7*clay + 1.239e-9*clay^2 + 8.33e-13*clay^3 -1.007e-14*clay^4; % Volumetric expansion coefficient of free water

    activation_energy_bound = 1467 + 2697e-2*clay - 980e-4*clay^2 + 1.368e-10*clay^3 - 8.61e-13*clay^4; %Activation energy of bound water over the universal gas constant
    entropy_activation_bound = 0.888 + 9.7-2*clay - 4.262e-4*clay^2 + 6.79e-21*clay^3 +4.263e-22*clay^4; %Zntropy of activation of bound water over the universal gas constant

    activation_energy_free = 2231 - 143.1*10^(-2)*clay + 223.2*10^(-4)*clay^2 - 142.1*10^(-6)*clay^3 + 27.14*10^(-8)*clay^4; %Activation energy of free water over the universal gas constant
    entropy_activation_free = 3.649 - 0.4894*10^(-2)*clay + 0.763*10^(-4)*clay^2 - 0.4859*10^(-6)*clay^3 + 0.0928*10^(-8)*clay^4; %Zntropy of activation of free water over the universal gas constant

    sigma_bound_starting = 0.3112 + 0.467*10^(-2)*clay; %Conductivity of bound water
    beta_temperature_bound = 0.0028 + 0.02094*10^(-2)*clay - 0.01229*10^(-4)*clay^2 - 5.03*10^(-22)*clay^3 + 4.163*10^(-24)*clay^4; %Temperature incrementation coefficient for conductivity of bound water

    sigma_free_starting = 0.05 + 14*(1-(1-clay*10^(-2))^(4.664)); %Conductivity of free water
    beta_temperature_free = 0.00108 + 0.1413e-2*clay - 0.2555e-4*clay^2 + 0.2147e-6*clay^3 - 0.0711e-8*clay^4; %Temperature incrementation coefficient for conductivity of free water

    sigma_bound = sigma_bound_starting + beta_temperature_bound * (temperature - temperature_starting); %Conductivity of bound water at a temperature t
    sigma_free = sigma_free_starting + beta_temperature_free * (temperature - temperature_starting); %Conductivity of free water at a temperature t 
    x1 = temperature + 273.15;

    relaxation_time_bound = (48*10^(-12))/(x1)*exp(activation_energy_bound/x1 - entropy_activation_bound)*10^(-12); %relaxation time of bound water
    relaxation_time_free = (48*10^(-12))/(x1)*exp(activation_energy_free/x1 - entropy_activation_free)*10^(-12); %relaxation time of free water

    f_starting_bound = log((epsilon_starting_bound - 1)/(epsilon_starting_bound +2)); %Computation of the function F in terms of the starting temperature of bound water
    f_starting_free = log((epsilon_starting_free - 1)/(epsilon_starting_free +2)); %Computation of the function F in terms of the starting temperature of free water

    %Calculation of the epsilon zero of bound water at the temperature t
    epsilon_zero_bound = (1 + 2*exp(f_starting_bound - beta_starting_bound*(temperature - temperature_starting)))/(1 - exp(f_starting_bound -  beta_starting_bound*(temperature - temperature_starting))); 

    %Calculation of the epsilon zero of free water at the temperature t
    epsilon_zero_free = (1 + 2*exp(f_starting_free - beta_starting_free*(temperature - temperature_starting)))/(1 - exp(f_starting_free -  beta_starting_free*(temperature - temperature_starting))); 
        
    epsilon_infinite = 4.9;
    epsilon_vaccum = 8.854187e-12;
    
    angular_velocity = 2*pi*frequency*epsilon_vaccum;
    
    %Computation of the CDC of bound water
    CDC_real_bound = epsilon_infinite + ((epsilon_zero_bound - epsilon_infinite)/(1 + (2*pi*frequency*relaxation_time_bound)^2));
    CDC_imaginary_bound = ((2*pi*frequency*relaxation_time_bound*(epsilon_zero_bound - epsilon_infinite))/(1 + (2*pi*frequency*relaxation_time_bound)^2))+((2*pi*frequency*sigma_bound*relaxation_time_bound*relaxation_time_bound)/(epsilon_vaccum*(1 + (2*pi*frequency*relaxation_time_bound)^2)))+((sigma_bound)/(angular_velocity*(1 + (2*pi*frequency*relaxation_time_bound)^2))) ;
    
    %Computation of the CDC of free water
    CDC_real_free = epsilon_infinite + ((epsilon_zero_free - epsilon_infinite)/(1 + (2*pi*frequency*relaxation_time_free)^2));
    CDC_imaginary_free = ((2*pi*frequency*relaxation_time_free*(epsilon_zero_free - epsilon_infinite))/(1 + (2*pi*frequency*relaxation_time_free)^2))+((2*pi*frequency*sigma_free*relaxation_time_free*relaxation_time_free)/(epsilon_vaccum*(1 + (2*pi*frequency*relaxation_time_free)^2)))+((sigma_free)/(angular_velocity*(1 + (2*pi*frequency*relaxation_time_free)^2))) ;
    
    %RI and NAC of bound and free water 
    RI_bound = sqrt(sqrt((CDC_real_bound)^2 + (CDC_imaginary_bound)^2) + CDC_real_bound)/sqrt(2); %RI_bound is the refractive index of bound water
    NAC_bound = sqrt(sqrt((CDC_real_bound)^2 + (CDC_imaginary_bound)^2) - CDC_real_bound)/sqrt(2); %NAC_bound is the Normalized Attenuation Coefficient known as k for bound water 

    RI_free = sqrt(sqrt((CDC_real_free)^2 + (CDC_imaginary_free)^2) + CDC_real_free)/sqrt(2); %RI_free is the refractive index of free water
    NAC_free = sqrt(sqrt((CDC_real_free)^2 + (CDC_imaginary_free)^2) - CDC_real_free)/sqrt(2); %NAC_free is the Normalized Attenuztion Coefficient known as k for free water 
    
   if moisture<MBWF
        RI = RI_dry + (RI_bound - 1)*moisture;
        NAC = NAC_dry + NAC_bound*moisture
    else
        RI = RI_dry + (RI_bound - 1)*MBWF + (RI_free - 1)*(moisture - MBWF);
        NAC = NAC_dry + NAC_bound*MBWF + NAC_free*(moisture - MBWF);
    end
    
    epsilon_real_CDC = RI^2 - NAC^2 ;
    epsilon_imaginary_CDC = 2*RI*NAC ;
    
    %Graphical representation;
    subplot (1,2,1);
    plot(frequency ,epsilon_real_CDC, '-o', 'MarkerIndices', 1:5:length(epsilon_real_CDC));
    ylabel('Dielectric Constant (DC)','FontWeight','bold')
    xlabel('Frequency (Hz)','FontWeight','bold')
    hold on
    
    subplot (1,2,2);
    plot(frequency ,epsilon_imaginary_CDC, '-o', 'MarkerIndices', 1:5:length(epsilon_imaginary_CDC));
    ylabel('Loss Factor (LF)','FontWeight','bold')
    xlabel('Frequency (Hz)','FontWeight','bold')
    hold on
    
    X(i) = epsilon_real_CDC;
    Y(i) = epsilon_imaginary_CDC;
    
    i = i+1;
end
DC = mean(X)
LF = max(Y)

end

