function [epsilon_real_CDC, epsilon_imaginary_CDC] = MBSDM(moisture,clay)
close all
i = 1;
for frequency = 300000000:1000000:900000000
    
    
    RI_dry = 1.634 - (0.539e-2)*clay + (0.2748e-4)*clay^2; %RI_dry is the refractive index of dry soil
    NAC_dry = 0.03952 - (0.04038e-2)*clay; %NAC_dry is the Normalized Attenuztion Coefficient known as k for dry soil 
    MBWF = 0.02863 + (0.30673e-2)*clay; %MBWF is the Maximum Bound Water Fraction
    epsilon_infinite = 4.9;
    epsilon_vaccum = 8.854187e-12;

    epsilon_starting_bound = 79.8 - (85.4e-2)*clay + (32.7e-4)*clay^2 ;
    relaxation_bound = 1.062e-11 + (3.450e-12)*10^(-2)*clay; % relaxation time of bound water
    sigma_bound = 0.3112 + (0.467e-2)*clay; %Conductivity of bound water

    epsilon_starting_free = 100;
    relaxation_free = 8.5e-12 ; % relaxation time of free water
    sigma_free = 0.3631 + (1.217e-2)*clay; %Conductivity of free water

    CDC_real_bound = epsilon_infinite + (epsilon_starting_bound - epsilon_infinite)/(1+(2*pi*frequency*relaxation_bound)^2);%real part of the CDC of bound water 
    CDC_imaginary_bound = (((epsilon_starting_bound - epsilon_infinite)*(2*pi*frequency*relaxation_bound))/(1+(2*pi*frequency*relaxation_bound)^2)) + ((sigma_bound)/(2*pi*frequency*epsilon_vaccum));%imaginary part of the CDC of bound water 

    CDC_real_free = epsilon_infinite + (epsilon_starting_free - epsilon_infinite)/(1+(2*pi*frequency*relaxation_free)^2);%real part of the CDC of free water 
    CDC_imaginary_free = (((epsilon_starting_free - epsilon_infinite)*(2*pi*frequency*relaxation_free))/(1+(2*pi*frequency*relaxation_free)^2)) + ((sigma_free)/(2*pi*frequency*epsilon_vaccum));%imaginary part of the CDC of free water 

    RI_bound = sqrt(sqrt((CDC_real_bound)^2 + (CDC_imaginary_bound)^2) + CDC_real_bound)/sqrt(2); %RI_bound is the refractive index of bound water
    NAC_bound = sqrt(sqrt((CDC_real_bound)^2 + (CDC_imaginary_bound)^2) - CDC_real_bound)/sqrt(2); %NAC_bound is the Normalized Attenuation Coefficient known as k for bound water 

    RI_free = sqrt(sqrt((CDC_real_free)^2 + (CDC_imaginary_free)^2) + CDC_real_free)/sqrt(2); %RI_free is the refractive index of free water
    NAC_free = sqrt(sqrt((CDC_real_free)^2 + (CDC_imaginary_free)^2) - CDC_real_free)/sqrt(2); %NAC_free is the Normalized Attenuztion Coefficient known as k for free water 

    if moisture<MBWF
        RI = RI_dry +(RI_bound - 1)*moisture ; %resulting Refractive index
        NAC = NAC_dry + NAC_bound*moisture ; %resulting Normalized attenuation 
    else
        RI = RI_dry +(RI_bound - 1)*MBWF + (RI_free - 1)*(moisture - MBWF); %resulting Refractive index
        NAC = NAC_dry + NAC_bound*MBWF + NAC_free*(moisture - MBWF); %resulting Normalized attenuation
    end

    epsilon_real_CDC = RI^2 - NAC^2  %Real part of the mixing CDC
    epsilon_imaginary_CDC = 2*RI*NAC %Imaginary part of the mixing CDC

    %Graphical representation;
       
    subplot (1,2,1);hold on
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
DC = max(X)
LF = max(Y)

end

