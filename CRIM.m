function [ Loss_Path ] = CRIM(DC, FL)

    distance = linspace(0,5);
    frequency = 868e+6;
    frequency2 = 433000000;
    frequency3 = 915e6;
    
    electrical_conduct = 2.32;
    epsilon_vacuum = 8.854e-12;
    
    alph_num = 8.68*60*pi*(2*pi*frequency2*epsilon_vacuum*DC+electrical_conduct);
    alph_denum = sqrt((DC/2)*(1+(sqrt(1+((FL+(electrical_conduct)/(2*pi*frequency2*epsilon_vacuum))/(DC))^2))));
    alph = alph_num / alph_denum
    
    alph_num2 = 8.68*60*pi*(2*pi*frequency*epsilon_vacuum*DC+electrical_conduct);
    alph_denum2 = sqrt((DC/2)*(1+(sqrt(1+((FL+(electrical_conduct)/(2*pi*frequency*epsilon_vacuum))/(DC))^2))));
    alph2 = alph_num2 / alph_denum2
    
    alph_num3 = 8.68*60*pi*(2*pi*frequency3*epsilon_vacuum*DC+electrical_conduct);
    alph_denum3 = sqrt((DC/2)*(1+(sqrt(1+((FL+(electrical_conduct)/(2*pi*frequency3*epsilon_vacuum))/(DC))^2))));
    alph3 = alph_num3 / alph_denum3
    
    reflect = ((1-sqrt(DC))/(1+sqrt(DC)))^2;
    reflect_attenuation = 10*log((2*reflect)/(1+reflect));
    
    Loss_Path = alph*distance + reflect_attenuation
    
    Loss_Path2 = alph2*distance + reflect_attenuation
    
    Loss_Path3 = alph3*distance + reflect_attenuation
     
    %plot(distance ,Loss_Path, '-o', 'MarkerIndices', 1:5:length(Loss_Path));
    p1 = plot(distance,Loss_Path,'-o','MarkerIndices',1:5:length(Loss_Path))
    hold on
    
    p2 = plot(distance,Loss_Path2,'-*','MarkerIndices',1:5:length(Loss_Path2))
    %hold on
    
    p3 = plot(distance,Loss_Path3,'-+','MarkerIndices',1:5:length(Loss_Path3))
    hold off

    legend([p1 p2 p3],{'433 MHz','868 MHz', '915 MHz'})
    xlabel('Distance (m)','FontWeight','bold')
    ylabel('Path Loss (dB)','FontWeight','bold')

    %plot(distance,Loss_Path,'g',distance,Loss_Path2,'b--o',distance,Loss_Path3,'c*')
    
end

