function [Bp,Ipp] = BNC_calculator(mVpp)
%BNC_CALCULATOR calculates the generated magnetic field in the middle of
%the loop for a given Vpp value.

    m = 0.0134;     %Value obtained in Multimeter_measurement_4.m
    h = 3.3685e-05; %Value obtained in Multimeter_measurement_4.m
    
    %for the phantom with the new resistor of 83,3 Ohm we get:
%     m = 0.0073;         %Values obtained in Multimeter_measurement_5.m with measurements on 21/02/2023
%     h = 5.5902e-06;

    Vpp_source = mVpp*1E-3; %[V]
    Vrms_source = (Vpp_source/(2*sqrt(2))); 
    Irms_fit = Vrms_source.*m + h; %[A]
    Ipp = Irms_fit.*(2*sqrt(2)); %[A]

    mw = 4*pi*1E-7;
    r = 8.5e-3; %[m]
    Bp = (((Ipp/2)*mw)/(2*r))*1E9;     % Magnetic field in the center of the loop [nT]. BNC = Bp.sin(....
end

