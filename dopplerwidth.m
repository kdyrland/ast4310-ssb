function [doppler] = dopplerwidth(wav,temp,v_t,m)
% Takes in central wavelength in cm,temperature in K, v_t in km/s ...
% and m in grams and returns dopplerwidth in cm
kerg = 1.380658e-16;            % Boltzmann's constant [erg/K]
c = 3e10;                       % Speed of light [cm/s]

doppler = wav./c.*sqrt(2.*kerg.*temp./m + v_t.*v_t.*1e10);
 
end