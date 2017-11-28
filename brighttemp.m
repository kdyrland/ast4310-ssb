function [ Tb ] = brighttemp( B,wav )
% Temperature brightness

kerg = 1.380658e-16;        % Boltzmann's constant [erg K]
h = 6.62607e-27;            % Planck's constant [erg s]
c = 2.99792e10;             % Speed of light [cm/s]
hc = h*c;

Tb = (hc./(wav.*kerg) ./ log((2*h*c^2) ./ (wav.^5 .*B.*1e10) + 1)); %[K]

end