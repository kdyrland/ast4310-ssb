function [ B ] = planck( temp,wav )
% Planck function

kerg = 1.380658e-16;        % Boltzmann's constant [erg K]
h = 6.62607e-27;            % Planck's constant [erg s]
kT = kerg*temp;
c = 2.99792e10;             % Speed of light [cm/s]

B = 2*h*c^2 ./ (wav.^5.*(exp(h*c./(wav*kT)) - 1));  %[erg/cm^2/s/cm/steradian]

end

