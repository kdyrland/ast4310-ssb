function [ relnrs ] = Boltz_Na( temp,r,s )
% Boltzmann distribution n_r,s/N_r
keV = 8.61734e-5;               % Boltzmann's constant [eV/K]
%kerg = 1.380658e-16;            % Boltzmann's constant [erg/K]
h = 6.62607e-27;                % Planck's constant [erg*s]
c = 3e10;                       % Speed of light [cm/s]
erg2eV = 6.242e11;              % 1 erg to .. eV

E_n1 = h.*c ./ 5895.94e-8 .* erg2eV;
E_n2 = h.*c ./ 5889.97e-8 .* erg2eV;
u = partfunc_Na(temp);
chi = [0 E_n1 E_n2];
g = [2 2 4];

relnrs = g(s+1)./u(r+1).*exp(-(chi(s+1))./(keV.*temp));

end

