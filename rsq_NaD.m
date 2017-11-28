function [ rsq ] = rsq_NaD( s )

h = 6.62607e-27;                % Planck's constant [erg*s]
c = 3e10;                       % Speed of light [cm/s]
erg2eV = 1/1.60219e-12;              % 1 erg to .. eV
E_ioniz = 5.139;

E_n = zeros([1 3]);
E_n(2) = h*c ./ 5895.94e-8 .* erg2eV;
E_n(3) = h*c ./ 5889.97e-8 .* erg2eV;
Z = 1;

Rydberg = 13.6;

l = [0 1 1];
nstar_sq = Rydberg .* Z.^2 ./ (E_ioniz - E_n(s));
rsq = nstar_sq ./ 2 ./ Z.^2 .* (5.*nstar_sq + 1 - 3.*l(s).*(l(s) + 1));

end
