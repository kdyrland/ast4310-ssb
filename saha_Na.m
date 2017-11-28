function [ nstagerel ] = saha_Na( temp,eldens,ionstage )

keV = 8.61734e-5;               % Boltzmann's constant [eV/K]
kerg = 1.380658e-16;            % Boltzmann's constant [erg/K]
h = 6.62607e-27;                % Planck's constant [erg*s]
m_e = 9.10939e-28;              % [g]
E_ioniz = 5.139;

keVT = keV*temp;
kergT =kerg*temp;

u = partfunc_Na(temp);
u(end+1) = 2;       % append element to array

sahaconst = (2.*pi .*m_e*kergT./(h.^2)).^(3./2).*2./eldens;
nstage = zeros([1 4]);
nstage(1) = 1;

for r = 1:3
    nstage(r+1) =nstage(r) .* sahaconst .*u(r+1)./u(r) .* exp(-E_ioniz./keVT);	
    ntotal = sum(nstage);
    nstagerel = nstage./ntotal;
end
nstagerel = nstagerel(ionstage+1);
end


