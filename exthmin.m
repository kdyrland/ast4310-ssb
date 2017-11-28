function [ kappa ] = exthmin( wav,temp,eldens )

k = 1.380658e-16;           % Boltzmann's constant [erg K]
h = 6.62607e-27;            % Planck's constant [erg s]
c = 2.99792e10;             % Speed of light [cm/s]
theta = 5040./temp;
elpress = eldens.*k.*temp;

sigmabf = (1.99654 - 1.18267e-5.*wav + 2.64243e-6.*wav.^2 ...
            - 4.40524e-10.*wav.^3 +3.23992e-14.*wav.^4 ...
            - 1.39568e-18.*wav.^5 + 2.78701e-23.*wav.^6);
sigmabf = sigmabf.*1e-18;       % cm^2 per H-min ion


if sum(wav) >= 16444       % H-min ionization limit 
    sigmabf(find(wav >= 16444)) = 0; %#ok<FNDSB>
end

graysaha = 4.158e-10 .* elpress .* theta.^(2.5) .* 10.^(0.754.*theta);

% disp([ num2str(size(sigmabf)), ',', num2str(size(graysaha)) ])

kappabf = sigmabf .* graysaha;
kappabf = kappabf .* (1 - exp(-h.*c./(wav.*1e-8.*k.*temp)));

% logratio = -0.1761 - log10(elpress) + log10(2) + 2.5 .* log10(temp) - theta .* 0.754;
% disp([ 'Hmin/H ratio = ', num2str(1/(10.^logratio)) ])

lwav = log10(wav);
f0 =  -2.2763 - 1.6850 .* lwav + 0.76661 .* lwav.^2 - 0.0533464 .* lwav.^3;
f1 =  15.2827 - 9.2846 .* lwav + 1.99381 .* lwav.^2 - 0.142631 .* lwav.^3;
f2 = - 197.789 + 190.266 .* lwav - 67.9775 .* lwav.^2 + 10.6913 .* lwav.^3 - 0.625151 .* lwav.^4;
ltheta = log10(theta);
kappaff = 1e-26 .* elpress .* 10.^(f0 + f1.*ltheta + f2.*ltheta.^2);

kappa = kappabf+kappaff;
end

