%% SSB

%% FALC
load('falc.mat')

h = h(:);               % Height above tau500=1     [km] = [1e5 cm]
tau500 = tau5(:);       % Optical depth at 500nm    [1]
cmass = colm(:);        % Mass of column            [g/cm^2]
temp = temp(:);         % Temperature               [K]
v_t = vturb(:);         % Microturbulent velocity   [km/s] = [1e5 cm/s]
n_H = nhyd(:);          % Hydrogen density          [1/cm^3]
n_p = nprot(:);         % Proton density            [1/cm^3]
n_e = nel(:);           % Electron density          [1/cm^3]
Ptot = ptot(:);         % Total pressure            [dyn/cm^2] = [g/cm/s^2]
Pratio = pgasptot(:);   % Pgas/Ptot ratio           [1]
rho = dens(:);          % Density                   [g/cm^3]

m_H = 1.67352e-24;      % Hydrogen mass             [g]
m_He = 3.97*m_H;        % Helium mass               [g]
rho_H = m_H.*n_H;       % Total H mass density      [g/cm^3]
Pgas = Ptot.*Pratio;    % Gas pressure              [dyn/cm^2] = [g/cm/s^2]
n_He = 0.1.*n_H;        % Helium density            [1/cm^3]
n_nonionz = n_e-n_H;

m_e = 9.10939e-28;      % [g]
% rho_e = n_e.*m_e;       % [g/cm^3]
% rho_nonionz = (n_e-n_H).*m_e;

rho_He = m_He.*n_He;
k = 1.38065e-16;        % Boltzmann's constant [erg/K] = [cm^2*g/s^2/K]
rho_ratH = rho_H./rho;  % Mass density ratio
rho_ratHHe = (rho_H + rho_He)./rho;
kT = k.*temp;
nHnekT = (n_H + n_e).*kT;               % [g/cm/s^2]
nHnHenekT = (n_H + n_e + n_He).*kT;     % [g/cm/s^2]



%% plot FALC
% temperature against height
figure
plot(h,temp)
ylim([3000 10000])
xlabel('h [km]','Fontsize',16)
ylabel('T [K]','Fontsize',16)

% total pressure against column mass
figure
plot(cmass,Ptot)
xlabel('m [g/cm^2]','Fontsize',16)
ylabel('P_{tot} [dyn/cm^2]','Fontsize',16)

figure
loglog(cmass,Ptot)
xlabel('m [g/cm^2]','Fontsize',16)
ylabel('P_{tot} [dyn/cm^2]','Fontsize',16)

% ratio of the hydrogen mass density to total mass density against height
figure
plot(h,rho_ratH)
xlabel('h [km]','Fontsize',16)
ylabel('\rho_H / \rho','Fontsize',16)
ylim([0.71 0.72])

figure
plot(h,rho_ratHHe)
xlabel('h [km]','Fontsize',16)
ylabel('\rho_H + \rho_{He} / \rho','Fontsize',16)
ylim([0.99 1])

% column mass against height
figure
plot(h,cmass)
xlabel('h [km]','Fontsize',16)
ylabel('m_{\sigma} [g/cm^2]','Fontsize',16)

figure
semilogy(h,cmass)
xlabel('h [km]','Fontsize',16)
ylabel('m_{\sigma} [g/cm^2]','Fontsize',16)

% density against height
figure
plot(h,rho)
xlabel('h [km]','Fontsize',16)
ylabel('\rho [g/cm^3]','Fontsize',16)

%scale height
ix = find(h == 0);
iy = find(h == 2017);
rho0 = rho(ix);
rhoN = rho(iy);
Hp = 2017./(log(rho0./rhoN));
rho_a = rho0.*exp(-h./Hp);

figure
semilogy(h,rho,'Displayname','\rho')
hold on
semilogy(h,rho_a,'-.','Displayname','\rho_{approx}')
xlabel('h [km]','Fontsize',16)
ylabel('\rho [g/cm^3]','Fontsize',16)
lgd = legend('show');
lgd.FontSize = 14;

%% gas pressure against height
figure
hold on
plot(h,Pgas,'Displayname','P_{gas}')
plot(h,nHnekT,'Displayname','(n_H+n_e)kT')
plot(h,nHnHenekT,'Displayname','(n_H+n_e+n_{He})kT')
xlabel('h [km]','Fontsize',16)
ylabel('Pressure [g/cm/s^2]','Fontsize',16)
lgd = legend('show');
lgd.FontSize = 14;

figure
hold on
plot(h,Pgas./nHnekT,'Displayname','P_{gas} / (n_H+n_e)kT')
plot(h,Pgas./nHnHenekT,'Displayname','P_{gas} / (n_H+n_{He}+n_e)kT')
xlabel('h [km]','Fontsize',16)
ylabel('Pressure ratios','Fontsize',16)
lgd = legend('show');
lgd.FontSize = 14;

%% Height vs densities
figure
hold on
plot(h,n_H,'Displayname','n_H')
plot(h,n_e,'Displayname','n_e')
plot(h,n_p,'Displayname','n_p')
plot(h,n_e-n_p,'Displayname','n_{e,nonionz}')
xlabel('h [km]','Fontsize',16)
ylabel('Number density [1/cm^3]','Fontsize',16)
% ylim([0 1e14])  % limit y scale to see details
    % at larger heights, H is ionized fully (high T)
lgd = legend('show');
lgd.FontSize = 14;

figure
semilogy(h,n_H,'Displayname','n_H')
hold on
semilogy(h,n_e,'Displayname','n_e')
semilogy(h,n_p,'Displayname','n_p')
semilogy(h,n_e-n_p,'Displayname','n_{e,nonionz}')
xlabel('h [km]','Fontsize',16)
ylabel('Number density [1/cm^3]','Fontsize',16)
% ylim([0 1e14])  % limit y scale to see details
    % at larger heights, H is ionized fully (high T)
lgd = legend('show');
lgd.FontSize = 14;


% ionization fraction of H vs height
figure
semilogy(h,n_p./n_H)
xlabel('h [km]','Fontsize',16)
ylabel('Ionization fracion of Hydrogen','Fontsize',16)

% photon density
hmin = find(h == min(h));
hmax = find(h == max(h));
t1 = temp(hmin);
t2 = temp(hmax);
teff = 5770;        % [K] effectiv temp

nphot1 = 20*t1^3;           % photon dens at lowest location
nphot2 = 20*teff^3/(2*pi);  % photon dens at highest location
disp(['Photon density: ',newline,'h_low = ',num2str(nphot1),newline,'h_high = ',num2str(nphot2),newline])
nH1 = n_H(hmin);
nH2 = n_H(hmax);
disp(['Hydrogen density: ',newline,'h_low = ',num2str(nH1),newline,'h_high = ',num2str(nH2),newline])



%% Earth data
load('earth.mat')

hE = h_E(:);
P_E = 10.^logP_E(:);
temp_E = T_E(:);
rho_E = 10.^logrho_E(:);
n_E = 10.^logN_E(:);

%% Earth plots
%temp
figure
plot(hE,temp_E)
xlabel('h [km]','Fontsize',16)
ylabel('T [K]','Fontsize',16)

%lin plots
% pressure 
figure
plot(hE,P_E)
xlabel('h [km]','Fontsize',16)
ylabel('P [g/cm/s^2]','Fontsize',16)

% particle density
figure
plot(hE,n_E)
xlabel('h [km]','Fontsize',16)
ylabel('Particle density [g/cm^3]','Fontsize',16)

% gas density
figure
plot(hE,rho_E)
xlabel('h [km]','Fontsize',16)
ylabel('Gas density [g/cm^3]','Fontsize',16)

% log plots - pressure, particle density, gas density
figure
semilogy(hE,P_E)
xlabel('h [km]','Fontsize',16)
ylabel('P [g/cm/s^2]','Fontsize',16)

figure
semilogy(hE,n_E)
xlabel('h [km]','Fontsize',16)
ylabel('Particle density [1/cm^3]','Fontsize',16)

figure
semilogy(hE,rho_E)
xlabel('h [km]','Fontsize',16)
ylabel('Gas density [g/cm^3]','Fontsize',16)

%% pressure + density
Pn = P_E./max(P_E);     % normalized
rhon = rho_E./max(rho_E);

figure
hold on
plot(hE,Pn,'Displayname','P/P_{max}')
plot(hE,rhon,'Displayname','\rho/\rho_{max}')
xlabel('h [km]','Fontsize',16)
ylabel('Normalized units','Fontsize',16)
lgd = legend('show');
lgd.FontSize = 14;

figure
semilogy(hE,Pn,'Displayname','P/P_{max}')
hold on
semilogy(hE,rhon,'Displayname','\rho/\rho_{max}')
xlabel('h [km]','Fontsize',16)
ylabel('Normalized units','Fontsize',16)
lgd = legend('show');
lgd.FontSize = 14;

% mean molecular weight vs height
mu = rho_E./(n_E.*m_H);
figure
plot(hE,mu)
xlabel('h [km]','Fontsize',16)
ylabel('\mu','Fontsize',20)

%% scale height


ii = find(hE == 0);
ij = find(hE == 120);
rho0_E = rho_E(ii);
rhoN_E = rho_E(ij);
HpE = 120./(log(rho0_E./rhoN_E));
rho_aE = rho0_E.*exp(-hE./HpE);

figure
semilogy(hE,rho_E,'Displayname','\rho')
hold on
semilogy(hE,rho_aE,'-.','Displayname','\rho_{approx}')
xlabel('h [km]','Fontsize',16)
ylabel('\rho [g/cm^3]','Fontsize',16)
lgd = legend('show');
lgd.FontSize = 14;

%% Chapter 2

load('solarspecs.mat')

F = flux;   % [1e10 erg/cm^2/s/microns/ster]
Fc = flux_cont;
I = irrad;  % [1e10 erg/cm^2/s/microns/ster]
Ic = irrad_cont;
w = lambda; % [microns]

maxIc = max(Ic);
maxat = find(Ic == max(Ic));
disp(['Ic max = ', num2str(maxIc),'e10 at w = ',num2str(w(maxat))])
    
figure
hold on
plot(w,F,'Displayname','F_{\lambda}')
plot(w,Fc,'Displayname','F_{\lambda}^c')
plot(w,I,'Displayname','I_{\lambda}')
plot(w,Ic,'Displayname','I_{\lambda}^c')
xlim([0 2])
xlabel('\lambda [\mum]','Fontsize',16)
ylabel('[10^{10} erg/cm^2/s/\mum/ster]','Fontsize',14)
lgd = legend('show');
lgd.FontSize = 14;

c = 2.99792e10;   % cm/s
conv = (1/c).*w.^2.*1e6; % dw = v^2/c dv * 1e6 microns..

figure
hold on
plot(w,F.*conv,'Displayname','F_{\nu}')
plot(w,Fc.*conv,'Displayname','F_{\nu}^c')
plot(w,I.*conv,'Displayname','I_{\nu}')
plot(w,Ic.*conv,'Displayname','I_{\nu}^c')
xlim([0 2])
xlabel('\lambda [\mum]','Fontsize',16)
ylabel('Irradiance [10^{10} erg/cm^2/s/Hz/ster]','Fontsize',14)
lgd = legend('show');
lgd.FontSize = 14;

%% planck fit

wav= w./1e4;   % [cm]
Bwav = planck(6460,wav)./1e14; % converting back to 10^10 erg/.../microns

% planck curve for T = 6460
figure
hold on
plot(w,Bwav,'Displayname','B_{\lambda}')        % wavelength
plot(w,Ic,'Displayname','I_{\lambda}^c')
xlim([0 2])
title('Plack function at T = 6460 K','Fontsize',14)
xlabel('\lambda [\mum]','Fontsize',16)
ylabel('Irradiance [10^{10} erg/cm^2/s/Hz/ster]','Fontsize',14)
lgd = legend('show');
lgd.FontSize = 14;

% loop - planck for all
% temperature = linspace(5900,6500,4);
% Btest = zeros([length(wav) 11]);
% 
% figure
% hold on
% for i = 1:length(temperature)
%     Btest(:,i) = planck(temperature(i),wav)./1e14;
%     plot(w,Btest(:,i))
% end
% legend show


%% Brightness temperature

Tb = brighttemp(Ic.*1e4,wav);   %both in cgs

figure
plot(w,Tb)
xlabel('\lambda [\mum]','Fontsize',16)
ylabel('Brightness temperature T_B [K]','Fontsize',14)

%% Continuous extinction
load('falc.mat')

wavA = w(1:38).*10000;      % [?]
temp1 = temp(1:38);         % [K]
eldens1 = n_e(1:38);        % [1/cm^3]

h0 = find(h == 0);
temp0 = temp(h0);
eldens0 = n_e(h0);

Hmin = exthmin(wavA,temp0,eldens0).*1e24;

figure
plot(w,Hmin)
xlim([0 2])
title('h = 0 km, T = 6520 K, n_e = 7.6970e13 /cm^3')
xlabel('\lambda [\mum]','Fontsize',16)
ylabel('H_{bf+ff}^- extinction [10^{-26} cm^2/H-atom/dyn/cm^2]','Fontsize',14)

% Invert to make look like Tb
figure
plot(w,1./Hmin)
xlim([0 5])
title('h = 0 km, T = 6520 K, n_e = 7.6970e13 /cm^3')
xlabel('\lambda [\mum]','Fontsize',16)
ylabel('Inverted H_{bf+ff}^- extinction [10^{-26} cm^2/H-atom/dyn/cm^2]^{-1}','Fontsize',12)


%% extintion vs height

nHN = n_H - n_p;         % [1/cm^3]
idx1 = find(w == 0.5);
w05 = w(idx1).*10000;    % [?]
temp05 = temp(idx1);
eldens05 = n_e(idx1);
sigmaT = 6.648e-25;     % [cm^2]

% alpha(H-)
aHmin = nHN.*exthmin(w05,temp05,eldens05.*sigmaT).*1e26;

% alpha(e)
aHminT = n_e.*sigmaT;

% add together
aHminC = aHmin + aHminT;

figure
semilogy(h,aHmin,'Displayname','\alpha(H^-)')
hold on
semilogy(h,aHminT,'Displayname','\alpha(e)')
semilogy(h,aHminC,'Displayname','\alpha(H^-) + \alpha(e)')
xlabel('h [km]','Fontsize',16)
ylabel('Extinction [10^{-26} /cm]','Fontsize',14)
lgd = legend('show');
lgd.FontSize = 14;

%% Optical depth
% wavelength in ?

tau = zeros([1 length(tau500)]);

for i = 2:length(tau500)
    tau(i) = tau(i-1) + 0.5.*((exthmin(5000,temp(i),n_e(i)) + n_e(i).*sigmaT) + ...
        (exthmin(5000,temp(i-1),n_e(i-1)) + n_e(i).*sigmaT)) .* (h(i-1)-h(i)).*1e5;
end

figure
semilogy(h,tau,'Displayname','tau_{500} Numerical')
hold on
semilogy(h,tau500,'Displayname','tau_{500} FALC')
xlabel('h [km]','Fontsize',16)
ylabel('\tau_{\lambda}','Fontsize',20)
lgd = legend('show');
lgd.FontSize = 14;


%% Emergent intensity and height of formation

sigmaT = 6.648e-25;     % [cm^2]

ext = zeros([1 length(tau500)]);
integrand = zeros([1 length(tau500)]);
contfunc = zeros([1 length(tau500)]);
intt = 0;
hint = 0;
wl = [0.5 1 1.6 5];
cf = zeros([4 length(tau500)]);
cfn = zeros([4 length(tau500)]);

for j = 1:4
    for i = 2:length(tau500)-1
        ext(i) = exthmin(wl(j)*1e4,temp(i), n_e(i)) .* (n_H(i) - n_p(i) + sigmaT.*n_e(i));
        tau(i) = tau(i-1) + 0.5 .* (exthmin(wl(j)*1e4,temp(i),n_e(i)) + ...
            exthmin(wl(j)*1e4,temp(i-1),n_e(i-1))) .* (h(i-1)-h(i)).*1e23;
%         tau(i) = tau(i-1) + 0.5.*((exthmin(5000,temp(i),n_e(i)) + n_e(i).*sigmaT) + ...
%            (exthmin(5000,temp(i-1),n_e(i-1)) + n_e(i).*sigmaT)) .* (h(i-1)-h(i)).*1e5;
        integrand(i) = planck(temp(i),wl(j)*1e-4).*exp(-tau(i));
        intt = intt + 0.5 * (integrand(i)+integrand(i-1)) * (tau(i)-tau(i-1));
        hint = hint + h(i) * 0.5 * (integrand(i) + integrand(i-1)) * (tau(i) - tau(i-1));
        contfunc(i) = integrand(i) .* ext(i);
        cf(j,i) = contfunc(i);
    end
hmean = hint/intt;
disp(intt)
% disp(Ic(idx1)*1e14)
disp(hmean)
cfn(j,:) = cf(j,:)./max(cf(j,:));   % peak-normalized
hm = [1.0106e+03 960.1713 909.7930 910.6937];

end

figure
hold on
plot(h,cfn(1,:),'Color',[50/255 205/255 50/255],'Displayname', '\lambda = 0.5 \mum, <h> = 1010.6')
plot(h,cfn(2,:),'Color',[255/255 140/255 0/255],'Displayname', '\lambda = 1.0 \mum, <h> = 960.17')
plot(h,cfn(3,:),'Color',[199/255 21/255 133/255],'Displayname', '\lambda = 1.6 \mum, <h> = 909.79')
plot(h,cfn(4,:),'Color',[70/255 130/255 180/255],'Displayname', '\lambda = 5.0 \mum, <h> = 910.69')
line([hm(1) hm(1)],[0 1],'Color',[50/255 205/255 50/255])
line([hm(2) hm(2)],[0 1],'Color',[255/255 140/255 0/255])
line([hm(3) hm(3)],[0 1],'Color',[199/255 21/255 133/255])
line([hm(4) hm(4)],[0 1],'Color',[70/255 130/255 180/255])

ylim([0 1])
xlabel('h [km]','Fontsize',16)
ylabel('Contribution function')
lgd = legend('show','\lambda = 0.5 \mum, <h> = 1010.6','\lambda = 1.0 \mum, <h> = 960.17', ...
    '\lambda = 1.6 \mum, <h> = 909.79','\lambda = 5.0 \mum, <h> = 910.69');
lgd.FontSize = 10;


%% Disk-center intensity
% Emergent continuum
load('solarspecs.mat')

ext = zeros([1 length(tau500)]);
integrand = zeros([1 length(tau500)]);
contfunc = zeros([1 length(tau500)]);
intt = 0;
hint = 0;
cf = zeros([length(w) length(tau500)]);
calcwav = zeros([1 length(w)]);

for j = 1:length(w)
    for i = 2:length(tau500)-1
        ext(i) = exthmin(w(j)*1e4,temp(i), (n_e(i) * (n_H(i) - n_p(i)) + sigmaT.*n_e(i)));
        tau(i) = tau(i-1) + 0.5 .* (exthmin(w(j)*1e4,temp(i),n_e(i)) + ...
            exthmin(w(j)*1e4,temp(i-1),n_e(i-1))) .* (h(i-1)-h(i)).*1e19;
        integrand(i) = planck(temp(i),w(j)*1e-4).*exp(-tau(i));
        intt = intt + 0.5 * (integrand(i)+integrand(i-1)) * (tau(i)-tau(i-1));
        hint = hint + h(i) * 0.5 * (integrand(i) + integrand(i-1)) * (tau(i) - tau(i-1));
        contfunc(i) = integrand(i) .* ext(i);
        cf(j,i) = contfunc(i);
    end
    calcwav(j) = intt/1e14;
end

figure
hold on
plot(w,Ic,'Displayname','FALC')
%plot(w,calcwav,'Displayname','Numerical')
plot(w,cf(:,41),'Displayname','Numerical')     % looks best
ylabel('I_C [10^{10} erg/cm^2/s/cm/ster]','Fontsize',14)
xlabel('\lambda [\mum]','Fontsize',16)
lgd = legend('show');
lgd.FontSize = 14;

%% Limb darknening

mudrk = 0.1:0.1:1.0;
intmu = zeros([length(mudrk) length(w)]);
intspec = zeros([1 length(w)]);
intfinal = zeros([length(mudrk) length(w)]);

for imu = 1:length(mudrk)
    mu = mudrk(imu);
    for iw = 1:length(w)
        wl = w(iw);
        ext = zeros([1 length(tau500)]);
        tau = zeros([1 length(tau500)]);
        integrand = zeros([1 length(tau500)]);
        intt = 0;
        for i = 2:length(w)
            ext(i) = exthmin(wl*1e4,temp(i), (n_e(i) * (n_H(i) - n_p(i)) + sigmaT.*n_e(i)));
            tau(i) = tau(i-1) + 0.5 .* (exthmin(wl*1e4,temp(i),n_e(i)) + ...
                exthmin(wl*1e4,temp(i-1),n_e(i-1))) .* (h(i-1)-h(i)).*1e5;
            integrand(i) = planck(temp(i),wl*1e-4).*exp(-tau(i)./mu);
            intt = intt + 0.5 * (integrand(i)+integrand(i-1)) * (tau(i)-tau(i-1))./mu;
        end
        intmu(imu,iw) = intt;
        intspec(iw) = intspec(iw) + intmu(imu,iw).*mu;
        intfinal(imu,iw) = intspec(iw);
    end
end

figure
hold on
plot(w,intfinal(1,:),'Displayname','\mu = 0.1')
plot(w,intfinal(4,:),'Displayname','\mu = 0.4')
plot(w,intfinal(7,:),'Displayname','\mu = 0.7')
plot(w,intfinal(10,:),'Displayname','\mu = 1.0')
ylabel('I_C(0,\mu) [10^{10} erg/cm^2/s/cm/ster]','Fontsize',14)
xlabel('\lambda [\mum]','Fontsize',16)
lgd = legend('show');
lgd.FontSize = 14;

%% Flux integration

%integrate intensity in 2.6 over mu in order to get flux.
xgauss = [-0.7745966692 0.0000000000 0.7745966692];
wgauss = [ 0.5555555555 0.8888888888 0.5555555555];
fluxspec = zeros([1 length(w)]);
intmu = zeros([3 length(w)]);


for imu = 1:3
	mu = 0.5 + xgauss(imu)./2;       % rescale xrange [-1,+1] to [0,1]
    wg = wgauss(imu)./2;             % weights add up to 2 on [-1,+1]
    for iw = 1:length(w)
        wl = w(iw);
        ext = zeros([1 length(tau500)]);
        tau = zeros([1 length(tau500)]);
        integrand = zeros([1 length(tau500)]);
        ints = 0;
        for i = 2:length(w)
            ext(i) = exthmin(wl*1e4,temp(i), (n_e(i) * (n_H(i) - n_p(i)) + sigmaT.*n_e(i)));
            tau(i) = tau(i-1) + 0.5 .* (exthmin(wl*1e4,temp(i),n_e(i)) + ...
                exthmin(wl*1e4,temp(i-1),n_e(i-1))) .* (h(i-1)-h(i)).*1e5;
            integrand(i) = planck(temp(i),wl*1e-4).*exp(-tau(i)./mu);
            ints = ints + 0.5 * (integrand(i)+integrand(i-1)) * (tau(i)-tau(i-1))./mu;
        end
        intmu(imu,iw) = ints;
        fluxspec(iw) = fluxspec(iw) + wg.*intmu(imu,iw).*mu;
    end
end

fluxspec = 2.*fluxspec;

% flux
figure
hold on
plot(w,fluxspec.*4e4,'Displayname','Numerical')
plot(w,Fc,'Displayname','FALC')
ylabel('F_C [10^{10} erg/cm^2/s/cm/ster]','Fontsize',14)
xlabel('\lambda [\mum]','Fontsize',16)
lgd = legend('show');
lgd.FontSize = 14;



%% Chapter 3
load('nadata.mat')
load('falc.mat')

figure
plot(wavnr,Sspec)
xlabel('Wavenumber [1/cm]','Fontsize',14)
title('Observed solar spetrum','Fontsize',14)

wav = 1./wavnr.*1e8;

figure
plot(wav,Sspec)
xlabel('Wavelength [?]','Fontsize',14)
title('Observed solar spetrum','Fontsize',14)

% Split into two regimes

Sspec1 = Sspec(1:3600);
Sspec2 = Sspec(3601:end);

idx1 = find(Sspec1 == min(Sspec1));
idx2 = find(Sspec2 == min(Sspec2));

NaID1minv = wav(idx1);
NaID2minv = wav(3600+idx2);
disp(NaID1minv)
disp(NaID2minv)      % [?]
disp(newline)

% Correct for air
wair = 0.99972683.*wav + 0.0107 - 196.25./wav;

NaID1min = wair(idx1);
NaID2min = wair(3600+idx2);
disp(NaID1min)
disp(NaID2min)      % [?]

figure
plot(wair,Sspec)
title('Observed Solar spectrum in air','Fontsize',14)
xlabel('Wavelength [?]','Fontsize',14)

%

b_l = 1;
b_u = 1;

A_Na = 1.8*1e-6;
f_lu = [0.318 0.631];
E_ioniz = [5.139 47.29 71.64];

boltz = zeros([3 length(temp)]);

for i = 1:length(temp)
	boltz(1,i) = Boltz_Na(temp(i), 0, 0);
	boltz(2,i) = Boltz_Na(temp(i), 0, 1);
	boltz(3,i) = Boltz_Na(temp(i), 0, 2);
end

% figure
% plot(h,boltz)

saha = zeros([2 length(temp)]);

for i = 1:length(temp)
	saha(1,i) = saha_Na(temp(i),nel(i),1,E_ioniz);
	saha(2,i) = saha_Na(temp(i),nel(i),2,E_ioniz);
end

% figure
% plot(h,saha(1,:))
% figure
% plot(h,saha(2,:))

sahaboltz = zeros([3 length(temp)]);
for i = 1:length(temp)
	sahaboltz(1,i) = sahabolt_Na(temp(i), nel(i), 0, 0, E_ioniz);
	sahaboltz(2,i) = sahabolt_Na(temp(i), nel(i), 0, 1, E_ioniz);
	sahaboltz(3,i) = sahabolt_Na(temp(i), nel(i), 0, 2, E_ioniz);
end
m_Na = 22.99*1.6605e-24;        % [g]
kerg = 1.380658e-16;

doppler_term1 = sqrt(2.*kerg.*temp./m_Na).*1e-5;
doppler_term2 = sqrt(vturb.*vturb);

figure
hold on
plot(h,doppler_term1.*1e8)
plot(h,doppler_term2.*1e8)

doppler = zeros([3 length(temp)]);
wav = [NaID1min NaID2min]*1e-8;                         %wavelengths of line in cm
doppler(1,:) = dopplerwidth(wav(1),temp,vturb,m_Na);    %values for NaID1
doppler(2,:) = dopplerwidth(wav(2),temp,vturb,m_Na);    %for NaID1


vdwbroadening = zeros([1 length(temp)]);

for i = 1:length(temp)
	vdwbroadening(i) =  gammavdw_NaD(temp(i), Pgas(i), 2, wav, E_ioniz);
end

figure
plot(h,vdwbroadening)

% works perfect so far
%% voigt
c = 3e10; 
wavelen = linspace(wav(1) - 2e-8, wav(1) + 2e-8, 1001);

a = zeros([length(h) length(wavelen)]);
v = zeros([length(h) length(wavelen)]);
voigtprofile = zeros([length(h) length(wavelen)]);

gamma = gammavdw_NaD(temp,Pgas,2,wav, E_ioniz);	% van der waal damping

for j = 1:length(h)
    for i = 1:length(wavelen)
		a(j,i) = wavelen(i).^2./(4 .* pi .*c) .*gamma(j)./doppler(1,j);
		v(j,i) = (wavelen(i)-wav(1))./doppler(1,j);
		voigtprofile(j,i) = voigt(a(j,i),v(j,i)).*doppler(1,j).*sqrt(pi);
    end
end

%
deltawav = (wavelen - wav(1)).*1e8;	%delta wav in angstrom


h1 = find(h == 0);
h2 = find(h == 200);
h3 = find(h == 400);

figure
hold on
plot(deltawav,real(voigtprofile(h1,1)))
plot(deltawav,real(voigtprofile(h2,1)))
plot(deltawav,real(voigtprofile(h3,1)))

% voigt profile does not work :(



