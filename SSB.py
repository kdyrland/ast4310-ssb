import numpy as np
import matplotlib.pyplot as plt
from scipy import special

# Chapter 3 of SSB

falc = np.loadtxt('data/falc.dat', unpack=True)
height = falc[0]
tau5 = falc[1]
colm = falc[2]
temp = falc[3]
vturb = falc[4]
nhyd = falc[5]
nprot = falc[6]
nel = falc[7]
ptot = falc[8]
pgasptot = falc[9]
dens = falc[10]

# Read in NaD data
Nadata = np.loadtxt('data/int_nad.dat', unpack=True)
wavenumber = Nadata[0]
E_spec = Nadata[1]
S_spec = Nadata[2]
S_spec_corr = Nadata[3]

k = 1.380658e-16  # Boltzmann constant [erg/K]
kerg = 1.380658e-16
h = 6.626076e-27  # Planck constant [erg s]
c = 2.997929e10  # velocity of light [cm/s]
# other parameters
theta = 5040. / temp
elpress = nel * k * temp

h = 6.62607e-27  # Planck constant (erg s)
c = 2.99792e10  # light speed [cm/s]
erg2eV = 1 / 1.60219e-12  # erg to eV conversion
E_ionization = 5.139
keV = 8.61734e-5
keVT = keV * temp;
kergT = kerg * temp;
m_e = 9.10939e-28;
m_Na = 22.99 * 1.6605e-24;

plt.figure(0)
plt.plot(wavenumber, S_spec)
plt.title('Observed Solar spectrum')
# plt.ylabel('')
plt.xlabel(r'Wavenumber[cm$^{-1}$]')
plt.xlim([np.min(wavenumber), np.max(wavenumber)])
# plt.show()

wavelength = 1. / wavenumber * 1e8

plt.figure(1)
plt.plot(wavelength, S_spec)
plt.title('Observed Solar spectrum')
# plt.ylabel('')
plt.xlabel(r'Wavelength[$\AA$]')
plt.xlim([np.min(wavelength), np.max(wavelength)])
plt.show()


# Split spectrum into two regimes to find local minimums
S_spec1 = S_spec[:3600]
S_spec2 = S_spec[3600:]
# print wavelength of for each minima
print('First line has minimum at : ', wavelength[3600 + np.argmin(S_spec2)], 'angstrom')
print('Second line has minimum at : ', wavelength[np.argmin(S_spec1)], 'angstrom')


def vac_to_air(wavelength):
    "takes in vacuum wavelength in ansgtrom and returns air wavelength in angstrom"""
    return 0.99972683 * wavelength + 0.0107 - 196.25 / wavelength


wavelength_air = vac_to_air(wavelength)
NaI_D2_minimum = vac_to_air(wavelength[3600 + np.argmin(S_spec2)])
NaI_D1_minimum = vac_to_air(wavelength[np.argmin(S_spec1)])
print('First line after correction has minimum at : ', NaI_D1_minimum, 'angstrom')
print('Second line after corection has minimum at : ', NaI_D2_minimum, 'angstrom')

plt.figure(2)
plt.plot(wavelength_air, S_spec)
plt.title('Observed Solar spectrum in air')
# plt.ylabel('')
plt.xlabel(r'Wavelength[$\AA$]')
plt.xlim([np.min(wavelength_air), np.max(wavelength_air)])
plt.show()


# 3.6
def planck(temp,wav):
    kT = kerg * temp
    B = 2*h*c**2. / (wav**5.*(np.exp(h*c / (wav*kT)) - 1));
    return B


def exthmin(wav, temp, nel):
    theta = 5040. / temp
    elpress = nel * k * temp
    sigmabf = (1.99654 - 1.18267E-5 * wav + 2.64243E-6 * wav ** 2
               - 4.40524E-10 * wav ** 3 + 3.23992E-14 * wav ** 4
               - 1.39568E-18 * wav ** 5 + 2.78701E-23 * wav ** 6)
    sigmabf *= 1E-18  # cm^2 per H-min ion
    if np.size(wav) > 1:
        sigmabf[np.where(wav > 16444)] = 0  # H-min ionization limit at lambda=
    elif (np.size(wav) == 1):
        if wav > 16444:
            sigmabf = 0
    graysaha = 4.158E-10 * elpress * theta ** 2.5 * 10. ** (0.754 * theta)  # Gray (8.12)
    kappabf = sigmabf * graysaha  # per neutral H atom
    kappabf = kappabf * (1. - np.exp(-h * c / (wav * 1E-8 * k * temp)))
    lwav = np.log10(wav)
    f0 = -2.2763 - 1.6850 * lwav + 0.76661 * lwav ** 2 - 0.0533464 * lwav ** 3
    f1 = 15.2827 - 9.2846 * lwav + 1.99381 * lwav ** 2 - 0.142631 * lwav ** 3
    f2 = (-197.789 + 190.266 * lwav - 67.9775 * lwav ** 2 + 10.6913 * lwav ** 3
          - 0.625151 * lwav ** 4)
    ltheta = np.log10(theta)
    kappaff = 1E-26 * elpress * 10 ** (f0 + f1 * ltheta + f2 * ltheta ** 2)  # Gray (8.13)
    return kappabf + kappaff

def partfunc_Na(temp):
    u = np.zeros(3)
    theta = 5040. / temp
    c0 = 0.30955
    c1 = -0.17778
    c2 = 1.10594
    c3 = -2.42847
    c4 = 1.70721
    logU1 = (c0 + c1 * np.log10(theta) + c2 * np.log10(theta) ** 2 + c3 * np.log10(theta) ** 3
             + c4 * np.log10(theta) ** 4)
    u[0] = 10 ** logU1
    u[1] = 1.
    u[2] = 6.
    return u


b_l = 1.
b_u = 1.

A_Na = 1.8 * 1e-6
f_lu = [0.318, 0.631]
E_ionization = np.array([5.139, 47.29, 71.64])


# Boltzmann
def boltz_Na(temp, r, s):
    "Boltzmann distribution n_r,s/N_r"
    E_n1 = h * c / 5895.94e-8 * erg2eV
    E_n2 = h * c / 5889.97e-8 * erg2eV
    u = partfunc_Na(temp)
    chi = [0, E_n1, E_n2]
    g = [2, 2, 4]
    relnrs = g[s] / u[r] * np.exp(-(chi[s]) / (keV * temp))
    return relnrs


boltz = np.zeros((3, len(temp)))
for i in range(len(temp)):
    boltz[0, i] = boltz_Na(temp[i], 0, 0)
    boltz[1, i] = boltz_Na(temp[i], 0, 1)
    boltz[2, i] = boltz_Na(temp[i], 0, 2)

plt.figure(0)
plt.plot(height, boltz[0], label=r's=1 groundstate')
plt.plot(height, boltz[1], '--', label=r's=2 upper level D1')
plt.plot(height, boltz[2], '-.', label=r's=3 upper level D2')
plt.title(r'Boltzmann distribution of Na I')
plt.xlabel(r'Height [km]')
plt.ylabel(r'Population fraction $n_{1,s}/N_1$')
plt.xlim([np.min(height), 2000])
plt.legend(loc='best')
plt.show()


# Saha
def saha_Na(temp, eldens, ionstage):
    keVT = keV * temp
    kergT = kerg * temp
    u = partfunc_Na(temp)
    u = np.append(u, 2)  # append element to array
    sahaconst = (2. * np.pi * m_e * kergT / (h * h)) ** (3. / 2) * 2. / eldens
    nstage = np.zeros(4)
    nstage[0] = 1.
    for r in range(3):
        nstage[r + 1] = nstage[r] * sahaconst * u[r + 1] / u[r] * np.exp(-E_ionization[r] / keVT)
    ntotal = np.sum(nstage)
    nstagerel = nstage / ntotal
    return nstagerel[ionstage - 1]


saha = np.zeros((2, len(temp)))
for i in range(len(temp)):
    saha[0, i] = saha_Na(temp[i], nel[i], 1)
    saha[1, i] = saha_Na(temp[i], nel[i], 2)

plt.figure(1)
plt.plot(height, saha[0], '-', label='Na I')
plt.plot(height, saha[1], '--', label='Na II')
plt.xlim([np.min(height), 2000])
plt.ylim([1e-4, 10])
plt.yscale('log')
plt.title(r'Saha distribution of Na ')
plt.xlabel(r'Height [km]')
plt.ylabel(r'Ionization state fraction $N_r/N_{total}$')
plt.legend(loc='best')
plt.show()



# Sahaboltzmann
def sahabolt_Na(temp, eldens, ionstage, level):
    return saha_Na(temp, eldens, ionstage) * boltz_Na(temp, ionstage, level)


sahaboltz = np.zeros((3, len(temp)))
for i in range(len(temp)):
    sahaboltz[0, i] = sahabolt_Na(temp[i], nel[i], 0, 0)
    sahaboltz[1, i] = sahabolt_Na(temp[i], nel[i], 0, 1)
    sahaboltz[2, i] = sahabolt_Na(temp[i], nel[i], 0, 2)
plt.figure(0)
plt.plot(height, sahaboltz[0])
plt.plot(height, sahaboltz[1])
plt.plot(height, sahaboltz[2])
plt.show()



# Dopplerwidth testing

def dopplerwidth(wav, temp, v_t, m):
    '''Takes in central wavelength in cm,temperature in K, v_t in km/s,
    and m in grams and returns dopplerwidth in cm.'''
    return wav / c * np.sqrt(2. * kerg * temp / m + v_t * v_t * 1e10)


doppler_term1 = np.sqrt(2. * kerg * temp / m_Na) * 1e-5
doppler_term2 = np.sqrt(vturb * vturb)

plt.figure(2)
plt.plot(height, doppler_term1, label=r'$\sqrt{2kT/m_{Na}}$')
plt.plot(height, doppler_term2, '--', label=r'$v_{turb}$')
plt.xlim(np.min(height), 2000)
plt.ylim(0, 5)
plt.ylabel('[km/s]')
plt.xlabel('Height [km]')
plt.title('Na Doppler broadening and microturbulence')
plt.legend(loc='best')
plt.show()

############################################################################################
doppler = np.zeros((2, len(temp)))
wav = np.array([NaI_D1_minimum, NaI_D2_minimum]) * 1e-8  # wavelengths of line in cm
doppler[0] = dopplerwidth(wav[0], temp, vturb, m_Na)  # values for NaID1
doppler[1] = dopplerwidth(wav[1], temp, vturb, m_Na)  # for NaID2

plt.figure(3)
plt.plot(height, doppler[0] * 1e8, label='Na I D1')
plt.plot(height, doppler[1] * 1e8, 'r--', label='Na I D2')
plt.xlim(np.min(height), 2000)
plt.ylabel(r'Doppler width [$\AA$]')
plt.xlabel(r'Height [km]')
plt.title(r'Doppler width')
plt.legend()
plt.show()



# VanderWaal broadening NaID1 line
def gammavdw_NaD(temp, pgas, s):
    rsq_u = rsq_NaD(s)
    rsq_l = rsq_NaD(1)
    loggvdw = 6.33 + 0.4 * np.log10(rsq_u - rsq_l) + np.log10(pgas) - 0.7 * np.log10(temp)
    return 10 ** loggvdw


def rsq_NaD(s):
    E_n = np.zeros(3, dtype=float)
    E_n[1] = h * c / wav[0] * erg2eV
    E_n[2] = h * c / wav[1] * erg2eV
    Z = 1.
    Rydberg = 13.6
    l = [0., 1., 1.]
    nstar_sq = Rydberg * Z ** 2 / (E_ionization[0] - E_n[s - 1])
    rsq = nstar_sq / 2. / Z ** 2 * (5 * nstar_sq + 1 - 3 * l[s - 1] * (l[s - 1] + 1))
    return rsq


pgas = pgasptot * ptot
vdwbroadening = np.zeros(len(temp))
for i in range(len(temp)):
    vdwbroadening[i] = gammavdw_NaD(temp[i], pgas[i], 2)

plt.figure(4)
plt.plot(height, vdwbroadening)
plt.yscale('log')
plt.xlim([np.min(height), 2000])
plt.title('Van der Waal broadening')
plt.xlabel('Height [km]')
plt.ylabel('$\gamma_{vdW}$ [s$^{-1}$]')
plt.show()



# Define voigt function
def voigt(a, u):
    "Calculates the voigt function for values u and a"
    z = u + 1.j * a
    return special.wofz(z).real


# find damping parameter
wavelen = np.linspace(wav[0] - 2e-8, wav[0] + 2e-8, 1001)  # range around center

# prepare arrays
a = np.zeros((len(height), len(wavelen)))
v = np.zeros((len(height), len(wavelen)))
voigtprofile = np.zeros((len(height), len(wavelen)))
gamma = np.zeros(len(height))

# calculate the voigtprofile for the different height
gamma = gammavdw_NaD(temp, pgas, 2)  # van der waal damping
for j in range(len(height)):
    for i in range(len(wavelen)):
        a[j][i] = wavelen[i] ** 2 / (4 * np.pi * c) * gamma[j] / doppler[0][j]
        v[j][i] = (wavelen[i] - wav[0]) / doppler[0][j]
        voigtprofile[j][i] = voigt(a[j][i], v[j][i]) * doppler[0][j] * np.sqrt(np.pi)

deltawav = (wavelen - wav[0]) * 1e8  # delta wav in angstrom

plt.figure(5)
plt.plot(deltawav, voigtprofile[np.where(height == 0)][0], '-', label='0km')
plt.plot(deltawav, voigtprofile[np.where(height == 200)][0], '--', label='200km')
plt.plot(deltawav, voigtprofile[np.where(height == 400)][0], '-.', label='400km')
plt.yscale('log')
plt.xlabel('$\Delta \lambda$[$\AA$]')
plt.title('Voigt function for Na I in FALC')
plt.xlim([np.min(deltawav), np.max(deltawav)])
plt.legend()
plt.show()

correction = 1 - b_u / b_l * np.exp(-h * c / wav[0] / kerg / temp)

plt.figure(5)
plt.plot(height, correction)
plt.xlim([np.min(height), 2000])
plt.title(r'Correction term')
plt.xlabel(r'Height [km]')
plt.ylabel(r'')
plt.show()


# NaID1 extiction

NaD1_extinction = np.zeros((len(wavelen), len(height)))
ionstage = 1
level = 2
exthminD1 = np.zeros((len(wavelen), len(height)))
extcont = np.zeros((len(wavelen), len(height)))
exttotal = np.zeros((len(wavelen), len(height)))
thomson = 6.648 * 1e-25  # cm**2
extincelec = thomson * nel
for j in range(len(wavelen)):
    for i in range(len(height)):
        # Calculate NaID1 extinction
        part1 = np.sqrt(np.pi) * np.exp(2) / m_e / c * wavelen[j] * wavelen[i] / c * b_l
        part2 = sahabolt_Na(temp[i], nel[i], ionstage, level)
        part3 = nhyd[i] * A_Na * f_lu[0]
        part4 = np.sqrt(np.pi) * voigtprofile[i, j]
        part5 = 1. - b_u / b_l * np.exp(-h * c / wav[0] / kerg / temp[i])
        NaD1_extinction[j][i] = part1 * part2 * part3 * part4 * part5
        # Calculate continuum extinction
        exthminD1[j][i] = exthmin(wavelen[j] * 1e8, temp[i], nel[i]) * (nhyd[i] - nprot[i])
        extcont[j][i] = exthminD1[j][i] + extincelec[i]
        # Sum to get the total
        exttotal[j][i] = NaD1_extinction[j][i] + extcont[j][i]

plt.figure(6)
plt.plot(height, NaD1_extinction[np.where(wavelen == wav[0])][0], '-', label='NaID1')
plt.plot(height, extcont[np.where(wavelen == wav[0])][0], '--', label='Continuum')
plt.plot(height,exttotal[np.where(wavelen==wav[0])][0],'-.',label='Sum')
plt.xlim([np.min(height), 2000])
plt.title(r'Extinction at line center')
plt.xlabel(r'Height [km]')
plt.ylabel(r'Extinction $\alpha_\lambda^l$ [cm$^{-1}$]')
plt.yscale('log')
plt.legend(loc='best')
plt.ylim([1e-14, 1e-2])
plt.show()


# 3.7 Computed Na D1 line profile


NaD1_extinction_swapped = np.zeros((len(height), len(deltawav)))
extcont_swapped = np.zeros((len(height), len(deltawav)))
exttotal_swapped = np.zeros((len(height), len(deltawav)))
for i in range(len(height)):
    for j in range(len(deltawav)):
        NaD1_extinction_swapped[i][j] = NaD1_extinction[j][i]
        extcont_swapped[i][j] = extcont[j][i]
        exttotal_swapped[i][j] = exttotal[j][i]

plt.figure(7)
plt.plot(deltawav, NaD1_extinction_swapped[np.where(height == 0)][0], '-.', label='NaI h=0')
plt.plot(deltawav, extcont_swapped[np.where(height == 0)][0], '--', label='Cont h=0')
plt.plot(deltawav, exttotal_swapped[np.where(height == 0)][0], label=r'Total h=0')
plt.plot(deltawav, NaD1_extinction_swapped[np.where(height == 560)][0], '-.', label='NaI h=560')
plt.plot(deltawav, extcont_swapped[np.where(height == 560)][0], '--', label='Cont h=560')
plt.plot(deltawav, exttotal_swapped[np.where(height == 560)][0], label=r'Total h=560')
plt.yscale('log')
plt.xlim([np.min(deltawav), np.max(deltawav)])
plt.title(r'Extinction')
plt.ylabel(r'Extinction $\alpha_\lambda^l$ [cm$^{-1}$]')
plt.xlabel(r'$\Delta\lambda$[$\AA$]')
plt.legend(loc='best')
plt.show()
############################################################################################

# Find the extinction for NaID1line
tau = np.zeros(np.size(tau5))
integrand = np.zeros(np.size(tau5))
contfunc = np.zeros(np.size(tau5))

intt = np.zeros(np.size(wavelen))
hint = np.zeros(np.size(wavelen))

for j in range(np.size(wavelen)):
    for i in range(1, len(tau)):  # the index zero is not accounted for
        tau[i] = tau[i - 1] + 0.5 * (exttotal[j][i] + exttotal[j][i - 1]) * (height[i - 1] - height[i]) * 1e5
        integrand[i] = planck(temp[i], wavelen[j] * 1e4) * np.exp(-tau[i])
        intt[j] += 0.5 * (integrand[i] + integrand[i - 1]) * (tau[i] - tau[i - 1])
        hint[j] += height[i] * 0.5 * (integrand[i] + integrand[i - 1]) * (tau[i] - tau[i - 1])

wavelen *= 1e8  # convert wavelen from cm to Angstrom

plt.figure()
plt.plot(wavelen, intt / np.max(intt), '-', label='NaI line profile')
plt.plot(wavelength_air, S_spec, '--', label='Observed intensity')
plt.xlim([wav[0] * 1e8 - 2, wav[0] * 1e8 + 2])
plt.title(r'NaID1 intensity line', Fontsize=12)
plt.xlabel(r'Wavelength[$\AA$]')
plt.ylabel(r'Normalized intensity')
ax = plt.gca()  # This disables the offset
ax.ticklabel_format(useOffset=False)  # disable offset
plt.legend(loc='best')
plt.show()
