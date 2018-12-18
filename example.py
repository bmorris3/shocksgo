import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from scipy.signal import periodogram
from astropy.constants import M_sun, L_sun

from shocksgo import generate_solar_fluxes, generate_stellar_fluxes

# Generate 1e7 fluxes at 60 second cadence: 
y, kernel = generate_solar_fluxes(size=5e7, cadence=60*u.s)
x = np.arange(len(y)) / 60 / 24

fig, ax = plt.subplots(3, 2, figsize=(14, 8))
ax[0, 0].plot(x, y)
ax[0, 0].set(xlabel='Time [day]', ylabel='Flux', title='Sun')

ftest, Ptest = periodogram(y, fs=1/60)

ax[0, 1].loglog(ftest * 1e6, Ptest, ',', label='Sample power')
ax[0, 1].loglog(ftest * 1e6, 2*np.pi*kernel.get_psd(2*np.pi*ftest), alpha=0.7, label='Kernel')
ax[0, 1].set_ylim([1e-10, 1e0])
ax[0, 1].set_xlim([1e-2, 5e4])
ax[0, 1].set(xlabel='Frequency [$\mu$Hz]', ylabel='Power', title='Sun')
ax[0, 1].legend()


# scale for Kepler-62
Mstar = 0.69 * M_sun
Teffstar = 4925 * u.K
Lstar = 0.21 * L_sun

y, kernel = generate_stellar_fluxes(5e7, Mstar, Teffstar, Lstar, cadence=60*u.s)
x = np.arange(len(y)) / 60 / 24

ax[1, 0].plot(x, y)
ax[1, 0].set(xlabel='Time [day]', ylabel='Flux', title='Kepler-62')

ftest, Ptest = periodogram(y, fs=1/60)

ax[1, 1].loglog(ftest * 1e6, Ptest, ',', label='Sample power')
ax[1, 1].loglog(ftest * 1e6, 2*np.pi*kernel.get_psd(2*np.pi*ftest), alpha=0.7, label='Kernel')
ax[1, 1].set_ylim([1e-10, 1e0])
ax[1, 1].set_xlim([1e-2, 5e4])
ax[1, 1].set(xlabel='Frequency [$\mu$Hz]', ylabel='Power', title='Kepler-62')
ax[1, 1].legend()


# scale for Kepler-296
Mstar = 0.498 * M_sun
Teffstar = 3740 * u.K
Lstar = 0.031 * L_sun # Not in barclay paper

y, kernel = generate_stellar_fluxes(1e7, Mstar, Teffstar, Lstar, cadence=1*u.s)
x = np.arange(len(y)) / 60 / 60 / 24

ax[2, 0].plot(x, y)
ax[2, 0].set(xlabel='Time [day]', ylabel='Flux', title='Kepler-296')

ftest, Ptest = periodogram(y, fs=1)

ax[2, 1].loglog(ftest * 1e6, Ptest, ',', label='Sample power')
ax[2, 1].loglog(ftest * 1e6, 2*np.pi*kernel.get_psd(2*np.pi*ftest), alpha=0.7, label='Kernel')
ax[2, 1].set_ylim([1e-10, 1e0])
ax[2, 1].set_xlim([1e-2, 5e4])
ax[2, 1].set(xlabel='Frequency [$\mu$Hz]', ylabel='Power', title='Kepler-296')
ax[2, 1].legend()

fig.tight_layout()
plt.show()