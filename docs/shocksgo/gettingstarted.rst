Getting Started
===============

This tutorial will show some examples of how to make solar and stellar light
curves using ``shocksgo``.

Generating a Solar Light Curve
------------------------------

To generate a sample of 1e5 solar fluxes at 60 second cadence, we can use
`~shocksgo.generate_solar_fluxes`::

    import matplotlib.pyplot as plt
    import astropy.units as u

    from shocksgo import generate_solar_fluxes

    fluxes, kernel = generate_solar_fluxes(size=1e5, cadence=60*u.s)

    times = np.arange(len(fluxes))

    plt.plot(times, 1e6 * fluxes)
    plt.gca().set(xlabel='Time [minutes]', ylabel='Relative Flux [ppm]')
    plt.show()

.. plot::

    import matplotlib.pyplot as plt
    import astropy.units as u

    from shocksgo import generate_solar_fluxes

    fluxes, kernel = generate_solar_fluxes(size=1e5, cadence=60*u.s)

    times = np.arange(len(fluxes))

    plt.plot(times, 1e6 * fluxes)
    plt.gca().set(xlabel='Time [minutes]', ylabel='Relative Flux [ppm]')
    plt.show()

We can check that the power spectrum of the fluxes that we've generated
reproduce the solar power spectrum::

    import matplotlib.pyplot as plt
    import numpy as np
    import astropy.units as u
    from scipy.signal import periodogram

    from shocksgo import generate_solar_fluxes

    fluxes, kernel = generate_solar_fluxes(size=1e7, cadence=60*u.s)

    freq, power = periodogram(fluxes, fs=1/60)

    plt.loglog(freq * 1e6, power, ',', label='Samples')
    plt.loglog(freq * 1e6, 2*np.pi*kernel.get_psd(2*np.pi*freq), alpha=0.7, label='Kernel')
    plt.ylim([1e-10, 1e0])
    plt.xlim([1e-2, 1e4])
    plt.gca().set(xlabel='Frequency [$\mu$Hz]', ylabel='Power')
    plt.show()


.. plot::

    import matplotlib.pyplot as plt
    import numpy as np
    import astropy.units as u
    from scipy.signal import periodogram

    from shocksgo import generate_solar_fluxes

    fluxes, kernel = generate_solar_fluxes(size=1e7, cadence=60*u.s)

    freq, power = periodogram(fluxes, fs=1/60)

    plt.loglog(freq * 1e6, power, ',', label='Samples')
    plt.loglog(freq * 1e6, 2*np.pi*kernel.get_psd(2*np.pi*freq), alpha=0.7, label='Kernel')
    plt.ylim([1e-10, 1e0])
    plt.xlim([1e-2, 1e4])
    plt.gca().set(xlabel='Frequency [$\mu$Hz]', ylabel='Power')
    plt.show()

Zooming into the p-mode oscillations, we can see the peaks are reproduced:

.. plot::

    import matplotlib.pyplot as plt
    import numpy as np
    import astropy.units as u
    from scipy.signal import periodogram

    from shocksgo import generate_solar_fluxes

    fluxes, kernel = generate_solar_fluxes(size=5e7, cadence=60*u.s)

    freq, power = periodogram(fluxes, fs=1/60)

    plt.semilogy(freq * 1e6, power, ',', label='Samples')
    plt.semilogy(freq * 1e6, 2*np.pi*kernel.get_psd(2*np.pi*freq), alpha=0.7, label='Kernel')
    plt.ylim([1e-8, 1e-4])
    plt.xlim([2000, 4000])
    plt.gca().set(xlabel='Frequency [$\mu$Hz]', ylabel='Power')
    plt.show()


Generating a Stellar Light Curve
--------------------------------


To generate a sample of 1e5 *steller* fluxes at 60 second cadence, we can use
`~shocksgo.generate_stellar_fluxes`::

    import matplotlib.pyplot as plt
    import astropy.units as u
    from astropy.constants import M_sun, L_sun

    from shocksgo import generate_stellar_fluxes

    # Stellar properties
    M = 0.9 * M_sun
    T_eff = 5340 * u.K
    L = 0.56 * L_sun

    fluxes, kernel = generate_stellar_fluxes(size=1e7, M=M, T_eff=T_eff, L=L, cadence=60*u.s)

    times = np.arange(len(fluxes)) / 60 / 24

    plt.plot(times, 1e6 * fluxes)
    plt.gca().set(xlabel='Time [days]', ylabel='Relative Flux [ppm]', title='G9V star')
    plt.show()

.. plot::

    import matplotlib.pyplot as plt
    import astropy.units as u
    from astropy.constants import M_sun, L_sun

    from shocksgo import generate_stellar_fluxes

    # Stellar properties
    M = 0.9 * M_sun
    T_eff = 5340 * u.K
    L = 0.56 * L_sun

    fluxes, kernel = generate_stellar_fluxes(size=1e7, M=M, T_eff=T_eff, L=L, cadence=60*u.s)

    times = np.arange(len(fluxes)) / 60 / 24

    plt.plot(times, 1e6 * fluxes)
    plt.gca().set(xlabel='Time [days]', ylabel='Relative Flux [ppm]', title='G9V star')
    plt.show()

We can see the shift in the p-mode oscillations relative to the solar ones above
if we plot the power spectrum::

    import matplotlib.pyplot as plt
    import numpy as np
    import astropy.units as u
    from scipy.signal import periodogram
    from astropy.constants import M_sun, L_sun

    from shocksgo import generate_stellar_fluxes

    # Stellar properties
    M = 0.9 * M_sun
    T_eff = 5340 * u.K
    L = 0.56 * L_sun

    fluxes, kernel = generate_stellar_fluxes(size=1e7, M=M, T_eff=T_eff, L=L, cadence=1*u.s)

    freq, power = periodogram(fluxes, fs=1)

    plt.semilogy(freq * 1e6, power, ',', label='Samples')
    plt.semilogy(freq * 1e6, 2*np.pi*kernel.get_psd(2*np.pi*freq), alpha=0.7, label='Kernel')
    plt.ylim([1e-8, 1e-4])
    plt.xlim([3000, 6000])
    plt.gca().set(xlabel='Frequency [$\mu$Hz]', ylabel='Power')
    plt.show()


.. plot::

    import matplotlib.pyplot as plt
    import numpy as np
    import astropy.units as u
    from scipy.signal import periodogram
    from astropy.constants import M_sun, L_sun

    from shocksgo import generate_stellar_fluxes

    # Stellar properties
    M = 0.9 * M_sun
    T_eff = 5340 * u.K
    L = 0.56 * L_sun

    fluxes, kernel = generate_stellar_fluxes(size=1e7, M=M, T_eff=T_eff, L=L, cadence=1*u.s)

    freq, power = periodogram(fluxes, fs=1)

    plt.semilogy(freq * 1e6, power, ',', label='Samples')
    plt.semilogy(freq * 1e6, 2*np.pi*kernel.get_psd(2*np.pi*freq), alpha=0.7, label='Kernel')
    plt.ylim([1e-8, 1e-4])
    plt.xlim([3000, 6000])
    plt.gca().set(xlabel='Frequency [$\mu$Hz]', ylabel='Power')
    plt.show()