Getting Started
===============

This tutorial will show some examples of how to make solar and stellar light
curves using ``shocksgo``.

Generating a Solar Light Curve
------------------------------

To generate a sample of ten hours of solar fluxes at 60 second cadence, we can use
`~shocksgo.generate_solar_fluxes`::

    import matplotlib.pyplot as plt
    import astropy.units as u

    from shocksgo import generate_solar_fluxes

    times, fluxes, kernel = generate_solar_fluxes(duration=10*u.hour, cadence=60*u.s)

    plt.plot(times.to(u.hour), 1e6 * fluxes)
    plt.gca().set(xlabel='Time [hours]', ylabel='Relative Flux [ppm]')
    plt.show()

.. plot::

    import matplotlib.pyplot as plt
    import astropy.units as u

    from shocksgo import generate_solar_fluxes

    times, fluxes, kernel = generate_solar_fluxes(duration=10*u.hour, cadence=60*u.s)

    plt.plot(times.to(u.hour), 1e6 * fluxes)
    plt.gca().set(xlabel='Time [hours]', ylabel='Relative Flux [ppm]')
    plt.show()

We can check that the power spectrum of the fluxes that we've generated
reproduce the solar power spectrum::

    import matplotlib.pyplot as plt
    import numpy as np
    import astropy.units as u

    from shocksgo import generate_solar_fluxes, power_spectrum

    times, fluxes, kernel = generate_solar_fluxes(duration=100*u.day, cadence=60*u.s)

    freq, power = power_spectrum(fluxes, d=60)

    plt.loglog(freq * 1e6, power, ',', label='Samples')
    plt.loglog(freq * 1e6, 1e6 * kernel.get_psd(2*np.pi*freq)/(2*np.pi), alpha=0.7, label='Kernel')
    plt.ylim([1e-5, 1e3])
    plt.xlim([1e-2, 1e4])
    plt.gca().set(xlabel='Frequency [$\mu$Hz]', ylabel='Power')
    plt.show()


.. plot::

    import matplotlib.pyplot as plt
    import numpy as np
    import astropy.units as u
    from shocksgo import power_spectrum

    from shocksgo import generate_solar_fluxes

    times, fluxes, kernel = generate_solar_fluxes(duration=100*u.day, cadence=60*u.s)

    freq, power = power_spectrum(fluxes, d=60)

    plt.loglog(freq * 1e6, power, ',', label='Samples')
    plt.loglog(freq * 1e6, 1e6 * kernel.get_psd(2*np.pi*freq)/(2*np.pi), alpha=0.7, label='Kernel')
    plt.ylim([1e-5, 1e3])
    plt.xlim([1e-2, 1e4])
    plt.gca().set(xlabel='Frequency [$\mu$Hz]', ylabel='Power')
    plt.show()

Zooming into the p-mode oscillations, we can see the peaks are reproduced:

.. plot::

    import matplotlib.pyplot as plt
    import numpy as np
    import astropy.units as u
    from shocksgo import power_spectrum

    from shocksgo import generate_solar_fluxes

    times, fluxes, kernel = generate_solar_fluxes(duration=100*u.day, cadence=60*u.s)

    freq, power = power_spectrum(fluxes, d=60)

    plt.semilogy(freq * 1e6, power, ',', label='Samples')
    plt.semilogy(freq * 1e6, 1e6 * kernel.get_psd(2*np.pi*freq)/(2*np.pi), alpha=0.7, label='Kernel')
    plt.ylim([1e-5, 1e-1])
    plt.xlim([2000, 4000])
    plt.gca().set(xlabel='Frequency [$\mu$Hz]', ylabel='Power')
    plt.show()


Generating a Stellar Light Curve
--------------------------------


To generate a sample of *steller* fluxes at 60 second cadence, we can use
`~shocksgo.generate_stellar_fluxes`::

    import matplotlib.pyplot as plt
    import astropy.units as u
    from astropy.constants import M_sun, L_sun, R_sun

    from shocksgo import generate_stellar_fluxes

    # Stellar properties
    M = 0.9 * M_sun
    T_eff = 5340 * u.K
    L = 0.56 * L_sun
    R = 0.7 * R_sun

    times, fluxes, kernel = generate_stellar_fluxes(duration=100*u.day, M=M, T_eff=T_eff, R=R, L=L, cadence=60*u.s)

    plt.plot(times.to(u.day), 1e6 * fluxes)
    plt.gca().set(xlabel='Time [days]', ylabel='Relative Flux [ppm]', title='G9V star')
    plt.show()

.. plot::

    import matplotlib.pyplot as plt
    import astropy.units as u
    from astropy.constants import M_sun, L_sun, R_sun

    from shocksgo import generate_stellar_fluxes

    # Stellar properties
    M = 0.9 * M_sun
    T_eff = 5340 * u.K
    L = 0.56 * L_sun
    R = 0.7 * R_sun

    times, fluxes, kernel = generate_stellar_fluxes(duration=100*u.day, M=M, T_eff=T_eff, R=R, L=L, cadence=60*u.s)

    plt.plot(times.to(u.day), 1e6 * fluxes)
    plt.gca().set(xlabel='Time [days]', ylabel='Relative Flux [ppm]', title='G9V star')
    plt.show()

We can see the shift in the p-mode oscillations relative to the solar ones above
if we plot the power spectrum::

    import matplotlib.pyplot as plt
    import numpy as np
    import astropy.units as u
    from astropy.constants import M_sun, L_sun, R_sun

    from shocksgo import generate_stellar_fluxes, power_spectrum

    # Stellar properties
    M = 0.9 * M_sun
    T_eff = 5340 * u.K
    L = 0.56 * L_sun
    R = 0.876 * R_sun

    times, fluxes, kernel = generate_stellar_fluxes(duration=10*u.day, M=M, T_eff=T_eff, R=R, L=L, cadence=60*u.s)

    freq, power = power_spectrum(fluxes, d=60)

    plt.semilogy(freq * 1e6, power, ',', label='Samples')
    plt.semilogy(freq * 1e6, 1e6 * kernel.get_psd(2*np.pi*freq)/(2*np.pi), alpha=0.7, label='Kernel')
    plt.ylim([1e-5, 1e-1])
    plt.xlim([2500, 5000])
    plt.gca().set(xlabel='Frequency [$\mu$Hz]', ylabel='Power')
    plt.show()


.. plot::

    import matplotlib.pyplot as plt
    import numpy as np
    import astropy.units as u
    from astropy.constants import M_sun, L_sun, R_sun

    from shocksgo import generate_stellar_fluxes, power_spectrum

    # Stellar properties
    M = 0.9 * M_sun
    T_eff = 5340 * u.K
    L = 0.56 * L_sun
    R = 0.876 * R_sun

    times, fluxes, kernel = generate_stellar_fluxes(duration=10*u.day, M=M, T_eff=T_eff, R=R, L=L, cadence=60*u.s)

    freq, power = power_spectrum(fluxes, d=60)

    plt.semilogy(freq * 1e6, power, ',', label='Samples')
    plt.semilogy(freq * 1e6, 1e6 * kernel.get_psd(2*np.pi*freq)/(2*np.pi), alpha=0.7, label='Kernel')
    plt.ylim([1e-5, 1e-1])
    plt.xlim([2500, 5000])
    plt.gca().set(xlabel='Frequency [$\mu$Hz]', ylabel='Power')
    plt.show()