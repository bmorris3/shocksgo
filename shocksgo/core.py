import numpy as np
import celerite
from celerite import terms
import astropy.units as u
from astropy.constants import L_sun, M_sun, R_sun
import os

__all__ = ['generate_solar_fluxes', 'generate_stellar_fluxes']

dirname = os.path.dirname(os.path.abspath(__file__))
parameter_vector = np.loadtxt(os.path.join(dirname, 'data',
                                           'parameter_vector.txt'))


@u.quantity_input(cadence=u.s, duration=u.s)
def generate_solar_fluxes(duration, cadence=60*u.s,
                          parameter_vector=parameter_vector):
    """
    Generate an array of fluxes with zero mean which mimic the power spectrum of
    the SOHO/VIRGO SPM observations.
    
    Parameters
    ----------
    size : int
        Number of fluxes to generate. Note: Assumes ``size``>>500.
    cadence : `~astropy.units.Quantity`
        Length of time between fluxes
    
    Returns
    -------
    times : `~astropy.units.Quantity`
        Array of times at cadence ``cadence`` of length ``size``
    fluxes : `~numpy.ndarray`
        Array of fluxes at cadence ``cadence`` of length ``size``.
    kernel : `~celerite.terms.TermSum`
        Celerite kernel used to approximate the solar power spectrum.
    """
    ##########################
    # Assemble celerite kernel
    ##########################
    parameter_vector = np.copy(parameter_vector)

    nterms = len(parameter_vector)//3

    kernel = terms.SHOTerm(log_S0=0, log_omega0=0, log_Q=0) 

    for term in range(nterms-1): 
        kernel += terms.SHOTerm(log_S0=0, log_omega0=0, log_Q=0)

    kernel.set_parameter_vector(parameter_vector)

    gp = celerite.GP(kernel)

    # For more efficient computation for large datasets, split into chunks:
    if duration.to(u.s).value >= 1e5:
        divisor = 500
        x = np.arange(0, duration.to(u.s).value//divisor, cadence.to(u.s).value)
        times = np.arange(0, duration.to(u.s).value,
                          cadence.to(u.s).value) * u.s
    else:
        divisor = 1
        times = np.arange(0, duration.to(u.s).value,
                          cadence.to(u.s).value) * u.s
        x = times.value

    gp.compute(x, check_sorted=False)

    ###################################
    # Get samples with the kernel's PSD
    ###################################

    y = gp.sample(size=divisor+1 if divisor != 1 else divisor)

    # Reassemble the chunks
    y_concatenated = []

    for i, yi in enumerate(y):
        xi = np.arange(len(yi))
        fit = np.polyval(np.polyfit(xi - xi.mean(), yi, 1),
                         xi-xi.mean())
        yi -= fit

        if i == 0: 
            y_concatenated.append(yi)
        else: 
            offset = yi[0] - y_concatenated[i-1][-1]
            y_concatenated.append(yi - offset)

    y_concatenated = np.hstack(y_concatenated)[:len(times)]
    
    x_c = np.arange(len(y_concatenated))
    
    y_concatenated -= np.polyval(np.polyfit(x_c - x_c.mean(),
                                            y_concatenated, 1),
                                 x_c - x_c.mean())
    
    return times, y_concatenated, kernel


@u.quantity_input(duration=u.s, cadence=u.s, M=u.kg, T_eff=u.K, L=u.W, R=u.m)
def generate_stellar_fluxes(duration, M, T_eff, R, L, cadence=60*u.s,
                            parameter_vector=parameter_vector):
    """
    Generate an array of fluxes with zero mean which mimic the power spectrum of
    the SOHO/VIRGO SPM observations, scaled for a star with a given mass,
    effective temperature and luminosity.
    
    Parameters
    ----------
    size : int
        Number of fluxes to generate. Note: Assumes ``size``>>500.
    M : `~astropy.units.Quantity`
        Stellar mass
    T_eff : `~astropy.units.Quantity`
        Stellar effective temperature
    R : `~astropy.units.Quantity`
        Stellar radius
    L : `~astropy.units.Quantity`
        Stellar luminosity
    cadence : `~astropy.units.Quantity`
        Length of time between fluxes
    
    Returns
    -------
    times : `~astropy.units.Quantity`
        Array of times at cadence ``cadence`` of length ``size``
    fluxes : `~numpy.ndarray`
        Array of fluxes at cadence ``cadence`` of length ``size``.
    kernel : `~celerite.terms.TermSum`
        Celerite kernel used to approximate the stellar power spectrum.
    """
    ##########################
    # Scale p-mode frequencies
    ##########################
    parameter_vector = np.copy(parameter_vector)

    # Scale frequencies
    tunable_amps = np.exp(parameter_vector[::3][2:])
    tunable_freqs = np.exp(parameter_vector[2::3][2:]) * 1e6/2/np.pi
    peak_ind = np.argmax(tunable_amps)
    peak_freq = tunable_freqs[peak_ind] # 3090 uHz in Huber 2011
    delta_freqs = tunable_freqs - peak_freq
    
    T_eff_solar = 5777 * u.K

    # Huber 2012 Eqn 3
    delta_nu_factor = (M/M_sun)**0.5 * (R/R_sun)**(-3/2)
    # Huber 2012 Eqn 4
    nu_factor = (M/M_sun) * (R/R_sun)**-2 * (T_eff/T_eff_solar)**-0.5

    new_peak_freq = nu_factor * peak_freq
    new_delta_freqs = delta_freqs * delta_nu_factor

    new_freqs = new_peak_freq + new_delta_freqs

    new_log_omegas = np.log(2*np.pi*new_freqs*1e-6).value

    parameter_vector[2::3][2:] = new_log_omegas

    # Scale amplitudes of p-mode oscillations following Huber 2011
    # Huber 2011 Eqn 8:
    c = (T_eff/(5934 * u.K))**0.8
    c_sun = ((5777 * u.K)/(5934 * u.K))**0.8
    s = 0.886
    r = 2
    t = 1.89

    # Huber 2011 Eqn 9:
    pmode_amp_star = (L/L_sun)**s / ((M/M_sun)**t * T_eff.value**(r-1) * c)
    pmode_amp_sun = (L_sun/L_sun)**s / ((M_sun/M_sun)**t * 5777**(r-1) * c_sun)
    pmode_amp_factor = pmode_amp_star / pmode_amp_sun

    new_pmode_amps = np.log(np.exp(parameter_vector[0::3][2:]) *
                            pmode_amp_factor)
    parameter_vector[0::3][2:] = new_pmode_amps

    #############################
    # Scale granulation frequency
    #############################

    # Kallinger 2014 pg 12:
    tau_eff_factor = (new_peak_freq/peak_freq)**-0.89
    # Kallinger 2014 pg 12:
    granulation_amplitude_factor = (new_peak_freq/peak_freq)**-0.56
    parameter_vector[5] = np.log(np.exp(parameter_vector[5]) / tau_eff_factor)
    parameter_vector[3] = np.log(np.exp(parameter_vector[3]) *
                                 granulation_amplitude_factor)

    ##########################
    # Assemble celerite kernel
    ##########################

    nterms = len(parameter_vector)//3

    kernel = terms.SHOTerm(log_S0=0, log_omega0=0, log_Q=0)

    for term in range(nterms-1):
        kernel += terms.SHOTerm(log_S0=0, log_omega0=0, log_Q=0)

    kernel.set_parameter_vector(parameter_vector)

    gp = celerite.GP(kernel)

    # For more efficient computation for large datasets, split into chunks:
    if duration.to(u.s).value >= 1e5:
        divisor = 500
        x = np.arange(0, duration.to(u.s).value//divisor, cadence.to(u.s).value)
        times = np.arange(0, duration.to(u.s).value,
                          cadence.to(u.s).value) * u.s
    else:
        divisor = 1
        times = np.arange(0, duration.to(u.s).value,
                          cadence.to(u.s).value) * u.s
        x = times.value

    gp.compute(x, check_sorted=False)

    ###################################
    # Get samples with the kernel's PSD
    ###################################

    y = gp.sample(size=divisor+1 if divisor != 1 else divisor)

    y_concatenated = []

    for i, yi in enumerate(y):
        xi = np.arange(len(yi))
        fit = np.polyval(np.polyfit(xi - xi.mean(), yi, 1), xi-xi.mean())
        yi -= fit

        if i == 0:
            y_concatenated.append(yi)
        else:
            offset = yi[0] - y_concatenated[i-1][-1]
            y_concatenated.append(yi - offset)

    y_concatenated = np.hstack(y_concatenated)[:len(times)]

    x_c = np.arange(len(y_concatenated))

    y_concatenated -= np.polyval(np.polyfit(x_c - x_c.mean(),
                                            y_concatenated, 1),
                                 x_c - x_c.mean())

    return times, y_concatenated, kernel
