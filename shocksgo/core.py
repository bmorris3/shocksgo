import numpy as np
import celerite
from celerite import terms
import astropy.units as u
from astropy.constants import L_sun, M_sun
import os

__all__ = ['generate_solar_fluxes', 'generate_stellar_fluxes']

dirname = os.path.dirname(os.path.abspath(__file__))
parameter_vector = np.loadtxt(os.path.join(dirname, 'data',
                                           'parameter_vector.txt'))


@u.quantity_input(cadence=u.s)
def generate_solar_fluxes(size, cadence=60*u.s,
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
    
    x = np.arange(0, size//500, cadence.to(u.s).value) 
    gp.compute(x, check_sorted=False)

    ###################################
    # Get samples with the kernel's PSD
    ###################################

    y = gp.sample(500)
    
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

    y_concatenated = np.hstack(y_concatenated)
    
    x_c = np.arange(len(y_concatenated))
    
    y_concatenated -= np.polyval(np.polyfit(x_c - x_c.mean(), y_concatenated, 1), 
                                 x_c - x_c.mean())
    
    return y_concatenated, kernel


@u.quantity_input(cadence=u.s, M=u.kg, T_eff=u.K, L=u.W)
def generate_stellar_fluxes(size, M, T_eff, L, cadence=60*u.s,
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
    L : `~astropy.units.Quantity`    
        Steller luminosity
    cadence : `~astropy.units.Quantity`
        Length of time between fluxes
    
    Returns
    -------
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
    peak_freq = tunable_freqs[peak_ind]
    delta_freqs = tunable_freqs - peak_freq
    
    T_eff_solar = 5777 * u.K
    
    # Huber 2011 Eqn 1
    nu_factor = ( (M/M_sun) * (T_eff/T_eff_solar)**3.5 / (L/L_sun) )
    # Huber 2011 Eqn 2
    delta_nu_factor = ( (M/M_sun)**0.5 * (T_eff/T_eff_solar)**3 /
                        (L/L_sun)**0.75 )
    
    new_peak_freq = nu_factor * peak_freq
    new_delta_freqs = delta_freqs * delta_nu_factor

    new_freqs = new_peak_freq + new_delta_freqs

    new_log_omegas = np.log(2*np.pi*new_freqs*1e-6).value

    parameter_vector[2::3][2:] = new_log_omegas

    # Scale amplitudes
    c = (T_eff/(5934 * u.K))**0.8
    c_sun = ((5777 * u.K)/(5934 * u.K))**0.8
    s = 0.886
    r = 2
    t = 1.89
    pmode_amp_factor = (L/L_sun)**s / ((M/M_sun)**t * T_eff.value**(r-1) * c)
    pmode_amp_factor_sun = (L_sun/L_sun)**s / ((M_sun/M_sun)**t * 5777**(r-1)
                                               * c_sun)

    new_pmode_amps = np.log(np.exp(parameter_vector[0::3][2:]) *
                            pmode_amp_factor/pmode_amp_factor_sun)
    parameter_vector[0::3][2:] = new_pmode_amps

    #############################
    # Scale granulation frequency
    #############################

    tau_eff_factor = (new_peak_freq/peak_freq)**-0.89
    granulation_amplitude_factor = (new_peak_freq/peak_freq)**-0.5
    parameter_vector[4] = np.log(np.exp(parameter_vector[4]) * tau_eff_factor)
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

    x = np.arange(0, size//500, cadence.to(u.s).value)
    gp.compute(x, check_sorted=False)


    ###################################
    # Get samples with the kernel's PSD
    ###################################

    y = gp.sample(500)

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

    y_concatenated = np.hstack(y_concatenated)

    x_c = np.arange(len(y_concatenated))

    y_concatenated -= np.polyval(np.polyfit(x_c - x_c.mean(),
                                            y_concatenated, 1),
                                 x_c - x_c.mean())

    return y_concatenated, kernel
