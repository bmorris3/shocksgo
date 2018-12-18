import numpy as np
import celerite
from celerite import terms
import astropy.units as u
from astropy.constants import L_sun, M_sun
import os

__all__ = ['generate_solar_fluxes', 'generate_stellar_fluxes']

dirname = os.path.dirname(os.path.abspath(__file__))
parameter_vector = np.loadtxt(os.path.join(dirname, 'data', 'parameter_vector.txt'))


@u.quantity_input(cadence=u.s)
def generate_solar_fluxes(size, cadence=60*u.s, parameter_vector=parameter_vector):
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
    y : `~numpy.ndarray`
        Array of fluxes at cadence ``cadence`` of length ``size``.
    """
    nterms = len(parameter_vector)//3

    kernel = terms.SHOTerm(log_S0=0, log_omega0=0, log_Q=0) 

    for term in range(nterms-1): 
        kernel += terms.SHOTerm(log_S0=0, log_omega0=0, log_Q=0)

    kernel.set_parameter_vector(parameter_vector)

    gp = celerite.GP(kernel)
    
    x = np.arange(0, size//500, cadence.to(u.s).value) 
    gp.compute(x, check_sorted=False)

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
def generate_stellar_fluxes(size, M, T_eff, L, cadence=60*u.s, parameter_vector=parameter_vector):
    """
    Generate an array of fluxes with zero mean which mimic the power spectrum of
    the SOHO/VIRGO SPM observations.
    
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
    y : `~numpy.ndarray`
        Array of fluxes at cadence ``cadence`` of length ``size``.
    """
    parameter_vector = parameter_vector.copy()
    
    tunable_amps = np.exp(parameter_vector[::3][2:])
    tunable_freqs = np.exp(parameter_vector[2::3][2:]) * 1e6/2/np.pi
    peak_ind = np.argmax(tunable_amps)
    peak_freq = tunable_freqs[peak_ind]
    delta_freqs = tunable_freqs - peak_freq
    
    T_eff_solar = 5777 * u.K 
    nu_max_sun = peak_freq * u.uHz
    delta_nu_sun = 135.1 * u.uHz
    
    # Huber 2011 Eqn 1
    nu_factor = ( (M/M_sun) * (T_eff/T_eff_solar)**3.5 / (L/L_sun) )
    # Huber 2011 Eqn 2
    delta_nu_factor = ( (M/M_sun)**0.5 * (T_eff/T_eff_solar)**3 / (L/L_sun)**0.75 )
    
    new_peak_freq = nu_factor * peak_freq
    new_delta_freqs = delta_freqs * delta_nu_factor

    new_peak_freq, new_delta_freqs
    new_freqs = new_peak_freq + new_delta_freqs

    new_log_omegas = np.log(2*np.pi*new_freqs*1e-6).value

    parameter_vector[2::3][2:] = new_log_omegas
    
    nterms = len(parameter_vector)//3

    kernel = terms.SHOTerm(log_S0=0, log_omega0=0, log_Q=0) 

    for term in range(nterms-1): 
        kernel += terms.SHOTerm(log_S0=0, log_omega0=0, log_Q=0)

    kernel.set_parameter_vector(parameter_vector)

    gp = celerite.GP(kernel)

    x = np.arange(0, size//500, cadence.to(u.s).value) 
    gp.compute(x, check_sorted=False)

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

