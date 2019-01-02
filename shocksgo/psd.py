import numpy as np

__all__ = ['power_spectrum']


def power_spectrum(samples, d=60):
    """
    Compute the power spectrum of ``samples``.

    Parameters
    ----------
    samples : `~numpy.ndarray`
        Samples
    d : float
        Time between samples [s].

    Returns
    -------
    freq : `~numpy.ndarray`
        Frequencies
    power : `~numpy.ndarray`
        Power at each frequency
    """
    fft = np.fft.rfft(samples)
    power = (fft * np.conj(fft)).real / len(samples) #**2
    freq = np.fft.rfftfreq(len(samples), d)
    return freq, power
