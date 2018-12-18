====================================================================================
Simple Harmonic Oscillator Celerite Kernels for Stellar Granulation and Oscillations
====================================================================================

This is the documentation for ``shocksgo``. The goal of ``shocksgo`` is to
generate light curves of stars accounting for the effects of granulation,
super-granulation and p-mode oscillations.

You can view the source code and/or contribute to ``shocksgo`` via
`GitHub <https://github.com/bmorris3/shocksgo>`_.

#############
Documentation
#############

.. toctree::
  :maxdepth: 2

  shocksgo/installation.rst
  shocksgo/gettingstarted.rst
  shocksgo/index.rst



########
Overview
########

Methods
=======

We compute these light curves efficiently by taking advantage of
`celerite <http://celerite.readthedocs.io>`_, a fast Gaussian process regression
package, which we use to approximate solar and stellar power spectrum densities
with sums of
`simple harmonic oscillator <https://celerite.readthedocs.io/en/stable/python/kernel/#celerite.terms.SHOTerm>`_
(SHO) kernels of the form:

.. math::

    S(\omega) = \sqrt{\frac{2}{\pi}} \frac{S_0\,\omega_0^4}
    {(\omega^2-{\omega_0}^2)^2 + {\omega_0}^2\,\omega^2/Q^2}

where :math:`\omega = 2\pi f` is the angular frequency. We use one
SHO kernel term for super/meso-granulation, another for ordinary granulation,
and about 100 terms for the comb of p-mode peaks.

Scaling relations for p-modes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For computation of stellar p-mode oscillation frequencies, we use the scaling
relations found in
`Huber et al. (2011) <http://adsabs.harvard.edu/abs/2011ApJ...743..143H>`_ and
references therein (e.g. Kjeldsen & Bedding 1995), namely:

.. math::

    \nu_\textrm{max} \approx \frac{M / M_\odot (T_\textrm{eff}/
    T_{\textrm{eff},\odot})^{3.5} }{L/L_\odot} \nu_{\textrm{max}, \odot}

(Equation 1), and

.. math::

    \Delta \nu_\textrm{max} \approx \frac{(M / M_\odot)^{0.5} (T_\textrm{eff}/
    T_{\textrm{eff},\odot})^{3} }{(L/L_\odot)^{0.75}} \Delta \nu_{\odot}

(Equation 2).

Scaling relations for granulation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For computation of stellar granulation frequency, we use the scaling relation
found in
`Kallinger et al. (2014) <http://adsabs.harvard.edu/abs/2014A%26A...570A..41K>`_:

.. math::

    \tau_\textrm{eff} \propto \nu^{-0.89}_\textrm{max},

where :math:`\tau_\textrm{eff}` is the characteristic granulation timescale. The
amplitude of the granulation scales as

.. math::

    a \propto \nu^{-0.5}_\textrm{max}.
