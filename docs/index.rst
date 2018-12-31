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
`simple harmonic oscillator
<https://celerite.readthedocs.io/en/stable/python/kernel/#celerite.terms.SHOTerm>`_
(SHO) kernels of the form:

.. math::

    S(\omega) = \sqrt{\frac{2}{\pi}} \frac{S_0\,\omega_0^4}
    {(\omega^2-{\omega_0}^2)^2 + {\omega_0}^2\,\omega^2/Q^2}

where :math:`\omega = 2\pi f` is the angular frequency. We use one
SHO kernel term for super/meso-granulation, another for ordinary granulation,
and about 50 terms for the comb of p-mode peaks.

Scaling relations for p-modes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For computation of stellar p-mode oscillation frequencies, we use the scaling
relations found in
`Huber et al. (2012) <http://adsabs.harvard.edu/abs/2012ApJ...760...32H>`_ and
references therein
(e.g. `Kjeldsen & Bedding 1995
<http://adsabs.harvard.edu/abs/1995A%26A...293...87K>`_ ),
namely Equation 4:

.. math::

     \nu_\textrm{max} \propto M R^{-2} T_{\rm eff}^{-1/2},

and Equation 3

.. math::

    \Delta \nu_\textrm{max} \propto M^{1/2} R^{-3/2}.


The amplitude scaling of the p-mode oscillations is given by Equation 9 of
`Huber et al. (2011) <http://adsabs.harvard.edu/abs/2011ApJ...743..143H>`_:

.. math::

    A \propto \frac{L^s}{M^t T_\textrm{eff}^{r-1} c(T_\textrm{eff})}

where :math:`r = 2`, :math:`s = 0.886`, :math:`t = 1.89` and

.. math::

    c(T_\textrm{eff}) = \left( \frac{T_\textrm{eff}}{5934 \textrm{K}}
    \right)^{0.8}.


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

    a \propto \nu^{-2}_\textrm{max}.

(`Kjeldsen & Bedding, 2011
<http://adsabs.harvard.edu/abs/2011A%26A...529L...8K>`_).