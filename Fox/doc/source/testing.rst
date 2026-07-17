*******
Testing
*******

Non-GUI test suite
==================

The repository includes a non-GUI test suite in the top-level ``test/`` directory.
It focuses on library and CLI workflows and avoids the wxWidgets GUI build/test path.

Run tests from ``ObjCryst/``
============================

.. code-block:: shell

  cd ObjCryst
  ln -sf rules-gnu.mak rules.mak
  make test wxcryst=0 opengl=0 fftw=0 cod=0

Run tests from repository root
==============================

.. code-block:: shell

  make -C test all

Current coverage
================

* unit tests for core classes and workflows (crystal creation, single-crystal and powder diffraction, indexing, CIF parsing, powder data import)
* quantitative ground-truth comparisons for simulated single-crystal and powder data (x-ray and neutron)
* powder profile model coverage (pseudo-Voigt Gaussian and Lorentzian limits, anisotropic profile)
* non-GUI integration tests for ``Fox-nogui`` including tutorial datasets
* regression checks for CLI output naming behavior
