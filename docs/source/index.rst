.. CubeSPA documentation master file, created by
   sphinx-quickstart on Thu Dec 28 19:25:16 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

*******
CubeSPA
*******

**CubeSPA** is a fully-featured set of utilities to process, analyze, and display useful information of datacubes, particularly
of astronomical data. It fully supports utilities such as moment map creation, spectral analysis, position-velocity diagrams, and
many others. 

To use **CubeSPA**, install it using 

``pip install cubespa`` (not yet, but soon)

To begin using CubeSPA, create a ``cubespa.CubeSPA`` object with the filename for your cube in the following way:

    >>> filename = "path/to/cube.fits"
    >>> c = cubespa.CubeSPA(cube_fn)

You can load in moment maps (assuming the convention from  `maskmoment <https://github.com/tonywong94/maskmoment>`_ ) with the following. If your maskmoment output is 
``path/to/maskmoment.mom0.fits.gz``, for the moment 0 map (.mom1, .mom2 for the others), these are loaded as follows below. With
moment maps loaded, you can also create a bounding box around "valid" data by calling the ``limits = "auto"`` feature.

    >>> filename = "path/to/cube.fits"
    >>> mommaps = "path/to/maskmoment"
    >>> c = cubespa.CubeSPA(cube_fn, mom_maps=mommaps, limits="auto")

See getting_started.rst for more information.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

Loading Data
============
.. toctree::
    :maxdepth: 1
    
    getting_started.rst

2D Analysis
==============

.. toctree::
    :maxdepth: 1
        
    momentmaps.rst
    overlays.rst

Cube Analysis
=============

.. toctree::
    :maxdepth: 1
        
    spectra.rst
    pvds.rst

Coming soon
===========

.. toctree::
   :maxdepth: 1

   modeling.rst
    
    

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
