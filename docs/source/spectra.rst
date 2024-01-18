Spectral Analysis
=================

To generate spectra, the user specifies the location of the aperture (elliptical by default), and the shape. 
Continuing with our ``blob`` example, we generate spectra like so:

>>> blob.create_spectra((15, 13), (6, 6), plot=True)

This creates the following plot, showing the location of the spectral features in either velocity or channel
space. CubeSPA uses sigma clipping to determine an RMS of the spectra, and then plots the 1,2 and 3-sigma levels 
as a shaded grey region.


.. image:: plots/blob_spectra.png
  :width: 800
  :alt: Alternative text


SNR Analysis
------------

To manually check the SNR of a given region (say, between channels 20 and 40 of the array ``spec``),
you can with the following method:

>>> chan_min, chan_max = 20, 40
>>> snr = cubespa.calc_snr(spec, chan_min, chan_max)

The signal to noise ratio between the two channels is given by:

.. math:: SNR = \frac{\Sigma F}{RMS_{spec} * \sqrt{n_{chan}}}





.. automodule:: cubespa.plotting.spectra_plots
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: cubespa.spectra
   :members:
   :undoc-members:
   :show-inheritance: