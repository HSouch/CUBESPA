Making Overlays 
===============

Overlays are useful features to display isophotal contours of one image onto another.

Once you have aligned datasets (see getting_started), you can create an overlay by using:

>>> cubespa.plotting.overlay_plot(c.mom_maps.mom0, c.additional_maps[0], 
...                               lims=c.limits, levels=[10, 20, 100, 150])

where in this instance, we are showing the moment 0 CO distribution with H-alpha contours overlaid in black. 
The output plot looks like this:

.. image:: plots/test_overlay.png
  :width: 600
  :alt: Alternative text


For an RGB image, use the following instead:

>>> norm = cubespa.normalized_rgb_image(hst_interp.data, sigma=(2, 1))
>>> cubespa.plotting.rgb_overlay(norm, c.mom_maps.mom0, 
...                              levels=[0.02, 0.05, 0.1, 1], lims=blob.limits, colors="white")

.. image:: plots/rgb_overlay.png
  :width: 600
  :alt: Alternative text


which plots CO contours on top of the HST RGB image, which for this region reveals
a small compact stellar feature at the head of the  "blob".


PSF Overlays for Sidelobe Checking
----------------------------------

If your cube/cleaning has a PSF that you have imported into cubespa with the ``psf=psf_fn`` feature, you can overplot it using the following
code:

>>> c = cubespa.CubeSPA(cube=cube_fn, psf=psf_fn, limits="auto")
>>> cubespa.plotting.plot_psf_overlay(c, y0=102, x0=25)

.. image:: plots/psf_overlay.png
  :width: 400
  :alt: Alternative text

Note that in this case, the psf has been shifted down by 102 pixels, and to the right by 25 pixels. This overlay is especially helpful when hunting 
for sidelobe features in a dataset, which can be revealed by looking at contours in the PSF.