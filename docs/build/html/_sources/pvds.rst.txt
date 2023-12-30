Position-Velocity diagrams
==========================

Position-velocity diagrams (PVDs) plot the strength of velocity features along some specified path, most often
along the major and minor axes. A disk with only circular motions will show the rotation curve along the 
major axis, and should be flat along the minor axis. Therefore, PVDs are useful tools to find non-circular
velocity features without explicit modeling of components.

To generate a PVD, do the following:

>>> pvd = cubespa.gen_pvd(blob, (15, 12), 25, 0, width=8)
>>> cubespa.plotting.pvd_plot(blob, pvd, vmin=-125)

Wich will spit out both the map of the region with the overlaid PVD path (with proper width displayed), as
well as the PVD itself. In the example, the region of the PVD cutting through the "blob" shows emission in
blue-shifted velocity channels, consistent with what is seen in the moment 1 mean velocity map.

.. image:: plots/example_pvd.png
  :width: 600
  :alt: Alternative text