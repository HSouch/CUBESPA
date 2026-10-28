# CubeSPA — CUBE SPectral Analysis

Fully-featured, aesthetically pleasing 2D and 3D analysis of astronomical datacubes.

CubeSPA bundles the common steps of working with a spectral cube: loading a cube with its moment maps, aligning ancillary images, making cutouts of regions, extracting spectra, building position-velocity diagrams, and producing clean overlay and RGB figures.

## Installation

```bash
pip install cubespa
```

Requires Python 3.10+. 

## Quick start

```python
import cubespa

# Load a cube, plus moment maps produced with maskmoment
# (e.g. path/to/maskmoment.mom0.fits.gz, .mom1, .mom2)
c = cubespa.CubeSPA("path/to/cube.fits", mom_maps="path/to/maskmoment", limits="auto")
c.plot_moment_maps()

# Cut out a region: centre (x, y) and size (int for square, tuple for rectangle)
blob = cubespa.gen_cutout(c, (345, 290), 15, show_bbox=True)

# Spectrum in an elliptical aperture, with 1-3 sigma RMS levels shaded
blob.create_spectra((15, 13), (6, 6), plot=True)

# Position-velocity diagram along a path through the cutout
pvd = cubespa.gen_pvd(blob, (15, 12), 25, 0, width=8)
cubespa.plotting.pvd_plot(blob, pvd, vmin=-125)
```

### Ancillary data and overlays

```python
# Load ancillary data and align it to parent image
halpha = cubespa.load_data("path/to/halpha.fits", label="HALPHA")
halpha_aligned = cubespa.align_image(c.mom_maps.mom0, halpha)
c.additional_maps.append(halpha_aligned)

# Overlay maps
cubespa.plotting.overlay_plot(c.mom_maps.mom0, halpha_aligned,
                              lims=c.limits, levels=[10, 20, 100, 150])

# Create RGB imaging
hst = cubespa.align_image(c.mom_maps.mom0, cubespa.load_data("path/to/hst.fits", rgb_index=None, label="HST"))
norm = cubespa.normalized_rgb_image(hst.data, sigma=(2, 5))
cubespa.plotting.rgb_overlay(norm, c.mom_maps.mom0, levels=[0.02, 0.05, 0.1, 1],
                             lims=c.limits, colors="white")
```

## Features

- **Moment maps:** load and plot maskmoment outputs, with automatic bounding of valid data (`limits="auto"`).
- **Cutouts:** trim the cube, moment maps and ancillary maps to a region in one call.
- **Spectra:** aperture spectra, sigma-clipped RMS, and channel-range SNR (`cubespa.calc_snr`).
- **PVDs:** position-velocity diagrams along arbitrary paths.
- **Alignment and overlays:** reproject ancillary images onto the cube, then draw contour and RGB overlays.
- **PSF overlays:** overplot a shifted PSF to hunt for sidelobe artefacts (`psf=...`, `cubespa.plotting.plot_psf_overlay`).

## Documentation

Full documentation, with example figures, is in [`docs/`](docs/source). 

## Contact

Harrison Souchereau — harrison.souchereau@yale.edu. 

Issues and pull requests are welcome at <https://codeberg.org/HSouch/CUBESPA>.
