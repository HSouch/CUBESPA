from photutils.aperture import CircularAnnulus, CircularAperture, EllipticalAnnulus, EllipticalAperture, aperture_photometry
import numpy as np

from . import plotting

from astropy.stats import sigma_clipped_stats


def create_aperture(cubespa_obj, position, shape, aper_type="elliptical", plot=False):

    """ Generate photutils aperture of desired type, position, and shape.

    Returns:
        photutils.aperture : Photutils aperture
    """

    if aper_type == "elliptical":
        aper = EllipticalAperture(positions=[position], a=shape[0], b=shape[1])
    elif aper_type == "circular":
        aper = CircularAperture(positions=[position], r=shape)
    
    if plot:
        plotting.plot_spectra(cubespa_obj.mom_maps.mom0.data, aper)

    return aper


def get_spectra(cube, aper):
    """ Get the spectra through a datacube at the position and size of a given aperture.

    Args:
        cube (ndarray): _description_
        aper (photutils aperture): Elliptical or circular aperture/annulus.

    Returns:
        _type_: _description_
    """
    reg_flux = []
    for frame in cube:
        phot = aperture_photometry(frame, aper)

        reg_flux.append(float(phot["aperture_sum"]))
    return np.array(reg_flux) / aper.area


def analyze_spectra(spec, sigma=2, cmin=None, cmax=None):

    if cmin is not None:
        spec[:cmin] = np.nan
    if cmax is not None:
        spec[cmax:] = np.nan

    spec_med, spec_std = sigma_clipped_stats(spec, sigma=sigma)[1:]

    return spec_med, spec_std


def multi_spec(cubespa_obj, spec_info):
    """ Generate a list of apertures and spectra for better diagnostics.

    Args:
        cubespa_obj (cubespa.CubeSPA): Input CubeSPA object to pull spectra from.
        spec_info (_type_): A list of spectra position and sizes, for example, to get two spectra, you would 
        have:
            spec_info = [[(p1, p2), (s1, s2)],
                         [(p3, p4), (s3, s4)]]
    Returns:
        tuple(array): A list of apertures and a list of spectra.
    """
    aper_list, spec_list = [], []

    for spec in spec_info:
        position, shape = spec
        # print(position, shape)

        aper = create_aperture(cubespa_obj, position, shape)
        spec = get_spectra(cubespa_obj.cube.data, aper)

        aper_list.append(aper)
        spec_list.append(spec)
    
    return aper_list, spec_list


def calc_snr(spec, chan_min, chan_max):
    stats = analyze_spectra(spec)
    n_chan = chan_max - chan_min
    spec_slice = spec[chan_min:chan_max]

    flux_sum = np.nansum(spec_slice)
    return flux_sum / (stats[1] * np.sqrt(n_chan))


def align_apertures(aper_list, wcs1, wcs2):
    """ Align and resize a set of apertures from wcs1 to wcs2

    Args:
        aper_list (list): List of apertures in the following format:
            [(p1, p2), (s1, s2)] where the first tuple is the (x,y) position and the second tuple is the (x,y) height.
        wcs1 (astropy.wcs.WCS): WCS that aper_list apertures are already aligned to
        wcs2 (astropy.wcs.WCS): WCS to align apertures to

    Returns:
        _type_: _description_
    """

    apers_out = []

    s_ratio = wcs2.wcs.cdelt / wcs1.wcs.cdelt

    for n in aper_list:
        (p1, p2), (s1, s2) = n
        
        x1, y1 = wcs1.wcs_pix2world(p1, p2, 1)
        
        x2, y2 = wcs2.wcs_world2pix(x1, y1, 1)

        apers_out.append([(int(x2), int(y2)), 
                          (np.round(s1 * s_ratio[0], 2), np.round(s2 * s_ratio[1], 2))])
        
    return apers_out 


def spectra_comparison(cubecomp, aper_list, outname=None, plot_ticks=True, chan_ranges=None):
    print(len(aper_list))
    apertures_aligned = align_apertures(aper_list, cubecomp.cube1.mom_maps.mom0.wcs, cubecomp.cube2.mom_maps.mom0.wcs)

    a1, s1 = multi_spec(cubecomp.cube1, aper_list)
    a2, s2 = multi_spec(cubecomp.cube2, apertures_aligned)

    plotting.spectra_comparison(cubecomp, a1, a2, s1, s2, cmap="rainbow", outname=outname, chan_ranges=chan_ranges, plot_ticks=plot_ticks)