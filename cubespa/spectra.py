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
    return np.array(reg_flux)


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
