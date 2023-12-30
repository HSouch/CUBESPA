import numpy as np
import logging
import copy

from . import utils

from astropy.wcs import WCS

from astropy.io import fits

from . import plotting


class DataSet:
    def __init__(self, data=None, wcs=None, header=None, label=None):
        self.data = data
        self.wcs = wcs
        self.header = header
        self.label = label


class RGBImage(DataSet):
    def __init__(self, data=None, wcs=None, header=None, label=None):
        super().__init__(data, wcs, header, label)
    


class MomentMaps:
    def __init__(self, mom0=None, mom1=None, mom2=None, data_index=0):
        if type(mom0) == str:
            self.mom0 = load_data(mom0, data_index=data_index)
        if type(mom1) == str:
            self.mom1 = load_data(mom1, data_index=data_index)
        if type(mom2) == str:
            self.mom2 = load_data(mom2, data_index=data_index)


def load_data(filename, data_index=0, rgb_index = None, label=None):
    try:
        with fits.open(filename) as HDUList:
            hdu = HDUList[data_index]
            
            if rgb_index is None:
                data = hdu.data
                wcs = WCS(hdu.header)
            else:
                data = hdu.data[rgb_index]
                wcs = WCS(hdu.header).dropaxis(2)

            out = DataSet(data=data, wcs=wcs, header=hdu.header, label=label)
            return out
    except Exception as e:
        print(f"Error loading data: {e}")
        return None


def load_moment_maps(topdir, data_index=0):
    try:
        mom_0_fn = topdir + ".mom0.fits.gz"
        mom_1_fn = topdir + ".mom1.fits.gz"
        mom_2_fn = topdir + ".mom2.fits.gz"

        moment_maps = MomentMaps(mom0=mom_0_fn, mom1=mom_1_fn, mom2=mom_2_fn, data_index=data_index)
        return moment_maps
    except Exception as e:
        print(f"Error loading Moment Maps: {e}")
        return None



def handle_data(data, handler, data_index=0):
    """ Handle incoming data for the CubeSPA object. 
    If not a string, CubeSPA will return the data enclosed in a dataset without additional info.
    Args:
        data (str or cubespa.DataSet): Incoming data. If str, CubeSPA will automatically load it in.
        handler (method): Method to handle data, either set to load_data or load_moment_maps
        data_index (int, optional): Index to find data in. Defaults to 0.
    Returns:
        utils.DataSet : The output DataSet object. 
    """
    if type(data) == str:
        return handler(data, data_index=data_index)
    elif type(data) == np.ndarray:
        logging.warn("Naked data loaded. No WCS or header included!")
        return utils.DataSet(data=data)


def gen_cutout(data, limits):
    """ Generates a cutout from an image based on a set of limits.

    Args:
        data (ndarray): Input image data array
        limits (arr): Limits, assuming [xmin, xmax, ymin, ymax] format

    Returns:
        _type_: _description_
    """
    
    xmin, xmax, ymin, ymax = limits

    return data[ymin:ymax, xmin:xmax]


def gen_cutout(cubespa_obj, cent, size, show_bbox=False):

    # TODO: change WCS coordinates to match new cutout location

    cutout = copy.deepcopy(cubespa_obj)

    if type(size) is int:
        size = (size, size)

    xmin, xmax = cent[1] - size[0], cent[1] + size[0]
    ymin, ymax = cent[0] - size[1], cent[0] + size[1]

    if show_bbox:
        plotting.plot_bbox(cubespa_obj, [ymin, ymax, xmin, xmax])

    new_cube = cubespa_obj.cube.data[:, xmin:xmax, ymin:ymax]
    if cubespa_obj.mom_maps is not None:
        new_mom0 = cubespa_obj.mom_maps.mom0.data[xmin:xmax, ymin:ymax]
        new_mom1 = cubespa_obj.mom_maps.mom1.data[xmin:xmax, ymin:ymax]
        new_mom2 = cubespa_obj.mom_maps.mom2.data[xmin:xmax, ymin:ymax]

    for i, additional in enumerate(cubespa_obj.additional_maps):
        cutout.additional_maps[i].data = additional.data[xmin:xmax, ymin:ymax]


    cutout.cube.data = new_cube

    cutout.limits = [ymin, ymax, xmin, xmax]

    cutout.mom_maps.mom0.data = new_mom0
    cutout.mom_maps.mom1.data = new_mom1
    cutout.mom_maps.mom2.data = new_mom2

    return cutout