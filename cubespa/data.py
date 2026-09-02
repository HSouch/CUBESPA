import numpy as np
import logging
import copy

from astropy.wcs import WCS
from astropy.convolution import Tophat2DKernel
from astropy.io import fits

from scipy.signal import convolve2d

from . import overlays, plotting

import os


class DataSet:
    def __init__(self, filename=None, data=None, wcs=None, header=None, label=None, dtype=None):
        self.filename = filename

        self.data = data
        self.wcs = wcs
        self.header = header
        self.label = label

        self.dtype = dtype
    
    def arcconv_from_wcs(self):
        wcs = self.wcs.celestial
        return np.sqrt(wcs.pixel_scale_matrix[0,0]**2 + wcs.pixel_scale_matrix[1,0]**2) * 3600
    

    def align_to(self, template, shape_out=None):
        return overlays.align_image(template, self, shape_out=None)


    def __repr__(self):
        dstring = str(self.data.shape) if self.data is not None else "NONE"

        return f"Data: {self.filename} | {dstring}"


class RGBImage(DataSet):
    def __init__(self, data=None, wcs=None, header=None, label=None):
        super().__init__(data, wcs, header, label)
    

class MomentMaps:
    def __init__(self, mom0=None, mom1=None, mom2=None,
                 emom0=None, emom1=None, emom2=None,
                   data_index=0):
        
        for attr, val in zip(["mom0", "mom1", "mom2", "emom0", "emom1", "emom2"],
                              [mom0, mom1, mom2, emom0, emom1, emom2]):
            if type(val) == str:
                self.__setattr__(attr, load_data(val, data_index=data_index))
            else:
                self.__setattr__(attr, val)


    def __repr__(self):
        outstring = ""
        for attr in ["mom0", "mom1", "mom2", "emom0", "emom1", "emom2"]:
            attr = self.__getattribute__(attr)
            if isinstance(attr, DataSet):
                outstring += f"{attr.__repr__()}\n"
        return outstring


def load_data(filename, data_index=0, rgb_index = None, label=None, verbose=True):
    try:
        with fits.open(filename) as HDUList:
            hdu = HDUList[data_index]
            
            if rgb_index is None:
                data = hdu.data
                if len(data.shape) == 4:        # Drop Stokes parameter (could remove this later maybe)
                    data = data[0]
                wcs = WCS(hdu.header)
            else:
                data = hdu.data[rgb_index]
                wcs = WCS(hdu.header).dropaxis(2)

            out = DataSet(filename=filename, data=data, wcs=wcs, header=hdu.header, label=label)
            return out
    except Exception as e:
        if verbose:
            print(f"Error loading data: {e}")
        return None


def load_moment_maps(topdir, data_index=0, load_error_maps=False, verbose=True):
    try:
        for stub in[".fits", ".fits.gz"]:
            mom_0_fn = topdir + f".mom0{stub}"
            if not os.path.exists(mom_0_fn):
                continue

            mom_0_fn = topdir + f".mom0{stub}"
            mom_1_fn = topdir + f".mom1{stub}"
            mom_2_fn = topdir + f".mom2{stub}"
            emom0_fn = topdir + f".emom0{stub}"
            emom1_fn = topdir + f".emom1{stub}"
            emom2_fn = topdir + f".emom2{stub}"
            try:
                moment_maps = MomentMaps(mom0=mom_0_fn, mom1=mom_1_fn, mom2=mom_2_fn,
                                        emom0=emom0_fn, emom1=emom1_fn, emom2=emom2_fn,
                                        data_index=data_index)
            except Exception as e:
                continue
            return moment_maps
    except Exception as e:
        if verbose:
            print(f"Error loading Moment Maps: {e}")
        return None
    
    return None



def handle_data(data, handler, data_index=0, verbose=True):
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
        return DataSet(data=data)
    else:
        return None


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


def gen_bg_mask(a, thresh_min=None, dilate=None, dilate_thresh = 0.7):
    """Generate mask from a map, assuming emission is either nan or below some threshold value.
    Mask can be dilated using the dilate parameter, and be further adjusted with dilate_thresh

    Args:
        a (ndarray): Input data array
        thresh_min (float, optional): Minimum threshold to keep unmasked. Defaults to None.
        dilate (int, optional): Side of tophat kernel to dilate mask. Defaults to None.
        dilate_thresh (float, optional): Threshold (between 0 and 1) to include "edge areas" in the dilation mask. 
        A higher value corresponds to a stronger dilation. Defaults to 0.7.

    Returns:
        ndarray: Boolean mask with same shape as input array.
    """

    mask = np.zeros(a.shape)
    mask[np.isnan(a)] = 1

    if thresh_min is not None:
        mask[a < thresh_min] = 1

    if dilate is not None:
        kernel = Tophat2DKernel(dilate)
        mask = convolve2d(mask, kernel, mode="same")
        mask[mask < dilate_thresh] = 0
        mask[mask >= dilate_thresh] = 1

    return mask.astype(bool)
