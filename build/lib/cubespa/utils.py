from astropy.io import fits
from astropy.wcs import WCS

import numpy as np

import os
    

def line_endpoints(x0, y0, L, theta):
    # Convert the angle from degrees to radians
    theta_radians = np.deg2rad(theta)

    # Calculate the endpoint coordinates
    xs = np.array([x0 + L * np.cos(theta_radians), x0 - L * np.cos(theta_radians)])
    ys = np.array([y0 + L * np.sin(theta_radians), y0 - L * np.sin(theta_radians)])

    return xs, ys


def line_corners(x0, y0, w, L, theta):
    theta = np.deg2rad(90 - theta)
    # Rotation matrix
    rotation_matrix = np.array([[np.cos(theta), -np.sin(theta)], 
                                [np.sin(theta),  np.cos(theta)]])
    
    # Half dimensions
    hw, hL = w / 2, L / 2
    
    # Define corners relative to center before rotation
    corners = np.array([(-hw, -hL), (hw, -hL), (hw, hL), (-hw, hL), (-hw, -hL)])
    
    # Rotate and translate corners
    rotated_corners = np.dot(corners, rotation_matrix) + np.array([x0, y0])
    
    return np.transpose(rotated_corners)


def centre_coords(input_wcs, ra, dec):
    cent_x, cent_y = input_wcs.wcs_world2pix(ra, dec, 0)

    return cent_x, cent_y



def match_wcs_axes(wcs1, wcs2):
    """ Match the axes in WCS axes (mostly for image alignment with cubes and images).
        This method assumes that the ra/dec axes are always at indices 0 and 1.

    Args:
        wcs1 (astropy.wcs): WCS A
        wcs2 (astropy.wcs): WCS B

    Returns:
        _type_: Both WCS objects.
    """
    
    naxis1, naxis2 = wcs1.wcs.naxis, wcs2.wcs.naxis

    if naxis1 == naxis2:
        return wcs1, wcs2
    
    to_reduce, target = np.argmax([naxis1, naxis2]), np.min([naxis1, naxis2])

    wcs_objs = [wcs1, wcs2]

    while wcs_objs[to_reduce].wcs.naxis > target:
        wcs_objs[to_reduce] = wcs_objs[to_reduce].dropaxis(-1)

    return wcs_objs


def bounds_from_moment_map(data, padding=0):
    
    non_nans = np.transpose(np.argwhere(~np.isnan(data)))
    ymin, ymax = np.min(non_nans[0]) - padding, np.max(non_nans[0]) + padding
    xmin, xmax = np.min(non_nans[1]) - padding, np.max(non_nans[1]) + padding

    return int(xmin), int(xmax), int(ymin), int(ymax)


def beam_area(bmaj, bmin):
    """ Get the area of the beam based on a standard elliptical Gaussian beam with major and 
        minor axes bmaj and bmin.

    Args:
        bmaj (float): Major axis of beam
        bmin (float): Minor axis of beam

    Returns:
        float: The beam area.
    """
    return np.pi * bmaj * bmin / (4 * np.log(2))


def RMS(a, sigclip=False):
    # todo: add sigma clipping
    a = np.asarray(a)
    return np.sqrt(np.nanmean(a * a))


def estimate_rms(cube, channel_min, channel_max):

    spec_empty = np.concatenate((cube[:channel_min], cube[-channel_max:]))
    
    return np.nanmean(spec_empty), RMS(spec_empty)


def normalize(a, clip_low = None, clip_high=None, stretch=None):
    a_norm = np.copy(a)

    if clip_low is not None:
        a_norm[a_norm < clip_low] = clip_low

    a_norm -= np.nanmin(a_norm)

    if clip_high is not None:
        a_norm[a_norm > clip_high] = clip_high

    a_norm /= np.nanmax(a_norm)

    if stretch == "square":
        a_norm = a_norm * a_norm


    return a_norm


def imstat(a):
    """ Simple tool to return statistics for a given array.

    Args:
        a (ndarray): Input array to get statistics.

    Returns:
        _type_: _description_
    """
    return {
        "MEDIAN": np.nanmedian(a),
        "STD": np.nanstd(a),
        "MAX": np.nanmax(a),
        "MIN": np.nanmin(a)
    }

def im_bounds(stats, sigma=1):
    if type(sigma) is not tuple:
        sigma = (sigma, sigma)

    low = stats["MEDIAN"] - sigma[0] * stats["STD"]
    high = stats["MEDIAN"] + sigma[1] * stats["STD"]

    return low, high


def normalized_rgb_image(image, sigma=1, stretch=None):
    r,g,b = image

    bounds_r = im_bounds(imstat(r), sigma=sigma)
    bounds_g = im_bounds(imstat(g), sigma=sigma)
    bounds_b = im_bounds(imstat(b), sigma=sigma)

    r = normalize(r, clip_low=bounds_r[0], clip_high=bounds_r[1], stretch=stretch)
    g = normalize(g, clip_low=bounds_g[0], clip_high=bounds_g[1], stretch=stretch)
    b = normalize(b, clip_low=bounds_b[0], clip_high=bounds_b[1], stretch=stretch)

    return np.asarray([r,g,b]).transpose(1,2,0)


def recommended_figsize(a, width=8):

    height = np.round(width * a.shape[0] / a.shape[1], 1)

    return (width, height)


def check_and_make_dir(directory):
    path = os.path.dirname(directory)
    if not os.path.isdir(path):
        os.makedirs(path, exist_ok=True)
    return directory



def check_kwarg(key, default, kwargs: dict):
    if key in kwargs.keys(): 
        return kwargs[key]
    else: 
        return default