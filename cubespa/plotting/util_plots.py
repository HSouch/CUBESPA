from matplotlib import pyplot as plt

from .. import utils

import numpy as np

def limit_plot(cubespa_obj):
    xmin, xmax, ymin, ymax = cubespa_obj.limits

    mom0 = cubespa_obj.mom_maps.mom0.data
    plt.figure()
    plt.imshow(mom0, origin="lower")

    plt.axvline(xmin)
    plt.axvline(xmax)

    plt.axhline(ymin)
    plt.axhline(ymax)
    plt.tight_layout()
    plt.show()


def plot_bbox(cubespa_obj, lims):

    xmin, xmax, ymin, ymax = lims

    if cubespa_obj.mom_maps is not None:
        data = cubespa_obj.mom_maps.mom0.data
    else:
        data = np.nansum(cubespa_obj.cube.data, axis=0)

    plt.imshow(data, origin="lower")
        
    plt.plot([xmin, xmax, xmax, xmin, xmin], [ymin, ymin, ymax, ymax, ymin])

    plt.xlim(cubespa_obj.limits[0], cubespa_obj.limits[1])
    plt.ylim(cubespa_obj.limits[2], cubespa_obj.limits[3])  

    plt.show()


def plot_rgb(rgb, lims=None, outname=None):
    """ Plot an RGB image using matplotlib

    Args:
        rgb (nxmx3 array): RGB image formatted for matplotlib
        lims (arr, optional): x and y limits for plotting. Defaults to None.
        outname (str, optional): Output filename. If not, show plot instead of
            save figure. Defaults to None.
    """
    figsize = utils.recommended_figsize(rgb)

    plt.figure(figsize=figsize)
    plt.imshow(rgb, origin="lower")

    if lims is not None:
        xmin, xmax, ymin, ymax = lims
        plt.xlim(xmin, xmax)
        plt.ylim(ymin, ymax)

    plt.tight_layout()

    if outname is not None:
        outname = utils.check_and_make_dir(outname)
        plt.savefig(outname, dpi=200)
    else:
        plt.show()