import numpy as np

from .. import utils

from scipy.interpolate import interp1d
from matplotlib import pyplot as plt


def pvd_plot(cubespa_obj, pvd, vmin=-200, vmax=200, vstep=25, filename=None, **kwargs):

    cmap = kwargs.get("cmap", "viridis")
    pvd_color = kwargs.get("pvd_color", "red")

    v_interp = interp1d(cubespa_obj.velocities - 43, np.arange(len(cubespa_obj.velocities)), 
                    bounds_error=False, fill_value="extrapolate")
    vs = np.arange(vmin, vmax + vstep, vstep)

    ys = v_interp(vs)
    
    fig, ax = plt.subplots(1,2, figsize=(8, 3), width_ratios=(1,3))

    # Use the moment 0 map for plotting, otherwise just take a nansum of the cube
    if cubespa_obj.mom_maps is not None:
        plot_data = cubespa_obj.mom_maps.mom0.data
    else:
        plot_data = np.nansum(cubespa_obj.cube.data, axis=0)

    ax[0].imshow(plot_data, origin="lower", cmap=cmap)

    ax[0].plot(pvd["CORNERS"][0], pvd["CORNERS"][1], lw=1, color=pvd_color, alpha=0.3)

    ax[1].imshow(pvd["PVD"], aspect="auto")
    ax[1].set_yticks(ys)
    ax[1].set_yticklabels(vs)
    ax[1].set_ylim(min(ys), max(ys))

    plt.tight_layout()
    if filename is None:
        plt.show()
        plt.close()
    else:
        plt.savefig(filename)
        plt.close()