from matplotlib import pyplot as plt

import numpy as np

from .. import utils


def moment_map_plot(cubespa_obj, filename = None, use_limits=True, **kwargs):
    """ Generate moment map plots.

    Args:
        cubespa_obj (cubespa.CubeSPA): The input CubeSPA object, with valid moment maps loaded.
        filename (str, optional): Output filename, in which the plot is just shown instead of saved. Defaults to None.
        use_limits (bool, optional): Whether or not to use limits from the CubeSPA object.
            
            It is a good idea to set to False for cutout objects, as their limits will be relative to the initial CubeSPA object,
            and their desired limits will be the entire array. Defaults to True.
    """
    if "kwargs" in kwargs.keys():
        kwargs.update(kwargs["kwargs"])


    xmin, xmax, ymin, ymax = cubespa_obj.limits
    dx, dy = xmax - xmin, ymax - ymin

    mom0, mom1, mom2 = cubespa_obj.mom_maps.mom0.data, cubespa_obj.mom_maps.mom1.data, cubespa_obj.mom_maps.mom2.data, 

    xs, ys = np.mgrid[:mom0.shape[0], :mom0.shape[1]]
    
    v_offset = utils.check_kwarg("vsys", 0, kwargs)
    mom0_lims = utils.check_kwarg("mom0_lims", [-3, 1.1], kwargs)
    divide_beam = utils.check_kwarg("divide_beam", False, kwargs)

    
    beam_area = cubespa_obj.beam_area if divide_beam else 1
    mom0_units = r'$\log_{10}($ I [Jy km/s ])' if divide_beam else r'$\log_{10}($ I [Jy beam$^{-1}$ km/s ])'

    fig, ax = plt.subplots(1,3, figsize=(12,7), sharey=True, facecolor="white")

    if use_limits:
        for i in range(3):
            ax[i].set_xlim(xmin, xmax)

        plt.ylim(ymin, ymax + 0.1 * dy)

    mom0_ax = ax[0].imshow(np.log10(mom0 / beam_area), origin="lower", cmap="Greys", vmin=mom0_lims[0], vmax=mom0_lims[1])
    mom0_cb = plt.colorbar(mappable=mom0_ax, ax=ax[0], location="top", label=mom0_units)

    ax[0].contour(ys, xs, mom0, origin="lower", colors="black", 
                levels=[0.04, 0.08, 0.1, 0.5, 1, 2, 2.5], linewidths=0.75)
    ax[0].contour(ys, xs, np.log10(mom0), origin="lower", colors="black", levels=[-5, -4, -3, -2, -1, 0, 1])
    ax[0].text(xmin + dx / 10, ymax, r'$F_{tot}=$ ' + str(np.round(np.nansum(mom0 / beam_area), 3)))

    mom1_ax = ax[1].imshow(mom1 - v_offset, origin="lower", cmap="rainbow", vmin= - 50, vmax= 50)
    # ax[1].contour(ys, xs, mom0, origin="lower", colors="white", 
    #               levels=[0.08, 0.1, 0.5, 1, 2, 2.5], linewidths=1.5)
    ax[1].contour(ys, xs, mom1 - v_offset, origin="lower", colors="black", linewidths=1.5, 
                levels=np.linspace(-100, 100, 10))
    mom1_cb = plt.colorbar(mappable=mom1_ax, ax=ax[1], location="top", label="Velocity [km/s]")


    mom2_ax = ax[2].imshow(mom2, origin="lower", cmap="magma", vmin=0, vmax=50)
    mom2_cb = plt.colorbar(mappable=mom2_ax, ax=ax[2], location="top", label="Velocity Dispersion [km/s]")

    plt.tight_layout()

    if filename is None:
        plt.show()
        plt.close()
    else:
        plt.savefig(filename, dpi=150)
        plt.close()