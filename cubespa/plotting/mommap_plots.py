from matplotlib import pyplot as plt

import numpy as np

from .. import utils


def moment_map_plot(cubespa_obj, outname = None, use_limits=True, **kwargs):
    """ Generate moment map plots.

    Args:
        cubespa_obj (cubespa.CubeSPA): The input CubeSPA object, with valid moment maps loaded.
        filename (str, optional): Output filename, in which the plot is just shown instead of saved. Defaults to None.
        use_limits (bool, optional): Whether or not to use limits from the CubeSPA object.
            
            It is a good idea to set to False for cutout objects, as their limits will be relative to the initial CubeSPA object,
            and their desired limits will be the entire array. Defaults to True.
    """
    # if "kwargs" in kwargs.keys():
    #     kwargs.update(kwargs["kwargs"])

    mom0, mom1, mom2 = cubespa_obj.mom_maps.mom0.data, cubespa_obj.mom_maps.mom1.data, cubespa_obj.mom_maps.mom2.data, 

    xs, ys = np.mgrid[:mom0.shape[0], :mom0.shape[1]]
    
    v_offset = kwargs.get("vsys", None)

    if v_offset is None:
        v_offset = cubespa_obj.vsys

    mom0_lims = kwargs.get("mom0_lims", [-3, 1.1])
    figsize = kwargs.get("figsize", (12, 7))
    top = kwargs.get("top", 0.95)
    name = kwargs.get("name", "")
    line = kwargs.get("line", "CO(2-1)")
    fontsize=kwargs.get("fontsize", 16)
    snr_levels = kwargs.get("snr_levels", [3, 5, 10, 20, 30, 50])
    
    beam_area = cubespa_obj.beam_area
    mom0_units = r'$\log_{10}($ I [Jy beam$^{-1}$ km s$^{-1}$ ])'

    wcs = kwargs.get("wcs", None)
    if wcs is None:
        fig, ax = plt.subplots(1,3, figsize=figsize, sharey=True, facecolor="white")
    else:
        fig = plt.figure(figsize=figsize, facecolor="white")
        ax = [fig.add_subplot(131, projection=wcs), 
              fig.add_subplot(132, projection=wcs), 
              fig.add_subplot(133, projection=wcs)]
    
    if use_limits:
        xmin, xmax, ymin, ymax = cubespa_obj.limits
        dx, dy = xmax - xmin, ymax - ymin
        for i in range(3):
            ax[i].set_xlim(xmin, xmax)

        for axis in ax:
            axis.set_ylim(ymin, ymax + 0.1 * dy)
    else:
        xmin, ymin = 0, 0
        ymax, xmax = mom0.shape             #TODO check if this is right (should it be xmax, ymax)
        dy, dx = mom0.shape

    mom0_ax = ax[0].imshow(np.log10(mom0), origin="lower", cmap="Greys", vmin=mom0_lims[0], vmax=mom0_lims[1])
    mom0_cb = plt.colorbar(mappable=mom0_ax, ax=ax[0], location="top", label=mom0_units)
    mom0_cb.set_label(mom0_units, fontsize=fontsize)
    
    mom0_levels = [cubespa_obj.cube_rms * np.sqrt(len(cubespa_obj.velocities)) * i for i in snr_levels]    
    ax[0].contour(ys, xs, mom0, origin="lower", colors="black", 
                levels=mom0_levels, linewidths=0.75)
    

    mom1_ax = ax[1].imshow(mom1, origin="lower", cmap="rainbow", 
                           vmin= kwargs.get("vmin", v_offset -50.), vmax=kwargs.get("vmax", v_offset + 50.))
    # ax[1].contour(ys, xs, mom0, origin="lower", colors="white", 
    #               levels=[0.08, 0.1, 0.5, 1, 2, 2.5], linewidths=1.5)
    ax[1].contour(ys, xs, mom1 - v_offset, origin="lower", colors="black", linewidths=1.5, 
                levels = kwargs.get("mom1_vcontours", np.linspace(-100, 100, 10)))
    mom1_cb = plt.colorbar(mappable=mom1_ax, ax=ax[1], location="top", label="Velocity [km/s]")
    mom1_cb.set_label("Velocity [km s$^{-1}$]", fontsize=fontsize)

    mom2_ax = ax[2].imshow(mom2, origin="lower", cmap="magma", vmin=0, vmax=50)
    mom2_cb = plt.colorbar(mappable=mom2_ax, ax=ax[2], location="top")
    mom2_cb.set_label("Velocity Dispersion [km s$^{-1}$]", fontsize=fontsize)
    
    if kwargs.get("write_info", True):
        ax[0].text(xmin + dx / 20, ymax + dy/20, f'{name} - {line}', fontsize=fontsize)
        ax[0].text(xmin + dx / 20, ymax, f'Moment 0 Intensity', fontsize=fontsize)
        if kwargs.get("plot_total_flux", False):
            ax[0].text(xmin + dx / 20, ymax - dy / 20, 
                    r'$S_{CO}=$ ' + str(np.round(np.nansum(mom0), 3)) + r" Jy km s$^{-1}$ beam$^{-1}$", 
                    fontsize=fontsize)

        ax[1].text(xmin + dx / 20, ymax + dy/20, f'{name} - {line}', fontsize=fontsize)
        ax[1].text(xmin + dx / 20, ymax, f'Moment 1 Velocity', fontsize=fontsize)

        ax[2].text(xmin + dx / 20, ymax + dy/20, f'{name} - {line}', fontsize=fontsize)
        ax[2].text(xmin + dx / 20, ymax, f'Moment 2 Velocity Dispersion', fontsize=fontsize)

    if wcs is None:
        if kwargs.get("hide_ticks", True):
            for axis in ax:
                axis.set_xticks([])
                axis.set_yticks([])
    else:
        for axis in ax[1:]:
            axis.coords[1].set_ticklabel_visible(False)
            axis.coords[1].set_axislabel("")
        for axis in ax:
            axis.coords[0].set_axislabel("Right Ascension", fontsize=kwargs.get("label_fontsize", 16))
        ax[0].coords[1].set_axislabel("Declination", fontsize=kwargs.get("label_fontsize", 16))            

    # PLot beam if possible
    if cubespa_obj.beam_pix is not None:
        bmaj, bmin, bpa = cubespa_obj.beam_pix
        beam_xs, beam_ys = utils.ellipse_coords(xmin + dx/12, ymin + dy/12, 
                                                bmaj / 2, bmin/2, np.deg2rad(bpa), num_points=30)
        for axis in ax:
            axis.plot(beam_xs, beam_ys, color="black", lw=1)

    plt.tight_layout()
    plt.subplots_adjust(top=top, wspace=kwargs.get("wspace", 0.02),
                        left=kwargs.get("left", 0.1), right=kwargs.get("right", 0.99),
                        bottom=kwargs.get("bottom", 0.05))

    if kwargs.get("return_fig", False):
        return fig, ax

    if outname is None:
        plt.show()
        plt.close()
    else:
        plt.savefig(outname, dpi=150, bbox_inches="tight")
        plt.close()
    
