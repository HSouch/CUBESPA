import numpy as np

from .. import data, utils

from matplotlib import pyplot as plt


def overlay_plot(img_obj, overlay_obj, lims=None,  **kwargs):
    if lims is None:
        xmin, xmax, ymin, ymax = 0, img_obj.data.shape[1], 0, img_obj.data.shape[0]
    else:
        xmin, xmax, ymin, ymax = lims

    ys, xs = np.mgrid[:img_obj.data.shape[0], :img_obj.data.shape[1]]    

    log_img = utils.check_kwarg("log_img", False, kwargs)
    colors = utils.check_kwarg("colors", "black", kwargs)
    levels = utils.check_kwarg("levels", None, kwargs)
    cmap = utils.check_kwarg("cmap", "Greys", kwargs)
    outname = utils.check_kwarg("outname", None, kwargs)
    img_levels = utils.check_kwarg("img_levels", [-3, 1] if log_img else [0, 0.1], kwargs)

    if log_img:
        plot_data = np.log10(img_obj.data)
    else:
        plot_data = np.copy(img_obj.data)

    figsize = utils.recommended_figsize(img_obj.data)
    plt.figure(figsize=figsize, facecolor="white")
    plt.imshow(plot_data, origin="lower", cmap=cmap, vmin=img_levels[0], vmax=img_levels[1])
    plt.contour(xs, ys, overlay_obj.data, levels=levels, colors=colors)

    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)

    plt.xticks([])
    plt.yticks([])

    # TODO: clean up the plotting functions for here and the rgb_overlay function
    if outname is None:
        plt.show()
        plt.close()
    else:
        plt.savefig(outname, dpi=150)
        plt.close()


def rgb_overlay(rgb_img, overlay_obj, lims=None, levels=None, colors=None, filename=None, **kwargs):
    # TODO: Remove the messiness and add it to kwargs
    colors = "black" if colors is None else colors

    if type(rgb_img) is data.DataSet:
        rgb_img = rgb_img.data

    if lims is None:
        xmin, xmax, ymin, ymax = 0, rgb_img.data.shape[1], 0, rgb_img.data.shape[0]
        figsize = utils.recommended_figsize(overlay_obj.data)
    else:
        xmin, xmax, ymin, ymax = lims
        figsize = utils.recommended_figsize(np.zeros([ymax - ymin, xmax-xmin]))

    ys, xs = np.mgrid[:rgb_img.data.shape[0], :rgb_img.data.shape[1]]    

    plt.figure(figsize=figsize, facecolor="white")
    plt.imshow(rgb_img, origin="lower")
    plt.contour(xs, ys, overlay_obj.data, levels=levels, colors=colors)

    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)

    plt.xticks([])
    plt.yticks([])

    plt.tight_layout()

    if filename is None:
        plt.show()
        plt.close()
    else:
        plt.savefig(filename, dpi=150)
        plt.close()


def plot_psf_overlay(cubespa_obj, psf_conv=None, x0=0, y0 = 0, **kwargs):

    figsize = utils.recommended_figsize(cubespa_obj.mom_maps.mom0.data)
    plt.figure(figsize=figsize, facecolor="white")

    show_ticks = utils.check_kwarg("show_xticks", True, kwargs)
    levels = np.array(utils.check_kwarg("levels", [0.02, 0.05, 0.1], kwargs))

    xmin, xmax, ymin, ymax = cubespa_obj.limits
    dx, dy = xmax-xmin, ymax-ymin

    mom0 = cubespa_obj.mom_maps.mom0.data
    psf = np.nanmean(cubespa_obj.psf.data, axis=0)


    ys, xs = np.mgrid[:psf.shape[0], :psf.shape[1]]
    xs -= x0
    ys -= y0


    plt.imshow(np.log10(mom0), origin="lower", cmap="Greys")
    plt.text(xmin + dx / 10, ymax - dy/18, f'PSF Contours Levels = {levels * 100} %', fontsize=15, color="Red")


    plt.contour(xs, ys, psf, levels=levels, colors="red", alpha=0.5)

    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)

    if not show_ticks:
        plt.xticks([])
        plt.yticks([])

    plt.tight_layout()
    plt.show()
