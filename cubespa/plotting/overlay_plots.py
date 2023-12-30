import numpy as np

from .. import data, utils

from matplotlib import pyplot as plt


def overlay_plot(img_obj, overlay_obj, lims=None, levels = None, colors=None, log_img=False, 
                 cmap="Greys", filename=None):
    if lims is None:
        xmin, xmax, ymin, ymax = 0, img_obj.data.shape[1], 0, img_obj.data.shape[0]
    else:
        xmin, xmax, ymin, ymax = lims

    ys, xs = np.mgrid[:img_obj.data.shape[0], :img_obj.data.shape[1]]    

    if log_img:
        plot_data = np.log10(img_obj.data)
    else:
        plot_data = np.copy(img_obj.data)

    figsize = utils.recommended_figsize(img_obj.data)
    plt.figure(figsize=figsize, facecolor="white")
    plt.imshow(plot_data, origin="lower", cmap=cmap)
    plt.contour(xs, ys, overlay_obj.data, levels=levels, colors="black")

    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)

    plt.xticks([])
    plt.yticks([])

    # TODO: clean up the plotting functions for here and the rgb_overlay function
    if filename is None:
        plt.show()
        plt.close()
    else:
        plt.savefig(filename, dpi=150)
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