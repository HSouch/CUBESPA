from matplotlib import pyplot as plt
import numpy as np

from .. import spectra, utils

def plot_spectra(data, aper):
    plt.imshow(data, origin="lower")
    aper.plot()
    plt.tight_layout()
    plt.show()


def spectra_plot(cubespa_obj, aper, spectrum):
    """ Create a plot showing both the image with overlaid spectra, as well as the spectrum with
        RMS levels shown.

    Args:
        cubespa_obj (cubespa.CubeSPA): CubeSPA object.
        aper (photutils.aperture): Input aperture, generated using cubespa.spectra
        spectrum (_type_): _description_
    """
    fig, ax = plt.subplots(1,2, figsize=(10, 4), width_ratios=(2, 3))

    ax[0].imshow(cubespa_obj.mom_maps.mom0.data, origin="lower")
    aper.plot(ax=ax[0])

    if cubespa_obj.velocities is None:
        xlabel = "Channel"
    else:
        xlabel = "Velocity"

    velocities = np.arange(len(spectrum)) if cubespa_obj.velocities is None else cubespa_obj.velocities

    spec_med, spec_std = spectra.analyze_spectra(spectrum)
    
    ax[1].stairs(spectrum[:-1], velocities, color="black")
    ax[1].axhline(spec_med, color="grey")
    ax[1].fill_between(velocities, spec_med - spec_std, spec_med + spec_std, color="black", alpha=0.1)
    ax[1].fill_between(velocities, spec_med - 2 * spec_std, spec_med + 2 * spec_std, color="black", alpha=0.1)
    ax[1].fill_between(velocities, spec_med - 3 * spec_std, spec_med + 3 * spec_std, color="black", alpha=0.1)
    
    ax[1].set_xlabel(xlabel)
    ax[1].set_ylabel("Intensity")

    ax[1].set_xlim(np.min(velocities), np.max(velocities))

    plt.show()


def multispec_plot(cubespa_obj, aper_list, spec_list, **kwargs):

    colors = ["blue", "red", "orange", "black", "green", "olive"]

    ypad, ysep = 0.02, 0.005

    if cubespa_obj.velocities is None:
        velocities = np.arange(len(cubespa_obj.data))
        xlabel = "Channel"
    else:
        velocities = cubespa_obj.velocities
        xlabel = "Velocity [km/s]"

    cmap = utils.check_kwarg("cmap", "viridis", kwargs)
    plot_ticks = utils.check_kwarg("plot_ticks", True, kwargs)


    fig = plt.figure(figsize=(11.5, 5), facecolor="white")
    
    img_ax = fig.add_axes([0.05, 0.05, 0.35, 0.90])
    spec_axes = []
    height = (1 - 2 * ypad - (len(spec_list) - 1) * (ysep)) / len(spec_list)
    for i, n in enumerate(spec_list):
        spec_axes.append(fig.add_axes([0.45, ypad +  i * ysep + i * height, 0.50, height]))

    img_ax.imshow(cubespa_obj.mom_maps.mom1.data, origin="lower", cmap=cmap)


    for i in range(len(spec_list)):
        spec_axes[i].stairs(spec_list[i][:-1], velocities, color=colors[i], lw=1)
        spec_axes[i].set_xlim(np.min(velocities), np.max(velocities))

        spec_med, spec_std = spectra.analyze_spectra(spec_list[i])
        
        spec_axes[i].axhline(spec_med, color=colors[i], alpha=0.3)
        for j in range(1, 4):
            spec_axes[i].fill_between(velocities, spec_med - j* spec_std, spec_med + j* spec_std, color=colors[i], alpha=0.15)

        aper_list[i].plot(ax=img_ax, color=colors[i], lw=2)

    for i in range(1, len(spec_axes)):
        spec_axes[i].set_xticks([])

    spec_axes[0].set_xlabel(xlabel)

    if not plot_ticks:
        img_ax.set_xticks([])
        img_ax.set_yticks([])
    
    plt.show()