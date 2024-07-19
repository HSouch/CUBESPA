from matplotlib import pyplot as plt
import numpy as np

from scipy.interpolate import interp1d

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

    cmap = kwargs.get("cmap", "viridis")
    plot_ticks = kwargs.get("plot_ticks", True)
    limits = kwargs.get("limits", None)

    fig = plt.figure(figsize=(11.5, 5), facecolor="white")
    
    img_ax = fig.add_axes([0.05, 0.05, 0.35, 0.90])
    spec_axes = []
    height = (1 - 2 * ypad - (len(spec_list) - 1) * (ysep)) / len(spec_list)
    for i, n in enumerate(spec_list):
        spec_axes.append(fig.add_axes([0.45, ypad +  i * ysep + i * height, 0.50, height]))

    img_ax.imshow(cubespa_obj.mom_maps.mom1.data, origin="lower", cmap=cmap)
    if limits is not None:
        xmin, xmax, ymin, ymax = limits
        img_ax.set_xlim(xmin, xmax)
        img_ax.set_ylim(ymin, ymax)

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


def spectra_comparison(cubecomp, a1, a2, s1, s2, chan_ranges=None, **kwargs):

    cube1, cube2 = cubecomp.cube1, cubecomp.cube2

    colors = ["blue", "red", "orange", "green", "olive", "indigo", "lime", "pink", "orangered", "slateblue"]

    ypad, ysep = 0.02, 0.005

    if cubecomp.cube1.velocities is None:
        v1, v2 = np.arange(len(cubecomp.cube1.data)), np.arange(len(cubecomp.cube2.data))
        xlabel = "Channel"
    else:
        v1, v2 = cubecomp.cube1.velocities, cubecomp.cube2.velocities

        xlabel = "Velocity [km/s]  (Matched to Vsys)"

    shape1, shape2 = cube1.mom_maps.mom1.data.shape, cube2.mom_maps.mom1.data.shape
    ys1, xs1 = np.mgrid[:shape1[0], :shape1[1]]
    ys2, xs2 = np.mgrid[:shape2[0], :shape2[1]]

    cmap = kwargs.get("cmap", "viridis")
    plot_ticks = kwargs.get("plot_ticks", True)

    align = kwargs.get("align", None)
    outname = kwargs.get("outname", None)
    limits = kwargs.get("limits", None)

    chan_ranges = [None for _ in range(len(s1))] if chan_ranges is None else chan_ranges


    fig = plt.figure(figsize=(16, 7), facecolor="white")
    
    img_axes = [fig.add_axes([0.05, 0.05, 0.30, 0.90]), 
                fig.add_axes([0.65, 0.05, 0.30, 0.90])]



    spec_axes = []
    height = (1 - 2 * ypad - (len(s1) - 1) * (ysep)) / len(s1)
    for i, n in enumerate(s1):
        spec_axes.append(fig.add_axes([0.39, 0.01 + ypad +  i * ysep + i * height, 0.24, height]))

    if limits is not None:
        xmin, xmax, ymin, ymax = limits
        img_axes[0].set_xlim(xmin, xmax)
        img_axes[0].set_ylim(ymin, ymax)
        img_axes[1].set_xlim(xmin, xmax)
        img_axes[1].set_ylim(ymin, ymax)
    else:
        if cube1.limits is not None:
            xmin, xmax, ymin, ymax = cube1.limits
            img_axes[0].set_xlim(xmin, xmax)
            img_axes[0].set_ylim(ymin, ymax)
        if cube2.limits is not None:
            xmin, xmax, ymin, ymax = cube2.limits
            img_axes[1].set_xlim(xmin, xmax)
            img_axes[1].set_ylim(ymin, ymax)

    levels = np.linspace(-100, 100, 10)

    im1 = img_axes[0].imshow(cube1.mom_maps.mom1.data, origin="lower", cmap=cmap)
    cont1 = img_axes[0].contour(xs1, ys1, cube1.mom_maps.mom1.data - cube1.vsys, levels=levels, colors="black", alpha=0.5)
    im2 = img_axes[1].imshow(cube2.mom_maps.mom1.data, origin="lower", cmap=cmap)
    cont2 = img_axes[1].contour(xs2, ys2, cube2.mom_maps.mom1.data - cube2.vsys, levels=levels, colors="black", alpha=0.5)

    cbar1 = plt.colorbar(mappable=im1, ax=img_axes[0], location="top", fraction=0.05, label="Relative Velocity [km/s]")
    cbar2 = plt.colorbar(mappable=im2, ax=img_axes[1], location="top", fraction=0.05, label="LSRK Velocity [km/s]")


    for i in range(len(s1)):
        spec_axes[i].set_xlim(np.min(v1), np.max(v1))

        vel_interp = interp1d(np.arange(len(v1)), v1, bounds_error=False, fill_value="extrapolate")

        spec_axes[i].stairs(s1[i][:-1], v1, color=colors[i], lw=1, alpha=0.8, zorder=5)
        spec_axes[i].stairs(s2[i][:-1], v2, color="black", lw=1, alpha=0.8)

        if chan_ranges[i] is not None:
            cmin, cmax = chan_ranges[i]
            spec_axes[i].axvline(vel_interp(cmin), color="black", ls="dashed")
            spec_axes[i].axvline(vel_interp(cmax), color="black", ls="dashed")

            snr = spectra.calc_snr(s1[i], cmin, cmax)
            
            spec_axes[i].text(min(v1) + 20, np.nanmax(s1[i]) * 0.8, f'SNR = {np.round(snr, 1)}')


        spec_med, spec_std = spectra.analyze_spectra(s1[i])
        spec_axes[i].axhline(spec_med, color=colors[i], alpha=0.3)
        spec_axes[i].fill_between([np.min(v1), np.max(v1)], spec_med - spec_std, spec_med + spec_std, color=colors[i], alpha=0.2)


        a1[i].plot(ax=img_axes[0], color=colors[i], lw=4)
        a2[i].plot(ax=img_axes[1], color=colors[i], lw=4)



    for i in range(1, len(spec_axes)):
        spec_axes[i].set_xticks([])

    spec_axes[0].set_xlabel(xlabel)

    if not plot_ticks:
        for ax in img_axes:
            ax.set_xticks([])
            ax.set_yticks([])
    
    if outname is not None:
        plt.savefig(outname, dpi=150)
    else:
        plt.show()
