from matplotlib import pyplot as plt
import numpy as np

from .. import spectra

def plot_spectra(data, aper):
    plt.imshow(data, origin="lower")
    aper.plot()
    plt.tight_layout()
    plt.show()


def spectra_plot(cubespa_obj, aper, spectrum):
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