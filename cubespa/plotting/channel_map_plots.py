from matplotlib import pyplot as plt


from astropy.stats import sigma_clipped_stats

from .. import utils

import numpy as np

def channel_maps(cubespa_obj, n_chan=25, limits = None, **kwargs):
    
    channels = np.linspace(50, 200, n_chan).astype(int)

    data = cubespa_obj.cube.data

    rows, cols = 5, 5

    def channel_maps(cubespa_obj, n_chan=25, limits=None, **kwargs):
        
        channels = np.linspace(50, 200, n_chan).astype(int)
        data = cubespa_obj.cube.data
        rows, cols = 5, 5

        cmap = kwargs.get("cmap", "Greys")
        show_ticks = kwargs.get("show_ticks", True)
        colors = kwargs.get("colors", "black")
        vmin = kwargs.get("vmin", 0)
        vmax = kwargs.get("vmax", 0.02)
        figsize = utils.recommended_figsize(data[0], width=10)
        
        fig, ax = plt.subplots(rows, cols, facecolor="white", figsize=figsize, sharex=True, sharey=True)
        count = 0
        for i in range(rows):
            for j in range(cols):
                
                if limits is not None:
                    xmin, xmax, ymin, ymax = limits
                    ax[i][j].set_xlim(xmin, xmax)
                    ax[i][j].set_ylim(ymin, ymax)
                if not show_ticks:
                    ax[i][j].set_xticks([])
                    ax[i][j].set_yticks([])
                this_data = data[channels[count]]
                data_median, data_std = sigma_clipped_stats(this_data)[1:]
                
                levels = [data_median + i * data_std for i in range(2, 6)]
                ys, xs = np.mgrid[:this_data.shape[0], :this_data.shape[1]]
                
                ax[i][j].imshow(this_data, cmap=cmap, origin="lower", aspect="auto", vmin=vmin, vmax=vmax)
                
                xlims, ylims = ax[i][j].get_xlim(), ax[i][j].get_ylim()
                ax[i][j].text(xlims[0] * 1.1, ylims[1], np.round(cubespa_obj.velocities[channels[count]], 0), 
                              fontsize=12, backgroundcolor="white")
                
                ax[i][j].contour(xs, ys, this_data, colors=colors, levels=levels, linewidths=0.8,)
                
                count += 1
        
        if limits is not None:
            xmin, xmax, ymin, ymax = limits
            plt.xlim(xmin, xmax)
            plt.ylim(ymin, ymax)
