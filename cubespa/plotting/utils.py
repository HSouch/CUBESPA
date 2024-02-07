import matplotlib as mpl

from matplotlib.colors import ListedColormap

import numpy as np

def reset_style():
    mpl.rcParams.update(mpl.rcParamsDefault)


def pyplot_style():
    mpl.rc('text', usetex=True)

    mpl.rcParams['xtick.major.size'] = 6
    mpl.rcParams['xtick.major.width'] = 2
    mpl.rcParams['xtick.minor.size'] = 4
    mpl.rcParams['xtick.minor.width'] = 1
    mpl.rcParams['ytick.major.size'] = 6
    mpl.rcParams['ytick.major.width'] = 3
    mpl.rcParams['ytick.minor.size'] = 4
    mpl.rcParams['ytick.minor.width'] = 1
    mpl.rcParams['axes.linewidth'] = 2

    mpl.rc('xtick', labelsize=12)
    mpl.rc('ytick', labelsize=12)

    #mpl.rcParams['text.latex.preamble'] = [r'\boldmath']

    font = {'family' : 'serif',
           'weight': 'bold',
           'size': 12}

    mpl.rc('font', **font)


def lavender_cmap(step_1=175):
    x, y = np.mgrid[:251, :251]

    lavender = [86, 82, 100]
    indigo = [29, 0, 51]
    yellow = [100, 100, 0]

    step_2 = 256 - step_1

    # ##########################################################33

    r, g, b = [], [], []

    r = np.array(np.append(np.linspace(lavender[0], indigo[0], step_1), 
                     np.linspace(indigo[0], yellow[0], step_2))) / 100
    g = np.array(np.append(np.linspace(lavender[1], indigo[1], step_1), 
                     np.linspace(indigo[1], yellow[1], step_2))) / 100
    b = np.array(np.append(np.linspace(lavender[2], indigo[2], step_1), 
                     np.linspace(indigo[2], yellow[2], step_2))) / 100

    a = np.array([1 for i in range(len(r))])

    full_arr = np.asarray([[r[i], g[i], b[i], a[i]] for i in range(0, len(r))])
    # full_arr = transpose(full_arr)

    new_cmap = ListedColormap(full_arr)
    return new_cmap
