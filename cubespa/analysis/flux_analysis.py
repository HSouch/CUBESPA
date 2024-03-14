import numpy as np
from .. import utils

from matplotlib import pyplot as plt

from photutils.aperture import EllipticalAnnulus


def quadrants(q1_angle, offset = 90, clockwise=False):
    mult = -1 if clockwise else 1
    return np.array([(q1_angle + i * 90 * mult + offset) for i in range(4)])


def halfline_points(cent, length, angle):

    p1 = (cent[0] , cent[0] + length * np.cos(angle))
    p2 = (cent[1], cent[1] + length * np.sin(angle))
    
    return (p1, p2)


def line_points(cent, length, angle):
    p1 = (cent[0] - length * np.cos(angle), cent[0] + length * np.cos(angle))
    p2 = (cent[1] - length * np.sin(angle), cent[1] + length * np.sin(angle))

    
    return (p1, p2)


def theta_list(theta_0, clockwise=False, n_samples=90, plot=False, **kwargs):
    outname = utils.check_kwarg("outname", None, kwargs)

    theta_max = theta_0 + 2*np.pi

    if clockwise:
        _ = theta_0
        theta_0, theta_max = theta_max, _
    
    print(theta_0, theta_max)

    angles = np.linspace(theta_0, theta_max, n_samples) % (2 * np.pi)

    if plot:
        thetas_norm = np.linspace(0, 2*np.pi, n_samples)
        xs, ys = np.cos(angles), np.sin(angles)
        plt.scatter(xs, ys, c=np.rad2deg(thetas_norm))
        plt.colorbar()
        plt.tight_layout()
        if outname is not None:
            plt.savefig(outname, dpi=100)
            plt.close()
        else:
            plt.show()
    
    return angles


def theta_range(shape, x0, y0, theta, dtheta):
    ys, xs = np.mgrid[:shape[0], :shape[1]].astype(float)
    xs -= x0
    ys -= y0
    
    thetas = (np.arctan2(ys, xs) ) % (2 * np.pi)

    theta_mask = np.ones(thetas.shape)

    theta_min = (theta - dtheta / 2) % (2 * np.pi)
    theta_max = (theta + dtheta / 2) % (2 * np.pi)
    
    # Adjust for cases where theta_max < theta_min due to wrapping around 2pi
    if theta_max < theta_min:
        # In this case, we need to handle the wrap around.
        theta_mask = np.logical_or(thetas >= theta_min, thetas <= theta_max)
    else:
        # Regular case where the range does not wrap around
        theta_mask = np.logical_and(thetas >= theta_min, thetas <= theta_max)

    return theta_mask.astype(int)


def process_azimuth(data, radius, pa, eps, thetas, stellar_x=0, stellar_y=0, annulus_width=50, trim_nans=False):

    
    ell_annulus = EllipticalAnnulus((stellar_x, stellar_y), radius, radius + annulus_width, 
                                    (radius + annulus_width) * (1-eps), theta = pa)
    mask = ell_annulus.to_mask(method='center').to_image(shape=data.shape)
    
    data_init_mask = np.copy(data)
    data_init_mask[mask == 0] = np.nan

    fluxes, flux_stds = [], []

    for i, theta in enumerate(thetas[:]):
    
        theta_mask = theta_range(data.shape, stellar_x, stellar_y, theta, 5 * np.pi / 180)
    
        this_data = np.copy(data_init_mask)
        this_data[theta_mask == 0] = np.nan
        try:
            med, std = np.nanmedian(this_data), np.nanstd(this_data)
        except RuntimeWarning:
            med, std = np.nan, np.nan
            
        fluxes.append(med)
        flux_stds.append(std)

    if trim_nans:
        fluxes = np.array(fluxes)
        fluxes[np.isnan(fluxes)] = 0

        flux_stds = np.array(flux_stds)
        flux_stds[np.isnan(flux_stds)] = 0
    
    return (np.array(fluxes), np.array(flux_stds), ell_annulus)


def segment_sum(cubespa_obj, data, radius, annulus_width, theta_0, plot=False, clockwise=False):
    stellar_x, stellar_y = cubespa_obj.center
    
    eps, pa = cubespa_obj.eps, cubespa_obj.position_angle

    # print(stellar_x, stellar_y, eps, pa)
    ell_annulus = EllipticalAnnulus((stellar_x, stellar_y), 
                                    radius, radius + annulus_width, 
                                    (radius + annulus_width) * (1-eps), 
                                    theta = pa)
    mask = ell_annulus.to_mask(method='center').to_image(shape=data.shape)
    
    data_init_mask = np.copy(data)
    data_init_mask[mask == 0] = np.nan

    theta_mask = theta_range(data.shape, stellar_x, stellar_y, np.deg2rad(theta_0) + np.pi/4, np.pi / 2)

    data_init_mask[theta_mask == 0] = np.nan
    pixel_sum = np.nansum(data_init_mask)

    if plot:
        plt.imshow(data, origin="lower", cmap="Greys")
        plt.imshow(data_init_mask, cmap="viridis", origin="lower", alpha=0.7)
        plt.imshow(theta_mask, alpha=0.2, cmap="Greys", origin="lower")
        ell_annulus.plot(color="black")
        plt.title(f'{radius:.1f} {theta_0}\n{pixel_sum:.2f}', fontsize=10)
        plt.show()

    return pixel_sum