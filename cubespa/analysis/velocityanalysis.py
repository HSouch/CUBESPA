import numpy as np
from astropy.io import fits

from ..utils import remap, norm

from astropy.cosmology import Planck18 as cosmo
from astropy.convolution import convolve, Gaussian2DKernel

from matplotlib import pyplot as plt


class ViewingGeometry:
    def __init__(self, d_obs, tail_angle_deg=78, north_vec=[0, 0, 1]):
        self.tail_angle_deg = tail_angle_deg
        self.tail_angle = np.radians(self.tail_angle_deg)
        self.tail_vector = np.array([0.0, np.sin(self.tail_angle), np.cos(self.tail_angle)])

        self.d_obs = -norm(np.array(d_obs))

    @staticmethod
    def from_header(header, tail_angle_deg=78, north_vec = [0, 0, 1]):
        d_obs = -np.array([header[f"LOS_{n}"] for n in ["X", "Y", "Z"]])
        return ViewingGeometry(d_obs, tail_angle_deg=tail_angle_deg, north_vec=north_vec)


class MinorAxisFit:
    def __init__(self, rs, ibins, means, stds, head = None):
        self.rs = rs
        self.ibins = ibins
        self.means = means
        self.stds = stds

        self.head = head

    def to_cartesian(self, x0=0, y0=0, mirror=True):
        
        if mirror:
            rs_full = np.concatenate([self.rs, self.rs])
            means_full = np.concatenate([self.means, self.means + np.pi])
        else:
            rs_full = np.copy(self.rs)
            means_full = np.copy(self.means)

        xs = rs_full * np.cos(means_full) + x0
        ys = rs_full * np.sin(means_full) + y0
        
        return xs, ys
    
    def bounds_cartesian(self, x0=0, y0=0, mirror=True):

        lower_bounds = self.means - self.stds
        upper_bounds = self.means + self.stds

        if mirror:
            rs_full = np.concatenate([self.rs, self.rs])
            lower_bounds_full = np.concatenate([lower_bounds, lower_bounds + np.pi])
            upper_bounds_full = np.concatenate([upper_bounds, upper_bounds + np.pi])
        else:
            rs_full = np.copy(self.rs)
            lower_bounds_full = np.copy(lower_bounds)
            upper_bounds_full = np.copy(upper_bounds)
            
        lower_xs = rs_full * np.cos(lower_bounds_full) + x0
        lower_ys = rs_full * np.sin(lower_bounds_full) + y0

        upper_xs = rs_full * np.cos(upper_bounds_full) + x0
        upper_ys = rs_full * np.sin(upper_bounds_full) + y0

        return (lower_xs, lower_ys), (upper_xs, upper_ys)
    
    @staticmethod
    def fit_minor_axis(mom0, mom1, vmax=10, imin=1e2, Rmin=1, Rmax=20, bin_width=2, extent=20,
        tail_angle_deg = 90, which="both", cent=None, head=None):

        if cent is None:
            cent = np.array(mom0.shape) / 2
        
        ys, xs = np.indices(mom0.shape)
        ys = ys - cent[1]
        xs = xs - cent[0]

        tail_angle = np.radians(tail_angle_deg)

        d_along = distance_along_tail(mom0.shape, cent, tail_angle)

        ys, xs = np.mgrid[0:mom1.shape[0], 0:mom1.shape[1]]
        ys = remap(ys, 0, mom1.shape[0], -extent, extent)
        xs = remap(xs, 0, mom1.shape[1], -extent, extent)

        rs = np.sqrt(xs ** 2 + ys ** 2)
        thetas = np.arctan2(ys, xs) % (2 * np.pi)

        mom1_mask = np.logical_and(np.abs(mom1) < vmax, np.abs(mom0) > imin)
        mom1_mask = np.logical_and(mom1_mask, rs > Rmin)

        if which == "lead" or which == "leading":
            mom1_mask = np.logical_and(mom1_mask, leading_side(d_along))
        if which == "trail" or which == "trailing":
            mom1_mask = np.logical_and(mom1_mask, trailing_side(d_along))

        zs = mom0[mom1_mask]
        rs_masked = rs[mom1_mask]
        thetas_masked = thetas[mom1_mask]

        r_bins = np.arange(Rmin, Rmax + bin_width, bin_width)

        res = []

        for i in range(len(r_bins) - 1):
            rmin, rmax = r_bins[i], r_bins[i + 1]
            rmean = (rmin + rmax) / 2

            bin_mask = (rs_masked > rmin) & (rs_masked <= rmax)

            thetas_submask = thetas_masked[bin_mask]
            z_submask = zs[bin_mask]

            ibin = np.nansum(z_submask) / np.nansum(bin_mask)
            if np.isnan(ibin):
                ibin = 0.0

            mean_theta = np.nanmedian(thetas_submask)
            std_theta = np.nanstd(thetas_submask)

            res.append((rmean, ibin, mean_theta, std_theta))
        
        rmeans, ibins, tmeans, tstds = np.array(res).T

        return MinorAxisFit(rs=rmeans, ibins=ibins, means=tmeans, stds=tstds, head=head)


def smooth_to_alma_beam(mom, extent, z=0.0231, beam=(1.1, 0.8, np.pi/4), plot=False):

    # Convert to angular size
    kpc_per_arcsec = 1/ cosmo.arcsec_per_kpc_proper(z).value

    bmaj, bmin, bpa = beam

    pix_per_kpc = mom.shape[0] / (extent[1] - extent[0])

    x_stddev = bmaj / (2 * np.sqrt(2 * np.log(2))) * kpc_per_arcsec * pix_per_kpc
    y_stddev = bmin / (2 * np.sqrt(2 * np.log(2))) * kpc_per_arcsec * pix_per_kpc

    kernel = Gaussian2DKernel(x_stddev=x_stddev, y_stddev=y_stddev, theta=bpa)

    convolved = convolve(mom, kernel)

    return convolved


def xy_map(shape, cent):
    ys, xs = np.indices(shape)
    ys = ys - cent[1]
    xs = xs - cent[0]

    return xs, ys


def tail_angle_vector(angle_rad):
    return np.array([np.cos(angle_rad), np.sin(angle_rad)])



def distance_along_tail(shape, cent, tail_angle_rad):
    xs, ys = xy_map(shape, cent)

    tail_vector = tail_angle_vector(tail_angle_rad)

    return xs * tail_vector[0] + ys * tail_vector[1]


def leading_side(distances_along_tail, offset=0.0):
        return distances_along_tail < offset


def trailing_side(distances_along_tail, offset=0.0):
    return distances_along_tail > offset
