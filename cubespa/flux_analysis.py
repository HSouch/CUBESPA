from photutils.aperture import EllipticalAperture
from .cubespa import CubeSPA
import numpy as np
from matplotlib import pyplot as plt


class FluxAnalysis:
    def __init__(self, c_obj, Rtrunc=10, stellar_ellip=0.0, stellar_pa=0.0, 
                clockwise=True, annulus_width=20):
        self.c_obj = c_obj
        self.center = np.copy(self.c_obj.center)
        self.shape = self.c_obj.mom_maps.mom0.data.shape
        
        self.Rtrunc=Rtrunc
        self.stellar_ellip = stellar_ellip
        self.stellar_pa = stellar_pa

        self.clockwise = clockwise

        self.annulus_width = annulus_width

    
    def tail_angle_vector(self):
        tail_angle = self.c_obj.tail_angle
        return [np.cos(tail_angle), np.sin(tail_angle)]

    def xy_map(self):
        ys, xs = np.indices(self.shape)
        ys = ys - self.center[1]
        xs = xs - self.center[0]
        return xs, ys

    def theta_map(self):
        xs, ys = self.xy_map()
        thetas = (np.arctan2(ys, xs) - self.c_obj.tail_angle) % (2 * np.pi)
        return thetas

    def r_ellip_pa(self, **kwargs):
        xs, ys = self.xy_map()

        pa = kwargs.get("pa", self.stellar_pa)
        ellip = kwargs.get("ellip", self.stellar_ellip)

        A = xs * np.cos(pa) + ys * np.sin(pa)
        B = -xs * np.sin(pa) + ys * np.cos(pa)
    
        return np.sqrt(A ** 2 + (B / (1 - ellip)) ** 2)
    
    
    def distance_along_tail(self):
        xs, ys = self.xy_map()

        tail_vector = self.tail_angle_vector()

        return xs * tail_vector[0] + ys * tail_vector[1]

    
    def distance_across_tail(self):
        xs, ys = self.xy_map()
        tail_vector = self.tail_angle_vector()

        return xs * tail_vector[1] - ys * tail_vector[0]

    def leading_side(self, offset=0.0):
        return self.distance_along_tail() < offset

    def trailing_side(self, offset=0.0):
        return self.distance_along_tail() > offset
    
    def rotating_into_wind(self):
        return self.distance_across_tail() > 0 if self.clockwise \
            else self.distance_across_tail() < 0

    def rotating_with_wind(self):
        return self.distance_across_tail() < 0 if self.clockwise \
            else self.distance_across_tail() > 0

    def quadrants(self):
        q1 = self.leading_side() & self.rotating_into_wind()
        q2 = self.leading_side() & self.rotating_with_wind()
        q3 = self.trailing_side() & self.rotating_with_wind()
        q4 = self.trailing_side() & self.rotating_into_wind()
        
        return q1, q2, q3, q4

    def annulus(self, rmin, rmax):
        return np.logical_and(self.r_ellip_pa() < rmax, self.r_ellip_pa() > rmin)

    def disk_region(self):
        return self.annulus(0, self.Rtrunc)

    def outer_annulus(self, **kwargs):
        return self.annulus(
            self.Rtrunc - kwargs.get("width", self.annulus_width),
            self.Rtrunc
        )

    def get_flux(self, mask, plot=False):
        mom0 = np.copy(self.c_obj.mom_maps.mom0.data)
        mom0[mask == 0] = 0

        if plot:
            plt.imshow(mom0, origin="lower")
            plt.show()

        return np.nansum(mom0)


    # ##############################################
    # ######### PLOTTING ROUTINES ##################
    # ##############################################

    def plot_cent(self, axis, **kwargs):
        cent = self.c_obj.center
        axis.scatter(cent[0], cent[1], **kwargs)

    
    def plot_disk(self, axis, **kwargs):
        aper = EllipticalAperture((self.center[0], self.center[1]), 
        self.Rtrunc, self.Rtrunc * (1 - self.stellar_ellip), self.stellar_pa,
            )
        
        aper.plot(ax=axis, **kwargs)


    def plot_tail_angle(self, axis, length=50, **kwargs):
        cent_x, cent_y = self.c_obj.center

        tail_vec = self.tail_angle_vector()

        axis.plot(
            [cent_x, cent_x + length * tail_vec[0]],
            [cent_y, cent_y + length * tail_vec[1]],
            **kwargs)

    def plot_distances_along_map(self, axis, contour=True, **kwargs):
        im = axis.imshow(self.distance_along_tail(), vmin=-100, vmax=100, cmap="RdBu_r")
        if contour:
            axis.contour(self.distance_along_tail(), levels=[0], colors="white", linestyles="dashed")

    def plot_distances_across_map(self, axis, contour=True, **kwargs):
        axis.imshow(self.distance_across_tail(), vmin=-100, vmax=100, cmap="RdBu_r")
        if contour:
            axis.contour(self.distance_across_tail(), levels=[0], colors="white", linestyles="dashed")

    def plot_rmap(self, axis, contour=True, **kwargs):
        axis.imshow(self.r_ellip_pa(**kwargs), origin="lower", cmap="inferno")
        if contour:
            axis.contour(self.r_ellip_pa(**kwargs), levels=[20, 50, 80, 110], colors="white", linestyles="dashed")

    def plot_leading_side(self, axis, **kwargs):
        axis.imshow(self.leading_side(), cmap=kwargs.pop("cmap", "Greys"), **kwargs)

    def overplot_mom0(self, axis, **kwargs):
        mom0 = np.copy(self.c_obj.mom_maps.mom0.data)
        mom0[mom0 <=0] = np.nan
        axis.imshow(np.log10(mom0), cmap=kwargs.pop("cmap", "inferno_r"), 
            vmin=-1, vmax=1.5, origin="lower", **kwargs)


class FluxMeasurements:
    def __init__(self, total=None, total_disk=None,
                 lead=None, trail=None,
                lead_disk=None, trail_disk=None,
                lead_annulus=None, trail_annulus=None,
                ):
        self.total = total
        self.total_disk = total_disk
        self.lead = lead
        self.trail = trail
        self.lead_disk = lead_disk
        self.trail_disk = trail_disk
        self.lead_annulus = lead_annulus
        self.trail_annulus = trail_annulus
    

    def results_dict(self):
        res = {}
        for name in [
            "TOTAL", "TOTAL_DISK",
            "LEAD", "TRAIL", 
            "LEAD_DISK", "TRAIL_DISK",
            "LEAD_ANNULUS", "TRAIL_ANNULUS"]:
            res[name] = float(getattr(self, name.lower()))
        return res
    
    @staticmethod
    def from_fluxanalysis(fa, debug_plot=False):
        return FluxMeasurements(
            total = fa.get_flux(mask = np.ones(fa.shape), plot=debug_plot),
            total_disk = fa.get_flux(mask=fa.disk_region(), plot=debug_plot),

            lead = fa.get_flux(mask=fa.leading_side(), plot=debug_plot),
            trail = fa.get_flux(mask=fa.trailing_side(), plot=debug_plot),

            lead_disk=fa.get_flux(mask=fa.leading_side() & fa.disk_region(), plot=debug_plot),
            trail_disk=fa.get_flux(mask=fa.trailing_side() & fa.disk_region(), plot=debug_plot),

            lead_annulus=fa.get_flux(mask=fa.leading_side() & fa.outer_annulus(), plot=debug_plot),
            trail_annulus=fa.get_flux(mask=fa.trailing_side() & fa.outer_annulus(), plot=debug_plot),
        )



def fluxanalysis_debug_plot(c: CubeSPA, isophote_info, 
    outdir="plots/Debug/FluxAnalysisDebug/", close=True):
    
    name = c.name
    med_ellip, med_pa = isophote_info["ELLIP_MED"], isophote_info["PA_MED"]

    fa = FluxAnalysis(c, stellar_ellip=med_ellip, stellar_pa=med_pa)
    cent_x, cent_y = fa.c_obj.center

    fig, ax = plt.subplots(2, 4, figsize=(12, 8))


    fa.plot_distances_along_map(ax[0, 0])
    fa.plot_distances_across_map(ax[0, 1])
    fa.plot_rmap(ax[0, 2])
    fa.plot_leading_side(ax[0, 3])


    for axis in ax.flatten():
        fa.overplot_mom0(axis)
        fa.c_obj.apply_limits(axis, pad=10)

        axis.scatter(cent_x, cent_y, color="black", marker="x", s=50)
        fa.plot_tail_angle(axis, color="black", lw=2)

        axis.set_xticks([])
        axis.set_yticks([])

    ax[0, 0].set_title("Distances Along Tail")
    ax[0, 1].set_title("Distances Across Tail")
    ax[0, 2].set_title("R(x,y) - PA/Ellip Adj.")
    ax[0, 3].set_title("Leading / Trailing Sides")

    plt.tight_layout()
    plt.savefig(f"{outdir}{name}.png", bbox_inches="tight")
    if close:
        plt.close()
