import os

import numpy as np

from photutils.aperture import EllipticalAperture

from scipy.signal import savgol_filter
from scipy.interpolate import interp1d

from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation


def R_ellip(r, ellip, theta):
    axial_ratio = 1 - ellip
    b = r * axial_ratio
    eccentricity = np.sqrt(1 - axial_ratio**2)

    return b / np.sqrt(1 - eccentricity**2 * np.cos(theta)**2)


def deviation_from_major_axis(theta, pa):
    """
    Calculate the deviation from the major axis in radians.
    """
    delta = (theta - pa) % np.pi 

    return np.min([delta, np.pi - delta], axis=0)  # Return the smallest angle to the major axis


class Profile:
    def __init__(self, rs, fluxes, rmax=None, **kwargs):
        self.rs = rs
        self.fluxes = fluxes
        self.rmax = rmax

        nanmask = np.isnan(fluxes)
        if np.all(nanmask):
            raise ValueError("All fluxes are NaN. Cannot create profile.")
        self.rs = rs[~nanmask]
        self.fluxes = fluxes[~nanmask]

        if rmax is not None:
            self.rs = self.rs[self.rs <= rmax]
            self.fluxes = self.fluxes[:len(self.rs)]

        self.prof_interp = interp1d(rs, fluxes, bounds_error=False, fill_value=np.nan)

        self.theta = kwargs.get("theta", 0)
        self.tmin = kwargs.get("tmin", 0)
        self.tmax = kwargs.get("tmax", 2 * np.pi)

        self.min_val = kwargs.get("min_val", 0.1)
        self.r_trunc = self.calc_rtrunc(self.min_val)


    def prof_smoothed(self, fluxes=None, **kwargs):
        if fluxes is None:
            fluxes = self.fluxes
        return savgol_filter(fluxes, kwargs.get("window_length", 7), kwargs.get("polyorder", 1), mode='nearest')


    def oversampled_fluxes(self, dsub=0.2, smooth=True, **kwargs):
        rs_fine = np.arange(self.rs.min(), self.rs.max(), dsub)
        fluxes_fine = self.prof_interp(rs_fine)
        if smooth:
            fluxes_fine = self.prof_smoothed(fluxes_fine, window_length=kwargs.get("window_length", 11), 
                                             polyorder=kwargs.get("polyorder", 3))
        return rs_fine, fluxes_fine


    def calc_rtrunc(self, min_val = 0.1):
        rs_fine = np.linspace(self.rs.min(), self.rs.max(), 1000)
        fluxes_smoothed = self.prof_smoothed(self.prof_interp(rs_fine))

        below_min = fluxes_smoothed < min_val
        if not np.any(below_min):
            return rs_fine[-1]
        else:
            return rs_fine[np.argmax(below_min)]
        

class RadiiMeasurer:
    def __init__(self, mom0, **kwargs):
        self.mom0 = mom0
        self.name = kwargs.get("name", "Unknown")
        self.tail_angle = kwargs.get("tail_angle", 0)
        self.cent_x, self.cent_y = kwargs.get("cent", (mom0.shape[1] // 2, mom0.shape[0] // 2))
        self.ellip = kwargs.get("ellip", 0.1)  # Default value if not provided
        self.pa = kwargs.get("pa", 0)  #

        self.ra = kwargs.get("ra", 0)
        self.dec = kwargs.get("dec", 0)
        

        self.leading_side_angle = (self.tail_angle + np.pi) % (2 * np.pi)
        
        ys, xs = np.indices(mom0.shape)
        self.rmap = np.sqrt((xs - self.cent_x)**2 + (ys - self.cent_y)**2)
        self.theta_map = np.arctan2(ys - self.cent_y, xs - self.cent_x) % (2 * np.pi)
        self.deviation_map = deviation_from_major_axis(self.theta_map, self.pa)
        self.correction_map = R_ellip(1, self.ellip, self.deviation_map)

        self.theta_range = kwargs.get("theta_range", np.pi / 2)
        self.tmin = self.leading_side_angle - self.theta_range / 2
        self.tmax = self.leading_side_angle + self.theta_range / 2

        self.truncation_value = kwargs.get("truncation_value", 0.5)

        self.profiles = []
        self.r_trunc = 0
    

    def gen_theta_mask(self, tmin, tmax):
        """
        Generate a mask for the theta map based on the provided angle range.
        """
        tmin = tmin % (2 * np.pi)
        tmax = tmax % (2 * np.pi)

        if tmin < tmax:
        # Simple case: no wraparound
            return (self.theta_map >= tmin) & (self.theta_map <= tmax)
        else:
            # Wraparound case: e.g. tmin=3π/4, tmax=-3π/4
            return (self.theta_map >= tmin) | (self.theta_map <= tmax)


    def obtain_profile(self, tmin, tmax, mask=None, ellip_correction=True, **kwargs):
        if mask is None:
            mask = self.gen_theta_mask(tmin, tmax)

        mom0_masked = np.copy(self.mom0)
        mom0_masked[~mask] = np.nan  # Mask out the regions

        rs_flat = self.rmap[mask].flatten()

        if ellip_correction:
            corr_map = np.copy(self.correction_map)
            rs_flat /= corr_map[mask].flatten()


        fluxes_flat = mom0_masked[mask].flatten()

        sort_indices = np.argsort(rs_flat)
        rs_flat = rs_flat[sort_indices]
        fluxes_flat = fluxes_flat[sort_indices]

        profile = Profile(
            rs=rs_flat,
            fluxes=fluxes_flat,
            theta=np.nanmean([tmin, tmax]),
            tmin =tmin,
            tmax=tmax,
            min_val=self.truncation_value,
            rmax = kwargs.get("rmax", None), 
        )

        return profile
    

    def profile_at_angle(self, angle, width = 15, **kwargs):
        """
        Obtain a profile at a specific angle with a given width.
        """
        tmin, tmax = angle - np.deg2rad(width) / 2, angle + np.deg2rad(width) / 2
        mask = self.gen_theta_mask(tmin, tmax)

        return self.obtain_profile(tmin, tmax, mask, **kwargs)
    
    def calc_truncs(self):
        rtruncs = [prof.r_trunc for prof in self.profiles if prof.r_trunc is not None]
        if not rtruncs:
            print(f"No valid truncation radii found for {self.name}.")
            self.r_trunc = 0
            self.r_trunc_scatter = 0
            return
        self.r_trunc = np.nanmedian(rtruncs)
        self.r_trunc_scatter = np.nanstd(rtruncs)


    def measure_rtrunc(self, n_slices=12):
        
        self.theta_bins = np.linspace(self.tmin, self.tmax, n_slices)

        self.theta_avgs = np.mean([self.theta_bins[:-1], self.theta_bins[1:]], axis=0)

        for i in range(0, len(self.theta_bins) - 1):
            try:
                theta_low, theta_high = self.theta_bins[i], self.theta_bins[i + 1]

                prof = self.obtain_profile(theta_low, theta_high)
                self.profiles.append(prof)
            except Exception as e:
                continue
        
        self.calc_truncs()

    
    def get_data_row(self, arcconv=1, kpc_per_arcsec=1):
        """
        Add the measured data to the results dictionary.
        """

        return  {
            "NAME": self.name,
            "R_TRUNC": self.r_trunc,
            "R_TRUNC_ARCSEC": self.r_trunc * arcconv,
            "R_TRUNC_KPC": self.r_trunc * arcconv * kpc_per_arcsec,
            "R_TRUNC_SCATTER": self.r_trunc_scatter,
            "R_TRUNC_SCATTER_ARCSEC": self.r_trunc_scatter * arcconv,
            "R_TRUNC_SCATTER_KPC": self.r_trunc_scatter * arcconv * kpc_per_arcsec,
            "TAIL_PA": self.tail_angle,
            "MAJOR_AXIS": self.pa,
            "ELLIP": self.ellip,
            "RA": float(self.ra),
            "DEC": float(self.dec),
            "CENT_X": float(self.cent_x),
            "CENT_Y": float(self.cent_y),
        }
    

    def diagnostic_plot(self, pad=100, save=True, **kwargs):
        """
        Create a diagnostic plot for the radii measurement.
        """
        plotter = RadiiPlotter(self)
        plotter.diagnostic_plot(pad=pad, save=save, **kwargs )

    
    def plot_single_profile(self, profile, ax=None, **kwargs):
        """
        Plot a single profile on the given axis.
        """
        plotter = RadiiPlotter(self)
        plotter.plot_single_profile(profile, ax=ax, **kwargs)



class RadiiPlotter:
    def __init__(self, rm):
        self.rm = rm
    
    def plot_tail_angle(self, axis, length=25, color='blue'):
        """
        Plot the tail angle on the given axis.
        """
        axis.plot([self.rm.cent_x, self.rm.cent_x + length * np.cos(self.rm.tail_angle)],
                  [self.rm.cent_y, self.rm.cent_y + length * np.sin(self.rm.tail_angle)],
                  color=color, lw=2, label='Tail Angle')

        axis.plot([self.rm.cent_x, self.rm.cent_x + length * np.cos(self.rm.leading_side_angle)],
                  [self.rm.cent_y, self.rm.cent_y + length * np.sin(self.rm.leading_side_angle)],
                  color=color, lw=2, ls="dashed")
        

    def plot_major_axis(self, axis, length=25):
        """
        Plot the major axis on the given axis.
        """
        axis.plot([self.rm.cent_x - length * np.cos(self.rm.pa), 
                   self.rm.cent_x + length * np.cos(self.rm.pa)],
                  [self.rm.cent_y - length * np.sin(self.rm.pa), 
                   self.rm.cent_y + length * np.sin(self.rm.pa)],
                  color='red', lw=1, label='Major Axis')
        
    
    def plot_angles(self, axis, length=25):
        """
        Plot both the tail angle and the major axis on the given axis.
        """
        self.plot_tail_angle(axis, length)
        self.plot_major_axis(axis, length)

    
    def plot_slices(self, axis, length=25, cmap="rainbow"):
        """
        Plot the slices defined by the tail angle and the major axis.
        """
        for pa in [self.rm.tmin, self.rm.tmax]:
            axis.plot([self.rm.cent_x, self.rm.cent_x + length * np.cos(pa)],
                      [self.rm.cent_y, self.rm.cent_y + length * np.sin(pa)],
                      color='black', lw=3, ls="dotted", label='Slice Boundaries')
            
        for i, pa in enumerate(self.rm.theta_avgs):
            #pa += self.rm.tail_angle
            tnorm = (i + 1) / len(self.rm.theta_avgs)

            color = plt.get_cmap(cmap)(tnorm)

            axis.plot([self.rm.cent_x, self.rm.cent_x + length * np.cos(pa)],
                      [self.rm.cent_y, self.rm.cent_y + length * np.sin(pa)],
                      color=color, lw=3, ls="dashed", alpha=0.5, label='Slice Averages')
        

    def plot_ellipse(self, ax, r, ellip, pa, alpha=1):
        ell = EllipticalAperture((int(self.rm.cent_x), int(self.rm.cent_y)), r, r * (1 - ellip), theta=pa)
        ell.plot(ax=ax, color='red', lw=2, alpha=alpha,  zorder=6)


    def plot_profiles(self, ax, arcconv=1, cmap="rainbow"):
        """
        Plot the radial profiles on the given axis.
        """
        for i, prof in enumerate(self.rm.profiles):
            tnorm = (i + 1) / len(self.rm.profiles)
            
            color = plt.get_cmap(cmap)(tnorm)
            r_trunc = prof.r_trunc
            ax.plot(prof.rs, prof.fluxes, color=color, lw=2, alpha=0.7)
            ax.axvline(r_trunc * arcconv, color="black", lw=0.5, ls="dotted")

        ax.axvline(self.rm.r_trunc * arcconv, color="black", lw=5, ls="dashed", label=r'$R_{CO}$')


    def diagnostic_plot(self, pad=100, save=True, **kwargs):
        fig, ax = plt.subplots(1, 4, figsize=(12, 3.5), facecolor='white')

        ax[0].imshow(self.rm.theta_map, cmap='rainbow', origin='lower', zorder=1, alpha=0.2)
        ax[0].scatter(self.rm.cent_x, self.rm.cent_y, marker='x', color='red', s=100, zorder=2)

        ax[1].imshow(self.rm.deviation_map, cmap="Blues_r", origin="lower", alpha=0.3, zorder=0)

        for i in range(3):
            ax[i].imshow(np.log10(self.rm.mom0), cmap='Greys', origin='lower',
                          vmin=np.log10(self.rm.truncation_value), 
                          vmax=np.nanmax(np.log10(self.rm.mom0)))
            ax[i].set(xlim=(self.rm.cent_x - pad, self.rm.cent_x + pad),
                        ylim=(self.rm.cent_y - pad, self.rm.cent_y + pad))

        self.plot_angles(ax[0], length=pad / 1.5)
        ax[0].legend(loc='upper left', fontsize=10, frameon=False)

        self.plot_tail_angle(ax[1], length=pad, color="black")
        self.plot_slices(ax[1], length=pad)



        for i, r in enumerate([self.rm.r_trunc,
                  self.rm.r_trunc - self.rm.r_trunc_scatter,
                  self.rm.r_trunc + self.rm.r_trunc_scatter]):
            if r <= 0:
                continue
            alpha = 1 if i == 0 else 0.5
            self.plot_ellipse(ax[2], r, self.rm.ellip, self.rm.pa, alpha=alpha)

        ax[2].text(0.02, 0.98, f"$R_{{CO}} = {self.rm.r_trunc:.2f}$", transform=ax[2].transAxes,
                   ha='left', va='top', fontsize=10)
        ax[2].text(0.02, 0.02, f"$sigma_{{R_{{CO}}}} = {self.rm.r_trunc_scatter:.2f}$", transform=ax[2].transAxes,
                   ha='left', va='bottom', fontsize=10)

        self.plot_profiles(ax[3], arcconv=1)

        ax[3].axhline(self.rm.truncation_value, color="blue", lw=2, alpha=0.5, ls="dotted", label=r"$I_{trunc}$")

        if kwargs.get("log_y", True):
            ax[3].set_yscale('log')
        if kwargs.get("log_x", False):
            ax[3].set_xscale('log')
        
        ax[3].set_xlim(0, self.rm.r_trunc * 1.5)

        ax[3].legend(loc="lower left", fontsize=10, frameon=False, ncols=1)
        ax[0].set_title("Theta Map and\nTail/Major Axis Angles")
        ax[1].set_title("Theta Slices")
        ax[2].set_title("Overplotted Ellipse")
        ax[3].set_title("Radial Profiles")

        for axis in ax[:3]:
            axis.set(xticks=[], yticks=[])
        
        ylims = ax[3].get_ylim()
        ax[3].fill_between([self.rm.r_trunc - 1 * self.rm.r_trunc_scatter,
                            self.rm.r_trunc + 1 * self.rm.r_trunc_scatter], 
                            *ylims, color="grey", alpha=0.5)
        ax[3].set_ylim(*ylims)
        plt.tight_layout()


        if save:
            outname = kwargs.get("outname", f"plots/TruncationRadii/{self.rm.name}_diagnostic.png")
            os.makedirs(os.path.dirname(outname), exist_ok=True)

            plt.savefig(outname, dpi=200, bbox_inches='tight')
            plt.close(fig)
        else:
            plt.show()
            
    
    def plot_single_profile(self, profile, ax=None, **kwargs):
        if ax is None:
            fig, ax = plt.subplots(1, 2, figsize=(6, 4))
        
        color = kwargs.get("color", "blue")
        
        theta_mask = self.rm.gen_theta_mask(profile.tmin, profile.tmax)
        theta_mask[theta_mask == 1] = np.nan

        ax[0].imshow(self.rm.mom0, cmap='Greys', origin='lower', vmin=-3, vmax=1, alpha=0.1)
        ax[0].imshow(self.rm.mom0, cmap='Greys', origin='lower', 
                     vmin=np.log10(self.rm.truncation_value), vmax=1)

        ax[0].imshow(theta_mask, cmap='Reds', origin='lower', alpha=0.5)

        size = kwargs.get("size", 100)
        ax[0].scatter(self.rm.cent_x, self.rm.cent_y, marker='x', color='black', s=size, zorder=2)
        self.plot_tail_angle(ax[0], length=size, color=kwargs.get("color", "blue"))

        ax[0].plot([self.rm.cent_x, self.rm.cent_x + size * np.cos(profile.theta)],
                  [self.rm.cent_y, self.rm.cent_y + size * np.sin(profile.theta)],
                  color=color, lw=1, label='Major Axis')

        ax[0].set(xlim=(self.rm.cent_x - size, self.rm.cent_x + size),
                  ylim=(self.rm.cent_y - size, self.rm.cent_y + size))

        ax[1].plot(profile.rs, profile.fluxes, color=kwargs.get("color", "blue"), lw=2, alpha=0.7)
        ax[1].axvline(profile.r_trunc, color="black", lw=2, ls="dashed")
        ax[1].axhline(self.rm.truncation_value, color="blue", lw=2, alpha=0.5, ls="dotted")

        ax[1].loglog()

        if kwargs.get("grid", True):
            ax[0].grid(True, which='both', linestyle='--', linewidth=0.5)
            ax[1].grid(True, which='both', linestyle='--', linewidth=0.5)
        
        self.plot_ellipse(ax[0], profile.r_trunc, self.rm.ellip, self.rm.pa)

        if ax is None:
            plt.tight_layout()
            plt.show()


    def animated_profile_plot(self, tmin, tmax, n_thetas=20,  width = np.deg2rad(15), ax_width = 100, ylim=[1e-3, 30],
                              outname = "animated_profile.mp4"):
        
        rm = self.rm

        name = rm.name.replace(" ", "_")

        fig, ax = plt.subplots(1, 2, figsize=(12, 6), facecolor='white')
        ax[0].imshow(rm.mom0, cmap='Greys', origin='lower', vmin=-3, vmax=1, alpha=0.1)
        ax[0].imshow(rm.mom0, cmap='Greys', origin='lower',
                    vmin=np.log10(rm.truncation_value), vmax=1)   
        ax[0].set(xlim=(rm.cent_x - ax_width, rm.cent_x + ax_width),
                ylim=(rm.cent_y - ax_width, rm.cent_y + ax_width))
        tmask = rm.gen_theta_mask(tmin - width / 2, tmin + width / 2)
        # tmask[tmask > 1e-10] = np.nan
        theta_mask = ax[0].imshow(tmask, cmap='rainbow',
                                origin='lower', zorder=1, alpha=0.2)
        
        theta_bins = np.linspace(tmin, tmax, n_thetas)


        prof = rm.profile_at_angle(tmin, width=np.rad2deg(width))

        line = ax[1].plot(prof.rs, prof.fluxes, color='blue', lw=2, alpha=0.7)[0]

        hline = ax[1].axhline(y=rm.truncation_value, color='blue', lw=2, alpha=0.5, ls='dotted')
        vline = ax[1].axvline(x=prof.r_trunc, color='blue', lw=2, ls='--')

        xlim = 1, ax_width / 2
        ax[1].set(xlim=xlim, ylim=ylim, xlabel='R [pix]', ylabel='Flux [Jy beam km s$^{-1}$]')
        ax[1].loglog()
        ax[1].grid(True, which='both', linestyle='--', linewidth=0.5, zorder=0, alpha=0.5)


        plotter = RadiiPlotter(rm)
        
        thetas = np.linspace(tmin, tmax, n_thetas)
        circle_xs = rm.cent_x + prof.r_trunc * np.cos(thetas)
        circle_ys = rm.cent_y + prof.r_trunc * np.sin(thetas)
        circle = ax[0].plot(circle_xs, circle_ys, color='red', lw=2, alpha=0.5)[0]

        def update(frame):
            theta = theta_bins[frame]
            tmask = rm.gen_theta_mask(theta - width / 2, theta + width / 2)

            theta_mask.set_data(tmask)
            
            prof = rm.profile_at_angle(theta, width=np.rad2deg(width), rmax=40)
            
            rs, fluxes = prof.oversampled_fluxes(dsub=0.2, smooth=True)
            line.set_data(prof.rs, prof.fluxes)

            vline.set_xdata(prof.r_trunc)
            
            circle_xs = rm.cent_x + prof.r_trunc * np.cos(thetas)
            circle_ys = rm.cent_y + prof.r_trunc * np.sin(thetas)
            circle.set_data(circle_xs, circle_ys)

            ax[1].plot(rs, fluxes, color='Grey', lw=1, alpha=0.2)
            ax[1].set(xlim=xlim, ylim=ylim)
            ax[1].loglog()

            plt.tight_layout()
            ani = FuncAnimation(fig, func=update, frames=len(theta_bins))

            ani.save(outname, writer='ffmpeg', fps=10, dpi=200)
