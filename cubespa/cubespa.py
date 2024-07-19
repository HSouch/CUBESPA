from . import data, utils, spectra

from . import plotting

import os, logging

import numpy as np


class CubeSPA:
    """ Base input class for a CubeSPA object
    """

    def __init__(self, cube, data_index=0, 
                 psf = None,
                 mom_maps=None, additional_maps = [],
                 center = None, position_angle = 0, eps=0,
                 vsys = 0,
                 limits = None, 
                 plot_dir = None,
                 **kwargs) -> None:
        """ CubeSPA object

        Args:
            cube (ndarray or str): Input cube or filename that points to cube.
            data_index (int, optional): If loading cube, index of cube data. Defaults to 0.
            mom_maps (str, optional): maskmoment-formatted moment map filename prefix. Defaults to None.
            additional_maps (list, optional): Additional data.DataSet objects. Defaults to [].
            center (tuple, optional): Central region, typically defined by stellar light isophotal
                center. Defaults to None.
            position_angle (float, optional): Position angle of disk. Defaults to None.
            eps (float, optional): Ellipticity (1 - b/a) of disk. Defaults to None.
            limits (array or str, optional): Bounding box containing relevant data. Defaults to None.
                If "auto", will try to automatically generate from moment maps.
            plot_dir (str, optional): Directory to place plots. Defaults to None.
        """
        
        self.cube = data.handle_data(cube, handler=data.load_data, data_index=data_index)
        
        if self.cube is not None:
            self.cube_noise_level, self.cube_rms = utils.estimate_rms(self.cube.data,
                                                                      kwargs.get("cmin", 5),
                                                                        kwargs.get("cmax", 5))

        self.psf = data.handle_data(psf, handler=data.load_data, data_index=data_index) if psf is not None else psf

        # TODO implement better handling of cx and cy. Right now assuming a 3D cube with Stokes. Can be better!
        self.center = center
        if self.center is not None:
            center_x, center_y, center_type = center
            if center_type == "radec":
                cx, cy = self.cube.wcs.wcs_world2pix(center_x, center_y, 0, 0, 0)[:2]
                self.center = (float(cx), float(cy))
            else:
                self.center = (center_x, center_y)

        # Position angle and ellipticty (1 - b/a) of the stellar disk, for flux analysis
        self.position_angle = np.deg2rad(position_angle)
        self.eps = eps

        self.velocities = self.velocities_from_wcs(vsys=vsys)
        self.vsys = vsys
        
        # TODO Add this to a separate function with error handling
        try:
            self.beam = [self.cube.header["BMIN"], self.cube.header["BMAJ"], self.cube.header["BPA"]]
            self.beam_pix = self.get_beam_pix()

            self.beam_area = self.get_beam_area()
            self.beam_area_arcsec = self.get_beam_area(in_pixels=False)
        except KeyError: 
            logging.warn(f"Unable to get beam for {cube}")
            self.beam = self.beam_pix = self.beam_area_arcsec = None

        try:
            self.mom_maps = data.handle_data(mom_maps, handler=data.load_moment_maps, data_index=data_index)
            self.additional_maps = additional_maps
        except Exception as e:
            logging.warn(f'Failed to load moment maps for {cube}\n{e}')

        if limits == "auto":
            if "padding" in kwargs.keys():
                padding = kwargs["padding"]
            else:
                padding = 10

            self.limits = utils.bounds_from_moment_map(self.mom_maps.mom0.data, padding = padding)
        
        if plot_dir is not None:
            self.plot_dir = plot_dir
            self.load_dir(plot_dir)
    

    def load_dir(path):
        if not os.path.isdir(path):
            os.makedirs(path, exist_ok=True)
        else:
            pass

    
    def velocities_from_wcs(self, vsys=0):
        vmin = self.cube.wcs.wcs.crval[2]
        vdelt = self.cube.wcs.wcs.cdelt[2]

        return (np.array([vmin + i * vdelt for i in range(len(self.cube.data))]) / 1000) - vsys
    

    def get_beam_area(self, in_pixels=True):
        try:
            header, w = self.cube.header, self.cube.wcs
            
            delt = np.mean(np.abs(w.wcs.cdelt[:2]))

            if in_pixels:
                bmin, bmaj = header["BMIN"] / delt , header["BMAJ"] / delt
            else:
                bmin, bmaj = header["BMIN"] * 3600 , header["BMAJ"] * 3600

            area = np.pi * bmaj * bmin / (4 * np.log(2))
        except Exception as e:
            print("Failed to get beam area:", e)
            area = 1
            
        return area

    def get_beam_pix(self):
        try:
            header, w = self.cube.header, self.cube.wcs
            
            delt = np.mean(np.abs(w.wcs.cdelt[:2]))
            bmaj, bmin, bpa = self.beam
            bmaj /= delt
            bmin /= delt
            return bmaj, bmin, bpa
        except:
            return 1, 1

    def get_beam_coords(self, x0=0, y0=0):
        try:
            header, w = self.cube.header, self.cube.wcs
            delt = np.mean(np.abs(w.wcs.cdelt[:2]))

            bmin, bmaj, bpa = header["BMIN"] / delt , header["BMAJ"] / delt, np.deg2rad(header["BPA"] + 90)
            return utils.ellipse_coords(x0, y0, bmaj / 2, bmin / 2, bpa)
        
        except Exception as e:
            print("Failed to get beam coordinates", e)
            
    # UTILITY FUNCTIONS
    def plot_moment_maps(self, use_limits=True, **kwargs):
        outname = kwargs.get("outname", None)

        plotting.moment_map_plot(self, use_limits=use_limits, **kwargs)

    
    def create_spectra(self, position, size, return_products=False, plot=False):
        aper = spectra.create_aperture(self, position, size)
        spectrum = spectra.get_spectra(self.cube.data, aper)
        if plot:
            plotting.spectra_plot(self, aper, spectrum)
        if return_products:
            return aper, spectrum


class CubeComparison:

    def __init__(self, cube1: CubeSPA, cube2: CubeSPA) -> None:
        self.cube1 = cube1
        self.cube2 = cube2


    def create_spectra(self, position, size, return_products=False, plot=False):
        aper = spectra.create_aperture(self, position, size)
        spectrum = spectra.get_spectra(self.cube1.data, aper)

        if plot:
            plotting.spectra_plot(self, aper, spectrum)
        if return_products:
            return aper, spectrum

    
    def readout(self):
        print(f'Cube 1: {self.cube1.cube.data.shape}    Cube 2: {self.cube2.cube.data.shape}')

