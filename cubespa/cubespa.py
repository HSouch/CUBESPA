from . import data, utils, spectra

from . import plotting

import os

import numpy as np


class CubeSPA:
    """ Base input class for a CubeSPA object
    """

    def __init__(self, cube, data_index=0, 
                 psf = None,
                 mom_maps=None, additional_maps = [],
                 center = None, position_angle = None, eps=None,
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
                                                                    utils.check_kwarg("cmin", 5, kwargs), 
                                                                    utils.check_kwarg("cmax", 5, kwargs))

        self.psf = data.handle_data(psf, handler=data.load_data, data_index=data_index) if psf is not None else psf


        self.center = center
        self.position_angle = position_angle
        self.eps = eps

        self.velocities = self.velocities_from_wcs(vsys=vsys)
        self.vsys = vsys
        self.beam_area = self.get_beam_area()

        self.mom_maps = data.handle_data(mom_maps, handler=data.load_moment_maps, data_index=data_index)
        self.additional_maps = additional_maps

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
    

    def get_beam_area(self):
        try:
            header, w = self.cube.header, self.cube.wcs
            
            delt = np.mean(np.abs(w.wcs.cdelt[:2]))

            bmin, bmaj = header["BMIN"] / delt , header["BMAJ"] / delt
            
            area = np.pi * bmaj * bmin / (4 * np.log(2))
        except Exception as e:
            print("Failed to get beam area:", e)
            area = 1
            
        return area
    
    # UTILITY FUNCTIONS
    def plot_moment_maps(self, use_limits=True, **kwargs):
        filename = utils.check_kwarg("filename", None, kwargs)

        plotting.moment_map_plot(self, use_limits=use_limits, filename=filename, kwargs=kwargs)

    
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

