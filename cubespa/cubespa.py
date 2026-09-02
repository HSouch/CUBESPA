from . import data, utils, spectra

from . import plotting

import os, logging

import numpy as np

import maskmoment
from spectral_cube import SpectralCube

from astropy.coordinates import SkyCoord


class CubeSPA:
    """ Base input class for a CubeSPA object
    """

    def __init__(self, cube,
                 psf = None,
                 mom_maps=None, additional_maps = [],
                 center = None, position_angle = 0, eps=0,
                 vsys = 0,
                 metadata: dict = {},
                 limits = None, 
                 **kwargs) -> None:
        """ Primary CubeSPA object for data handling

        Args:
            cube (DataSet): CubeSPA dataset containing the primary cube
            psf (_type_, optional): _description_. Defaults to None.
            mom_maps (_type_, optional): _description_. Defaults to None.
            additional_maps (list, optional): _description_. Defaults to [].
            center (_type_, optional): _description_. Defaults to None.
            position_angle (int, optional): _description_. Defaults to 0.
            eps (int, optional): _description_. Defaults to 0.
            vsys (int, optional): _description_. Defaults to 0.
            metadata (dict, optional): _description_. Defaults to {}.
            limits (_type_, optional): _description_. Defaults to None.
            plot_dir (_type_, optional): _description_. Defaults to None.
        """

        self.name = kwargs.get("name", "")
        self.metadata = metadata
        data_index = kwargs.get("data_index", 0)
        self.cube = data.handle_data(cube, handler=data.load_data, data_index=data_index)

        self.vsys = vsys
        self.velocities = self.velocities_from_wcs()
        self.arcconv = self.arcconv_from_wcs()

        self.tail_angle = kwargs.get("tail_angle", None)
        
        if self.cube:
            self.cube_noise_level, self.cube_rms = utils.estimate_rms(self.cube.data,
                                                                      kwargs.get("cmin", 5),
                                                                        kwargs.get("cmax", 5))
        self.psf = data.handle_data(psf, handler=data.load_data, 
                                    data_index=data_index) if psf is not None else psf
        
        try:
            self.mom_maps = data.handle_data(mom_maps, handler=data.load_moment_maps, 
                                            data_index=data_index, 
                                            verbose=kwargs.get("verbose", True))
            self.additional_maps = additional_maps
        except Exception as e:
            logging.warn(f'Failed to load moment maps for {cube}\n{e}')

        if not self.mom_maps and kwargs.get("gen_from_cube", False):
            self.gen_mom_maps_from_cube()

        # Position angle and ellipticty (1 - b/a) of the stellar disk, for flux analysis
        self.position_angle = np.deg2rad(position_angle)
        self.eps = eps

        self.handle_center(center)
        self.handle_beam()
        self.handle_limits()
        
        self.plot_dir = kwargs.get("plot_dir", None)
        if self.plot_dir:
            os.makedirs(self.plot_dir, exist_ok=True)


    def handle_center(self, center, center_type="radec"):
        if center:
            center_x, center_y = center
            if center_type == "radec":
                cx, cy = self.cube.wcs.celestial.wcs_world2pix(center_x, center_y, 0, 0)
                self.center = (float(cx), float(cy))
            else:
                self.center = (center_x, center_y)
        else:
            try:
                sc = SkyCoord.from_name(name=self.name)
                center_x, center_y = self.cube.wcs.celestial.wcs_world2pix(sc.ra, sc.dec, 0)
                self.center = (float(center_x), float(center_y))
            except:
                print(f"Failed to get center for {self.name}.")
                self.center = (None, None)


    def handle_limits(self, limits="auto", padding=10):
        if limits == "auto":
            if self.mom_maps:
                self.limits = utils.bounds_from_moment_map(self.mom_maps.mom0.data, 
                                                        padding=padding)
            else:
                self.limits = [0, self.cube.data[0].shape[0],
                               0, self.cube.data[0].shape[1]]
        else:
            self.limits = limits


    def handle_beam(self):

        # TODO Add this to a separate function with error handling
        try:
            self.beam = [self.cube.header["BMIN"], self.cube.header["BMAJ"], self.cube.header["BPA"]]
            self.beam_pix = self.get_beam_pix()

            self.beam_area = self.get_beam_area()
            self.beam_area_arcsec = self.get_beam_area(in_pixels=False)
        except KeyError: 
            logging.warn(f"Unable to get beam for {self.cube}")
            self.beam = self.beam_pix = self.beam_area_arcsec = None


    def gen_mom_maps_from_cube(self, **kwargs):
        import warnings
        warnings.filterwarnings("ignore")

        rmscube = np.ones(self.cube.data.shape) * self.cube_rms

        maps, emaps = maskmoment.momfuncs.calc_moments(SpectralCube(self.cube.data, wcs=self.cube.wcs),
                                                    SpectralCube(rmscube, wcs=self.cube.wcs))
        mom0, mom1, mom2 = maps
        emom0, emom1, emom2 = emaps

        rms_lim = self.cube_rms * np.sqrt(len(self.velocities)) * kwargs.get("rms_sigma", 2)
        mom0_mask = mom0 < rms_lim

        for arr in [mom0, mom1, mom2, emom0, emom1, emom2]:
            arr[mom0_mask] = np.nan

        moment_maps = data.MomentMaps(mom0=mom0, mom1=mom1, mom2=mom2,
                                    emom0=emom0, emom1=emom1, emom2=emom2)
        self.mom_maps = moment_maps

        warnings.filterwarnings("always")
    

    def velocities_from_wcs(self, vsys=0):
        vmin = self.cube.wcs.wcs.crval[2]
        vdelt = self.cube.wcs.wcs.cdelt[2]

        return (np.array([vmin + i * vdelt for i in range(len(self.cube.data))]) / 1000) - vsys
    

    def arcconv_from_wcs(self):
        return self.cube.arcconv_from_wcs()
    

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

        return(plotting.moment_map_plot(self, use_limits=use_limits, **kwargs))

    
    def create_spectra(self, position, size, return_products=False, plot=False):
        aper = spectra.create_aperture(self, position, size)
        spectrum = spectra.get_spectra(self.cube.data, aper)
        if plot:
            plotting.spectra_plot(self, aper, spectrum)
        if return_products:
            return aper, spectrum
    
    def apply_limits(self, axis, pad=0):
        xmin, xmax, ymin, ymax = utils.pad_limits(self.limits, padding=pad)

        axis.set_xlim(xmin, xmax)
        axis.set_ylim(ymin, ymax)


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



class CubeSpaCollection:
    """ 
    Class to represent a collection of CubeSPA objects.
    Useful for doing bulk analyses.
    """

    def __init__(self, metadata={}):
        self.members: list[CubeSPA] = []

        self.metadata = metadata

    def names(self):
        return [self.members[i].name for i in range(len(self.members))]

    def get(self, name):
        mask = np.array(self.names()) == name

        if np.sum(mask) == 0:
            return None
        else:
            return self.members[np.argmax(mask)]
    
    @staticmethod
    def load(CO_DIR, MOM_MAP_DIR=None, ADDL_DIRS ={}, 
             ISOLIST_DIR=None, ATLAS=None, **kwargs):
        """ 
        Uses the utils.load_cubespa_objects method to generate a list of Cubespa objects dynamically

        Args:
            CO_DIR (_type_): _description_
            MOM_MAP_DIR (_type_, optional): _description_. Defaults to None.
            LEGACY_DIR (_type_, optional): _description_. Defaults to None.
            ISOLIST_DIR (_type_, optional): _description_. Defaults to None.
            ATLAS (_type_, optional): _description_. Defaults to None.

        Returns:
            CubeSpaCollection: _description_
        """
        csc = CubeSpaCollection()
        csc.members = utils.load_cubespa_objects(
            CO_DIR,
            MOM_MAP_DIR=MOM_MAP_DIR,
            ADDL_DIRS=ADDL_DIRS,
            ISOLIST_DIR=ISOLIST_DIR, 
            ATLAS=ATLAS,
            **kwargs,
        )

        csc.metadata["CO_DIR"] = CO_DIR
        csc.metadata["MOM_MAP_DIR"] = MOM_MAP_DIR

        for name, path in ADDL_DIRS.items():
            csc.metadata[name] = path
       
        csc.metadata["ISOLIST_DIR"] = ISOLIST_DIR
        csc.metadata["ATLAS"] = ATLAS

        return csc
