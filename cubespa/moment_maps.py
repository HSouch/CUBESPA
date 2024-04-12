import os, time, warnings
from maskmoment import maskmoment

from . import utils

from numpy import round


def mask_cube(cube_fn, snr_lo=3., snr_hi=4., gain_fn=None, 
              minbeam=2., vmin=None, vmax=None,
              snr_lo_minch=3, snr_hi_minch=2,
              edgech = 5, smooth_fwhm=None,
              outdir="", prefix="", **kwargs):
    """ Mask a cube using Tony Wong's MASKMOMENT routine.
        Please note that using the current version of MASKMOMENT on PyPi will not be successful, as the velocity ranges are not specified.
        To use this iteration, you need to download at least version (1.2.1) from the Github page:
        
        https://github.com/tonywong94/maskmoment

    Args:
        cube_fn (str): Filename of the input datacube.
        snr_lo (float, optional): Low threshold bound that the mask is expanded to. Defaults to 3..
        snr_hi (float, optional): High threshold to find valid signal pixels. Defaults to 4..
        gain_fn (str, optional): Filename for the gain map (likely a .pb from CASA). Defaults to None.
        minbeam (float, optional): Number of beam areas required for valid signal. Defaults to 2.
        vmin (int, optional): Minimum velocity channel for masking. Defaults to None.
        vmax (int, optional): Maximum velocity channel for masking. Defaults to None.
        snr_lo_minch (int, optional): Number of channels the lower threshold must span. Defaults to 3.
        snr_hi_minch (int, optional): Number of channels the higher threshold must span.. Defaults to 2.
        edgech (int, optional): Number of channels at each end of velocity axis to estimate RMS from. Defaults to 5.
        smooth_fwhm (float, optional): Size of kernel to smooth mask after it is created. Defaults to None (no smoothing.)
        outdir (str, optional): Output directory to place all files. Defaults to "".
        prefix (str, optional): Prefix for output files, for example if prefix="foo", then the moment 0 map will be "foo.mom0.fits.gz". Defaults to "".
    """

    snr_lo = [snr_lo] if type(snr_lo) in (int, float) else snr_lo
    snr_hi = [snr_hi] if type(snr_hi) in (int, float) else snr_hi

    keep_uncertainties = kwargs.get("keep_uncertainties", True)
    keep_fluxcsv = kwargs.get("keep_fluxcsv", False)

    t_total_init = time.time()

    for lo in snr_lo:
        for hi in snr_hi:
            t_this = time.time()
            outname = f'{prefix}.{hi}_{lo}'
            
            print(f'\n{outname}')

            warnings.filterwarnings('ignore')
            maskmoment(img_fits=cube_fn, gain_fits=gain_fn,
                       snr_hi=hi, snr_lo=lo,
                       snr_hi_minch=snr_hi_minch, snr_lo_minch=snr_lo_minch,
                       edgech=edgech,
                       vmin=vmin, vmax=vmax,
                       fwhm=smooth_fwhm,
                       outdir=outdir, outname=outname,
                       to_kelvin=False,
                       minbeam=minbeam)
            warnings.filterwarnings('default')

            # Delete unneeded files (if requested)
            if not keep_uncertainties:
                for suffix in ['.ecube.fits.gz', '.emom0.fits.gz', '.emom1.fits.gz', '.emom2.fits.gz']:
                    try:
                        os.remove(outdir + outname + suffix)
                    except Exception as e:
                        print(f'Failed to delete uncertainty file: {suffix}', e)
            if not keep_fluxcsv:
                for suffix in ['.flux.csv']:
                    try:
                        os.remove(outdir + outname + suffix)
                    except Exception as e:
                        print(f'Failed to delete uncertainty file: {suffix}', e)
    

            print(f'\nTime Elapsed:  {round((time.time() - t_this) / 60, 2)} minutes (this masking)')
            print(f'Time Elapsed:  {round((time.time() -  t_total_init) / 60, 2)} minutes (cumulative)')


    print(f'\n Total Time Elapsed:  {round((time.time() -  t_total_init) / 60, 2)} minutes ')
