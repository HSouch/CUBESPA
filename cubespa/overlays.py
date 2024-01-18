from astropy.io import fits
from reproject import reproject_interp

from astropy.nddata import NDData

from . import data, utils


def align_image(input_data: data.DataSet, data_to_align: data.DataSet, shape_out=None):
    """ Align two dataset objects together.

    Args:
        input_data (data.DataSet): The DataSet object that serves as the template alignment.
        data_to_align (data.DataSet): The DataSet to align to input_data

    Returns:
        DataSet: A new DataSet object where the data is aligned with input_data
    """

    input_wcs, align_wcs = utils.match_wcs_axes(input_data.wcs, data_to_align.wcs)

    nddata = NDData(data=data_to_align.data, wcs=align_wcs)

    if shape_out is None:
        shape_out = input_data.data.shape[:2]

    interp = reproject_interp(nddata, input_wcs, shape_out=shape_out)[0]
    
    return data.DataSet(data=interp, wcs=input_wcs, header=data_to_align.header, label=data_to_align.label)

