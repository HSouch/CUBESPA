from pvextractor import Path, extract_pv_slice

from . import utils

import numpy as np




def gen_pvd(cubespa_obj, center, length=50, pa=0, width=10):

    endpoints = utils.line_endpoints(center[0], center[1], length/2, pa)
    corners = utils.line_corners(center[0], center[1], width, length, pa)

    coords = np.transpose(endpoints)

    path = Path([coords[0], coords[1]], width=width)
    pvslice = extract_pv_slice(cubespa_obj.cube.data, path).data

    return {"ENDPOINTS": endpoints, "CORNERS": corners, "PVD": pvslice}
