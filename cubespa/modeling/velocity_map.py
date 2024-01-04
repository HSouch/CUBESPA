import numpy as np
from .. import utils

def brandt_curve(rs, kwargs):
    Rmax = utils.check_kwarg("rmax", 50, kwargs)
    Vmax = utils.check_kwarg("vmax", 50, kwargs)
    n = utils.check_kwarg("n", 1, kwargs)

    R_norm = rs/Rmax
    
    return Vmax * R_norm * (1/3 + (2/3)*(R_norm ** n)) ** (-3 / (2 * n))


def df_to_of(x, y, alpha, i): 
    """ Convert the disk frame to the observer's frame.

    Args:
        x (_type_): _description_
        y (_type_): _description_
        alpha (_type_): _description_
        i (_type_): _description_

    Returns:
        _type_: _description_
    """
    X = x * np.cos(alpha) - y * np.sin(alpha) * np.cos(i)
    Y = x * np.sin(alpha) + y * np.cos(alpha) * np.cos(i)
    return X, Y

def of_to_df(X, Y, alpha, i):
    """ Convert the observer's frame to the disk frame.

    Args:
        X (_type_): _description_
        Y (_type_): _description_
        alpha (_type_): _description_
        i (_type_): _description_

    Returns:
        _type_: _description_
    """
    x = X * np.cos(alpha) + Y * np.sin(alpha)
    y = -X * np.sin(alpha) * np.cos(i) + Y * np.cos(alpha) * np.cos(i)
    
    return x, y


class VelocityModel:
    def __init__(self, shape, x0=None, y0=None, pa=None, inc=None, vsys=0, rmax=50):
        
        if type(shape) == int:
            shape = (shape, shape)

        self.shape = shape

        self.x0 = x0 if x0 is not None else int(shape[1] / 2)
        self.y0 = y0 if y0 is not None else int(shape[0] / 2)
        self.pa = np.deg2rad(pa) if pa is not None else 0
        self.inc = np.deg2rad(inc) if inc is not None else np.deg2rad(30)
        self.vsys = vsys if vsys is not None else 0
        self.rmax = rmax



    def gen_model(self, velocity_model, fill_value=np.nan, **vmod_kwargs):
        xs, ys = np.arange(0, self.shape[0]).astype(float), np.arange(0, self.shape[1]).astype(float)
        X, Y = np.meshgrid(ys, xs)

        X -= self.x0
        Y -= self.y0
    
        x, y = of_to_df(X, Y, self.pa, self.inc)
    
        r = np.sqrt(x ** 2 + y ** 2)

        vcirc = velocity_model(r, vmod_kwargs)
        theta = np.arctan2(x, y)
        
        vr = self.vsys + vcirc * np.sin(self.inc) * np.cos(theta)
        vr[r > self.rmax] = fill_value

        return vr
    
