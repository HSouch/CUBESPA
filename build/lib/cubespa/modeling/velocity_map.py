import numpy as np


def brandt_curve(R, Rmax=100, Vmax=140, n=1):
    R_norm = R/Rmax
    
    return Vmax * R_norm * (1/3 + (2/3)*(R_norm ** n)) ** (-3 / (2 * n))


def df_to_of(x, y, alpha, i): 
    X = x * np.cos(alpha) - y * np.sin(alpha) * np.cos(i)
    Y = x * np.sin(alpha) + y * np.cos(alpha) * np.cos(i)
    return X, Y

def of_to_df(X, Y, alpha, i):
    x = X * np.cos(alpha) + Y * np.sin(alpha)
    y = -X * np.sin(alpha) * np.cos(i) + Y * np.cos(alpha) * np.cos(i)
    
    return x, y


def vrad(X, Y, x0, y0, alpha, i, Vrot=brandt_curve, vsys=45, n=1, vmax=110, rmax=220, makeplots = False):

    X, Y = np.copy(X), np.copy(Y)
    
    X -= x0
    Y -= y0
    
    x, y = of_to_df(X, Y, alpha, i)
    
    r = np.sqrt(x ** 2 + y ** 2)
    
    #r[r < 0.5] = np.nan
    vcirc = Vrot(r, Rmax=rmax, Vmax=vmax, n=n)
    
    theta = np.arctan2(x, y)

    
    vr = vsys + vcirc * np.sin(i) * np.cos(theta)
    vr[r > rmax] = np.nan

        
    return vr