import sys
sys.path.append("../")

import cubespa
from numpy import array, min, max


def test_H2_mass():
    """ Test the H2 mass calculations using the values from Cramer et al, 2020 (on NGC 4402)
    """
    cloud_SCOs = [0.53, 0.53, 0.64, 0.99, 0.41, 21.4]

    cloud_masses_paper = array([6.8e5, 6.8e5, 8.1e5, 1.3e6, 5.2e5, 2.7e7])

    cloud_masses_calced = []
    d_virgo = 16.5

    for SCO in cloud_SCOs:
        cloud_masses_calced.append(cubespa.H2_Mass(SCO, D_L=d_virgo, 
                                                   z_gal=0.0008, freq=223.538,
                                                   a_CO=4.3, R_21=0.8))
    
    ratios = array(cloud_masses_calced) / cloud_masses_paper

    # Return a boolean that sees if the ratios agree to within 5%
    return (min(ratios) > 0.95) and (max(ratios) < 1.05)
