import numpy as np
from scipy.interpolate import spline
from pylandau import langau

class PDGLookup():
    def __init__(self):
        self.load("./data/mutable.txt")

    def from_file(filename):
        self.rr_v   = None
        self.dedx_v = None
        seld.ke_v   = None

    def fit_spline():
       self.pdg_rr_fit_v   = np.linspace(1, 1000, 1000)
       self.pdg_dedx_fit_v = spline(pdg_rr_v, pdg_dedx_v, pdg_rr_fit_v)
       self.pdg_ke_fit_v   = spline(pdg_rr_v, pdg_ke_v, pdg_rr_fit_v)

pdg_rr_fit_v   = np.linspace(1, 1000, 1000)
pdg_dedx_fit_v = spline(pdg_rr_v, pdg_dedx_v, pdg_rr_fit_v)
pdg_ke_fit_v   = spline(pdg_rr_v, pdg_ke_v, pdg_rr_fit_v)

def ERange(RR):
    for i in xrange(len(pdg_rr_fit_v)):
        if (RR <= pdg_rr_fit_v[i]):
            return pdg_ke_fit_v[i]
    return -1

def ModBoxInverse(dqdx, fModBoxA = 0.93, fModBoxB = 0.562):
    Wion = 23.6 * (10 ** (-6))
    dedx = (np.exp(fModBoxB * Wion * dqdx) - fModBoxA) / fModBoxB
    return dedx

def dEdx(x, elecgain):
    dqdx = x["dqdx"] * elecgain
    return ModBoxInverse(dqdx)

def GL(x_v, mpv, sL, sG, A):
    return A * langau(x_v * 100.0, mpv * 100, sL * 100, sG * 100)
