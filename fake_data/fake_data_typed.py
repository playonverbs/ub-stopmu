#!/usr/bin/env python2

import numpy as np
import pandas as pd
import pylandau
from scipy.optimize import curve_fit

import typing
from collections.abc import Callable
from collections import namedtuple

class LandauParams(namedtuple('LandauParams', ['mpv', 'eta'])):
    @property
    def scale(self, a):
        self.mpv = self.mpv * a - (2 * a * self.eta * np.log(a))/np.pi
        self.eta = self.eta * a
    def shift(self, m):
        self.mpv = self.mpv + m
    def sum(self, L_2):
        return LandauParams(mpv = self.mpv + L_2.mpv, eta = self.eta + L_2.eta)

class GLParams(namedtuple('GLParams', LandauParams._fields + ('sigma',))):
    @property
    def scale(self, a):
        self.mpv = self.mpv * a - (2 * a * self.eta * np.log(a))/np.pi
        self.eta = self.eta * a

# LandauParams = namedtuple('LandauParams', ['mpv', 'scale'])

# invert Recombination Modified Box Model to get dE/dx from dQ/dx
def ModBoxInverse(dqdx: float) -> float:
    # argon density [g/cm^3]
    rho = 1.396
    # electric field [kV/cm]
    efield = 0.273
    # ionization energy [MeV/e]
    Wion = 23.6 * (10 ** (-6))
    fModBoxA = 0.93
    fModBoxB = 0.562

    dedx = (np.exp(fModBoxB * Wion * dqdx) - fModBoxA) / fModBoxB
    return dedx

def dEdx(x, elecgain):
    dqdx = x["dqdx"] * elecgain
    return ModBoxInverse(dqdx)

def modulate_theory_mpv(
    func:        Callable[[float], float],
    # theory_mpvs: list[float],
) -> dict[str, list[float]]:
    rr_ranges = np.linspace(100, 150, 10)
    BINS = np.linspace(1, 6.0, 100)

    rr_xs         = []
    mpv_theory_xs = []
    mpv_modded_xs = []

    for n in xrange(len(rr_ranges) - 1):
        rrmin = rr_ranges[n]
        rrmax = rr_ranges[n+1]

        if (n==1) or (n==3):
            continue

        rravg = 0.5 * (rrmax + rrmin)

        rr_xs.append(rravg)
        mpv_theory_xs.append(DPDX(ERange(rravg), 0.35, 105.6))

    mpv_modded_xs = map(func, mpv_theory_xs)

    return {"rr": rr_xs, "mpv_theory": mpv_theory_xs, "mpv_modulated": mpv_modded_xs}
    
def minimise_mpvs(
    rr_vals:     list[float],
    mpv_theory:  list[float],
    mpv_data:    list[float],
    egain_range: tuple[float, float],
    egain_bins:  int
) -> tuple[float, float]:
    elec_gain = np.linspace(egain_range[0], egain_range[1], egain_bins)

    for egain in elec_gain:
        chilocal = 0.0

        # a*X ~ Landau(a*mu * log)

        for rr in rr_vals:
            dftmp = datafit.query('rr > %i and rr < %i'%(rrmin,rrmax))
        
            dedx_v = dftmp['dedx'].values
        
            vals,bine = np.histogram(dedx_v,bins=BINS)
            binc = 0.5*(bine[1:]+bine[:-1])
            guess = [1.8,0.1,0.1,30.]
    #         popt,popv = curve_fit(GL,binc,vals,p0=guess,bounds=([1.2,0.02,0.02,1],[4.0,0.3,0.3,200]))#,sigma=np.sqrt(vals),absolute_sigma=True)
            popt,popv = curve_fit(GL,binc[0:FIT_TRUNC],vals[0:FIT_TRUNC],p0=guess,bounds=([1.2,0.02,0.02,1],[4.0,0.3,0.3,200]))#,sigma=np.sqrt(vals),absolute_sigma=True)

            mpv_data = popt[0]
            mpv_err  = np.sqrt(np.diag(popv))[0]

            rr_xs.append(rravg)
            mpv_theory_xs.append(DPDX(ERange(rravg), 0.35, 105.6))

    return (1, -1)

f: Callable[[float], float] = lambda x: x**2
