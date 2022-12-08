import numpy as np
import matplotlib.pyplot as plt
import pylandau
from ELOSS import *

from collections import namedtuple

import unittest

"""
We know how the landau and gaussian distributions transform under scaling and shifting

Landau scale:
    L(m, c)
    a*L ~ L(m*a - (2 * a * c * log(a))/pi, c*a)

Gaussian scale:
    N(m, s^2),
    a*N ~ N(a*m, a^2*s^2)

But what about the Landau-Gaussian convolution?
Test this by comparing scaled GL points to parameter shifted GL:
a*GL(m, c, s)   to  GL(m*a - (2 * a * c * log(a))/pi, c*a, a^2*s^2)
"""

def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    """
    implement math.isclose() from Python 3.5, given by
    PEP 485 <https://peps.python.org/pep-0485/#proposed-implementation>
    """
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

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


sf = 2.3 # test scale factor fo

# x = np.arange(0, 100, 0.01)
# # Sole landau test
# for A, eta, mpv in (
#     (1, 1, 10),
#     (1, 2, 30),
#     (0.5, 5, 50)
# ): 
#     y_lan = pylandau.landau(x, mpv, eta, A)

def test_landau_to_landau_scaling():
    """produce ratio plot of"""

    x = np.arange(0, 100, 0.01)
    sf = 2.3

    for A, eta, mpv in (
        (1, 1, 10),
        (1, 2, 30),
        (0.5, 5, 50)
    ): 
        y1 = pylandau.landau(x, mpv, eta, A)
        y2 = pylandau.landau(x, mpv*sf - (2 * eta * np.log(sf))/np.pi, eta * sf, A)

        # plt.plot(x, y1/y2, label='A=%d, mpv=%d, eta=%d' % (A, mpv, eta))
        # plt.plot(x, y1 * sf, label='A=%d, mpv=%d, eta=%d' % (A, mpv, eta))
        # plt.plot(x, y2, label='A=%d, mpv=%d, eta=%d' % (A, mpv, eta))
        plt.plot(x, np.divide(y1, y2, out=np.zeros_like(y1), where=y2!=0), label='A=%d, mpv=%d, eta=%d' % (A, mpv, eta))

    plt.legend(loc=0, fontsize=11)
    plt.xlabel("$\Delta$ (eV)")
    plt.ylabel("$f(\Delta)$")
    plt.show()

if __name__ == '__main__':
    test_landau_to_landau_scaling()
