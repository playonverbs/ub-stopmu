import numpy as np
import pandas as pd
from pylandau import langau
from scipy.optimize import curve_fit

import helpers

class ExtractorConfig:
    def __init__(self):
        self.rr_fit_limits = (100, 150)
        self.sample_version = "v08_00_00_63"
        self.sample_name = "bnb_run4a_ext_v08_00_00_63_reco2_reco2_all"
        self.GL_params = {
                'guess':      [1.8, 0.1, 0.1, 30.],
                'bound_low':  [1.2, 0.02, 0.02, 1],
                'bound_high': [4.0, 0.3, 0.3 , 200]
            }

class GainExtractor:
    """
    Class to minimise MPV chi2 values and abstract the minimisation
    functionality

    Currently intended to be run under Python 2.x
    """
    def __init__(self):
        self.rr_fit_limits = (100, 150)
        self.df = None
        self.config = ExtractorConfig()

    def load_dataframe(filename="stopmu_data.csv"):
        """
        take a filename corresponding to a tagged track truncated tldqdx
        output, read it conditionally on its filetype and load it as a
        dataframe
        """
        ftype = filename.split('.')[-1]
        if ftype == "csv":
            self.df = pd.read_csv(filename)
        elif ftype == "hdf5":
            self.df = pd.read_hdf5(filename)
        else:
            self.df = None

    def minimise_chi2(steps=20):
        """
        iterate through a range of electronics gains, compute and compare the
        data to theory MPV and minimise the comparison chi2 between the two.
        """
        elec_gain = np.linspace(230,250,steps)
        rr_ranges = np.linspace(self.rr_fit_limits[0], self.rr_fit_limits[1], 10)
        chisq_v = []
        BINS = np.linspace(1,6.0,100)
        datafit = self.df.query( 'rr > 100 and rr < 150 and pitch > 0.3 and\
                pitch < 0.4 and px < 0.2 and px >-0.2')
        print(datafit.shape)
        
        for elecgain in elec_gain:
            datafit['dedx'] = datafit.apply(lambda x: helpers.dEdx(x, elecgain), axis=1)

            chilocal = 0.

            for n in xrange(len(rr_ranges)-1):
                rrmin = rr_ranges[n]
                rrmax = rr_ranges[n+1]

                if (n == 1) or (n == 3):
                    #print "skipping"
                    continue
                rravg = 0.5*(rrmax + rrmin)

                dftmp = datafit.query('rr > %i and rr < %i'%(rrmin, rrmax))
                dedx_v = dftmp['dedx'].values

                vals, bine = np.histogram(dedx_v, bins=BINS)
                binc = 0.5*(bine[1:]+bine[:-1])
                popt, popv = curve_fit(helpers.GL, binc, vals,
                        p0=self.config.GL_params['guess'],
                        bounds=(self.config.GL_params['bound_low'],self.config.GL_params['bound_high']))

                mpv_data, mpv_err = popt[0], np.sqrt(np.diag(popv))[0]
                mpv_theory = DPDX(helpers.ERange(rravg), 0.35, 105.6)
                
                chilocal += ( (mpv_data - mpv_theory)**2 ) / (mpv_err**2 + 0.015**2 + 0.01**2)

            print('Gain = %.02f -> chisq = %.02f' % (elecgain, chilocal))
            chisq_v.append(chilocal / float(len(rr_ranges)-1.))

        chi_min = 100000.
        gain_min = 0.

        for i, chi in enumerate(chisq_v):
            if chi < chi_min:
                chi_min = chi
                gain_min = elec_gain[i]

