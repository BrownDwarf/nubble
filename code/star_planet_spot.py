#!/usr/bin/env python

# All of the argument parsing is done in the `parallel.py` module.

import multiprocessing
import time
import numpy as np
import Starfish
#from Starfish.model import ThetaParam, PhiParam
from model_planet import ThetaParamPlanet as ThetaParam
from model_planet import PhiParam

import argparse
parser = argparse.ArgumentParser(prog="star_planet_spot.py", description="Run Starfish fitting for both ambient, starspot, and planet properties")
parser.add_argument("--samples", type=int, default=5, help="How many samples to run?")
parser.add_argument("--incremental_save", type=int, default=100, help="How often to save incremental progress of MCMC samples.")
parser.add_argument("--resume", action="store_true", help="Continue from the last sample. If this is left off, the chain will start from your initial guess specified in config.yaml.")
parser.add_argument("--plot", action="store_true", help="Don't sample, just plot the config parameters")
args = parser.parse_args()

import os

import Starfish.grid_tools
from Starfish.spectrum import DataSpectrum, Mask, ChebyshevSpectrum
from Starfish.emulator import Emulator
import Starfish.constants as C
from Starfish.covariance import get_dense_C, make_k_func, make_k_func_region

from scipy.special import j1
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.linalg import cho_factor, cho_solve
from numpy.linalg import slogdet
from astropy.stats import sigma_clip

import gc
import logging

from itertools import chain
#from collections import deque
from operator import itemgetter
import yaml
import shutil
import json

from star_base import Order as OrderBase
from star_base import SampleThetaPhi as SampleThetaPhiBase

Starfish.routdir = ""

# list of keys from 0 to (norders - 1)
order_keys = np.arange(1)
DataSpectra = [DataSpectrum.open(os.path.expandvars(file), orders=Starfish.data["orders"]) for file in Starfish.data["files"]]
# list of keys from 0 to (nspectra - 1) Used for indexing purposes.
spectra_keys = np.arange(len(DataSpectra))

#Instruments are provided as one per dataset
Instruments = [eval("Starfish.grid_tools." + inst)() for inst in Starfish.data["instruments"]]


logging.basicConfig(format="%(asctime)s - %(levelname)s - %(name)s -  %(message)s", filename="{}log.log".format(
    Starfish.routdir), level=logging.DEBUG, filemode="w", datefmt='%m/%d/%Y %I:%M:%S %p')

class Order(OrderBase):

    def initialize(self, key):
        OrderBase.initialize(self, key)
        self.flux_scalar2 = None
        self.mus2, self.C_GP2 = None, None
        self.Omega2 = None
        self.emulator2 = Emulator.open(filename=Starfish.PCA["path2"])
        self.emulator2.determine_chunk_log(self.wl)
        self.pca2 = self.emulator2.pca

        # The raw eigenspectra and mean flux components
        self.EIGENSPECTRA2 = np.vstack((self.pca2.flux_mean[np.newaxis,:], self.pca2.flux_std[np.newaxis,:], self.pca2.eigenspectra))

        self.eigenspectra2 = np.empty((self.pca2.m, self.ndata))
        self.flux_mean2 = np.empty((self.ndata,))
        self.flux_std2 = np.empty((self.ndata,))
        self.flux_scalar3 = None

        self.mus3, self.C_GP3, self.data_mat3 = None, None, None


    def evaluate(self):
        '''
        Return the lnprob using the current version of the C_GP matrix, data matrix,
        and other intermediate products.
        '''

        self.lnprob_last = self.lnprob

        X = (self.chebyshevSpectrum.k * self.flux_std * np.eye(self.ndata)).dot(self.eigenspectra.T)
        X2 = (self.flux_std2 * self.eigenspectra2).T

        part1 = self.Omega**2 * self.flux_scalar**2 * X.dot(self.C_GP.dot(X.T))
        part2 = self.Omega2**2 * self.flux_scalar2**2 * X.dot(self.C_GP2.dot(X.T))
        part3 = self.data_mat

        #CC = X.dot(self.C_GP.dot(X.T)) + self.data_mat
        CC = part1 + part2 + part3

        try:
            factor, flag = cho_factor(CC)
        except np.linalg.linalg.LinAlgError:
            print("Spectrum:", self.spectrum_id, "Order:", self.order)
            self.CC_debugger(CC)
            raise

        try:
            #hard code distance to GJ1214!
            planet_solid_angle = 2.0*np.pi*(self.flux_scalar3 * (self.flux_mean2 + X2.dot(self.mus3)) /  4.5051e14)**2
            model1 = self.flux_scalar *(self.chebyshevSpectrum.k * self.flux_mean + X.dot(self.mus))
            model2 = self.flux_scalar2 * (self.chebyshevSpectrum.k * self.flux_mean + X.dot(self.mus2))

            net_model = (self.Omega - planet_solid_angle)*model1 + self.Omega2 * model2
            R = self.fl - net_model

            logdet = np.sum(2 * np.log((np.diag(factor))))
            self.lnprob = -0.5 * (np.dot(R, cho_solve((factor, flag), R)) + logdet)

            self.logger.debug("Evaluating lnprob={}".format(self.lnprob))
            return self.lnprob

        # To give us some debugging information about what went wrong.
        except np.linalg.linalg.LinAlgError:
            print("Spectrum:", self.spectrum_id, "Order:", self.order)
            raise


    def update_Theta(self, p):
        OrderBase.update_Theta(self, p)

        eigenspectra_full2 = self.EIGENSPECTRA2.copy()
        for lres2, hres2 in zip(chain([self.flux_mean2, self.flux_std2], self.eigenspectra2), eigenspectra_full2):
            interp2 = InterpolatedUnivariateSpline(self.wl_FFT_shift, hres2, k=5)
            lres2[:] = interp2(self.wl)
            del interp2

        self.emulator.params = np.append(p.teff2, p.grid[1:])
        self.mus2, self.C_GP2 = self.emulator.matrix
        self.flux_scalar2 = self.emulator.absolute_flux
        self.Omega2 = 10**p.logOmega2

        self.emulator2.params = p.grid2
        self.mus3, self.C_GP3 = self.emulator2.matrix
        self.flux_scalar3 = self.emulator2.absolute_flux
        gc.collect()

    def draw_save(self):
        '''
        Return the lnprob using the current version of the C_GP matrix, data matrix,
        and other intermediate products.
        '''

        self.lnprob_last = self.lnprob

        X = (self.chebyshevSpectrum.k * self.flux_std * np.eye(self.ndata)).dot(self.eigenspectra.T)
        X2 = (self.flux_std2 * self.eigenspectra2).T

        part1 = self.Omega**2 * self.flux_scalar**2 * X.dot(self.C_GP.dot(X.T))
        part2 = self.Omega2**2 * self.flux_scalar2**2 * X.dot(self.C_GP2.dot(X.T))
        part3 = self.data_mat

        #CC = X.dot(self.C_GP.dot(X.T)) + self.data_mat
        CC = part1 + part2 + part3

        planet_solid_angle = 2.0*np.pi*(self.flux_scalar3 * (self.flux_mean2 + X2.dot(self.mus3)) /  4.5051e14)**2
        model1 = self.flux_scalar *(self.chebyshevSpectrum.k * self.flux_mean + X.dot(self.mus))
        model2 = self.flux_scalar2 * (self.chebyshevSpectrum.k * self.flux_mean + X.dot(self.mus2))

        model_out = (self.Omega - planet_solid_angle)*model1 + self.Omega2 * model2
        if np.log10(self.Omega) > -30:
            np.save('model1.npy', model1)
            np.save('planet_solid_angle.npy', planet_solid_angle)
            np.save('Omega.npy', self.Omega)
        if np.log10(self.Omega2) > -30:
            np.save('model2.npy', model2)
            np.save('Omega2.npy', self.Omega2)
        return model_out

class SampleThetaPhi(Order, SampleThetaPhiBase):
    pass #put custom behavior here


# Run the program.

model = SampleThetaPhi(debug=True)

model.initialize((0,0))

def lnlike(p):
    # Now we can proceed with the model
    try:
        #pars1 = ThetaParam(grid=p[0:3], vz=p[3], vsini=p[4], logOmega=p[5], teff2=p[6], logOmega2=p[7])
        pars1 = ThetaParam(grid=p[0:3], grid2=p[14:], vz=p[3], vsini=p[4], logOmega=p[5], teff2=p[6], logOmega2=p[7])
        model.update_Theta(pars1)
        # hard code npoly=3 (for fixc0 = True with npoly=4)
        #pars2 = PhiParam(0, 0, True, p[6:9], p[9], p[10], p[11])
        pars2 = PhiParam(0, 0, True, p[8:11], p[11], p[12], p[13])
        model.update_Phi(pars2)
        lnp = model.evaluate()
        return lnp
    except C.ModelError:
        model.logger.debug("ModelError in stellar parameters, sending back -np.inf {}".format(p))
        return -np.inf


# Must load a user-defined prior
try:
    sourcepath_env = Starfish.config['Theta_priors']
    sourcepath = os.path.expandvars(sourcepath_env)
    with open(sourcepath, 'r') as f:
        sourcecode = f.read()
    code = compile(sourcecode, sourcepath, 'exec')
    exec(code)
    lnprior = user_defined_lnprior
    print("Using the user defined prior in {}".format(sourcepath_env))
except:
    print("Don't you want to use a user defined prior??")
    raise

# Insert the prior here
def lnprob(p):
    lp = lnprior(p)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(p)

def lnprob_all(p):
    '''Return a model draw with the input parameters'''
    pars1 = ThetaParam(grid=p[0:3], grid2=p[14:], vz=p[3], vsini=p[4], logOmega=p[5], teff2=p[6], logOmega2=p[7])
    model.update_Theta(pars1)
    # hard code npoly=3 (for fixc0 = True with npoly=4)
    #pars2 = PhiParam(0, 0, True, p[6:9], p[9], p[10], p[11])
    pars2 = PhiParam(0, 0, True, p[8:11], p[11], p[12], p[13])
    model.update_Phi(pars2)
    junk = model.evaluate()
    draw = model.draw_save()
    return draw


import emcee

start = Starfish.config["Theta"]
fname = Starfish.specfmt.format(model.spectrum_id, model.order) + "phi.json"
phi0 = PhiParam.load(fname)

ndim, nwalkers = 16, 40

p0 = np.array(start["grid"] + [start["vz"], start["vsini"], start["logOmega"], start["teff2"], start["logOmega2"]] +
             phi0.cheb.tolist() + [phi0.sigAmp, phi0.logAmp, phi0.l] + start["grid2"])

p0_std = [5, 0.02, 0.005, 0.5, 0.5, 0.01, 5, 0.01, 0.005, 0.005, 0.005, 0.01, 0.001, 0.5, 5, 0.05]

if args.plot:
    wl = model.wl
    data = model.fl
    import pandas as pd
    import json

    df_out = pd.DataFrame({'wl':wl, 'data':data})

    with open('s0_o0phi.json') as f:
        s0phi = json.load(f)

    psl = (Starfish.config['Theta']['grid']+
      [Starfish.config['Theta'][key] for key in ['vz', 'vsini', 'logOmega', 'teff2', 'logOmega2']] +
      s0phi['cheb'] +
      [s0phi['sigAmp']] + [s0phi['logAmp']] + [s0phi['l']]+
      Starfish.config['Theta']['grid2'])

    ps = np.array(psl)
    df_out['model_composite'] = lnprob_all(ps)

    pset2 = ps.copy()
    pset2[7] = -35.0
    df_out['model_hot50'] = lnprob_all(pset2)
    pset1 = ps.copy()
    pset1[5] = -35.0
    df_out['model_cool50'] = lnprob_all(pset1)

    df_out.to_csv('spec_config.csv', index=False)
    import sys
    print('Model complete, look for a file spec_config.csv to plot.  Exiting.')
    sys.exit(0)


if args.resume:
    try:
        p0_ball = np.load("emcee_chain.npy")[:,-1,:]
    except:
        final_samples = np.load("temp_emcee_chain.npy")
        max_obs = final_samples.any(axis=(0,2)).sum()
        p0_ball = final_samples[:,max_obs-1,:]
else:
    p0_ball = emcee.utils.sample_ball(p0, p0_std, size=nwalkers)

n_threads = multiprocessing.cpu_count()
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, threads=n_threads)

nsteps = args.samples
ninc = args.incremental_save
for i, (pos, lnp, state) in enumerate(sampler.sample(p0_ball, iterations=nsteps)):
    if (i+1) % ninc == 0:
        time.ctime()
        t_out = time.strftime('%Y %b %d,%l:%M %p')
        print("{0}: {1:}/{2:} = {3:.1f}%".format(t_out, i, nsteps, 100 * float(i) / nsteps))
        np.save('temp_emcee_chain.npy',sampler.chain)

np.save('emcee_chain.npy',sampler.chain)

print("The end.")