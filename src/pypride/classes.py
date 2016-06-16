# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 17:05:04 2013

Definitions of classes used in pypride

@author: Dmitry A. Duev
"""

import ConfigParser
import datetime
from astropy.time import Time

#import sys
import numpy as np
import scipy as sp
#from matplotlib.cbook import flatten
from math import *
#from math import sin, cos, sqrt, floor
from sklearn.base import BaseEstimator
from sklearn.grid_search import GridSearchCV
from sklearn.linear_model import LinearRegression

import struct
import os
import inspect
from copy import deepcopy

#from pypride.vintflib import lagint
try:
    from pypride.vintflib import lagint, pleph#, iau_xys00a_fort, admint2
except:
    # compile the Fortran code if necessary
    from numpy import f2py
    abs_path = os.path.dirname(inspect.getfile(inspect.currentframe()))
    fid = open(os.path.join(abs_path, 'vintflib.f'))
    source = fid.read()
    fid.close()
    f2py.compile(source, modulename='vintflib')
    from pypride.vintflib import lagint, pleph#, iau_xys00a_fort, admint2

#from time import time as _time
#from numba import jit

cheb = np.polynomial.chebyshev
norm = np.linalg.norm

'''
#==============================================================================
# Polynomial regression
#==============================================================================
'''
class PolynomialRegression(BaseEstimator):
    def __init__(self, deg=None):
        self.deg = deg
        self.model = LinearRegression(fit_intercept=False)
    
    def fit(self, X, y):
        # self.model.fit(np.vander(X, N=self.deg + 1), y, n_jobs=-1)
        self.model.fit(np.vander(X, N=self.deg + 1), y)
    
    def predict(self, x):
        try:
            len(x)
            x = np.array(x)
        except:
            x = np.array([x])
        
        return self.model.predict(np.vander(x, N=self.deg + 1))
    
    @property
    def coef_(self):
        return self.model.coef_

'''
#==============================================================================
# Chebyshev regression
#==============================================================================
'''
class ChebyshevRegression(BaseEstimator):
    def __init__(self, deg=None):
        self.deg = deg
    
    def fit(self, X, y):
        self.chefit = cheb.chebfit(X, y, self.deg)
    
    def predict(self, x):
        return cheb.chebval(x, self.chefit)
    
    @property
    def coef_(self):
        return self.chefit


'''
#==============================================================================
# Optimal fit to data
#==============================================================================
'''
def optimalFit(x, y, min_order=0, max_order=8, fit_type='poly'):
    # initialise optimal estimator:
    if fit_type == 'poly':
        estimator = PolynomialRegression()
    elif fit_type == 'cheb':
        estimator = ChebyshevRegression()
    else:
        raise Exception('unknown fit type')
    degrees = np.arange(min_order, max_order)
    cv_model = GridSearchCV(estimator,
                            param_grid={'deg': degrees},
                            scoring='mean_squared_error')
    cv_model.fit(x, y)
    # use as: cv_model.predict(x_new)
    return cv_model


'''
#==============================================================================
# Input settings
#==============================================================================
'''
class inp_set(object):
    """
    input settings: catalogs, directories, etc
    """
    def __init__(self, inp_file):
        self.inp_file = inp_file

        self.config = ConfigParser.RawConfigParser()
        self.config.read(self.inp_file)

        # absolute path
#        self.abs_path = self.config.get('Catalogues', 'abs_path')
#        self.abs_path = os.getcwd()
#        self.abs_path = os.path.dirname(os.path.abspath(__file__))
        self.abs_path = os.path.dirname(inspect.getfile(inspect.currentframe()))
        
        # paths
        self.sta_xyz = self.config.get('Catalogues', 'sta_xyz')
        self.sta_vxvyvz = self.config.get('Catalogues', 'sta_vxvyvz')
        self.sta_axo = self.config.get('Catalogues', 'sta_axo')
        self.oc_load = self.config.get('Catalogues', 'oc_load')
        self.atm_load = self.config.get('Catalogues', 'atm_load')
        self.sta_nam = self.config.get('Catalogues', 'sta_nam')
        self.cat_eop = self.config.get('Catalogues', 'cat_eop')
        self.sta_thermdef = self.config.get('Catalogues', 'sta_thermdef')
        self.source_cat = self.config.get('Catalogues', 'source_cat')
        self.source_nam = self.config.get('Catalogues', 'source_nam')
        self.shnames_cat = self.config.get('Catalogues', 'shnames_cat')
        self.shnames_cat_igs = self.config.get('Catalogues', 'shnames_cat_igs')
        self.meteo_cat = self.config.get('Catalogues', 'meteo_cat')
        self.ion_cat = self.config.get('Catalogues', 'ion_cat')
        self.f_ramp = self.config.get('Catalogues', 'f_ramp')
        self.f_ramp1w =  self.config.get('Catalogues', 'f_ramp1w')
        self.f_gc = self.config.get('Catalogues', 'f_gc')
        self.jpl_eph = self.config.get('Ephemerides', 'jpl_eph')
        self.sc_eph_cat = self.config.get('Ephemerides', 'sc_eph_cat')
        self.obs_path = self.config.get('Directories', 'obs_path')
        self.out_path = self.config.get('Directories', 'out_path')
        
        # relative paths? make them absolute:
        if self.sta_xyz[0]!='/':
            self.sta_xyz = os.path.join(self.abs_path, self.sta_xyz)
        if self.sta_vxvyvz[0]!='/':
            self.sta_vxvyvz = os.path.join(self.abs_path, self.sta_vxvyvz)
        if self.sta_axo[0]!='/':
            self.sta_axo = os.path.join(self.abs_path, self.sta_axo)
        if self.oc_load[0]!='/':
            self.oc_load = os.path.join(self.abs_path, self.oc_load)
        if self.atm_load[0]!='/':
            self.atm_load = os.path.join(self.abs_path, self.atm_load)
        if self.sta_nam[0]!='/':
            self.sta_nam = os.path.join(self.abs_path, self.sta_nam)
        if self.cat_eop[0]!='/':
            self.cat_eop = os.path.join(self.abs_path, self.cat_eop)
        if self.sta_thermdef[0]!='/':
            self.sta_thermdef = os.path.join(self.abs_path, self.sta_thermdef)
        if self.source_cat[0]!='/':
            self.source_cat = os.path.join(self.abs_path, self.source_cat)
        if self.source_nam[0]!='/':
            self.source_nam = os.path.join(self.abs_path, self.source_nam)
        if self.shnames_cat[0]!='/':
            self.shnames_cat = os.path.join(self.abs_path, self.shnames_cat)
        if self.shnames_cat_igs[0]!='/':
            self.shnames_cat_igs = os.path.join(self.abs_path, self.shnames_cat_igs)
        if self.meteo_cat[0]!='/':
            self.meteo_cat = os.path.join(self.abs_path, self.meteo_cat)
        if self.ion_cat[0]!='/':
            self.ion_cat = os.path.join(self.abs_path, self.ion_cat)
        if self.f_ramp[0]!='/':
            self.f_ramp = os.path.join(self.abs_path, self.f_ramp)
        if self.f_ramp1w[0]!='/':
            self.f_ramp1w = os.path.join(self.abs_path, self.f_ramp1w)
        if self.f_gc[0]!='/':
            self.f_gc = os.path.join(self.abs_path, self.f_gc)
        if self.jpl_eph[0]!='/':
            self.jpl_eph = os.path.join(self.abs_path, self.jpl_eph)
        if self.sc_eph_cat[0]!='/':
            self.sc_eph_cat = os.path.join(self.abs_path, self.sc_eph_cat)
        if self.obs_path[0]!='/':
            self.obs_path = os.path.join(self.abs_path, self.obs_path)
        if self.out_path[0]!='/':
            self.out_path = os.path.join(self.abs_path, self.out_path)
        
        # phase center:
        self.phase_center = self.config.get('Models', 'phase_center')
        # Near-field model:
        self.nf_model = self.config.get('Models', 'nf_model')
        #tropo and iono:
        self.do_trp_calc = self.config.getboolean('Switches', 'do_trp_calc')
        self.tropo_model = self.config.get('Models', 'tropo_model')
        self.do_trp_grad_calc = self.config.getboolean('Switches', 'do_trp_grad_calc')
        self.do_ion_calc = self.config.getboolean('Switches', 'do_ion_calc')
        self.iono_model = self.config.get('Models', 'iono_model')
        # Calculate delay prediction?
        self.delay_calc = self.config.getboolean('Switches', 'delay_calc')
        # Calculate uvws?
        self.uvw_calc = self.config.getboolean('Switches', 'uvw_calc')
        # note that 2nd station is the transmitter if 2(3)-way:
        self.doppler_calc = self.config.getboolean('Switches', 'doppler_calc')
        # Doppler prediction model
        self.dop_model = self.config.get('Models', 'dop_model')
        # Doppler mode parametres for 2(3)-way
#        self.uplink_sta = self.config.get('Models', 'uplink_sta')
#        self.freq_type = self.config.get('Models', 'freq_type')
#        self.freq = self.config.getfloat('Models', 'freq')
        self.tr = self.config.getfloat('Models', 'tr')
        # Jacobians
        # generate additional ephs and calc delays for S/C position correction
        self.sc_rhophitheta = self.config.getboolean('Switches', 
                                                     'sc_rhophitheta')
        self.mas_step = self.config.getfloat('Models', 'mas_step')
        self.m_step = self.config.getfloat('Models', 'm_step')
        # generate additional ephs and calc delays for RA/GNSS XYZ pos corr
        self.sc_xyz = self.config.getboolean('Switches', 'sc_xyz')
        self.m_step_xyz = self.config.getfloat('Models', 'm_step_xyz')
        # force update s/c ephemetis
        self.sc_eph_force_update = self.config.getboolean('Switches', \
                                                         'sc_eph_force_update')
        
    def get_section(self, section='all'):
        # returns dictionary containing section (without '__name__' key)
        if section!='all':
            if section in ('Catalogues', 'Ephemerides', 'Directories'):
                out = dict((k,os.path.join(self.abs_path, v)) for k, v in \
                          self.config._sections[section].iteritems() if k!='__name__')
            else:
                out = dict((k,v) for k, v in \
                      self.config._sections[section].iteritems() if k!='__name__')
        else:
            out = {}
            for section in self.config._sections.keys():
                if section in ('Catalogues', 'Ephemerides', 'Directories'):
                    o = dict((k,os.path.join(self.abs_path, v)) for k, v in \
                              self.config._sections[section].iteritems() if k!='__name__')
                else:
                    o = dict((k,v) for k, v in \
                          self.config._sections[section].iteritems() if k!='__name__')
                out.update(o)
        # replace stupidness:
        for k, v in out.iteritems():
            if v=='False': out[k] = False
            if v=='True': out[k] = True
            if v=='None': out[k] = None
            if k in ('tr', 'mas_step', 'm_step', 'm_step_xyz'):
                out[k] = float(out[k])
        return out

'''
#==============================================================================
# Obs object - contains info about obs to be processed
#==============================================================================
'''
class obs(object):
    """
    ([sta], source, sou_type, start, step, stop)
    the class contains data previously stored in .obs-files
    for each source on each baseline back in Matlab days
    """
    def __init__(self, sta, source, sou_type, exp_name='', sou_radec=None, inp=None):
        #self.sta = [sta1, sta2]
        self.sta = sta
        self.source = source
        if sou_type not in ('C', 'S', 'R', 'G'):
            raise Exception('Unrecognised source type.')
        self.sou_type = sou_type
        self.sou_radec = sou_radec
        self.exp_name = exp_name
        self.upSta = [] # list of uplink stations in case of n-way Doppler, n>1
        self.tstamps = []   # list of datetime objs [datetime(2013,9,18,14,28),.]
        self.scanStartTimes = []  # might be useful.. integrate the phase?
        self.scanStopTimes = []  # to put trailing zeros in proper places (.del)
        self.freqs = []  # list of obs freqs as defined in $MODE
        self.dude = dude()  # Delays, Uvws, Doppler, Etc. = DUDE
        # pointings
        self.pointingsJ2000 = []
        self.pointingsDate = []
        self.azels = []
        # input switches (like what to calculate and what not)
        # default values
        self.inp = {'do_trp_calc':False, 'do_ion_calc':False,
                    'delay_calc':False, 'uvw_calc':False,
                    'doppler_calc':False, 'sc_rhophitheta':False,
                    'sc_xyz':False}
        # set things that have been passed in inp
        if inp!=None:
            for k,v in inp.iteritems():
                self.inp[k] = v
            self.inp = inp

    def __str__(self):
        if self.inp['delay_calc']:
            echo = 'baseline: {:s}-{:s}, source: {:s}'.\
                    format(self.sta[0], self.sta[1], self.source)
            if len(self.scanStartTimes)>0:
                echo += ', time range: {:s} - {:s}, n_obs: {:d}'.\
                        format(str(self.scanStartTimes[0]),
                               str(self.scanStopTimes[-1]), len(self.tstamps))
            return echo

        elif self.inp['doppler_calc']:
            if len(self.sta)==2:
                echo = 'RX station: {:s}, TX station: {:s}, source: {:s}'.\
                        format(self.sta[0], self.sta[1], self.source)
            elif len(self.sta)==1:
                echo = 'RX station: {:s}, source: {:s}'.\
                        format(self.sta[0], self.source)
            if len(self.scanStartTimes)>0:
                echo += ', time range: {:s} - {:s}, n_obs: {:d}'.\
                        format(str(self.scanStartTimes[0]),
                               str(self.scanStopTimes[-1]), len(self.tstamps))
            return echo
        
        else:
            if len(self.scanStartTimes)>0:
                echo = 'time range: {:s} - {:s}, n_obs: {:d}'.\
                        format(str(self.scanStartTimes[0]),
                               str(self.scanStopTimes[-1]), len(self.tstamps))
            return 'probably a fake obs object for one reason or another\n' + echo

#    def __str__(self):
#        if self.tstamps == []:
#            if self.exp_name == '':
#                return "baseline: %s-%s\nsource: %s\nsouce type: %s" % (self.sta[0], 
#                    self.sta[1], self.source, self.sou_type)
#            else:
#                return "exp name: %s\nbaseline: %s-%s\nsource: %s\nsouce type: %s" \
#                 % (self.exp_name, self.sta[0], self.sta[1], self.source, self.sou_type)
#        else:
#            tf = []
#            for jj in range(len(self.tstamps)):
#                if len(self.freqs)!=0:
#                    tf.append(str(self.tstamps[jj]) + ' f=' + str(self.freqs[jj]) + ' Hz')
#                else:
#                    tf.append(str(self.tstamps[jj]))
#            if self.exp_name == '':
#                return "baseline: %s-%s\nsource: %s\nsouce type: %s\nscans:\n%s" \
#                        % (self.sta[0], self.sta[1], self.source, self.sou_type,
#                           '\n'.join(map(str, tf )))
#            else:
#                return "exp name: %s\nbaseline: %s-%s\nsource: %s\nsouce type: %s\nscans:\n%s" \
#                      % (self.exp_name, self.sta[0], self.sta[1], self.source, self.sou_type,
#                           '\n'.join(map(str, tf )))

    @staticmethod
    def factors(n):
        # factorise a number
        facs = set(reduce(list.__add__, 
                    ([i, n//i] for i in range(1, int(n**0.5) + 1) if n % i == 0)))
        return sorted(list(facs))

    def resample(self, tstep=1):
        # resample scan time stamps
        tstamps = []
#        freqs = []
        for start, stop in zip(self.scanStartTimes, self.scanStopTimes):
            stopMinusStart = stop - start
            # step must be a multiple of stopMinusStart.total_seconds()
            # also allow for non-integer step
#            if (stopMinusStart.total_seconds())%step == 0:
            if (stopMinusStart.total_seconds())%tstep < 1e-9:
                nobs = stopMinusStart.total_seconds()/tstep + 1
            else:
                print 'bad t_step: not a multiple of N_sec. set to nearest possible.'
                facs = self.factors(stopMinusStart.total_seconds())
                # use smaller value (to the left). max for if pos is [0]
                pos = max(0, np.searchsorted(facs, tstep) - 1)
                tstep = facs[pos]
                nobs = stopMinusStart.total_seconds()/tstep + 1
            for ii in range(int(nobs)):
                tstamps.append(start + ii*datetime.timedelta(seconds=tstep))
#            if len(self.freqs)>0:
#                f = [ f for t, f in zip(self.tstamps, self.freqs) \
#                    if start <= t <= stop ][0]
#                for ii in range(int(nobs)):
#                    freqs.append(f)
        self.tstamps = tstamps
#        self.freqs = freqs

#    @staticmethod
#    def freqRamp(cat_dir='cats/ramp.', sc=None):
#        # read frequency ramping parameters
#        if sc==None:
#            raise Exception('Can\'t load ramp params: S/C not specified.')
#        else:
#            rampFile = ''.join((cat_dir, sc.lower()))
#            try:
#                with open(rampFile, 'r') as f:
#                    f_lines = f.readlines()
#                    ramp = []
#                    for line in f_lines:
#                        line = line.split()
#                        t_start = ''.join((line[0], ' ', line[1]))
#                        t_stop = ''.join((line[2], ' ', line[3]))
#                        # [t_start t_stop f_0 df uplink_sta]
#                        ramp.append([datetime.datetime.strptime(t_start, \
#                                                      "%Y-%m-%d %H:%M:%S"), \
#                                    datetime.datetime.strptime(t_stop, \
#                                                      "%Y-%m-%d %H:%M:%S"), \
#                                    float(line[4]), float(line[5]), line[6] ])
#                return ramp
#            except Exception, err:
#                print str(err)
#                print 'Could not load ramp params for ' + sc + '.'
    
    def addScan(self, start, step, stop=None, nobs=None, freq=None,
                force_step=False):
        # start and stop must be a datetime object
        # step is an integer time step in seconds

        if stop is None and nobs is None:
            raise Exception('You should specify either "stop" time \
                            or number of time stamps "nobs" to add')
        elif stop is None:
            # nobs is set
            self.scanStartTimes.append(start)
            nobs = int(nobs)  # this is needed for the range() function
            stop = start + (nobs-1)*datetime.timedelta(seconds=step)
            self.scanStopTimes.append(stop)                            
            for ii in range(nobs):
                tmp = start + ii*datetime.timedelta(seconds=step)
                # avoid duplicates! this might become too slow, so drop it...
#                if tmp not in self.tstamps:
                self.tstamps.append(tmp)
        elif nobs is None:
            # stop is set
            self.scanStartTimes.append(start)
            self.scanStopTimes.append(stop)
            stopMinusStart = stop - start
            # step must be a multiple of stopMinusStart.total_seconds()
            # also allow for non-integer step
#            if (stopMinusStart.total_seconds())%step == 0:
#            if (stopMinusStart.total_seconds())%step < 1e-9:
            if modf(stopMinusStart.total_seconds()/step)[0] < 1e-9:
                nobs = int(stopMinusStart.total_seconds()/step + 1)
            else:
                if not force_step:
                    print 'bad t_step: not a multiple of N_sec. set to nearest possible.'
                    facs = self.factors(stopMinusStart.total_seconds())
                    # use smaller value (to the left). max for if pos is [0]
                    pos = max(0, np.searchsorted(facs, step) - 1)
                    step = facs[pos]
                    nobs = int(stopMinusStart.total_seconds() / step + 1)
                else:
                    # force time step?
                    nobs = int(stopMinusStart.total_seconds() // step + 1)

            for ii in range(int(nobs)):
                tmp = start + ii*datetime.timedelta(seconds=step)
                # avoid duplicates! this might become too slow, so drop it...
#                if tmp not in self.tstamps:
                self.tstamps.append(tmp)

            # fix the case with the forced step:
            if force_step:
                self.tstamps.append(stop)
        
        # frequency for ionospheric delay calculation
        if freq is not None:
            self.freqs.append([start, stop, freq, 0, None])

    def splitScans(self, duration=60):
        """
            Split longer scans into shorter ones.
            useful for running together with handy.py (on very long scan) and self.smoothDude:
            increase t_step for a faster computation, then split into shorter subscans,
            and run self.smoothDude()
        Args:
            duration: subscan duration in seconds. if last subscan is shorter,
                      it get appended to the last but one

        Returns:

        """
        scanStartTimes = []
        scanStopTimes = []
        freqs = []
        s = 0
        for start, stop in zip(self.scanStartTimes, self.scanStopTimes):
            stopMinusStart = stop - start
            # skip scan if it's already shorter than duration
            if stopMinusStart < datetime.timedelta(seconds=duration):
                continue
            n_subscans = int(stopMinusStart.total_seconds() // duration)
            for ii in range(n_subscans):
                subScanStartNominal = start + ii * datetime.timedelta(seconds=duration)
                subScanStopNominal = start + (ii+1) * datetime.timedelta(seconds=duration)
                tstamps_cut = [t for t in self.tstamps
                               if subScanStartNominal <= t <= subScanStopNominal]
                scanStartTimes.append(tstamps_cut[0])
                scanStopTimes.append(tstamps_cut[-1])
                if len(self.freqs) > 0:
                    freqs.append([start, stop, self.freqs[s][2], 0, None])
            # fix last subscan end time
            scanStopTimes[-1] = stop
            s += 1
        # update self:
        self.scanStartTimes = scanStartTimes
        self.scanStopTimes = scanStopTimes
        self.freqs = freqs

    def smoothPointing(self, tstep=1, method='cheb'):
        """
        Scan-based smoothing of pointing data
        """
        if method == 'poly':
            # initialise optimal polinomial estimator:
            estimator = PolynomialRegression()
            _degrees = np.arange(0, 8)
            cv_model = GridSearchCV(estimator,
                                    param_grid={'deg': _degrees},
                                    scoring='mean_squared_error')
        elif method == 'cheb':
            # initialise optimal polinomial estimator:
            estimator = ChebyshevRegression()
            _degrees = np.arange(0, 10)
            cv_model = GridSearchCV(estimator,
                                    param_grid={'deg': _degrees},
                                    scoring='mean_squared_error')
        else:
            raise Exception('Unknown smoothing method. Use \'poly\' or \'cheb\'')

        # pointingsJ2000 = np.empty(shape=(len(self.sta), 0, 2))
        # pointingsDate = np.empty(shape=(len(self.sta), 0, 2))
        # azels = np.empty(shape=(len(self.sta), 0, 2))
        pointingsJ2000 = [[]*len(self.sta)]
        pointingsDate = [[]*len(self.sta)]
        azels = [[]*len(self.sta)]

        # iterate over scans, as fits must be scan-based
        for start, stop in zip(self.scanStartTimes, self.scanStopTimes):
            dd0 = datetime.datetime(start.year, start.month, start.day)
            # time scale and proper indices to make fit
            time = np.array([
                    ( ii, t.hour*3600 + t.minute*60.0 + t.second +
                      (t-dd0).days*86400.0 )
                    for ii, t in enumerate(self.tstamps)
                    if start <= t <= stop ])

            # time scale used for smoothing:
            t_dense = np.arange(time[0, 1], time[-1, 1]+tstep, tstep)

            # renorm
            time[:, 1] = 24.0*time[:, 1]/86400.0
            t_dense = 24.0*t_dense/86400.0

            # iterate over stations:
            for jj, _ in enumerate(self.sta):
                if len(self.pointingsJ2000) > 0:
                    pointingsJ2000_scan = []
                    for ii in range(self.pointingsJ2000.shape[2]):
                        ind = map(int, time[:, 0])
                        cv_model.fit(time[:, 1], np.unwrap(self.pointingsJ2000[ind, jj, ii]))
                        if ii == 1:
                            # fix RA if necessary
                            # wrap back if necessary:
                            wrap = np.angle(np.exp(1j*cv_model.predict(t_dense)))
                            negative = wrap < 0
                            wrap[negative] += 2*np.pi
                            pointingsJ2000_scan.append(wrap)
                        else:
                            pointingsJ2000_scan.append(cv_model.predict(t_dense))
                    try:
                        pointingsJ2000[jj] = np.vstack((pointingsJ2000[jj],
                                                              np.array(pointingsJ2000_scan).T))
                    except:
                        pointingsJ2000[jj] = np.array(pointingsJ2000_scan).T

                if len(self.pointingsDate) > 0:
                    pointingsDate_scan = []
                    for ii in range(self.pointingsDate.shape[2]):
                        ind = map(int, time[:, 0])
                        cv_model.fit(time[:, 1], np.unwrap(self.pointingsDate[ind, jj, ii]))
                        if ii == 1:
                            # fix RA if necessary
                            wrap = np.angle(np.exp(1j * cv_model.predict(t_dense)))
                            negative = wrap < 0
                            wrap[negative] += 2 * np.pi
                            pointingsDate_scan.append(wrap)
                        else:
                            pointingsDate_scan.append(cv_model.predict(t_dense))
                    try:
                        pointingsDate[jj] = np.vstack((pointingsDate[jj],
                                                             np.array(pointingsDate_scan).T))
                    except:
                        pointingsDate[jj] = np.array(pointingsDate_scan).T

                if len(self.azels) > 0:
                    azels_scan = []
                    for ii in range(self.azels.shape[2]):
                        ind = map(int, time[:, 0])
                        cv_model.fit(time[:, 1], np.unwrap(self.azels[ind, jj, ii]*np.pi/180.0))
                        if ii == 1:
                            wrap = np.angle(np.exp(1j * cv_model.predict(t_dense)))
                            negative = wrap < 0
                            wrap[negative] += 2 * np.pi
                            azels_scan.append(wrap)
                        else:
                            azels_scan.append(cv_model.predict(t_dense))
                    try:
                        azels[jj] = np.vstack((azels[jj], np.array(azels_scan).T*180.0/pi))
                    except:
                        azels[jj] = np.array(azels_scan).T*180.0/pi

        # update
        self.pointingsJ2000 = np.swapaxes(np.array(pointingsJ2000), 0, 1)
        self.pointingsDate = np.swapaxes(np.array(pointingsDate), 0, 1)
        self.azels = np.swapaxes(np.array(azels), 0, 1)

        # resample tstamps and freqs:
        self.resample(tstep=tstep)

        def smoothDude(self, tstep=1, method='cheb'):
            """
            Scan-based smoothing of DUDE data
            """
            if method == 'poly':
                # initialise optimal polinomial estimator:
                estimator = PolynomialRegression()
                _degrees = np.arange(0, 8)
                cv_model = GridSearchCV(estimator,
                                        param_grid={'deg': _degrees},
                                        scoring='mean_squared_error')
            elif method == 'cheb':
                # initialise optimal polinomial estimator:
                estimator = ChebyshevRegression()
                _degrees = np.arange(0, 10)
                cv_model = GridSearchCV(estimator,
                                        param_grid={'deg': _degrees},
                                        scoring='mean_squared_error')
            else:
                raise Exception('Unknown smoothing method. Use \'poly\' or \'cheb\'')
            # iterate over scans, as fits must be scan-based
            delay_smooth = []
            uvw_smooth = []
            doppler_smooth = []

            for start, stop in zip(self.scanStartTimes, self.scanStopTimes):
                dd0 = datetime.datetime(start.year, start.month, start.day)
                # time scale and proper indices to make fit
                time = np.array([
                                    (ii, t.hour * 3600 + t.minute * 60.0 + t.second +
                                     (t - dd0).days * 86400.0)
                                    for ii, t in enumerate(self.tstamps)
                                    if start <= t <= stop])

                # time scale used for smoothing:
                t_dense = np.arange(time[0, 1], time[-1, 1] + tstep, tstep)

                # renorm
                time[:, 1] = 24.0 * time[:, 1] / 86400.0
                t_dense = 24.0 * t_dense / 86400.0

                if len(self.dude.delay) > 0:
                    # make optimal fit to each of delay 'components'
                    delay_smooth_scan = []
                    for ii in range(self.dude.delay.shape[1]):
                        # time[:,0] - indices of current scan in full-len tstamps
                        # time[:,1] - time stamps in the flesh
                        ind = map(int, time[:, 0])
                        cv_model.fit(time[:, 1], self.dude.delay[ind, ii])
                        # print [grid_score.mean_validation_score for \
                        #        grid_score in cv_model.grid_scores_]
                        delay_smooth_scan.append(cv_model.predict(t_dense))
                    try:
                        delay_smooth = np.vstack((delay_smooth,
                                                  np.array(delay_smooth_scan).T))
                    except:
                        delay_smooth = np.array(delay_smooth_scan).T

                if len(self.dude.uvw) > 0:
                    uvw_smooth_scan = []
                    for ii in range(self.dude.uvw.shape[1]):
                        ind = map(int, time[:, 0])
                        cv_model.fit(time[:, 1], self.dude.uvw[ind, ii])
                        uvw_smooth_scan.append(cv_model.predict(t_dense))
                    try:
                        uvw_smooth = np.vstack((uvw_smooth, np.array(uvw_smooth_scan).T))
                    except:
                        uvw_smooth = np.array(uvw_smooth_scan).T

                if len(self.dude.doppler) > 0:
                    doppler_smooth_scan = []
                    for ii in range(self.dude.doppler.shape[1]):
                        ind = map(int, time[:, 0])
                        cv_model.fit(time[:, 1], self.dude.doppler[ind, ii])
                        doppler_smooth_scan.append(cv_model.predict(t_dense))
                    try:
                        doppler_smooth = np.vstack((doppler_smooth, np.array(doppler_smooth_scan).T))
                    except:
                        doppler_smooth = np.array(doppler_smooth_scan).T
                        #            print '___'
            self.dude.delay = delay_smooth
            self.dude.uvw = uvw_smooth
            self.dude.doppler = doppler_smooth

            # resample tstamps and freqs:
            self.resample(tstep=tstep)

'''
#==============================================================================
# Container for Delays, Uvws, Doppler and Etc.
#==============================================================================
'''
class dude(object):
    """
    the class contains values of Delays, Uvws, Doppler and Etc. 
    calculated for a given obs-object
    """
    def __init__(self):
        # initialise lists for delays, uvws, dopplers
        self.delay = []
        self.uvw = []
        self.doppler = []
#    def __init__(self, delay_calc=False, uvw_calc=False, doppler_calc=False):
#        self.delay_calc = delay_calc
#        self.uvw_calc = uvw_calc
#        self.doppler_calc = doppler_calc
#        # initialise lists for delays, uvws, dopplers
#        if self.delay_calc:
#            self.delay = []
#        if self.uvw_calc:
#            self.uvw = []
#        if self.doppler_calc:
#            self.doppler = []
            
'''
#==============================================================================
# Constants
#==============================================================================
'''
class constants(object):
    """
    Physical constants used in the code.
    
    IMPORTANT!!! TDB->TCB->TT
     to convert to TCB-compatible values, the appropriate TDB-compatible mass 
     value has to be multiplied by (1+L_B)
     then to get the TT-compatible values, the TCB-compatible mass should be
     multiplied by (1-L_G)
    """
    def __init__(self, jpl_eph='421'):
        # jpl_eph used can be 405 or 421, the latter is the default
        # math consts
        self.CDEGRAD = 1.7453292519943296e-02
        self.CARCRAD = 4.8481368110953599e-06
        self.CTIMRAD = 7.2722052166430399e-05
        self.SECDAY = 86400.0
        self.JUL_CENT = 36525.0
        # JULIAN DATE OF STANDARD EPOCH J2000.0
        self.JD2000 = 2451545.0
        # Algemeene Physical constants
        self.C = 2.99792458e8
        self.C_km = 2.99792458e5
#        self.F = 298.25765 # the tide-free value
        self.F = 298.25642 # the tide-free value, IERS2010
        self.AE = 6378136.3 # the tide-free value
        self.J_2 = 1.0826e-3;
        self.AU = 149597870691.0
        self.TAUA = 499.0047838061
        self.G = 6.67428e-11
        
        self.L_B = 1.550519768e-8 #tn36
        self.L_C = 1.48082686741e-8
        self.L_G = 6.969290134e-10
        
        ## DE/LE405 Header. TDB-compatible!!
        if '403' in jpl_eph:
            AU_DE405 = 1.49597870691000015e+11 #m
            self.GSUN = 0.295912208285591095e-03*(AU_DE405)**3/(86400.0)**2
            self.MU = (0.813005600000000044e+02)**(-1)
            self.GEARTH = 0.899701134671249882e-09*(AU_DE405)**3/(86400.0)**2 / (1+self.MU)
            self.GMOON = self.GEARTH*self.MU
            self.GMPlanet = [0.491254745145081187e-10*(AU_DE405)**3/(86400.0)**2,
                             0.724345248616270270e-09*(AU_DE405)**3/(86400.0)**2,
                             0.954953510577925806e-10*(AU_DE405)**3/(86400.0)**2,
                             0.282534590952422643e-06*(AU_DE405)**3/(86400.0)**2,
                             0.845971518568065874e-07*(AU_DE405)**3/(86400.0)**2,
                             0.129202491678196939e-07*(AU_DE405)**3/(86400.0)**2,
                             0.152435890078427628e-07*(AU_DE405)**3/(86400.0)**2,
                             0.218869976542596968e-11*(AU_DE405)**3/(86400.0)**2]
        if '405' in jpl_eph:
            AU_DE405 = 1.49597870691000015e+11 #m
            self.GSUN = 0.295912208285591095e-03*(AU_DE405)**3/(86400.0)**2
            self.MU = (0.813005600000000044e+02)**(-1)
            self.GEARTH = 0.899701134671249882e-09*(AU_DE405)**3/(86400.0)**2 / (1+self.MU)
            self.GMOON = self.GEARTH*self.MU
            self.GMPlanet = [0.491254745145081187e-10*(AU_DE405)**3/(86400.0)**2,
                             0.724345248616270270e-09*(AU_DE405)**3/(86400.0)**2,
                             0.954953510577925806e-10*(AU_DE405)**3/(86400.0)**2,
                             0.282534590952422643e-06*(AU_DE405)**3/(86400.0)**2,
                             0.845971518568065874e-07*(AU_DE405)**3/(86400.0)**2,
                             0.129202491678196939e-07*(AU_DE405)**3/(86400.0)**2,
                             0.152435890078427628e-07*(AU_DE405)**3/(86400.0)**2,
                             0.218869976542596968e-11*(AU_DE405)**3/(86400.0)**2]
        ## DE/LE421 Header. TDB-compatible!!
        if '421' in jpl_eph:
            AU_DE421 = 1.49597870699626200e+11 #m
            self.GSUN = 0.295912208285591100e-03*(AU_DE421)**3/(86400.0)**2
            self.MU = (0.813005690699153000e+02)**(-1)
            self.GEARTH = 0.899701140826804900e-09*(AU_DE421)**3/(86400.0)**2 / (1+self.MU)
            self.GMOON = self.GEARTH*self.MU
            self.GMPlanet = [0.491254957186794000e-10*(AU_DE421)**3/(86400.0)**2,
                             0.724345233269844100e-09*(AU_DE421)**3/(86400.0)**2,
                             0.954954869562239000e-10*(AU_DE421)**3/(86400.0)**2,
                             0.282534584085505000e-06*(AU_DE421)**3/(86400.0)**2,
                             0.845970607330847800e-07*(AU_DE421)**3/(86400.0)**2,
                             0.129202482579265000e-07*(AU_DE421)**3/(86400.0)**2,
                             0.152435910924974000e-07*(AU_DE421)**3/(86400.0)**2,
                             0.217844105199052000e-11*(AU_DE421)**3/(86400.0)**2]
        ## INPOP13c Header. TDB-compatible!!
        if '13c' in jpl_eph:
            AU_13c = 1.495978707000000e+11
            self.GSUN = 0.2959122082912712e-03*(AU_13c)**3/(86400.0)**2
            self.MU = (0.8130056945994197e+02)**(-1)
            self.GEARTH = 0.8997011572788968e-09*(AU_13c)**3/(86400.0)**2 / (1+self.MU)
            self.GMOON = self.GEARTH*self.MU
            self.GMPlanet = [0.4912497173300158e-10*(AU_13c)**3/(86400.0)**2,
                             0.7243452327305554e-09*(AU_13c)**3/(86400.0)**2,
                             0.9549548697395966e-10*(AU_13c)**3/(86400.0)**2,
                             0.2825345791109909e-06*(AU_13c)**3/(86400.0)**2,
                             0.8459705996177680e-07*(AU_13c)**3/(86400.0)**2,
                             0.1292024916910406e-07*(AU_13c)**3/(86400.0)**2,
                             0.1524357330444817e-07*(AU_13c)**3/(86400.0)**2,
                             0.2166807318808926e-11*(AU_13c)**3/(86400.0)**2]

        
        self.TDB_TCB = (1.0+self.L_B) # F^-1
        # G*masses in TCB-frame!
        self.GM_TCB = np.hstack([self.GMPlanet[0:2], self.GEARTH,\
                        self.GMPlanet[2:], self.GMOON, self.GSUN]) * self.TDB_TCB
        # G*masses in TDB-frame!
        self.GM = np.hstack([self.GMPlanet[0:2], self.GEARTH,\
                        self.GMPlanet[2:], self.GMOON, self.GSUN])

'''
#==============================================================================
# 
#==============================================================================
'''
class site(object):
    '''
    Class containing site-specific data for a station:
        - station name
        - position
        - velocity
        - thermal def coeffs and antenna axis offsets
        - ocean loading parameters
        - atmospheric loading parameters
        - azimuth/elevation angles in case of a S/C obs
    '''
    def __init__(self, name):
        self.name = name.strip()
        # crds in terrestrial RF at J2000.0:
        self.r_GTRS = np.zeros(3)
        self.v_GTRS = np.zeros(3)
        # crds in terrestrial RF at the epoch of obs (acc for tectonic plate motion):
        self.r_GTRS_date = np.zeros(3)
        # crds in celestial geocentric RF at J2000.0:
        self.r_GCRS = np.zeros(3)
        self.v_GCRS = np.zeros(3)
        self.a_GCRS = np.zeros(3)
        # crds in celestial barycentric RF at J2000.0:
        self.r_BCRS = np.zeros(3)
        self.v_BCRS = np.zeros(3)
        self.a_BCRS = np.zeros(3)
        # geodetic data:
        self.sph_rad = 0.0 # site spherical radius
        self.lat_gcen = 0.0 # site spherical radius
        self.lon_gcen = 0.0 # site spherical radius
        self.lat_geod = 0.0 # site spherical radius
        self.h_geod = 0.0 # site spherical radius
        # distance from the Earth spin axis and from the equatorial plane (in km):
        self.u = 0.0
        self.v = 0.0
        # VEN-to-crust-fixed rotation matrix:
        self.vw = np.identity(3)
        self.eta = np.zeros(3)
        self.theta = 0.0
        self.R_E = np.zeros(3) # Pierce point with the Earth ellipsoid
        # antennae Thermal deformations coefficients + some auxiliary info, 
        # e.g. mount type:
        self.ivs_name = ''
        self.focus_type = ''
        self.mount_type = ''
        self.radome = ''
        self.meas_type = ''
        self.T0 = 0.0
        self.sin_T = 0.0
        self.cos_T = 0.0
        self.h0 = 0.0
        self.ant_diam = 0.0
        self.hf = 0.0
        self.df = 0.0
        self.gamma_hf = 0.0
        self.hp = 0.0
        self.gamma_hp = 0.0
        self.AO = 0.0
        self.gamma_AO = 0.0
        self.hv = 0.0
        self.gamma_hv = 0.0
        self.hs = 0.0
        self.gamma_hs = 0.0
        # ocean loading coeffitients:
        self.amp_ocean = np.zeros((11,3))
        self.phs_ocean = np.zeros((11,3))
        # Atmospheric loading displacement:
        # !!! not implemented!
        if 1==0:
            self.amp_atm = np.zeros((9,3))
        # azimuth/elevation of a S/C or a source
#        self.azel = np.zeros(2)
        self.azel_interp = []
        # meteo data, stored in a dictionary:
        self.met = {'mjd':[],'ahz':[], 'awz':[], 'zhz':[], 'zwz':[], 'TC':[], \
                   'pres':[], 'hum':[], 'gnh':[], 'geh':[], 'gnw':[], 'gew':[]}
#        print self.met
        # interpolants for meteo data:
        self.fMet = {'fT':[], 'fP':[], 'fH':[],\
                     'fAhz':[], 'fAwz':[], 'fZhz':[], 'fZwz':[],\
                     'fGnh':[], 'fGeh':[], 'fGnw':[], 'fGew':[]}
        ## displacements in GCRS due to different phenomena:
        # solid Earth tides:
        self.dr_tide = np.zeros(3)
        self.dv_tide = np.zeros(3)
        # ocean loading:
        self.dr_oclo = np.zeros(3)
        self.dv_oclo = np.zeros(3)
        # pole tide:
        self.dr_poltide = np.zeros(3)
        self.dv_poltide = np.zeros(3)
        # atmospheric loading:
        self.dr_atlo = np.zeros(3)
        self.dv_atlo = np.zeros(3)
        ## delays due to instrumental and propagation effects
        # thermal deformation of telescope:
        self.dtau_therm = 0.0
        # axis offset:
        self.dtau_ao = 0.0
        # troposphere:
        self.dtau_tropo = 0.0
        # ionosphere:
        self.dtau_iono = 0.0


    def geodetic(self, const):
        
        '''
        Calculate the site position in geodetic coordinate systems
        for the stations participating in the current observation.
        Transformation VW from local geodetic coordinate system (Vertical,East,North)
        to the Earth-fixed coordinate system is calculated for each cite.
         
        For transformation to geodetic coordinate system the 
        Reference Ellipsoid from Table 1.1, IERS 2010 Conventions is used.
        The geophysical values are those for the "zero-frequency" tide system.
        NOTE !!! The ITRF2000 is "conventional tide free crust" system.
      
         Input  - 
           1. sta               - list of site-objects
           2. const             - physical constants
           
         Output to each site-object:
           lat_geod          - The geodetic latitude at each site (RAD)
           h_geod            - The geodetic height of each site. (M)
           lat_gcen          - The geocentric latitude at each site. (RAD)
           lon_gcen          - The geocentric east longitude at each site. (RAD)
                                  (0 <= lon_gcen <= 2*pi)
           sph_rad           - The site spherical radius (M)
           u                 - The stations distance from the Earth spin axis (KM)
           v                 - The stations distance from the equatorial plane(KM)
           vw                - The transformation matrix of the VEN geodetic
                                  system to the Earth-fixed coordinate system
        '''
        # IERS 2010
        AE = 6378136.3 # m
        F = 298.25642
        
        if self.name == 'RA' or self.name == 'GEOCENTR':
            return

        # Compute the site spherical radius
        self.sph_rad = sqrt( sum(x*x for x in self.r_GTRS_date) )
        # Compute geocentric latitude
        self.lat_gcen = asin( self.r_GTRS_date[2] / self.sph_rad )
        # Compute geocentric longitude       
        self.lon_gcen = atan2( self.r_GTRS_date[1], self.r_GTRS_date[0] )
        if self.lon_gcen < 0.0:
            self.lon_gcen = self.lon_gcen + 2.0*pi
        # Compute the stations distance from the Earth spin axis and
        # from the equatorial plane (in KM)
        req = sqrt( sum(x*x for x in self.r_GTRS_date[:-1]) )
        self.u = req*1e-3
        self.v = self.r_GTRS_date[2]*1e-3
        # Compute geodetic latitude and height.
        # The geodetic longitude is equal to the geocentric longitude
        self.lat_geod, self.h_geod = self.geoid( req, self.r_GTRS_date[2], \
                                                 AE, F )
        # Compute the local VEN-to-crust-fixed rotation matrices by rotating
        # about the geodetic latitude and the longitude. 
        # w - rotation matrix by an angle lat_geod around the y axis
        w = self.R_123(2, self.lat_geod)
        # v - rotation matrix by an angle -lon_gcen around the z axis
        v = self.R_123(3, -self.lon_gcen)
        # product of the two matrices:        
        self.vw = np.dot(v, w)    
        
        
        # WGS84 Earth ellipsoid:
#        a = 6378.1370 # km semi-major axis
#        b = a*(1 - 1/298.257223563) # km semi-minor axis
        # f = (a-b)/a = 1/298.257223563
        # IERS 2010 Ellipsoid:
        a = AE
        b = a*(1-1/F)
        e = np.sqrt(a**2-b**2)/a
        
        R = self.r_GTRS_date
        # initial approximation. Formula (43a), p. 501
        sq = np.sqrt( R[0]**2 + R[1]**2 + (1-e**2)*(R[2]**2) )
        R_E = np.array([a*R[0]/sq, a*R[1]/sq, a*R[2]*(1-e**2)/sq])

        # iterative process to solve for the Earth diameter at the site
        R_E_old = np.zeros(3) # initialise
        while norm(R_E-R_E_old) > 1e-9:
            R_E_old = np.copy(R_E) #for comparison
            # (21a), p.498
            r_s_shtrih = np.array([R[0]-R_E[0]*(e**2), R[1]-R_E[1]*(e**2), R[2]])
            # (21a), p.498
            R_s_shtrih = np.array([R_E[0]*(1-e**2), R_E[1]*(1-e**2), R_E[2]])
            # (22a), p.498
            eta = r_s_shtrih / norm(r_s_shtrih)
            # (32), p.499
            theta = R_s_shtrih[0]/r_s_shtrih[0]

            # new XeYeZe, (33a), p.499
            R_E = np.array([theta*(R[0] - (e**2)*R_E[0])/(1-e**2),\
                            theta*(R[1] - (e**2)*R_E[1])/(1-e**2), theta*R[2]])
        self.eta = eta
        self.theta = theta
        self.R_E = R_E

    @staticmethod
    def geoid(r, z, a, fr):
        '''
        Transform Cartesian to geodetic coordinates
        based on the exact solution (Borkowski,1989)
         
               Input variables :  
                     r, z = equatorial [m] and polar [m] components
               Output variables: 
                     fi, h = geodetic coord's (latitude [rad], height [m])
        
        IERS ellipsoid: semimajor axis (a) and inverse flattening (fr)
        '''
    
        if z>=0.0:
            b = abs(a - a/fr)
        else:
            b = -abs(a - a/fr)
        E = ((z + b)*b/a - a)/r
        F = ((z - b)*b/a + a)/r
        # Find solution to: t**4 + 2*E*t**3 + 2*F*t - 1 = 0
        P = (E*F + 1.0)*4.0/3.0;
        Q = (E*E - F*F)*2.0
        D = P*P*P + Q*Q
        if D >= 0.0:
            s = sqrt(D) + Q
            if s>=0:
                s = abs(exp(log(abs(s))/3.0))
            else:
                s = -abs(exp(log(abs(s))/3.0))
            v = P/s - s
            # Improve the accuracy of numeric values of v
            v = -(Q + Q + v*v*v)/(3.0*P)
        else:
            v = 2.0*sqrt(-P)*cos(acos(Q/P/sqrt(-P))/3.0)
    
        G = 0.5*(E + sqrt(E*E + v))
        t = sqrt(G*G + (F - v*G)/(G + G - E)) - G
        fi = atan((1.0 - t*t)*a/(2.0*b*t))
        h = (r - a*t)*cos(fi) + (z - b)*sin(fi)
    
        return fi, h
    
    @staticmethod
    def R_123 (i, theta):
        '''
        function R_123 creates a matrix 'R_i' which describes
        a right rotation by an angle 'THETA' about coordinate axis 'I'
      
         Input variables:
            1. I      -  The number which determines the rotation axis.
                         (I = 1, 2, 3 corresponds X, Y, and Z axes respectfully)
            2. THETA  -  The rotation angle. (RAD)
      
         Output variables:
            1. R(3,3) -  The 3x3 rotation matrix. (UNITLESS)
        '''
        r = np.zeros((3,3))
        c = cos(theta)
        s = sin(theta)
    
        if i==1:
        #Rotation around the X-axis:
        #        ( 1  0  0 )
        # R(X) = ( 0  C  S )
        #        ( 0 -S  C )
            r[0,0] = 1.0
            r[1,0] = 0.0
            r[2,0] = 0.0
            r[0,1] = 0.0
            r[1,1] = +c
            r[2,1] = -s
            r[0,2] = 0.0
            r[1,2] = +s
            r[2,2] = +c     
        elif i==2:
        #Rotation around the Y-axis:
        #        ( C  0 -S )
        # R(Y) = ( 0  1  0 )
        #        ( S  0  C )
            r[0,0] = +c
            r[1,0] = 0.0
            r[2,0] = +s
            r[0,1] = 0.0
            r[1,1] = 1.0
            r[2,1] = 0.0
            r[0,2] = -s
            r[1,2] = 0.0
            r[2,2] = +c
        elif i==3:
        #Rotation around the Z-axis:
        #        ( C  S  0 )
        # R(Z) = (-S  C  0 )
        #        ( 0  0  1 ) 
            r[0,0] = +c
            r[1,0] = -s
            r[2,0] = 0.0
            r[0,1] = +s
            r[1,1] = +c
            r[2,1] = 0.0
            r[0,2] = 0.0
            r[1,2] = 0.0
            r[2,2] = 1.0
    
        return r
    
    def AzEl2(self, gcrs, utc, JD, t_1, r2000, jpl_eph):
        '''
        Calculate topocentric [Elevation, Azimuth] of the S/C
        
        input: 
            gcrs - eph.gcrs
            utc  - eph.UT in decimal days
            t_1  - obs epoch in decimal days
        '''
        const = constants()
        C = const.C # m/s
        GM = const.GM
        
        precision = 1e-12
        n_max = 3
        lag_order = 5
        
        # initial approximation:
        nn = 0
        lt_01_tmp = 0.0
    
        astropy_t_1 = Time(JD + t_1/86400.0, format='jd', scale='utc', precision=9)
        eph_t_0 = datetime.datetime(*map(int, gcrs[0,:3]))
        # first time stamp negative?
        if utc[0]<0:
            eph_t_0 += datetime.timedelta(days=1)
        # cut eph and it's tomorrow?
        if utc[0]//1 > 0:
            eph_t_0 -= datetime.timedelta(days=utc[0]//1)
        dd = (astropy_t_1.datetime - eph_t_0).days
        
#        print utc, JD, t_1, dd
        
        x, _ = lagint(lag_order, utc, gcrs[:,6], t_1+dd)
        y, _ = lagint(lag_order, utc, gcrs[:,7], t_1+dd)
        z, _ = lagint(lag_order, utc, gcrs[:,8], t_1+dd)
        R_0_0 = np.hstack((x,y,z))

        lt_01 = norm(R_0_0 - self.r_GCRS)/C
        
        
        mjd = JD - 2400000.5
        astropy_t_0 = Time(mjd + t_1 - lt_01/86400.0, \
                            format='mjd', scale='utc', precision=9)
        t_0_0 = astropy_t_0.tdb.jd2

        ''' BCRS state vectors of celestial bodies at JD+CT, [m, m/s]: '''
        ## Earth:
        rrd = pleph(JD+t_0_0, 3, 12, jpl_eph)
        earth = np.reshape(np.asarray(rrd), (3,2), 'F') * 1e3
        ## Sun:
        rrd = pleph(JD+t_0_0, 11, 12, jpl_eph)
        sun = np.reshape(np.asarray(rrd), (3,2), 'F') * 1e3
        ## Moon:
        rrd = pleph(JD+t_0_0, 10, 12, jpl_eph)
        moon = np.reshape(np.asarray(rrd), (3,2), 'F') * 1e3
    
        state_ss = []
        for jj in (1,2,4,5,6,7,8,9):
            rrd = pleph(JD+t_0_0, jj, 12, jpl_eph)
            state_ss.append(np.reshape(np.asarray(rrd), (3,2), 'F') * 1e3)
        state_ss.insert(2, earth)
        state_ss.append(moon)
        state_ss.append(sun)

        while (abs(lt_01 - lt_01_tmp) > precision) and (nn < n_max):
            lt_01_tmp = lt_01
            t_0 = t_1 - lt_01/86400.0
#            print 't_0 =', t_0
            x, _ = lagint(lag_order, utc, gcrs[:,6], t_0+dd)
            y, _ = lagint(lag_order, utc, gcrs[:,7], t_0+dd)
            z, _ = lagint(lag_order, utc, gcrs[:,8], t_0+dd)
            vx, _ = lagint(lag_order, utc, gcrs[:,9], t_0+dd)
            vy, _ = lagint(lag_order, utc, gcrs[:,10], t_0+dd)
            vz, _ = lagint(lag_order, utc, gcrs[:,11], t_0+dd)
            R_0 = np.hstack((x,y,z))
            V_0 = np.hstack((vx,vy,vz))
            
            # vector needed for RLT calculation
#            R_01 = -(self.r_GCRS - R_0)  ## WTF, Mityaj???
            R_01 = self.r_GCRS - R_0
            
            RLT = 0.0
            for ii, state in enumerate(state_ss):
                if ii==2 and norm(self.r_GCRS)==0.0: continue
                rb = state[:,0]
                vb = state[:,1]
                R_0_B  = R_0 - (rb - (t_0-t_0_0)*86400.0*vb)
                R_1_B  = self.r_GCRS - rb
                R_01_B = R_1_B - R_0_B

                RLT += (2.0*GM[ii]/C**3) * \
                      log(  ( norm(R_0_B) + norm(R_1_B) + norm(R_01_B) + \
                              2.0*GM[ii]/C**2 ) / \
                            ( norm(R_0_B) + norm(R_1_B) - norm(R_01_B) + \
                              2.0*GM[ii]/C**2 ) )

            lt_01 = lt_01 - (lt_01 - norm(R_01)/C - RLT) / \
                        ( 1.0 - np.dot(R_01, V_0)/(C*norm(R_01)) )
            
#            RLT = 0.0
#
#            lt_01 = lt_01 - (lt_01 - norm(R_01)/C - RLT) / \
#                        ( 1.0 - np.dot(R_01, V_0)/(C*norm(R_01)) )

            t_0 = t_1 - lt_01/86400.0

            nn += 1

        x, _ = lagint(lag_order, utc, gcrs[:,6], t_0+dd)
        y, _ = lagint(lag_order, utc, gcrs[:,7], t_0+dd)
        z, _ = lagint(lag_order, utc, gcrs[:,8], t_0+dd)
        # it's still station-centric!
        r = -(self.r_GCRS - np.hstack((x,y,z)))

        # S/C position is given at a moment LT seconds ago, which
        # means r is abberated in the far-field case sense

        # normalise aberrated vector:
        K_unit_aber = r / norm(r)
        
        # compute the rotation matrix which rotates from the geocentric 
        # crust fixed system to the VEN system
        
        # Rotate the aberrated vector to the crust fixed system:
        Crust_star = np.dot(r2000.T, K_unit_aber)
        ven_star = np.dot(self.vw.T, Crust_star)
        
        el = np.arcsin( ven_star[0] )
        az = np.arctan2(ven_star[1], ven_star[2])
        if az < 0.0:
            az += 2.0*np.pi
        
        return az, el, lt_01
    
    
    def AzEl(self, gtrs, JD, UTC, t_obs, jpl_eph, interpolants=False):
        '''
        Calculate topocentric [Elevation, Azimuth] of the S/C
        on the basis of lt corrected GC SC ephemeris
         input:
             gtrs/UT - target ITRF ephemeris
             [t_obs] - if this is set, method returns Az/El at epoch
                       otherwise, constructs Az/El series corresponding
                       to eph.gtrs
                       t_obs should be in decimal days
        '''
        raise Exception('Depricated. Use AzEl2 instead.')
        
        if self.name=='GEOCENTR' or self.name=='RA':
            return

        const = constants()
        C = const.C
        GM = const.GM
        R_1 = self.r_GTRS_date
        UT = np.copy(UTC)
        UT *= 86400.0
        
        # check input, create list with epochs:
        try:
            len(t_obs) # this will fail if t_obs is a number
            epochs = 86400.0*np.array(t_obs)
        except:
            epochs = np.array([86400.0*t_obs])
        
        # compute LTs from the target to the station,
        # then correct the ITRF eph for these LTs using linear interpolation
        n_max = 3
        lag_order = 5
        
        r_lt = [] # lt-to-station-corrected s/c positions
        
        mjd = JD - 2400000.5
        
        for jj, t_1 in enumerate(epochs):
#            tic = _time()
#            R_0_0 = gtrs[jj,6:9]
#            V_0_0 = gtrs[jj,9:12]
#            A_0_0 = gtrs[jj,12:15]
            # dummy acceleration=0 for RA/GNSS (lt is too short to bother..):
#            if len(A_0_0)==0: A_0_0 = np.zeros(3)

            x, _ = lagint(3, UT, gtrs[:,6], t_1)
            y, _ = lagint(3, UT, gtrs[:,7], t_1)
            z, _ = lagint(3, UT, gtrs[:,8], t_1)
            R_0_0 = np.hstack((x,y,z))            
            lt_01 = norm(R_1 - R_0_0)/C
            
            astropy_t_0 = Time(mjd + (t_1 - lt_01)/86400.0, \
                format='mjd', scale='utc', precision=9)
            t_0_0 = astropy_t_0.tdb.jd2
            
            ''' BCRS state vectors of celestial bodies at JD+CT, [m, m/s]: '''
            ## Earth:
            rrd = pleph(JD+t_0_0, 3, 12, jpl_eph)
            earth = np.reshape(np.asarray(rrd), (3,2), 'F') * 1e3
            ## Sun:
            rrd = pleph(JD+t_0_0, 11, 12, jpl_eph)
            sun = np.reshape(np.asarray(rrd), (3,2), 'F') * 1e3
            ## Moon:
            rrd = pleph(JD+t_0_0, 10, 12, jpl_eph)
            moon = np.reshape(np.asarray(rrd), (3,2), 'F') * 1e3
        
            state_ss = []
            for jj in (1,2,4,5,6,7,8,9):
                rrd = pleph(JD+t_0_0, jj, 12, jpl_eph)
                state_ss.append(np.reshape(np.asarray(rrd), (3,2), 'F') * 1e3)
            state_ss.insert(2, earth)
            state_ss.append(moon)
            state_ss.append(sun)
            
            numlt = max(15.0, 2.1*lt_01)
            mini = max(np.searchsorted(UT, t_1-numlt)-1, 0)
            maxi = np.searchsorted(UT, t_1+numlt)

            # make n_max iterations
            for nn in range(n_max):
                t_0 = t_1 - lt_01 # in sec

#                R_0 = R_0_0 - lt_01*V_0_0 + (lt_01)**2*A_0_0/2.0
#                V_0 = V_0_0 - lt_01*A_0_0
                x, _ = lagint(lag_order, UT[mini:maxi], gtrs[mini:maxi,6], t_0)
                y, _ = lagint(lag_order, UT[mini:maxi], gtrs[mini:maxi,7], t_0)
                z, _ = lagint(lag_order, UT[mini:maxi], gtrs[mini:maxi,8], t_0)
                R_0 = np.hstack((x,y,z))
                vx, _ = lagint(lag_order, UT[mini:maxi], gtrs[mini:maxi,9], t_0)
                vy, _ = lagint(lag_order, UT[mini:maxi], gtrs[mini:maxi,10], t_0)
                vz, _ = lagint(lag_order, UT[mini:maxi], gtrs[mini:maxi,11], t_0)
                V_0 = np.hstack((vx,vy,vz))

                # vector needed for RLT calculation
                R_01 = R_1 - R_0
                
                RLT = 0.0
                for ii, state in enumerate(state_ss):
                    if ii==2 and norm(R_1)==0.0: continue
                    rb = state[:,0]
                    vb = state[:,1]
                    R_0_B  = R_0 - (rb - (t_0-t_0_0*86400.0)*vb)
                    R_1_B  = R_1 - rb
                    R_01_B = R_1_B - R_0_B

                    RLT += (2.0*GM[ii]/C**3) * \
                          log(  ( norm(R_0_B) + norm(R_1_B) + norm(R_01_B) + \
                                  2.0*GM[ii]/C**2 ) / \
                                ( norm(R_0_B) + norm(R_1_B) - norm(R_01_B) + \
                                  2.0*GM[ii]/C**2 ) )

                lt_01 = lt_01 - (lt_01 - norm(R_01)/C - RLT) / \
                            ( 1.0 - np.dot(R_01, V_0)/(C*norm(R_01)) )
                t_0 = t_1 - lt_01

                nn += 1
            
#            R_0 = R_0_0 - lt_01*V_0_0 + (lt_01**2)*A_0_0/2.0
            x, _ = lagint(lag_order, UT[mini:maxi], gtrs[mini:maxi,6], t_0)
            y, _ = lagint(lag_order, UT[mini:maxi], gtrs[mini:maxi,7], t_0)
            z, _ = lagint(lag_order, UT[mini:maxi], gtrs[mini:maxi,8], t_0)
            R_0 = np.hstack((x,y,z))
            r_lt.append(R_0)
#            print lt_01
#            print _time()-tic
#            
#            print '__'
        
        r_lt = np.array(r_lt)
        
#        np.set_printoptions(precision=16)
#        print r_lt[0]
#        print R_1

        ## iterative eta + (phi,lambda,h) determination
#        print 'uuu'
#        tic = _time()
        eta = self.eta
        
        ## (H, A) PROPER CALCULATION
        # elevation
        rho = r_lt - R_1
        H = []
        for rh in rho:
            H.append(np.arcsin( np.dot(rh, eta) / ( norm(rh) ) ))
        H = np.array(H).T
        # azimuth from North to East
        Z = np.array([0.0, 0.0, 1.0])
        E = np.cross(Z, R_1) / norm(np.cross(Z, R_1))
        N = np.cross(eta, E) / norm(np.cross(eta, E))
        A = []
        for rh in rho:
            A.append(np.arctan2(np.dot(rh, E), np.dot(rh, N)))
            if A[-1]<0.0: A[-1] += 2.0*np.pi
        A = np.array(A).T
#        self.azel = np.array(np.vstack((A, H))).T
        
#        print A, H
#        raw_input()
        if not interpolants:
            if len(epochs)==1:
                return A[0], H[0]
            else:
                return A, H
        else:
            # store interpolant of these - for faster access later
            fA = sp.interpolate.interp1d(epochs/86400.0, A)
            fH = sp.interpolate.interp1d(epochs/86400.0, H)
            self.azel_interp = [fA, fH]
#        print _time()-tic
        
    def LT_radec_bc(self, bcrs, tdb, JD_UTC, t_1_UTC, jpl_eph):
        """
        Calculate station-centric LTs and LT-corrected ra/decs
        
        input: 
            bcrs - eph.bcrs
            tdb  - eph.CT in decimal days
            JD_UTC, t_1_UTC  - obs epoch in decimal days, TDB scale
        """
        const = constants()
        C = const.C  # m/s
        GM = const.GM
        
        precision = 1e-12
        n_max = 3
        lag_order = 5

        # zero point
        astropy_t_1 = Time(JD_UTC, t_1_UTC, format='jd', scale='utc', precision=9).tdb
        mjd_full = astropy_t_1.mjd
        mjd = np.floor(mjd_full)
        # t_1_coarse = mjd_full - mjd
        JD = mjd + 2400000.5
        astropy_t_1_day_start = Time(JD, format='jd', scale='tdb', precision=9)
        t_1 = (astropy_t_1-astropy_t_1_day_start).jd
        t_start_day = datetime.datetime(*map(int, bcrs[0, :3]))
        dd = (astropy_t_1.datetime - t_start_day).total_seconds() // 86400

        # print astropy_t_1.tdb.datetime, t_start_day, dd, JD, t_1
        # initial approximation:
        nn = 0
        lt_01_tmp = 0.0

        x, _ = lagint(lag_order, tdb, bcrs[:, 6], t_1+dd)
        y, _ = lagint(lag_order, tdb, bcrs[:, 7], t_1+dd)
        z, _ = lagint(lag_order, tdb, bcrs[:, 8], t_1+dd)
        R_0_0 = np.hstack((x, y, z))
        
        # check if self.r_BCRS is set
        ## Earth:
        rrd = pleph(JD + t_1, 3, 12, jpl_eph)
        earth = np.reshape(np.asarray(rrd), (3,2), 'F') * 1e3
        # low accuracy is good enough for this application:
        self.r_BCRS = earth[:, 0] + self.r_GCRS

        lt_01 = norm(R_0_0 - self.r_BCRS)/C

        ''' BCRS state vectors of celestial bodies at JD+CT, [m, m/s]: '''
        state_ss = []
        for jj in (1, 2, 4, 5, 6, 7, 8, 9, 10, 11):
            rrd = pleph(astropy_t_1.jd, jj, 12, jpl_eph)
            state_ss.append(np.reshape(np.asarray(rrd), (3, 2), 'F') * 1e3)
        
        while (abs(lt_01 - lt_01_tmp) > precision) and (nn < n_max):
            lt_01_tmp = lt_01
            t_0 = t_1 - lt_01/86400.0
            # print t_0

            x, _ = lagint(lag_order, tdb, bcrs[:, 6], t_0+dd)
            y, _ = lagint(lag_order, tdb, bcrs[:, 7], t_0+dd)
            z, _ = lagint(lag_order, tdb, bcrs[:, 8], t_0+dd)
            vx, _ = lagint(lag_order, tdb, bcrs[:, 9], t_0+dd)
            vy, _ = lagint(lag_order, tdb, bcrs[:, 10], t_0+dd)
            vz, _ = lagint(lag_order, tdb, bcrs[:, 11], t_0+dd)
            R_0 = np.hstack((x, y, z))
            V_0 = np.hstack((vx, vy, vz))
            
            # vector needed for RLT calculation
#            R_01 = -(self.r_GCRS - R_0)  ## WTF, Mityaj???
            R_01 = self.r_BCRS - R_0
            
            RLT = 0.0
            
            for j, ii in enumerate((1, 2, 4, 5, 6, 7, 8, 9, 10, 11)):
                rrd = pleph(JD + t_0, ii, 12, jpl_eph)
                state = np.reshape(np.asarray(rrd), (3, 2), 'F') * 1e3
                R_B = state[:,0]

                R_0_B  = R_B - R_0
                R_1_B  = state_ss[j][:,0] - self.r_BCRS
                R_01_B = R_1_B - R_0_B
#                print ii, R_0, rb, vb, t_0, t_0_0
                RLT += (2.0*GM[ii-1]/C**3) * \
                    log(  ( norm(R_0_B) + norm(R_1_B) + norm(R_01_B) + 2.0*GM[ii-1]/C**2 ) /
                        ( norm(R_0_B) + norm(R_1_B) - norm(R_01_B) + 2.0*GM[ii-1]/C**2 ) )

            lt_01 = lt_01 - (lt_01 - norm(R_01)/C - RLT) / \
                        ( 1.0 - np.dot(R_01, V_0)/(C*norm(R_01)) )

            t_0 = t_1 - lt_01/86400.0

            nn += 1

        x, _ = lagint(lag_order, tdb, bcrs[:, 6], t_0+dd)
        y, _ = lagint(lag_order, tdb, bcrs[:, 7], t_0+dd)
        z, _ = lagint(lag_order, tdb, bcrs[:, 8], t_0+dd)
        # it's still station-centric!
        r = -(self.r_BCRS - np.hstack((x, y, z)))

        ra = np.arctan2(r[1], r[0])  # right ascension
        dec = np.arctan(r[2]/np.sqrt(r[0]**2+r[1]**2))  # declination
        if ra < 0:
            ra += 2.0*np.pi
#        print ra, dec
#        print ra*12/np.pi, dec*180/np.pi
#        raw_input()
#        self.lt = lt_01
#        self.ra = ra
#        self.dec = dec

        return lt_01, ra, dec


    def addMet(self, date_start, date_stop, inp):
        '''
        Load VMF1 site/gridded meteo data. Files with these should
        be pre-downloaded with doup() function.
        
         input:
             date_start - datetime object with start date
             date_stop - datetime object with start date
             inp        - input settings (cats etc.)
        '''
        if self.name=='GEOCENTR' or self.name=='RA':
            return
        
        day_start = datetime.datetime(date_start.year,date_start.month,\
                                      date_start.day)
        day_stop = datetime.datetime(date_stop.year,date_stop.month,\
                                    date_stop.day) + datetime.timedelta(days=1)
        dd = (day_stop - day_start).days
        
        # make list with datetime objects ranging from 1st day to last+1
        dates = [day_start]
        for d in range(1,dd+1):
            dates.append(day_start + datetime.timedelta(days=d))

        for day in dates:
            year = day.year
            doy = day.timetuple().tm_yday # day of year
            vmf_file = '{:4d}{:03d}.vmf1_r'.format(year, doy)
            # vmf1 site files:
            with open(inp['meteo_cat']+'/'+vmf_file,'r') as f:
                f_lines = f.readlines()
            
            # if met data are present in the vmf1 site file, load 'em
            if any(self.name in s.split() for s in f_lines):
                noSiteData = False
                for f_line in f_lines:                    
#                    if self.name in f_line.split(' '):
                    if self.name in f_line.split():
                        tmp = [float(i) for i in f_line[10:].split()]
                        self.met['mjd'].append(tmp[0])
                        self.met['ahz'].append(tmp[1])
                        self.met['awz'].append(tmp[2])
                        self.met['zhz'].append(tmp[3])
                        self.met['zwz'].append(tmp[4])
                        self.met['pres'].append(tmp[6])
                        self.met['TC'].append(tmp[7])
                        wvp = tmp[8]
                        TK = tmp[7] + 273.15
                        ew = 10**(10.79574*(1-273.16/TK) - \
                             5.028*np.log10(TK/273.16) + \
                             1.50475e-4*(1-10**(8.2969*(1-TK/273.16))) + \
                             0.42873e-3*(10**(4.76955*(1-TK/273.16))-1) + 0.78614)
                        self.met['hum'].append(100.0*wvp/ew)

            # if met data are not there, load gridded data
            else:
                noSiteData = True
                
        if noSiteData: # if met data are not there, load gridded data
            print 'Meteo data for '+self.name+\
                  ' not found. Using VMF1 grids instead.'

            lat = np.arange(90,-92,-2)
            lon = np.arange(0,359,2.5)
#                grid_lat, grid_lon = np.meshgrid(lat,lon)
            
            # for each day load data
            for day in dates:
                # mjd:
                yy = day.year
                mm = day.month
                dd = day.day
                if mm <= 2: #January & February
                    yy  = yy - 1.0
                    mm = mm + 12.0
                jd = floor( 365.25*(yy + 4716.0)) + \
                     floor( 30.6001*( mm + 1.0)) + 2.0 - \
                     floor( yy/100.0 ) + \
                     floor( floor( yy/100.0 )/4.0 ) + dd - 1524.5
                mjd = jd - 2400000.5
                
                for hh in (0.0, 6.0, 12.0, 18.0):
                    self.met['mjd'].append(mjd + hh/24.0)

                vmf_grid_file = 'VMFG_{:4d}{:02d}{:02d}.H'.\
                                 format(day.year, day.month, day.day)
                
                # load gridded data
                for hh in ('00','06','12','18'):
                    # vmf1 grid files:
                    with open(inp['meteo_cat']+'/'+vmf_grid_file+hh,'r') as f:
                        f_lines = f.readlines()
                        f_lines = [line for line in f_lines if line[0]!='!']
                    
                    ahz_grid = []
                    awz_grid = []
                    zhz_grid = []
                    zwz_grid = []
                    grid = []
                    
                    for f_line in f_lines:
                        tmp = [float(i) for i in f_line.split()]
                        grid.append([tmp[0], tmp[1]])
                        ahz_grid.append(tmp[2])
                        awz_grid.append(tmp[3])
                        zhz_grid.append(tmp[4])
                        zwz_grid.append(tmp[5])
                    
                    ahz_grid = np.array(ahz_grid)
                    awz_grid = np.array(awz_grid)
                    zhz_grid = np.array(zhz_grid)
                    zwz_grid = np.array(zwz_grid)
                    
                    f = sp.interpolate.interp2d(lon, lat, \
                                                ahz_grid.reshape((91,144)))
                    ahz = f(self.lon_gcen, self.lat_geod)[0]
                    f = sp.interpolate.interp2d(lon, lat, \
                                                awz_grid.reshape((91,144)))
                    awz = f(self.lon_gcen, self.lat_geod)[0]
                    f = sp.interpolate.interp2d(lon, lat, \
                                                zhz_grid.reshape((91,144)))
                    zhz = f(self.lon_gcen, self.lat_geod)[0]
                    f = sp.interpolate.interp2d(lon, lat, \
                                                zwz_grid.reshape((91,144)))
                    zwz = f(self.lon_gcen, self.lat_geod)[0]
                    self.met['ahz'].append(ahz)
                    self.met['awz'].append(awz)
                    self.met['zhz'].append(zhz)
                    self.met['zwz'].append(zwz)
            # convert lists of one-valued array into arrays
            self.met['ahz'] = np.array(self.met['ahz'])
            self.met['awz'] = np.array(self.met['awz'])
            self.met['zhz'] = np.array(self.met['zhz'])
            self.met['zwz'] = np.array(self.met['zwz'])
            
        ## lhg files with precomputed site-specific data (tropo gradients):
        for day in dates:
            year = day.year
            doy = day.timetuple().tm_yday # day of year
            # lhg site files:
            lhg_file = '{:4d}{:03d}.lhg_r'.format(year, doy)
            with open(inp['meteo_cat']+'/'+lhg_file,'r') as f:
                f_lines = f.readlines()
            
            # if met data are present in the lhg site file, load 'em
            if any(self.name in s for s in f_lines):
                for f_line in f_lines:
#                    if self.name in f_line.split(' '):
                    if self.name in f_line.split():
                        tmp = [float(i) for i in f_line[10:].split()]
                        self.met['gnh'].append(tmp[1])
                        self.met['geh'].append(tmp[2])
                        self.met['gnw'].append(tmp[3])
                        self.met['gew'].append(tmp[4])

        ## now make interpolants for a faster access:
        self.fMet['fAhz'] = sp.interpolate.interp1d(self.met['mjd'],\
                                            self.met['ahz'], kind='cubic')
        self.fMet['fAwz'] = sp.interpolate.interp1d(self.met['mjd'],\
                                            self.met['awz'], kind='cubic')
        self.fMet['fZhz'] = sp.interpolate.interp1d(self.met['mjd'],\
                                            self.met['zhz'], kind='cubic')
        self.fMet['fZwz'] = sp.interpolate.interp1d(self.met['mjd'],\
                                            self.met['zwz'], kind='cubic')
        if len(self.met['TC'])>0:
            self.fMet['fT'] = sp.interpolate.interp1d(self.met['mjd'],\
                                                      self.met['TC'])
            self.fMet['fP'] = sp.interpolate.interp1d(self.met['mjd'],\
                                                      self.met['pres'])
            self.fMet['fH'] = sp.interpolate.interp1d(self.met['mjd'],\
                                                      self.met['hum'])
        if len(self.met['gnh'])>0:
            self.fMet['fGnh'] = sp.interpolate.interp1d(self.met['mjd'],\
                                            self.met['gnh'], kind='cubic')
            self.fMet['fGeh'] = sp.interpolate.interp1d(self.met['mjd'],\
                                            self.met['geh'], kind='cubic')
            self.fMet['fGnw'] = sp.interpolate.interp1d(self.met['mjd'],\
                                            self.met['gnw'], kind='cubic')
            self.fMet['fGew'] = sp.interpolate.interp1d(self.met['mjd'],\
                                            self.met['gew'], kind='cubic')
#        print self.met
        
    def j2000gp(self, r2000, gcrs=None, utc=None, t=None):
        '''
        Compute station GCRS state.
        Add up displacements due to geophysical effects to it.
        '''
        if self.name=='GEOCENTR':
            # keep the zeros
            return
        elif self.name=='RA':
            # interpolate RA's GCRS orbit to date
#            for jj in range(0,3):
#                self.r_GCRS[jj] = sp.interpolate.splev(t,fGcrs[jj],der=0)
#            for jj in range(3,6):
#                self.v_GCRS[jj-3] = sp.interpolate.splev(t,fGcrs[jj],der=0)
#            for jj in range(6,9):
#                self.a_GCRS[jj-6] = sp.interpolate.splev(t,fGcrs[jj],der=0)
            lag_order = 9
            x, _ = lagint(lag_order, utc, gcrs[:,6], t)
            y, _ = lagint(lag_order, utc, gcrs[:,7], t)
            z, _ = lagint(lag_order, utc, gcrs[:,8], t)
            vx, _ = lagint(lag_order, utc, gcrs[:,9], t)
            vy, _ = lagint(lag_order, utc, gcrs[:,10], t)
            vz, _ = lagint(lag_order, utc, gcrs[:,11], t)
            self.r_GCRS = np.hstack((x,y,z))
            self.v_GCRS = np.hstack((vx,vy,vz))
            try:
                ax, _ = lagint(lag_order, utc, gcrs[:,12], t)
                ay, _ = lagint(lag_order, utc, gcrs[:,13], t)
                az, _ = lagint(lag_order, utc, gcrs[:,14], t)
                self.a_GCRS = np.hstack((ax,ay,az))
            except:
                self.a_GCRS = np.zeros(3)
        else:
            # transformation GTRS -> GCRS
            self.r_GCRS = np.dot(r2000[:,:,0], self.r_GTRS_date)
            self.v_GCRS = np.dot(r2000[:,:,1], self.r_GTRS_date)
            self.a_GCRS = np.dot(r2000[:,:,2], self.r_GTRS_date)
            # geophysics:
#            print self.dr_tide
#            print self.dr_oclo
#            print self.dr_poltide
#            print self.dv_tide
#            print self.dv_oclo
#            print self.dv_poltide
#            raw_input()
#            print self.lat_geod, self.lon_gcen, self.h_geod
            self.r_GCRS += self.dr_tide + self.dr_oclo +\
                           self.dr_poltide + self.dr_atlo
            self.v_GCRS += self.dv_tide + self.dv_oclo +\
                           self.dv_poltide + self.dv_atlo

#    def thermal_def(self, alpha, delta, elv, T, C):
#        '''
#        thermal_def computes delta_tau due to the 
#        thermal deformation effect of the telescope
#    
#        input:
#           self - site object
#           alpha, delta, elv, T - right ascention, declination, elevation, 
#                                   air temperature in C
#           const
#        output delay:
#            dt_thermal
#        '''
#        
#        dl = -0.12*np.pi/180.0
#        phi0 = 39.06*np.pi/180.0
#        
#        # Antenna focus factor
#        if self.focus_type == 'FO_PRIM':
#            Fa = 0.9
#        else:
#            Fa = 1.8
#        
#        # Alt-azimuth
#        if self.mount_type == 'MO_AZEL':
#            dt_thermal = ( self.gamma_hf * (T - self.T0) * (self.hf * sin(elv)) + \
#                           self.gamma_hp * (T - self.T0) * (self.hp * sin(elv) + \
#                           self.AO * cos(elv) + self.hv - Fa * self.hs) ) / C
#        # Equatorial
#        elif self.mount_type == 'MO_EQUA':
#            dt_thermal = ( self.gamma_hf * (T - self.T0) * (self.hf * sin(elv)) + \
#                           self.gamma_hp * (T - self.T0) * (self.hp * sin(elv) + \
#                           self.AO * cos(delta) + self.hv - Fa * self.hs) ) / C
#        # XY north
#        elif self.mount_type == 'MO_XYNO':
#            dt_thermal = ( self.gamma_hf * (T - self.T0) * (self.hf * sin(elv)) + \
#                           self.gamma_hp * (T - self.T0) * (self.hp * sin(elv) + \
#    		 self.AO * sqrt( 1.0 - cos(elv)*cos(elv)*cos(alpha)*cos(alpha) ) + \
#    		 self.hv - Fa * self.hs) ) / C
#        # XY east
#        elif self.mount_type == 'MO_XYEA':
#            dt_thermal = ( self.gamma_hf * (T - self.T0) * (self.hf * sin(elv)) + \
#                           self.gamma_hp * (T - self.T0) * (self.hp * sin(elv) + \
#    		 self.AO * sqrt( 1.0 - cos(elv)*cos(elv)*cos(alpha)*cos(alpha) ) + \
#    		 self.hv - Fa * self.hs) ) / C
#        # misplaced equatorial RICHMOND
#        elif self.mount_type == 'MO_RICH':
#            dt_thermal = ( self.gamma_hf * (T - self.T0) * (self.hf * sin(elv)) + \
#                           self.gamma_hp * (T - self.T0) * (self.hp * sin(elv) + \
#                           self.AO * sqrt( 1.0 - ( sin(elv)*sin(phi0) + \
#    		 cos(elv)*cos(phi0)*(cos(alpha)*cos(dl) + sin(alpha)*sin(dl)) )**2 ) + \
#    		 self.hv - Fa * self.hs) ) / C    
#    
#        self.dtau_therm = dt_thermal

'''
#==============================================================================
#
#==============================================================================
'''
class source(object):
    '''
    Class containing a far-field source-specific data:
        - source name
        - source IVS-name
        - ra/dec
        - source type (C - calibrator, far-field
                       R - RadioAstron, near-field
                       G - GNSS, near-field
                       S - ESA's deep space s/c, near-field)
    '''
    def __init__(self, name, sou_type):
        self.name = name
        self.ivsname = name
        self.sou_type = sou_type
        self.radec = [] # in h m s
        self.ra = 0.0 # in rad
        self.dec = 0.0 # in rad
        self.K_s = np.zeros(3) # J2000.0 source unit vector
        
'''
#==============================================================================
# 
#==============================================================================
'''
class ephem(object):
    """
    Class containing spacecraft ephemerides
    in GTRS, GCRS and BCRS
    """
    def __init__(self, sc_name):
        self.sc_name = sc_name
        self.gtrs = np.array([]) # empty numpy array
        self.gcrs = np.array([]) # empty numpy array
#        self.bcrs = np.array([]) # empty numpy array
        self.bcrs = [] # empty list
        self.fGtrs = [] # spline interpolant
        self.fGcrs = [] # spline interpolant
        self.fBcrs = [] # spline interpolant, it will be a list of lists
        self.UT = np.array([]) # decimal time stamps for gtrs and gcrs
        self.CT = np.array([]) # decimal time stamps for bcrs
        self.CT_sec = np.array([]) # [s] time stamps for bcrs
#        self.lt_gc = np.array([]) # geocentric LTs
#        self.radec = np.array([]) # geocentric LT-corrected alpha/delta
#        self.fRadec = []
#        self.fLt_gc = 0.0
        
    def RaDec_bc_sec(self, jd, T_1, jpl_eph, return_R10=False):
        """
        Calculate geocentric [Ra, Dec] of the S/C
        
        input: 
            JD   - Julian Date
            t_1  - obs epoch in decimal days, TDB
            return_R10 - return vector R_GC - R_SC in BCRS or not
        output:
            ra, dec [rad] 
            lt_01 [seconds]
        """
        JD = deepcopy(jd)
        t_1 = deepcopy(T_1)
        
        const = constants()
        C = const.C  # m/s
        GM = const.GM

        # must go below 1 ps to stop iterating:
        precision = 1e-13
        # but should do so in no more than n_max iterations:
        n_max = 3
        # lagrange poly order for interpolation
        lag_order = 5
        
        # initial approximation:
        nn = 0
        lt_01_tmp = 0.0
        
        bcrs = self.bcrs[0]
        tdb = self.CT*86400.0
        
        # correct 'overnighter'
        if t_1 >= 86400.0:
            JD += 1
            t_1 -= 86400.0
        
        astropy_t_1 = Time(JD + t_1/86400.0, format='jd', scale='tdb', precision=9)
        eph_t_0 = datetime.datetime(*map(int, bcrs[0, :3]))
        # first time stamp negative?
        if tdb[0] < 0:
            eph_t_0 += datetime.timedelta(days=1)
        # cut eph and it's tomorrow?
        if tdb[0]//86400 > 0:
            eph_t_0 -= datetime.timedelta(days=tdb[0]//86400)
#        print astropy_t_1.datetime, eph_t_0
        dd = (astropy_t_1.datetime - eph_t_0).days
#        print 'dd = ', dd
        
#        print tdb, JD, t_1, dd, '\n'
#        print 't_1 = {:.18f}'.format(t_1)
        
        x, _ = lagint(lag_order, tdb, bcrs[:, 6], t_1+dd*86400)
        y, _ = lagint(lag_order, tdb, bcrs[:, 7], t_1+dd*86400)
        z, _ = lagint(lag_order, tdb, bcrs[:, 8], t_1+dd*86400)
        R_0 = np.hstack((x,y,z))
#        print 'R_0 = {:.18f} {:.18f} {:.18f}'.format(*R_0)
        
        ## Earth:
        rrd = pleph(astropy_t_1.jd, 3, 12, jpl_eph)
        earth = np.reshape(np.asarray(rrd), (3,2), 'F') * 1e3
        R_1 = earth[:,0]
#        print 'R_1 = {:.18f} {:.18f} {:.18f}'.format(*R_1)
        
        R_0_0 = R_0 - R_1
#        print 'R_0_0 =', R_0_0
        lt_01 = norm(R_0_0)/C
#        print 'lt_01 =', lt_01
        
        mjd = JD - 2400000.5
#        print mjd, t_1, lt_01
        astropy_t_0 = Time(mjd + (t_1 - lt_01)/86400.0, format='mjd', scale='tdb', precision=9)
#        print astropy_t_0.datetime
        ''' BCRS! state vectors of celestial bodies at JD+CT, [m, m/s]: '''
        state_ss = []
        for jj in (1,2,4,5,6,7,8,9,10,11):
            rrd = pleph(astropy_t_1.jd, jj, 12, jpl_eph)
            state_ss.append(np.reshape(np.asarray(rrd), (3,2), 'F') * 1e3)
        # back to UTC!:
#        t_0_0 = t_1 - lt_01
        
        while (abs(lt_01 - lt_01_tmp) > precision) and (nn < n_max):
            lt_01_tmp = lt_01
            t_0 = t_1 - lt_01
            astropy_t_0 = Time(mjd + t_0/86400.0, format='mjd', scale='tdb', precision=9)
#            print astropy_t_0.datetime
#            dd = (astropy_t_0.datetime - eph_t_0).days
#            print t_0, dd
            x, _ = lagint(lag_order, tdb, bcrs[:,6], t_0+dd*86400)
            y, _ = lagint(lag_order, tdb, bcrs[:,7], t_0+dd*86400)
            z, _ = lagint(lag_order, tdb, bcrs[:,8], t_0+dd*86400)
            vx, _ = lagint(lag_order, tdb, bcrs[:,9], t_0+dd*86400)
            vy, _ = lagint(lag_order, tdb, bcrs[:,10], t_0+dd*86400)
            vz, _ = lagint(lag_order, tdb, bcrs[:,11], t_0+dd*86400)
            R_0 = np.hstack((x,y,z))
            V_0 = np.hstack((vx,vy,vz))
            
            # vector needed for RLT calculation
            R_01 = R_1 - R_0
#            print '(t_0-t_0_0)=', (t_0-t_0_0)
            
            RLT = 0.0
            
            for j, ii in enumerate((1,2,4,5,6,7,8,9,10,11)):
                rrd = pleph(astropy_t_0.jd, ii, 12, jpl_eph)
                state = np.reshape(np.asarray(rrd), (3,2), 'F') * 1e3
                R_B = state[:,0]

                R_0_B  = R_B - R_0
                R_1_B  = state_ss[j][:,0] - R_1
                R_01_B = R_1_B - R_0_B
#                print ii, R_0, rb, vb, t_0, t_0_0
                RLT += (2.0*GM[ii-1]/C**3) * \
                      log(  ( norm(R_0_B) + norm(R_1_B) + norm(R_01_B) + \
                              2.0*GM[ii-1]/C**2 ) / \
                            ( norm(R_0_B) + norm(R_1_B) - norm(R_01_B) + \
                              2.0*GM[ii-1]/C**2 ) )
#            print 'RLT=', RLT
            lt_01 = lt_01 - (lt_01 - norm(R_01)/C - RLT) / \
                        ( 1.0 - np.dot(R_01, V_0)/(C*norm(R_01)) )
#            print 'lt_01=', lt_01
            t_0 = t_1 - lt_01

            nn += 1

        x, _ = lagint(lag_order, tdb, bcrs[:, 6], t_0+dd*86400)
        y, _ = lagint(lag_order, tdb, bcrs[:, 7], t_0+dd*86400)
        z, _ = lagint(lag_order, tdb, bcrs[:, 8], t_0+dd*86400)
        R_0 = np.hstack((x, y, z))
        r = R_0 - R_1
#        print 'R_0 = {:.18f} {:.18f} {:.18f}'.format(*R_0)
#        print 'R_0-R_1 = {:.18f} {:.18f} {:.18f}'.format(*r)
#        print 'lt_01 = {:.18f}'.format(lt_01)

#        raw_input()

        # S/C position is given at a moment LT seconds ago, which
        # means r is abberated in the far-field case sense

        ra = np.arctan2(r[1], r[0])  # right ascension
        dec = np.arctan(r[2]/np.sqrt(r[0]**2+r[1]**2))  # declination
        if ra < 0:
            ra += 2.0*np.pi

#        if dec > -1*np.pi/180: raise Exception()

        if not return_R10:
            return ra, dec, lt_01
        else:
            return ra, dec, lt_01, r
        
'''  
#==============================================================================
# 
#==============================================================================
'''
class ion(object):
    '''
    Class containing ionospheric TEC-maps
    '''    
    def __init__(self, date_start, date_stop, inp):
        self.lat = np.arange(87.5,-90.0,-2.5) # len = 71
        self.lon = np.arange(-180.0,185.0,5.0) # len = 73
        self.date_tec = []
        vTEC_grid = []
        
        # set up dates
        day_start = datetime.datetime(date_start.year,date_start.month,\
                                  date_start.day,0,0,0)
        day_stop = datetime.datetime(date_stop.year,date_stop.month,\
                          date_stop.day,0,0,0) + datetime.timedelta(days=1)

        dd = (day_stop - day_start).days

        # make list with datetime objects ranging from 1st day to last+1
        dates = [day_start]
        for d in range(1,dd+1):
            dates.append(day_start+ datetime.timedelta(days=d))
        
        for day in dates:
            yy = str(day.year)[2:]
            doy = day.timetuple().tm_yday # day of year
            ionex = inp['iono_model'] + 'g{:03d}0.'.format(doy) + yy + 'i'
            try:
                with open(inp['ion_cat']+'/'+ionex,'r') as f:
                    f_lines = f.readlines()
            except Exception, err:
                print str(err)
#                print 'Failed to load TEC data, no iono delay will be computed.'
                raise Exception('Failed to load TEC data, no iono delay will be computed.')
                    
            for jj in range(len(f_lines)):
                if 'START OF TEC MAP' in f_lines[jj]:
                    tmp = [int(i) for i in f_lines[jj+1][0:50].split()]
                    # end of day1==start of day2 => don't include duplicate:
                    if len(self.date_tec)==0 or \
                            datetime.datetime(*tmp)!=self.date_tec[-1]:
                        self.date_tec.append(datetime.datetime(*tmp))
                    else:
                        continue

                    grid = []
                    for ii in range(jj+2,(jj+2)+71*6,6):
#                        lat = float(f_lines[ii][0:8])
                        # you can't transpose a 1D-array, have to first 
                        # convert it to a matrix:
#                        tmp1 = np.array(np.matrix(lat*np.ones(73)).T)
#                        tmp2 = np.array(np.matrix(np.arange(-180.0,185.0,5.0)).T)
#                        tmp = np.hstack((tmp1,tmp2))
                        tecs = []
                        for kk in range(1,6):
                            tmpflt = [float(i) for i in f_lines[ii+kk].split()]
                            tecs.append(tmpflt)
                        tecs = [item for sublist in tecs for item in sublist]
#                        tecs = np.array(np.matrix(tecs).T)
                        tecs = np.array(tecs)
#                        tmp = np.hstack((tmp,tecs))
                        if len(grid)==0:
#                            grid = tmp
                            grid = tecs
                        else:
#                            grid = np.vstack((grid,tmp))
                            grid = np.hstack((grid,tecs))
                    # reshape for interpolation:                            
                    vTEC_grid.append(grid.reshape((71,73)))
                    jj = ii # shift current index
        vTEC_grid = np.array(vTEC_grid)
        # the resulting array has a shape (dd*12+1,71x73): dd*12+1 epochs 
        # each containing TEC values on a 2D-grid
#        import matplotlib.pyplot as plt
#        import seaborn as sns
#        sns.set_style('whitegrid') # plot em niice!
#        plt.close('all')
#        print dd
#        for i in range(vTEC_grid.shape[0]):
#            plt.figure()
#            plt.imshow(vTEC_grid[i,:,:])
#            plt.show()
        # build a decimal time scale for a more convenient interpolation
        # in the future
        self.UT_tec = []
        for t in self.date_tec:
            self.UT_tec.append( (t.hour+t.minute/60.0+t.second/3600.0)/24.0 + \
                                (t-day_start).days )
#        print self.UT_tec
        # produce interpolants:
        self.fVTEC = []
        for vTEC in vTEC_grid:
            f = sp.interpolate.interp2d(self.lon, self.lat, vTEC, kind='cubic')
            self.fVTEC.append(f)

#%%
'''
#==============================================================================
# Load binary delays
#==============================================================================
'''
class bindel(object):
    '''
    Parse binary delay files in SFXC format
    '''
    def __init__(self, fname, fdir='.'):
        self.fname = fname
        self.fdir = fdir
        
        self.parse()

    def parse(self):
        self.scans = {}
        with open(os.path.join(self.fdir, self.fname),'rb') as f_del:
#            content = f_del.read()
#            line = struct.unpack('<i2sx', content[:7])

            # binary header (header_size, sta_name):
            self.header_size = struct.unpack('<i', f_del.read(4))[0] 
            # (unpack always returns a tuple)
            self.sta = struct.unpack('<2sx', f_del.read(self.header_size))[0]
            #print self.sta
            
            # unpack scan by scan
            sn = 1 # scan number
            while True:
                try:
                    source = struct.unpack('<80sx', f_del.read(81))[0]
                    mjd = struct.unpack('<i', f_del.read(4))[0]
                    # init dict for scan
                    self.scans[sn] = {'source':None, 'mjd':None, 'time':None, \
                                      'uvw':None, 'delay':None, \
                                      'phase':None, 'amp':None}
                    self.scans[sn]['source'] = source
                    self.scans[sn]['mjd'] = mjd
                    
                    time, uvw, delay, phase, amp = [], [], [], [], []
                    
                    while True:                
                        t,u,v,w,d,p,a = struct.unpack('<7d', f_del.read(8*7))
                        if t==0 and d==0:
                            break
                        else:
                            time.append(t)
                            uvw.append((u,v,w))
                            delay.append(d)
                            phase.append(p)
                            amp.append(a)
                            
                    #print np.array(time)

                    self.scans[sn]['time'] = np.array(time)
                    self.scans[sn]['uvw'] = np.array(uvw)
                    self.scans[sn]['delay'] = np.array(delay)
                    self.scans[sn]['phase'] = np.array(phase)
                    self.scans[sn]['amp'] = np.array(amp)
                    
                    sn += 1

                except:
                    break
                
    def getSources(self):
        '''
        Return list of sources
        '''
        return list(set([self.scans[sn]['source'].strip() \
                         for sn in self.scans.keys()]))
                             
    def phaseCor(self, source, day, t, f):
        '''
        Integrate f_gc to get a phase correction for a S/C
        input:
            source - source name
            t - seconds from beginning of the first day of experiment
            f - freqs in Hz
        '''
        # set f_0 - will integrate difference (f-f_0) to reduce phase dyn.range
        f_0 = f[0]
        for sn in self.scans.keys():
            if self.scans[sn]['source'].strip() == source:
                at = Time(self.scans[sn]['mjd'], format='mjd', scale='utc')
                dt = at.datetime
                # same time scale as t:
                ts = 86400.0*(dt-day).days + self.scans[sn]['time']
                # cut proper piece of (f-f_0) to integrate:
                bounds = np.searchsorted(t, (ts[0]+1, ts[-1]))
#                print bounds
                ti = t[bounds[0]-1 : bounds[1]+1]
                fi = f[bounds[0]-1 : bounds[1]+1] - f_0
#                if sn==191:
#                    print ti
#                    print fi
#                    raw_input('lala')
#                print len(fi)
                try:
#                    # make an optimal polyfit:
#                    ffit = optimalFit(ti, fi, min_order=3, \
#                                      max_order=7, fit_type='poly')
#                    # integrate it to a phase poly
#                    pint = np.polyint(ffit.best_estimator_.coef_)
#                    # evaluate it at scan start
#                    i0 = np.polyval(pint, ts[0])
#                    # the phase is then int(ts[i])-int(ts[0])
#                    phase = np.array([2.0*np.pi*(np.polyval(pint, to) - i0) \
#                                      for to in ts])
#                    print phase
                    # scale ti for a more robust fit:
                    mu_ti = np.array([np.mean(ti), np.std(ti)])
                    ti = (ti-mu_ti[0])/mu_ti[1]
                    # make an optimal polyfit:
                    ffit = optimalFit(ti, fi, min_order=3, \
                                      max_order=7, fit_type='poly')
#                    print len(ffit.best_estimator_.coef_)
                    # integrate it to a phase poly
                    pint = np.polyint(ffit.best_estimator_.coef_)
                    # evaluate it at scan start
                    i0 = np.polyval(pint, (ts[0] - mu_ti[0]) / mu_ti[1])
                    # the phase is then int(ts[i])-int(ts[0])
                    phase = np.array([2.0*np.pi*(np.polyval(pint, to) - i0) \
                                      for to in (ts-mu_ti[0])/mu_ti[1]])*mu_ti[1]
#                    print phase
#                    raw_input('ololo?')
                    # store it:
                    self.scans[sn]['phase'] = phase
                except:
                    # too few points to make a fit - skip this scan then
                    continue
                
    def dump(self, binary=True, txt=False, out_name=None, out_dir=None):
        '''
        Dump parsed (and processed) data back to a binary .del-file
        '''
        if out_name is None:
            dot = self.fname.index('.')
            out_name = self.fname[:dot] + 'i.del'
#            out_name = self.fname
        if txt:
            dot = self.fname.index('.')
            out_name_txt = self.fname[:dot] + '.txt'
        if out_dir is None:
            out_dir = self.fdir
        # create output dir if it doesn't exist:
        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)
        
        if binary:
            with open(os.path.join(out_dir, out_name),'wb') as f_del:
                # binary header (header_size, sta_name):
                line = struct.pack('<i2sx', self.header_size, self.sta)
                f_del.write(line)
                
                for sn in self.scans.keys():
                    # source name and mjd
                    line = struct.pack('<80sxi', self.scans[sn]['source'],
                                                 self.scans[sn]['mjd'])
                    f_del.write(line)
                    for t, (u,v,w), d, p, a in zip(self.scans[sn]['time'],
                          self.scans[sn]['uvw'], self.scans[sn]['delay'],
                          -self.scans[sn]['phase'], self.scans[sn]['amp']):
    #                      np.zeros_like(self.scans[sn]['delay']), self.scans[sn]['amp']):
                        line = struct.pack('<7d', t,u,v,w,d,p,a)
                        f_del.write(line)
                    # trailing zeros at scan end:
                    line = struct.pack('<7d', *list(np.zeros(7)))
                    f_del.write(line)

        # dump txt:
        if txt:
            with open(os.path.join(out_dir, out_name_txt),'w') as f_txt:
                # binary header (header_size, sta_name):
                line = '{:s}\n'.format(self.sta)
                f_txt.write(line)
                
                for sn in self.scans.keys():
                    # source name and mjd
                    line = '{:s} {:f}\n'.format(self.scans[sn]['source'].strip(),
                                                 self.scans[sn]['mjd'])
                    f_txt.write(line)
                    for t, (u,v,w), d, p, a in zip(self.scans[sn]['time'],
                          self.scans[sn]['uvw'], self.scans[sn]['delay'],
                          -self.scans[sn]['phase'], self.scans[sn]['amp']):
                         # np.zeros_like(self.scans[sn]['delay']), self.scans[sn]['amp'])
                        line = '{:f} {:f} {:f} {:f} {:.15e} {:.15e} {:.15e}\n'.format(t,
                                                                   u,v,w,d,p,a)
                        f_txt.write(line)

