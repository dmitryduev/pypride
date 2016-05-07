# -*- coding: utf-8 -*-
"""
A library containing a handful of useful functions.

Created on Wed Oct  2 16:23:02 2013

@author: Dr. Dmitry A. Duev
"""

from math import *
import numpy as np
import scipy as sp
#from scipy import interpolate
from scipy.interpolate import BarycentricInterpolator as bi
# import matplotlib.pyplot as plt

import datetime
from astropy.time import Time, TimeDelta
import astropy.units as units
from astropy.coordinates import EarthLocation, Angle
import os
import re # regular expressions
import urllib2
import requests
from ftplib import FTP
import paramiko # ssh client
import gzip
import collections
from copy import deepcopy

# import all class declarations
from pypride.classes import *

## fortran stuff
#from lagint import lagint
#from pleph import pleph
#from iau_xys00a import iau_xys00a_fort
#from admint2 import admint2
try:
    from pypride.vintflib import lagint, pleph, iau_xys00a_fort, admint2
except:
    # compile the Fortran code if necessary
    import numpy.f2py as f2py
    fid = open('vintflib.f')
    source = fid.read()
    fid.close()
    f2py.compile(source, modulename='vintflib')
    from pypride.vintflib import lagint, pleph, iau_xys00a_fort, admint2
#from pypride.vintflib import lagint, pleph, iau_xys00a_fort, admint2#, lagintt

## parallelism
import multiprocessing as mp

#from decimal import *

## mostly used numpy stuff:
from numpy.linalg import norm
from numpy import dot
np.set_printoptions(precision=18)
from numpy.polynomial import chebyshev as cheb
#cheb = np.polynomial.chebyshev

#from numba import jit, double#, autojit
import numba

import inspect

#from time import time as _time

# go off and set up missing folders if git-cloned:
#if 'vispy' in os.getcwd():
#    if not os.path.isdir('ion'):
#        os.makedirs('ion')
#    if not os.path.isdir('meteo'):
#        os.makedirs('meteo')
#    if not os.path.isdir('sc_eph'):
#        os.makedirs('sc_eph')
#    for sc_name in ('mex', 'vex', 'radioastron', 'gnss', 'gaia'):
#        pth = os.path.join('sc_eph','raw_'+sc_name)
#        if not os.path.isdir(pth):
#            os.makedirs(pth)
#abs_path = os.path.dirname(os.path.abspath(__file__))
abs_path = os.path.dirname(inspect.getfile(inspect.currentframe()))
if not os.path.isdir(os.path.join(abs_path, 'ion')):
    os.makedirs(os.path.join(abs_path, 'ion'))
if not os.path.isdir(os.path.join(abs_path, 'meteo')):
    os.makedirs(os.path.join(abs_path, 'meteo'))
if not os.path.isdir(os.path.join(abs_path, 'sc_eph')):
    os.makedirs(os.path.join(abs_path, 'sc_eph'))
for sc_name in ('mex', 'vex', 'rosetta', 'radioastron', 'gnss', 'gaia'):
    pth = os.path.join(abs_path, 'sc_eph','raw_'+sc_name)
    if not os.path.isdir(pth):
        os.makedirs(pth)

'''
#==============================================================================
# Flatten (irregular) list of lists
#==============================================================================
'''
def flatten_generator(l):
    for el in l:
        if isinstance(el, collections.Iterable) and not isinstance(el, basestring):
            for sub in flatten(el):
                yield sub
        else:
            yield el
            
def flatten(x):
    result = []
    for el in x:
        if hasattr(el, "__iter__") and not isinstance(el, basestring):
            result.extend(flatten(el))
        else:
            result.append(el)
    return result

#flatten = lambda x: [y for l in x for y in flatten(l)] if type(x) is list else [x]

'''
#==============================================================================
# Factorise a number
#==============================================================================
'''
def factors(n):    
    facs = set(reduce(list.__add__, 
                ([i, n//i] for i in range(1, int(n**0.5) + 1) if n % i == 0)))
    return sorted(list(facs))

'''
#==============================================================================
# Use interpolator as exptrapolator
#==============================================================================
'''
def extrap1d(interpolator):
    xs = interpolator.x
    ys = interpolator.y

    def pointwise(x):
        if x < xs[0]:
            return ys[0]+(x-xs[0])*(ys[1]-ys[0])/(xs[1]-xs[0])
        elif x > xs[-1]:
            return ys[-1]+(x-xs[-1])*(ys[-1]-ys[-2])/(xs[-1]-xs[-2])
        else:
            return interpolator(x)

    def ufunclike(xs):
        return np.array(map(pointwise, np.array(xs)))

    return ufunclike

'''
#==============================================================================
# Linear extrapolator
#==============================================================================
'''
def extrap(x, xp, yp, interp_type='linear'):
    """np.interp function with linear extrapolation"""
    # convert to numpy arrays if necessary:
    if type(x).__module__ != np.__name__:
        x = np.array(x)
    if type(xp).__module__ != np.__name__:
        xp = np.array(xp)
    if type(yp).__module__ != np.__name__:
        yp = np.array(yp)
    if interp_type=='linear':
        y = np.interp(x, xp, yp)
    elif interp_type=='lag':
        y, _ = lagint(min(len(xp),15), xp, yp, x)
    elif interp_type=='cheb':
        p = cheb.chebfit(xp, yp, min(len(xp), 15))
        y = cheb.chebval(x, p)
    # linear fit outside:
    y[x < xp[0]] = yp[0] + (x[x<xp[0]]-xp[0]) * (yp[0]-yp[1]) / (xp[0]-xp[1])
    y[x > xp[-1]]= yp[-1] + (x[x>xp[-1]]-xp[-1])*(yp[-1]-yp[-2])/(xp[-1]-xp[-2])
    return y

'''
#==============================================================================
# Memoize function values on previously used arguments
#==============================================================================
'''
def memoize(f):
    '''
    Memoize function values on previously used arguments
    '''
    def memf(*x):
        if x not in memf.cache:
            memf.cache[x] = f(*x)
        return memf.cache[x]
    
    memf.cache = {}
    return memf

'''
#==============================================================================
# Check internet connection
#==============================================================================
'''
def internet_on():
    ''' Check internet connection '''
    try:
        urllib2.urlopen('http://google.com', timeout=1)
        return True
    except urllib2.URLError:
        pass
    return False

'''
#==============================================================================
# Clever Lagrange interpolation
#==============================================================================
'''
@memoize
def make_lag_bi(x, y):
    '''
    construct Lagrange interpolator
    input should be two tupils (immutable, thus hushable)
    '''
    L = bi(x)
    L.set_yi(y)
    return L

'''
#==============================================================================
# Optimal fit to data
#==============================================================================
'''
def optimalFit(x, y, min_order=0, max_order=8, fit_type='poly'):
    # initialise optimal estimator:
    if fit_type=='poly':
        estimator = PolynomialRegression()
    elif fit_type=='cheb':
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
# Smooth the data using a window with requested size
#==============================================================================
'''
def smooth(x, window_len=11, window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TO DO: the window parameter could be the window itself if an array 
          instead of a string
    NOTE: length(output) != length(input), to correct this: 
       return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."

    if window_len<3:
        return x

    if window not in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, \
        "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"


    s=np.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')

#    return y
    return y[(window_len/2-1):-(window_len/2)]


'''
#==============================================================================
# Convert time stamp in seconds to HMS
#==============================================================================
'''
def sec2hms(sec):
    sec = int(sec)
    days = sec/86400
    hours = sec/3600 - days * 24
    minutes = sec/60 - days * 1440 - hours * 60
    seconds = fmod(sec, 60)
    return [hours, minutes, seconds]

'''
#==============================================================================
# Convert from spectral to lag domain
#==============================================================================
'''
def Sp2Lag(spec):
    spf = np.zeros_like(spec)
    for n, v in enumerate(spec):
        spf[n] = v*(-1)**n
#        spf[jj] = v
    pad = np.zeros_like(spec)
    spfp = np.hstack((spf, pad))
    return np.fft.fft(spfp)

'''
#==============================================================================
# Get Fringe SNR
#==============================================================================
'''
def weighted_avg_and_var(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, weights=weights)   
    variance = np.average((values-average)**2, weights=weights)  # Fast and numerically precise
#    std = sqrt(variance)
    return (average, variance)

def GetSNR(data, lmax, spread, void):

    dx = data[lmax-spread:lmax+spread+1]
    kx = np.arange(lmax-spread,lmax+spread+1)
    
    weights = (np.abs(kx - lmax)) > void
    
#    fit = cheb.chebfit(kx, dx, 6, w=weights)
#    print fit
#    dx -= cheb.chebval(kx, fit)
    
#    np.ma.masked_inside() - instead of using weights? hm..
    
#    fit = np.polyfit(kx, dx, 6, w=weights)
    fitC = np.polynomial.Chebyshev.fit(kx, dx, 6, w=weights)
#    figa = plt.figure()
#    ax = figa.add_subplot(211)
#    ax.plot(kx, dx)
#    ax.plot(kx, np.polyval(fit, kx))
#    ax.plot(kx, fitC(kx))
#    print np.polyval(fit, kx)
#    dx -= np.polyval(fit, kx)
    dx -= fitC(kx)
#    ax2 = figa.add_subplot(212)
#    ax2.plot(kx, dx)
    
#    print weighted_avg_and_var(dx, weights)
    rd = np.sqrt(weighted_avg_and_var(dx, weights)[1])
    pa = max(dx)

#    print pa, rd
    return pa/rd

'''
#==============================================================================
# Find function maximum using a cubic fi
#==============================================================================
'''
def FindMax(Spec, Fmin, Fmax):

    Fmin = int(Fmin)
    Fmax = int(Fmax)
    # create index mask: 1 means mask and don't use, 0 - don't mask and use
    mask = np.hstack((np.ones(Fmin), \
                      np.zeros(Fmax-Fmin), np.ones(len(Spec)-Fmax)))
    mask = np.array(mask, dtype=bool)

    sp_masked = np.ma.array(Spec, mask=mask)
    mx = np.max(sp_masked)
    jmax = np.argmax(sp_masked)
    
    a2 = 0.5*(Spec[jmax-1] + Spec[jmax+1] - 2.0*Spec[jmax])
    a1 = 0.5*(Spec[jmax+1] - Spec[jmax-1])
    djx = -a1/(2.0*a2)
    xmax = jmax + djx
    
    return jmax, xmax, mx

'''
#==============================================================================
# Factorise a number
#==============================================================================
'''
def factorise(n):
    # factorise a number
    facs = set(reduce(list.__add__, 
                ([i, n//i] for i in range(1, int(n**0.5) + 1) if n % i == 0)))
    return sorted(list(facs))

'''
#==============================================================================
# Calculate derivative using polyfit
#==============================================================================
'''
def derivative(x, y, points=5, poly=2):
    '''
    Calculate a derivative dy/dx at each x
    Don't forget to normalise output
    '''
    dydx = np.zeros_like(y)
    for n0, xi in enumerate(x):
        # number of points to cut from the left-hand side
        nl = int(np.floor(points/2.0))
        # number of points to cut from the right-hand side
        nr = int(np.ceil(points/2.0))
        # check/correct bounds:
        if len(x[:n0]) < nl:
            nr = nr + nl - len(x[:n0])
            nl = len(x[:n0])
        if len(x[n0:]) < nr:
            nl = nl + nr - len(x[n0:])
            nr = len(x[n0:])
        
        # make a fit
        yfit = np.polyfit(x[n0-nl:n0+nr], y[n0-nl:n0+nr], poly)
        dydx[n0] = np.polyval(np.polyder(yfit), xi)
    
    return dydx
    
'''
#==============================================================================
# Interpolate using polyfit
#==============================================================================
'''
def interpolate(x, y, xn, points=5, poly=2):
    '''
    interpolate at each xn
    Don't forget to normalise output
    '''
    yn = np.zeros_like(xn)
    for nn, xi in enumerate(xn):
        n0 = np.searchsorted(x, xi)
        # number of points to cut from the left-hand side
        nl = int(np.floor(points/2.0))
        # number of points to cut from the right-hand side
        nr = int(np.ceil(points/2.0))
        # check/correct bounds:
        if len(x[:n0]) < nl:
            nr = nr + nl - len(x[:n0])
            nl = len(x[:n0])
        if len(x[n0:]) < nr:
            nl = nl + nr - len(x[n0:])
            nr = len(x[n0:])
        
        # make a fit
        yfit = np.polyfit(x[n0-nl:n0+nr], y[n0-nl:n0+nr], poly)
        yn[nn] = np.polyval(yfit, xi)
    
    return yn
    
'''
#==============================================================================
# leap seconds
#==============================================================================
'''
#@jit
def nsec(mjd):
    ''' 
    This routine determines the number of leap seconds or
    difference TAI - UTC_iers starting from 1972 January 1 
    '''
    idelt = 0
    if (mjd >= 41317.0 and mjd < 41499.0): idelt = 10 # 1972 JAN 1
    elif (mjd >= 41499.0 and mjd < 41683.0): idelt = 11 # 1972 JUL 1
    elif (mjd >= 41683.0 and mjd < 42048.0): idelt = 12 # 1973 JAN 1
    elif (mjd >= 42048.0 and mjd < 42413.0): idelt = 13 # 1974 JAN 1
    elif (mjd >= 42413.0 and mjd < 42778.0): idelt = 14 # 1975 JAN 1
    elif (mjd >= 42778.0 and mjd < 43144.0): idelt = 15 # 1976 JAN 1
    elif (mjd >= 43144.0 and mjd < 43509.0): idelt = 16 # 1977 JAN 1
    elif (mjd >= 43509.0 and mjd < 43874.0): idelt = 17 # 1978 JAN 1
    elif (mjd >= 43874.0 and mjd < 44239.0): idelt = 18 # 1979 JAN 1     
    elif (mjd >= 44239.0 and mjd < 44786.0): idelt = 19 # 1980 JAN 1
    elif (mjd >= 44786.0 and mjd < 45151.0): idelt = 20 # 1981 JUL 1
    elif (mjd >= 45151.0 and mjd < 45516.0): idelt = 21 # 1982 JUL 1
    elif (mjd >= 45516.0 and mjd < 46247.0): idelt = 22 # 1983 JUL 1
    elif (mjd >= 46247.0 and mjd < 47161.0): idelt = 23 # 1985 JUL 1
    elif (mjd >= 47161.0 and mjd < 47892.0): idelt = 24 # 1988 JAN 1
    elif (mjd >= 47892.0 and mjd < 48257.0): idelt = 25 # 1990 JAN 1
    elif (mjd >= 48257.0 and mjd < 48804.0): idelt = 26 # 1991 JAN 1
    elif (mjd >= 48804.0 and mjd < 49169.0): idelt = 27 # 1992 JUL 1
    elif (mjd >= 49169.0 and mjd < 49534.0): idelt = 28 # 1993 JUL 1
    elif (mjd >= 49534.0 and mjd < 50083.0): idelt = 29 # 1994 JUL 1
    elif (mjd >= 50083.0 and mjd < 50630.0): idelt = 30 # 1996 JAN 1
    elif (mjd >= 50630.0 and mjd < 51179.0): idelt = 31 # 1997 JUL 1
    elif (mjd >= 51179.0 and mjd < 53736.0): idelt = 32 # 1999 JAN 1
    elif (mjd >= 53736.0 and mjd < 54832.0): idelt = 33 # 2006 JAN 1
    elif (mjd >= 54832.0 and mjd < 56109.0): idelt = 34 # 2009 JAN 1
    elif (mjd >= 56109.0 and mjd < 57204.0): idelt = 35 # 2012 JUL 1
    elif (mjd >= 57204.0):                   idelt = 36 # 2015 JUL 1
    
    return float(idelt)

'''
#==============================================================================
#     Calculate Modified Julian Date
#==============================================================================
'''
def mjuliandate(year,month,day,hour=0.0,minu=0.0,sec=0.0):
    '''
    Calculate Modified Julian Date
    '''
    if month <= 2: #January & February
        year  = year - 1.0
        month = month + 12.0
    jd = floor( 365.25*(year + 4716.0)) + floor( 30.6001*( month + 1.0)) + 2.0 - \
         floor( year/100.0 ) + floor( floor( year/100.0 )/4.0 ) + day - 1524.5 + \
         (hour + minu/60.0 + sec/3600.0)/24.0;
    mjd = jd - 2400000.5;
    return mjd


'''
#==============================================================================
# Load frequency ramping parameters for a S/C for 2(3)-way Doppler
#==============================================================================
'''
def freqRamp(cat_dir=None, sc=None, tx_type=None):
    '''
        tx_type: ['1way', '3way']
    '''
    # read frequency ramping parameters
    if sc==None or tx_type==None:
        raise Exception('Can\'t load ramp params: S/C/Tx_type not specified.')
    else:
        if cat_dir==None:
            cat_dir = 'cats/ramp.' if tx_type=='3way' else 'cats/ramp1w.'
        rampFile = ''.join((cat_dir, sc.lower()))
        try:
            if tx_type=='3way':
                with open(rampFile, 'r') as f:
                    f_lines = f.readlines()
                    # skip comments
                    f_lines = [l for l in f_lines if l[0]!='#' and len(l)>10]
                    ramp = []
                    for line in f_lines:
                        line = line.split()
                        t_start = ''.join((line[0], ' ', line[1]))
                        t_stop = ''.join((line[2], ' ', line[3]))
                        # [t_start t_stop f_0 df uplink_sta]
                        ramp.append([datetime.datetime.strptime(t_start, \
                                                      "%Y-%m-%d %H:%M:%S.%f"), \
                                    datetime.datetime.strptime(t_stop, \
                                                      "%Y-%m-%d %H:%M:%S.%f"), \
                                    float(line[4]), float(line[5]), line[6] ])
                return ramp
                
            elif tx_type=='1way':
                # for ESA s/c time is in TDB!!
                with open(rampFile, 'r') as f:
                    f_lines = f.readlines()
                    # skip comments
                    f_lines = [l for l in f_lines if l[0]!='#' and len(l)>10]
                    ramp = []
                    for line in f_lines:
                        line = line.split()
                        t_start = ''.join((line[0], ' ', line[1]))
                        t_stop = ''.join((line[2], ' ', line[3]))
                        # [t_start t_stop f_0 df uplink_sta]
                        ramp.append([datetime.datetime.strptime(t_start, \
                                                      "%Y-%m-%d %H:%M:%S.%f"), \
                                    datetime.datetime.strptime(t_stop, \
                                                      "%Y-%m-%d %H:%M:%S.%f"), \
                                    float(line[4]), float(line[5]) ])
                return ramp
        except Exception, err:
            print str(err)
            raise Exception('Could not load ramp params for ' + sc + '.')
            
'''
#==============================================================================
# Load frequency value for a S/C for 1-way Doppler
#==============================================================================
'''
def freqConst(cat_file='cats/sc.freq', sc=None):
    # read frequency ramping parameters
    if sc==None:
        raise Exception('Can\'t load freq value: S/C not specified.')
    else:
        try:
            with open(cat_file, 'r') as f:
                f_lines = f.readlines()
                
            return [(float(x.split()[1]), x.split()[2]) for x in f_lines \
                    if x[0]!='#' and x.split()[0].lower()==sc.lower()][0]
                        
        except Exception, err:
            print str(err)
            raise Exception('Could not load frequency settings for '+sc+'.')


'''
#==============================================================================
# Download raw ephemeride files for ESA spacecraft
#==============================================================================
'''
def getESAraweph(source, orb_path='.', replaceDE=True):
    ''' 
    Download raw ephemeride files for VEX, MEX, and ROSETTA
    '''
    if internet_on():
        eph_url = 'http://tasc.esa.int/data/' + source[0:3].upper()
        response = urllib2.urlopen(eph_url)
        html = response.read()
        html = html.split('\n')
        p = re.compile('(?<=HREF=")(V|M|R)ORB_\d.*gz(?=">)', flags=re.IGNORECASE)
        orbFiles = []
        # get list of orb files
        for line in html:
            match = p.search(line)
            if match!=None: # if there's a match
                orbFiles.append(match.group(0))
        
        # check contents of orb_path folder (final raw ephs:)
        orbFilesLocal = [ f for f in os.listdir(orb_path) \
                            if os.path.isfile(os.path.join(orb_path, f)) and \
                            ('.gz' not in f) and ('ORB_0' in f)]
        # find updates
        new = [f for f in orbFiles if f[:-3] not in orbFilesLocal]
        
        # check contents of orb_path folder (planned raw ephs:)
        if source[0:3].upper()=='MEX' and 1==0:
            p = re.compile('(?<=HREF=")MORB_pln\d.*gz(?=">)', flags=re.IGNORECASE)
            plnOrbFiles = []
            # get list of orb files
            for line in html:
                match = p.search(line)
                if match!=None: # if there's a match
                    plnOrbFiles.append(match.group(0))
            if len(plnOrbFiles)>0:
                # get only the latest:
                plnOrbFile = sorted(plnOrbFiles)[-1]
            else:
                plnOrbFile = None
#            print plnOrbFile
            
            plnOrbFilesLocal = [ f for f in os.listdir(orb_path) \
                            if os.path.isfile(os.path.join(orb_path, f)) and \
                            ('.gz' not in f) and ('ORB_pln' in f)]
            if plnOrbFile is not None and plnOrbFile[:-3] not in plnOrbFilesLocal:
                # remove old planned data:
                for fOld in plnOrbFilesLocal:
                    os.remove(os.path.join(orb_path, fOld))
                # fetch the new file:
                new.append(plnOrbFile)
        
        # get 'em. if they exist
        if len(new)>0:
            for fNew in new:
                try:
                    url = eph_url+'/'+fNew
                    fu = urllib2.urlopen(url)
                    print "downloading " + fNew
                    # Save zipped file
                    with open(os.path.join(orb_path,fNew), "wb") as local_file:
                        local_file.write(fu.read())
                    # unzip:
                    print "unzipping " + fNew
                    with gzip.open(os.path.join(orb_path,fNew), 'rb') as fIn, \
                         open(os.path.join(orb_path,fNew[:-3]), "w") as fOut:
                        fOut.write(fIn.read())
                    # replace all D+/- with E+/-
                    if replaceDE:
                        with open(os.path.join(orb_path,fNew[:-3]), 'r') as fOut:
                            data = fOut.read()
                        data = data.replace('D+','E+')
                        data = data.replace('D-','E-')
                        with open(os.path.join(orb_path,fNew[:-3]), 'w') as fOut:
                            fOut.write(data)
                    # delete zipped file:
                    os.remove(os.path.join(orb_path,fNew))
                    # find and delete older local files:
                    # (eph numbers are the same, dates are different)
                    [os.remove(os.path.join(orb_path, x)) \
                        for x in orbFilesLocal if x[:9] in fNew]
                #handle errors
                except urllib2.HTTPError, e:
                    print "HTTP Error:", e.code, url
                except urllib2.URLError, e:
                    print "URL Error:", e.reason, url
            
            return True
    else:
        print 'No Internet connection, unable to get ESA\'s ephs.'
        return False

'''
#==============================================================================
# Check time boundaries of orb files
#==============================================================================
'''
def checkbound(source, orb_path='.', n=7):
    '''
    Check time boundaries of orb files
    '''
    # There must be a file 'boundaries.raw'
    # If it's >1 week old, recalculate it
    # first check that it actually exist:
    boundFile = 'boundaries.raw'
    path = os.path.join(orb_path, boundFile)
#    if os.path.isfile(path):
#        age = datetime.datetime.now() - \
#              datetime.datetime.utcfromtimestamp(os.path.getmtime(path))
    
    updates = getESAraweph(source, orb_path=orb_path, replaceDE=True)
    
    bounds = []
#    if (os.path.isfile(path) and age.days > n) or (not os.path.isfile(path)):
    if updates or (not os.path.isfile(path)):
        r = re.compile('[-T:]+')
        for root, subFolders, files in os.walk(orb_path):
            for fle in files:
                if 'ORB_' in fle and 'gz' not in fle and 'svn' not in fle:
                    with open(os.path.join(root, fle), 'r') as f:
                        data = f.readlines()
                    starts = [i for i,v in enumerate(data) if 'START_TIME' in v]
                    stops = [i for i,v in enumerate(data) if 'STOP_TIME' in v]
                    eq = data[starts[0]].index('=')            
                    t = r.split(data[starts[0]][eq+1:].strip())
                    t_start = map(int,t[:-1])
                    t_start.append(float(t[-1]))
                    eq = data[stops[-1]].index('=')            
                    t = r.split(data[stops[-1]][eq+1:].strip())
                    t_stop = map(int,t[:-1])
                    t_stop.append(float(t[-1]))
#                    print fle, t_start, t_stop
                    # return file_name, t_start, t_stop
                    bounds.append([fle, t_start, t_stop])

        with open(path, 'w') as f:
            for line in bounds:
#                print line[1]
                out = '{:s}  '.format(line[0])
                out += '{:5d} {:02d} {:02d} {:02d} {:02d} {:011.8f}   '.\
                        format(*line[1])
                out += '{:5d} {:02d} {:02d} {:02d} {:02d} {:011.8f}\n'.\
                        format(*line[2])
                f.write(out)
    else:
        # load file with bounds
        with open(path, 'r') as f:
            f_lines = f.readlines()
        for line in f_lines:
            lns = line.split()
            lns[1:] = map(float, lns[1:])
            bounds.append([lns[0], list(lns[1:7]), list(lns[7:])])
#        print bounds
    return bounds
 

'''
#==============================================================================
# Load time slots needed from ESA orb file
#==============================================================================
'''
def load_slots(orb_file, t_start, t_stop, t_step, source, inp=None):
    '''
    The orbital data proper are just lines providing at discrete time steps the
    epoch of the state, the state (position in km, velocity in km/s) and, if
    applicable, the state derivative (w.r.t time scale in days).
    
    cbody = 'venus' for VEX, 'mars' for MEX
    '''
    with open(orb_file,'r') as f:
        f_lines = f.readlines()
    
    slots = [ii for ii,v in enumerate(f_lines) if 'START_TIME' in v]
    
    slotData = []
    r = re.compile('[-T:]+')
    
    for slot in slots:
        out = [slot]
        for jj in (0,1):
            eq = f_lines[slot+jj].index('=')            
            t = r.split(f_lines[slot+jj][eq+1:].strip())
            ttag = map(int,t[:-1])
            ttag.append(float(t[-1]))
            out.append(ttag)
#            slot_start = datetime(*map(int,ttag))
        slotData.append(out)
    
    # line number of 'start_time', start_time, stop_time
    s = [x for x in slotData if (mjuliandate(*t_start)>=mjuliandate(*x[1]) and\
                                mjuliandate(*t_start)<=mjuliandate(*x[2])) or\
                                (mjuliandate(*t_start)<=mjuliandate(*x[1]) and\
                                (mjuliandate(*t_stop)<=mjuliandate(*x[2]) and\
                                 mjuliandate(*t_stop)>=mjuliandate(*x[1]))) or\
                                 (mjuliandate(*t_start)<=mjuliandate(*x[1]) and\
                                  mjuliandate(*t_stop)>=mjuliandate(*x[2]))]

    # load blocks of ephm corresponding to s'es
    eph_blocks = []
    # count not from block start (commented below), but t_start
    date_start = datetime.datetime(*t_start[0:3])
    for ii, sta, sto in s:
        kk = ii + 8
        block = []
        while len(f_lines[kk])>5 and ('META_START' not in f_lines[kk]):
            # skip state derivatives:
            if f_lines[kk][0]==' ': 
                kk += 1
                continue
            line_split = f_lines[kk].split()
            ttag = map(float, r.split(line_split[0]))
            dd = (datetime.datetime(*map(int,ttag[0:3])) - date_start).days
            ct = ttag[3] + ttag[4]/60.0 + ttag[5]/3600.0 + dd*24.0
            block.append(flatten([ct, map(float,line_split[1:])]))
            kk += 1
        eph_blocks.append( np.array(block) )


    ''' make ephemerides for output '''
    # fake obs
    ob = obs(['DUMMY'],'DUMMY','C')
    date_t_start = datetime.datetime(*t_start)
    date_t_stop = datetime.datetime(*t_stop)
    ob.addScan(date_t_start, t_step, stop=date_t_stop)
    
    ''' load input sittings: '''
    if inp==None:
#        inp = inp_set('inp.cfg')
        raise Exception('inp-file not provided')
    jpl_eph = inp['jpl_eph']
    cat_eop = inp['cat_eop']
#    jpl_eph = '/Users/oasis/_jive/python/vispy/jpl_eph/JPLEPH.421'
#    jpl_eph = '/Users/oasis/_jive/python/vispy/jpl_eph/JPLEPH.405'
#    cat_eop = '/Users/oasis/_jive/python/vispy/cats/eopc04.cat'
    
    #load eops
    ''' get the relevant eop entries from the catalogue: '''
    mjd_start = mjuliandate(date_t_start.year, date_t_start.month,\
                            date_t_start.day, date_t_start.hour, \
                            date_t_start.minute, date_t_start.second)
    with open(cat_eop, 'r') as fc:
        fc_lines = fc.readlines()
    eops = np.zeros((7,7)) # +/- 3 days
    for jj in range(len(fc_lines)):
        if fc_lines[jj][0]!=' ' and fc_lines[jj][0]!='*':
            entry = [float(x) for x in fc_lines[jj].split()]
            if len(entry) > 0 and entry[3] == np.floor(mjd_start) - 3:
                for kk in range(7):
                    eops[kk,0] = entry[3] # mjd
                    eops[kk,1] = entry[6] # UT1-UTC
                    eops[kk,2] = entry[6] - nsec(entry[3]) # UT1 - TAI
                    eops[kk,3] = entry[4] # x_p
                    eops[kk,4] = entry[5] # y_p
                    eops[kk,5] = entry[8] # dX
                    eops[kk,6] = entry[9] # dY
                    entry = [float(x) for x in fc_lines[jj+kk+1].split()]
                break #exit loop
    
    ''' stick relevant eph blocks together and cut relevant peace '''
    # raw cut
    mjd_end = mjuliandate(date_t_stop.year, date_t_stop.month,\
                            date_t_stop.day, date_t_stop.hour, \
                            date_t_stop.minute, date_t_stop.second)
    ii_start = [ii for ii,v in enumerate(s) if \
                        mjuliandate(*v[1]) <= mjd_start < mjuliandate(*v[2])][0]
    ii_end = [ii for ii,v in enumerate(s) if \
                        mjuliandate(*v[1]) <= mjd_end < mjuliandate(*v[2])][0]
    eph_blocks = eph_blocks[ii_start:ii_end+1]
    s = s[ii_start:ii_end+1]
    
    # stick
#    tmp = eph_blocks[0]
#    for b in eph_blocks[1:]:
#        # last entry of block n = first entry of block n+1. skip duplicates
#        tmp = np.vstack((tmp, b[1:]))
#    eph_blocks = tmp
    
    # tailored cut of each of the first and the last blocks
    # start and stop t in h to compare with eph_blocks
    th_start = date_t_start.hour + date_t_start.minute/60.0 +\
                date_t_start.second/3600.0
    th_end = date_t_stop.hour + date_t_stop.minute/60.0 +\
                date_t_stop.second/3600.0
    dd = (datetime.datetime(*t_stop[0:3]) -\
          datetime.datetime(*t_start[0:3])).days

    th_start_ii = np.searchsorted(eph_blocks[0][:,0], th_start)
    th_end_ii = np.searchsorted(eph_blocks[-1][:,0], th_end + dd*24.0)
    eph_blocks[0] = eph_blocks[0][th_start_ii-1:, :]
    eph_blocks[-1] = eph_blocks[-1][:th_end_ii+1, :]
    
#    print eph_blocks, len(eph_blocks)
#    print s
#    raw_input('bugagaga')
    '''
    make ephs only for the entries in eph_blocks, then interp 'em to 1s grid
    '''
    mjd0 = mjuliandate(date_t_start.year, date_t_start.month,date_t_start.day)
#    tic = _time()
    eph_bcrs_blocks = []
    eph_gcrs_blocks = []
    eph_gtrs_blocks = []

    for block in eph_blocks:
        eph_bcrs = []
        eph_gcrs = []
        eph_gtrs = []
        for state in block:
#            tic2 = _time()
            # full mjd:
            mjd = mjd0 + state[0]/24.0
            jd = mjd + 2400000.5
            
            rv = state[1:]
            
            # add barycentric r/v of cbody at CT:
            if source.lower()=='vex':
                cbodyState = pleph(jd, 2, 12, jpl_eph)
            elif source.lower() in ('her', 'rosetta'):
                cbodyState = pleph(jd, 4, 12, jpl_eph)
            elif source.lower()=='mex':
                cbodyState = pleph(jd, 4, 12, jpl_eph)
            else:
                raise Exception('Dunno the central body for ' + source)
            # append line to output SS barycentric ephemeris of VEX:
            rv_bcrs = rv + cbodyState
            eph_bcrs.append(rv_bcrs)
            
            # calculate geocentric orbits
            # GCRS first:
            earth = pleph(jd, 3, 12, jpl_eph)
            rv_gcrs = rv_bcrs - earth
            eph_gcrs.append(rv_gcrs)
            
            # then GTRS
            astrotime = Time(mjd, format='mjd', scale='tdb', precision=9)
            UTCdatetime = astrotime.utc.datetime
            UTC = (UTCdatetime.hour*3600.0 + \
                   UTCdatetime.minute*60.0 + \
                   UTCdatetime.second + UTCdatetime.microsecond*1e-6)/86400.0

            ''' interpolate eops to t_obs '''
#            UT1, eop_int = eop_iers(mjd, UTC, eops)
            UT1, eop_int = eop_iers(np.floor(mjd0+state[0]), UTC, eops)
            
            ''' rotation matrix IERS '''
            r2000 = ter2cel(UTCdatetime, eop_int, mode='der')
            rv_gtrs = np.zeros(6)
            rv_gtrs[0:3] = dot(r2000[:,:,0].T, rv_gcrs[0:3].T)
            rv_gtrs[3:6] = dot(r2000[:,:,0].T, rv_gcrs[3:6].T) + \
                           dot(r2000[:,:,1].T, rv_gcrs[0:3].T)
            eph_gtrs.append(rv_gtrs)
#            print _time() - tic2
        # turn into numpy array and append to block lists:
        eph_bcrs_blocks.append(np.array(eph_bcrs))
        eph_gcrs_blocks.append(np.array(eph_gcrs))
        eph_gtrs_blocks.append(np.array(eph_gtrs))

#    print _time()-tic
#    raw_input('bugagaga')
    
    # interpolate to 1 s grid
    t0 = datetime.datetime(*t_start[0:3])
    th = []
    for t in ob.tstamps:
        th.append(t.hour + t.minute/60.0 + t.second/3600.0 + 24.0*(t-t0).days)
    th = np.array(th)
    
#    print eph_bcrs_blocks
#    print len(eph_bcrs_blocks)
#    print th
    
#    print eph_blocks
#    raw_input('bugagaga')
#    tic = _time()
    bcrs = np.zeros((len(th), 6))
    gcrs = np.zeros((len(th), 6))
    gtrs = np.zeros((len(th), 6))
    # which eph_block to use with each time stamp
    ii_cur = []
#    print eph_blocks
    for thi in th:
        ii_cur.append([i for i,v in enumerate(eph_blocks) if \
                        v[0,0] <= thi < v[-1,0]][0])
    ii_cur = flatten(ii_cur)

    # choose proper eph block and cut array in a clever way:
    for ii in range(len(eph_blocks)):
        # clever: start and end indices for the current block
        i_s = ii_cur.index(ii)
        i_e = len(th) - ii_cur[::-1].index(ii)
#        print 'ii =',ii
#        print i_s, i_e
        for jj in range(6):
#            bcrs[i_s:i_e, jj], _ = lagint(9, eph_blocks[ii][:, 0], \
#                                       eph_bcrs_blocks[ii][:, jj], \
#                                       th[i_s:i_e])
#            gcrs[i_s:i_e, jj], _ = lagint(9, eph_blocks[ii][:, 0], \
#                                       eph_gcrs_blocks[ii][:, jj], \
#                                       th[i_s:i_e])
#            gtrs[i_s:i_e, jj], _ = lagint(9, eph_blocks[ii][:, 0], \
#                                       eph_gtrs_blocks[ii][:, jj], \
#                                       th[i_s:i_e])
            for kk in range(i_s, i_e):
                bcrs[kk, jj], _ = lagint(9, eph_blocks[ii][:, 0], \
                                       eph_bcrs_blocks[ii][:, jj], th[kk])
                gcrs[kk, jj], _ = lagint(9, eph_blocks[ii][:, 0], \
                                       eph_gcrs_blocks[ii][:, jj], th[kk])
                gtrs[kk, jj], _ = lagint(9, eph_blocks[ii][:, 0], \
                                       eph_gtrs_blocks[ii][:, jj], th[kk])

#            gtrs[:,jj], _ = lagint(9, eph_blocks[:,0], eph_gtrs[:,jj], th)

    astrotime_out = Time(map(str, ob.tstamps), format='iso', scale='tdb')
#    print _time()-tic  
#    print bcrs
#    raw_input('bugagaga')
    return astrotime_out, bcrs, gcrs, gtrs 
 
 
'''
#==============================================================================
# Load time slots needed from ESA orb file
#==============================================================================
'''
#def load_slots_obsolete(orb_file, t_start, t_stop, t_step, source, inp=None):
#    '''
#    The orbital data proper are just lines providing at discrete time steps the
#    epoch of the state, the state (position in km, velocity in km/s) and, if
#    applicable, the state derivative (w.r.t time scale in days).
#    
#    cbody = 'venus' for VEX, 'mars' for MEX
#    '''
#    with open(orb_file,'r') as f:
#        f_lines = f.readlines()
#    
#    slots = [ii for ii,v in enumerate(f_lines) if 'START_TIME' in v]
#    
#    slotData = []
#    r = re.compile('[-T:]+')
#    
#    for slot in slots:
#        out = [slot]
#        for jj in (0,1):
#            eq = f_lines[slot+jj].index('=')            
#            t = r.split(f_lines[slot+jj][eq+1:].strip())
#            ttag = map(int,t[:-1])
#            ttag.append(float(t[-1]))
#            out.append(ttag)
##            slot_start = datetime(*map(int,ttag))
#        slotData.append(out)
#    
#    # line number of 'start_time', start_time, stop_time
#    s = [x for x in slotData if (mjuliandate(*t_start)>=mjuliandate(*x[1]) and\
#                                mjuliandate(*t_start)<=mjuliandate(*x[2])) or\
#                                (mjuliandate(*t_start)<=mjuliandate(*x[1]) and\
#                                (mjuliandate(*t_stop)<=mjuliandate(*x[2]) and\
#                                 mjuliandate(*t_stop)>=mjuliandate(*x[1]))) or\
#                                 (mjuliandate(*t_start)<=mjuliandate(*x[1]) and\
#                                  mjuliandate(*t_stop)>=mjuliandate(*x[2]))]
##    print s
#    # load blocks of ephm corresponding to s'es
#    eph_blocks = []
#    # count not from block start (commented below), but t_start
#    date_start = datetime.datetime(*t_start[0:3])
#    for ii, sta, sto in s:
##        date_start = datetime.datetime(*sta[0:3])
##        print date_start
#        kk = ii + 8
#        block = []
#        while len(f_lines[kk])>5 and ('META_START' not in f_lines[kk]):
#            # skip state derivatives:
#            if f_lines[kk][0]==' ': 
#                kk += 1
#                continue
#            line_split = f_lines[kk].split()
#            ttag = map(float, r.split(line_split[0]))
#            dd = (datetime.datetime(*map(int,ttag[0:3])) - date_start).days
#            ct = ttag[3] + ttag[4]/60.0 + ttag[5]/3600.0 + dd*24.0
##            print flatten([ct, map(float,line_split[1:])])
##            raw_input()
#            block.append(flatten([ct, map(float,line_split[1:])]))
#            kk += 1
#        eph_blocks.append( np.array(block) )
#
#
#    ''' make the ephemeris for output '''
#    # fake obs
#    ob = obs(['DUMMY'],'DUMMY','C')
#    date_t_start = datetime.datetime(*t_start)
#    date_t_stop = datetime.datetime(*t_stop)
#    ob.addScan(date_t_start, t_step, stop=date_t_stop)
#    
#    ''' load input sittings: '''
#    if inp==None:
#        inp = inp_set('inp.cfg')
#    jpl_eph = inp.jpl_eph
#    cat_eop = inp.cat_eop
##    jpl_eph = '/Users/oasis/_jive/python/vispy/jpl_eph/JPLEPH.421'
##    jpl_eph = '/Users/oasis/_jive/python/vispy/jpl_eph/JPLEPH.405'
##    cat_eop = '/Users/oasis/_jive/python/vispy/cats/eopc04.cat'
#    
#    #load eops
#    ''' get the relevant eop entries from the catalogue: '''
#    mjd_start = mjuliandate(date_t_start.year, date_t_start.month,\
#                            date_t_start.day, date_t_start.hour, \
#                            date_t_start.minute,date_t_start.second)
#    with open(cat_eop, 'r') as fc:
#        fc_lines = fc.readlines()
#    eops = np.zeros((7,7)) # +/- 3 days
#    for jj in range(len(fc_lines)):
#        if fc_lines[jj][0]!=' ' and fc_lines[jj][0]!='*':
#            entry = [float(x) for x in fc_lines[jj].split()]
#            if len(entry) > 0 and entry[3] == np.floor(mjd_start) - 3:
#                for kk in range(7):
#                    eops[kk,0] = entry[3] # mjd
#                    eops[kk,1] = entry[6] # UT1-UTC
#                    eops[kk,2] = entry[6] - nsec(entry[3]) # UT1 - TAI
#                    eops[kk,3] = entry[4] # x_p
#                    eops[kk,4] = entry[5] # y_p
#                    eops[kk,5] = entry[8] # dX
#                    eops[kk,6] = entry[9] # dY
#                    entry = [float(x) for x in fc_lines[jj+kk+1].split()]
#                break #exit loop
#    
#    eph_bcrs = []
#    eph_gcrs = []
#    eph_gtrs = []
#    mjd0 = mjuliandate(ob.tstamps[0].year,\
#                        ob.tstamps[0].month, ob.tstamps[0].day)
#    for t in ob.tstamps:
#        # dd is counted from 0h of date_t_start
#        mjd = mjuliandate(t.year, t.month, t.day)
#        dd = mjd - mjd0
#        # CT in hours from mjd
#        ct_h = t.hour + t.minute/60.0 + t.second/3600.0 + dd*24.0
#        mjd_full = mjuliandate(t.year,t.month,t.day,t.hour,t.minute,t.second)
#        
#        jd = mjd + 2400000.5
#        ct = (t.hour*3600.0 + t.minute*60.0 + t.second)/86400.0
#        
#        # choose proper eph block
#        ii_cur = [ii for ii,v in enumerate(s) if \
#                        mjuliandate(*v[1]) <= mjd_full < mjuliandate(*v[2])][0]
#        
##        n = 8 # Lagrange poly order
##        # cut n+1 points around ct_h
##        ctz = tuple(eph_blocks[ii_cur][:,0])
##        position = np.searchsorted(ctz, ct_h)
###        print position
##        # number of points to cut from the left-hand side
##        n_left = int(floor((n+1)/2.0))
##        # number of points to cut from the right-hand side
##        n_right = int(ceil((n+1)/2.0))
##    
###        print n_left, n_right
##        
##        # check/correct bounds:
###        print len(ctz[:position])
##        if len(ctz[:position]) < n_left:
##            n_right += n_left - len(ctz[:position])
##            n_left = len(ctz[:position])
###        print len(ctz[position:])
##        if len(ctz[position:]) < n_right:
##            n_left += n_right - len(ctz[position:])
##            n_right = len(ctz[position:])
##
###        print n_left, n_right
##        
##        # cut the proper piece:
##        ctz = ctz[position-n_left:position+n_right]
##        # cut proper pieces of r/v:
##        rv = []
##        for jj in range(1,7):
##            yp = tuple(eph_blocks[ii_cur][:,jj])
##            yp = yp[position-n_left:position+n_right]
##            # interpolate to current CT
##            L = make_lag_bi(ctz, yp)
##            rv.append(L(ct_h))
##        rv = np.array(rv)
#        
#        # this is an order of magnitude faster with the same precision:
#        rv = []
#        for jj in range(1,7):
#            x, _ = lagint(9, eph_blocks[ii_cur][:,0],\
#                             eph_blocks[ii_cur][:,jj], ct_h)
#            rv.append(x)
#        rv = np.array(flatten(rv))
#        
#        # add barycentric r/v of cbody at CT:
#        if source.lower()=='vex':
#            cbodyState = pleph(jd+ct, 2, 12, jpl_eph)
#        elif source.lower()=='mex':
#            cbodyState = pleph(jd+ct, 4, 12, jpl_eph)
#        # append line to output SS barycentric ephemeris of VEX:
#        rv_bcrs = rv + cbodyState
#        eph_bcrs.append(rv_bcrs)
#        
#        # calculate geocentric orbits
##        dtdbtt = dtdb( jd, 0.0, 0.5, 0.0, 0.0, 0.0 )
#        # GCRS first:
#        earth = pleph(jd+ct, 3, 12, jpl_eph)
#        rv_gcrs = rv_bcrs - earth
#        eph_gcrs.append(rv_gcrs)
#        
#        # then GTRS
#        astrotime = Time(str(t), format='iso', scale='tdb', precision=9)
#        UTCdatetime = astrotime.utc.datetime
#        UTC = (UTCdatetime.hour*3600.0 + \
#               UTCdatetime.minute*60.0 + \
#               UTCdatetime.second + UTCdatetime.microsecond*1e-6)/86400.0
#
#        ''' interpolate eops to t_obs '''
#        UT1, eop_int = eop_iers(mjd, UTC, eops)
#        
#        ''' rotation matrix IERS '''
#        r2000 = ter2cel(UTCdatetime, eop_int, mode='der')
#        rv_gtrs = np.zeros(6)
#        rv_gtrs[0:3] = dot(r2000[:,:,0].T, rv_gcrs[0:3].T)
#        rv_gtrs[3:6] = dot(r2000[:,:,0].T, rv_gcrs[3:6].T) + \
#                       dot(r2000[:,:,1].T, rv_gcrs[0:3].T)
#        eph_gtrs.append(rv_gtrs)
#        
#        
#    # turn into numpy array:
#    eph_bcrs = np.array(eph_bcrs)
#    eph_gcrs = np.array(eph_gcrs)
#    eph_gtrs = np.array(eph_gtrs)
#    
#    astrotime_out = Time(map(str, ob.tstamps), format='iso', scale='tdb')
#
#    return astrotime_out, eph_bcrs, eph_gcrs, eph_gtrs

'''
#==============================================================================
# Create vispy eph-files for ESA spacecraft
# [download from https://www.fd-tasc.info]
#==============================================================================
'''
def esa_eph_download_helper(args):
    sc_name, seg_start, seg_stop, ref_object, frame, scale = args
    # print(sc_name, seg_start, seg_stop, ref_object, frame, scale)
    known_spacecraft = {'vex': 'Venus-Express', 'mex': 'Mars-Express', 'gai': 'Gaia'}
    eph_url = 'https://www.fd-tasc.info/{:s}cmd-cgi-bin/seqgenExec.pl?'.format(sc_name) + \
              'queryType=run&ops=ops_Routine&' + \
              'job={:s}_job_tmpl_fdtool_crttab_state.dat&'.format(sc_name[0]) + \
              'source=seqgenExec&opsConfig=opsSelect&jobConfig=jobTable&TIME_SCALE={:s}&'.format(scale) + \
              'INITIAL_TIME={:s}&'.format(datetime.datetime.strftime(seg_start,
                                                                     '%Y/%m/%d+%H:%M:%S.000')) + \
              'FINAL_TIME={:s}&'.format(datetime.datetime.strftime(seg_stop,
                                                                   '%Y/%m/%d+%H:%M:%S.000')) + \
              'TIME_STEP=000+00:00:01.000&' + \
              'OBJECT={:s}&'.format(known_spacecraft[sc_name]) + \
              'REF_OBJECT={:s}&'.format(ref_object) + \
              'FRAME={:s}&'.format(frame) + \
              'LTC_FLAG=NO' + \
              '&scenarioId=NEW_ID&inputCase=DEFAULT_case'

    # these guys are just unable to unify things somehow...
    if sc_name == 'vex':
        eph_url = eph_url.replace('_Routine', '')
    try:
        # get and parse response
        # response = requests.get(eph_url, timeout=5)
        response = requests.get(eph_url)
        html_seg = response.text
        # print html
        html_seg = html_seg.split('\n')
        html_seg = [l for l in html_seg if len(l) > 5 and (l.strip()[0] != '#' and l.strip()[0] != '<')]
    # handle errors
    except urllib2.HTTPError, e:
        print "HTTP Error:", e.code, eph_url
        html_seg = []
    except urllib2.URLError, e:
        print "URL Error:", e.reason, eph_url
        html_seg = []

    # save current segment
    return html_seg


def esa_sc_eph_make(sc_name, start, stop, inp, paddLeft=30, paddRight=2, parallel=True):
    """
    Let ESOC do all the transformations for us.

    Note that they're using different JPL ephemerides for different spacecraft (e.g. DE405 for MEX and VEX)

    This should be somehow addressed in the code revisions to come,
    just as it was done when ESOC used to provide raw orbits with the actual
    (varying) t_int steps and in the system where integration was carried out

    :param sc_name:
    :param start:
    :param stop:
    :return:
    """

    sc_name = sc_name.lower()

    known_spacecraft = {'vex': 'Venus-Express', 'mex': 'Mars-Express', 'gai': 'Gaia'}
    # 'ros':'Rosetta', 'lpf':'', 'pla':'', 'her':''
    if sc_name not in known_spacecraft:
        raise Exception('Spacecraft unknown')

    # beidseitig padding:
    t_start = start - datetime.timedelta(minutes=paddLeft)
    t_stop = stop + datetime.timedelta(minutes=paddRight)

    # date string:
    date_string = start.strftime("%y%m%d")

    # file names
    sc_bcrs_eph = ''.join((sc_name, '.bcrs.tdb.', date_string, '.eph'))
    sc_gcrs_eph = ''.join((sc_name, '.gcrs.utc.', date_string, '.eph'))
    sc_gtrs_eph = ''.join((sc_name, '.gtrs.utc.', date_string, '.eph'))

    # check whether the file with BCRS eph exists
    if os.path.isfile(os.path.join(inp['sc_eph_cat'], sc_bcrs_eph)):
        eph_bc_exist = True
    else:
        eph_bc_exist = False
    # if it does, check that T_obs is within the time slot of the existing ephemeris
    t_obs_out_of_eph_boundary = False
    if eph_bc_exist:
        try:
            # load it:
            with open(os.path.join(inp['sc_eph_cat'], sc_bcrs_eph), 'r') as f:
                tmp = f.readlines()
            firstLine = map(int, [float(x) for x in tmp[0].split()[0:6]])
            lastLine = map(int, [float(x) for x in tmp[-1].split()[0:6]])
            eph_start = datetime.datetime(*firstLine)
            eph_stop = datetime.datetime(*lastLine)
            if (t_start < eph_start) or (t_stop > eph_stop):
                t_obs_out_of_eph_boundary = True
        except:
            t_obs_out_of_eph_boundary = True

    # force update requested?
    if eph_bc_exist and inp['sc_eph_force_update']:
        print('Precalculated S/C ephs forced update requested.')
        print('removing {:s}'.format(sc_bcrs_eph))
        os.remove(os.path.join(inp['sc_eph_cat'], sc_bcrs_eph))
        print('removing {:s}'.format(sc_gcrs_eph))
        os.remove(os.path.join(inp['sc_eph_cat'], sc_gcrs_eph))
        print('removing {:s}'.format(sc_gtrs_eph))
        os.remove(os.path.join(inp['sc_eph_cat'], sc_gtrs_eph))

    # now make/update the eph
    if not eph_bc_exist or t_obs_out_of_eph_boundary:
        if not eph_bc_exist:
            print('S/C BCRS ephemeris file ' + sc_bcrs_eph + ' not found, creating...')
        if t_obs_out_of_eph_boundary:
            #            print t_start, t_end, eph_start, eph_end
            print('T_obs not within existing BCRS ephemeris time range. Updating {:s}...'.format(sc_bcrs_eph))

        # start time is in the future? notify the user!
        if t_start > datetime.datetime.now():
            print('Note that start date is in the future! ' +
                  'Using planning ephs. \n' + 'Force update computed ephs when final version is available!')

        '''create inputs for the helper function'''
        inps_bcrs = []
        inps_gcrs = []
        inps_gtrs = []
        for seg in range(1, 13):
            seg_start = t_start + datetime.timedelta(hours=2 * (seg - 1))
            seg_stop = t_start + datetime.timedelta(hours=2 * seg) - datetime.timedelta(seconds=1)
            ref_object = 'Solar+System+Barycentre'
            frame = 'mean+equatorial+J2000'
            scale = 'TDB'
            inps_bcrs.append([sc_name, seg_start, seg_stop, ref_object, frame, scale])
            ref_object = 'Earth'
            scale = 'UTC'
            inps_gcrs.append([sc_name, seg_start, seg_stop, ref_object, frame, scale])
            frame = 'Earth+fixed'
            inps_gtrs.append([sc_name, seg_start, seg_stop, ref_object, frame, scale])
        inps = inps_bcrs + inps_gcrs + inps_gtrs

        if parallel:  # Parallel way
            n_cpu = mp.cpu_count()
            # create a Thread pool pool. a Process pool silently fails for some reason
            pool = mp.pool.ThreadPool(np.min((n_cpu, len(inps))))  #
            # asynchronously apply helper to each of inps
            '''bcrs, gcrs, and gtrs'''
            result = pool.map_async(esa_eph_download_helper, inps)
            # close bassejn
            pool.close()
            # tell it to wait until all threads are done before going on
            pool.join()
            # get the ordered results
            allstacked = np.hstack(np.array(result.get()))
            bcrs = allstacked[:86400]
            gcrs = allstacked[86400:86400 * 2]
            gtrs = allstacked[86400 * 2:]
        else:  # Serial way
            bcrs = []
            gcrs = []
            gtrs = []
            for _inp in inps_bcrs:
                bcrs.append(esa_eph_download_helper(_inp))
            for _inp in inps_gcrs:
                gcrs.append(esa_eph_download_helper(_inp))
            for _inp in inps_gtrs:
                gtrs.append(esa_eph_download_helper(_inp))
            bcrs = np.hstack(np.array(bcrs))
            gcrs = np.hstack(np.array(gcrs))
            gtrs = np.hstack(np.array(gtrs))

        ''' dump to files '''
        tstamps = [datetime.datetime.strptime(t.split()[0], '%Y/%m/%dT%H:%M:%S.%f') for t in bcrs]
        eph_bcrs = np.array([map(float, s.split()[1:]) for s in bcrs])
        eph_gcrs = np.array([map(float, s.split()[1:]) for s in gcrs])
        eph_gtrs = np.array([map(float, s.split()[1:4]) + [0, 0, 0] for s in gtrs])
        ''' BCRS '''
        with open(os.path.join(inp['sc_eph_cat'], sc_bcrs_eph), 'w') as f:
            for ii, t in enumerate(tstamps):
                #            line = '{:4d}-{:02d}-{:02d}T{:02d}:{:02d}:{:06.3f}'.\
                line = '{:4d} {:02d} {:02d} {:02d} {:02d} {:06.3f}'. \
                    format(t.year, t.month, t.day, t.hour, t.minute, t.second)
                line += '{:20.6f} {:18.6f} {:18.6f} {:16.10f} {:15.10f} {:15.10f}\n' \
                    .format(*eph_bcrs[ii, 0:6])
                f.write(line)

        ''' GCRS '''
        with open(os.path.join(inp['sc_eph_cat'], sc_gcrs_eph), 'w') as f:
            for ii, t in enumerate(tstamps):
                line = '{:4d} {:02d} {:02d} {:02d} {:02d} {:06.3f}' \
                    .format(t.year, t.month, t.day, t.hour, t.minute, t.second)
                line += '{:20.6f} {:18.6f} {:18.6f} {:16.10f} {:15.10f} {:15.10f} ' \
                    .format(*eph_gcrs[ii, 0:6])
                line += '{:16.10f} {:15.10f} {:15.10f}\n' \
                    .format(*np.zeros(3))
                f.write(line)

        ''' GTRS '''
        with open(os.path.join(inp['sc_eph_cat'], sc_gtrs_eph), 'w') as f:
            for ii, t in enumerate(tstamps):
                line = '{:4d} {:02d} {:02d} {:02d} {:02d} {:06.3f}'. \
                    format(t.year, t.month, t.day, t.hour, t.minute, t.second)
                line += '{:20.6f} {:18.6f} {:18.6f} {:16.10f} {:15.10f} {:15.10f}\n' \
                    .format(*eph_gtrs[ii, 0:6])
                f.write(line)

        return sc_bcrs_eph, sc_gtrs_eph, sc_gcrs_eph


'''
#==============================================================================
# Create vispy eph-files from ESA orb file parsed/loaded with load_slots(*args)
#==============================================================================
'''
def esa_sc_eph_make_from_raw(source, date_t_start, date_t_stop, inp,
                             paddLeft=30, paddRight=0):
    '''
    Make vispy eph-files from ESA orb file parsed/loaded with load_slots(*args)
    
    padding - in minutes for each side

    NOTE from 03/05/2016: these raw files are not available anymore
    since the old ESA TASC site was turned down, and the new one does not have those
    '''
    # 30 min beidseitig padding:
    t_start = date_t_start - datetime.timedelta(minutes=paddLeft)
    t_end = date_t_stop + datetime.timedelta(minutes=paddRight)

    # date string:
    date_string = date_t_start.strftime("%y%m%d")

    # file names    
    sc_bcrs_eph = ''.join((source.lower(),'.bcrs.tdb.', date_string, '.eph'))
    sc_gcrs_eph = ''.join((source.lower(),'.gcrs.utc.', date_string, '.eph'))
    sc_gtrs_eph = ''.join((source.lower(),'.gtrs.utc.', date_string, '.eph'))
    
    # check whether the file with BCRS eph exists
    if os.path.isfile(inp['sc_eph_cat']+'/'+sc_bcrs_eph):
        eph_bc_exist = True
    else:
        eph_bc_exist = False
    # if it does, check that T_obs is within the time slot of the existing ephemeris
    t_obs_out_of_eph_boundary = False
    if eph_bc_exist:
        # load it:
        with open(inp['sc_eph_cat']+'/'+sc_bcrs_eph,'r') as f:
            tmp = f.readlines()
        firstLine = map(int, [float(x) for x in tmp[0].split()[0:6]])
        lastLine  = map(int, [float(x) for x in tmp[-1].split()[0:6]])
        eph_start = datetime.datetime(*firstLine)
        eph_end = datetime.datetime(*lastLine)
        if (t_start<eph_start) or (t_end>eph_end):
            t_obs_out_of_eph_boundary = True            
    
    # force update requested?
    if eph_bc_exist and inp['sc_eph_force_update']:
        print 'Precalculated S/C ephs forced update requested.'
        print 'removing {:s}'.format(sc_bcrs_eph)
        os.remove(os.path.join(inp['sc_eph_cat'], sc_bcrs_eph))
        print 'removing {:s}'.format(sc_gcrs_eph)
        os.remove(os.path.join(inp['sc_eph_cat'], sc_gcrs_eph))
        print 'removing {:s}'.format(sc_gtrs_eph)
        os.remove(os.path.join(inp['sc_eph_cat'], sc_gtrs_eph))
        
    
    #now make/update the eph
    if not eph_bc_exist or t_obs_out_of_eph_boundary:
        if not eph_bc_exist:
            print 'S/C BCRS ephemeris file '+sc_bcrs_eph+' not found, creating...'
        if t_obs_out_of_eph_boundary:
#            print t_start, t_end, eph_start, eph_end
            print 'T_obs not within existing BCRS ephemeris time range. Updating '\
                  +sc_bcrs_eph+'...'
        
        # start time is in the future? notify the user!
        if t_start > datetime.datetime.now():
            print 'Note that start date is in the future! '+\
                  'Using planning ephs. \n'+\
                  'Force update computed ephs when final version is available!'
        
        # time slot
        t_start = [t_start.year, t_start.month, t_start.day,\
                   t_start.hour, t_start.minute, t_start.second]
        t_stop  = [t_end.year, t_end.month, t_end.day,\
                   t_end.hour, t_end.minute, t_end.second]
        
        t_step = 1 # seconds
    
        path = os.path.join(inp['sc_eph_cat'], 'raw_'+source.lower())
        # downlaod fresh ephs from tasc
        boundaries = checkbound(source, orb_path=path)
        try:
            eph_file = os.path.join(path, [x[0] for x in boundaries if \
                        mjuliandate(*x[1]) < mjuliandate(*t_start) < \
                        mjuliandate(*t_stop) < mjuliandate(*x[2]) ][-1] )
        except:
            raise Exception('No suitable raw ephemeris found. Fail!')

        astrotime, eph_bcrs, eph_gcrs, eph_gtrs = \
                 load_slots(eph_file, t_start, t_stop, t_step, source, inp)

        ''' interpolate geocentric ephs to utc time stamps'''
        utz = astrotime.utc.datetime
        t_s = date_t_start - datetime.timedelta(minutes=paddLeft)
        t0 = datetime.datetime(t_s.year, t_s.month, t_s.day)
        ut = np.array([(float((tu-t0).days) + (tu.hour*3600.0 + tu.minute*60.0 + \
              tu.second + tu.microsecond*1e-6)/86400.0) for tu in utz])

        ctz = astrotime.tdb.datetime
        ct = np.array([(float((tc-t0).days) + (tc.hour*3600.0 + tc.minute*60.0 + \
              tc.second + tc.microsecond*1e-6)/86400.0) for tc in ctz])
        for jj in range(6):
            eph_gcrs[:,jj] = extrap(ct, ut, eph_gcrs[:,jj])
            eph_gtrs[:,jj] = extrap(ct, ut, eph_gtrs[:,jj])
#            eph_gcrs[:,jj] = np.interp(ct, ut, eph_gcrs[:,jj])
#            eph_gtrs[:,jj] = np.interp(ct, ut, eph_gtrs[:,jj])
#            eph_gcrs[:,jj], _ = lagint(9, ut, eph_gcrs[:,jj], ct)
#            eph_gtrs[:,jj], _ = lagint(9, ut, eph_gtrs[:,jj], ct)
        
        ''' output to files '''
        tstamps = astrotime.tdb.datetime # t stamps look the same in utc and tdb
    
        ''' BCRS '''
        with open(os.path.join(inp['sc_eph_cat'], sc_bcrs_eph),'w') as f:
            for ii, t in enumerate(tstamps):
    #            line = '{:4d}-{:02d}-{:02d}T{:02d}:{:02d}:{:06.3f}'.\
                line = '{:4d} {:02d} {:02d} {:02d} {:02d} {:06.3f}'.\
                        format(t.year,t.month,t.day,t.hour,t.minute,t.second)
                line += '{:20.6f} {:18.6f} {:18.6f} {:16.10f} {:15.10f} {:15.10f}\n'\
                        .format(*eph_bcrs[ii,0:6])
                f.write(line)
                
        ''' GCRS '''
        with open(os.path.join(inp['sc_eph_cat'], sc_gcrs_eph),'w') as f:
            for ii, t in enumerate(tstamps):
                line = '{:4d} {:02d} {:02d} {:02d} {:02d} {:06.3f}'\
                         .format(t.year,t.month,t.day,t.hour,t.minute,t.second)
                line += '{:20.6f} {:18.6f} {:18.6f} {:16.10f} {:15.10f} {:15.10f} '\
                         .format(*eph_gcrs[ii,0:6])
                line += '{:16.10f} {:15.10f} {:15.10f}\n'\
                         .format(*np.zeros(3))
                f.write(line)
        
        ''' GTRS '''
        with open(os.path.join(inp['sc_eph_cat'], sc_gtrs_eph),'w') as f:
            for ii, t in enumerate(tstamps):
                line = '{:4d} {:02d} {:02d} {:02d} {:02d} {:06.3f}'.\
                        format(t.year,t.month,t.day,t.hour,t.minute,t.second)
                line += '{:20.6f} {:18.6f} {:18.6f} {:16.10f} {:15.10f} {:15.10f}\n'\
                        .format(*eph_gtrs[ii,0:6])
                f.write(line)

    return (sc_bcrs_eph, sc_gtrs_eph, sc_gcrs_eph)

'''
#==============================================================================
# Download ESA s/c ephemerides from tasc.esa.int
#==============================================================================
'''
def esa_sc_eph_down(source, date_t_start, date_t_end, inp):
    orbtyp = 'ops'

    if source=='VEX': staobj='Venus-Express'
    if source=='MEX': staobj='Mars-Express'
    if source=='HER': staobj='Herschel'

    sc_eph_cat = inp['sc_eph_cat'] # location of the dir with s/c eph cat
    
    ''' Barycentric ephemeris '''    
    # 30 min beidseitig padding:
    t_start = date_t_start - datetime.timedelta(seconds=30*60)
    t_end = date_t_end + datetime.timedelta(seconds=30*60)
    
    # because of a stupid limitation on the website of 10000 maximum epochs
    # in one go...:
    if source=='MEX' or source=='HER': 
        t_step_bc=5
    else:
        t_step_bc=1
    
    sc_bcrs_eph = source.lower()+'.bcrs.tdb.'+str(t_start.year)[2:]+\
                '%02d'%date_t_start.month+'%02d'%date_t_start.day+'.eph'
    # check whether the file with eph exists
    if os.path.isfile(sc_eph_cat+'/'+sc_bcrs_eph):
        eph_bc_exist = 1
    else:
        eph_bc_exist = 0
    # if it does, check that T_obs is within the time slot of the existing ephemeris
    t_obs_out_of_eph_boundary = 0
    if eph_bc_exist:
        # load it:
        with open(sc_eph_cat+'/'+sc_bcrs_eph,'r') as f:
            tmp = f.readlines()
        eph = []
        # do what the matlab function load does (convert strings to floats):
        for line in tmp:
            line = [float(x) for x in line.split()]
            line[0:6] = map(int,line[0:6])
            eph.append(line)
        eph_start = datetime.datetime(*eph[0][0:6])
        eph_end = datetime.datetime(*eph[-1][0:6])
        if (t_start<eph_start) or (t_end>eph_end):
            t_obs_out_of_eph_boundary = 1            
            
    #now download/update the eph
    if eph_bc_exist!=1 or t_obs_out_of_eph_boundary==1:
        if eph_bc_exist!=1:
            print 'S/C bc ephemeris file: '+sc_bcrs_eph+' not found, downloading...'
        if t_obs_out_of_eph_boundary==1:
            print 'T_obs is not within the existing bc ephemeris time range. Updating '+sc_bcrs_eph+'...'
        
        eph_url = 'http://tasc.esa.int/cgi-bin/query.html?'+\
                'mission='+source+'&'+\
                'querytyp=data&'+\
                'qtyp=sta&'+\
                'staobj='+staobj+'&'+\
                'starefobj=Solar+System+Barycentre&'+\
                'staltc=NO&'+\
                'stafrm=mean+equatorial+J2000&'+\
                'matfrm=mean+equatorial+J2000&'+\
                'quafrm=mean+equatorial+J2000&'+\
                'angobj1=Earth&'+\
                'angobj2=Sun&'+\
                'angrefobj='+staobj+'&'+\
                'tscl=TDB&'+\
                'tstart='+str(t_start.year)+\
                            '%2F'+'%02d'%t_start.month+\
                            '%2F'+'%02d'%t_start.day+'+'+\
                            '%02d'%t_start.hour+'%3A'+\
                            '%02d'%t_start.minute+'%3A'+\
                            '%06.3f'%t_start.second+'&'+\
                'tend='+str(t_end.year)+\
                            '%2F'+'%02d'%t_end.month+\
                            '%2F'+'%02d'%t_end.day+'+'+\
                            '%02d'%t_end.hour+'%3A'+\
                            '%02d'%t_end.minute+'%3A'+\
                            '%06.3f'%t_end.second+'&'+\
                'tstp=000+00%3A00%3A'+'%06.3f'%t_step_bc+'&orbtyp='+\
                orbtyp
        #print eph_url                
        response = urllib2.urlopen(eph_url)
        html = response.read()
        #print html
        html = html.split('\n')
        html_copy = []
        for line in html:
            if len(line)>5 and (line.strip()[0]!='#' and line.strip()[0]!='<'):
                for char in (':','T','/'):             
                    line = line.replace(char,' ')
                html_copy.append(line)
        # print it to file:
        with open(os.path.join(sc_eph_cat,sc_bcrs_eph),'w') as out:
            for line in html_copy:
                out.write(line+'\n')
        
        
    ''' Geocentric ephemeris without lt-correction '''
    # if no need for padding:
    #t_start = date_t_start
    #t_end = date_t_end
    # if 30 min beidseitig padding - keep the same t_start and t_end
    
    #t_step_gc = 10 #seconds
    
    # because of a stupid limitation on the website of 10000 maximum epochs
    # in one go...:
    if source=='MEX' or source=='HER': 
        t_step_gc=5
    else:
        t_step_gc=1
    
    #sc_gtrs_eph = source.lower()+'.gc.ut.ltcorr'+str(t_start.year)[2:]+\
    sc_gtrs_eph = source.lower()+'.gtrs.utc.'+str(t_start.year)[2:]+\
                '%02d'%date_t_start.month+'%02d'%date_t_start.day+'.eph'
    # check whether the file with eph exists
    if os.path.isfile(sc_eph_cat+'/'+sc_gtrs_eph):
        eph_gc_exist = 1
    else:
        eph_gc_exist = 0
    # if it does, check that T_obs is within the time slot of the existing ephemeris
    t_obs_out_of_eph_boundary = 0
    if eph_gc_exist:
        # load it:
        f = open(sc_eph_cat+'/'+sc_gtrs_eph,'r')
        tmp = f.readlines()
        eph = []
        # do what the matlab function load does (convert strings to floats):
        for line in tmp:
            line = [float(x) for x in line.split()]
            line[0:6] = map(int,line[0:6])
            eph.append(line)
        eph_start = datetime.datetime(eph[0][0],eph[0][1],\
                                      eph[0][2],eph[0][3],eph[0][4],eph[0][5])
        eph_end = datetime.datetime(eph[-1][0],eph[-1][1],\
                                    eph[-1][2],eph[-1][3],eph[-1][4],eph[-1][5])
        if (t_start<eph_start) or (t_end>eph_end):
            t_obs_out_of_eph_boundary = 1            
            
    #now download/update the GTRS eph
    if eph_gc_exist!=1 or t_obs_out_of_eph_boundary==1:
        if eph_gc_exist!=1:
            print 'S/C gtrs ephemeris file: '+sc_gtrs_eph+' not found, downloading...'
        if t_obs_out_of_eph_boundary==1:
            print 'T_obs is not within the existing gtrs ephemeris time range. Updating '+sc_gtrs_eph+'...'

        # DO NOT LT-CORRECT!!
        eph_url = 'http://tasc.esa.int/cgi-bin/query.html?'+\
                'mission='+source+'&'+\
                'querytyp=data&'+\
                'qtyp=sta&'+\
                'staobj='+staobj+'&'+\
                'starefobj=Earth&'+\
                'staltc=NO&'+\
                'stafrm=Earth+fixed&'+\
                'matfrm=mean+equatorial+J2000&'+\
                'quafrm=mean+equatorial+J2000&'+\
                'angobj1=Earth&'+\
                'angobj2=Sun&'+\
                'angrefobj='+staobj+'&'+\
                'tscl=UTC&'+\
                'tstart='+str(t_start.year)+\
                            '%2F'+'%02d'%t_start.month+\
                            '%2F'+'%02d'%t_start.day+'+'+\
                            '%02d'%t_start.hour+'%3A'+\
                            '%02d'%t_start.minute+'%3A'+\
                            '%06.3f'%t_start.second+'&'+\
                'tend='+str(t_end.year)+\
                            '%2F'+'%02d'%t_end.month+\
                            '%2F'+'%02d'%t_end.day+'+'+\
                            '%02d'%t_end.hour+'%3A'+\
                            '%02d'%t_end.minute+'%3A'+\
                            '%06.3f'%t_end.second+'&'+\
                'tstp=000+00%3A00%3A'+'%06.3f'%t_step_gc+'&orbtyp='+\
                orbtyp
        #print eph_url                
        response = urllib2.urlopen(eph_url)
        html = response.read()
        #print html
        html = html.split('\n')
        html_copy = []
        for line in html:
#            print len(line)
            if len(line)>2 and (line.strip()[0]!='#' and line.strip()[0]!='<'):
                for char in (':','T','/'):             
                    line = line.replace(char,' ')
                html_copy.append(line[0:75])
        
        eph = []
        # do what the matlab function load does (convert strings to floats):
        for line in html_copy:
            line = [float(x) for x in line.split()]
            line[0:5] = map(int,line[0:5])
            eph.append(line)
        
        UT_eph = [(ep[3] + ep[4]/60.0 + ep[5]/3600.0)/24.0 for ep in eph]
        # allow for 'overnighters':
        dd = mjuliandate(*eph[-1][0:3])-mjuliandate(*eph[0][0:3])
        for cr in range(1,int(dd)+1):
            for nn in range(1,len(UT_eph)):
                 if UT_eph[nn]<UT_eph[nn-1]:
                      UT_eph[nn] = UT_eph[nn] + 1        
        
        ''' calculate v and a using a nth-order Chebyshёv poly '''
        eph = np.array(eph)
        nth = 7
        v = np.zeros((len(UT_eph),3))
        while nth<12:
            try:
                for jj in range(9,12):
                    p = cheb.chebfit(UT_eph,eph[:,jj-3], nth)
                    v[:,jj-9] = cheb.chebval(UT_eph,
                                cheb.chebder(p)) / 86400.0
            except:
                # didn't work? try a higher order then
                nth += 1
            else:
                # went well? exit loop then
                break
        eph = np.hstack((eph,v))
        
        nth = 7
        a = np.zeros((len(UT_eph),3))
        while nth<12:
            try:
                for jj in range(12,15):
                    p = cheb.chebfit(UT_eph,eph[:,jj-3], nth)
                    a[:,jj-12] = cheb.chebval(UT_eph,
                                 cheb.chebder(p)) / 86400.0
            except:
                # didn't? try a higher order then
                nth += 1
            else:
                # went well? exit loop then
                break
        eph = np.hstack((eph,a))
        
        
        # print it to file:
        f = open(sc_eph_cat+'/'+sc_gtrs_eph,'w')
        for ep in eph:
            datestr = map(int,ep[0:5])
            datestr.append(ep[5])
            s = '{:5d} {:02d} {:02d} {:02d} {:02d} {:06.3f} '\
                 .format(*datestr)+\
                '{:20.6f} {:18.6f} {:18.6f} {:16.10f} {:15.10f} {:15.10f} '\
                .format(*ep[6:12])+\
                '{:16.10f} {:15.10f} {:15.10f}\n'\
                .format(*ep[12:15])
            f.write(s)
        f.close()
        
    # make a GCRS ephemeris:
    sc_gcrs_eph = source.lower()+'.gcrs.utc.'+str(t_start.year)[2:]+\
                  '%02d'%date_t_start.month+'%02d'%date_t_start.day+'.eph'
    # check whether the file with eph exists
    if os.path.isfile(sc_eph_cat+'/'+sc_gcrs_eph):
        eph_gcrs_exist = 1
    else:
        eph_gcrs_exist = 0
    # if it does, check that T_obs is within the time slot of the existing ephemeris
    t_obs_out_of_eph_boundary = 0
    if eph_gcrs_exist:
        # load it:
        f = open(sc_eph_cat+'/'+sc_gcrs_eph,'r')
        tmp = f.readlines()
        eph = []
        # do what the matlab function load does (convert strings to floats):
        for line in tmp:
            line = [float(x) for x in line.split()]
            line[0:6] = map(int,line[0:6])
            eph.append(line)
        eph_start = datetime.datetime(*eph[0][0:6])
        eph_end = datetime.datetime(*eph[-1][0:6])
        if (t_start<eph_start) or (t_end>eph_end):
            t_obs_out_of_eph_boundary = 1            
            
    # now make/update the GCRS eph in UTC time stamps
    # this is used to calculate ra/dec's
    if eph_gcrs_exist!=1 or t_obs_out_of_eph_boundary==1:
        if eph_gcrs_exist!=1:
            print 'S/C gcrs ephemeris file: '+sc_gcrs_eph+' not found, calculating...'
        if t_obs_out_of_eph_boundary==1:
            print 'T_obs is not within the existing gtrs ephemeris time range. Updating '+sc_gcrs_eph+'...'
            
        # load BCRS eph in UTC time stamps:
        eph_url = 'http://tasc.esa.int/cgi-bin/query.html?'+\
                'mission='+source+'&'+\
                'querytyp=data&'+\
                'qtyp=sta&'+\
                'staobj='+staobj+'&'+\
                'starefobj=Earth&'+\
                'staltc=NO&'+\
                'stafrm=mean+equatorial+J2000&'+\
                'matfrm=mean+equatorial+J2000&'+\
                'quafrm=mean+equatorial+J2000&'+\
                'angobj1=Earth&'+\
                'angobj2=Sun&'+\
                'angrefobj='+staobj+'&'+\
                'tscl=UTC&'+\
                'tstart='+str(t_start.year)+\
                            '%2F'+'%02d'%t_start.month+\
                            '%2F'+'%02d'%t_start.day+'+'+\
                            '%02d'%t_start.hour+'%3A'+\
                            '%02d'%t_start.minute+'%3A'+\
                            '%06.3f'%t_start.second+'&'+\
                'tend='+str(t_end.year)+\
                            '%2F'+'%02d'%t_end.month+\
                            '%2F'+'%02d'%t_end.day+'+'+\
                            '%02d'%t_end.hour+'%3A'+\
                            '%02d'%t_end.minute+'%3A'+\
                            '%06.3f'%t_end.second+'&'+\
                'tstp=000+00%3A00%3A'+'%06.3f'%t_step_bc+'&orbtyp='+\
                orbtyp
        #print eph_url                
        response = urllib2.urlopen(eph_url)
        html = response.read()
        #print html
        html = html.split('\n')
        eph = []
        for line in html:
            if len(line)>2 and (line.strip()[0]!='#' and line.strip()[0]!='<'):
                for char in (':','T','/'):             
                    line = line.replace(char,' ')
                eph.append(line)
        
        # save to file:
        f = open(sc_eph_cat+'/'+sc_gcrs_eph,'w')
        for ep in eph:
            ep = [float(x) for x in ep.split()]
            ep[0:5] = map(int,ep[0:5])                                             
#            mjd = mjuliandate(ep[0],ep[1],ep[2])
#            JD = mjd + 2400000.5
#            CT = (ep[3] + ep[4]/60.0 + ep[5]/3600.0)/24.0
            ''' calculate state vector of the Earth at JD+CT, in [km, km/s] '''
#            earth = pleph(JD+CT, 3, 12, inp['jpl_eph'])
#                .format(*(ep[6:12]-earth))+\ 
            s = '{:5d} {:02d} {:02d} {:02d} {:02d} {:06.3f} '\
                 .format(*ep[0:6])+\
                '{:20.6f} {:18.6f} {:18.6f} {:16.10f} {:15.10f} {:15.10f} '\
                .format(*ep[6:12])+\
                '{:16.10f} {:15.10f} {:15.10f}\n'\
                .format(*np.zeros(3))
            f.write(s)
        f.close()
        
    return (sc_bcrs_eph, sc_gtrs_eph, sc_gcrs_eph)


'''
#==============================================================================
# Load time slots from Gaia txt ephemeris file
#==============================================================================
'''
def load_slots_gaia(eph_file, t_start, t_stop, t_step):
    with open(eph_file,'r') as f:
        f_lines = f.readlines()

    slots = [ii for ii,v in enumerate(f_lines) if 'START_TIME' in f_lines[ii]]
    
    slotData = []
    r = re.compile('[-T:]+')
    
    for slot in slots:
        out = [slot]
        for jj in (0,1):
            eq = f_lines[slot+jj].index('=')            
            t = r.split(f_lines[slot+jj][eq+1:].strip())
            ttag = map(int,t[:-1])
            ttag.append(float(t[-1]))
            out.append(ttag)
#            slot_start = datetime(*map(int,ttag))
        slotData.append(out)
    
    # line number of 'start_time', start_time, stop_time
    s = [x for x in slotData if (mjuliandate(*t_start)>=mjuliandate(*x[1]) and\
                                mjuliandate(*t_start)<=mjuliandate(*x[2])) or\
                                (mjuliandate(*t_start)<=mjuliandate(*x[1]) and\
                                (mjuliandate(*t_stop)<=mjuliandate(*x[2]) and\
                                 mjuliandate(*t_stop)>=mjuliandate(*x[1]))) or\
                                 (mjuliandate(*t_start)<=mjuliandate(*x[1]) and\
                                  mjuliandate(*t_stop)>=mjuliandate(*x[2]))]
#    print s
    # load blocks of ephm corresponding to s'es
    eph_blocks = []
    # count not from block start (commented below), but t_start
    date_start = datetime.datetime(*t_start[0:3])
    for ii, sta, sto in s:
#        date_start = datetime.datetime(*sta[0:3])
#        print date_start
        kk = ii + 4
        block = []
        while len(f_lines[kk])>5 and ('META_START' not in f_lines[kk]):
            line_split = f_lines[kk].split()
            ttag = map(float, r.split(line_split[0]))
            dd = (datetime.datetime(*map(int,ttag[0:3])) - date_start).days
            ut = ttag[3] + ttag[4]/60.0 + ttag[5]/3600.0 + dd*24.0
            block.append(flatten([ut, map(float,line_split[1:])]))
            kk += 1
        eph_blocks.append( np.array(block) )


    ''' make the ephemeris for output '''
    # fake obs
    ob = obs(['DUMMY'],'DUMMY','C')
    date_t_start = datetime.datetime(*t_start)
    date_t_stop = datetime.datetime(*t_stop)
    ob.addScan(date_t_start, t_step, stop=date_t_stop)
    
    eph = []
    for t in ob.tstamps:
        # dd is counted from 0h of date_t_start
        dd = (t - date_t_start).days
        ut = t.hour + t.minute/60.0 + t.second/3600.0 + dd*24.0
        mjd = mjuliandate(t.year,t.month,t.day,t.hour,t.minute,t.second)
        
        # choose proper eph block
        ii_cur = [ii for ii,v in enumerate(s) if \
                            mjuliandate(*v[1]) <= mjd < mjuliandate(*v[2])][0]

        n = 8 # Lagrange poly order
        # cut n+1 points around ut
        utz = tuple(eph_blocks[ii_cur][:,0])
        position = np.searchsorted(utz, ut)
#        print position
        # number of points to cut from the left-hand side
        n_left = int(floor((n+1)/2.0))
        # number of points to cut from the right-hand side
        n_right = int(ceil((n+1)/2.0))
    
#        print n_left, n_right
        
        # check/correct bounds:
#        print len(utz[:position])
        if len(utz[:position]) < n_left:
            n_right += n_left - len(utz[:position])
            n_left = len(utz[:position])
#        print len(utz[position:])
        if len(utz[position:]) < n_right:
            n_left += n_right - len(utz[position:])
            n_right = len(utz[position:])

#        print n_left, n_right
        
        # cut the proper piece:
        utz = utz[position-n_left:position+n_right]
        # cut proper pieces of r/v/a:
        rva = []
        for jj in range(1,10):
            yp = tuple(eph_blocks[ii_cur][:,jj])
            yp = yp[position-n_left:position+n_right]
            L = make_lag_bi(utz, yp)
            rva.append(L(ut))
        # append line to output ephemeris
        eph.append(rva)
    # turn it into numpy array
    eph = np.array(eph)
    
    tstamps = [[x.year, x.month, x.day, x.hour, x.minute, x.second] \
                for x in ob.tstamps]
    tstamps = np.array(tstamps)
    
    # return the same thing as load_scf
    return np.hstack((tstamps, eph))
    
'''
#==============================================================================
# Freshest Gaia ephemeris
#==============================================================================
'''
def gaia_fresh(cat_gaia):
    # was a more fresh eph downloaded? should update pyp-ephs, if they exist
    updated = False
    # find 'freshest' txt eph file of Gaia
    for root, subFolders, files in os.walk(cat_gaia):
        txt = [x for x in files if x[-3:] == 'txt']
        fresh = sorted(txt)[-1]
    # if there's internet connection, check if there's newer one on ESA server
    if internet_on():
        try:
            host = 'ssh.sciops.esa.int'
            user = 'rgbot'
            pwd = ''
            port = 22
            
            eph_path = '/home/dpce/dpceInterface/OUT/ORBIT/current/'
            
            client = paramiko.SSHClient()
            client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
            client.connect(hostname=host, username=user, password=pwd, port=port)
            # contents of current dir
            stdin, stdout, stderr = client.exec_command('ls -l ' + eph_path)
            # one line output
            data = stdout.read() + stderr.read()
            # parse 'data' and get the freshest eph file name
            eph_file = [x for x in data.split() if '.txt' in x][-1]
            # close connection
            client.close()
            # check if it's already downloaded, scp if necessary:
            if eph_file not in txt:
                print ''.join(('Downloading ', eph_file, '...'))
                os.system(''.join(('scp rgbot@ssh.sciops.esa.int:',\
                            eph_path, eph_file, ' ', cat_gaia)))
                # refresh:)
                fresh = eph_file
                updated = True
        except Exception, err:
            print str(err)
            print 'Could not check fresh Gaia ephemeris.'
        
    return os.path.join(cat_gaia, fresh), updated
    
'''
#==============================================================================
# Make ephemeris for RA and GAIA and CE3
#==============================================================================
'''
def ra_eph_down(source, date_t_start, date_t_end, inp):
    '''
        Make ephemeris for RA and GAIA and CE3
    '''
    sc_eph_cat = inp['sc_eph_cat']
    cat_eop = inp['cat_eop']
    jpl_eph = inp['jpl_eph']
    
    const = constants()
    
    date_str = str(date_t_start.year)[2:]+\
                '%02d'%date_t_start.month+'%02d'%date_t_start.day
    sc_bcrs_eph = source.lower()+'.bcrs.tdb.'+date_str+'.eph'
    sc_gtrs_eph = source.lower()+'.gtrs.utc.'+date_str+'.eph'
    sc_gcrs_eph = source.lower()+'.gcrs.utc.'+date_str+'.eph'
    # check whether the file with eph exists
    if os.path.isfile(os.path.join(sc_eph_cat, sc_bcrs_eph)):
        eph_bc_exist = True
    else:
        eph_bc_exist = False
    # if it does, check that T_obs is within the time slot of the existing eph
    t_obs_out_of_eph_boundary = False
    if eph_bc_exist:
        # load it:
        with open(os.path.join(sc_eph_cat, sc_bcrs_eph),'r') as f:
            tmp = f.readlines()
        eph = []
        # check boundaries
        start = [float(x) for x in tmp[0].split()]
        eph_start = datetime.datetime(*map(int,start[0:6]))
        stop = [float(x) for x in tmp[-1].split()]
        eph_end = datetime.datetime(*map(int,stop[0:6]))
        if (date_t_start<eph_start) or (date_t_end>eph_end):
            t_obs_out_of_eph_boundary = True            
    
    # force update requested?
    if eph_bc_exist and inp['sc_eph_force_update']:
        print 'Precalculated S/C ephs forced update requested.'
        print 'removing {:s}'.format(sc_bcrs_eph)
        os.remove(os.path.join(inp['sc_eph_cat'], sc_bcrs_eph))
        print 'removing {:s}'.format(sc_gcrs_eph)
        os.remove(os.path.join(inp['sc_eph_cat'], sc_gcrs_eph))
        print 'removing {:s}'.format(sc_gtrs_eph)
        os.remove(os.path.join(inp['sc_eph_cat'], sc_gtrs_eph))
        eph_bc_exist = False
    
    # check Gaia eph updates from ESA
    gaia_updated = False
    if source.lower()=='gaia':
        eph_file, gaia_updated = gaia_fresh(os.path.join(sc_eph_cat, 'raw_gaia'))

    # now download/update the eph
    if not eph_bc_exist or t_obs_out_of_eph_boundary or gaia_updated:
        if not eph_bc_exist:
            print 'S/C bc ephemeris file '+sc_bcrs_eph+' not found'
        if t_obs_out_of_eph_boundary:
            print 'T_obs not within existing bc ephemeris time range '\
                    +sc_bcrs_eph+'...'
        if gaia_updated:
            print 'Found Gaia orbit update'
        
        # set eph file manually:
        #org_eph = load_scf(sc_eph_cat+'/raw_radioastron/RA121031-121101jj.scf')
        
        # default orbit storage:
        if 'org_eph' not in locals(): # var org_eph doesn't exist:
            if source.lower()=='ra':
                try:
#                    org_eph = load_scf(sc_eph_cat+'/raw_radioastron/RA'+date_str+'j.scf')
                    org_eph = load_scf('{:s}/raw_radioastron/RA{:s}j.scf'.\
                                            format(sc_eph_cat, date_str))
                    kk=1;
                    # if observation span across more than 1 day, load more eph(s) and
                    # append it(them) to the org_eph
                    while mjuliandate(*org_eph[-1,0:6]) < \
                          mjuliandate(date_t_end.year,date_t_end.month,date_t_end.day,\
                                      date_t_end.hour,date_t_end.minute,date_t_end.second):
                        new_date = date_t_start + datetime.timedelta(kk)
                        date_str = str(new_date.year)[2:]+\
                                '%02d'%new_date.month+'%02d'%new_date.day
                        org_eph = np.vstack( (org_eph, \
                                  load_scf(sc_eph_cat+'/raw_radioastron/RA'+\
                                  date_str+'j.scf')) )
                        kk+=1
                except IOError:
                    raise Exception('Input orbit for RA was not specified by hand and not'+ \
                           ' found in the /raw_radioastron directory. Fail!')

            if source.lower()[0:2]=='pr' or source.lower()[0:2]=='pg':
                org_eph = load_sp3(sc_eph_cat, source, date_t_start)

            elif source.lower()=='gaia':
                try:
                    #eph_file, _ = gaia_fresh(os.path.join(sc_eph_cat, 'raw_gaia'))
#                    eph_file = sc_eph_cat+'/raw_gaia/ORB1_20140205_000001.txt'
                    t_step = 10.0 #seconds
#                    padd = 10.0 # minutes
                    padd = 0.0 # minutes
                    t_start = date_t_start - datetime.timedelta(minutes=padd)
                    t_stop = date_t_end + datetime.timedelta(minutes=padd)
                    t_start = [t_start.year, t_start.month, t_start.day, \
                               t_start.hour, t_start.minute, t_start.second]
                    t_stop = [t_stop.year, t_stop.month, t_stop.day, \
                               t_stop.hour, t_stop.minute, t_stop.second]
                    org_eph = load_slots_gaia(eph_file, t_start, t_stop, t_step)
#                    org_eph = load_scf(sc_eph_cat+'/raw_gaia/GAIA'\
#                                       +date_str+'j.scf')
                except IOError:
                    raise Exception('Input orbit for Gaia was not specified by hand and not \
                           found in the /raw_gaia directory. Fail!')
            elif source.lower()=='ce3':
                try:
                    org_eph = load_scf(sc_eph_cat+'/raw_ce3/ce3'+date_str+'j.scf')
                except IOError:
                    raise Exception('Input orbit for CE3 was not specified by hand and not \
                           found in the /raw_ce3 directory. Fail!')
        
        # now, if there were no errors, org_eph exists for sure:
        org_eph[:,6:] *= 1e3
        #org_eph[:,6:12] = [x*1e3 for x in org_eph[:,6:12]] # convert to meters
        #org_eph[:,6:12] = map(lambda x: x*1e3, org_eph[:,6:12]) # convert to meters
        t_eph_utc = (org_eph[:,3]*3600.0 + org_eph[:,4]*60.0 + org_eph[:,5])/86400.0
        # allow for 'overnighters':
        dd = mjuliandate(*org_eph[-1,0:3]) - mjuliandate(*org_eph[0,0:3])
        for cr in range(1,int(dd)+1):
            for nn in range(1,len(t_eph_utc)):
                if t_eph_utc[nn]<t_eph_utc[nn-1]:
                    t_eph_utc[nn] = t_eph_utc[nn] + 1
        # for GNSS, first point is always in the previous day, fix this:
        if source.lower()[0]=='p':
            t_eph_utc -= 1
        
        t_step = 10.0 #seconds
        #+- 5 min (600 sec) at the beginnig and the end
#        padding = 5 #minutes
        # load_eph should take care of this
        padding = 0 #minutes
        mjd_start = mjuliandate(date_t_start.year,date_t_start.month,date_t_start.day,\
                                date_t_start.hour,date_t_start.minute,date_t_start.second)
        mjd_end = mjuliandate(date_t_end.year,date_t_end.month,date_t_end.day,\
                                date_t_end.hour,date_t_end.minute,date_t_end.second)
        N_obs = round(86400.0*(mjd_end-mjd_start)/t_step) + \
                60.0*2*padding/t_step + 1.0
        N_obs = int(N_obs)
        ob = obs(['DUMMY'],'DUMMY','C')
        # for L_2 Sun-Earth Lagrange point LT=1500000e5/3e10= 5 sec:
        dt = datetime.timedelta(minutes=padding)
        ob.addScan(date_t_start - dt, t_step, nobs=N_obs)
        
        # if obss go overnight, t_end will be 23:59:59 (see load_eph), which
        # will result in last t_stamp being 00:00:00, given the 10 s step
        # shouldn't be so, as the next day's file will have the same
        # first t stamp
        # to prevent this:
        if date_t_end < ob.tstamps[-1]:
            ob.tstamps[-1] = date_t_end
            # another way: (not so good, involves extrapolation)
#            ob.tstamps = ob.tstamps[:-1]
#            N_obs -= 1
        
        #allocate memory:
        #r2000 = np.zeros((3,3,N_obs)) # allocate memory for the rotation matrix
        sc_bcrs = np.zeros((N_obs,12)) # allocate BCRS state
        sc_gcrs = np.zeros((N_obs,15)) # allocate GCRS state + acceleration
        sc_gtrs = np.zeros((N_obs,12)) # allocate GTRS state

        UT_eph = [(UT.hour + UT.minute/60.0 + UT.second/3600.0)/24.0 \
                  for UT in ob.tstamps]
        CT_eph = []
        # allow for 'overnighters':
        dd = mjuliandate(date_t_end.year,date_t_end.month,date_t_end.day)-\
             mjuliandate(date_t_start.year,date_t_start.month,date_t_start.day)
        for cr in range(1,int(dd)+1):
            for nn in range(1,len(UT_eph)):
                 if UT_eph[nn]<UT_eph[nn-1]:
                      UT_eph[nn] = UT_eph[nn] + 1

        #load eops
        ''' get the relevant eop entries from the catalogue: '''
        with open(cat_eop, 'r') as fc:
            fc_lines = fc.readlines()
        eops = np.zeros((7,7)) # +/- 3 days
        for jj in range(len(fc_lines)):
            if fc_lines[jj][0]!=' ' and fc_lines[jj][0]!='*':
                entry = [float(x) for x in fc_lines[jj].split()]
                if len(entry) > 0 and entry[3] == np.floor(mjd_start) - 3:
                    for kk in range(7):
                        eops[kk,0] = entry[3] # mjd
                        eops[kk,1] = entry[6] # UT1-UTC
                        eops[kk,2] = entry[6] - nsec(entry[3]) # UT1 - TAI
                        eops[kk,3] = entry[4] # x_p
                        eops[kk,4] = entry[5] # y_p
                        eops[kk,5] = entry[8] # dX
                        eops[kk,6] = entry[9] # dY
                        entry = [float(x) for x in fc_lines[jj+kk+1].split()]
                    break #exit loop

        ''' go-go-go: '''
        # interpolants for org_eph:
#        peez = []
##        plt.close('all')
##        plt.figure()
#        for jj in range(6,12):
#            p_regr = optimalFit(t_eph_utc, org_eph[:,jj], \
#                               min_order=3, max_order=15, fit_type='cheb')
#            p_cheb = cheb.chebfit(t_eph_utc, org_eph[:,jj], 15)
##            plt.subplot(2,3,jj-5)
##            print p_regr.predict(t_eph_utc)-cheb.chebval(t_eph_utc, p_cheb)
##            plt.plot(t_eph_utc,p_regr.predict(t_eph_utc)-org_eph[:,jj],'.')
##            plt.plot(t_eph_utc,cheb.chebval(t_eph_utc, p_cheb)-org_eph[:,jj],'.')
#            peez.append([p_cheb, p_regr])
##            peez.append(p_cheb)
        
        mjd0 = mjuliandate(ob.tstamps[0].year,\
                           ob.tstamps[0].month,ob.tstamps[0].day)
        if source.lower()[0] == 'p':
            r2000_keep = []
            
        for ii, tstamp in enumerate(ob.tstamps):
            ''' set dates: '''
            mjd = mjuliandate(tstamp.year, tstamp.month, tstamp.day)
            dd = mjd - mjd0
            UTC = (tstamp.hour*3600.0 + tstamp.minute*60.0 + tstamp.second)/86400.0
            JD = mjd + 2400000.5

            ''' compute tai & tt '''
            TAI, TT = taitime(mjd, UTC)
            
            ''' interpolate eops to t_obs '''
            UT1, eop_int = eop_iers(mjd, UTC, eops)

            ''' interp GC eph to the UTC of teh obs '''
            if source.lower()[0] == 'p': # GNSS case - initial orbit in GTRS/UTC
                for jj in range(6,9):
                    sc_gtrs[ii,jj], _ = lagint(11, t_eph_utc, org_eph[:,jj], UTC+dd)
            else: # RA/Gaia case - initial orbit in GCRS/UTC
                for jj in range(6,12):
    #                sc_gcrs[ii,jj], _ = lagintt(9, t_eph_utc, org_eph[:,jj], UTC+dd)
                    sc_gcrs[ii,jj], _ = lagint(9, t_eph_utc, org_eph[:,jj], UTC+dd)
    #                sc_gcrs[ii,jj] = peez[jj-6][1].predict(UTC+dd)
    #                sc_gcrs[ii,jj] = cheb.chebval(UTC+dd, peez[jj-6][0])
    #                sc_gcrs[ii,jj] = np.interp(UTC+dd, t_eph_utc, org_eph[:,jj])

            ''' rotation matrix IERS '''
            r2000 = ter2cel(tstamp, eop_int, mode='der')
            if source.lower()[0] == 'p': # keep it for GTRS vel calc after loop
                r2000_keep.append(r2000)

            ''' GCRS position of GNSS / GTRS position of RA/Gaia '''
            if source.lower()[0] == 'p':
                sc_gcrs[ii,6:9]  = dot(r2000[:,:,0], sc_gtrs[ii,6:9].T)
            else:
                # rotate RA/Gaia GCRS position to convert it to ITRF:
                sc_gtrs[ii,6:9]  = dot(r2000[:,:,0].T, sc_gcrs[ii,6:9].T)
                sc_gtrs[ii,9:12] = dot(r2000[:,:,0].T, sc_gcrs[ii,9:12].T) + \
                                   dot(r2000[:,:,1].T, sc_gcrs[ii,6:9].T)
                                     
            ''' compute the coordinate time fraction of the CT day at S/C position '''
            # Compute geocentric longitude of RadioAstron
            lon_gcen = atan2 ( sc_gtrs[ii,7], sc_gtrs[ii,6] )
            if ( lon_gcen < 0.0 ): lon_gcen = lon_gcen + 2.0*pi
            # Compute the stations distance from the Earth spin axis and
            # from the equatorial plane (in KM)
            u_site = sqrt(sc_gtrs[ii,6]**2 + sc_gtrs[ii,7]**2 ) *1e-3
            v_site = sc_gtrs[ii,8] *1e-3
            # subj
            CT, dTAIdCT = t_eph(JD, UT1, TT, lon_gcen, u_site, v_site)
            
            ''' calculate state vector of the Earth at JD+CT, in [m, m/s] '''
            rrd = pleph(JD+CT, 3, 12, jpl_eph)
            earth = np.reshape(np.asarray(rrd), (3,2), 'F') * 1e3

            ''' BCRS position'''
            if source.lower()[0]=='p':
                ## Sun:
                rrd = pleph(JD+CT, 11, 12, inp['jpl_eph'])
                sun = np.reshape(np.asarray(rrd), (3,2), 'F') * 1e3
                U = const.GSUN / norm(sun[:,0]-earth[:,0]);
#                sc_bcrs[ii,6:9] = dot(r2000[:,:,0], sc_gtrs[ii,6:9].T) \
#                                   + earth[:,0]
#                print sc_bcrs[ii,6:9]
                sc_bcrs[ii,6:9] = dot(r2000[:,:,0], sc_gtrs[ii,6:9].T)*\
                                    (1.0-const.L_C-U/(const.C**2)) \
                                - earth[:,1]*dot(earth[:,1], sc_gtrs[ii,6:9])/\
                                    (2.0*const.C**2) + earth[:,0]
#                print sc_bcrs[ii,6:9]
                # >> add Earth bc velocity, as sc velocities are 
                # calculated afterwards
                sc_bcrs[ii,9:12] = earth[:,1]
            else:
                # >> add Earth bc position to make SC pos barycentric
                sc_bcrs[ii,6:9] = sc_gcrs[ii,6:9] + earth[:,0]
                # >> add Earth bc velocity to make SC vel barycentric
                sc_bcrs[ii,9:12] = sc_gcrs[ii,9:12] + earth[:,1]

            ''' make CT_eph '''
            CT_eph.append(CT)


        ''' correct CT_eph if needed '''
        for cr in range(1,int(dd)+1):
            for nn in range(1,len(CT_eph)):
                 if CT_eph[nn]<CT_eph[nn-1]:
                      CT_eph[nn] = CT_eph[nn] + 1
        
        ''' calculate velosities if GNSS '''
        if source.lower()[0] == 'p':
            for jj in range(9,12):
#                p = optimalFit(UT_eph, sc_gtrs[:,jj-3], \
#                                   min_order=3, max_order=15, fit_type='cheb')
#                sc_gtrs[:,jj] = cheb.chebval(UT_eph,
#                                cheb.chebder(p.best_estimator_.coef_)) / 86400.0
                sc_gtrs[:,jj] = derivative(UT_eph, sc_gtrs[:,jj-3], \
                                            points=5, poly=2) / 86400.0
            for ii, _ in enumerate(ob.tstamps):
                sc_gcrs[ii,9:12] = dot(r2000_keep[ii][:,:,0], \
                                        sc_gtrs[ii,9:12].T) + \
                                   dot(r2000_keep[ii][:,:,1], \
                                        sc_gtrs[ii,6:9].T)
                sc_bcrs[ii,9:12] += sc_gcrs[ii,9:12]
                
        ''' calculate acceleration using an optimal Chebyshёv poly '''
        for jj in range(12,15):
            p = optimalFit(UT_eph, sc_gcrs[:,jj-3], \
                               min_order=3, max_order=15, fit_type='cheb')
            sc_gcrs[:,jj] = cheb.chebval(UT_eph,
                            cheb.chebder(p.best_estimator_.coef_)) / 86400.0
                            
        ''' resample to 1 sec grid '''
        if t_step!=1:
            N_obs = int((date_t_end - date_t_start + 2*dt).total_seconds()) + 1
            date_t_stamps = [date_t_start - dt + datetime.timedelta(seconds=x) \
                             for x in range(N_obs)]
#            print np.array(date_t_stamps)
            t0 = ((date_t_start-dt).hour*3600.0 + (date_t_start-dt).minute*60.0 + \
                         (date_t_start-dt).second)/86400.0
            t_stamps = [t0 + (x - (date_t_start-dt)).total_seconds()/86400.0 \
                        for x in date_t_stamps]
            t_stamps = np.array(t_stamps)
#            print t_stamps
            sc_gcrs_save = np.zeros((N_obs, 15))
            sc_gtrs_save = np.zeros((N_obs, 12))
            for jj in range(6,12):
#                p_regr = optimalFit(UT_eph, sc_gtrs[:,jj], \
#                                   min_order=3, max_order=15, fit_type='cheb')
#                sc_gtrs_save[:,jj] = p_regr.predict(t_stamps)
#                p_regr = optimalFit(UT_eph, sc_gcrs[:,jj], \
#                                   min_order=3, max_order=15, fit_type='cheb')
#                sc_gcrs_save[:,jj] = p_regr.predict(t_stamps)
                sc_gcrs_save[:,jj],_ = lagint(9, UT_eph, sc_gcrs[:,jj], t_stamps)
                sc_gtrs_save[:,jj],_ = lagint(9, UT_eph, sc_gtrs[:,jj], t_stamps)
            ob.resample()
        else:
            t_stamps = UT_eph
            sc_gcrs_save = sc_gcrs
            sc_gtrs_save = sc_gtrs
        
        ''' interp RA bc eph to proper time stamps '''
        sc_bcrs_save = np.zeros((N_obs, 12))
        for jj in range(6,12):
#            p_regr = optimalFit(CT_eph, sc_bcrs[:,jj], \
#                               min_order=3, max_order=15, fit_type='cheb')
#            sc_bcrs_save[:,jj] = p_regr.predict(t_stamps)
            # this works well only if UT_eph is within CT_eph...:
            sc_bcrs_save[:,jj], _ = lagint(9, CT_eph, sc_bcrs[:,jj], t_stamps)
            # do not use the smoothing (s=0)
            #tck = interpolate.splrep(CT_eph,sc_bcrs[:,jj],s=0)
            # derivative is undesired:
            #sc_bcrs_save[:,jj] = interpolate.splev(UT_eph,tck,der=0)
            # turn interpolator into extrapolator (slow on long series)
            # [best so far]
#            f = interpolate.interp1d(CT_eph, sc_bcrs[:,jj])
#            f = extrap1d(f)
#            sc_bcrs_save[:,jj] = f(UT_eph)
#            # linear extrapolator/ lagrange interpolator:
#            sc_bcrs_save[:,jj] = extrap(UT_eph, CT_eph, sc_bcrs[:,jj], \
#                                        interp_type='linear')

        #print '%20.12f'%sc_bcrs_save[-1,6]
        #plt.plot(UT_eph,sc_bcrs_save[:,6])
        #plt.show()
        
        ''' save to files '''
        # GTRS
        with open(sc_eph_cat+'/'+sc_gtrs_eph,'w') as f:
            for jj in range(int(N_obs)):
                s = '{:5d} {:02d} {:02d} {:02d} {:02d} {:06.3f} '\
                     .format(ob.tstamps[jj].year,ob.tstamps[jj].month,
                             ob.tstamps[jj].day,ob.tstamps[jj].hour,
                             ob.tstamps[jj].minute,ob.tstamps[jj].second)+\
                    '{:20.6f} {:18.6f} {:18.6f} {:16.10f} {:15.10f} {:15.10f}\n'\
                    .format(*sc_gtrs_save[jj,6:12]/1e3)
                f.write(s)
        # GCRS
        with open(sc_eph_cat+'/'+sc_gcrs_eph,'w') as f:
            for jj in range(int(N_obs)):
                s = '{:5d} {:02d} {:02d} {:02d} {:02d} {:06.3f} '\
                     .format(ob.tstamps[jj].year,ob.tstamps[jj].month,
                             ob.tstamps[jj].day,ob.tstamps[jj].hour,
                             ob.tstamps[jj].minute,ob.tstamps[jj].second)+\
                    '{:20.6f} {:18.6f} {:18.6f} {:16.10f} {:15.10f} {:15.10f} '\
                    .format(*sc_gcrs_save[jj,6:12]/1e3)+\
                    '{:16.10f} {:15.10f} {:15.10f}\n'\
                    .format(*sc_gcrs_save[jj,12:15]/1e3)
                f.write(s)
        # BCRS        
        with open(sc_eph_cat+'/'+sc_bcrs_eph,'w') as f:
            for jj in range(int(N_obs)):
                s = '{:5d} {:02d} {:02d} {:02d} {:02d} {:06.3f} '\
                     .format(ob.tstamps[jj].year,ob.tstamps[jj].month,
                             ob.tstamps[jj].day,ob.tstamps[jj].hour,
                             ob.tstamps[jj].minute,ob.tstamps[jj].second)+\
                    '{:20.6f} {:18.6f} {:18.6f} {:16.10f} {:15.10f} {:15.10f}\n'\
                    .format(*sc_bcrs_save[jj,6:12]/1e3)
                f.write(s)
        
    return (sc_bcrs_eph, sc_gtrs_eph, sc_gcrs_eph)
            
'''
#==============================================================================
#  Load eph in IPM-style scf format       
#==============================================================================
'''
def load_scf(scf_file):
    with open(scf_file,'r') as f:
        f_lines = f.readlines()
    # rewrite the source file?
#    out = open(scf_file,'w')
#    for line in f_lines:
#        out.write(line+'\n')
#    out.close()
    eph = []
    # do what the matlab function 'load' does (convert strings to floats):
    for line in f_lines:
        # cut trailing spaces
        line = line.strip()
        if len(line)>0 and (line[0]=='2' or line[0]=='1'):
            for char in (':','T','/','-'):
                # this determines time stamp accuracy, actually:
                firstSpacePos = line.index(' ')
                line = line[0:firstSpacePos].replace(char,' ') + \
                       line[firstSpacePos:]
            line = [float(x) for x in line.split()]
#            line[0:6] = map(int,line[0:6])
            line[0:5] = map(int,line[0:5])
            eph.append(line)
    #convert output to a numpy array
    return np.asarray(eph)

'''
#==============================================================================
#  Download and load eph in sp3-format for GNSS
#==============================================================================
'''
def load_sp3(sc_eph_cat, source, date_t_start, load=True):
    '''
        Download and load sp3-file for GNSS
    '''
    # set up dates and names:
    mjd = mjuliandate(date_t_start.year, date_t_start.month, date_t_start.day)
    jd = mjd + 2400000.5
    dow = np.fmod(jd + 0.5, 7) + 1
    gnss_week = np.floor( (mjd - mjuliandate(1980,1,7))/7.0 )
    if dow==7:
        dow=0
        gnss_week += 1
    
    gnss_week = int(gnss_week)
    dow = int(dow)
    
    if source[0:2].lower() == 'pr': # GLONASS
        gnss_sp3 = os.path.join(sc_eph_cat,'raw_gnss/igl{:04d}{:01d}.sp3'.\
                                format(gnss_week, dow))
        sp3_name = 'igl{:04d}{:01d}.sp3'.format(gnss_week, dow)
    elif source[0:2].lower() == 'pg': # GPS
        gnss_sp3 = os.path.join(sc_eph_cat,'raw_gnss/igs{:04d}{:01d}.sp3'.\
                                format(gnss_week, dow))
        sp3_name = 'igs{:04d}{:01d}.sp3'.format(gnss_week, dow)
    else:
        raise Exception('Wrong GNSS source name')
    
    # doesn't exist? download first then:
    if not os.path.isfile(gnss_sp3):
        print 'GNSS sp3-file: {:s}.Z not found, fetching...'.format(sp3_name)
        try:
            ftp = FTP('cddis.nasa.gov')
            ftp.login() # user anonymous, passwd anonymous
            if 'igl' in sp3_name:
                folder = 'glonass'
            elif 'igs' in sp3_name:
                folder = 'gps'
            ftp.cwd(folder+'/products/{:04d}'.format(gnss_week))
            if '{:s}.Z'.format(sp3_name) in ftp.nlst(): # check that the file is there
                ftp.retrbinary('RETR {:s}.Z'.format(sp3_name), \
                           open('{:s}.Z'.format(gnss_sp3), 'wb').write)
                # uncompress:
                print 'uncompressing: ' + '{:s}.Z'.format(sp3_name) + '...'
                os.system('uncompress -f {:s}'.format(gnss_sp3))
            else:
                raise Exception('file {:s} not found on the server. fail!'.\
                        format(sp3_name+'.Z'))
            ftp.quit()
            
        except Exception, err:
            print str(err)
            print 'Failed to download {:s} from cddis.nasa.gov'.format(sp3_name)

    # load:
    if load:
        with open(gnss_sp3) as f:
            f_lines = f.readlines()
        
        # remove comments:
        f_lines = [l for l in f_lines if l[0] not in ['#','+','%','/']]
        # extract time in GPS time scale and convert to UTC:
        utc_gps = nsec(mjd) - nsec(mjuliandate(1980,1,6))
        time = [map(float, t[1:].split()) for t in f_lines if t[0]=='*']
        time = [datetime.datetime(*map(int, t)) - \
                datetime.timedelta(seconds=utc_gps) for t in time]
        time = [[t.year, t.month, t.day, t.hour, t.minute, t.second] for t in time]
    
        xyz = [map(float, x.split()[1:4]) for x in f_lines \
                if source.lower() in x.lower()]
    
        # output
        eph = []
        for t, r in zip(time, xyz):
            eph.append(t+r)
    
        #convert output to a numpy array
        return np.asarray(eph)


'''
#==============================================================================
#     
#==============================================================================
'''
#@jit
def taitime(mjd, UTC):
    '''
    TAITIME computes the atomic (TAI) and terrestrial time (TT) from the UTC time
    Input variables:
        1.   mjd - THE MODIFIED JULIAN DATE AT 0:00 UTC OF THE DATE IN QUESTION.
        2.   UTC - THE UTC FRACTION OF THE UTC DAY. (HOUR/HOUR=SEC/SEC)
    Output variables:
        1.   TAI - THE ATOMIC TIME FRACTION OF THE ATOMIC TIME DAY. (DAYS)
        2.   TT  - TIME measured in days of TT. (DAYS)
    '''
    if UTC<0:
        # Compute the atomic time fraction of the atomic time day.
        TAI = UTC+1 + nsec(float(mjd-1)) / 86400.0
    else:
        TAI = UTC + nsec(float(mjd)) / 86400.0
    # Calculation of TT (in fraction of day)
    TT = TAI + 32.184/86400.0
    return TAI, TT
    
'''
#==============================================================================
# 
#==============================================================================
'''
def eop_iers(mjd, UTC, eops):
    '''
    This function takes a series of x, y, UT1-UTC, dX and dY values
    and interpolates them to an epoch of choice. This routine
    assumes that the values of x and y are in seconds of
    arc and that UT1-UTC is in seconds of time. At least
    one point before and one point after the epoch of the
    interpolation point are necessary in order for the
    interpolation scheme to work. 

  parameters are :
    RJD     - array of the epochs of data (given in mjd)
    X       - array of x polar motion (arcsec)
    Y       - array of y polar motion (arcsec)
    UT1     - array of UT1-UTC (sec)
    n       - number of points in arrays
    rjd_int - epoch for the interpolated value
  eop_int[0:6]:
    ut1_int - interpolated value of ut1-utc
    x_int   - interpolated value of x
    y_int   - interpolated value of y
    dX_int  - interpolated value of dX
    dY_int  - interpolated value of dY
  
  CALLED SUBROUTINE : LAGINT (Lagrange interpolation)
                         PMUT1_OCEANS (Diurnal and semidiurnal oceanic effects)
                         PM_GRAVI (Diurnal and semidiurnal lunisolar effects)
     coded by Ch. BIZOUARD (Observatoire de Paris) : November 2002
                                         Corrected : September 2007
     matlab version by D. DUEV (JIVE): November 2011
     python version by D. DUEV (JIVE): October 2013
     '''

    RJD = eops[:,0]
    X = eops[:,3]
    Y = eops[:,4]
    UT1 = eops[:,1]
    dX = eops[:,5]
    dY = eops[:,6]
#    n = len(eops[:,0])
    
    rjd_int = mjd + UTC

#    x_int = LAGINT (RJD,X,n,rjd_int)
#    y_int = LAGINT (RJD,Y,n,rjd_int)          
#    ut1_int = LAGINT (RJD,UT1,n,rjd_int)          
#    dX_int = LAGINT (RJD,dX,n,rjd_int)          
#    dY_int = LAGINT (RJD,dY,n,rjd_int)
#    print x_int, y_int, ut1_int, dX_int, dY_int

    x_int = lagint(4,RJD,X,rjd_int)[0][0]
    y_int = lagint(4,RJD,Y,rjd_int)[0][0]
    ut1_int = lagint(4,RJD,UT1,rjd_int)[0][0]
    dX_int = lagint(4,RJD,dX,rjd_int)[0][0]
    dY_int = lagint(4,RJD,dY,rjd_int)[0][0]

    # --------------
    # Oceanic effect      
    # --------------
    cor_x, cor_y, cor_ut1 = PMUT1_OCEANS (rjd_int)

    x_int = x_int + cor_x
    y_int = y_int + cor_y
    ut1_int = ut1_int + cor_ut1
    
    # Lunisolar effect 
    cor_x, cor_y = PM_GRAVI(rjd_int)

    x_int = x_int + cor_x
    y_int = y_int + cor_y

    # output
    eop_int = [ut1_int, x_int, y_int, dX_int, dY_int]
          
    # ut1 in seconds of the day
    ut1 = eop_int[0] + UTC * 86400.0
    
    return ut1, eop_int
     
## ----------------------------------------------------------------
#@numba.jit('f8(f8[:], f8[:], i8, f8)')
def LAGINT (X,Y,n,xint):
    ''' 
     This subroutine performs lagrangian interpolation
     within a set of (X,Y) pairs to give the y
     value corresponding to xint. This program uses a
     window of 4 data points to perform the interpolation.
     if the window size needs to be changed, this can be
     done by changing the indices in the do loops for
     variables m and j.

     PARAMETERS ARE :
     X     - array of values of the independent variable
     Y     - array of function values corresponding to x
     n     - number of points
     xint  - the x-value for which estimate of y is desired
     yout  - the y value returned to caller
     '''

    yout = 0.0
    for i in range(n-1):
        if ( xint >= X[i-1] and xint < X[i] ): k = i
    
    if ( k < 2 ): k = 2
    if ( k > n-2 ): k = n-2

    for m in range(k-1,k+3):
        term = Y[m-1]
        for j in range(k-1,k+3):
            if ( m != j ):
                term = term * (xint - X[j-1])/(X[m-1] - X[j-1]);
        yout = yout + term
    return yout


#@numba.jit('f8(f8)')
def PMUT1_OCEANS(rjd):
    '''
     This subroutine provides, in time domain, the diurnal/subdiurnal
     tidal effets on polar motion ("), UT1 (s) and LOD (s). The tidal terms,
     listed in the program above, have been extracted from the procedure   
     ortho_eop.f coed by Eanes in 1997.
     
     N.B.:  The fundamental lunisolar arguments are those of Simon et al.  
%
     These corrections should be added to "average"
     EOP values to get estimates of the instantaneous values.
%
     PARAMETERS ARE :
     rjd      - epoch of interest given in mjd
     cor_x    - tidal correction in x (sec. of arc)
     cor_y    - tidal correction in y (sec. of arc)
     cor_ut1  - tidal correction in UT1-UTC (sec. of time)
     cor_lod  - tidal correction in length of day (sec. of time)
%
     coded by Ch. Bizouard (2002), initially coded by McCarthy and 
     D.Gambis(1997) for the 8 prominent tidal waves.  
      
     arg(6),    % Array of the tidal arguments   
     Darg(6)    % Array of their time derivative 
    '''      
    halfpi = 1.5707963267948966
    secrad = 2.0*halfpi/(180.0*3600.0)

#  Oceanic tidal terms present in x (microas),y(microas),ut1(microseconds)       
#  narg(j,6) : Multipliers of GMST+pi and Delaunay arguments. 

    narg = np.asarray([\
    [1,-1, 0,-2,-2,-2],\
    [1,-2, 0,-2, 0,-1],\
    [1,-2, 0,-2, 0,-2],\
    [1, 0, 0,-2,-2,-1],\
    [1, 0, 0,-2,-2,-2],\
    [1,-1, 0,-2, 0,-1],\
    [1,-1, 0,-2, 0,-2],\
    [1, 1, 0,-2,-2,-1],\
    [1, 1, 0,-2,-2,-2],\
    [1, 0, 0,-2, 0, 0],\
    [1, 0, 0,-2, 0,-1],\
    [1, 0, 0,-2, 0,-2],\
    [1,-2, 0, 0, 0, 0],\
    [1, 0, 0, 0,-2, 0],\
    [1,-1, 0,-2, 2,-2],\
    [1, 1, 0,-2, 0,-1],\
    [1, 1, 0,-2, 0,-2],\
    [1,-1, 0, 0, 0, 0],\
    [1,-1, 0, 0, 0,-1],\
    [1, 1, 0, 0,-2, 0],\
    [1, 0,-1,-2, 2,-2],\
    [1, 0, 0,-2, 2,-1],\
    [1, 0, 0,-2, 2,-2],\
    [1, 0, 1,-2, 2,-2],\
    [1, 0,-1, 0, 0, 0],\
    [1, 0, 0, 0, 0, 1],\
    [1, 0, 0, 0, 0, 0],\
    [1, 0, 0, 0, 0,-1],\
    [1, 0, 0, 0, 0,-2],\
    [1, 0, 1, 0, 0, 0],\
    [1, 0, 0, 2,-2, 2],\
    [1,-1, 0, 0, 2, 0],\
    [1, 1, 0, 0, 0, 0],\
    [1, 1, 0, 0, 0,-1],\
    [1, 0, 0, 0, 2, 0],\
    [1, 2, 0, 0, 0, 0],\
    [1, 0, 0, 2, 0, 2],\
    [1, 0, 0, 2, 0, 1],\
    [1, 0, 0, 2, 0, 0],\
    [1, 1, 0, 2, 0, 2],\
    [1, 1, 0, 2, 0, 1],\
    [2,-3, 0,-2, 0,-2],\
    [2,-1, 0,-2,-2,-2],\
    [2,-2, 0,-2, 0,-2],\
    [2, 0, 0,-2,-2,-2],\
    [2, 0, 1,-2,-2,-2],\
    [2,-1,-1,-2, 0,-2],\
    [2,-1, 0,-2, 0,-1],\
    [2,-1, 0,-2, 0,-2],\
    [2,-1, 1,-2, 0,-2],\
    [2, 1, 0,-2,-2,-2],\
    [2, 1, 1,-2,-2,-2],\
    [2,-2, 0,-2, 2,-2],\
    [2, 0,-1,-2, 0,-2],\
    [2, 0, 0,-2, 0,-1],\
    [2, 0, 0,-2, 0,-2],\
    [2, 0, 1,-2, 0,-2],\
    [2,-1, 0,-2, 2,-2],\
    [2, 1, 0,-2, 0,-2],\
    [2,-1, 0, 0, 0, 0],\
    [2,-1, 0, 0, 0,-1],\
    [2, 0,-1,-2, 2,-2],\
    [2, 0, 0,-2, 2,-2],\
    [2, 0, 1,-2, 2,-2],\
    [2, 0, 0, 0, 0, 1],\
    [2, 0, 0, 0, 0, 0],\
    [2, 0, 0, 0, 0,-1],\
    [2, 0, 0, 0, 0,-2],\
    [2, 1, 0, 0, 0, 0],\
    [2, 1, 0, 0, 0,-1],\
    [2, 0, 0, 2, 0, 2]])
    
     
    data = np.asarray([\
           [ -0.05,   0.94,  -0.94,  -0.05,  0.396, -0.078],\
           [  0.06,   0.64,  -0.64,   0.06,  0.195, -0.059],\
           [  0.30,   3.42,  -3.42,   0.30,  1.034, -0.314],\
           [  0.08,   0.78,  -0.78,   0.08,  0.224, -0.073],\
           [  0.46,   4.15,  -4.15,   0.45,  1.187, -0.387],\
           [  1.19,   4.96,  -4.96,   1.19,  0.966, -0.474],\
           [  6.24,  26.31, -26.31,   6.23,  5.118, -2.499],\
           [  0.24,   0.94,  -0.94,   0.24,  0.172, -0.090],\
           [  1.28,   4.99,  -4.99,   1.28,  0.911, -0.475],\
           [ -0.28,  -0.77,   0.77,  -0.28, -0.093,  0.070],\
           [  9.22,  25.06, -25.06,   9.22,  3.025, -2.280],\
           [ 48.82, 132.91,-132.90,  48.82, 16.020,-12.069],\
           [ -0.32,  -0.86,   0.86,  -0.32, -0.103,  0.078],\
           [ -0.66,  -1.72,   1.72,  -0.66, -0.194,  0.154],\
           [ -0.42,  -0.92,   0.92,  -0.42, -0.083,  0.074],\
           [ -0.30,  -0.64,   0.64,  -0.30, -0.057,  0.050],\
           [ -1.61,  -3.46,   3.46,  -1.61, -0.308,  0.271],\
           [ -4.48,  -9.61,   9.61,  -4.48, -0.856,  0.751],\
           [ -0.90,  -1.93,   1.93,  -0.90, -0.172,  0.151],\
           [ -0.86,  -1.81,   1.81,  -0.86, -0.161,  0.137],\
           [  1.54,   3.03,  -3.03,   1.54,  0.315, -0.189],\
           [ -0.29,  -0.58,   0.58,  -0.29, -0.062,  0.035],\
           [ 26.13,  51.25, -51.25,  26.13,  5.512, -3.095],\
           [ -0.22,  -0.42,   0.42,  -0.22, -0.047,  0.025],\
           [ -0.61,  -1.20,   1.20,  -0.61, -0.134,  0.070],\
           [  1.54,   3.00,  -3.00,   1.54,  0.348, -0.171],\
           [-77.48,-151.74, 151.74, -77.48,-17.620,  8.548],\
           [-10.52, -20.56,  20.56, -10.52, -2.392,  1.159],\
           [  0.23,   0.44,  -0.44,   0.23,  0.052, -0.025],\
           [ -0.61,  -1.19,   1.19,  -0.61, -0.144,  0.065],\
           [ -1.09,  -2.11,   2.11,  -1.09, -0.267,  0.111],\
           [ -0.69,  -1.43,   1.43,  -0.69, -0.288,  0.043],\
           [ -3.46,  -7.28,   7.28,  -3.46, -1.610,  0.187],\
           [ -0.69,  -1.44,   1.44,  -0.69, -0.320,  0.037],\
           [ -0.37,  -1.06,   1.06,  -0.37, -0.407, -0.005],\
           [ -0.17,  -0.51,   0.51,  -0.17, -0.213, -0.005],\
           [ -1.10,  -3.42,   3.42,  -1.09, -1.436, -0.037],\
           [ -0.70,  -2.19,   2.19,  -0.70, -0.921, -0.023],\
           [ -0.15,  -0.46,   0.46,  -0.15, -0.193, -0.005],\
           [ -0.03,  -0.59,   0.59,  -0.03, -0.396, -0.024],\
           [ -0.02,  -0.38,   0.38,  -0.02, -0.253, -0.015],\
           [ -0.49,  -0.04,   0.63,   0.24, -0.089, -0.011],\
           [ -1.33,  -0.17,   1.53,   0.68, -0.224, -0.032],\
           [ -6.08,  -1.61,   3.13,   3.35, -0.637, -0.177],\
           [ -7.59,  -2.05,   3.44,   4.23, -0.745, -0.222],\
           [ -0.52,  -0.14,   0.22,   0.29, -0.049, -0.015],\
           [  0.47,   0.11,  -0.10,  -0.27,  0.033,  0.013],\
           [  2.12,   0.49,  -0.41,  -1.23,  0.141,  0.058],\
           [-56.87, -12.93,  11.15,  32.88, -3.795, -1.556],\
           [ -0.54,  -0.12,   0.10,   0.31, -0.035, -0.015],\
           [-11.01,  -2.40,   1.89,   6.41, -0.698, -0.298],\
           [ -0.51,  -0.11,   0.08,   0.30, -0.032, -0.014],\
           [  0.98,   0.11,  -0.11,  -0.58,  0.050,  0.022],\
           [  1.13,   0.11,  -0.13,  -0.67,  0.056,  0.025],\
           [ 12.32,   1.00,  -1.41,  -7.31,  0.605,  0.266],\
         [ -330.15, -26.96,  37.58, 195.92,-16.195, -7.140],\
           [ -1.01,  -0.07,   0.11,   0.60, -0.049, -0.021],\
           [  2.47,  -0.28,  -0.44,  -1.48,  0.111,  0.034],\
           [  9.40,  -1.44,  -1.88,  -5.65,  0.425,  0.117],\
           [ -2.35,   0.37,   0.47,   1.41, -0.106, -0.029],\
           [ -1.04,   0.17,   0.21,   0.62, -0.047, -0.013],\
           [ -8.51,   3.50,   3.29,   5.11, -0.437, -0.019],\
         [ -144.13,  63.56,  59.23,  86.56, -7.547, -0.159],\
           [  1.19,  -0.56,  -0.52,  -0.72,  0.064,  0.000],\
           [  0.49,  -0.25,  -0.23,  -0.29,  0.027, -0.001],\
           [-38.48,  19.14,  17.72,  23.11, -2.104,  0.041],\
           [-11.44,   5.75,   5.32,   6.87, -0.627,  0.015],\
           [ -1.24,   0.63,   0.58,   0.75, -0.068,  0.002],\
           [ -1.77,   1.79,   1.71,   1.04, -0.146,  0.037],\
           [ -0.77,   0.78,   0.75,   0.45, -0.064,  0.017],\
           [ -0.33,   0.62,   0.65,   0.19, -0.049,  0.018 ]])
 
    XSIN = data[:,0]
    XCOS = data[:,1]
    YSIN = data[:,2]
    YCOS = data[:,3]
    UTSIN = data[:,4]
    UTCOS = data[:,5]

    T = (rjd - 51544.5)/36525.0  # julian century

    # arguments in the following order : chi=GMST+pi,l,lp,F,D,Omega
    # et leur derivee temporelle 

    arg = np.zeros(6)
#    Darg = np.zeros(6)
    
    arg[0] = (67310.54841 +\
               (876600.0*3600.0 + 8640184.812866)*T +\
               0.093104*T**2 -\
               6.2e-6*T**3)*15.0 + 648000.0;

    arg[0] = fmod(arg[0],1296000.0)*secrad

#    Darg[0] = (876600.0*3600.0 + 8640184.812866 \
#                 + 2.0 * 0.093104 * T - 3.0 * 6.2e-6*T**2)*15.0
#    Darg[0] = Darg[0]* secrad / 36525.0   # rad/day


    arg[1] = -0.00024470*T**4 + 0.051635*T**3 + 31.8792*T**2 \
               + 1717915923.2178*T + 485868.249036
    arg[1] = fmod(arg[1],1296000.0)*secrad
      
#    Darg[1] = -4.*0.00024470*T**3 + 3.*0.051635*T**2 \
#               + 2.*31.8792*T + 1717915923.2178
#    Darg[1] = Darg[1]* secrad / 36525.0 # rad/day

    arg[2] = -0.00001149*T**4 - 0.000136*T**3 \
              -  0.5532*T**2 + 129596581.0481*T + 1287104.79305;
    arg[2] = fmod(arg[2],1296000.0)*secrad

#    Darg[2] = -4.0*0.00001149*T**3 - 3.*0.000136*T**2 \
#                -  2.0*0.5532*T + 129596581.0481
#    Darg[2] = Darg[2]* secrad / 36525.0 # rad/day
          
    arg[3] = 0.00000417*T**4 - 0.001037*T**3 - 12.7512*T**2 \
               + 1739527262.8478*T + 335779.526232
    arg[3] = fmod(arg[3],1296000.0)*secrad

#    Darg[3] = 4.0*0.00000417*T**3 - 3.0*0.001037*T**2 \
#              - 2.0 * 12.7512*T + 1739527262.8478
#    Darg[3] = Darg[3]* secrad / 36525.0   # rad/day
    
    arg[4] = -0.00003169*T**4 + 0.006593*T**3 - 6.3706*T**2 \
               + 1602961601.2090*T + 1072260.70369
    arg[4] = fmod(arg[4],1296000.0)*secrad

#    Darg[4] = -4.0*0.00003169*T**3 + 3.0*0.006593*T**2 \
#                - 2.0 * 6.3706*T + 1602961601.2090
#    Darg[4] = Darg[4]* secrad / 36525.0   # rad/day

    arg[5] = -0.00005939*T**4 + 0.007702*T**3 \
              + 7.4722*T**2 - 6962890.2665*T + 450160.398036
    arg[5] = fmod(arg[5],1296000.0)*secrad

#    Darg[5] = -4.0*0.00005939*T**3 + 3.0 * 0.007702*T**2 \
#               + 2.0 * 7.4722*T - 6962890.2665
#    Darg[5] = Darg[5]* secrad / 36525.0   # rad/day

    # CORRECTIONS
    ag = np.sum(narg[:,]*arg,1)
    ag = np.fmod(ag,4.0*halfpi)
#    dag = np.sum(narg[:,]*Darg,1)
   
    cor_x = np.sum(XCOS*np.cos(ag) + XSIN*np.sin(ag),0)
    cor_y = np.sum(YCOS*np.cos(ag) + YSIN*np.sin(ag),0)
    cor_ut1 = np.sum(UTCOS*np.cos(ag) + UTSIN*np.sin(ag),0)
#    cor_lod = np.sum((UTCOS*np.cos(ag) - UTSIN*np.sin(ag))*dag,0)

    cor_x   = cor_x * 1.0e-6   # arcseconds (")
    cor_y   = cor_y * 1.0e-6   # arcseconds (")
    cor_ut1 = cor_ut1 * 1.0e-6 # seconds (s)
#    cor_lod = cor_lod * 1.0e-6 # seconds (s)
    
#    return cor_x, cor_y, cor_ut1, cor_lod
#    return cor_x, cor_y, cor_ut1
    return np.array([cor_x, cor_y, cor_ut1])


def PM_GRAVI(rjd):
    '''
     This subroutine provides, in time domain, the diurnal
     lunisolar effet on polar motion (")
     
     N.B.:  The fundamental lunisolar arguments are those of Simon et al.  
%
     These corrections should be added to "average"
     EOP values to get estimates of the instantaneous values.
%
     PARAMETERS ARE :
     rjd      - epoch of interest given in mjd
     cor_x    - tidal correction in x (sec. of arc)
     cor_y    - tidal correction in y (sec. of arc)
%
     coded by Ch. Bizouard (2002)
    '''      
    halfpi = 1.5707963267948966
    secrad = 2.0*halfpi/(180.0*3600.0)

#  Diurnal lunisolar tidal terms present in x (microas),y(microas)      
#  narg(j,6) : Multipliers of GMST+pi and Delaunay arguments. 
	

    narg = np.asarray([\
         [1,-1, 0,-2, 0,-1],
         [1,-1, 0,-2, 0,-2],
         [1, 1, 0,-2,-2,-2],
         [1, 0, 0,-2, 0,-1],
         [1, 0, 0,-2, 0,-2],
         [1,-1, 0, 0, 0, 0],
         [1, 0, 0,-2, 2,-2],
         [1, 0, 0, 0, 0, 0],
         [1, 0, 0, 0, 0,-1],
         [1, 1, 0, 0, 0, 0]])
     
    XSIN = [ -.44, -2.31, -.44, -2.14, -11.36, .84, -4.76, 14.27, 1.93, .76 ]
     
    XCOS = [ .25, 1.32, .25, 1.23, 6.52, -.48, 2.73, -8.19, -1.11, -.43 ]
    
    YSIN = [ -.25, -1.32, -.25, -1.23, -6.52, .48, -2.73, 8.19, 1.11, .43 ]
    
    YCOS = [ -.44, -2.31, -.44, -2.14, -11.36, .84, -4.76, 14.27, 1.93, .76 ]
 
    T = (rjd - 51544.5)/36525.0  # julian century

# arguments in the following order : chi=GMST+pi,l,lp,F,D,Omega
# et leur derivee temporelle 

    arg = np.zeros(6)
    
    arg[0] = (67310.54841 + \
              (876600.0*3600.0 + 8640184.812866)*T + \
               0.093104*T**2 - \
               6.2e-6*T**3)*15.0 + 648000.0
    arg[0] = fmod(arg[0],1296000.0)*secrad
   
    arg[1] = -0.00024470*T**4 + 0.051635*T**3 + 31.8792*T**2 \
               + 1717915923.2178*T + 485868.249036
    arg[1] = fmod(arg[1],1296000)*secrad
     
    arg[2] = -0.00001149*T**4 - 0.000136*T**3 \
               -  0.5532*T**2 + 129596581.0481*T + 1287104.79305
    arg[2] = fmod(arg[2],1296000.0)*secrad
    
    arg[3] = 0.00000417*T**4 - 0.001037*T**3 - 12.7512*T**2 \
               + 1739527262.8478*T + 335779.526232
    arg[3] = fmod(arg[3],1296000.0)*secrad

    arg[4] = -0.00003169*T**4 + 0.006593*T**3 - 6.3706*T**2 \
              + 1602961601.2090*T + 1072260.70369
    arg[4] = fmod(arg[4],1296000.0)*secrad
  
    arg[5] = -0.00005939*T**4 + 0.007702*T**3 \
              + 7.4722*T**2 - 6962890.2665*T + 450160.398036
    arg[5] = fmod(arg[5],1296000.0)*secrad

    # CORRECTIONS    
    ag = np.sum(narg[:,]*arg,1)
    ag = np.fmod(ag,4.0*halfpi)
    
    cor_x = np.sum(XCOS*np.cos(ag) + XSIN*np.sin(ag),0)
    cor_y = np.sum(YCOS*np.cos(ag) + YSIN*np.sin(ag),0)
    
    cor_x   = cor_x * 1.0e-6   # arcseconds (")
    cor_y   = cor_y * 1.0e-6   # arcseconds (")
    
    return cor_x, cor_y
    

#==============================================================================
# IAU TDB-TT
#==============================================================================
def dtdb( DATE1, DATE2, UT, ELONG, U, V ):
    '''
   +
     - - - - - - - - -
      i a u _ D T D B
     - - - - - - - - -
   
     An approximation to TDB-TT, the difference between barycentric
     dynamical time and terrestrial time, for an observer on the Earth.
   
     The different time scales - proper, coordinate and realized - are
     related to each other
   
               TAI             <-  physically realized
                
             offset            <-  observed (nominally +32.184s)
                
               TT              <-  terrestrial time
                
       rate adjustment (L_G)   <-  definition of TT
                
               TCG             <-  time scale for GCRS
                
         "periodic" terms      <-  iau_DTDB is an implementation
                
       rate adjustment (L_C)   <-  function of solar-system ephemeris
                
               TCB             <-  time scale for BCRS
                
       rate adjustment (-L_B)  <-  definition of TDB
                
               TDB             <-  TCB scaled to track TT
                
         "periodic" terms      <-  -iau_DTDB is an approximation
                
               TT              <-  terrestrial time
   
     Adopted values for the various constants can be found in the IERS
     Conventions (McCarthy & Petit 2003).
   
     This routine is part of the International Astronomical Union's
     SOFA (Standards of Fundamental Astronomy) software collection.
   
     Status  support routine.
   
     Given
        DATE1,DATE2     d    date, TDB (Notes 1-3)
        UT              d    universal time (UT1, fraction of one day)
        ELONG           d    longitude (east positive, radians)
        U               d    distance from Earth spin axis (km)
        V               d    distance north of equatorial plane (km)
   
     Returned
       iau_DTDB         d    TDB-TT (seconds)
   
     Notes
   
     1) The date DATE1+DATE2 is a Julian Date, apportioned in any
        convenient way between the arguments DATE1 and DATE2.  For
        example, JD(TDB)=2450123.7 could be expressed in any of these
        ways, among others
   
               DATE1          DATE2
   
            2450123.7D0        0D0        (JD method)
             2451545D0      -1421.3D0     (J2000 method)
            2400000.5D0     50123.2D0     (MJD method)
            2450123.5D0       0.2D0       (date & time method)
   
        The JD method is the most natural and convenient to use in cases
        where the loss of several decimal digits of resolution is
        acceptable.  The J2000 method is best matched to the way the
        argument is handled internally and will deliver the optimum
        resolution.  The MJD method and the date & time methods are both
        good compromises between resolution and convenience.
   
        Although the date is, formally, barycentric dynamical time (TDB),
        the terrestrial dynamical time (TT) can be used with no practical
        effect on the accuracy of the prediction.
   
     2) TT can be regarded as a coordinate time that is realized as an
        offset of 32.184s from International Atomic Time, TAI.  TT is a
        specific linear transformation of geocentric coordinate time TCG,
        which is the time scale for the Geocentric Celestial Reference
        System, GCRS.
   
     3) TDB is a coordinate time, and is a specific linear transformation
        of barycentric coordinate time TCB, which is the time scale for
        the Barycentric Celestial Reference System, BCRS.
   
     4) The difference TCG-TCB depends on the masses and positions of the
        bodies of the solar system and the velocity of the Earth.  It is
        dominated by a rate difference, the residual being of a periodic
        character.  The latter, which is remeled by the present routine,
        comprises a main (annual) sinusoidal term of amplitude
        approximately 0.00166 seconds, plus planetary terms up to about
        20 microseconds, and lunar and diurnal terms up to 2 microseconds.
        These effects come from the changing transverse Doppler effect
        and gravitational ree-shift as the observer (on the Earth's
        surface) experiences variations in speed (with respect to the
        BCRS) and gravitational potential.
   
     5) TDB can be regarded as the same as TCB but with a rate adjustment
        to keep it close to TT, which is convenient for many applications.
        The history of successive attempts to define TDB is set out in
        Resolution 3 adopted by the IAU General Assembly in 2006, which
        defines a fixed TDB(TCB) transformation that is consistent with
        contemporary solar-system ephemerides.  Future ephemerides will
        imply slightly changed transformations between TCG and TCB, which
        could introduce a linear drift between TDB and TT;  however, any
        such drift is unlikely to exceed 1 nanosecond per century.
   
     6) The geocentric TDB-TT remel used in the present routine is that of
        Fairhead & Bretagnon (1990), in its full form.  It was originally
        supplied by Fairhead (private communications with P.T.Wallace,
        1990) as a Fortran subroutine.  The present routine contains an
        adaptation of the Fairhead code.  The numerical results are
        essentially unaffected by the changes, the differences with
        respect to the Fairhead & Bretagnon original being at the 1e-20 s
        level.
   
        The topocentric part of the remel is from Moyer (1981) and
        Murray (1983), with fundamental arguments adapted from
        Simon et al. 1994.  It is an approximation to the expression
        ( v / c ) . ( r / c ), where v is the barycentric velocity of
        the Earth, r is the geocentric position of the observer and
        c is the speed of light.
   
        By supplying zeroes for U and V, the topocentric part of the
        remel can be nullified, and the routine will return the Fairhead
        & Bretagnon result alone.
   
     7) During the interval 1950-2050, the absolute accuracy is better
        than +/- 3 nanoseconds relative to time ephemerides obtained by
        direct numerical integrations based on the JPL DE405 solar system
        ephemeris.
   
     8) It must be stressed that the present routine is merely a remel,
        and that numerical integration of solar-system ephemerides is the
        definitive method for predicting the relationship between TCG and
        TCB and hence between TT and TDB.
   
     References
   
        Fairhead, L., & Bretagnon, P., Astron.Astrophys., 229, 240-247
        (1990).
   
        IAU 2006 Resolution 3.
   
        McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
        IERS Technical Note No. 32, BKG (2004)
   
        Moyer, T.D., Cel.Mech., 23, 33 (1981).
   
        Murray, C.A., Vectorial Astrometry, Adam Hilger (1983).
   
        Seidelmann, P.K. et al., Explanatory Supplement to the
        Astronomical Almanac, Chapter 2, University Science Books (1992).
   
        Simon, J.L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
        Francou, G. & Laskar, J., Astron.Astrophys., 282, 663-683 (1994).
   
     This revision  2010 July 29
   
     SOFA release 2010-12-01
   
     Copyright (C) 2010 IAU SOFA Board.  See notes at end.
   
   -----------------------------------------------------------------------
    '''

    #  2Pi
    D2PI = 6.283185307179586476925287
    
    #  Degrees to radians
    DD2R = 1.745329251994329576923691e-2
    
    #  Reference epoch (J2000.0), JD
    DJ00 = 2451545.0
    
    #  Days per Julian millennium
    DJM = 365250.0
    
    #
    #  =====================
    #  Fairhead et al. remel
    #  =====================
    #
    #  787 sets of three coefficients.
    #
    #  Each set is amplitude (microseconds)
    #              frequency (radians per Julian millennium since J2000.0),
    #              phase (radians).
    #
    #  Sets   1-474 are the T**0 terms,
    #   "   475-679  "   "  T**1   "
    #   "   680-764  "   "  T**2   "
    #   "   765-784  "   "  T**3   "
    #   "   785-787  "   "  T**4   "  .
    #
    #FAIRHD = np.zeros((3,787))
    
    FAIRHD = np.array([\
     1656.674564e-6,    6283.075849991, 6.240054195,\
       22.417471e-6,    5753.384884897, 4.296977442,\
       13.839792e-6,   12566.151699983, 6.196904410,\
        4.770086e-6,     529.690965095, 0.444401603,\
        4.676740e-6,    6069.776754553, 4.021195093,\
        2.256707e-6,     213.299095438, 5.543113262,\
        1.694205e-6,      -3.523118349, 5.025132748,\
        1.554905e-6,   77713.771467920, 5.198467090,\
        1.276839e-6,    7860.419392439, 5.988822341,\
        1.193379e-6,    5223.693919802, 3.649823730,\
         1.115322e-6,    3930.209696220, 1.422745069,\
         0.794185e-6,   11506.769769794, 2.322313077,\
         0.447061e-6,      26.298319800, 3.615796498,\
         0.435206e-6,    -398.149003408, 4.349338347,\
         0.600309e-6,    1577.343542448, 2.678271909,\
         0.496817e-6,    6208.294251424, 5.696701824,\
         0.486306e-6,    5884.926846583, 0.520007179,\
         0.432392e-6,      74.781598567, 2.435898309,\
         0.468597e-6,    6244.942814354, 5.866398759,\
         0.375510e-6,    5507.553238667, 4.103476804,\
         0.243085e-6,    -775.522611324, 3.651837925,\
         0.173435e-6,   18849.227549974, 6.153743485,\
         0.230685e-6,    5856.477659115, 4.773852582,\
         0.203747e-6,   12036.460734888, 4.333987818,\
         0.143935e-6,    -796.298006816, 5.957517795,\
         0.159080e-6,   10977.078804699, 1.890075226,\
         0.119979e-6,      38.133035638, 4.551585768,\
         0.118971e-6,    5486.777843175, 1.914547226,\
         0.116120e-6,    1059.381930189, 0.873504123,\
         0.137927e-6,   11790.629088659, 1.135934669,\
         0.098358e-6,    2544.314419883, 0.092793886,\
         0.101868e-6,   -5573.142801634, 5.984503847,\
         0.080164e-6,     206.185548437, 2.095377709,\
         0.079645e-6,    4694.002954708, 2.949233637,\
         0.062617e-6,      20.775395492, 2.654394814,\
         0.075019e-6,    2942.463423292, 4.980931759,\
         0.064397e-6,    5746.271337896, 1.280308748,\
         0.063814e-6,    5760.498431898, 4.167901731,\
         0.048042e-6,    2146.165416475, 1.495846011,\
         0.048373e-6,     155.420399434, 2.251573730,\
         0.058844e-6,     426.598190876, 4.839650148,\
         0.046551e-6,      -0.980321068, 0.921573539,\
         0.054139e-6,   17260.154654690, 3.411091093,\
         0.042411e-6,    6275.962302991, 2.869567043,\
         0.040184e-6,      -7.113547001, 3.565975565,\
         0.036564e-6,    5088.628839767, 3.324679049,\
         0.040759e-6,   12352.852604545, 3.981496998,\
         0.036507e-6,     801.820931124, 6.248866009,\
         0.036955e-6,    3154.687084896, 5.071801441,\
         0.042732e-6,     632.783739313, 5.720622217,\
         0.042560e-6,  161000.685737473, 1.270837679,\
         0.040480e-6,   15720.838784878, 2.546610123,\
         0.028244e-6,   -6286.598968340, 5.069663519,\
         0.033477e-6,    6062.663207553, 4.144987272,\
         0.034867e-6,     522.577418094, 5.210064075,\
         0.032438e-6,    6076.890301554, 0.749317412,\
         0.030215e-6,    7084.896781115, 3.389610345,\
         0.029247e-6,  -71430.695617928, 4.183178762,\
         0.033529e-6,    9437.762934887, 2.404714239,\
         0.032423e-6,    8827.390269875, 5.541473556,\
         0.027567e-6,    6279.552731642, 5.040846034,\
         0.029862e-6,   12139.553509107, 1.770181024,\
         0.022509e-6,   10447.387839604, 1.460726241,\
         0.020937e-6,    8429.241266467, 0.652303414,\
         0.020322e-6,     419.484643875, 3.735430632,\
         0.024816e-6,   -1194.447010225, 1.087136918,\
         0.025196e-6,    1748.016413067, 2.901883301,\
         0.021691e-6,   14143.495242431, 5.952658009,\
         0.017673e-6,    6812.766815086, 3.186129845,\
         0.022567e-6,    6133.512652857, 3.307984806,\
         0.016155e-6,   10213.285546211, 1.331103168,\
         0.014751e-6,    1349.867409659, 4.308933301,\
         0.015949e-6,    -220.412642439, 4.005298270,\
         0.015974e-6,   -2352.866153772, 6.145309371,\
         0.014223e-6,   17789.845619785, 2.104551349,\
         0.017806e-6,      73.297125859, 3.475975097,\
         0.013671e-6,    -536.804512095, 5.971672571,\
         0.011942e-6,    8031.092263058, 2.053414715,\
         0.014318e-6,   16730.463689596, 3.016058075,\
         0.012462e-6,     103.092774219, 1.737438797,\
         0.010962e-6,       3.590428652, 2.196567739,\
         0.015078e-6,   19651.048481098, 3.969480770,\
         0.010396e-6,     951.718406251, 5.717799605,\
         0.011707e-6,   -4705.732307544, 2.654125618,\
         0.010453e-6,    5863.591206116, 1.913704550,\
         0.012420e-6,    4690.479836359, 4.734090399,\
         0.011847e-6,    5643.178563677, 5.489005403,\
         0.008610e-6,    3340.612426700, 3.661698944,\
         0.011622e-6,    5120.601145584, 4.863931876,\
         0.010825e-6,     553.569402842, 0.842715011,\
         0.008666e-6,    -135.065080035, 3.293406547,\
         0.009963e-6,     149.563197135, 4.870690598,\
         0.009858e-6,    6309.374169791, 1.061816410,\
         0.007959e-6,     316.391869657, 2.465042647,\
         0.010099e-6,     283.859318865, 1.942176992,\
         0.007147e-6,    -242.728603974, 3.661486981,\
         0.007505e-6,    5230.807466803, 4.920937029,\
         0.008323e-6,   11769.853693166, 1.229392026,\
         0.007490e-6,   -6256.777530192, 3.658444681,\
         0.009370e-6,  149854.400134205, 0.673880395,\
         0.007117e-6,      38.027672636, 5.294249518,\
         0.007857e-6,   12168.002696575, 0.525733528,\
         0.007019e-6,    6206.809778716, 0.837688810,\
         0.006056e-6,     955.599741609, 4.194535082,\
         0.008107e-6,   13367.972631107, 3.793235253,\
         0.006731e-6,    5650.292110678, 5.639906583,\
         0.007332e-6,      36.648562930, 0.114858677,\
         0.006366e-6,    4164.311989613, 2.262081818,\
         0.006858e-6,    5216.580372801, 0.642063318,\
         0.006919e-6,    6681.224853400, 6.018501522,\
         0.006826e-6,    7632.943259650, 3.458654112,\
         0.005308e-6,   -1592.596013633, 2.500382359,\
         0.005096e-6,   11371.704689758, 2.547107806,\
         0.004841e-6,    5333.900241022, 0.437078094,\
         0.005582e-6,    5966.683980335, 2.246174308,\
         0.006304e-6,   11926.254413669, 2.512929171,\
         0.006603e-6,   23581.258177318, 5.393136889,\
         0.005123e-6,      -1.484472708, 2.999641028,\
         0.004648e-6,    1589.072895284, 1.275847090,\
         0.005119e-6,    6438.496249426, 1.486539246,\
         0.004521e-6,    4292.330832950, 6.140635794,\
         0.005680e-6,   23013.539539587, 4.557814849,\
         0.005488e-6,      -3.455808046, 0.090675389,\
         0.004193e-6,    7234.794256242, 4.869091389,\
         0.003742e-6,    7238.675591600, 4.691976180,\
         0.004148e-6,    -110.206321219, 3.016173439,\
         0.004553e-6,   11499.656222793, 5.554998314,\
         0.004892e-6,    5436.993015240, 1.475415597,\
         0.004044e-6,    4732.030627343, 1.398784824,\
         0.004164e-6,   12491.370101415, 5.650931916,\
         0.004349e-6,   11513.883316794, 2.181745369,\
         0.003919e-6,   12528.018664345, 5.823319737,\
         0.003129e-6,    6836.645252834, 0.003844094,\
         0.004080e-6,   -7058.598461315, 3.690360123,\
         0.003270e-6,      76.266071276, 1.517189902,\
         0.002954e-6,    6283.143160294, 4.447203799,\
         0.002872e-6,      28.449187468, 1.158692983,\
         0.002881e-6,     735.876513532, 0.349250250,\
         0.003279e-6,    5849.364112115, 4.893384368,\
         0.003625e-6,    6209.778724132, 1.473760578,\
         0.003074e-6,     949.175608970, 5.185878737,\
         0.002775e-6,    9917.696874510, 1.030026325,\
         0.002646e-6,   10973.555686350, 3.918259169,\
         0.002575e-6,   25132.303399966, 6.109659023,\
         0.003500e-6,     263.083923373, 1.892100742,\
         0.002740e-6,   18319.536584880, 4.320519510,\
         0.002464e-6,     202.253395174, 4.698203059,\
         0.002409e-6,       2.542797281, 5.325009315,\
         0.003354e-6,  -90955.551694697, 1.942656623,\
         0.002296e-6,    6496.374945429, 5.061810696,\
         0.003002e-6,    6172.869528772, 2.797822767,\
         0.003202e-6,   27511.467873537, 0.531673101,\
         0.002954e-6,   -6283.008539689, 4.533471191,\
         0.002353e-6,     639.897286314, 3.734548088,\
         0.002401e-6,   16200.772724501, 2.605547070,\
         0.003053e-6,  233141.314403759, 3.029030662,\
         0.003024e-6,   83286.914269554, 2.355556099,\
         0.002863e-6,   17298.182327326, 5.240963796,\
         0.002103e-6,   -7079.373856808, 5.756641637,\
         0.002303e-6,   83996.847317911, 2.013686814,\
         0.002303e-6,   18073.704938650, 1.089100410,\
         0.002381e-6,      63.735898303, 0.759188178,\
         0.002493e-6,    6386.168624210, 0.645026535,\
         0.002366e-6,       3.932153263, 6.215885448,\
         0.002169e-6,   11015.106477335, 4.845297676,\
         0.002397e-6,    6243.458341645, 3.809290043,\
         0.002183e-6,    1162.474704408, 6.179611691,\
         0.002353e-6,    6246.427287062, 4.781719760,\
         0.002199e-6,    -245.831646229, 5.956152284,\
         0.001729e-6,    3894.181829542, 1.264976635,\
         0.001896e-6,   -3128.388765096, 4.914231596,\
         0.002085e-6,      35.164090221, 1.405158503,\
         0.002024e-6,   14712.317116458, 2.752035928,\
         0.001737e-6,    6290.189396992, 5.280820144,\
         0.002229e-6,     491.557929457, 1.571007057,\
         0.001602e-6,   14314.168113050, 4.203664806,\
         0.002186e-6,     454.909366527, 1.402101526,\
         0.001897e-6,   22483.848574493, 4.167932508,\
         0.001825e-6,   -3738.761430108, 0.545828785,\
         0.001894e-6,    1052.268383188, 5.817167450,\
         0.001421e-6,      20.355319399, 2.419886601,\
         0.001408e-6,   10984.192351700, 2.732084787,\
         0.001847e-6,   10873.986030480, 2.903477885,\
         0.001391e-6,   -8635.942003763, 0.593891500,\
         0.001388e-6,      -7.046236698, 1.166145902,\
         0.001810e-6,  -88860.057071188, 0.487355242,\
         0.001288e-6,   -1990.745017041, 3.913022880,\
         0.001297e-6,   23543.230504682, 3.063805171,\
         0.001335e-6,    -266.607041722, 3.995764039,\
         0.001376e-6,   10969.965257698, 5.152914309,\
         0.001745e-6,  244287.600007027, 3.626395673,\
         0.001649e-6,   31441.677569757, 1.952049260,\
         0.001416e-6,    9225.539273283, 4.996408389,\
         0.001238e-6,    4804.209275927, 5.503379738,\
         0.001472e-6,    4590.910180489, 4.164913291,\
         0.001169e-6,    6040.347246017, 5.841719038,\
         0.001039e-6,    5540.085789459, 2.769753519,\
         0.001004e-6,    -170.672870619, 0.755008103,\
         0.001284e-6,   10575.406682942, 5.306538209,\
         0.001278e-6,      71.812653151, 4.713486491,\
         0.001321e-6,   18209.330263660, 2.624866359,\
         0.001297e-6,   21228.392023546, 0.382603541,\
         0.000954e-6,    6282.095528923, 0.882213514,\
         0.001145e-6,    6058.731054289, 1.169483931,\
         0.000979e-6,    5547.199336460, 5.448375984,\
         0.000987e-6,   -6262.300454499, 2.656486959,\
         0.001070e-6, -154717.609887482, 1.827624012,\
         0.000991e-6,    4701.116501708, 4.387001801,\
         0.001155e-6,     -14.227094002, 3.042700750,\
         0.001176e-6,     277.034993741, 3.335519004,\
         0.000890e-6,   13916.019109642, 5.601498297,\
         0.000884e-6,   -1551.045222648, 1.088831705,\
         0.000876e-6,    5017.508371365, 3.969902609,\
         0.000806e-6,   15110.466119866, 5.142876744,\
         0.000773e-6,   -4136.910433516, 0.022067765,\
         0.001077e-6,     175.166059800, 1.844913056,\
         0.000954e-6,   -6284.056171060, 0.968480906,\
         0.000737e-6,    5326.786694021, 4.923831588,\
         0.000845e-6,    -433.711737877, 4.749245231,\
         0.000819e-6,    8662.240323563, 5.991247817,\
         0.000852e-6,     199.072001436, 2.189604979,\
         0.000723e-6,   17256.631536341, 6.068719637,\
         0.000940e-6,    6037.244203762, 6.197428148,\
         0.000885e-6,   11712.955318231, 3.280414875,\
         0.000706e-6,   12559.038152982, 2.824848947,\
         0.000732e-6,    2379.164473572, 2.501813417,\
         0.000764e-6,   -6127.655450557, 2.236346329,\
         0.000908e-6,     131.541961686, 2.521257490,\
         0.000907e-6,   35371.887265976, 3.370195967,\
         0.000673e-6,    1066.495477190, 3.876512374,\
         0.000814e-6,   17654.780539750, 4.627122566,\
         0.000630e-6,      36.027866677, 0.156368499,\
         0.000798e-6,     515.463871093, 5.151962502,\
         0.000798e-6,     148.078724426, 5.909225055,\
         0.000806e-6,     309.278322656, 6.054064447,\
         0.000607e-6,     -39.617508346, 2.839021623,\
         0.000601e-6,     412.371096874, 3.984225404,\
         0.000646e-6,   11403.676995575, 3.852959484,\
         0.000704e-6,   13521.751441591, 2.300991267,\
         0.000603e-6,  -65147.619767937, 4.140083146,\
         0.000609e-6,   10177.257679534, 0.437122327,\
         0.000631e-6,    5767.611978898, 4.026532329,\
         0.000576e-6,   11087.285125918, 4.760293101,\
         0.000674e-6,   14945.316173554, 6.270510511,\
         0.000726e-6,    5429.879468239, 6.039606892,\
         0.000710e-6,   28766.924424484, 5.672617711,\
         0.000647e-6,   11856.218651625, 3.397132627,\
         0.000678e-6,   -5481.254918868, 6.249666675,\
         0.000618e-6,   22003.914634870, 2.466427018,\
         0.000738e-6,    6134.997125565, 2.242668890,\
         0.000660e-6,     625.670192312, 5.864091907,\
         0.000694e-6,    3496.032826134, 2.668309141,\
         0.000531e-6,    6489.261398429, 1.681888780,\
         0.000611e-6, -143571.324284214, 2.424978312,\
         0.000575e-6,   12043.574281889, 4.216492400,\
         0.000553e-6,   12416.588502848, 4.772158039,\
         0.000689e-6,    4686.889407707, 6.224271088,\
         0.000495e-6,    7342.457780181, 3.817285811,\
         0.000567e-6,    3634.621024518, 1.649264690,\
         0.000515e-6,   18635.928454536, 3.945345892,\
         0.000486e-6,    -323.505416657, 4.061673868,\
         0.000662e-6,   25158.601719765, 1.794058369,\
         0.000509e-6,     846.082834751, 3.053874588,\
         0.000472e-6,  -12569.674818332, 5.112133338,\
         0.000461e-6,    6179.983075773, 0.513669325,\
         0.000641e-6,   83467.156352816, 3.210727723,\
         0.000520e-6,   10344.295065386, 2.445597761,\
         0.000493e-6,   18422.629359098, 1.676939306,\
         0.000478e-6,    1265.567478626, 5.487314569,\
         0.000472e-6,     -18.159247265, 1.999707589,\
         0.000559e-6,   11190.377900137, 5.783236356,\
         0.000494e-6,    9623.688276691, 3.022645053,\
         0.000463e-6,    5739.157790895, 1.411223013,\
         0.000432e-6,   16858.482532933, 1.179256434,\
         0.000574e-6,   72140.628666286, 1.758191830,\
         0.000484e-6,   17267.268201691, 3.290589143,\
         0.000550e-6,    4907.302050146, 0.864024298,\
         0.000399e-6,      14.977853527, 2.094441910,\
         0.000491e-6,     224.344795702, 0.878372791,\
         0.000432e-6,   20426.571092422, 6.003829241,\
         0.000481e-6,    5749.452731634, 4.309591964,\
         0.000480e-6,    5757.317038160, 1.142348571,\
         0.000485e-6,    6702.560493867, 0.210580917,\
         0.000426e-6,    6055.549660552, 4.274476529,\
         0.000480e-6,    5959.570433334, 5.031351030,\
         0.000466e-6,   12562.628581634, 4.959581597,\
         0.000520e-6,   39302.096962196, 4.788002889,\
         0.000458e-6,   12132.439962106, 1.880103788,\
         0.000470e-6,   12029.347187887, 1.405611197,\
         0.000416e-6,   -7477.522860216, 1.082356330,\
         0.000449e-6,   11609.862544012, 4.179989585,\
         0.000465e-6,   17253.041107690, 0.353496295,\
         0.000362e-6,   -4535.059436924, 1.583849576,\
         0.000383e-6,   21954.157609398, 3.747376371,\
         0.000389e-6,      17.252277143, 1.395753179,\
         0.000331e-6,   18052.929543158, 0.566790582,\
         0.000430e-6,   13517.870106233, 0.685827538,\
         0.000368e-6,   -5756.908003246, 0.731374317,\
         0.000330e-6,   10557.594160824, 3.710043680,\
         0.000332e-6,   20199.094959633, 1.652901407,\
         0.000384e-6,   11933.367960670, 5.827781531,\
         0.000387e-6,   10454.501386605, 2.541182564,\
         0.000325e-6,   15671.081759407, 2.178850542,\
         0.000318e-6,     138.517496871, 2.253253037,\
         0.000305e-6,    9388.005909415, 0.578340206,\
         0.000352e-6,    5749.861766548, 3.000297967,\
         0.000311e-6,    6915.859589305, 1.693574249,\
         0.000297e-6,   24072.921469776, 1.997249392,\
         0.000363e-6,    -640.877607382, 5.071820966,\
         0.000323e-6,   12592.450019783, 1.072262823,\
         0.000341e-6,   12146.667056108, 4.700657997,\
         0.000290e-6,    9779.108676125, 1.812320441,\
         0.000342e-6,    6132.028180148, 4.322238614,\
         0.000329e-6,    6268.848755990, 3.033827743,\
         0.000374e-6,   17996.031168222, 3.388716544,\
         0.000285e-6,    -533.214083444, 4.687313233,\
         0.000338e-6,    6065.844601290, 0.877776108,\
         0.000276e-6,      24.298513841, 0.770299429,\
         0.000336e-6,   -2388.894020449, 5.353796034,\
         0.000290e-6,    3097.883822726, 4.075291557,\
         0.000318e-6,     709.933048357, 5.941207518,\
         0.000271e-6,   13095.842665077, 3.208912203,\
         0.000331e-6,    6073.708907816, 4.007881169,\
         0.000292e-6,     742.990060533, 2.714333592,\
         0.000362e-6,   29088.811415985, 3.215977013,\
         0.000280e-6,   12359.966151546, 0.710872502,\
         0.000267e-6,   10440.274292604, 4.730108488,\
         0.000262e-6,     838.969287750, 1.327720272,\
         0.000250e-6,   16496.361396202, 0.898769761,\
         0.000325e-6,   20597.243963041, 0.180044365,\
         0.000268e-6,    6148.010769956, 5.152666276,\
         0.000284e-6,    5636.065016677, 5.655385808,\
         0.000301e-6,    6080.822454817, 2.135396205,\
         0.000294e-6,    -377.373607916, 3.708784168,\
         0.000236e-6,    2118.763860378, 1.733578756,\
         0.000234e-6,    5867.523359379, 5.575209112,\
         0.000268e-6, -226858.238553767, 0.069432392,\
         0.000265e-6,  167283.761587465, 4.369302826,\
         0.000280e-6,   28237.233459389, 5.304829118,\
         0.000292e-6,   12345.739057544, 4.096094132,\
         0.000223e-6,   19800.945956225, 3.069327406,\
         0.000301e-6,   43232.306658416, 6.205311188,\
         0.000264e-6,   18875.525869774, 1.417263408,\
         0.000304e-6,   -1823.175188677, 3.409035232,\
         0.000301e-6,     109.945688789, 0.510922054,\
         0.000260e-6,     813.550283960, 2.389438934,\
         0.000299e-6,  316428.228673312, 5.384595078,\
         0.000211e-6,    5756.566278634, 3.789392838,\
         0.000209e-6,    5750.203491159, 1.661943545,\
         0.000240e-6,   12489.885628707, 5.684549045,\
         0.000216e-6,    6303.851245484, 3.862942261,\
         0.000203e-6,    1581.959348283, 5.549853589,\
         0.000200e-6,    5642.198242609, 1.016115785,\
         0.000197e-6,     -70.849445304, 4.690702525,\
         0.000227e-6,    6287.008003254, 2.911891613,\
         0.000197e-6,     533.623118358, 1.048982898,\
         0.000205e-6,   -6279.485421340, 1.829362730,\
         0.000209e-6,  -10988.808157535, 2.636140084,\
         0.000208e-6,    -227.526189440, 4.127883842,\
         0.000191e-6,     415.552490612, 4.401165650,\
         0.000190e-6,   29296.615389579, 4.175658539,\
         0.000264e-6,   66567.485864652, 4.601102551,\
         0.000256e-6,   -3646.350377354, 0.506364778,\
         0.000188e-6,   13119.721102825, 2.032195842,\
         0.000185e-6,    -209.366942175, 4.694756586,\
         0.000198e-6,   25934.124331089, 3.832703118,\
         0.000195e-6,    4061.219215394, 3.308463427,\
         0.000234e-6,    5113.487598583, 1.716090661,\
         0.000188e-6,    1478.866574064, 5.686865780,\
         0.000222e-6,   11823.161639450, 1.942386641,\
         0.000181e-6,   10770.893256262, 1.999482059,\
         0.000171e-6,    6546.159773364, 1.182807992,\
         0.000206e-6,      70.328180442, 5.934076062,\
         0.000169e-6,   20995.392966449, 2.169080622,\
         0.000191e-6,   10660.686935042, 5.405515999,\
         0.000228e-6,   33019.021112205, 4.656985514,\
         0.000184e-6,   -4933.208440333, 3.327476868,\
         0.000220e-6,    -135.625325010, 1.765430262,\
         0.000166e-6,   23141.558382925, 3.454132746,\
         0.000191e-6,    6144.558353121, 5.020393445,\
         0.000180e-6,    6084.003848555, 0.602182191,\
         0.000163e-6,   17782.732072784, 4.960593133,\
         0.000225e-6,   16460.333529525, 2.596451817,\
         0.000222e-6,    5905.702242076, 3.731990323,\
         0.000204e-6,     227.476132789, 5.636192701,\
         0.000159e-6,   16737.577236597, 3.600691544,\
         0.000200e-6,    6805.653268085, 0.868220961,\
         0.000187e-6,   11919.140866668, 2.629456641,\
         0.000161e-6,     127.471796607, 2.862574720,\
         0.000205e-6,    6286.666278643, 1.742882331,\
         0.000189e-6,     153.778810485, 4.812372643,\
         0.000168e-6,   16723.350142595, 0.027860588,\
         0.000149e-6,   11720.068865232, 0.659721876,\
         0.000189e-6,    5237.921013804, 5.245313000,\
         0.000143e-6,    6709.674040867, 4.317625647,\
         0.000146e-6,    4487.817406270, 4.815297007,\
         0.000144e-6,    -664.756045130, 5.381366880,\
         0.000175e-6,    5127.714692584, 4.728443327,\
         0.000162e-6,    6254.626662524, 1.435132069,\
         0.000187e-6,   47162.516354635, 1.354371923,\
         0.000146e-6,   11080.171578918, 3.369695406,\
         0.000180e-6,    -348.924420448, 2.490902145,\
         0.000148e-6,     151.047669843, 3.799109588,\
         0.000157e-6,    6197.248551160, 1.284375887,\
         0.000167e-6,     146.594251718, 0.759969109,\
         0.000133e-6,   -5331.357443741, 5.409701889,\
         0.000154e-6,      95.979227218, 3.366890614,\
         0.000148e-6,   -6418.140930027, 3.384104996,\
         0.000128e-6,   -6525.804453965, 3.803419985,\
         0.000130e-6,   11293.470674356, 0.939039445,\
         0.000152e-6,   -5729.506447149, 0.734117523,\
         0.000138e-6,     210.117701700, 2.564216078,\
         0.000123e-6,    6066.595360816, 4.517099537,\
         0.000140e-6,   18451.078546566, 0.642049130,\
         0.000126e-6,   11300.584221356, 3.485280663,\
         0.000119e-6,   10027.903195729, 3.217431161,\
         0.000151e-6,    4274.518310832, 4.404359108,\
         0.000117e-6,    6072.958148291, 0.366324650,\
         0.000165e-6,   -7668.637425143, 4.298212528,\
         0.000117e-6,   -6245.048177356, 5.379518958,\
         0.000130e-6,   -5888.449964932, 4.527681115,\
         0.000121e-6,    -543.918059096, 6.109429504,\
         0.000162e-6,    9683.594581116, 5.720092446,\
         0.000141e-6,    6219.339951688, 0.679068671,\
         0.000118e-6,   22743.409379516, 4.881123092,\
         0.000129e-6,    1692.165669502, 0.351407289,\
         0.000126e-6,    5657.405657679, 5.146592349,\
         0.000114e-6,     728.762966531, 0.520791814,\
         0.000120e-6,      52.596639600, 0.948516300,\
         0.000115e-6,      65.220371012, 3.504914846,\
         0.000126e-6,    5881.403728234, 5.577502482,\
         0.000158e-6,  163096.180360983, 2.957128968,\
         0.000134e-6,   12341.806904281, 2.598576764,\
         0.000151e-6,   16627.370915377, 3.985702050,\
         0.000109e-6,    1368.660252845, 0.014730471,\
         0.000131e-6,    6211.263196841, 0.085077024,\
         0.000146e-6,    5792.741760812, 0.708426604,\
         0.000146e-6,     -77.750543984, 3.121576600,\
         0.000107e-6,    5341.013788022, 0.288231904,\
         0.000138e-6,    6281.591377283, 2.797450317,\
         0.000113e-6,   -6277.552925684, 2.788904128,\
         0.000115e-6,    -525.758811831, 5.895222200,\
         0.000138e-6,    6016.468808270, 6.096188999,\
         0.000139e-6,   23539.707386333, 2.028195445,\
         0.000146e-6,   -4176.041342449, 4.660008502,\
         0.000107e-6,   16062.184526117, 4.066520001,\
         0.000142e-6,   83783.548222473, 2.936315115,\
         0.000128e-6,    9380.959672717, 3.223844306,\
         0.000135e-6,    6205.325306007, 1.638054048,\
         0.000101e-6,    2699.734819318, 5.481603249,\
         0.000104e-6,    -568.821874027, 2.205734493,\
         0.000103e-6,    6321.103522627, 2.440421099,\
         0.000119e-6,    6321.208885629, 2.547496264,\
         0.000138e-6,    1975.492545856, 2.314608466,\
         0.000121e-6,     137.033024162, 4.539108237,\
         0.000123e-6,   19402.796952817, 4.538074405,\
         0.000119e-6,   22805.735565994, 2.869040566,\
         0.000133e-6,   64471.991241142, 6.056405489,\
         0.000129e-6,     -85.827298831, 2.540635083,\
         0.000131e-6,   13613.804277336, 4.005732868,\
         0.000104e-6,    9814.604100291, 1.959967212,\
         0.000112e-6,   16097.679950283, 3.589026260,\
         0.000123e-6,    2107.034507542, 1.728627253,\
         0.000121e-6,   36949.230808424, 6.072332087,\
         0.000108e-6,  -12539.853380183, 3.716133846,\
         0.000113e-6,   -7875.671863624, 2.725771122,\
         0.000109e-6,    4171.425536614, 4.033338079,\
         0.000101e-6,    6247.911759770, 3.441347021,\
         0.000113e-6,    7330.728427345, 0.656372122,\
         0.000113e-6,   51092.726050855, 2.791483066,\
         0.000106e-6,    5621.842923210, 1.815323326,\
         0.000101e-6,     111.430161497, 5.711033677,\
         0.000103e-6,     909.818733055, 2.812745443,\
         0.000101e-6,    1790.642637886, 1.965746028,\
#  T
       102.156724e-6,    6283.075849991, 4.249032005,\
         1.706807e-6,   12566.151699983, 4.205904248,\
         0.269668e-6,     213.299095438, 3.400290479,\
         0.265919e-6,     529.690965095, 5.836047367,\
         0.210568e-6,      -3.523118349, 6.262738348,\
         0.077996e-6,    5223.693919802, 4.670344204,\
         0.054764e-6,    1577.343542448, 4.534800170,\
         0.059146e-6,      26.298319800, 1.083044735,\
         0.034420e-6,    -398.149003408, 5.980077351,\
         0.032088e-6,   18849.227549974, 4.162913471,\
         0.033595e-6,    5507.553238667, 5.980162321,\
         0.029198e-6,    5856.477659115, 0.623811863,\
         0.027764e-6,     155.420399434, 3.745318113,\
         0.025190e-6,    5746.271337896, 2.980330535,\
         0.022997e-6,    -796.298006816, 1.174411803,\
         0.024976e-6,    5760.498431898, 2.467913690,\
         0.021774e-6,     206.185548437, 3.854787540,\
         0.017925e-6,    -775.522611324, 1.092065955,\
         0.013794e-6,     426.598190876, 2.699831988,\
         0.013276e-6,    6062.663207553, 5.845801920,\
         0.011774e-6,   12036.460734888, 2.292832062,\
         0.012869e-6,    6076.890301554, 5.333425680,\
         0.012152e-6,    1059.381930189, 6.222874454,\
         0.011081e-6,      -7.113547001, 5.154724984,\
         0.010143e-6,    4694.002954708, 4.044013795,\
         0.009357e-6,    5486.777843175, 3.416081409,\
         0.010084e-6,     522.577418094, 0.749320262,\
         0.008587e-6,   10977.078804699, 2.777152598,\
         0.008628e-6,    6275.962302991, 4.562060226,\
         0.008158e-6,    -220.412642439, 5.806891533,\
         0.007746e-6,    2544.314419883, 1.603197066,\
         0.007670e-6,    2146.165416475, 3.000200440,\
         0.007098e-6,      74.781598567, 0.443725817,\
         0.006180e-6,    -536.804512095, 1.302642751,\
         0.005818e-6,    5088.628839767, 4.827723531,\
         0.004945e-6,   -6286.598968340, 0.268305170,\
         0.004774e-6,    1349.867409659, 5.808636673,\
         0.004687e-6,    -242.728603974, 5.154890570,\
         0.006089e-6,    1748.016413067, 4.403765209,\
         0.005975e-6,   -1194.447010225, 2.583472591,\
         0.004229e-6,     951.718406251, 0.931172179,\
         0.005264e-6,     553.569402842, 2.336107252,\
         0.003049e-6,    5643.178563677, 1.362634430,\
         0.002974e-6,    6812.766815086, 1.583012668,\
         0.003403e-6,   -2352.866153772, 2.552189886,\
         0.003030e-6,     419.484643875, 5.286473844,\
         0.003210e-6,      -7.046236698, 1.863796539,\
         0.003058e-6,    9437.762934887, 4.226420633,\
         0.002589e-6,   12352.852604545, 1.991935820,\
         0.002927e-6,    5216.580372801, 2.319951253,\
         0.002425e-6,    5230.807466803, 3.084752833,\
         0.002656e-6,    3154.687084896, 2.487447866,\
         0.002445e-6,   10447.387839604, 2.347139160,\
         0.002990e-6,    4690.479836359, 6.235872050,\
         0.002890e-6,    5863.591206116, 0.095197563,\
         0.002498e-6,    6438.496249426, 2.994779800,\
         0.001889e-6,    8031.092263058, 3.569003717,\
         0.002567e-6,     801.820931124, 3.425611498,\
         0.001803e-6,  -71430.695617928, 2.192295512,\
         0.001782e-6,       3.932153263, 5.180433689,\
         0.001694e-6,   -4705.732307544, 4.641779174,\
         0.001704e-6,   -1592.596013633, 3.997097652,\
         0.001735e-6,    5849.364112115, 0.417558428,\
         0.001643e-6,    8429.241266467, 2.180619584,\
         0.001680e-6,      38.133035638, 4.164529426,\
         0.002045e-6,    7084.896781115, 0.526323854,\
         0.001458e-6,    4292.330832950, 1.356098141,\
         0.001437e-6,      20.355319399, 3.895439360,\
         0.001738e-6,    6279.552731642, 0.087484036,\
         0.001367e-6,   14143.495242431, 3.987576591,\
         0.001344e-6,    7234.794256242, 0.090454338,\
         0.001438e-6,   11499.656222793, 0.974387904,\
         0.001257e-6,    6836.645252834, 1.509069366,\
         0.001358e-6,   11513.883316794, 0.495572260,\
         0.001628e-6,    7632.943259650, 4.968445721,\
         0.001169e-6,     103.092774219, 2.838496795,\
         0.001162e-6,    4164.311989613, 3.408387778,\
         0.001092e-6,    6069.776754553, 3.617942651,\
         0.001008e-6,   17789.845619785, 0.286350174,\
         0.001008e-6,     639.897286314, 1.610762073,\
         0.000918e-6,   10213.285546211, 5.532798067,\
         0.001011e-6,   -6256.777530192, 0.661826484,\
         0.000753e-6,   16730.463689596, 3.905030235,\
         0.000737e-6,   11926.254413669, 4.641956361,\
         0.000694e-6,    3340.612426700, 2.111120332,\
         0.000701e-6,    3894.181829542, 2.760823491,\
         0.000689e-6,    -135.065080035, 4.768800780,\
         0.000700e-6,   13367.972631107, 5.760439898,\
         0.000664e-6,    6040.347246017, 1.051215840,\
         0.000654e-6,    5650.292110678, 4.911332503,\
         0.000788e-6,    6681.224853400, 4.699648011,\
         0.000628e-6,    5333.900241022, 5.024608847,\
         0.000755e-6,    -110.206321219, 4.370971253,\
         0.000628e-6,    6290.189396992, 3.660478857,\
         0.000635e-6,   25132.303399966, 4.121051532,\
         0.000534e-6,    5966.683980335, 1.173284524,\
         0.000543e-6,    -433.711737877, 0.345585464,\
         0.000517e-6,   -1990.745017041, 5.414571768,\
         0.000504e-6,    5767.611978898, 2.328281115,\
         0.000485e-6,    5753.384884897, 1.685874771,\
         0.000463e-6,    7860.419392439, 5.297703006,\
         0.000604e-6,     515.463871093, 0.591998446,\
         0.000443e-6,   12168.002696575, 4.830881244,\
         0.000570e-6,     199.072001436, 3.899190272,\
         0.000465e-6,   10969.965257698, 0.476681802,\
         0.000424e-6,   -7079.373856808, 1.112242763,\
         0.000427e-6,     735.876513532, 1.994214480,\
         0.000478e-6,   -6127.655450557, 3.778025483,\
         0.000414e-6,   10973.555686350, 5.441088327,\
         0.000512e-6,    1589.072895284, 0.107123853,\
         0.000378e-6,   10984.192351700, 0.915087231,\
         0.000402e-6,   11371.704689758, 4.107281715,\
         0.000453e-6,    9917.696874510, 1.917490952,\
         0.000395e-6,     149.563197135, 2.763124165,\
         0.000371e-6,    5739.157790895, 3.112111866,\
         0.000350e-6,   11790.629088659, 0.440639857,\
         0.000356e-6,    6133.512652857, 5.444568842,\
         0.000344e-6,     412.371096874, 5.676832684,\
         0.000383e-6,     955.599741609, 5.559734846,\
         0.000333e-6,    6496.374945429, 0.261537984,\
         0.000340e-6,    6055.549660552, 5.975534987,\
         0.000334e-6,    1066.495477190, 2.335063907,\
         0.000399e-6,   11506.769769794, 5.321230910,\
         0.000314e-6,   18319.536584880, 2.313312404,\
         0.000424e-6,    1052.268383188, 1.211961766,\
         0.000307e-6,      63.735898303, 3.169551388,\
         0.000329e-6,      29.821438149, 6.106912080,\
         0.000357e-6,    6309.374169791, 4.223760346,\
         0.000312e-6,   -3738.761430108, 2.180556645,\
         0.000301e-6,     309.278322656, 1.499984572,\
         0.000268e-6,   12043.574281889, 2.447520648,\
         0.000257e-6,   12491.370101415, 3.662331761,\
         0.000290e-6,     625.670192312, 1.272834584,\
         0.000256e-6,    5429.879468239, 1.913426912,\
         0.000339e-6,    3496.032826134, 4.165930011,\
         0.000283e-6,    3930.209696220, 4.325565754,\
         0.000241e-6,   12528.018664345, 3.832324536,\
         0.000304e-6,    4686.889407707, 1.612348468,\
         0.000259e-6,   16200.772724501, 3.470173146,\
         0.000238e-6,   12139.553509107, 1.147977842,\
         0.000236e-6,    6172.869528772, 3.776271728,\
         0.000296e-6,   -7058.598461315, 0.460368852,\
         0.000306e-6,   10575.406682942, 0.554749016,\
         0.000251e-6,   17298.182327326, 0.834332510,\
         0.000290e-6,    4732.030627343, 4.759564091,\
         0.000261e-6,    5884.926846583, 0.298259862,\
         0.000249e-6,    5547.199336460, 3.749366406,\
         0.000213e-6,   11712.955318231, 5.415666119,\
         0.000223e-6,    4701.116501708, 2.703203558,\
         0.000268e-6,    -640.877607382, 0.283670793,\
         0.000209e-6,    5636.065016677, 1.238477199,\
         0.000193e-6,   10177.257679534, 1.943251340,\
         0.000182e-6,    6283.143160294, 2.456157599,\
         0.000184e-6,    -227.526189440, 5.888038582,\
         0.000182e-6,   -6283.008539689, 0.241332086,\
         0.000228e-6,   -6284.056171060, 2.657323816,\
         0.000166e-6,    7238.675591600, 5.930629110,\
         0.000167e-6,    3097.883822726, 5.570955333,\
         0.000159e-6,    -323.505416657, 5.786670700,\
         0.000154e-6,   -4136.910433516, 1.517805532,\
         0.000176e-6,   12029.347187887, 3.139266834,\
         0.000167e-6,   12132.439962106, 3.556352289,\
         0.000153e-6,     202.253395174, 1.463313961,\
         0.000157e-6,   17267.268201691, 1.586837396,\
         0.000142e-6,   83996.847317911, 0.022670115,\
         0.000152e-6,   17260.154654690, 0.708528947,\
         0.000144e-6,    6084.003848555, 5.187075177,\
         0.000135e-6,    5756.566278634, 1.993229262,\
         0.000134e-6,    5750.203491159, 3.457197134,\
         0.000144e-6,    5326.786694021, 6.066193291,\
         0.000160e-6,   11015.106477335, 1.710431974,\
         0.000133e-6,    3634.621024518, 2.836451652,\
         0.000134e-6,   18073.704938650, 5.453106665,\
         0.000134e-6,    1162.474704408, 5.326898811,\
         0.000128e-6,    5642.198242609, 2.511652591,\
         0.000160e-6,     632.783739313, 5.628785365,\
         0.000132e-6,   13916.019109642, 0.819294053,\
         0.000122e-6,   14314.168113050, 5.677408071,\
         0.000125e-6,   12359.966151546, 5.251984735,\
         0.000121e-6,    5749.452731634, 2.210924603,\
         0.000136e-6,    -245.831646229, 1.646502367,\
         0.000120e-6,    5757.317038160, 3.240883049,\
         0.000134e-6,   12146.667056108, 3.059480037,\
         0.000137e-6,    6206.809778716, 1.867105418,\
         0.000141e-6,   17253.041107690, 2.069217456,\
         0.000129e-6,   -7477.522860216, 2.781469314,\
         0.000116e-6,    5540.085789459, 4.281176991,\
         0.000116e-6,    9779.108676125, 3.320925381,\
         0.000129e-6,    5237.921013804, 3.497704076,\
         0.000113e-6,    5959.570433334, 0.983210840,\
         0.000122e-6,    6282.095528923, 2.674938860,\
         0.000140e-6,     -11.045700264, 4.957936982,\
         0.000108e-6,   23543.230504682, 1.390113589,\
         0.000106e-6,  -12569.674818332, 0.429631317,\
         0.000110e-6,    -266.607041722, 5.501340197,\
         0.000115e-6,   12559.038152982, 4.691456618,\
         0.000134e-6,   -2388.894020449, 0.577313584,\
         0.000109e-6,   10440.274292604, 6.218148717,\
         0.000102e-6,    -543.918059096, 1.477842615,\
         0.000108e-6,   21228.392023546, 2.237753948,\
         0.000101e-6,   -4535.059436924, 3.100492232,\
         0.000103e-6,      76.266071276, 5.594294322,\
         0.000104e-6,     949.175608970, 5.674287810,\
         0.000101e-6,   13517.870106233, 2.196632348,\
         0.000100e-6,   11933.367960670, 4.056084160,\
#  T^2
         4.322990e-6,    6283.075849991, 2.642893748,\
         0.406495e-6,       0.000000000, 4.712388980,\
         0.122605e-6,   12566.151699983, 2.438140634,\
         0.019476e-6,     213.299095438, 1.642186981,\
         0.016916e-6,     529.690965095, 4.510959344,\
         0.013374e-6,      -3.523118349, 1.502210314,\
         0.008042e-6,      26.298319800, 0.478549024,\
         0.007824e-6,     155.420399434, 5.254710405,\
         0.004894e-6,    5746.271337896, 4.683210850,\
         0.004875e-6,    5760.498431898, 0.759507698,\
         0.004416e-6,    5223.693919802, 6.028853166,\
         0.004088e-6,      -7.113547001, 0.060926389,\
         0.004433e-6,   77713.771467920, 3.627734103,\
         0.003277e-6,   18849.227549974, 2.327912542,\
         0.002703e-6,    6062.663207553, 1.271941729,\
         0.003435e-6,    -775.522611324, 0.747446224,\
         0.002618e-6,    6076.890301554, 3.633715689,\
         0.003146e-6,     206.185548437, 5.647874613,\
         0.002544e-6,    1577.343542448, 6.232904270,\
         0.002218e-6,    -220.412642439, 1.309509946,\
         0.002197e-6,    5856.477659115, 2.407212349,\
         0.002897e-6,    5753.384884897, 5.863842246,\
         0.001766e-6,     426.598190876, 0.754113147,\
         0.001738e-6,    -796.298006816, 2.714942671,\
         0.001695e-6,     522.577418094, 2.629369842,\
         0.001584e-6,    5507.553238667, 1.341138229,\
         0.001503e-6,    -242.728603974, 0.377699736,\
         0.001552e-6,    -536.804512095, 2.904684667,\
         0.001370e-6,    -398.149003408, 1.265599125,\
         0.001889e-6,   -5573.142801634, 4.413514859,\
         0.001722e-6,    6069.776754553, 2.445966339,\
         0.001124e-6,    1059.381930189, 5.041799657,\
         0.001258e-6,     553.569402842, 3.849557278,\
         0.000831e-6,     951.718406251, 2.471094709,\
         0.000767e-6,    4694.002954708, 5.363125422,\
         0.000756e-6,    1349.867409659, 1.046195744,\
         0.000775e-6,     -11.045700264, 0.245548001,\
         0.000597e-6,    2146.165416475, 4.543268798,\
         0.000568e-6,    5216.580372801, 4.178853144,\
         0.000711e-6,    1748.016413067, 5.934271972,\
         0.000499e-6,   12036.460734888, 0.624434410,\
         0.000671e-6,   -1194.447010225, 4.136047594,\
         0.000488e-6,    5849.364112115, 2.209679987,\
         0.000621e-6,    6438.496249426, 4.518860804,\
         0.000495e-6,   -6286.598968340, 1.868201275,\
         0.000456e-6,    5230.807466803, 1.271231591,\
         0.000451e-6,    5088.628839767, 0.084060889,\
         0.000435e-6,    5643.178563677, 3.324456609,\
         0.000387e-6,   10977.078804699, 4.052488477,\
         0.000547e-6,  161000.685737473, 2.841633844,\
         0.000522e-6,    3154.687084896, 2.171979966,\
         0.000375e-6,    5486.777843175, 4.983027306,\
         0.000421e-6,    5863.591206116, 4.546432249,\
         0.000439e-6,    7084.896781115, 0.522967921,\
         0.000309e-6,    2544.314419883, 3.172606705,\
         0.000347e-6,    4690.479836359, 1.479586566,\
         0.000317e-6,     801.820931124, 3.553088096,\
         0.000262e-6,     419.484643875, 0.606635550,\
         0.000248e-6,    6836.645252834, 3.014082064,\
         0.000245e-6,   -1592.596013633, 5.519526220,\
         0.000225e-6,    4292.330832950, 2.877956536,\
         0.000214e-6,    7234.794256242, 1.605227587,\
         0.000205e-6,    5767.611978898, 0.625804796,\
         0.000180e-6,   10447.387839604, 3.499954526,\
         0.000229e-6,     199.072001436, 5.632304604,\
         0.000214e-6,     639.897286314, 5.960227667,\
         0.000175e-6,    -433.711737877, 2.162417992,\
         0.000209e-6,     515.463871093, 2.322150893,\
         0.000173e-6,    6040.347246017, 2.556183691,\
         0.000184e-6,    6309.374169791, 4.732296790,\
         0.000227e-6,  149854.400134205, 5.385812217,\
         0.000154e-6,    8031.092263058, 5.120720920,\
         0.000151e-6,    5739.157790895, 4.815000443,\
         0.000197e-6,    7632.943259650, 0.222827271,\
         0.000197e-6,      74.781598567, 3.910456770,\
         0.000138e-6,    6055.549660552, 1.397484253,\
         0.000149e-6,   -6127.655450557, 5.333727496,\
         0.000137e-6,    3894.181829542, 4.281749907,\
         0.000135e-6,    9437.762934887, 5.979971885,\
         0.000139e-6,   -2352.866153772, 4.715630782,\
         0.000142e-6,    6812.766815086, 0.513330157,\
         0.000120e-6,   -4705.732307544, 0.194160689,\
         0.000131e-6,  -71430.695617928, 0.000379226,\
         0.000124e-6,    6279.552731642, 2.122264908,\
         0.000108e-6,   -6256.777530192, 0.883445696,\
#  T^3
         0.143388e-6,    6283.075849991, 1.131453581,\
         0.006671e-6,   12566.151699983, 0.775148887,\
         0.001480e-6,     155.420399434, 0.480016880,\
         0.000934e-6,     213.299095438, 6.144453084,\
         0.000795e-6,     529.690965095, 2.941595619,\
         0.000673e-6,    5746.271337896, 0.120415406,\
         0.000672e-6,    5760.498431898, 5.317009738,\
         0.000389e-6,    -220.412642439, 3.090323467,\
         0.000373e-6,    6062.663207553, 3.003551964,\
         0.000360e-6,    6076.890301554, 1.918913041,\
         0.000316e-6,     -21.340641002, 5.545798121,\
         0.000315e-6,    -242.728603974, 1.884932563,\
         0.000278e-6,     206.185548437, 1.266254859,\
         0.000238e-6,    -536.804512095, 4.532664830,\
         0.000185e-6,     522.577418094, 4.578313856,\
         0.000245e-6,   18849.227549974, 0.587467082,\
         0.000180e-6,     426.598190876, 5.151178553,\
         0.000200e-6,     553.569402842, 5.355983739,\
         0.000141e-6,    5223.693919802, 1.336556009,\
         0.000104e-6,    5856.477659115, 4.239842759,\

#  T^4
         0.003826e-6,    6283.075849991, 5.705257275,\
         0.000303e-6,   12566.151699983, 5.407132842,\
         0.000209e-6,     155.420399434, 1.989815753 ])

    FAIRHD = np.reshape(FAIRHD, (3,787), order='F')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    #  Time since J2000.0 in Julian millennia.
    T = ( ( DATE1-DJ00 ) + DATE2 ) / DJM
    
    #  =================
    #  Topocentric terms
    #  =================
    
    #  Convert UT to local solar time in radians.
    TSOL = fmod(UT, 1.0) * D2PI + ELONG

    #  FUNDAMENTAL ARGUMENTS  Simon et al. 1994.

    #  Combine time argument (millennia) with deg/arcsec factor.
    W = T / 3600.0

    #  Sun Mean Longitude.
    ELSUN = fmod(280.46645683 + 1296027711.03429 * W, 360.0) * DD2R

    #  Sun Mean Anomaly.
    EMSUN = fmod(357.52910918 + 1295965810.481 * W, 360.0) * DD2R

    #  Mean Elongation of Moon from Sun.
    D = fmod(297.85019547 + 16029616012.090 * W, 360.0) * DD2R

    #  Mean Longitude of Jupiter.
    ELJ = fmod(34.35151874 + 109306899.89453 * W, 360.0) * DD2R

    #  Mean Longitude of Saturn.
    ELS = fmod(50.07744430 + 44046398.47038 * W, 360.0) * DD2R

    #  TOPOCENTRIC TERMS  Moyer 1981 and Murray 1983.
    WT =  + 0.00029e-10 * U * sin(TSOL + ELSUN - ELS) \
           + 0.00100e-10 * U * sin(TSOL - 2.0*EMSUN) \
           + 0.00133e-10 * U * sin(TSOL - D) \
           + 0.00133e-10 * U * sin(TSOL + ELSUN - ELJ) \
           - 0.00229e-10 * U * sin(TSOL + 2.0*ELSUN + EMSUN) \
           - 0.02200e-10 * V * cos(ELSUN + EMSUN) \
           + 0.05312e-10 * U * sin(TSOL - EMSUN) \
           - 0.13677e-10 * U * sin(TSOL + 2.0*ELSUN) \
           - 1.31840e-10 * V * cos(ELSUN) \
           + 3.17679e-10 * U * sin(TSOL)

    #  =====================
    #  Fairhead et al. remel
    #  =====================
    
    #  T**0
    W0 = 0.0
    for J in range(474,0,-1):
        W0 = W0 + FAIRHD[0,J-1] * sin(FAIRHD[1,J-1]*T + FAIRHD[2,J-1])
    #  T**1
    W1 = 0.0
    for J in range(679,474,-1):
        W1 = W1 + FAIRHD[0,J-1] * sin(FAIRHD[1,J-1]*T + FAIRHD[2,J-1])

    #  T**2
    W2 = 0.0
    for J in range(764,679,-1):
        W2 = W2 + FAIRHD[0,J-1] * sin(FAIRHD[1,J-1]*T + FAIRHD[2,J-1])

    #  T**3
    W3 = 0.0
    for J in range(784,764,-1):
        W3 = W3 + FAIRHD[0,J-1] * sin(FAIRHD[1,J-1]*T + FAIRHD[2,J-1])

    #  T**4
    W4 = 0.0
    for J in range(787,784,-1):
        W4 = W4 + FAIRHD[0,J-1] * sin(FAIRHD[1,J-1]*T + FAIRHD[2,J-1])
    
    #  Multiply by powers of T and combine.
    WF = T * ( T * ( T * ( T * W4 + W3 ) + W2 ) + W1 ) + W0

    #  Adjustments to use JPL planetary masses instead of IAU.
    WJ =    0.00065e-6 * sin(6069.776754*T + 4.021194) + \
            0.00033e-6 * sin( 213.299095*T + 5.543132) + \
         ( -0.00196e-6 * sin(6208.294251*T + 5.696701) ) + \
         ( -0.00173e-6 * sin(  74.781599*T + 2.435900) ) + \
            0.03638e-6 * T * T

    #  ============
    #  Final result
    #  ============
    
    #  TDB-TT in seconds.
    iau_DTDB = WT + WF + WJ

    return iau_DTDB

#==============================================================================
# 
#==============================================================================
def t_eph( JD, UT1, TT, lon_gcen, U, V ):
    ''' 
    #    SUBROUTINE T_EPH is the utility routine which computes the coordinate time
    #    "T_eph" which is the independent variable in the JPL ephemerids DE405/LE405.
    # 
    #    T_eph is approximately equal to TDB, but not exactly.
    #    Difference T_eph-TDB is ignored in SUBROUTINE T_EPH %%%
    # 
    #    In the routine difference T_eph-TT ~ TDB-TT is calculated on base
    #    of model of Fairhead  Bretagnon (1990), in its full form, where TDB is
    #    barycentric dynamical time and TT is terrestrial time.
    #    The present routine contains an adaptation of the Fairhead subroutine.
    # 
    #    The result is calculation of periodic part that has a main (annual)
    #    sinusoidal term of amplitude approximately 0.00166 seconds, plus
    #    planetary terms up to about 20 microseconds, and quasiperiodic
    #    topocentric part up to 2 microseconds and the T_eph fraction of
    #    the coordinate time day.   (DAYS)
    # 
    #    The topocentric part of the model is from Moyer (1981) and
    #    Murray (1983), with fundamental arguments adapted from
    #    Simon et al. 1994.  The topocentric part that is equal to
    #    V_E / c**2 . r', where V_E is the barycentric velocity of
    #    the Earth, r' is the geocentric position of the observer and
    #    c is the speed of light, in routine is approximated
    #    by Moyer's equation.
    # 
    #   References:
    # 
    #      1  Fairhead,L.,  Bretagnon,P., Astron.Astrophys., 229, 240-247
    #         (1990).
    # 
    #      2  Moyer,T.D., Cel.Mech., 23, 33 (1981).
    # 
    #      3  Murray,C.A., Vectorial Astrometry, Adam Hilger (1983).
    # 
    #      4  Simon J.L., Bretagnon P., Chapront J., Chapront-Touze M.,
    #         Francou G.  Laskar J., 1994, Astron.Astrophys., 282, 663-683.
    # 
    # 
    #    Input variables:
    #        1.        JD - The JULIAN DATE at zero hours UTC of the date.   (DAYS)
    #        2.     tAI - The atomic time fraction of the atomic time day. (DAYS)
    #        3.       UT1 - The UT1 time of the day.                         (Sec)
    #        4.     tT  - Time measured in days of TT.                     (DAYS)
    #        5.  lon_gcen - The geocentric east longitude of the station.    (RAD)
    #                              (0 <= lon_gcen <= 2*pi)
    #        6.        U  - Station 1 distance from the Earth spin axis.     (KM)
    #        7.        V  - Station 1 distance north of equatorial plane.    (KM)
    # 
    #    Output variables:
    #        1.        CT - The T_eph fraction of the coordinate time day.   (DAYS)
    #        2.   dTAIdCT - The partial derivative of the atomic time with
    #                        respect to the coordinate time.                 (SEC/SEC)
    # 
    #    The algorithm goes as follows:
    # 
    #   1) Calculation of MJD measured in TAI (in DAYS) (it was done in TAITIME)
    #        tAI = dfloat(mjd) + utc + dfloat(idelt)/86400.
    #      where idelt = TAI - UTC = n [sec]
    #      Value of n (=idelt) is computed by SUBROUTINE NSEC
    #      Fraction of the atomic time day at observing site #1 (=TAI)
    #      is computed by SUBROUTINE TAITIME
    # 
    #   2) Calculation of TT=TDT - Terrestrial Time (it was done in TAITIME)
    #        tT = TAI + 32.184
    # 
    #   3) Calculation of the Coordinate Time T_eph
    #       t_eph = TT + 1/c**2 * V_E * r' + P,
    #      where P are periodic term that are calculated
    #      in SUBROUTINE PERIOD_TERM
    # 
    # -----------------------------------------------------------------------
    #    SUBROUTINE T_EPH was written by V.Zharov 
    # 
    #    VERSION 1.0   21 May 2005
    # -----------------------------------------------------------------------
    # 
    '''
    
    # Sun Mean Anomaly
    f2 = [ -0.00001149,  -0.000136,  -0.5532, 129596581.0481,  1287104.793048 ]
    # Mean Elongation of Moon from Sun
    f4 = [ -0.00003169, 0.006593,  -6.3706, 1602961601.2090,  1072260.703692 ]
    # Sun Mean Longitude
    LSun = [ -0.00002353, 0.000072,   1.09158, 129602771.103429,  1009679.244588 ]
    # Mean Longitude of Jupiter
    LJup = [ -0.00001885, 0.000133,   0.80387,  10930689.989453,   123665.467464 ]
    # Mean Longitude of Saturn
    LSat = [ -0.00003500,    -0.000107,   1.86868, 4404639.847038,   180278.79948 ]
    
    sec360 = 1296000.0 # arcseconds in one turn
    
    # Calculation of periodic terms
    per_terms, dperdct = PERIOD_TERM (JD, TT)
    
    # Calculation of CT (in fraction of day)
    CT0 = TT + per_terms/86400.0
    
    # Calculation of the partial derivative of the atomic time
    # with respect to the coordinate time.  The derivatives of
    # the long period terms (without diurnal topocentric terms) are
    # computed.
    dTAIdCT = 1.0 - dperdct
    
    # Correct the CT time for the topocentric terms
    # 
    # NOTE: t in formulae of mean elements is TDB measured in
    # CENTURES of Julian years from J2000.0
    # In Simon et al. 1994 T is measured in ThOUSANDS of years.
    # Coefficients of Simon et al. 1994 were changed to agree with
    # adopted fundamental arguments
    t = ( JD + CT0 - 2451545.0 ) / 36525.0
    t2 = t * t
    t3 = t * t2
    t4 = t * t3
    
    # Convert UT1 to local solar time in radians.
    TSOL = UT1/86400.0 * (2.0*pi) + lon_gcen
    # Arguments in radians
    # Sun Mean Anomaly
    w1 = f2[0]*t4 + f2[1]*t3 + f2[2]*t2 + f2[3]*t + f2[4]
    EMSUN = fmod( w1, sec360 ) * 4.8481368110953599e-06
    # Mean Elongation of Moon from Sun
    w2 = f4[0]*t4 + f4[1]*t3 + f4[2]*t2 + f4[3]*t + f4[4]
    D = fmod( w2, sec360 ) * 4.8481368110953599e-06
    # Sun Mean Longitude
    w3 = LSun[0]*t4 + LSun[1]*t3 + LSun[2]*t2 + LSun[3]*t + LSun[4]
    ELSUN = fmod( w3, sec360 ) * 4.8481368110953599e-06
    # Mean Longitude of Jupiter
    w4 = LJup[0]*t4 + LJup[1]*t3 + LJup[2]*t2 + LJup[3]*t + LJup[4]
    ELJ = fmod( w4, sec360 ) * 4.8481368110953599e-06
    # Mean Longitude of Saturn
    w5 = LSat[0]*t4 + LSat[1]*t3 + LSat[2]*t2 + LSat[3]*t + LSat[4]
    ELS = fmod( w5, sec360 ) * 4.8481368110953599e-06
    
    # TOPOCENTRIC TERMS (in Sec):  from Moyer 1981
    CT_top = + 3.17679e-10 * U * sin(TSOL) \
             + 5.312e-12   * U * sin(TSOL - EMSUN) \
             + 1.00e-13    * U * sin(TSOL - 2.0*EMSUN) \
             - 1.3677e-11  * U * sin(TSOL + 2.0*ELSUN) \
             - 2.29e-13    * U * sin(TSOL + 2.0*ELSUN + EMSUN) \
             + 1.33e-13    * U * sin(TSOL - D) \
             - 1.3184e-10  * V * cos(ELSUN) \
             + 1.33e-13    * U * sin(TSOL + ELSUN - ELJ) \
             + 2.9e-14     * U * sin(TSOL + ELSUN - ELS) \
             - 2.20e-12    * V * cos(ELSUN + EMSUN)   # from Murray 1983.
    # Calculation of CT (in fraction of day)
    CT = CT_top/86400.0 + CT0
    return CT, dTAIdCT

#@jit
def PERIOD_TERM(jd, TT):
    '''
    Provides TDB-TT by Fairhead, Bretagnon  Lestrade
    Extended Version transmitted by A. Irwin
    Original subroutine FBL.FOR was re-written by V.Zharov.
    Derivatives of the periodic terms with respect to
    coordinate time are calculated with accuracy of
    order of 10**{-16} sec/sec (not all terms are
    included for calculation of dt3, dt5, etc)
    
    Returned value is per_terms = TDB - TT in sec.
    '''
    T   = ( (jd  - 2451545.0) + TT )/36525.0
    T_2 = T*T
    T_3 = T*T_2
    T_4 = T*T_3
    # Amplitudes of terms are given in microseconds
    # T**0
    t1 = 1656.674564 * sin(   628.3075849991 *T + 6.240054195 )   + \
       22.417471 * sin(   575.3384884897 *T + 4.296977442 )   + \
       13.839792 * sin(  1256.6151699983 *T + 6.196904410 )   + \
        4.770086 * sin(    52.9690965095 *T + 0.444401603 )   + \
        4.676740 * sin(   606.9776754553 *T + 4.021195093 )   + \
        2.256707 * sin(    21.3299095438 *T + 5.543113262 )   + \
        1.694205 * sin(    -0.3523118349 *T + 5.025132748 )   + \
        1.554905 * sin(  7771.3771467920 *T + 5.198467090 )   + \
        1.276839 * sin(   786.0419392439 *T + 5.988822341 )   + \
        1.193379 * sin(   522.3693919802 *T + 3.649823730 )   + \
        1.115322 * sin(   393.0209696220 *T + 1.422745069 )   + \
        0.794185 * sin(  1150.6769769794 *T + 2.322313077 )   + \
        0.447061 * sin(     2.6298319800 *T + 3.615796498 )   + \
        0.435206 * sin(   -39.8149003408 *T + 4.349338347 )   + \
        0.600309 * sin(   157.7343542448 *T + 2.678271909 )   + \
        0.496817 * sin(   620.8294251424 *T + 5.696701824 )   + \
        0.486306 * sin(   588.4926846583 *T + 0.520007179 )   + \
        0.432392 * sin(     7.4781598567 *T + 2.435898309 )   + \
        0.468597 * sin(   624.4942814354 *T + 5.866398759 )   + \
        0.375510 * sin(   550.7553238667 *T + 4.103476804 )

    t2 = 0.243085 * sin(   -77.5522611324 *T + 3.651837925 )   + \
        0.173435 * sin(  1884.9227549974 *T + 6.153743485 )   + \
        0.230685 * sin(   585.6477659115 *T + 4.773852582 )   + \
        0.203747 * sin(  1203.6460734888 *T + 4.333987818 )   + \
        0.143935 * sin(   -79.6298006816 *T + 5.957517795 )   + \
        0.159080 * sin(  1097.7078804699 *T + 1.890075226 )   + \
        0.119979 * sin(     3.8133035638 *T + 4.551585768 )   + \
        0.118971 * sin(   548.6777843175 *T + 1.914547226 )   + \
        0.116120 * sin(   105.9381930189 *T + 0.873504123 )   + \
        0.137927 * sin(  1179.0629088659 *T + 1.135934669 )   + \
        0.098358 * sin(   254.4314419883 *T + 0.092793886 )   + \
        0.101868 * sin(  -557.3142801634 *T + 5.984503847 )   + \
        0.080164 * sin(    20.6185548437 *T + 2.095377709 )   + \
        0.079645 * sin(   469.4002954708 *T + 2.949233637 )   + \
        0.062617 * sin(     2.0775395492 *T + 2.654394814 )   + \
        0.075019 * sin(   294.2463423292 *T + 4.980931759 )   + \
        0.064397 * sin(   574.6271337896 *T + 1.280308748 )   + \
        0.063814 * sin(   576.0498431898 *T + 4.167901731 )   + \
        0.048042 * sin(   214.6165416475 *T + 1.495846011 )   + \
        0.048373 * sin(    15.5420399434 *T + 2.251573730 )

    t3 = 0.058844 * sin(    42.6598190876 *T + 4.839650148 )   + \
        0.046551 * sin(    -0.0980321068 *T + 0.921573539 )   + \
        0.054139 * sin(  1726.0154654690 *T + 3.411091093 )   + \
        0.042411 * sin(   627.5962302991 *T + 2.869567043 )   + \
        0.040184 * sin(    -0.7113547001 *T + 3.565975565 )   + \
        0.036564 * sin(   508.8628839767 *T + 3.324679049 )   + \
        0.040759 * sin(  1235.2852604545 *T + 3.981496998 )   + \
        0.036507 * sin(    80.1820931124 *T + 6.248866009 )   + \
        0.036955 * sin(   315.4687084896 *T + 5.071801441 )   + \
        0.042732 * sin(    63.2783739313 *T + 5.720622217 )   + \
        0.042560 * sin( 16100.0685737473 *T + 1.270837679 )   + \
        0.040480 * sin(  1572.0838784878 *T + 2.546610123 )   + \
        0.028244 * sin(  -628.6598968340 *T + 5.069663519 )   + \
        0.033477 * sin(   606.2663207553 *T + 4.144987272 )   + \
        0.034867 * sin(    52.2577418094 *T + 5.210064075 )   + \
        0.032438 * sin(   607.6890301554 *T + 0.749317412 )   + \
        0.030215 * sin(   708.4896781115 *T + 3.389610345 )   + \
        0.029247 * sin( -7143.0695617928 *T + 4.183178762 )   + \
        0.033529 * sin(   943.7762934887 *T + 2.404714239 )   + \
        0.032423 * sin(   882.7390269875 *T + 5.541473556 )

    t4 = 0.027567 * sin(   627.9552731642 *T + 5.040846034 )   + \
        0.029862 * sin(  1213.9553509107 *T + 1.770181024 )   + \
        0.022509 * sin(  1044.7387839604 *T + 1.460726241 )   + \
        0.020937 * sin(   842.9241266467 *T + 0.652303414 )   + \
        0.020322 * sin(    41.9484643875 *T + 3.735430632 )   + \
        0.024816 * sin(  -119.4447010225 *T + 1.087136918 )   + \
        0.025196 * sin(   174.8016413067 *T + 2.901883301 )   + \
        0.021691 * sin(  1414.3495242431 *T + 5.952658009 )   + \
        0.017673 * sin(   681.2766815086 *T + 3.186129845 )   + \
        0.022567 * sin(   613.3512652857 *T + 3.307984806 )   + \
        0.016155 * sin(  1021.3285546211 *T + 1.331103168 )   + \
        0.014751 * sin(   134.9867409659 *T + 4.308933301 )   + \
        0.015949 * sin(   -22.0412642439 *T + 4.005298270 )   + \
        0.015974 * sin(  -235.2866153772 *T + 6.145309371 )   + \
        0.014223 * sin(  1778.9845619785 *T + 2.104551349 )   + \
        0.017806 * sin(     7.3297125859 *T + 3.475975097 )   + \
        0.013671 * sin(   -53.6804512095 *T + 5.971672571 )   + \
        0.011942 * sin(   803.1092263058 *T + 2.053414715 )   + \
        0.014318 * sin(  1673.0463689596 *T + 3.016058075 )   + \
        0.012462 * sin(    10.3092774219 *T + 1.737438797 )

    t5 = 0.010962 * sin(     0.3590428652 *T + 2.196567739 )   + \
        0.015078 * sin(  1965.1048481098 *T + 3.969480770 )   + \
        0.010396 * sin(    95.1718406251 *T + 5.717799605 )   + \
        0.011707 * sin(  -470.5732307544 *T + 2.654125618 )   + \
        0.010453 * sin(   586.3591206116 *T + 1.913704550 )   + \
        0.012420 * sin(   469.0479836359 *T + 4.734090399 )   + \
        0.011847 * sin(   564.3178563677 *T + 5.489005403 )   + \
        0.008610 * sin(   334.0612426700 *T + 3.661698944 )   + \
        0.011622 * sin(   512.0601145584 *T + 4.863931876 )   + \
        0.010825 * sin(    55.3569402842 *T + 0.842715011 )   + \
        0.008666 * sin(   -13.5065080035 *T + 3.293406547 )   + \
        0.009963 * sin(    14.9563197135 *T + 4.870690598 )   + \
        0.009858 * sin(   630.9374169791 *T + 1.061816410 )   + \
        0.007959 * sin(    31.6391869657 *T + 2.465042647 )   + \
        0.010099 * sin(    28.3859318865 *T + 1.942176992 )   + \
        0.007147 * sin(   -24.2728603974 *T + 3.661486981 )   + \
        0.007505 * sin(   523.0807466803 *T + 4.920937029 )   + \
        0.008323 * sin(  1176.9853693166 *T + 1.229392026 )   + \
        0.007490 * sin(  -625.6777530192 *T + 3.658444681 )   + \
        0.009370 * sin( 14985.4400134205 *T + 0.673880395 )

    t6 = 0.007117 * sin(     3.8027672636 *T + 5.294249518 )   + \
        0.007857 * sin(  1216.8002696575 *T + 0.525733528 )   + \
        0.007019 * sin(   620.6809778716 *T + 0.837688810 )   + \
        0.006056 * sin(    95.5599741609 *T + 4.194535082 )   + \
        0.008107 * sin(  1336.7972631107 *T + 3.793235253 )   + \
        0.006731 * sin(   565.0292110678 *T + 5.639906583 )   + \
        0.007332 * sin(     3.6648562930 *T + 0.114858677 )   + \
        0.006366 * sin(   416.4311989613 *T + 2.262081818 )   + \
        0.006858 * sin(   521.6580372801 *T + 0.642063318 )   + \
        0.006919 * sin(   668.1224853400 *T + 6.018501522 )   + \
        0.006826 * sin(   763.2943259650 *T + 3.458654112 )   + \
        0.005308 * sin(  -159.2596013633 *T + 2.500382359 )   + \
        0.005096 * sin(  1137.1704689758 *T + 2.547107806 )   + \
        0.004841 * sin(   533.3900241022 *T + 0.437078094 )   + \
        0.005582 * sin(   596.6683980335 *T + 2.246174308 )   + \
        0.006304 * sin(  1192.6254413669 *T + 2.512929171 )   + \
        0.006603 * sin(  2358.1258177318 *T + 5.393136889 )   + \
        0.005123 * sin(    -0.1484472708 *T + 2.999641028 )   + \
        0.004648 * sin(   158.9072895284 *T + 1.275847090 )   + \
        0.005119 * sin(   643.8496249426 *T + 1.486539246 )

    t7 = 0.004521 * sin(   429.2330832950 *T + 6.140635794 )   + \
        0.005680 * sin(  2301.3539539587 *T + 4.557814849 )   + \
        0.005488 * sin(    -0.3455808046 *T + 0.090675389 )   + \
        0.004193 * sin(   723.4794256242 *T + 4.869091389 )   + \
        0.003742 * sin(   723.8675591600 *T + 4.691976180 )   + \
        0.004148 * sin(   -11.0206321219 *T + 3.016173439 )   + \
        0.004553 * sin(  1149.9656222793 *T + 5.554998314 )   + \
        0.004892 * sin(   543.6993015240 *T + 1.475415597 )   + \
        0.004044 * sin(   473.2030627343 *T + 1.398784824 )   + \
        0.004164 * sin(  1249.1370101415 *T + 5.650931916 )   + \
        0.004349 * sin(  1151.3883316794 *T + 2.181745369 )   + \
        0.003919 * sin(  1252.8018664345 *T + 5.823319737 )   + \
        0.003129 * sin(   683.6645252834 *T + 0.003844094 )   + \
        0.004080 * sin(  -705.8598461315 *T + 3.690360123 )   + \
        0.003270 * sin(     7.6266071276 *T + 1.517189902 )   + \
        0.002954 * sin(   628.3143160294 *T + 4.447203799 )   + \
        0.002872 * sin(     2.8449187468 *T + 1.158692983 )   + \
        0.002881 * sin(    73.5876513532 *T + 0.349250250 )   + \
        0.003279 * sin(   584.9364112115 *T + 4.893384368 )   + \
        0.003625 * sin(   620.9778724132 *T + 1.473760578 )

    t8 = 0.003074 * sin(    94.9175608970 *T + 5.185878737 )   + \
        0.002775 * sin(   991.7696874510 *T + 1.030026325 )   + \
        0.002646 * sin(  1097.3555686350 *T + 3.918259169 )   + \
        0.002575 * sin(  2513.2303399966 *T + 6.109659023 )   + \
        0.003500 * sin(    26.3083923373 *T + 1.892100742 )   + \
        0.002740 * sin(  1831.9536584880 *T + 4.320519510 )   + \
        0.002464 * sin(    20.2253395174 *T + 4.698203059 )   + \
        0.002409 * sin(     0.2542797281 *T + 5.325009315 )   + \
        0.003354 * sin( -9095.5551694697 *T + 1.942656623 )   + \
        0.002296 * sin(   649.6374945429 *T + 5.061810696 )   + \
        0.003002 * sin(   617.2869528772 *T + 2.797822767 )   + \
        0.003202 * sin(  2751.1467873537 *T + 0.531673101 )   + \
        0.002954 * sin(  -628.3008539689 *T + 4.533471191 )   + \
        0.002353 * sin(    63.9897286314 *T + 3.734548088 )   + \
        0.002401 * sin(  1620.0772724501 *T + 2.605547070 )   + \
        0.003053 * sin( 23314.1314403759 *T + 3.029030662 )   + \
        0.003024 * sin(  8328.6914269554 *T + 2.355556099 )   + \
        0.002863 * sin(  1729.8182327326 *T + 5.240963796 )   + \
        0.002103 * sin(  -707.9373856808 *T + 5.756641637 )   + \
        0.002303 * sin(  8399.6847317911 *T + 2.013686814 )

    t9 = 0.002303 * sin(  1807.3704938650 *T + 1.089100410 )   + \
        0.002381 * sin(     6.3735898303 *T + 0.759188178 )   + \
        0.002493 * sin(   638.6168624210 *T + 0.645026535 )   + \
        0.002366 * sin(     0.3932153263 *T + 6.215885448 )   + \
        0.002169 * sin(  1101.5106477335 *T + 4.845297676 )   + \
        0.002397 * sin(   624.3458341645 *T + 3.809290043 )   + \
        0.002183 * sin(   116.2474704408 *T + 6.179611691 )   + \
        0.002353 * sin(   624.6427287062 *T + 4.781719760 )   + \
        0.002199 * sin(   -24.5831646229 *T + 5.956152284 )   + \
        0.001729 * sin(   389.4181829542 *T + 1.264976635 )   + \
        0.001896 * sin(  -312.8388765096 *T + 4.914231596 )   + \
        0.002085 * sin(     3.5164090221 *T + 1.405158503 )   + \
        0.002024 * sin(  1471.2317116458 *T + 2.752035928 )   + \
        0.001737 * sin(   629.0189396992 *T + 5.280820144 )   + \
        0.002229 * sin(    49.1557929457 *T + 1.571007057 )   + \
        0.001602 * sin(  1431.4168113050 *T + 4.203664806 )   + \
        0.002186 * sin(    45.4909366527 *T + 1.402101526 )   + \
        0.001897 * sin(  2248.3848574493 *T + 4.167932508 )   + \
        0.001825 * sin(  -373.8761430108 *T + 0.545828785 )   + \
        0.001894 * sin(   105.2268383188 *T + 5.817167450 )

    t10 = 0.001421 * sin(     2.0355319399 *T + 2.419886601 )   + \
        0.001408 * sin(  1098.4192351700 *T + 2.732084787 )   + \
        0.001847 * sin(  1087.3986030480 *T + 2.903477885 )   + \
        0.001391 * sin(  -863.5942003763 *T + 0.593891500 )   + \
        0.001388 * sin(    -0.7046236698 *T + 1.166145902 )   + \
        0.001810 * sin( -8886.0057071188 *T + 0.487355242 )   + \
        0.001288 * sin(  -199.0745017041 *T + 3.913022880 )   + \
        0.001297 * sin(  2354.3230504682 *T + 3.063805171 )   + \
        0.001335 * sin(   -26.6607041722 *T + 3.995764039 )   + \
        0.001376 * sin(  1096.9965257698 *T + 5.152914309 )   + \
        0.001745 * sin( 24428.7600007027 *T + 3.626395673 )   + \
        0.001649 * sin(  3144.1677569757 *T + 1.952049260 )   + \
        0.001416 * sin(   922.5539273283 *T + 4.996408389 )   + \
        0.001238 * sin(   480.4209275927 *T + 5.503379738 )   + \
        0.001472 * sin(   459.0910180489 *T + 4.164913291 )   + \
        0.001169 * sin(   604.0347246017 *T + 5.841719038 )   + \
        0.001039 * sin(   554.0085789459 *T + 2.769753519 )   + \
        0.001004 * sin(   -17.0672870619 *T + 0.755008103 )   + \
        0.001284 * sin(  1057.5406682942 *T + 5.306538209 )   + \
        0.001278 * sin(     7.1812653151 *T + 4.713486491 )

    t11 = 0.001321 * sin(  1820.9330263660 *T + 2.624866359 )   + \
        0.001297 * sin(  2122.8392023546 *T + 0.382603541 )   + \
        0.000954 * sin(   628.2095528923 *T + 0.882213514 )   + \
        0.001145 * sin(   605.8731054289 *T + 1.169483931 )   + \
        0.000979 * sin(   554.7199336460 *T + 5.448375984 )   + \
        0.000987 * sin(  -626.2300454499 *T + 2.656486959 )   + \
        0.001070 * sin(-15471.7609887482 *T + 1.827624012 )   + \
        0.000991 * sin(   470.1116501708 *T + 4.387001801 )   + \
        0.001155 * sin(    -1.4227094002 *T + 3.042700750 )   + \
        0.001176 * sin(    27.7034993741 *T + 3.335519004 )   + \
        0.000890 * sin(  1391.6019109642 *T + 5.601498297 )   + \
        0.000884 * sin(  -155.1045222648 *T + 1.088831705 )   + \
        0.000876 * sin(   501.7508371365 *T + 3.969902609 )   + \
        0.000806 * sin(  1511.0466119866 *T + 5.142876744 )   + \
        0.000773 * sin(  -413.6910433516 *T + 0.022067765 )   + \
        0.001077 * sin(    17.5166059800 *T + 1.844913056 )   + \
        0.000954 * sin(  -628.4056171060 *T + 0.968480906 )   + \
        0.000737 * sin(   532.6786694021 *T + 4.923831588 )   + \
        0.000845 * sin(   -43.3711737877 *T + 4.749245231 )   + \
        0.000819 * sin(   866.2240323563 *T + 5.991247817 )

    t12 = 0.000852 * sin(    19.9072001436 *T + 2.189604979 )   + \
        0.000723 * sin(  1725.6631536341 *T + 6.068719637 )   + \
        0.000940 * sin(   603.7244203762 *T + 6.197428148 )   + \
        0.000885 * sin(  1171.2955318231 *T + 3.280414875 )   + \
        0.000706 * sin(  1255.9038152982 *T + 2.824848947 )   + \
        0.000732 * sin(   237.9164473572 *T + 2.501813417 )   + \
        0.000764 * sin(  -612.7655450557 *T + 2.236346329 )   + \
        0.000908 * sin(    13.1541961686 *T + 2.521257490 )   + \
        0.000907 * sin(  3537.1887265976 *T + 3.370195967 )   + \
        0.000673 * sin(   106.6495477190 *T + 3.876512374 )   + \
        0.000814 * sin(  1765.4780539750 *T + 4.627122566 )   + \
        0.000630 * sin(     3.6027866677 *T + 0.156368499 )   + \
        0.000798 * sin(    51.5463871093 *T + 5.151962502 )   + \
        0.000798 * sin(    14.8078724426 *T + 5.909225055 )   + \
        0.000806 * sin(    30.9278322656 *T + 6.054064447 )   + \
        0.000607 * sin(    -3.9617508346 *T + 2.839021623 )   + \
        0.000601 * sin(    41.2371096874 *T + 3.984225404 )   + \
        0.000646 * sin(  1140.3676995575 *T + 3.852959484 )   + \
        0.000704 * sin(  1352.1751441591 *T + 2.300991267 )   + \
        0.000603 * sin( -6514.7619767937 *T + 4.140083146 )

    t13 = 0.000609 * sin(  1017.7257679534 *T + 0.437122327 )   + \
        0.000631 * sin(   576.7611978898 *T + 4.026532329 )   + \
        0.000576 * sin(  1108.7285125918 *T + 4.760293101 )   + \
        0.000674 * sin(  1494.5316173554 *T + 6.270510511 )   + \
        0.000726 * sin(   542.9879468239 *T + 6.039606892 )   + \
        0.000710 * sin(  2876.6924424484 *T + 5.672617711 )   + \
        0.000647 * sin(  1185.6218651625 *T + 3.397132627 )   + \
        0.000678 * sin(  -548.1254918868 *T + 6.249666675 )   + \
        0.000618 * sin(  2200.3914634870 *T + 2.466427018 )   + \
        0.000738 * sin(   613.4997125565 *T + 2.242668890 )   + \
        0.000660 * sin(    62.5670192312 *T + 5.864091907 )   + \
        0.000694 * sin(   349.6032826134 *T + 2.668309141 )   + \
        0.000531 * sin(   648.9261398429 *T + 1.681888780 )   + \
        0.000611 * sin(-14357.1324284214 *T + 2.424978312 )   + \
        0.000575 * sin(  1204.3574281889 *T + 4.216492400 )   + \
        0.000553 * sin(  1241.6588502848 *T + 4.772158039 )   + \
        0.000689 * sin(   468.6889407707 *T + 6.224271088 )   + \
        0.000495 * sin(   734.2457780181 *T + 3.817285811 )   + \
        0.000567 * sin(   363.4621024518 *T + 1.649264690 )   + \
        0.000515 * sin(  1863.5928454536 *T + 3.945345892 )

    t14 = 0.000486 * sin(   -32.3505416657 *T + 4.061673868 )   + \
        0.000662 * sin(  2515.8601719765 *T + 1.794058369 )   + \
        0.000509 * sin(    84.6082834751 *T + 3.053874588 )   + \
        0.000472 * sin( -1256.9674818332 *T + 5.112133338 )   + \
        0.000461 * sin(   617.9983075773 *T + 0.513669325 )   + \
        0.000641 * sin(  8346.7156352816 *T + 3.210727723 )   + \
        0.000520 * sin(  1034.4295065386 *T + 2.445597761 )   + \
        0.000493 * sin(  1842.2629359098 *T + 1.676939306 )   + \
        0.000478 * sin(   126.5567478626 *T + 5.487314569 )   + \
        0.000472 * sin(    -1.8159247265 *T + 1.999707589 )   + \
        0.000559 * sin(  1119.0377900137 *T + 5.783236356 )   + \
        0.000494 * sin(   962.3688276691 *T + 3.022645053 )   + \
        0.000463 * sin(   573.9157790895 *T + 1.411223013 )   + \
        0.000432 * sin(  1685.8482532933 *T + 1.179256434 )   + \
        0.000574 * sin(  7214.0628666286 *T + 1.758191830 )   + \
        0.000484 * sin(  1726.7268201691 *T + 3.290589143 )   + \
        0.000550 * sin(   490.7302050146 *T + 0.864024298 )   + \
        0.000399 * sin(     1.4977853527 *T + 2.094441910 )   + \
        0.000491 * sin(    22.4344795702 *T + 0.878372791 )   + \
        0.000432 * sin(  2042.6571092422 *T + 6.003829241 )

    t15 = 0.000481 * sin(   574.9452731634 *T + 4.309591964 )   + \
        0.000480 * sin(   575.7317038160 *T + 1.142348571 )   + \
        0.000485 * sin(   670.2560493867 *T + 0.210580917 )   + \
        0.000426 * sin(   605.5549660552 *T + 4.274476529 )   + \
        0.000480 * sin(   595.9570433334 *T + 5.031351030 )   + \
        0.000466 * sin(  1256.2628581634 *T + 4.959581597 )   + \
        0.000520 * sin(  3930.2096962196 *T + 4.788002889 )   + \
        0.000458 * sin(  1213.2439962106 *T + 1.880103788 )   + \
        0.000470 * sin(  1202.9347187887 *T + 1.405611197 )   + \
        0.000416 * sin(  -747.7522860216 *T + 1.082356330 )   + \
        0.000449 * sin(  1160.9862544012 *T + 4.179989585 )   + \
        0.000465 * sin(  1725.3041107690 *T + 0.353496295 )   + \
        0.000362 * sin(  -453.5059436924 *T + 1.583849576 )   + \
        0.000383 * sin(  2195.4157609398 *T + 3.747376371 )   + \
        0.000389 * sin(     1.7252277143 *T + 1.395753179 )   + \
        0.000331 * sin(  1805.2929543158 *T + 0.566790582 )   + \
        0.000430 * sin(  1351.7870106233 *T + 0.685827538 )   + \
        0.000368 * sin(  -575.6908003246 *T + 0.731374317 )   + \
        0.000330 * sin(  1055.7594160824 *T + 3.710043680 )   + \
        0.000332 * sin(  2019.9094959633 *T + 1.652901407 )

    t16 = 0.000384 * sin(  1193.3367960670 *T + 5.827781531 )   + \
        0.000387 * sin(  1045.4501386605 *T + 2.541182564 )   + \
        0.000325 * sin(  1567.1081759407 *T + 2.178850542 )   + \
        0.000318 * sin(    13.8517496871 *T + 2.253253037 )   + \
        0.000305 * sin(   938.8005909415 *T + 0.578340206 )   + \
        0.000352 * sin(   574.9861766548 *T + 3.000297967 )   + \
        0.000311 * sin(   691.5859589305 *T + 1.693574249 )   + \
        0.000297 * sin(  2407.2921469776 *T + 1.997249392 )   + \
        0.000363 * sin(   -64.0877607382 *T + 5.071820966 )   + \
        0.000323 * sin(  1259.2450019783 *T + 1.072262823 )   + \
        0.000341 * sin(  1214.6667056108 *T + 4.700657997 )   + \
        0.000290 * sin(   977.9108676125 *T + 1.812320441 )   + \
        0.000342 * sin(   613.2028180148 *T + 4.322238614 )   + \
        0.000329 * sin(   626.8848755990 *T + 3.033827743 )   + \
        0.000374 * sin(  1799.6031168222 *T + 3.388716544 )   + \
        0.000285 * sin(   -53.3214083444 *T + 4.687313233 )   + \
        0.000338 * sin(   606.5844601290 *T + 0.877776108 )   + \
        0.000276 * sin(     2.4298513841 *T + 0.770299429 )   + \
        0.000336 * sin(  -238.8894020449 *T + 5.353796034 )   + \
        0.000290 * sin(   309.7883822726 *T + 4.075291557 )

    t17 = 0.000318 * sin(    70.9933048357 *T + 5.941207518 )   + \
        0.000271 * sin(  1309.5842665077 *T + 3.208912203 )   + \
        0.000331 * sin(   607.3708907816 *T + 4.007881169 )   + \
        0.000292 * sin(    74.2990060533 *T + 2.714333592 )   + \
        0.000362 * sin(  2908.8811415985 *T + 3.215977013 )   + \
        0.000280 * sin(  1235.9966151546 *T + 0.710872502 )   + \
        0.000267 * sin(  1044.0274292604 *T + 4.730108488 )   + \
        0.000262 * sin(    83.8969287750 *T + 1.327720272 )   + \
        0.000250 * sin(  1649.6361396202 *T + 0.898769761 )   + \
        0.000325 * sin(  2059.7243963041 *T + 0.180044365 )   + \
        0.000268 * sin(   614.8010769956 *T + 5.152666276 )   + \
        0.000284 * sin(   563.6065016677 *T + 5.655385808 )   + \
        0.000301 * sin(   608.0822454817 *T + 2.135396205 )   + \
        0.000294 * sin(   -37.7373607916 *T + 3.708784168 )   + \
        0.000236 * sin(   211.8763860378 *T + 1.733578756 )   + \
        0.000234 * sin(   586.7523359379 *T + 5.575209112 )   + \
        0.000268 * sin(-22685.8238553767 *T + 0.069432392 )   + \
        0.000265 * sin( 16728.3761587465 *T + 4.369302826 )   + \
        0.000280 * sin(  2823.7233459389 *T + 5.304829118 )   + \
        0.000292 * sin(  1234.5739057544 *T + 4.096094132 )

    t18 =  0.000223 * sin(  1980.0945956225 *T + 3.069327406 )   + \
        0.000301 * sin(  4323.2306658416 *T + 6.205311188 )   + \
        0.000264 * sin(  1887.5525869774 *T + 1.417263408 )   + \
        0.000304 * sin(  -182.3175188677 *T + 3.409035232 )   + \
        0.000301 * sin(    10.9945688789 *T + 0.510922054 )   + \
        0.000260 * sin(    81.3550283960 *T + 2.389438934 )   + \
        0.000299 * sin( 31642.8228673312 *T + 5.384595078 )   + \
        0.000211 * sin(   575.6566278634 *T + 3.789392838 )   + \
        0.000209 * sin(   575.0203491159 *T + 1.661943545 )   + \
        0.000240 * sin(  1248.9885628707 *T + 5.684549045 )   + \
        0.000216 * sin(   630.3851245484 *T + 3.862942261 )   + \
        0.000203 * sin(   158.1959348283 *T + 5.549853589 )   + \
        0.000200 * sin(   564.2198242609 *T + 1.016115785 )   + \
        0.000197 * sin(    -7.0849445304 *T + 4.690702525 )   + \
        0.000227 * sin(   628.7008003254 *T + 2.911891613 )   + \
        0.000197 * sin(    53.3623118358 *T + 1.048982898 )   + \
        0.000205 * sin(  -627.9485421340 *T + 1.829362730 )   + \
        0.000209 * sin( -1098.8808157535 *T + 2.636140084 )   + \
        0.000208 * sin(   -22.7526189440 *T + 4.127883842 )   + \
        0.000191 * sin(    41.5552490612 *T + 4.401165650 )

    t19 =  0.000190 * sin(  2929.6615389579 *T + 4.175658539 )   + \
        0.000264 * sin(  6656.7485864652 *T + 4.601102551 )   + \
        0.000256 * sin(  -364.6350377354 *T + 0.506364778 )   + \
        0.000188 * sin(  1311.9721102825 *T + 2.032195842 )   + \
        0.000185 * sin(   -20.9366942175 *T + 4.694756586 )   + \
        0.000198 * sin(  2593.4124331089 *T + 3.832703118 )   + \
        0.000195 * sin(   406.1219215394 *T + 3.308463427 )   + \
        0.000234 * sin(   511.3487598583 *T + 1.716090661 )   + \
        0.000188 * sin(   147.8866574064 *T + 5.686865780 )   + \
        0.000222 * sin(  1182.3161639450 *T + 1.942386641 )   + \
        0.000181 * sin(  1077.0893256262 *T + 1.999482059 )   + \
        0.000171 * sin(   654.6159773364 *T + 1.182807992 )   + \
        0.000206 * sin(     7.0328180442 *T + 5.934076062 )   + \
        0.000169 * sin(  2099.5392966449 *T + 2.169080622 )   + \
        0.000191 * sin(  1066.0686935042 *T + 5.405515999 )   + \
        0.000228 * sin(  3301.9021112205 *T + 4.656985514 )   + \
        0.000184 * sin(  -493.3208440333 *T + 3.327476868 )   + \
        0.000220 * sin(   -13.5625325010 *T + 1.765430262 )   + \
        0.000166 * sin(  2314.1558382925 *T + 3.454132746 )   + \
        0.000191 * sin(   614.4558353121 *T + 5.020393445 )

    t20 =  0.000180 * sin(   608.4003848555 *T + 0.602182191 )   + \
        0.000163 * sin(  1778.2732072784 *T + 4.960593133 )   + \
        0.000225 * sin(  1646.0333529525 *T + 2.596451817 )   + \
        0.000222 * sin(   590.5702242076 *T + 3.731990323 )   + \
        0.000204 * sin(    22.7476132789 *T + 5.636192701 )   + \
        0.000159 * sin(  1673.7577236597 *T + 3.600691544 )   + \
        0.000200 * sin(   680.5653268085 *T + 0.868220961 )   + \
        0.000187 * sin(  1191.9140866668 *T + 2.629456641 )   + \
        0.000161 * sin(    12.7471796607 *T + 2.862574720 )   + \
        0.000205 * sin(   628.6666278643 *T + 1.742882331 )   + \
        0.000189 * sin(    15.3778810485 *T + 4.812372643 )   + \
        0.000168 * sin(  1672.3350142595 *T + 0.027860588 )   + \
        0.000149 * sin(  1172.0068865232 *T + 0.659721876 )   + \
        0.000189 * sin(   523.7921013804 *T + 5.245313000 )   + \
        0.000143 * sin(   670.9674040867 *T + 4.317625647 )   + \
        0.000146 * sin(   448.7817406270 *T + 4.815297007 )   + \
        0.000144 * sin(   -66.4756045130 *T + 5.381366880 )   + \
        0.000175 * sin(   512.7714692584 *T + 4.728443327 )   + \
        0.000162 * sin(   625.4626662524 *T + 1.435132069 )   + \
        0.000187 * sin(  4716.2516354635 *T + 1.354371923 )

    t21 = 0.000146 * sin(  1108.0171578918 *T + 3.369695406 )   + \
        0.000180 * sin(   -34.8924420448 *T + 2.490902145 )   + \
        0.000148 * sin(    15.1047669843 *T + 3.799109588 )   + \
        0.000157 * sin(   619.7248551160 *T + 1.284375887 )   + \
        0.000167 * sin(    14.6594251718 *T + 0.759969109 )   + \
        0.000133 * sin(  -533.1357443741 *T + 5.409701889 )   + \
        0.000154 * sin(     9.5979227218 *T + 3.366890614 )   + \
        0.000148 * sin(  -641.8140930027 *T + 3.384104996 )   + \
        0.000128 * sin(  -652.5804453965 *T + 3.803419985 )   + \
        0.000130 * sin(  1129.3470674356 *T + 0.939039445 )   + \
        0.000152 * sin(  -572.9506447149 *T + 0.734117523 )   + \
        0.000138 * sin(    21.0117701700 *T + 2.564216078 )   + \
        0.000123 * sin(   606.6595360816 *T + 4.517099537 )   + \
        0.000140 * sin(  1845.1078546566 *T + 0.642049130 )   + \
        0.000126 * sin(  1130.0584221356 *T + 3.485280663 )   + \
        0.000119 * sin(  1002.7903195729 *T + 3.217431161 )   + \
        0.000151 * sin(   427.4518310832 *T + 4.404359108 )   + \
        0.000117 * sin(   607.2958148291 *T + 0.366324650 )   + \
        0.000165 * sin(  -766.8637425143 *T + 4.298212528 )   + \
        0.000117 * sin(  -624.5048177356 *T + 5.379518958 )

    t22 =  0.000130 * sin(  -588.8449964932 *T + 4.527681115 )   + \
        0.000121 * sin(   -54.3918059096 *T + 6.109429504 )   + \
        0.000162 * sin(   968.3594581116 *T + 5.720092446 )   + \
        0.000141 * sin(   621.9339951688 *T + 0.679068671 )   + \
        0.000118 * sin(  2274.3409379516 *T + 4.881123092 )   + \
        0.000129 * sin(   169.2165669502 *T + 0.351407289 )   + \
        0.000126 * sin(   565.7405657679 *T + 5.146592349 )   + \
        0.000114 * sin(    72.8762966531 *T + 0.520791814 )   + \
        0.000120 * sin(     5.2596639600 *T + 0.948516300 )   + \
        0.000115 * sin(     6.5220371012 *T + 3.504914846 )   + \
        0.000126 * sin(   588.1403728234 *T + 5.577502482 )   + \
        0.000158 * sin( 16309.6180360983 *T + 2.957128968 )   + \
        0.000134 * sin(  1234.1806904281 *T + 2.598576764 )   + \
        0.000151 * sin(  1662.7370915377 *T + 3.985702050 )   + \
        0.000109 * sin(   136.8660252845 *T + 0.014730471 )   + \
        0.000131 * sin(   621.1263196841 *T + 0.085077024 )   + \
        0.000146 * sin(   579.2741760812 *T + 0.708426604 )   + \
        0.000146 * sin(    -7.7750543984 *T + 3.121576600 )   + \
        0.000107 * sin(   534.1013788022 *T + 0.288231904 )   + \
        0.000138 * sin(   628.1591377283 *T + 2.797450317 )

    t23 = 0.000113 * sin(  -627.7552925684 *T + 2.788904128 )   + \
        0.000115 * sin(   -52.5758811831 *T + 5.895222200 )   + \
        0.000138 * sin(   601.6468808270 *T + 6.096188999 )   + \
        0.000139 * sin(  2353.9707386333 *T + 2.028195445 )   + \
        0.000146 * sin(  -417.6041342449 *T + 4.660008502 )   + \
        0.000107 * sin(  1606.2184526117 *T + 4.066520001 )   + \
        0.000142 * sin(  8378.3548222473 *T + 2.936315115 )   + \
        0.000128 * sin(   938.0959672717 *T + 3.223844306 )   + \
        0.000135 * sin(   620.5325306007 *T + 1.638054048 )   + \
        0.000101 * sin(   269.9734819318 *T + 5.481603249 )   + \
        0.000104 * sin(   -56.8821874027 *T + 2.205734493 )   + \
        0.000103 * sin(   632.1103522627 *T + 2.440421099 )   + \
        0.000119 * sin(   632.1208885629 *T + 2.547496264 )   + \
        0.000138 * sin(   197.5492545856 *T + 2.314608466 )   + \
        0.000121 * sin(    13.7033024162 *T + 4.539108237 )   + \
        0.000123 * sin(  1940.2796952817 *T + 4.538074405 )   + \
        0.000119 * sin(  2280.5735565994 *T + 2.869040566 )   + \
        0.000133 * sin(  6447.1991241142 *T + 6.056405489 )   + \
        0.000129 * sin(    -8.5827298831 *T + 2.540635083 )   + \
        0.000131 * sin(  1361.3804277336 *T + 4.005732868 )

    t24 =  0.000104 * sin(   981.4604100291 *T + 1.959967212 )   + \
        0.000112 * sin(  1609.7679950283 *T + 3.589026260 )   + \
        0.000123 * sin(   210.7034507542 *T + 1.728627253 )   + \
        0.000121 * sin(  3694.9230808424 *T + 6.072332087 )   + \
        0.000108 * sin( -1253.9853380183 *T + 3.716133846 )   + \
        0.000113 * sin(  -787.5671863624 *T + 2.725771122 )   + \
        0.000109 * sin(   417.1425536614 *T + 4.033338079 )   + \
        0.000101 * sin(   624.7911759770 *T + 3.441347021 )   + \
        0.000113 * sin(   733.0728427345 *T + 0.656372122 )   + \
        0.000113 * sin(  5109.2726050855 *T + 2.791483066 )   + \
        0.000106 * sin(   562.1842923210 *T + 1.815323326 )   + \
        0.000101 * sin(    11.1430161497 *T + 5.711033677 )   + \
        0.000103 * sin(    90.9818733055 *T + 2.812745443 )   + \
        0.000101 * sin(   179.0642637886 *T + 1.965746028 )
    # --------------------------------------------------------------------
    # T**1
    t24s = \
      10.2156724 * sin(  628.3075849991 *T + 4.249032005 )   + \
       0.1706807 * sin( 1256.6151699983 *T + 4.205904248 )   + \
       0.0269668 * sin(   21.3299095438 *T + 3.400290479 )   + \
       0.0265919 * sin(   52.9690965095 *T + 5.836047367 )   + \
       0.0210568 * sin(   -0.3523118349 *T + 6.262738348 )   + \
       0.0077996 * sin(  522.3693919802 *T + 4.670344204 )

    t25 = \
       0.0054764 * sin(  157.7343542448 *T + 4.534800170 )   + \
       0.0059146 * sin(    2.6298319800 *T + 1.083044735 )   + \
       0.0034420 * sin(  -39.8149003408 *T + 5.980077351 )   + \
       0.0032088 * sin( 1884.9227549974 *T + 4.162913471 )   + \
       0.0033595 * sin(  550.7553238667 *T + 5.980162321 )   + \
       0.0029198 * sin(  585.6477659115 *T + 0.623811863 )   + \
       0.0027764 * sin(   15.5420399434 *T + 3.745318113 )   + \
       0.0025190 * sin(  574.6271337896 *T + 2.980330535 )   + \
       0.0022997 * sin(  -79.6298006816 *T + 1.174411803 )   + \
       0.0024976 * sin(  576.0498431898 *T + 2.467913690 )   + \
       0.0021774 * sin(   20.6185548437 *T + 3.854787540 )   + \
       0.0017925 * sin(  -77.5522611324 *T + 1.092065955 )   + \
       0.0013794 * sin(   42.6598190876 *T + 2.699831988 )   + \
       0.0013276 * sin(  606.2663207553 *T + 5.845801920 )   + \
       0.0011774 * sin( 1203.6460734888 *T + 2.292832062 )   + \
       0.0012869 * sin(  607.6890301554 *T + 5.333425680 )   + \
       0.0012152 * sin(  105.9381930189 *T + 6.222874454 )   + \
       0.0011081 * sin(   -0.7113547001 *T + 5.154724984 )   + \
       0.0010143 * sin(  469.4002954708 *T + 4.044013795 )   + \
       0.0009357 * sin(  548.6777843175 *T + 3.416081409 )

    t26 = \
       0.0010084 * sin(   52.2577418094 *T + 0.749320262 )   + \
       0.0008587 * sin( 1097.7078804699 *T + 2.777152598 )   + \
       0.0008628 * sin(  627.5962302991 *T + 4.562060226 )   + \
       0.0008158 * sin(  -22.0412642439 *T + 5.806891533 )   + \
       0.0007746 * sin(  254.4314419883 *T + 1.603197066 )   + \
       0.0007670 * sin(  214.6165416475 *T + 3.000200440 )   + \
       0.0007098 * sin(    7.4781598567 *T + 0.443725817 )   + \
       0.0006180 * sin(  -53.6804512095 *T + 1.302642751 )   + \
       0.0005818 * sin(  508.8628839767 *T + 4.827723531 )   + \
       0.0004945 * sin( -628.6598968340 *T + 0.268305170 )   + \
       0.0004774 * sin(  134.9867409659 *T + 5.808636673 )   + \
       0.0004687 * sin(  -24.2728603974 *T + 5.154890570 )   + \
       0.0006089 * sin(  174.8016413067 *T + 4.403765209 )   + \
       0.0005975 * sin( -119.4447010225 *T + 2.583472591 )   + \
       0.0004229 * sin(   95.1718406251 *T + 0.931172179 )   + \
       0.0005264 * sin(   55.3569402842 *T + 2.336107252 )   + \
       0.0003049 * sin(  564.3178563677 *T + 1.362634430 )   + \
       0.0002974 * sin(  681.2766815086 *T + 1.583012668 )   + \
       0.0003403 * sin( -235.2866153772 *T + 2.552189886 )   + \
       0.0003030 * sin(   41.9484643875 *T + 5.286473844 )

    t27 = \
       0.0003210 * sin(   -0.7046236698 *T + 1.863796539 )   + \
       0.0003058 * sin(  943.7762934887 *T + 4.226420633 )   + \
       0.0002589 * sin( 1235.2852604545 *T + 1.991935820 )   + \
       0.0002927 * sin(  521.6580372801 *T + 2.319951253 )   + \
       0.0002425 * sin(  523.0807466803 *T + 3.084752833 )   + \
       0.0002656 * sin(  315.4687084896 *T + 2.487447866 )   + \
       0.0002445 * sin( 1044.7387839604 *T + 2.347139160 )   + \
       0.0002990 * sin(  469.0479836359 *T + 6.235872050 )   + \
       0.0002890 * sin(  586.3591206116 *T + 0.095197563 )   + \
       0.0002498 * sin(  643.8496249426 *T + 2.994779800 )   + \
       0.0001889 * sin(  803.1092263058 *T + 3.569003717 )   + \
       0.0002567 * sin(   80.1820931124 *T + 3.425611498 )   + \
       0.0001803 * sin(-7143.0695617928 *T + 2.192295512 )   + \
       0.0001782 * sin(    0.3932153263 *T + 5.180433689 )   + \
       0.0001694 * sin( -470.5732307544 *T + 4.641779174 )   + \
       0.0001704 * sin( -159.2596013633 *T + 3.997097652 )   + \
       0.0001735 * sin(  584.9364112115 *T + 0.417558428 )   + \
       0.0001643 * sin(  842.9241266467 *T + 2.180619584 )   + \
       0.0001680 * sin(    3.8133035638 *T + 4.164529426 )   + \
       0.0002045 * sin(  708.4896781115 *T + 0.526323854 )

    t28 = \
       0.0001458 * sin(  429.2330832950 *T + 1.356098141 )   + \
       0.0001437 * sin(    2.0355319399 *T + 3.895439360 )   + \
       0.0001738 * sin(  627.9552731642 *T + 0.087484036 )   + \
       0.0001367 * sin( 1414.3495242431 *T + 3.987576591 )   + \
       0.0001344 * sin(  723.4794256242 *T + 0.090454338 )   + \
       0.0001438 * sin( 1149.9656222793 *T + 0.974387904 )   + \
       0.0001257 * sin(  683.6645252834 *T + 1.509069366 )   + \
       0.0001358 * sin( 1151.3883316794 *T + 0.495572260 )   + \
       0.0001628 * sin(  763.2943259650 *T + 4.968445721 )   + \
       0.0001169 * sin(   10.3092774219 *T + 2.838496795 )   + \
       0.0001162 * sin(  416.4311989613 *T + 3.408387778 )   + \
       0.0001092 * sin(  606.9776754553 *T + 3.617942651 )   + \
       0.0001008 * sin( 1778.9845619785 *T + 0.286350174 )   + \
       0.0001008 * sin(   63.9897286314 *T + 1.610762073 )   + \
       0.0000918 * sin( 1021.3285546211 *T + 5.532798067 )   + \
       0.0001011 * sin( -625.6777530192 *T + 0.661826484 )   + \
       0.0000753 * sin( 1673.0463689596 *T + 3.905030235 )   + \
       0.0000737 * sin( 1192.6254413669 *T + 4.641956361 )   + \
       0.0000694 * sin(  334.0612426700 *T + 2.111120332 )   + \
       0.0000701 * sin(  389.4181829542 *T + 2.760823491 )

    t29 = \
       0.0000689 * sin(  -13.5065080035 *T + 4.768800780 )   + \
       0.0000700 * sin( 1336.7972631107 *T + 5.760439898 )   + \
       0.0000664 * sin(  604.0347246017 *T + 1.051215840 )   + \
       0.0000654 * sin(  565.0292110678 *T + 4.911332503 )   + \
       0.0000788 * sin(  668.1224853400 *T + 4.699648011 )   + \
       0.0000628 * sin(  533.3900241022 *T + 5.024608847 )   + \
       0.0000755 * sin(  -11.0206321219 *T + 4.370971253 )   + \
       0.0000628 * sin(  629.0189396992 *T + 3.660478857 )   + \
       0.0000635 * sin( 2513.2303399966 *T + 4.121051532 )   + \
       0.0000534 * sin(  596.6683980335 *T + 1.173284524 )   + \
       0.0000543 * sin(  -43.3711737877 *T + 0.345585464 )   + \
       0.0000517 * sin( -199.0745017041 *T + 5.414571768 )   + \
       0.0000504 * sin(  576.7611978898 *T + 2.328281115 )   + \
       0.0000485 * sin(  575.3384884897 *T + 1.685874771 )   + \
       0.0000463 * sin(  786.0419392439 *T + 5.297703006 )   + \
       0.0000604 * sin(   51.5463871093 *T + 0.591998446 )   + \
       0.0000443 * sin( 1216.8002696575 *T + 4.830881244 )   + \
       0.0000570 * sin(   19.9072001436 *T + 3.899190272 )   + \
       0.0000465 * sin( 1096.9965257698 *T + 0.476681802 )   + \
       0.0000424 * sin( -707.9373856808 *T + 1.112242763 )

    t30 = \
       0.0000427 * sin(   73.5876513532 *T + 1.994214480 )   + \
       0.0000478 * sin( -612.7655450557 *T + 3.778025483 )   + \
       0.0000414 * sin( 1097.3555686350 *T + 5.441088327 )   + \
       0.0000512 * sin(  158.9072895284 *T + 0.107123853 )   + \
       0.0000378 * sin( 1098.4192351700 *T + 0.915087231 )   + \
       0.0000402 * sin( 1137.1704689758 *T + 4.107281715 )   + \
       0.0000453 * sin(  991.7696874510 *T + 1.917490952 )   + \
       0.0000395 * sin(   14.9563197135 *T + 2.763124165 )   + \
       0.0000371 * sin(  573.9157790895 *T + 3.112111866 )   + \
       0.0000350 * sin( 1179.0629088659 *T + 0.440639857 )   + \
       0.0000356 * sin(  613.3512652857 *T + 5.444568842 )   + \
       0.0000344 * sin(   41.2371096874 *T + 5.676832684 )   + \
       0.0000383 * sin(   95.5599741609 *T + 5.559734846 )   + \
       0.0000333 * sin(  649.6374945429 *T + 0.261537984 )   + \
       0.0000340 * sin(  605.5549660552 *T + 5.975534987 )   + \
       0.0000334 * sin(  106.6495477190 *T + 2.335063907 )   + \
       0.0000399 * sin( 1150.6769769794 *T + 5.321230910 )   + \
       0.0000314 * sin( 1831.9536584880 *T + 2.313312404 )   + \
       0.0000424 * sin(  105.2268383188 *T + 1.211961766 )   + \
       0.0000307 * sin(    6.3735898303 *T + 3.169551388 )

    t31 = \
       0.0000329 * sin(    2.9821438149 *T + 6.106912080 )   + \
       0.0000357 * sin(  630.9374169791 *T + 4.223760346 )   + \
       0.0000312 * sin( -373.8761430108 *T + 2.180556645 )   + \
       0.0000301 * sin(   30.9278322656 *T + 1.499984572 )   + \
       0.0000268 * sin( 1204.3574281889 *T + 2.447520648 )   + \
       0.0000257 * sin( 1249.1370101415 *T + 3.662331761 )   + \
       0.0000290 * sin(   62.5670192312 *T + 1.272834584 )   + \
       0.0000256 * sin(  542.9879468239 *T + 1.913426912 )   + \
       0.0000339 * sin(  349.6032826134 *T + 4.165930011 )   + \
       0.0000283 * sin(  393.0209696220 *T + 4.325565754 )   + \
       0.0000241 * sin( 1252.8018664345 *T + 3.832324536 )   + \
       0.0000304 * sin(  468.6889407707 *T + 1.612348468 )   + \
       0.0000259 * sin( 1620.0772724501 *T + 3.470173146 )   + \
       0.0000238 * sin( 1213.9553509107 *T + 1.147977842 )   + \
       0.0000236 * sin(  617.2869528772 *T + 3.776271728 )   + \
       0.0000296 * sin( -705.8598461315 *T + 0.460368852 )   + \
       0.0000306 * sin( 1057.5406682942 *T + 0.554749016 )   + \
       0.0000251 * sin( 1729.8182327326 *T + 0.834332510 )   + \
       0.0000290 * sin(  473.2030627343 *T + 4.759564091 )   + \
       0.0000261 * sin(  588.4926846583 *T + 0.298259862 )

    t32 = \
       0.0000249 * sin(  554.7199336460 *T + 3.749366406 )   + \
       0.0000213 * sin( 1171.2955318231 *T + 5.415666119 )   + \
       0.0000223 * sin(  470.1116501708 *T + 2.703203558 )   + \
       0.0000268 * sin(  -64.0877607382 *T + 0.283670793 )   + \
       0.0000209 * sin(  563.6065016677 *T + 1.238477199 )   + \
       0.0000193 * sin( 1017.7257679534 *T + 1.943251340 )   + \
       0.0000182 * sin(  628.3143160294 *T + 2.456157599 )   + \
       0.0000184 * sin(  -22.7526189440 *T + 5.888038582 )   + \
       0.0000182 * sin( -628.3008539689 *T + 0.241332086 )   + \
       0.0000228 * sin( -628.4056171060 *T + 2.657323816 )   + \
       0.0000166 * sin(  723.8675591600 *T + 5.930629110 )   + \
       0.0000167 * sin(  309.7883822726 *T + 5.570955333 )   + \
       0.0000159 * sin(  -32.3505416657 *T + 5.786670700 )   + \
       0.0000154 * sin( -413.6910433516 *T + 1.517805532 )   + \
       0.0000176 * sin( 1202.9347187887 *T + 3.139266834 )   + \
       0.0000167 * sin( 1213.2439962106 *T + 3.556352289 )   + \
       0.0000153 * sin(   20.2253395174 *T + 1.463313961 )   + \
       0.0000157 * sin( 1726.7268201691 *T + 1.586837396 )   + \
       0.0000142 * sin( 8399.6847317911 *T + 0.022670115 )   + \
       0.0000152 * sin( 1726.0154654690 *T + 0.708528947 )

    t33 = \
       0.0000144 * sin(  608.4003848555 *T + 5.187075177 )   + \
       0.0000135 * sin(  575.6566278634 *T + 1.993229262 )   + \
       0.0000134 * sin(  575.0203491159 *T + 3.457197134 )   + \
       0.0000144 * sin(  532.6786694021 *T + 6.066193291 )   + \
       0.0000160 * sin( 1101.5106477335 *T + 1.710431974 )   + \
       0.0000133 * sin(  363.4621024518 *T + 2.836451652 )   + \
       0.0000134 * sin( 1807.3704938650 *T + 5.453106665 )   + \
       0.0000134 * sin(  116.2474704408 *T + 5.326898811 )   + \
       0.0000128 * sin(  564.2198242609 *T + 2.511652591 )   + \
       0.0000160 * sin(   63.2783739313 *T + 5.628785365 )   + \
       0.0000132 * sin( 1391.6019109642 *T + 0.819294053 )   + \
       0.0000122 * sin( 1431.4168113050 *T + 5.677408071 )   + \
       0.0000125 * sin( 1235.9966151546 *T + 5.251984735 )   + \
       0.0000121 * sin(  574.9452731634 *T + 2.210924603 )   + \
       0.0000136 * sin(  -24.5831646229 *T + 1.646502367 )   + \
       0.0000120 * sin(  575.7317038160 *T + 3.240883049 )   + \
       0.0000134 * sin( 1214.6667056108 *T + 3.059480037 )   + \
       0.0000137 * sin(  620.6809778716 *T + 1.867105418 )   + \
       0.0000141 * sin( 1725.3041107690 *T + 2.069217456 )   + \
       0.0000129 * sin( -747.7522860216 *T + 2.781469314 )

    t34 = \
       0.0000116 * sin(  554.0085789459 *T + 4.281176991 )   + \
       0.0000116 * sin(  977.9108676125 *T + 3.320925381 )   + \
       0.0000129 * sin(  523.7921013804 *T + 3.497704076 )   + \
       0.0000113 * sin(  595.9570433334 *T + 0.983210840 )   + \
       0.0000122 * sin(  628.2095528923 *T + 2.674938860 )   + \
       0.0000140 * sin(   -1.1045700264 *T + 4.957936982 )   + \
       0.0000108 * sin( 2354.3230504682 *T + 1.390113589 )   + \
       0.0000106 * sin(-1256.9674818332 *T + 0.429631317 )   + \
       0.0000110 * sin(  -26.6607041722 *T + 5.501340197 )   + \
       0.0000115 * sin( 1255.9038152982 *T + 4.691456618 )   + \
       0.0000134 * sin( -238.8894020449 *T + 0.577313584 )   + \
       0.0000109 * sin( 1044.0274292604 *T + 6.218148717 )   + \
       0.0000102 * sin(  -54.3918059096 *T + 1.477842615 )   + \
       0.0000108 * sin( 2122.8392023546 *T + 2.237753948 )   + \
       0.0000101 * sin( -453.5059436924 *T + 3.100492232 )   + \
       0.0000103 * sin(    7.6266071276 *T + 5.594294322 )   + \
       0.0000104 * sin(   94.9175608970 *T + 5.674287810 )   + \
       0.0000101 * sin( 1351.7870106233 *T + 2.196632348 )   + \
       0.0000100 * sin( 1193.3367960670 *T + 4.056084160 )   
    # --------------------------------------------------------------------
    # T**2

    t35 = \
      0.04322990 * sin(  628.3075849991 *T + 2.642893748 )  + \
      0.00406495 * sin(    0.0000000000 *T + 4.712388980 )  + \
      0.00122605 * sin( 1256.6151699983 *T + 2.438140634 )  + \
      0.00019476 * sin(   21.3299095438 *T + 1.642186981 )  + \
      0.00016916 * sin(   52.9690965095 *T + 4.510959344 )  + \
      0.00013374 * sin(   -0.3523118349 *T + 1.502210314 )  + \
      0.00008042 * sin(    2.6298319800 *T + 0.478549024 )  + \
      0.00007824 * sin(   15.5420399434 *T + 5.254710405 )  + \
      0.00004894 * sin(  574.6271337896 *T + 4.683210850 )  + \
      0.00004875 * sin(  576.0498431898 *T + 0.759507698 )  + \
      0.00004416 * sin(  522.3693919802 *T + 6.028853166 )  + \
      0.00004088 * sin(   -0.7113547001 *T + 0.060926389 )  + \
      0.00004433 * sin( 7771.3771467920 *T + 3.627734103 )  + \
      0.00003277 * sin( 1884.9227549974 *T + 2.327912542 )  + \
      0.00002703 * sin(  606.2663207553 *T + 1.271941729 )  + \
      0.00003435 * sin(  -77.5522611324 *T + 0.747446224 )  + \
      0.00002618 * sin(  607.6890301554 *T + 3.633715689 )  + \
      0.00003146 * sin(   20.6185548437 *T + 5.647874613 )  + \
      0.00002544 * sin(  157.7343542448 *T + 6.232904270 )  + \
      0.00002218 * sin(  -22.0412642439 *T + 1.309509946 )  + \
      0.00002197 * sin(  585.6477659115 *T + 2.407212349 )

    t36 = \
      0.00002897 * sin(  575.3384884897 *T + 5.863842246 )  + \
      0.00001766 * sin(   42.6598190876 *T + 0.754113147 )  + \
      0.00001738 * sin(  -79.6298006816 *T + 2.714942671 )  + \
      0.00001695 * sin(   52.2577418094 *T + 2.629369842 )  + \
      0.00001584 * sin(  550.7553238667 *T + 1.341138229 )  + \
      0.00001503 * sin(  -24.2728603974 *T + 0.377699736 )  + \
      0.00001552 * sin(  -53.6804512095 *T + 2.904684667 )  + \
      0.00001370 * sin(  -39.8149003408 *T + 1.265599125 )  + \
      0.00001889 * sin( -557.3142801634 *T + 4.413514859 )  + \
      0.00001722 * sin(  606.9776754553 *T + 2.445966339 )  + \
      0.00001124 * sin(  105.9381930189 *T + 5.041799657 )  + \
      0.00001258 * sin(   55.3569402842 *T + 3.849557278 )  + \
      0.00000831 * sin(   95.1718406251 *T + 2.471094709 )  + \
      0.00000767 * sin(  469.4002954708 *T + 5.363125422 )  + \
      0.00000756 * sin(  134.9867409659 *T + 1.046195744 )  + \
      0.00000775 * sin(   -1.1045700264 *T + 0.245548001 )  + \
      0.00000597 * sin(  214.6165416475 *T + 4.543268798 )  + \
      0.00000568 * sin(  521.6580372801 *T + 4.178853144 )  + \
      0.00000711 * sin(  174.8016413067 *T + 5.934271972 )  + \
      0.00000499 * sin( 1203.6460734888 *T + 0.624434410 )

    t37 = \
      0.00000671 * sin( -119.4447010225 *T + 4.136047594 )  + \
      0.00000488 * sin(  584.9364112115 *T + 2.209679987 )  + \
      0.00000621 * sin(  643.8496249426 *T + 4.518860804 )  + \
      0.00000495 * sin( -628.6598968340 *T + 1.868201275 )  + \
      0.00000456 * sin(  523.0807466803 *T + 1.271231591 )  + \
      0.00000451 * sin(  508.8628839767 *T + 0.084060889 )  + \
      0.00000435 * sin(  564.3178563677 *T + 3.324456609 )  + \
      0.00000387 * sin( 1097.7078804699 *T + 4.052488477 )  + \
      0.00000547 * sin(16100.0685737473 *T + 2.841633844 )  + \
      0.00000522 * sin(  315.4687084896 *T + 2.171979966 )  + \
      0.00000375 * sin(  548.6777843175 *T + 4.983027306 )  + \
      0.00000421 * sin(  586.3591206116 *T + 4.546432249 )  + \
      0.00000439 * sin(  708.4896781115 *T + 0.522967921 )  + \
      0.00000309 * sin(  254.4314419883 *T + 3.172606705 )  + \
      0.00000347 * sin(  469.0479836359 *T + 1.479586566 )  + \
      0.00000317 * sin(   80.1820931124 *T + 3.553088096 )  + \
      0.00000262 * sin(   41.9484643875 *T + 0.606635550 )  + \
      0.00000248 * sin(  683.6645252834 *T + 3.014082064 )  + \
      0.00000245 * sin( -159.2596013633 *T + 5.519526220 )  + \
      0.00000225 * sin(  429.2330832950 *T + 2.877956536 )

    t38 = \
      0.00000214 * sin(  723.4794256242 *T + 1.605227587 )   + \
      0.00000205 * sin(  576.7611978898 *T + 0.625804796 )   + \
      0.00000180 * sin( 1044.7387839604 *T + 3.499954526 )   + \
      0.00000229 * sin(   19.9072001436 *T + 5.632304604 )   + \
      0.00000214 * sin(   63.9897286314 *T + 5.960227667 )   + \
      0.00000175 * sin(  -43.3711737877 *T + 2.162417992 )   + \
      0.00000209 * sin(   51.5463871093 *T + 2.322150893 )   + \
      0.00000173 * sin(  604.0347246017 *T + 2.556183691 )   + \
      0.00000184 * sin(  630.9374169791 *T + 4.732296790 )   + \
      0.00000227 * sin(14985.4400134205 *T + 5.385812217 )   + \
      0.00000154 * sin(  803.1092263058 *T + 5.120720920 )   + \
      0.00000151 * sin(  573.9157790895 *T + 4.815000443 )   + \
      0.00000197 * sin(  763.2943259650 *T + 0.222827271 )   + \
      0.00000197 * sin(    7.4781598567 *T + 3.910456770 )   + \
      0.00000138 * sin(  605.5549660552 *T + 1.397484253 )   + \
      0.00000149 * sin( -612.7655450557 *T + 5.333727496 )   + \
      0.00000137 * sin(  389.4181829542 *T + 4.281749907 )   + \
      0.00000135 * sin(  943.7762934887 *T + 5.979971885 )   + \
      0.00000139 * sin( -235.2866153772 *T + 4.715630782 )   + \
      0.00000142 * sin(  681.2766815086 *T + 0.513330157 )

    t39 = \
      0.00000120 * sin( -470.5732307544 *T + 0.194160689 )  + \
      0.00000131 * sin(-7143.0695617928 *T + 0.000379226 )  + \
      0.00000124 * sin(  627.9552731642 *T + 2.122264908 )  + \
      0.00000108 * sin( -625.6777530192 *T + 0.883445696 )
    # --------------------------------------------------------------------
    # T**3
    t39s = \
      0.143388e-3  * sin(  628.3075849991 *T + 1.131453581 )  + \
      0.006671e-3  * sin( 1256.6151699983 *T + 0.775148887 )  + \
      0.001480e-3  * sin(   15.5420399434 *T + 0.480016880 )  + \
      0.000934e-3  * sin(   21.3299095438 *T + 6.144453084 )  + \
      0.000795e-3  * sin(   52.9690965095 *T + 2.941595619 )  + \
      0.000673e-3  * sin(  574.6271337896 *T + 0.120415406 )  + \
      0.000672e-3  * sin(  576.0498431898 *T + 5.317009738 )  + \
      0.000389e-3  * sin(  -22.0412642439 *T + 3.090323467 )  + \
      0.000373e-3  * sin(  606.2663207553 *T + 3.003551964 )  + \
      0.000360e-3  * sin(  607.6890301554 *T + 1.918913041 )  + \
      0.000316e-3  * sin(   -2.1340641002 *T + 5.545798121 )  + \
      0.000315e-3  * sin(  -24.2728603974 *T + 1.884932563 )  + \
      0.000278e-3  * sin(   20.6185548437 *T + 1.266254859 )  + \
      0.000238e-3  * sin(  -53.6804512095 *T + 4.532664830 )  + \
      0.000185e-3  * sin(   52.2577418094 *T + 4.578313856 )  + \
      0.000245e-3  * sin( 1884.9227549974 *T + 0.587467082 )

    t40 = \
      0.000180e-3  * sin(   42.6598190876 *T + 5.151178553 )  + \
      0.000200e-3  * sin(   55.3569402842 *T + 5.355983739 )  + \
      0.000141e-3  * sin(  522.3693919802 *T + 1.336556009 )  + \
      0.000104e-3  * sin(  585.6477659115 *T + 4.239842759 )
    # -------------------------------------------------------------------
    # T**4
    t40s = \
      0.003826e-4  * sin(  628.3075849991 *T + 5.705257275 )  + \
      0.000303e-4  * sin( 1256.6151699983 *T + 5.407132842 )  + \
      0.000209e-4  * sin(   15.5420399434 *T + 1.989815753 )

    # Correction to use JPL planetary masses instead of IAU.
    t41 = \
         0.00065   * sin(  606.9776754    *T + 4.021194  )    + \
         0.00033   * sin(   21.3299095    *T + 5.543132  )    - \
         0.00196   * sin(  620.8294251    *T + 5.696701  )    - \
         0.00173   * sin(    7.4781599    *T + 2.435900  )    + \
         0.0003638 * T_2

    per_terms = t41 + t40s*T_4 + (t39s+t40)*T_3 + (t39+t38+t37+t36+t35)*T_2 + \
                (t34+t33+t32+t31+t30+t29+t28+t27+t26+t25+t24s)*T + \
                t24+t23+t22+t21+t20+t19+t18+t17+t16+t15+t14+t13 + \
                t12+t11+t10 +t9 +t8 +t7 +t6 +t5 +t4 +t3 +t2 +t1

    per_terms = per_terms * 1e-6;  # Periodic terms  in Sec

    #   Accuracy of calculation of periodic terms is of order of 2x10**{-14} sec
    #   for one Julian century. Accuracy of calculation of derivatives with
    #   respect to coordinate time was limited by a value of 3x10**{-17).
    # 
    #   Derivatives [microsec/Julian Century] were calculated by MAPLE software.
    #   Order of terms can differ from the original order.
    # 
    #   Conversion factor from   [microsec/Julian Century] to [Sec/Sec] is equal
    #   to 10**(-6)/36525/86400 = 3x10**(-16). Only terms with amplitudes more
    #   than 0.1 [microsec/Julian Century] were kept for calculation.

    dt1 = \
      1040901.1944362769 *cos(   628.3075849991 *T+6.240054195 )+ \
       12897.6338809017 *cos(   575.3384884897 *T+4.296977442 )+ \
       17391.2925768211 *cos(  1256.6151699983 *T+6.196904410 )+ \
         252.6671456926 *cos(    52.9690965095 *T+0.444401603 )+ \
        2838.6767739088 *cos(   606.9776754553 *T+4.021195093 )+ \
          48.1353561769 *cos(    21.3299095438 *T+5.543113262 )- \
           0.5968884722 *cos(     0.3523118349 *T-5.025132748 )+ \
       12083.7531824326 *cos(  7771.3771467920 *T+5.198467090 )+ \
        1003.6490036622 *cos(   786.0419392439 *T+5.988822341 )+ \
         623.3846626319 *cos(   522.3693919802 *T+3.649823730 )+ \
         438.3449338807 *cos(   393.0209696220 *T+1.422745069 )+ \
         913.8503949624 *cos(  1150.6769769794 *T+2.322313077 )+ \
           1.1756953148 *cos(     2.6298319800 *T+3.615796498 )- \
          17.3276835177 *cos(    39.8149003408 *T-4.349338347 )+ \
          94.6893524623 *cos(   157.7343542448 *T+2.678271909 )+ \
         308.4386125110 *cos(   620.8294251424 *T+5.696701824 )+ \
         286.1875235054 *cos(   588.4926846583 *T+0.520007179 )+ \
           3.2334964968 *cos(     7.4781598567 *T+2.435898309 )+ \
         292.6361467978 *cos(   624.4942814354 *T+5.866398759 )+ \
         206.8141316652 *cos(   550.7553238667 *T+4.103476804 )

    dt2 = \
         -18.8517913974 *cos(    77.5522611324 *T-3.651837925 )+ \
         326.9115780130 *cos(  1884.9227549974 *T+6.153743485 )+ \
         135.1001548793 *cos(   585.6477659115 *T+4.773852582 )+ \
         245.2392765351 *cos(  1203.6460734888 *T+4.333987818 )- \
          11.4615153611 *cos(    79.6298006816 *T-5.957517795 )+ \
         174.6233696252 *cos(  1097.7078804699 *T+1.890075226 )+ \
           0.4575163483 *cos(     3.8133035638 *T+4.551585768 )+ \
          65.2767446780 *cos(   548.6777843175 *T+1.914547226 )+ \
          12.3015429734 *cos(   105.9381930189 *T+0.873504123 )+ \
         162.6246098311 *cos(  1179.0629088659 *T+1.135934669 )+ \
          25.0253677711 *cos(   254.4314419883 *T+0.092793886 )- \
          56.7724910917 *cos(   557.3142801634 *T-5.984503847 )+ \
           1.6528658305 *cos(    20.6185548437 *T+2.095377709 )+ \
          37.3853865328 *cos(   469.4002954708 *T+2.949233637 )+ \
           0.1300892940 *cos(     2.0775395492 *T+2.654394814 )+ \
          22.0740663552 *cos(   294.2463423292 *T+4.980931759 )+ \
          37.0042635346 *cos(   574.6271337896 *T+1.280308748 )+ \
          36.7600446933 *cos(   576.0498431898 *T+4.167901731 )+ \
          10.3106078938 *cos(   214.6165416475 *T+1.495846011 )+ \
           0.7518150982 *cos(    15.5420399434 *T+2.251573730 )

    dt3 = \
           2.5102743944 *cos(    42.6598190876 *T+4.839650148 )+ \
          93.4447512850 *cos(  1726.0154654690 *T+3.411091093 )+ \
          26.6169837232 *cos(   627.5962302991 *T+2.869567043 )+ \
          18.6060624897 *cos(   508.8628839767 *T+3.324679049 )+ \
          50.3489919309 *cos(  1235.2852604545 *T+3.981496998 )+ \
           2.9272076733 *cos(    80.1820931124 *T+6.248866009 )+ \
          11.6581461222 *cos(   315.4687084896 *T+5.071801441 )+ \
           2.7040114748 *cos(    63.2783739313 *T+5.720622217 )+ \
         685.2189184987 *cos( 16100.0685737473 *T+1.270837679 )+ \
          63.6379554012 *cos(  1572.0838784878 *T+2.546610123 )- \
          17.7558701262 *cos(   628.6598968340 *T-5.069663519 )+ \
          20.2959776199 *cos(   606.2663207553 *T+4.144987272 )+ \
           1.8220706837 *cos(    52.2577418094 *T+5.210064075 )+ \
          19.7122167602 *cos(   607.6890301554 *T+0.749317412 )+ \
          21.4070156241 *cos(   708.4896781115 *T+3.389610345 )- \
         208.9133554738 *cos(  7143.0695617928 *T-4.183178762 )+ \
          31.6438753444 *cos(   943.7762934887 *T+2.404714239 )+ \
          28.6210474720 *cos(   882.7390269875 *T+5.541473556 )

    dt4 = \
          17.3108430153 *cos(   627.9552731642 *T+5.040846034 )+ \
          36.2511346889 *cos(  1213.9553509107 *T+1.770181024 )+ \
          23.5160252882 *cos(  1044.7387839604 *T+1.460726241 )+ \
          17.6483024396 *cos(   842.9241266467 *T+0.652303414 )+ \
           0.8524766933 *cos(    41.9484643875 *T+3.735430632 )- \
           2.9641397006 *cos(   119.4447010225 *T-1.087136918 )+ \
           4.4043021544 *cos(   174.8016413067 *T+2.901883301 )+ \
          30.6786555304 *cos(  1414.3495242431 *T+5.952658009 )+ \
          12.0402027923 *cos(   681.2766815086 *T+3.186129845 )+ \
          13.8414980037 *cos(   613.3512652857 *T+3.307984806 )+ \
          16.4995627999 *cos(  1021.3285546211 *T+1.331103168 )+ \
           1.9911894160 *cos(   134.9867409659 *T+4.308933301 )- \
           0.3515361234 *cos(    22.0412642439 *T-4.005298270 )- \
           3.7584683940 *cos(   235.2866153772 *T-6.145309371 )+ \
          25.3024974250 *cos(  1778.9845619785 *T+2.104551349 )+ \
           0.1305128623 *cos(     7.3297125859 *T+3.475975097 )- \
           0.7338654485 *cos(    53.6804512095 *T-5.971672571 )+ \
           9.5907303805 *cos(   803.1092263058 *T+2.053414715 )+ \
          23.9546779108 *cos(  1673.0463689596 *T+3.016058075 )+ \
           0.1284742152 *cos(    10.3092774219 *T+1.737438797 )

    dt5 = \
          29.6298508998 *cos(  1965.1048481098 *T+3.969480770 )+ \
           0.9894064551 *cos(    95.1718406251 *T+5.717799605 )- \
           5.5090008124 *cos(   470.5732307544 *T-2.654125618 )+ \
           6.1292118878 *cos(   586.3591206116 *T+1.913704550 )+ \
           5.8255759568 *cos(   469.0479836359 *T+4.734090399 )+ \
           6.6854736444 *cos(   564.3178563677 *T+5.489005403 )+ \
           2.8762672994 *cos(   334.0612426700 *T+3.661698944 )+ \
           5.9511626514 *cos(   512.0601145584 *T+4.863931876 )+ \
           0.5992388786 *cos(    55.3569402842 *T+0.842715011 )- \
           0.1170473984 *cos(    13.5065080035 *T-3.293406547 )+ \
           0.1490098133 *cos(    14.9563197135 *T+4.870690598 )+ \
           6.2197810566 *cos(   630.9374169791 *T+1.061816410 )+ \
           0.2518162891 *cos(    31.6391869657 *T+2.465042647 )+ \
           0.2866695261 *cos(    28.3859318865 *T+1.942176992 )- \
           0.1734781333 *cos(    24.2728603974 *T-3.661486981 )+ \
           3.9257210038 *cos(   523.0807466803 *T+4.920937029 )+ \
           9.7960492288 *cos(  1176.9853693166 *T+1.229392026 )- \
           4.6863263701 *cos(   625.6777530192 *T-3.658444681 )+ \
         140.4135729258 *cos( 14985.4400134205 *T+0.673880395 )

    dt6 = \
           9.5603997187 *cos(  1216.8002696575 *T+0.525733528 )+ \
           4.3565597837 *cos(   620.6809778716 *T+0.837688810 )+ \
           0.5787112035 *cos(    95.5599741609 *T+4.194535082 )+ \
          10.8374154120 *cos(  1336.7972631107 *T+3.793235253 )+ \
           3.8032116197 *cos(   565.0292110678 *T+5.639906583 )+ \
           2.6510010126 *cos(   416.4311989613 *T+2.262081818 )+ \
           3.5775308197 *cos(   521.6580372801 *T+0.642063318 )+ \
           4.6227394761 *cos(   668.1224853400 *T+6.018501522 )+ \
           5.2102470690 *cos(   763.2943259650 *T+3.458654112 )- \
           0.8453499640 *cos(   159.2596013633 *T-2.500382359 )+ \
           5.7950207099 *cos(  1137.1704689758 *T+2.547107806 )+ \
           2.5821411067 *cos(   533.3900241022 *T+0.437078094 )+ \
           3.3306029978 *cos(   596.6683980335 *T+2.246174308 )+ \
           7.5183107824 *cos(  1192.6254413669 *T+2.512929171 )+ \
          15.5707047745 *cos(  2358.1258177318 *T+5.393136889 )+ \
           0.7386010817 *cos(   158.9072895284 *T+1.275847090 )+ \
           3.2958662301 *cos(   643.8496249426 *T+1.486539246 )

    dt7 = \
           1.9405627696 *cos(   429.2330832950 *T+6.140635794 )+ \
          13.0716904585 *cos(  2301.3539539587 *T+4.557814849 )+ \
           3.0335492316 *cos(   723.4794256242 *T+4.869091389 )+ \
           2.7087124064 *cos(   723.8675591600 *T+4.691976180 )+ \
           5.2357934782 *cos(  1149.9656222793 *T+5.554998314 )+ \
           2.6597769831 *cos(   543.6993015240 *T+1.475415597 )+ \
           1.9136331857 *cos(   473.2030627343 *T+1.398784824 )+ \
           5.2014065102 *cos(  1249.1370101415 *T+5.650931916 )+ \
           5.0073878545 *cos(  1151.3883316794 *T+2.181745369 )+ \
           4.9097305146 *cos(  1252.8018664345 *T+5.823319737 )+ \
           2.1391862996 *cos(   683.6645252834 *T+0.003844094 )- \
           2.8799081722 *cos(   705.8598461315 *T-3.690360123 )+ \
           1.8560404896 *cos(   628.3143160294 *T+4.447203799 )+ \
           0.2120060235 *cos(    73.5876513532 *T+0.349250250 )+ \
           1.9180064924 *cos(   584.9364112115 *T+4.893384368 )+ \
           2.2510447875 *cos(   620.9778724132 *T+1.473760578 )

    dt8 = \
           0.2917765822 *cos(    94.9175608970 *T+5.185878737 )+ \
           2.7521608827 *cos(   991.7696874510 *T+1.030026325 )+ \
           2.9036028346 *cos(  1097.3555686350 *T+3.918259169 )+ \
           6.4715681255 *cos(  2513.2303399966 *T+6.109659023 )+ \
           5.0195530243 *cos(  1831.9536584880 *T+4.320519510 )- \
          30.5064920384 *cos(  9095.5551694697 *T-1.942656623 )+ \
           1.4915676875 *cos(   649.6374945429 *T+5.061810696 )+ \
           1.8530954325 *cos(   617.2869528772 *T+2.797822767 )+ \
           8.8091720131 *cos(  2751.1467873537 *T+0.531673101 )- \
           1.8560007226 *cos(   628.3008539689 *T-4.533471191 )+ \
           0.1505678315 *cos(    63.9897286314 *T+3.734548088 )+ \
           3.8898055312 *cos(  1620.0772724501 *T+2.605547070 )+ \
          71.1780432875 *cos( 23314.1314403759 *T+3.029030662 )+ \
          25.1859628751 *cos(  8328.6914269554 *T+2.355556099 )+ \
           4.9524696003 *cos(  1729.8182327326 *T+5.240963796 )- \
           1.4887923221 *cos(   707.9373856808 *T-5.756641637 )+ \
          19.3444739373 *cos(  8399.6847317911 *T+2.013686814 )

    dt9 = \
           4.1623742474 *cos(  1807.3704938650 *T+1.089100410 )+ \
           1.5920718380 *cos(   638.6168624210 *T+0.645026535 )+ \
           2.3891765949 *cos(  1101.5106477335 *T+4.845297676 )+ \
           1.4965569645 *cos(   624.3458341645 *T+3.809290043 )+ \
           0.2537682280 *cos(   116.2474704408 *T+6.179611691 )+ \
           1.4697843406 *cos(   624.6427287062 *T+4.781719760 )+ \
           0.6733040383 *cos(   389.4181829542 *T+1.264976635 )- \
           0.5931425099 *cos(   312.8388765096 *T-4.914231596 )+ \
           2.9777729844 *cos(  1471.2317116458 *T+2.752035928 )+ \
           1.0926058983 *cos(   629.0189396992 *T+5.280820144 )+ \
           0.1095682625 *cos(    49.1557929457 *T+1.571007057 )+ \
           2.2931297317 *cos(  1431.4168113050 *T+4.203664806 )+ \
           4.2651860746 *cos(  2248.3848574493 *T+4.167932508 )- \
           0.6823239610 *cos(   373.8761430108 *T-0.545828785 )+ \
           0.1992996318 *cos(   105.2268383188 *T+5.817167450 )

    dt10 = \
           1.5465742831 *cos(  1098.4192351700 *T+2.732084787 )+ \
           2.0084252198 *cos(  1087.3986030480 *T+2.903477885 )- \
           1.2012595327 *cos(   863.5942003763 *T-0.593891500 )- \
          16.0836703298 *cos(  8886.0057071188 *T-0.487355242 )- \
           0.2564079581 *cos(   199.0745017041 *T-3.913022880 )+ \
           3.0535569964 *cos(  2354.3230504682 *T+3.063805171 )+ \
           1.5094672194 *cos(  1096.9965257698 *T+5.152914309 )+ \
          42.6281862012 *cos( 24428.7600007027 *T+3.626395673 )+ \
           5.1847326312 *cos(  3144.1677569757 *T+1.952049260 )+ \
           1.3063363610 *cos(   922.5539273283 *T+4.996408389 )+ \
           0.5947611083 *cos(   480.4209275927 *T+5.503379738 )+ \
           0.6757819785 *cos(   459.0910180489 *T+4.164913291 )+ \
           0.7061165930 *cos(   604.0347246017 *T+5.841719038 )+ \
           0.5756149135 *cos(   554.0085789459 *T+2.769753519 )+ \
           1.3578822180 *cos(  1057.5406682942 *T+5.306538209 )

    dt11 = \
           2.4054525278 *cos(  1820.9330263660 *T+2.624866359 )+ \
           2.7533224455 *cos(  2122.8392023546 *T+0.382603541 )+ \
           0.5993119135 *cos(   628.2095528923 *T+0.882213514 )+ \
           0.6937247057 *cos(   605.8731054289 *T+1.169483931 )+ \
           0.5430708150 *cos(   554.7199336460 *T+5.448375984 )- \
           0.6180890549 *cos(   626.2300454499 *T-2.656486959 )- \
          16.5547842580 *cos( 15471.7609887482 *T-1.827624012 )+ \
           0.4658806453 *cos(   470.1116501708 *T+4.387001801 )+ \
           1.2385257008 *cos(  1391.6019109642 *T+5.601498297 )- \
           0.1371123977 *cos(   155.1045222648 *T-1.088831705 )+ \
           0.4395337333 *cos(   501.7508371365 *T+3.969902609 )+ \
           1.2179035693 *cos(  1511.0466119866 *T+5.142876744 )- \
           0.3197831765 *cos(   413.6910433516 *T-0.022067765 )- \
           0.5994989587 *cos(   628.4056171060 *T-0.968480906 )+ \
           0.3925841793 *cos(   532.6786694021 *T+4.923831588 )+ \
           0.7094374825 *cos(   866.2240323563 *T+5.991247817 )

    dt12 = \
           1.2476544601 *cos(  1725.6631536341 *T+6.068719637 )+ \
           0.5675009552 *cos(   603.7244203762 *T+6.197428148 )+ \
           1.0365965457 *cos(  1171.2955318231 *T+3.280414875 )+ \
           0.8866680936 *cos(  1255.9038152982 *T+2.824848947 )+ \
           0.1741548395 *cos(   237.9164473572 *T+2.501813417 )- \
           0.4681528764 *cos(   612.7655450557 *T-2.236346329 )+ \
           3.2082301750 *cos(  3537.1887265976 *T+3.370195967 )+ \
           1.4370991359 *cos(  1765.4780539750 *T+4.627122566 )+ \
           0.7366775339 *cos(  1140.3676995575 *T+3.852959484 )+ \
           0.9519313015 *cos(  1352.1751441591 *T+2.300991267 )- \
           3.9284014720 *cos(  6514.7619767937 *T-4.140083146 )

    dt13 = \
           0.6197949927 *cos(  1017.7257679534 *T+0.437122327 )+ \
           0.3639363159 *cos(   576.7611978898 *T+4.026532329 )+ \
           0.6386276233 *cos(  1108.7285125918 *T+4.760293101 )+ \
           1.0073143101 *cos(  1494.5316173554 *T+6.270510511 )+ \
           0.3942092494 *cos(   542.9879468239 *T+6.039606892 )+ \
           2.0424516341 *cos(  2876.6924424484 *T+5.672617711 )+ \
           0.7670973468 *cos(  1185.6218651625 *T+3.397132627 )- \
           0.3716290835 *cos(   548.1254918868 *T-6.249666675 )+ \
           1.3598419244 *cos(  2200.3914634870 *T+2.466427018 )+ \
           0.4527627879 *cos(   613.4997125565 *T+2.242668890 )+ \
           0.2426246781 *cos(   349.6032826134 *T+2.668309141 )+ \
           0.3445797803 *cos(   648.9261398429 *T+1.681888780 )- \
           8.7722079138 *cos( 14357.1324284214 *T-2.424978312 )+ \
           0.6925055212 *cos(  1204.3574281889 *T+4.216492400 )+ \
           0.6866373442 *cos(  1241.6588502848 *T+4.772158039 )+ \
           0.3229266802 *cos(   468.6889407707 *T+6.224271088 )+ \
           0.3634516601 *cos(   734.2457780181 *T+3.817285811 )+ \
           0.2060830121 *cos(   363.4621024518 *T+1.649264690 )+ \
           0.9597503154 *cos(  1863.5928454536 *T+3.945345892 )

    dt14 = \
           1.6654994338 *cos(  2515.8601719765 *T+1.794058369 )- \
           0.5932886514 *cos(  1256.9674818332 *T-5.112133338 )+ \
           0.2848972198 *cos(   617.9983075773 *T+0.513669325 )+ \
           5.3502447222 *cos(  8346.7156352816 *T+3.210727723 )+ \
           0.5379033434 *cos(  1034.4295065386 *T+2.445597761 )+ \
           0.9082356274 *cos(  1842.2629359098 *T+1.676939306 )+ \
           0.6255421246 *cos(  1119.0377900137 *T+5.783236356 )+ \
           0.4754102009 *cos(   962.3688276691 *T+3.022645053 )+ \
           0.2657230057 *cos(   573.9157790895 *T+1.411223013 )+ \
           0.7282864454 *cos(  1685.8482532933 *T+1.179256434 )+ \
           4.1408720854 *cos(  7214.0628666286 *T+1.758191830 )+ \
           0.8357357810 *cos(  1726.7268201691 *T+3.290589143 )+ \
           0.2699016127 *cos(   490.7302050146 *T+0.864024298 )+ \
           0.8824278712 *cos(  2042.6571092422 *T+6.003829241 )

    dt15 = \
          0.27654867639 *cos(   574.9452731634 *T+4.309591964 )+ \
          0.27635121783 *cos(   575.7317038160 *T+1.142348571 )+ \
          0.32507418395 *cos(   670.2560493867 *T+0.210580917 )+ \
          0.25796641554 *cos(   605.5549660552 *T+4.274476529 )+ \
          0.28605938080 *cos(   595.9570433334 *T+5.031351030 )+ \
          0.58541849190 *cos(  1256.2628581634 *T+4.959581597 )+ \
          2.04370904203 *cos(  3930.2096962196 *T+4.788002889 )+ \
          0.55566575026 *cos(  1213.2439962106 *T+1.880103788 )+ \
          0.56537931783 *cos(  1202.9347187887 *T+1.405611197 )- \
          0.31106495098 *cos(   747.7522860216 *T-1.082356330 )+ \
          0.52128282823 *cos(  1160.9862544012 *T+4.179989585 )+ \
          0.80226641151 *cos(  1725.3041107690 *T+0.353496295 )- \
          0.16416915162 *cos(   453.5059436924 *T-1.583849576 )+ \
          0.84084423644 *cos(  2195.4157609398 *T+3.747376371 )+ \
          0.59755196788 *cos(  1805.2929543158 *T+0.566790582 )+ \
          0.58126841457 *cos(  1351.7870106233 *T+0.685827538 )- \
          0.21185421452 *cos(   575.6908003246 *T-0.731374317 )+ \
          0.34840060731 *cos(  1055.7594160824 *T+3.710043680 )+ \
          0.67060995266 *cos(  2019.9094959633 *T+1.652901407 )

    dt16 = \
          0.45824132969 *cos(  1193.3367960670 *T+5.827781531 )+ \
          0.40458920366 *cos(  1045.4501386605 *T+2.541182564 )+ \
          0.50931015718 *cos(  1567.1081759407 *T+2.178850542 )+ \
          0.28633418024 *cos(   938.8005909415 *T+0.578340206 )+ \
          0.20239513418 *cos(   574.9861766548 *T+3.000297967 )+ \
          0.21508323323 *cos(   691.5859589305 *T+1.693574249 )+ \
          0.71496576765 *cos(  2407.2921469776 *T+1.997249392 )+ \
          0.40673613564 *cos(  1259.2450019783 *T+1.072262823 )+ \
          0.41420134661 *cos(  1214.6667056108 *T+4.700657997 )+ \
          0.28359415161 *cos(   977.9108676125 *T+1.812320441 )+ \
          0.20971536376 *cos(   613.2028180148 *T+4.322238614 )+ \
          0.20624512407 *cos(   626.8848755990 *T+3.033827743 )+ \
          0.67305156569 *cos(  1799.6031168222 *T+3.388716544 )+ \
          0.20502554752 *cos(   606.5844601290 *T+0.877776108 )

    dt17 = \
          0.35489733622 *cos(  1309.5842665077 *T+3.208912203 )+ \
          0.20103976485 *cos(   607.3708907816 *T+4.007881169 )+ \
          1.05301497326 *cos(  2908.8811415985 *T+3.215977013 )+ \
          0.34607905224 *cos(  1235.9966151546 *T+0.710872502 )+ \
          0.27875532361 *cos(  1044.0274292604 *T+4.730108488 )+ \
          0.41240903491 *cos(  1649.6361396202 *T+0.898769761 )+ \
          0.66941042880 *cos(  2059.7243963041 *T+0.180044365 )+ \
          0.16476668863 *cos(   614.8010769956 *T+5.152666276 )+ \
          0.16006424647 *cos(   563.6065016677 *T+5.655385808 )+ \
          0.18303275589 *cos(   608.0822454817 *T+2.135396205 )+ \
          0.13730004661 *cos(   586.7523359379 *T+5.575209112 )- \
          6.07980079324 *cos( 22685.8238553767 *T-0.069432392 )+ \
          4.43301968207 *cos( 16728.3761587465 *T+4.369302826 )+ \
          0.79064253686 *cos(  2823.7233459389 *T+5.304829118 )+ \
          0.36049558048 *cos(  1234.5739057544 *T+4.096094132 )

    dt18 = \
          0.44156109482 *cos(  1980.0945956225 *T+3.069327406 )+ \
          1.30129243042 *cos(  4323.2306658416 *T+6.205311188 )+ \
          0.49831388296 *cos(  1887.5525869774 *T+1.417263408 )+ \
          9.46120403733 *cos( 31642.8228673312 *T+5.384595078 )+ \
          0.12146354848 *cos(   575.6566278634 *T+3.789392838 )+ \
          0.12017925296 *cos(   575.0203491159 *T+1.661943545 )+ \
          0.29975725509 *cos(  1248.9885628707 *T+5.684549045 )+ \
          0.13616318690 *cos(   630.3851245484 *T+3.862942261 )+ \
          0.11284396485 *cos(   564.2198242609 *T+1.016115785 )+ \
          0.14271508167 *cos(   628.7008003254 *T+2.911891613 )- \
          0.12872945114 *cos(   627.9485421340 *T-1.829362730 )- \
          0.22966609049 *cos(  1098.8808157535 *T-2.636140084 )

    dt19 = \
          0.55663569240 *cos(  2929.6615389579 *T+4.175658539 )+ \
          1.75738162683 *cos(  6656.7485864652 *T+4.601102551 )+ \
          0.24665075673 *cos(  1311.9721102825 *T+2.032195842 )+ \
          0.51349566176 *cos(  2593.4124331089 *T+3.832703118 )+ \
          0.11965560981 *cos(   511.3487598583 *T+1.716090661 )+ \
          0.26247418840 *cos(  1182.3161639450 *T+1.942386641 )+ \
          0.19495316794 *cos(  1077.0893256262 *T+1.999482059 )+ \
          0.11193933212 *cos(   654.6159773364 *T+1.182807992 )+ \
          0.35482214113 *cos(  2099.5392966449 *T+2.169080622 )+ \
          0.20361912046 *cos(  1066.0686935042 *T+5.405515999 )+ \
          0.75283368136 *cos(  3301.9021112205 *T+4.656985514 )+ \
          0.38414986916 *cos(  2314.1558382925 *T+3.454132746 )+ \
          0.11736106454 *cos(   614.4558353121 *T+5.020393445 )

    dt20 = \
          0.10951206927 *cos(   608.4003848555 *T+0.602182191 )+ \
          0.28985853279 *cos(  1778.2732072784 *T+4.960593133 )+ \
          0.37035750441 *cos(  1646.0333529525 *T+2.596451817 )+ \
          0.13110658977 *cos(   590.5702242076 *T+3.731990323 )+ \
          0.26612747806 *cos(  1673.7577236597 *T+3.600691544 )+ \
          0.13611306536 *cos(   680.5653268085 *T+0.868220961 )+ \
          0.22288793421 *cos(  1191.9140866668 *T+2.629456641 )+ \
          0.12887665871 *cos(   628.6666278643 *T+1.742882331 )+ \
          0.28095228240 *cos(  1672.3350142595 *T+0.027860588 )+ \
          0.17462902609 *cos(  1172.0068865232 *T+0.659721876 )+ \
          0.10132495193 *cos(   625.4626662524 *T+1.435132069 )+ \
          0.88193905583 *cos(  4716.2516354635 *T+1.354371923 )

    dt21 = \
          0.16177050505 *cos(  1108.0171578918 *T+3.369695406 )+ \
          0.14681511877 *cos(  1129.3470674356 *T+0.939039445 )+ \
          0.25831509965 *cos(  1845.1078546566 *T+0.642049130 )+ \
          0.14238736119 *cos(  1130.0584221356 *T+3.485280663 )+ \
          0.11933204803 *cos(  1002.7903195729 *T+3.217431161 )- \
          0.12653251751 *cos(   766.8637425143 *T-4.298212528 )

    dt22 = \
          0.15687423221*cos(    968.3594581116 *T+5.720092446 )+ \
          0.26837223068*cos(   2274.3409379516 *T+4.881123092 )+ \
          2.57691964970*cos(  16309.6180360983 *T+2.957128968 )+ \
          0.16538021252*cos(   1234.1806904281 *T+2.598576764 )+ \
          0.25107330082*cos(   1662.7370915377 *T+3.985702050 )

    dt23 = \
          0.32720193267*cos(   2353.9707386333 *T+2.028195445 )+ \
          0.17186537443*cos(   1606.2184526117 *T+4.066520001 )+ \
          1.18972638476*cos(   8378.3548222473 *T+2.936315115 )+ \
          0.12007628381*cos(    938.0959672717 *T+3.223844306 )+ \
          0.23865440252*cos(   1940.2796952817 *T+4.538074405 )+ \
          0.27138825324*cos(   2280.5735565994 *T+2.869040566 )+ \
          0.85747748351*cos(   6447.1991241142 *T+6.056405489 )+ \
          0.17834083603*cos(   1361.3804277336 *T+4.005732868 )

    dt24 = \
          0.57734780437  *cos( 5109.2726050855 *T+2.791483066 )- \
          0.13543041651  *cos( 1253.9853380183 *T-3.716133846 )+ \
          0.44708569278  *cos( 3694.9230808424 *T+6.072332087 )+ \
          0.18029401544  *cos( 1609.7679950283 *T+3.589026260 )+ \
          0.10207188264  *cos(  981.4604100291 *T+1.959967212 )+ \
         10.2156724      *sin(  628.3075849991 *T+4.249032005 )+ \
          0.1706807      *sin( 1256.6151699983 *T+4.205904248 )
    # -----------------------------------------------------------------------
    # *T
    dt24s = \
        214.47995684593  *cos( 1256.6151699983 *T+4.205904248 )+ \
       6418.58445478596  *cos(  628.3075849991 *T+4.249032005 )+ \
          4.07427230969  *cos(  522.3693919802 *T+4.670344204 )+ \
          1.40854891747  *cos(   52.9690965095 *T+5.836047367 )+ \
          0.57519940469  *cos(   21.3299095438 *T+3.400290479 )

    dt25 = \
          1.85026251053  *cos(  550.7553238667 *T+5.980162321 )+ \
          6.04834013624  *cos( 1884.9227549974 *T+4.162913471 )+ \
          1.44748575002  *cos(  574.6271337896 *T+2.980330535 )+ \
          0.47611271970  *cos(  469.4002954708 *T+4.044013795 )+ \
          1.70997434691  *cos(  585.6477659115 *T+0.623811863 )- \
          0.13901242808  *cos(   77.5522611324 *T-1.092065955 )+ \
          1.43874208835  *cos(  576.0498431898 *T+2.467913690 )- \
          0.18312465263  *cos(   79.6298006816 *T-1.174411803 )+ \
          0.78203501291  *cos(  607.6890301554 *T+5.333425680 )+ \
          1.41717288693  *cos( 1203.6460734888 *T+2.292832062 )+ \
          0.80487916743  *cos(  606.2663207553 *T+5.845801920 )+ \
          0.51339780279  *cos(  548.6777843175 *T+3.416081409 )+ \
          0.12873609216  *cos(  105.9381930189 *T+6.222874454 )- \
          0.13704288697  *cos(   39.8149003408 *T-5.980077351 )+ \
          0.86381641759  *cos(  157.7343542448 *T+4.534800170 )

    dt26 = \
          0.94260175696  *cos( 1097.7078804699 *T+2.777152598 )+ \
          0.16461088744  *cos(  214.6165416475 *T+3.000200440 )+ \
          0.19708259496  *cos(  254.4314419883 *T+1.603197066 )+ \
          0.54149002750  *cos(  627.5962302991 *T+4.562060226 )+ \
          0.10643671939  *cos(  174.8016413067 *T+4.403765209 )- \
          0.31087231898  *cos(  628.6598968340 *T-0.268305170 )+ \
          0.29605642590  *cos(  508.8628839767 *T+4.827723531 )+ \
          0.20261168508  *cos(  681.2766815086 *T+1.583012668 )+ \
          0.17206051440  *cos(  564.3178563677 *T+1.362634430 )

    dt27 = \
          0.28860679054  *cos(  943.7762934887 *T+4.226420633 )- \
          1.28789544199  *cos( 7143.0695617928 *T-2.192295512 )+ \
          0.15170733284  *cos(  803.1092263058 *T+3.569003717 )+ \
          0.14024534710  *cos(  469.0479836359 *T+6.235872050 )+ \
          0.16945778585  *cos(  586.3591206116 *T+0.095197563 )+ \
          0.16083363631  *cos(  643.8496249426 *T+2.994779800 )+ \
          0.25543863267  *cos( 1044.7387839604 *T+2.347139160 )+ \
          0.15268930751  *cos(  521.6580372801 *T+2.319951253 )+ \
          0.12684708106  *cos(  523.0807466803 *T+3.084752833 )+ \
          0.14488613917  *cos(  708.4896781115 *T+0.526323854 )+ \
          0.10148646734  *cos(  584.9364112115 *T+0.417558428 )+ \
          0.13849243400  *cos(  842.9241266467 *T+2.180619584 )+ \
          0.31981535393  *cos( 1235.2852604545 *T+1.991935820 )

    dt28 = \
          0.19334157996  *cos( 1414.3495242431 *T+3.987576591 )+ \
          0.16536505648  *cos( 1149.9656222793 *T+0.974387904 )+ \
          0.12426431626  *cos(  763.2943259650 *T+4.968445721 )+ \
          0.15635853544  *cos( 1151.3883316794 *T+0.495572260 )+ \
          0.12598039158  *cos( 1673.0463689596 *T+3.905030235 )+ \
          0.17932164384  *cos( 1778.9845619785 *T+0.286350174 )+ \
          0.10913862647  *cos(  627.9552731642 *T+0.087484036 )

    dt29 = 0.15959012659  *cos( 2513.2303399966 *T+4.121051532 )

    dt32 = 0.1192755232   *cos( 8399.6847317911 *T+0.022670115 )
    # ----------------------------------------------------------------------
    # *T**2  
    dt34 = 27.16167407    *cos(  628.3075849991 *T+2.642893748 )

    # Correction to use JPL planetary masses instead of IAU.
    dt41 = 0.394535489 * cos( 606.9776754 *T + 4.021194) - \
         1.216825673 * cos( 620.8294251 *T + 5.696701) + \
         0.0007276   * T
    # Other terms are very small

    dperdct =  dt34*T_2 + dt41+ \
               (dt32+dt29+dt28+dt27+dt26+dt25+dt24s)*T+ \
               dt24+dt23+dt22+dt21+dt20+dt19+dt18+ \
               dt17+dt16+dt15+dt14+dt13+dt12+dt11+ \
               dt10+dt9 +dt8 +dt7 +dt6 +dt5 +dt4 +dt3 +dt2+dt1

    # Conversion of derivative to [Sec/Sec]
    dperdct = dperdct * 1e-6/86400.0/36525.0
    return per_terms, dperdct
    
    
#==============================================================================
#     
#==============================================================================
def ter2cel(date_utc, eop_int, dTAIdCT=1.0, theory='iau2000', mode='der'):
    ''' 
    Compute celestial-to-terrestrial matrix for a given UTC date/time
    
    mode='der' - do calculate derivatives
    mode='noDer' - do not calculate derivatives
    '''
    if mode=='der':
        r2000 = np.zeros((3,3,3))
    elif mode=='noDer':
        r2000 = np.zeros((3,3))
    sec = datetime.timedelta(seconds=1)
    # following IAU 2000A:
    if theory=='iau2000':
        cel2ter = cel2ter00(date_utc, eop_int)
        # numerically calculate derivative:
        if mode=='der':
#            print date_utc - 2*sec
            cel2ter_m2 = cel2ter00(date_utc - 2*sec, eop_int) # -2s TAI
            cel2ter_m1 = cel2ter00(date_utc -   sec, eop_int) # -1s TAI
            cel2ter_p1 = cel2ter00(date_utc +   sec, eop_int) # +1s TAI
            cel2ter_p2 = cel2ter00(date_utc + 2*sec, eop_int) # +2s TAI
    # following IAU 2006/2000A:   
    if theory=='iau2006':
        raise NotImplemented
#        cel2ter = cel2ter06(date_utc, eop_int);
#        if mode=='der':
#            cel2ter_m2 = cel2ter06(date_utc - 2*sec, eop_int) # -2s TAI
#            cel2ter_m1 = cel2ter06(date_utc -   sec, eop_int) # -1s TAI
#            cel2ter_p1 = cel2ter06(date_utc +   sec, eop_int) # +1s TAI
#            cel2ter_p2 = cel2ter06(date_utc + 2*sec, eop_int) # +2s TAI

    if mode=='der':
    # transpose it to get terrestrial-to-celestial matrix and its CT derivative
        r2000[:,:,0] = cel2ter.T
    # (dr/dTAI)*(dTAI/dCT):
        r2000[:,:,1] = dTAIdCT * (cel2ter_p1.T - cel2ter_m1.T) / 2.0 
    # (dv/dTAI)*(dTAI/dCT):
        r2000[:,:,2] = dTAIdCT * (cel2ter_p2.T - 
                       2.0*cel2ter.T + cel2ter_m2.T) / 4.0
#        print '{:.18e}'.format(cel2ter_m2[0,0])
#        print '{:.18e}'.format(-2.0*cel2ter[0,0])
#        print '{:.18e}'.format(cel2ter_p2[0,0])
#        print '{:.18e}'.format(-np.float64(cel2ter[0,0])-np.float64(cel2ter[0,0])+\
#        np.float64(cel2ter_m2[0,0])+np.float64(cel2ter_p2[0,0]))
#        print cel2ter_p2[0,0]
#        print cel2ter_m2.T
#        print cel2ter_m2.T - 2.0*cel2ter.T + cel2ter_p2.T
#        print cel2ter_p2.T - 2.0*cel2ter.T + cel2ter_m2.T
#        print r2000[:,:,2]
    elif mode=='noDer':
        r2000 = cel2ter.T
    else:
        raise Exception('unrecognised mode')
    
    return r2000
    
    
#==============================================================================
#     
#==============================================================================
#def cel2ter00(date_utc, eop_int):
#    date_eop = (date_utc.year, date_utc.month, date_utc.day,\
#                date_utc.hour, date_utc.minute, date_utc.second,\
#                eop_int[0],eop_int[1],eop_int[2],eop_int[3],eop_int[4])
#    return cel2ter00_mem(date_eop)
#
#@memoize
#def cel2ter00_mem(date_eop):
def cel2ter00(date_utc, eop_int):
    '''
    # Function vint_cel2ter calculates transformation matrix ICRF -> ITRF
    # according to the IERS 2010 Conventions. IAU 2000A CIO-based remel 
    # with classiacal angles is used.
    #
    # input: date_utc = [2007 04 05 12 00 0.0]  : array with the date in UTC
    #    eop_int[0:5]:
    #        dut1 : UT1-UTC in seconds
    #        xp, yp : celestial pole coordinates in arcseconds
    #        dx06, dy06 : corrections to the CIP position, IAU2000A in rad
    #
    # output: rc2it : transformation matrix from ICRF to ITRF
    #
    # matlab version by D. Duev (JIVE) November 2011
    '''
    ## initialize
#    date_utc = datetime.datetime(*date_eop[0:6])
#    eop_int = date_eop[6:]
    
    # * Arcseconds to radians
    as2r = 4.848136811095359935899141e-6
    
    # UT1-UTC
    dut1 = eop_int[0]
    # * Polar motion (arcsec->radians).
    xp = eop_int[1] * as2r
    yp = eop_int[2] * as2r
    # * CIP offsets wrt IAU 2000A (as->radians)
    dx00 = eop_int[3] * as2r
    dy00 = eop_int[4] * as2r
    
    # date
    iy = date_utc.year
    im = date_utc.month
    iday = date_utc.day
    ih = date_utc.hour
    mn = date_utc.minute
    sc = date_utc.second + date_utc.microsecond*1e-6
    ## IAU 2000A remel
    
    # * TT (MJD).
    date = mjuliandate ( iy, im, iday )
    DJMJD0 =  2400000.5
    time = ( 60.0*(60.0*float(ih) + float(mn)) + sc ) / 86400.0
    UTC = date + time
#    DAT = iau_DAT ( iy, im, iday, time )
    DAT = nsec(date)
    TAI = UTC + DAT/86400.0
    TT = TAI + 32.184/86400.0
    # * UT1.
    TUT = time + dut1/86400.0
    # UT1 = date + TUT
    
    
    ## * =========================
    ## * IAU 2006/2000A, CIO based
    ## * =========================
    # * CIP and CIO, IAU 2006/2000A.
    try:
        X, Y, S = iau_xys00a_fort ( DJMJD0, TT ) # fortran version
    except:
        print 'f2pyed version of iau_xys00a did not work,\
               so it is gonna take a wwhhhiiillle. u betah fix dat!'
        # this is 2 orders of magnitude slower!!
        X, Y, S = iau_XYS00A ( DJMJD0, TT ) # python version

    # * Add CIP corrections.
    X = X + dx00
    Y = Y + dy00
    # * GCRS to CIRS matrix.
    RC2I = iau_C2IXYS ( X, Y, S )
    # * Earth rotation angle.
    ERA = iau_ERA00 ( DJMJD0+date, TUT )

    # * Form celestial-terrestrial matrix (no polar motion yet).
    RC2TI = iau_RZ ( ERA, RC2I )

    # * Polar motion matrix (TIRS->ITRS, IERS 2003).
    rpom = iau_POM00 ( xp, yp, iau_SP00(DJMJD0,TT) )
#    print '------'
#    print '{:.18e} {:.18e} {:.18e}'.format(X, Y, S)
#    print '------'

    # * Form celestial-terrestrial matrix (including polar motion).
    rc2it = dot(rpom, RC2TI)
#    print '{:.18e}'.format(rc2it[0,0])
#    print rc2it    
    return rc2it


def iau_ANP ( A ):
    '''
    # *+
    # *  - - - - - - - -
    # *   i a u _ A N P
    # *  - - - - - - - -
    # *
    # *  Normalize angle into the range 0 <= A < 2pi.
    # *
    # *  This routine is part of the International Astronomical Union's
    # *  SOFA (Standards of Fundamental Astronomy) software collection.
    # *
    # *  Status:  vector/matrix support routine.
    # *
    # *  Given:
    # *     A          d       angle (radians)
    # *
    # *  Returned:
    # *     iau_ANP    d       angle in range 0-2pi
    # *
    # *  This revision:  2000 December 15
    # *
    # *  SOFA release 2010-12-01
    # *
    # *  Copyright (C) 2010 IAU SOFA Board.  See notes at end.
    # *
    # *-----------------------------------------------------------------------
    '''

    # *  2Pi
    D2PI = 6.283185307179586476925287
    W = fmod(A, D2PI)
    if ( W < 0.0 ): W = W + D2PI
    return W

    # *  Finished.

def iau_C2IXYS ( X, Y, S ):
    '''
    # *+
    # *  - - - - - - - - - - -
    # *   i a u _ C 2 I X Y S
    # *  - - - - - - - - - - -
    # *
    # *  Form the celestial to intermediate-frame-of-date matrix given the CIP
    # *  X,Y and the CIO locator s.
    # *
    # *  This routine is part of the International Astronomical Union's
    # *  SOFA (Standards of Fundamental Astronomy) software collection.
    # *
    # *  Status:  support routine.
    # *
    # *  Given:
    # *     X,Y        d       Celestial Intermediate Pole (Note 1)
    # *     S          d       the CIO locator s (Note 2)
    # *
    # *  Returned:
    # *     RC2I     d(3,3)    celestial-to-intermediate matrix (Note 3)
    # *
    # *  Notes:
    # *
    # *  1) The Celestial Intermediate Pole coordinates are the x,y components
    # *     of the unit vector in the Geocentric Celestial Reference System.
    # *
    # *  2) The CIO locator s (in radians) positions the Celestial
    # *     Intermediate Origin on the equator of the CIP.
    # *
    # *  3) The matrix RC2I is the first stage in the transformation from
    # *     celestial to terrestrial coordinates:
    # *
    # *        [TRS]  =  RPOM * R_3(ERA) * RC2I * [CRS]
    # *
    # *               =  RC2T * [CRS]
    # *
    # *     where [CRS] is a vector in the Geocentric Celestial Reference
    # *     System and [TRS] is a vector in the International Terrestrial
    # *     Reference System (see IERS Conventions 2003), ERA is the Earth
    # *     Rotation Angle and RPOM is the polar motion matrix.
    # *
    # *  Called:
    # *     iau_IR       initialize r-matrix to identity
    # *     iau_RZ       rotate around Z-axis
    # *     iau_RY       rotate around Y-axis
    # *
    # *  Reference:
    # *
    # *     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    # *     IERS Technical Note No. 32, BKG (2004)
    # *
    # *  This revision:  2006 October 10
    # *
    # *  SOFA release 2010-12-01
    # *
    # *  Copyright (C) 2010 IAU SOFA Board.  See notes at end.
    # *
    # *-----------------------------------------------------------------------
    '''

    # *  Obtain the spherical angles E and d.
    R2 = X*X+Y*Y
#    print 'Y = {:.18e}'.format(Y)
#    print 'X = {:.18e}'.format(X)
    if ( R2!=0.0 ):
        E = atan2 ( Y, X )
    else:
        E = 0.0
    D = atan ( sqrt ( R2 / (1.0-R2) ) )
#    print 'E = {:.18e}'.format(E)
#    print 'D = {:.18e}'.format(D)
    # *  Form the matrix.
    RC2I = np.identity(3)
    RC2I = iau_RZ ( E, RC2I )
    RC2I = iau_RY ( D, RC2I )
    RC2I = iau_RZ ( -(E+S), RC2I )
    
    return RC2I
    # *  Finished.


def iau_DAT ( IY, IM, ID, FD ):
    '''
    # *+
    # *  - - - - - - - -
    # *   i a u _ D A T
    # *  - - - - - - - -
    # *
    # *  For a given UTC date, calculate delta(AT) = TAI-UTC.
    # *
    # *     :------------------------------------------:
    # *     :                                          :
    # *     :                 IMPORTANT                :
    # *     :                                          :
    # *     :  A new version of this routine must be   :
    # *     :  produced whenever a new leap second is  :
    # *     :  announced.  There are five items to     :
    # *     :  change on each such occasion:           :
    # *     :                                          :
    # *     :  1) The parameter NDAT must be           :
    # *     :     increased by 1.                      :
    # *     :                                          :
    # *     :  2) A new line must be added to the set  :
    # *     :     of DATA statements that initialize   :
    # *     :     the arrays IDATE and DATS.           :
    # *     :                                          :
    # *     :  3) The parameter IYV must be set to     :
    # *     :     the current year.                    :
    # *     :                                          :
    # *     :  4) The "Latest leap second" comment     :
    # *     :     below must be set to the new leap    :
    # *     :     second date.                         :
    # *     :                                          :
    # *     :  5) The "This revision" comment, later,  :
    # *     :     must be set to the current date.     :
    # *     :                                          :
    # *     :  Change (3) must also be carried out     :
    # *     :  whenever the routine is re-issued,      :
    # *     :  even if no leap seconds have been       :
    # *     :  added.                                  :
    # *     :                                          :
    # *     :  Latest leap second:  2015 Jun 30        :
    # *     :                                          :
    # *     :__________________________________________:
    # *
    # *  This routine is part of the International Astronomical Union's
    # *  SOFA (Standards of Fundamental Astronomy) software collection.
    # *
    # *  Status:  support routine.
    # *
    # *  Given:
    # *     IY       i     UTC:  year (Notes 1 and 2)
    # *     IM       i           month (Note 2)
    # *     ID       i           day (Notes 2 and 3)
    # *     FD       d           fraction of day (Note 4)
    # *
    # *  Returned:
    # *     DELTAT   d     TAI minus UTC, seconds
    # *     J        i     status (Note 5):
    # *                       1 = dubious year (Note 1)
    # *                       0 = OK
    # *                      -1 = bad year
    # *                      -2 = bad month
    # *                      -3 = bad day (Note 3)
    # *                      -4 = bad fraction (Note 4)
    # *
    # *  Notes:
    # *
    # *  1) UTC began at 1960 January 1.0 (JD 2436934.5) and it is improper
    # *     to call the routine with an earlier date.  if this is attempted,
    # *     zero is returned together with a warning status.
    # *
    # *     Because leap seconds cannot, in principle, be predicted in
    # *     advance, a reliable check for dates beyond the valid range is
    # *     impossible.  To guard against gross errors, a year five or more
    # *     after the release year of the present routine (see parameter IYV)
    # *     is considered dubious.  In this case a warning status is returned
    # *     but the result is computed in the normal way.
    # *
    # *     For both too-early and too-late years, the warning status is J=+1.
    # *     This is distinct from the error status J=-1, which signifies a
    # *     year so early that JD could not be computed.
    # *
    # *  2) if the specified date is for a day which ends with a leap second,
    # *     the UTC-TAI value returned is for the period leading up to the
    # *     leap second.  if the date is for a day which begins as a leap
    # *     second ends, the UTC-TAI returned is for the period following the
    # *     leap second.
    # *
    # *  3) The day number must be in the normal calendar range, for example
    # *     1 through 30 for April.  The "almanac" convention of allowing
    # *     such dates as January 0 and December 32 is not supported in this
    # *     routine, in order to avoid confusion near leap seconds.
    # *
    # *  4) The fraction of day is used only for dates before the introduction
    # *     of leap seconds, the first of which occurred at the end of 1971.
    # *     It is tested for validity (zero to less than 1 is the valid range)
    # *     even if not used;  if invalid, zero is used and status J=-4 is
    # *     returned.  For many applications, setting FD to zero is
    # *     acceptable  the resulting error is always less than 3 ms (and
    # *     occurs only pre-1972).
    # *
    # *  5) The status value returned in the case where there are multiple
    # *     errors refers to the first error detected.  For example, if the
    # *     month and day are 13 and 32 respectively, J=-2 (bad month) will be
    # *     returned.
    # *
    # *  6) In cases where a valid result is not available, zero is returned.
    # *
    # *  References:
    # *
    # *  1) For dates from 1961 January 1 onwards, the expressions from the
    # *     file ftp://maia.usno.navy.mil/ser7/tai-utc.dat are used.
    # *
    # *  2) The 5ms timestep at 1961 January 1 is taken from 2.58.1 (p87) of
    # *     the 1992 Explanatory Supplement.
    # *
    # *  Called:
    # *     iau_CAL2JD   Gregorian calendar to Julian Day number
    # *
    # *  This revision:  2009 November 12
    # *
    # *  SOFA release 2010-12-01
    # *
    # *  Copyright (C) 2010 IAU SOFA Board.  See notes at end.
    # *
    # *-----------------------------------------------------------------------
    '''
    # *  Release year for this version of iau_DAT
    # IYV = 2013

    # *  Number of Delta(AT) changes (increase by 1 for each new leap second)
    NDAT = 40

    # *  Number of Delta(AT) expressions before leap seconds were introduced
    NERA1 = 14

    # *  Dates (year, month) on which new Delta(AT) came into force
    IDATE = np.asarray([[1960,  1], [1961,  1], [1961,  8], [1962,  1], [1963, 11], [1964,  1],
             [1964,  4], [1964,  9], [1965,  1], [1965,  3], [1965,  7], [1965,  9],
             [1966,  1], [1968,  2], [1972,  1], [1972,  7], [1973,  1], [1974,  1],
             [1975,  1], [1976,  1], [1977,  1], [1978,  1], [1979,  1], [1980,  1],
             [1981,  7], [1982,  7], [1983,  7], [1985,  7], [1988,  1], [1990,  1],
             [1991,  1], [1992,  7], [1993,  7], [1994,  7], [1996,  1], [1997,  7],
             [1999,  1], [2006,  1], [2009,  1], [2012,  7], [2015,  7] ])
    # *  New Delta(AT) which came into force on the given dates
    DATS = np.asarray([ 1.4178180, 1.4228180, 1.3728180,
              1.8458580, 1.9458580, 3.2401300 ,
              3.3401300, 3.4401300, 3.5401300 ,
              3.6401300, 3.7401300, 3.8401300 ,
              4.3131700, 4.2131700, 10, 11, 12, 13, 14, 15, 16,
              17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
              30, 31, 32, 33, 34, 35, 36 ])

    # *  Reference dates (MJD) and DRIFT rates (s/day), pre leap seconds
    DRIFT = [ [37300.0, 0.001296],
              [37300.0, 0.001296],
              [37300.0, 0.001296],
              [37665.0, 0.0011232],
              [37665.0, 0.0011232],
              [38761.0, 0.001296],
              [38761.0, 0.001296],
              [38761.0, 0.001296],
              [38761.0, 0.001296],
              [38761.0, 0.001296],
              [38761.0, 0.001296],
              [38761.0, 0.001296],
              [39126.0, 0.002592],
              [39126.0, 0.002592] ]

    # * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # *  if invalid fraction of a day, set error status and give up.
    if ( FD<0.0 or FD>=1.0 ):
        raise Exception('Invalid fraction of a day. KARAUL!!')

    # *  Convert the date into an MJD.
    DJM = mjuliandate ( IY, IM, ID )

    # *  Combine year and month.
    M = 12.0*IY+IM

    # *  Find the most recent table entry.
    MORE = True
    for I in range(NDAT-1,-1,-1):
        if ( MORE ):
            IS = I
            MORE = M < ( 12*IDATE[I,0] + IDATE[I,1] )

    # *  Get the Delta(AT).
    DA = DATS[IS]

    # *  if pre-1972, adjust for DRIFT.
    if ( IS <= NERA1 ):
        DA = DA + ( DJM + FD - DRIFT[IS,0] ) * DRIFT[IS,1]
      
    DELTAT = DA
    return DELTAT
    # *  Finished.


def iau_ERA00 ( DJ1, DJ2 ):
    '''
    # *+
    # *  - - - - - - - - - -
    # *   i a u _ E R A 0 0
    # *  - - - - - - - - - -
    # *
    # *  Earth rotation angle (IAU 2000 remel).
    # *
    # *  This routine is part of the International Astronomical Union's
    # *  SOFA (Standards of Fundamental Astronomy) software collection.
    # *
    # *  Status:  canonical remel.
    # *
    # *  Given:
    # *     DJ1,DJ2     d      UT1 as a 2-part Julian Date (see note)
    # *
    # *  The result is the Earth rotation angle (radians), in the range 0 to
    # *  2pi.
    # *
    # *  Notes:
    # *
    # *  1) The UT1 date DJ1+DJ2 is a Julian Date, apportioned in any
    # *     convenient way between the arguments DJ1 and DJ2.  For example,
    # *     JD(UT1)=2450123.7 could be expressed in any of these ways,
    # *     among others:
    # *
    # *             DJ1            DJ2
    # *
    # *         2450123.7.0        0.0        (JD method)
    # *          2451545.0      -1421.3.0     (J2000 method)
    # *         2400000.5.0     50123.2.0     (MJD method)
    # *         2450123.5.0       0.2.0       (date & time method)
    # *
    # *     The JD method is the most natural and convenient to use in
    # *     cases where the loss of several decimal digits of resolution
    # *     is acceptable.  The J2000 and MJD methods are good compromises
    # *     between resolution and convenience.  The date & time method is
    # *     best matched to the algorithm used:  maximum accuracy (or, at
    # *     least, minimum noise) is delivered when the DJ1 argument is for
    # *     0hrs UT1 on the day in question and the DJ2 argument lies in the
    # *     range 0 to 1, or vice versa.
    # *
    # *  2) The algorithm is adapted from Expression 22 of Capitaine et al.
    # *     2000.  The time argument has been expressed in days directly,
    # *     and, to retain precision, integer contributions have been
    # *     eliminated.  The same formulation is given in IERS Conventions
    # *     (2003), Chap. 5, Eq. 14.
    # *
    # *  Called:
    # *     iau_ANP      normalize angle into range 0 to 2pi
    # *
    # *  References:
    # *
    # *     Capitaine N., Guinot B. and McCarthy D.D, 2000, Astron.
    # *     Astrophys., 355, 398-405.
    # *
    # *     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    # *     IERS Technical Note No. 32, BKG (2004)
    # *
    # *  This revision:  2009 December 15
    # *
    # *  SOFA release 2010-12-01
    # *
    # *  Copyright (C) 2010 IAU SOFA Board.  See notes at end.
    # *
    # *-----------------------------------------------------------------------
    '''
    # *  2Pi
    D2PI = 6.283185307179586476925287

    # *  Reference epoch (J2000.0), JD
    DJ00 = 2451545.0

    # * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # *  Days since fundamental epoch.
    if ( DJ1 < DJ2 ):
        D1 = DJ1
        D2 = DJ2
    else:
        D1 = DJ2
        D2 = DJ1

    T = D1 + ( D2-DJ00 )

    # *  Fractional part of T (days).
    F = fmod ( D1, 1.0 ) + fmod ( D2, 1.0 )

    # *  Earth rotation angle at this UT1.
    return iau_ANP ( D2PI * ( F + 0.7790572732640 + 0.00273781191135448 * T ) )

    # *  Finished.


def iau_FAD03 ( T ):
    '''
    # *+
    # *  - - - - - - - - - -
    # *   i a u _ F A D 0 3
    # *  - - - - - - - - - -
    # *
    # *  Fundamental argument, IERS Conventions (2003):
    # *  mean elongation of the Moon from the Sun.
    # *
    # *  This routine is part of the International Astronomical Union's
    # *  SOFA (Standards of Fundamental Astronomy) software collection.
    # *
    # *  Status:  canonical remel.
    # *
    # *  Given:
    # *     T           d    TDB, Julian centuries since J2000.0 (Note 1)
    # *
    # *  Returned:
    # *     iau_FA.03   d    D, radians (Note 2)
    # *
    # *  Notes:
    # *
    # *  1) Though T is strictly TDB, it is usually more convenient to use TT,
    # *     which makes no significant difference.
    # *
    # *  2) The expression used is as adopted in IERS Conventions (2003) and
    # *     is from Simon et al. (1994).
    # *
    # *  References:
    # *
    # *     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    # *     IERS Technical Note No. 32, BKG (2004)
    # *
    # *     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
    # *     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
    # *
    # *  This revision:  2009 December 15
    # *
    # *  SOFA release 2010-12-01
    # *
    # *  Copyright (C) 2010 IAU SOFA Board.  See notes at end.
    # *
    # *-----------------------------------------------------------------------
    '''
    
    # *  Arcseconds to radians.
    DAS2R = 4.848136811095359935899141e-6

    # *  Arcseconds in a full circle.
    TURNAS = 1296000.0

    # * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # *  Mean elongation of the Moon from the Sun (IERS Conventions 2003).
    return fmod (      1072260.703692 + 
                        T*( 1602961601.2090 + 
                        T*(        - 6.3706 + 
                        T*(          0.006593 + 
                        T*(        - 0.00003169 )))), TURNAS ) * DAS2R

    # *  Finished.



def iau_FAE03 ( T ):
    '''
    # *+
    # *  - - - - - - - - - - -
    # *   i a u _ F A E 0 3
    # *  - - - - - - - - - - -
    # *
    # *  Fundamental argument, IERS Conventions (2003):
    # *  mean longitude of Earth.
    # *
    # *  This routine is part of the International Astronomical Union's
    # *  SOFA (Standards of Fundamental Astronomy) software collection.
    # *
    # *  Status:  canonical remel.
    # *
    # *  Given:
    # *     T           d    TDB, Julian centuries since J2000.0 (Note 1)
    # *
    # *  Returned:
    # *     iau_FAE03   d    mean longitude of Earth, radians (Note 2)
    # *
    # *  Notes:
    # *
    # *  1) Though T is strictly TDB, it is usually more convenient to use TT,
    # *     which makes no significant difference.
    # *
    # *  2) The expression used is as adopted in IERS Conventions (2003) and
    # *     comes from Souchay et al. (1999) after Simon et al. (1994).
    # *
    # *  References:
    # *
    # *     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    # *     IERS Technical Note No. 32, BKG (2004)
    # *
    # *     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
    # *     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
    # *
    # *     Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
    # *     Astron.Astrophys.Supp.Ser. 135, 111
    # *
    # *  This revision:  2009 December 15
    # *
    # *  SOFA release 2010-12-01
    # *
    # *  Copyright (C) 2010 IAU SOFA Board.  See notes at end.
    # *
    # *-----------------------------------------------------------------------
    '''
    
    # *  2Pi.
    D2PI = 6.283185307179586476925287

    # * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # *  Mean longitude of Earth (IERS Conventions 2003).
    return fmod( 1.753470314 + 628.3075849991 * T, D2PI )

    # *  Finished.


def iau_FAF03 ( T ):
    '''
    # *+
    # *  - - - - - - - - - -
    # *   i a u _ F A F 0 3
    # *  - - - - - - - - - -
    # *
    # *  Fundamental argument, IERS Conventions (2003):
    # *  mean longitude of the Moon minus mean longitude of the ascending
    # *  node.
    # *
    # *  This routine is part of the International Astronomical Union's
    # *  SOFA (Standards of Fundamental Astronomy) software collection.
    # *
    # *  Status:  canonical remel.
    # *
    # *  Given:
    # *     T           d    TDB, Julian centuries since J2000.0 (Note 1)
    # *
    # *  Returned:
    # *     iau_FAF03   d    F, radians (Note 2)
    # *
    # *  Notes:
    # *
    # *  1) Though T is strictly TDB, it is usually more convenient to use TT,
    # *     which makes no significant difference.
    # *
    # *  2) The expression used is as adopted in IERS Conventions (2003) and
    # *     is from Simon et al. (1994).
    # *
    # *  References:
    # *
    # *     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    # *     IERS Technical Note No. 32, BKG (2004)
    # *
    # *     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
    # *     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
    # *
    # *  This revision:  2009 December 15
    # *
    # *  SOFA release 2010-12-01
    # *
    # *  Copyright (C) 2010 IAU SOFA Board.  See notes at end.
    # *
    # *-----------------------------------------------------------------------
    '''
    
    # *  Arcseconds to radians.
    DAS2R = 4.848136811095359935899141e-6

    # *  Arcseconds in a full circle.
    TURNAS = 1296000.0

    # * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # *  Mean longitude of the Moon minus that of the ascending node
    # *  (IERS Conventions 2003).
    return fmod(       335779.526232 + 
                        T*( 1739527262.8478 + 
                        T*(       - 12.7512 + 
                        T*(       -  0.001037 + 
                        T*(          0.00000417 )))), TURNAS ) * DAS2R

# *  Finished.


def iau_FAJU03 ( T ):
    '''
    # *+
    # *  - - - - - - - - - - -
    # *   i a u _ F A J U 0 3
    # *  - - - - - - - - - - -
    # *
    # *  Fundamental argument, IERS Conventions (2003):
    # *  mean longitude of Jupiter.
    # *
    # *  This routine is part of the International Astronomical Union's
    # *  SOFA (Standards of Fundamental Astronomy) software collection.
    # *
    # *  Status:  canonical remel.
    # *
    # *  Given:
    # *     T           d    TDB, Julian centuries since J2000.0 (Note 1)
    # *
    # *  Returned:
    # *     iau_FAJU03  d    mean longitude of Jupiter, radians (Note 2)
    # *
    # *  Notes:
    # *
    # *  1) Though T is strictly TDB, it is usually more convenient to use TT,
    # *     which makes no significant difference.
    # *
    # *  2) The expression used is as adopted in IERS Conventions (2003) and
    # *     comes from Souchay et al. (1999) after Simon et al. (1994).
    # *
    # *  References:
    # *
    # *     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    # *     IERS Technical Note No. 32, BKG (2004)
    # *
    # *     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
    # *     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
    # *
    # *     Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
    # *     Astron.Astrophys.Supp.Ser. 135, 111
    # *
    # *  This revision:  2009 December 15
    # *
    # *  SOFA release 2010-12-01
    # *
    # *  Copyright (C) 2010 IAU SOFA Board.  See notes at end.
    # *
    # *-----------------------------------------------------------------------
    '''

    # *  2Pi.
    D2PI = 6.283185307179586476925287

    # * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # *  Mean longitude of Jupiter (IERS Conventions 2003).
    return fmod( 0.599546497 + 52.9690962641 * T, D2PI )

    # *  Finished.


def iau_FAL03 ( T ):
    '''
    # *+
    # *  - - - - - - - - - -
    # *   i a u _ F A L 0 3
    # *  - - - - - - - - - -
    # *
    # *  Fundamental argument, IERS Conventions (2003):
    # *  mean anomaly of the Moon.
    # *
    # *  This routine is part of the International Astronomical Union's
    # *  SOFA (Standards of Fundamental Astronomy) software collection.
    # *
    # *  Status:  canonical remel.
    # *
    # *  Given:
    # *     T           d    TDB, Julian centuries since J2000.0 (Note 1)
    # *
    # *  Returned:
    # *     iau_FAL03   d    l, radians (Note 2)
    # *
    # *  Notes:
    # *
    # *  1) Though T is strictly TDB, it is usually more convenient to use TT,
    # *     which makes no significant difference.
    # *
    # *  2) The expression used is as adopted in IERS Conventions (2003) and
    # *     is from Simon et al. (1994).
    # *
    # *  References:
    # *
    # *     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    # *     IERS Technical Note No. 32, BKG (2004)
    # *
    # *     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
    # *     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
    # *
    # *  This revision:  2009 December 15
    # *
    # *  SOFA release 2010-12-01
    # *
    # *  Copyright (C) 2010 IAU SOFA Board.  See notes at end.
    # *
    # *-----------------------------------------------------------------------
    '''
    # *  Arcseconds to radians.
    DAS2R = 4.848136811095359935899141e-6

    # *  Arcseconds in a full circle.
    TURNAS = 1296000.0

    # * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # *  Mean anomaly of the Moon (IERS Conventions 2003).
    return fmod(       485868.249036 + 
                        T*( 1717915923.2178 + 
                        T*(         31.8792 + 
                        T*(          0.051635 + 
                        T*(        - 0.00024470 )))), TURNAS ) * DAS2R
                
    # *  Finished.


def iau_FALP03 ( T ):
    '''
    # *+
    # *  - - - - - - - - - - -
    # *   i a u _ F A L P 0 3
    # *  - - - - - - - - - - -
    # *
    # *  Fundamental argument, IERS Conventions (2003):
    # *  mean anomaly of the Sun.
    # *
    # *  This routine is part of the International Astronomical Union's
    # *  SOFA (Standards of Fundamental Astronomy) software collection.
    # *
    # *  Status:  canonical remel.
    # *
    # *  Given:
    # *     T           d    TDB, Julian centuries since J2000.0 (Note 1)
    # *
    # *  Returned:
    # *     iau_FALP03  d    l', radians (Note 2)
    # *
    # *  Notes:
    # *
    # *  1) Though T is strictly TDB, it is usually more convenient to use TT,
    # *     which makes no significant difference.
    # *
    # *  2) The expression used is as adopted in IERS Conventions (2003) and
    # *     is from Simon et al. (1994).
    # *
    # *  References:
    # *
    # *     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    # *     IERS Technical Note No. 32, BKG (2004)
    # *
    # *     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
    # *     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
    # *
    # *  This revision:  2009 December 15
    # *
    # *  SOFA release 2010-12-01
    # *
    # *  Copyright (C) 2010 IAU SOFA Board.  See notes at end.
    # *
    # *-----------------------------------------------------------------------
    '''

    # *  Arcseconds to radians.
    DAS2R = 4.848136811095359935899141e-6

    # *  Arcseconds in a full circle.
    TURNAS = 1296000.0

    # * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # *  Mean anomaly of the Sun (IERS Conventions 2003).
    return fmod(     1287104.793048 + 
                         T*( 129596581.0481 + 
                         T*(       - 0.5532 + 
                         T*(         0.000136 + 
                         T*(       - 0.00001149 )))), TURNAS ) * DAS2R

    # *  Finished.


def iau_FAMA03 ( T ):
    '''
    # *+
    # *  - - - - - - - - - - -
    # *   i a u _ F A M A 0 3
    # *  - - - - - - - - - - -
    # *
    # *  Fundamental argument, IERS Conventions (2003):
    # *  mean longitude of Mars.
    # *
    # *  This routine is part of the International Astronomical Union's
    # *  SOFA (Standards of Fundamental Astronomy) software collection.
    # *
    # *  Status:  canonical remel.
    # *
    # *  Given:
    # *     T           d    TDB, Julian centuries since J2000.0 (Note 1)
    # *
    # *  Returned:
    # *     iau_FAMA03  d    mean longitude of Mars, radians (Note 2)
    # *
    # *  Notes:
    # *
    # *  1) Though T is strictly TDB, it is usually more convenient to use TT,
    # *     which makes no significant difference.
    # *
    # *  2) The expression used is as adopted in IERS Conventions (2003) and
    # *     comes from Souchay et al. (1999) after Simon et al. (1994).
    # *
    # *  References:
    # *
    # *     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    # *     IERS Technical Note No. 32, BKG (2004)
    # *
    # *     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
    # *     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
    # *
    # *     Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
    # *     Astron.Astrophys.Supp.Ser. 135, 111
    # *
    # *  This revision:  2009 December 15
    # *
    # *  SOFA release 2010-12-01
    # *
    # *  Copyright (C) 2010 IAU SOFA Board.  See notes at end.
    # *
    # *-----------------------------------------------------------------------
    '''

    # *  2Pi.
    D2PI = 6.283185307179586476925287

    # * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # *  Mean longitude of Mars (IERS Conventions 2003).
    return fmod( 6.203480913 + 334.0612426700 * T, D2PI )

    # *  Finished.


def iau_FAME03 ( T ):
    '''
    # *+
    # *  - - - - - - - - - - -
    # *   i a u _ F A M E 0 3
    # *  - - - - - - - - - - -
    # *
    # *  Fundamental argument, IERS Conventions (2003):
    # *  mean longitude of Mercury.
    # *
    # *  This routine is part of the International Astronomical Union's
    # *  SOFA (Standards of Fundamental Astronomy) software collection.
    # *
    # *  Status:  canonical remel.
    # *
    # *  Given:
    # *     T           d    TDB, Julian centuries since J2000.0 (Note 1)
    # *
    # *  Returned:
    # *     iau_FAME03  d    mean longitude of Mercury, radians (Note 2)
    # *
    # *  Notes:
    # *
    # *  1) Though T is strictly TDB, it is usually more convenient to use TT,
    # *     which makes no significant difference.
    # *
    # *  2) The expression used is as adopted in IERS Conventions (2003) and
    # *     comes from Souchay et al. (1999) after Simon et al. (1994).
    # *
    # *  References:
    # *
    # *     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    # *     IERS Technical Note No. 32, BKG (2004)
    # *
    # *     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
    # *     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
    # *
    # *     Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
    # *     Astron.Astrophys.Supp.Ser. 135, 111
    # *
    # *  This revision:  2009 December 15
    # *
    # *  SOFA release 2010-12-01
    # *
    # *  Copyright (C) 2010 IAU SOFA Board.  See notes at end.
    # *
    # *-----------------------------------------------------------------------
    '''

    # *  2Pi.
    D2PI = 6.283185307179586476925287

    # * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # *  Mean longitude of Mercury (IERS Conventions 2003).
    return fmod( 4.402608842 + 2608.7903141574 * T, D2PI )

    # *  Finished.


def iau_FANE03 ( T ):
    '''
    # *+
    # *  - - - - - - - - - - -
    # *   i a u _ F A N E 0 3
    # *  - - - - - - - - - - -
    # *
    # *  Fundamental argument, IERS Conventions (2003):
    # *  mean longitude of Neptune.
    # *
    # *  This routine is part of the International Astronomical Union's
    # *  SOFA (Standards of Fundamental Astronomy) software collection.
    # *
    # *  Status:  canonical remel.
    # *
    # *  Given:
    # *     T           d    TDB, Julian centuries since J2000.0 (Note 1)
    # *
    # *  Returned:
    # *     iau_FANE03  d    mean longitude of Neptune, radians (Note 2)
    # *
    # *  Notes:
    # *
    # *  1) Though T is strictly TDB, it is usually more convenient to use TT,
    # *     which makes no significant difference.
    # *
    # *  2) The expression used is as adopted in IERS Conventions (2003) and
    # *     is adapted from Simon et al. (1994).
    # *
    # *  References:
    # *
    # *     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    # *     IERS Technical Note No. 32, BKG (2004)
    # *
    # *     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
    # *     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
    # *
    # *  This revision:  2009 December 15
    # *
    # *  SOFA release 2010-12-01
    # *
    # *  Copyright (C) 2010 IAU SOFA Board.  See notes at end.
    # *
    # * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    '''

    # *  2Pi.
    D2PI = 6.283185307179586476925287

    # *  Mean longitude of Neptune (IERS Conventions 2003).
    return fmod( 5.311886287 + 3.8133035638 * T, D2PI )

    # *  Finished.


def iau_FAOM03 ( T ):
    '''
    # *+
    # *  - - - - - - - - - - -
    # *   i a u _ F A O M 0 3
    # *  - - - - - - - - - - -
    # *
    # *  Fundamental argument, IERS Conventions (2003):
    # *  mean longitude of the Moon's ascending node.
    # *
    # *  This routine is part of the International Astronomical Union's
    # *  SOFA (Standards of Fundamental Astronomy) software collection.
    # *
    # *  Status:  canonical remel.
    # *
    # *  Given:
    # *     T           d    TDB, Julian centuries since J2000.0 (Note 1)
    # *
    # *  Returned:
    # *     iau_FAOM03  d    Omega, radians (Note 2)
    # *
    # *  Notes:
    # *
    # *  1) Though T is strictly TDB, it is usually more convenient to use TT,
    # *     which makes no significant difference.
    # *
    # *  2) The expression used is as adopted in IERS Conventions (2003) and
    # *     is from Simon et al. (1994).
    # *
    # *  References:
    # *
    # *     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    # *     IERS Technical Note No. 32, BKG (2004)
    # *
    # *     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
    # *     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
    # *
    # *  This revision:  2009 December 15
    # *
    # *  SOFA release 2010-12-01
    # *
    # *  Copyright (C) 2010 IAU SOFA Board.  See notes at end.
    # *
    # *-----------------------------------------------------------------------
    '''

    # *  Arcseconds to radians.
    DAS2R = 4.848136811095359935899141e-6

    # *  Arcseconds in a full circle.
    TURNAS = 1296000.0

    # * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # *  Mean longitude of the Moon's ascending node (IERS Conventions 2003).
    return fmod(      450160.398036 + 
                         T*( - 6962890.5431 + 
                         T*(         7.4722 + 
                         T*(         0.007702 + 
                         T*(       - 0.00005939 )))), TURNAS ) * DAS2R

    # *  Finished.


def iau_FAPA03 ( T ):
    '''
    # *+
    # *  - - - - - - - - - - -
    # *   i a u _ F A P A 0 3
    # *  - - - - - - - - - - -
    # *
    # *  Fundamental argument, IERS Conventions (2003):
    # *  general accumulated precession in longitude.
    # *
    # *  This routine is part of the International Astronomical Union's
    # *  SOFA (Standards of Fundamental Astronomy) software collection.
    # *
    # *  Status:  canonical remel.
    # *
    # *  Given:
    # *     T           d    TDB, Julian centuries since J2000.0 (Note 1)
    # *
    # *  Returned:
    # *     iau_FAPA03  d    general precession in longitude, radians (Note 2)
    # *
    # *  Notes:
    # *
    # *  1) Though T is strictly TDB, it is usually more convenient to use TT,
    # *     which makes no significant difference.
    # *
    # *  2) The expression used is as adopted in IERS Conventions (2003).  It
    # *     is taken from Kinoshita & Souchay (1990) and comes originally from
    # *     Lieske et al. (1977).
    # *
    # *  References:
    # *
    # *     Kinoshita, H. and Souchay J. 1990, Celest.Mech. and Dyn.Astron.
    # *     48, 187
    # *
    # *     Lieske, J.H., Lederle, T., Fricke, W. & Morando, B. 1977,
    # *     Astron.Astrophys. 58, 1-16
    # *
    # *     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    # *     IERS Technical Note No. 32, BKG (2004)
    # *
    # *  This revision:  2009 December 15
    # *
    # *  SOFA release 2010-12-01
    # *
    # *  Copyright (C) 2010 IAU SOFA Board.  See notes at end.
    # *
    # *-----------------------------------------------------------------------
    '''

    # *  General accumulated precession in longitude.
    return ( 0.024381750 + 0.00000538691 * T ) * T

    # *  Finished.



def iau_FASA03 ( T ):
    '''
    # *+
    # *  - - - - - - - - - - -
    # *   i a u _ F A S A 0 3
    # *  - - - - - - - - - - -
    # *
    # *  Fundamental argument, IERS Conventions (2003):
    # *  mean longitude of Saturn.
    # *
    # *  This routine is part of the International Astronomical Union's
    # *  SOFA (Standards of Fundamental Astronomy) software collection.
    # *
    # *  Status:  canonical remel.
    # *
    # *  Given:
    # *     T           d    TDB, Julian centuries since J2000.0 (Note 1)
    # *
    # *  Returned:
    # *     iau_FASA03  d    mean longitude of Saturn, radians (Note 2)
    # *
    # *  Notes:
    # *
    # *  1) Though T is strictly TDB, it is usually more convenient to use TT,
    # *     which makes no significant difference.
    # *
    # *  2) The expression used is as adopted in IERS Conventions (2003) and
    # *     comes from Souchay et al. (1999) after Simon et al. (1994).
    # *
    # *  References:
    # *
    # *     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    # *     IERS Technical Note No. 32, BKG (2004)
    # *
    # *     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
    # *     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
    # *
    # *     Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
    # *     Astron.Astrophys.Supp.Ser. 135, 111
    # *
    # *  This revision:  2009 December 15
    # *
    # *  SOFA release 2010-12-01
    # *
    # *  Copyright (C) 2010 IAU SOFA Board.  See notes at end.
    # *
    # *-----------------------------------------------------------------------
    '''

    # *  2Pi.
    D2PI = 6.283185307179586476925287

    # *  Mean longitude of Saturn (IERS Conventions 2003).
    return fmod( 0.874016757 + 21.3299104960 * T, D2PI )

    # *  Finished.



def iau_FAUR03 ( T ):
    '''
    # *+
    # *  - - - - - - - - - - -
    # *   i a u _ F A U R 0 3
    # *  - - - - - - - - - - -
    # *
    # *  Fundamental argument, IERS Conventions (2003):
    # *  mean longitude of Uranus.
    # *
    # *  This routine is part of the International Astronomical Union's
    # *  SOFA (Standards of Fundamental Astronomy) software collection.
    # *
    # *  Status:  canonical remel.
    # *
    # *  Given:
    # *     T           d    TDB, Julian centuries since J2000.0 (Note 1)
    # *
    # *  Returned:
    # *     iau_FAUR03  d    mean longitude of Uranus, radians (Note 2)
    # *
    # *  Notes:
    # *
    # *  1) Though T is strictly TDB, it is usually more convenient to use TT,
    # *     which makes no significant difference.
    # *
    # *  2) The expression used is as adopted in IERS Conventions (2003) and
    # *     is adapted from Simon et al. (1994).
    # *
    # *  References:
    # *
    # *     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    # *     IERS Technical Note No. 32, BKG (2004)
    # *
    # *     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
    # *     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
    # *
    # *  This revision:  2009 December 15
    # *
    # *  SOFA release 2010-12-01
    # *
    # *  Copyright (C) 2010 IAU SOFA Board.  See notes at end.
    # *
    # *-----------------------------------------------------------------------
    '''

    # *  2Pi.
    D2PI = 6.283185307179586476925287

    # * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # *  Mean longitude of Uranus (IERS Conventions 2003).
    return fmod( 5.481293872 + 7.4781598567 * T, D2PI )

    # *  Finished.



def iau_FAVE03 ( T ):
    '''
    # *+
    # *  - - - - - - - - - - -
    # *   i a u _ F A V E 0 3
    # *  - - - - - - - - - - -
    # *
    # *  Fundamental argument, IERS Conventions (2003):
    # *  mean longitude of Venus.
    # *
    # *  This routine is part of the International Astronomical Union's
    # *  SOFA (Standards of Fundamental Astronomy) software collection.
    # *
    # *  Status:  canonical remel.
    # *
    # *  Given:
    # *     T           d    TDB, Julian centuries since J2000.0 (Note 1)
    # *
    # *  Returned:
    # *     iau_FAVE03  d    mean longitude of Venus, radians (Note 2)
    # *
    # *  Notes:
    # *
    # *  1) Though T is strictly TDB, it is usually more convenient to use TT,
    # *     which makes no significant difference.
    # *
    # *  2) The expression used is as adopted in IERS Conventions (2003) and
    # *     comes from Souchay et al. (1999) after Simon et al. (1994).
    # *
    # *  References:
    # *
    # *     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    # *     IERS Technical Note No. 32, BKG (2004)
    # *
    # *     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
    # *     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
    # *
    # *     Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
    # *     Astron.Astrophys.Supp.Ser. 135, 111
    # *
    # *  This revision:  2009 December 15
    # *
    # *  SOFA release 2010-12-01
    # *
    # *  Copyright (C) 2010 IAU SOFA Board.  See notes at end.
    # *
    # *-----------------------------------------------------------------------
    '''

    # *  2Pi.
    D2PI = 6.283185307179586476925287

    # * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # *  Mean longitude of Venus (IERS Conventions 2003).
    return fmod( 3.176146697 + 1021.3285546211 * T, D2PI )

    # *  Finished.


#@jit
def iau_NUT00A ( DATE1, DATE2 ):
    '''
    # *+
    # *  - - - - - - - - - - -
    # *   i a u _ N U T 0 0 A
    # *  - - - - - - - - - - -
    # *
    # *  Nutation, IAU 2000A remel (MHB2000 luni-solar and planetary nutation
    # *  with free core nutation omitted).
    # *
    # *  This routine is part of the International Astronomical Union's
    # *  SOFA (Standards of Fundamental Astronomy) software collection.
    # *
    # *  Status   canonical remel.
    # *
    # *  Given 
    # *     DATE1,DATE2   d    TT as a 2-part Julian Date (Note 1)
    # *
    # *  Returned 
    # *     DPSI,DEPS     d    nutation, luni-solar + planetary (Note 2)
    # *
    # *  Notes 
    # *
    # *  1) The TT date DATE1+DATE2 is a Julian Date, apportioned in any
    # *     convenient way between the two arguments.  For example,
    # *     JD(TT)=2450123.7 could be expressed in any of these ways,
    # *     among others 
    # *
    # *            DATE1          DATE2
    # *
    # *         2450123.7        0.0        (JD method)
    # *          2451545.0     -1421.3     (J2000 method)
    # *         2400000.5     50123.2    (MJD method)
    # *         2450123.5       0.2       (date & time method)
    # *
    # *     The JD method is the most natural and convenient to use in
    # *     cases where the loss of several decimal digits of resolution
    # *     is acceptable.  The J2000 method is best matched to the way
    # *     the argument is handled internally and will deliver the
    # *     optimum resolution.  The MJD method and the date & time methods
    # *     are both good compromises between resolution and convenience.
    # *
    # *  2) The nutation components in longitude and obliquity are in radians
    # *     and with respect to the equinox and ecliptic of date.  The
    # *     obliquity at J2000.0 is assumed to be the Lieske et al. (1977)
    # *     value of 84381.448 arcsec.
    # *
    # *     Both the luni-solar and planetary nutations are included.  The
    # *     latter are due to direct planetary nutations and the perturbations
    # *     of the lunar and terrestrial orbits.
    # *
    # *  3) The routine computes the MHB2000 nutation series with the
    # *     associated corrections for planetary nutations.  It is an
    # *     implementation of the nutation part of the IAU 2000A precession-
    # *     nutation remel, formally adopted by the IAU General Assembly in
    # *     2000, namely MHB2000 (Mathews et al. 2002), but with the free core
    # *     nutation (FCN - see Note 4) omitted.
    # *
    # *  4) The full MHB2000 remel also contains contributions to the
    # *     nutations in longitude and obliquity due to the free-excitation of
    # *     the free-core-nutation during the period 1979-2000.  These FCN
    # *     terms, which are time-dependent and unpredictable, are NOT
    # *     included in the present routine and, if required, must be
    # *     independently computed.  With the FCN corrections included, the
    # *     present routine delivers a pole which is at current epochs
    # *     accurate to a few hundred microarcseconds.  The omission of FCN
    # *     introduces further errors of about that size.
    # *
    # *  5) The present routine provides classical nutation.  The MHB2000
    # *     algorithm, from which it is adapted, deals also with (i) the
    # *     offsets between the GCRS and mean poles and (ii) the adjustments
    # *     in longitude and obliquity due to the changed precession rates.
    # *     These additional functions, namely frame bias and precession
    # *     adjustments, are supported by the SOFA routines iau_BI00 and
    # *     iau_PR00.
    # *
    # *  6) The MHB2000 algorithm also provides "total" nutations, comprising
    # *     the arithmetic sum of the frame bias, precession adjustments,
    # *     luni-solar nutation and planetary nutation.  These total nutations
    # *     can be used in combination with an existing IAU 1976 precession
    # *     implementation, such as iau_PMAT76, to deliver GCRS-to-true
    # *     predictions of sub-mas accuracy at current epochs.  However, there
    # *     are three shortcomings in the MHB2000 remel that must be taken
    # *     into account if more accurate or definitive results are required
    # *     (see Wallace 2002) 
    # *
    # *       (i) The MHB2000 total nutations are simply arithmetic sums,
    # *           yet in reality the various components are successive Euler
    # *           rotations.  This slight lack of rigor leads to cross terms
    # *           that exceed 1 mas after a century.  The rigorous procedure
    # *           is to form the GCRS-to-true rotation matrix by applying the
    # *           bias, precession and nutation in that order.
    # *
    # *      (ii) Although the precession adjustments are stated to be with
    # *           respect to Lieske et al. (1977), the MHB2000 remel does
    # *           not specify which set of Euler angles are to be used and
    # *           how the adjustments are to be applied.  The most literal and
    # *           straightforward procedure is to adopt the 4-rotation
    # *           epsilon_0, psi_A, omega_A, xi_A option, and to add DPSIPR to
    # *           psi_A and DEPSPR to both omega_A and eps_A.
    # *
    # *     (iii) The MHB2000 remel predates the determination by Chapront
    # *           et al. (2002) of a 14.6 mas displacement between the J2000.0
    # *           mean equinox and the origin of the ICRS frame.  It should,
    # *           however, be noted that neglecting this displacement when
    # *           calculating star coordinates does not lead to a 14.6 mas
    # *           change in right ascension, only a small second-order
    # *           distortion in the pattern of the precession-nutation effect.
    # *
    # *     For these reasons, the SOFA routines do not generate the "total
    # *     nutations" directly, though they can of course easily be generated
    # *     by calling iau_BI00, iau_PR00 and the present routine and adding
    # *     the results.
    # *
    # *  7) The MHB2000 remel contains 41 instances where the same frequency
    # *     appears multiple times, of which 38 are duplicates and three are
    # *     triplicates.  To keep the present code close to the original MHB
    # *     algorithm, this small inefficiency has not been corrected.
    # *
    # *  Called 
    # *     iau_FAL03    mean anomaly of the Moon
    # *     iau_FAF03    mean argument of the latitude of the Moon
    # *     iau_FAOM03   mean longitude of the Moon's ascending node
    # *     iau_FAME03   mean longitude of Mercury
    # *     iau_FAVE03   mean longitude of Venus
    # *     iau_FAE03    mean longitude of Earth
    # *     iau_FAMA03   mean longitude of Mars
    # *     iau_FAJU03   mean longitude of Jupiter
    # *     iau_FASA03   mean longitude of Saturn
    # *     iau_FAUR03   mean longitude of Uranus
    # *     iau_FAPA03   general accumulated precession in longitude
    # *
    # *  References 
    # *
    # *     Chapront, J., Chapront-Touze, M. & Francou, G. 2002,
    # *     Astron.Astrophys. 387, 700
    # *
    # *     Lieske, J.H., Lederle, T., Fricke, W. & Morando, B. 1977,
    # *     Astron.Astrophys. 58, 1-16
    # *
    # *     Mathews, P.M., Herring, T.A., Buffet, B.A. 2002, J.Geophys.Res.
    # *     107, B4.  The MHB_2000 code itself was obtained on 9th September
    # *     2002 from ftp//maia.usno.navy.mil/conv2000/chapter5/IAU2000A.
    # *
    # *     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
    # *     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
    # *
    # *     Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
    # *     Astron.Astrophys.Supp.Ser. 135, 111
    # *
    # *     Wallace, P.T., "Software for Implementing the IAU 2000
    # *     Resolutions", in IERS Workshop 5.1 (2002)
    # *
    # *  This revision   2009 December 15
    # *
    # *  SOFA release 2010-12-01
    # *
    # *  Copyright (C) 2010 IAU SOFA Board.  See notes at end.
    # *
    # *-----------------------------------------------------------------------
    '''
    
    # *  Arcseconds to radians
    DAS2R = 4.848136811095359935899141e-6

    # *  Arcseconds in a full circle
    TURNAS = 1296000.0

    # *  2Pi
    D2PI = 6.283185307179586476925287

    # *  Units of 0.1 microarcsecond to radians
    U2R = DAS2R/1e7

    # *  Reference epoch (J2000.0), JD
    DJ00 = 2451545.0

    # *  Days per Julian century
    DJC = 36525.0

    # *  -------------------------
    # *  Luni-Solar nutation remel
    # *  -------------------------

    # *  Number of terms in the luni-solar nutation remel
    NLS = 678

    # *  Coefficients for fundamental arguments
    #NALS = np.zeros((NLS,5))

    # *  Longitude and obliquity coefficients
    #CLS = np.zeros((NLS,6))

    # *  ---------------
    # *  Planetary terms
    # *  ---------------

    # *  Number of terms in the planetary nutation remel
    NPL = 687
    
    #NAPL = np.zeros((NPL,14))
    
    # *  Longitude and obliquity coefficients
    #ICPL = np.zeros((NPL,4))
    
    # *  ----------------------------------------
    # *  Tables of argument and term coefficients
    # *  ----------------------------------------
    
    #
    #  Luni-Solar argument multipliers
    #               L     L'    F     D     Om
    NALS = np.asarray([
              [ 0,    0,    0,    0,    1],
              [ 0,    0,    2,   -2,    2],
              [ 0,    0,    2,    0,    2],
              [ 0,    0,    0,    0,    2],
              [ 0,    1,    0,    0,    0],
              [ 0,    1,    2,   -2,    2],
              [ 1,    0,    0,    0,    0],
              [ 0,    0,    2,    0,    1],
              [ 1,    0,    2,    0,    2],
              [ 0,   -1,    2,   -2,    2],
              [ 0,    0,    2,   -2,    1],
              [-1,    0,    2,    0,    2],
              [-1,    0,    0,    2,    0],
              [ 1,    0,    0,    0,    1],
              [-1,    0,    0,    0,    1],
              [-1,    0,    2,    2,    2],
              [ 1,    0,    2,    0,    1],
              [-2,    0,    2,    0,    1],
              [ 0,    0,    0,    2,    0],
              [ 0,    0,    2,    2,    2],
              [ 0,   -2,    2,   -2,    2],
              [-2,    0,    0,    2,    0],
              [ 2,    0,    2,    0,    2],
              [ 1,    0,    2,   -2,    2],
              [-1,    0,    2,    0,    1],
              [ 2,    0,    0,    0,    0],
              [ 0,    0,    2,    0,    0],
              [ 0,    1,    0,    0,    1],
              [-1,    0,    0,    2,    1],
              [ 0,    2,    2,   -2,    2],
              [ 0,    0,   -2,    2,    0],
              [ 1,    0,    0,   -2,    1],
              [ 0,   -1,    0,    0,    1],
              [-1,    0,    2,    2,    1],
              [ 0,    2,    0,    0,    0],
              [ 1,    0,    2,    2,    2],
              [-2,    0,    2,    0,    0],
              [ 0,    1,    2,    0,    2],
              [ 0,    0,    2,    2,    1],
              [ 0,   -1,    2,    0,    2],
              [ 0,    0,    0,    2,    1],
              [ 1,    0,    2,   -2,    1],
              [ 2,    0,    2,   -2,    2],
              [-2,    0,    0,    2,    1],
              [ 2,    0,    2,    0,    1],
              [ 0,   -1,    2,   -2,    1],
              [ 0,    0,    0,   -2,    1],
              [-1,   -1,    0,    2,    0],
              [ 2,    0,    0,   -2,    1],
              [ 1,    0,    0,    2,    0],
              [ 0,    1,    2,   -2,    1],
              [ 1,   -1,    0,    0,    0],
              [-2,    0,    2,    0,    2],
              [ 3,    0,    2,    0,    2],
              [ 0,   -1,    0,    2,    0],
              [ 1,   -1,    2,    0,    2],
              [ 0,    0,    0,    1,    0],
              [-1,   -1,    2,    2,    2],
              [-1,    0,    2,    0,    0],
              [ 0,   -1,    2,    2,    2],
              [-2,    0,    0,    0,    1],
              [ 1,    1,    2,    0,    2],
              [ 2,    0,    0,    0,    1],
              [-1,    1,    0,    1,    0],
              [ 1,    1,    0,    0,    0],
              [ 1,    0,    2,    0,    0],
              [-1,    0,    2,   -2,    1],
              [ 1,    0,    0,    0,    2],
              [-1,    0,    0,    1,    0],
              [ 0,    0,    2,    1,    2],
              [-1,    0,    2,    4,    2],
              [-1,    1,    0,    1,    1],
              [ 0,   -2,    2,   -2,    1],
              [ 1,    0,    2,    2,    1],
              [-2,    0,    2,    2,    2],
              [-1,    0,    0,    0,    2],
              [ 1,    1,    2,   -2,    2],
              [-2,    0,    2,    4,    2],
              [-1,    0,    4,    0,    2],
              [ 2,    0,    2,   -2,    1],
              [ 2,    0,    2,    2,    2],
              [ 1,    0,    0,    2,    1],
              [ 3,    0,    0,    0,    0],
              [ 3,    0,    2,   -2,    2],
              [ 0,    0,    4,   -2,    2],
              [ 0,    1,    2,    0,    1],
              [ 0,    0,   -2,    2,    1],
              [ 0,    0,    2,   -2,    3],
              [-1,    0,    0,    4,    0],
              [ 2,    0,   -2,    0,    1],
              [-2,    0,    0,    4,    0],
              [-1,   -1,    0,    2,    1],
              [-1,    0,    0,    1,    1],
              [ 0,    1,    0,    0,    2],
              [ 0,    0,   -2,    0,    1],
              [ 0,   -1,    2,    0,    1],
              [ 0,    0,    2,   -1,    2],
              [ 0,    0,    2,    4,    2],
              [-2,   -1,    0,    2,    0],
              [ 1,    1,    0,   -2,    1],
              [-1,    1,    0,    2,    0],
              [-1,    1,    0,    1,    2],
              [ 1,   -1,    0,    0,    1],
              [ 1,   -1,    2,    2,    2],
              [-1,    1,    2,    2,    2],
              [ 3,    0,    2,    0,    1],
              [ 0,    1,   -2,    2,    0],
              [-1,    0,    0,   -2,    1],
              [ 0,    1,    2,    2,    2],
              [-1,   -1,    2,    2,    1],
              [ 0,   -1,    0,    0,    2],
              [ 1,    0,    2,   -4,    1],
              [-1,    0,   -2,    2,    0],
              [ 0,   -1,    2,    2,    1],
              [ 2,   -1,    2,    0,    2],
              [ 0,    0,    0,    2,    2],
              [ 1,   -1,    2,    0,    1],
              [-1,    1,    2,    0,    2],
              [ 0,    1,    0,    2,    0],
              [ 0,   -1,   -2,    2,    0],
              [ 0,    3,    2,   -2,    2],
              [ 0,    0,    0,    1,    1],
              [-1,    0,    2,    2,    0],
              [ 2,    1,    2,    0,    2],
              [ 1,    1,    0,    0,    1],
              [ 1,    1,    2,    0,    1],
              [ 2,    0,    0,    2,    0],
              [ 1,    0,   -2,    2,    0],
              [-1,    0,    0,    2,    2],
              [ 0,    1,    0,    1,    0],
              [ 0,    1,    0,   -2,    1],
              [-1,    0,    2,   -2,    2],
              [ 0,    0,    0,   -1,    1],
              [-1,    1,    0,    0,    1],
              [ 1,    0,    2,   -1,    2],
              [ 1,   -1,    0,    2,    0],
              [ 0,    0,    0,    4,    0],
              [ 1,    0,    2,    1,    2],
              [ 0,    0,    2,    1,    1],
              [ 1,    0,    0,   -2,    2],
              [-1,    0,    2,    4,    1],
              [ 1,    0,   -2,    0,    1],
              [ 1,    1,    2,   -2,    1],
              [ 0,    0,    2,    2,    0],
              [-1,    0,    2,   -1,    1],
              [-2,    0,    2,    2,    1],
              [ 4,    0,    2,    0,    2],
              [ 2,   -1,    0,    0,    0],
              [ 2,    1,    2,   -2,    2],
              [ 0,    1,    2,    1,    2],
              [ 1,    0,    4,   -2,    2],
              [-1,   -1,    0,    0,    1],
              [ 0,    1,    0,    2,    1],
              [-2,    0,    2,    4,    1],
              [ 2,    0,    2,    0,    0],
              [ 1,    0,    0,    1,    0],
              [-1,    0,    0,    4,    1],
              [-1,    0,    4,    0,    1],
              [ 2,    0,    2,    2,    1],
              [ 0,    0,    2,   -3,    2],
              [-1,   -2,    0,    2,    0],
              [ 2,    1,    0,    0,    0],
              [ 0,    0,    4,    0,    2],
              [ 0,    0,    0,    0,    3],
              [ 0,    3,    0,    0,    0],
              [ 0,    0,    2,   -4,    1],
              [ 0,   -1,    0,    2,    1],
              [ 0,    0,    0,    4,    1],
              [-1,   -1,    2,    4,    2],
              [ 1,    0,    2,    4,    2],
              [-2,    2,    0,    2,    0],
              [-2,   -1,    2,    0,    1],
              [-2,    0,    0,    2,    2],
              [-1,   -1,    2,    0,    2],
              [ 0,    0,    4,   -2,    1],
              [ 3,    0,    2,   -2,    1],
              [-2,   -1,    0,    2,    1],
              [ 1,    0,    0,   -1,    1],
              [ 0,   -2,    0,    2,    0],
              [-2,    0,    0,    4,    1],
              [-3,    0,    0,    0,    1],
              [ 1,    1,    2,    2,    2],
              [ 0,    0,    2,    4,    1],
              [ 3,    0,    2,    2,    2],
              [-1,    1,    2,   -2,    1],
              [ 2,    0,    0,   -4,    1],
              [ 0,    0,    0,   -2,    2],
              [ 2,    0,    2,   -4,    1],
              [-1,    1,    0,    2,    1],
              [ 0,    0,    2,   -1,    1],
              [ 0,   -2,    2,    2,    2],
              [ 2,    0,    0,    2,    1],
              [ 4,    0,    2,   -2,    2],
              [ 2,    0,    0,   -2,    2],
              [ 0,    2,    0,    0,    1],
              [ 1,    0,    0,   -4,    1],
              [ 0,    2,    2,   -2,    1],
              [-3,    0,    0,    4,    0],
              [-1,    1,    2,    0,    1],
              [-1,   -1,    0,    4,    0],
              [-1,   -2,    2,    2,    2],
              [-2,   -1,    2,    4,    2],
              [ 1,   -1,    2,    2,    1],
              [-2,    1,    0,    2,    0],
              [-2,    1,    2,    0,    1],
              [ 2,    1,    0,   -2,    1],
              [-3,    0,    2,    0,    1],
              [-2,    0,    2,   -2,    1],
              [-1,    1,    0,    2,    2],
              [ 0,   -1,    2,   -1,    2],
              [-1,    0,    4,   -2,    2],
              [ 0,   -2,    2,    0,    2],
              [-1,    0,    2,    1,    2],
              [ 2,    0,    0,    0,    2],
              [ 0,    0,    2,    0,    3],
              [-2,    0,    4,    0,    2],
              [-1,    0,   -2,    0,    1],
              [-1,    1,    2,    2,    1],
              [ 3,    0,    0,    0,    1],
              [-1,    0,    2,    3,    2],
              [ 2,   -1,    2,    0,    1],
              [ 0,    1,    2,    2,    1],
              [ 0,   -1,    2,    4,    2],
              [ 2,   -1,    2,    2,    2],
              [ 0,    2,   -2,    2,    0],
              [-1,   -1,    2,   -1,    1],
              [ 0,   -2,    0,    0,    1],
              [ 1,    0,    2,   -4,    2],
              [ 1,   -1,    0,   -2,    1],
              [-1,   -1,    2,    0,    1],
              [ 1,   -1,    2,   -2,    2],
              [-2,   -1,    0,    4,    0],
              [-1,    0,    0,    3,    0],
              [-2,   -1,    2,    2,    2],
              [ 0,    2,    2,    0,    2],
              [ 1,    1,    0,    2,    0],
              [ 2,    0,    2,   -1,    2],
              [ 1,    0,    2,    1,    1],
              [ 4,    0,    0,    0,    0],
              [ 2,    1,    2,    0,    1],
              [ 3,   -1,    2,    0,    2],
              [-2,    2,    0,    2,    1],
              [ 1,    0,    2,   -3,    1],
              [ 1,    1,    2,   -4,    1],
              [-1,   -1,    2,   -2,    1],
              [ 0,   -1,    0,   -1,    1],
              [ 0,   -1,    0,   -2,    1],
              [-2,    0,    0,    0,    2],
              [-2,    0,   -2,    2,    0],
              [-1,    0,   -2,    4,    0],
              [ 1,   -2,    0,    0,    0],
              [ 0,    1,    0,    1,    1],
              [-1,    2,    0,    2,    0],
              [ 1,   -1,    2,   -2,    1],
              [ 1,    2,    2,   -2,    2],
              [ 2,   -1,    2,   -2,    2],
              [ 1,    0,    2,   -1,    1],
              [ 2,    1,    2,   -2,    1],
              [-2,    0,    0,   -2,    1],
              [ 1,   -2,    2,    0,    2],
              [ 0,    1,    2,    1,    1],
              [ 1,    0,    4,   -2,    1],
              [-2,    0,    4,    2,    2],
              [ 1,    1,    2,    1,    2],
              [ 1,    0,    0,    4,    0],
              [ 1,    0,    2,    2,    0],
              [ 2,    0,    2,    1,    2],
              [ 3,    1,    2,    0,    2],
              [ 4,    0,    2,    0,    1],
              [-2,   -1,    2,    0,    0],
              [ 0,    1,   -2,    2,    1],
              [ 1,    0,   -2,    1,    0],
              [ 0,   -1,   -2,    2,    1],
              [ 2,   -1,    0,   -2,    1],
              [-1,    0,    2,   -1,    2],
              [ 1,    0,    2,   -3,    2],
              [ 0,    1,    2,   -2,    3],
              [ 0,    0,    2,   -3,    1],
              [-1,    0,   -2,    2,    1],
              [ 0,    0,    2,   -4,    2],
              [-2,    1,    0,    0,    1],
              [-1,    0,    0,   -1,    1],
              [ 2,    0,    2,   -4,    2],
              [ 0,    0,    4,   -4,    4],
              [ 0,    0,    4,   -4,    2],
              [-1,   -2,    0,    2,    1],
              [-2,    0,    0,    3,    0],
              [ 1,    0,   -2,    2,    1],
              [-3,    0,    2,    2,    2],
              [-3,    0,    2,    2,    1],
              [-2,    0,    2,    2,    0],
              [ 2,   -1,    0,    0,    1],
              [-2,    1,    2,    2,    2],
              [ 1,    1,    0,    1,    0],
              [ 0,    1,    4,   -2,    2],
              [-1,    1,    0,   -2,    1],
              [ 0,    0,    0,   -4,    1],
              [ 1,   -1,    0,    2,    1],
              [ 1,    1,    0,    2,    1],
              [-1,    2,    2,    2,    2],
              [ 3,    1,    2,   -2,    2],
              [ 0,   -1,    0,    4,    0],
              [ 2,   -1,    0,    2,    0],
              [ 0,    0,    4,    0,    1],
              [ 2,    0,    4,   -2,    2],
              [-1,   -1,    2,    4,    1],
              [ 1,    0,    0,    4,    1],
              [ 1,   -2,    2,    2,    2],
              [ 0,    0,    2,    3,    2],
              [-1,    1,    2,    4,    2],
              [ 3,    0,    0,    2,    0],
              [-1,    0,    4,    2,    2],
              [ 1,    1,    2,    2,    1],
              [-2,    0,    2,    6,    2],
              [ 2,    1,    2,    2,    2],
              [-1,    0,    2,    6,    2],
              [ 1,    0,    2,    4,    1],
              [ 2,    0,    2,    4,    2],
              [ 1,    1,   -2,    1,    0],
              [-3,    1,    2,    1,    2],
              [ 2,    0,   -2,    0,    2],
              [-1,    0,    0,    1,    2],
              [-4,    0,    2,    2,    1],
              [-1,   -1,    0,    1,    0],
              [ 0,    0,   -2,    2,    2],
              [ 1,    0,    0,   -1,    2],
              [ 0,   -1,    2,   -2,    3],
              [-2,    1,    2,    0,    0],
              [ 0,    0,    2,   -2,    4],
              [-2,   -2,    0,    2,    0],
              [-2,    0,   -2,    4,    0],
              [ 0,   -2,   -2,    2,    0],
              [ 1,    2,    0,   -2,    1],
              [ 3,    0,    0,   -4,    1],
              [-1,    1,    2,   -2,    2],
              [ 1,   -1,    2,   -4,    1],
              [ 1,    1,    0,   -2,    2],
              [-3,    0,    2,    0,    0],
              [-3,    0,    2,    0,    2],
              [-2,    0,    0,    1,    0],
              [ 0,    0,   -2,    1,    0],
              [-3,    0,    0,    2,    1],
              [-1,   -1,   -2,    2,    0],
              [ 0,    1,    2,   -4,    1],
              [ 2,    1,    0,   -4,    1],
              [ 0,    2,    0,   -2,    1],
              [ 1,    0,    0,   -3,    1],
              [-2,    0,    2,   -2,    2],
              [-2,   -1,    0,    0,    1],
              [-4,    0,    0,    2,    0],
              [ 1,    1,    0,   -4,    1],
              [-1,    0,    2,   -4,    1],
              [ 0,    0,    4,   -4,    1],
              [ 0,    3,    2,   -2,    2],
              [-3,   -1,    0,    4,    0],
              [-3,    0,    0,    4,    1],
              [ 1,   -1,   -2,    2,    0],
              [-1,   -1,    0,    2,    2],
              [ 1,   -2,    0,    0,    1],
              [ 1,   -1,    0,    0,    2],
              [ 0,    0,    0,    1,    2],
              [-1,   -1,    2,    0,    0],
              [ 1,   -2,    2,   -2,    2],
              [ 0,   -1,    2,   -1,    1],
              [-1,    0,    2,    0,    3],
              [ 1,    1,    0,    0,    2],
              [-1,    1,    2,    0,    0],
              [ 1,    2,    0,    0,    0],
              [-1,    2,    2,    0,    2],
              [-1,    0,    4,   -2,    1],
              [ 3,    0,    2,   -4,    2],
              [ 1,    2,    2,   -2,    1],
              [ 1,    0,    4,   -4,    2],
              [-2,   -1,    0,    4,    1],
              [ 0,   -1,    0,    2,    2],
              [-2,    1,    0,    4,    0],
              [-2,   -1,    2,    2,    1],
              [ 2,    0,   -2,    2,    0],
              [ 1,    0,    0,    1,    1],
              [ 0,    1,    0,    2,    2],
              [ 1,   -1,    2,   -1,    2],
              [-2,    0,    4,    0,    1],
              [ 2,    1,    0,    0,    1],
              [ 0,    1,    2,    0,    0],
              [ 0,   -1,    4,   -2,    2],
              [ 0,    0,    4,   -2,    4],
              [ 0,    2,    2,    0,    1],
              [-3,    0,    0,    6,    0],
              [-1,   -1,    0,    4,    1],
              [ 1,   -2,    0,    2,    0],
              [-1,    0,    0,    4,    2],
              [-1,   -2,    2,    2,    1],
              [-1,    0,    0,   -2,    2],
              [ 1,    0,   -2,   -2,    1],
              [ 0,    0,   -2,   -2,    1],
              [-2,    0,   -2,    0,    1],
              [ 0,    0,    0,    3,    1],
              [ 0,    0,    0,    3,    0],
              [-1,    1,    0,    4,    0],
              [-1,   -1,    2,    2,    0],
              [-2,    0,    2,    3,    2],
              [ 1,    0,    0,    2,    2],
              [ 0,   -1,    2,    1,    2],
              [ 3,   -1,    0,    0,    0],
              [ 2,    0,    0,    1,    0],
              [ 1,   -1,    2,    0,    0],
              [ 0,    0,    2,    1,    0],
              [ 1,    0,    2,    0,    3],
              [ 3,    1,    0,    0,    0],
              [ 3,   -1,    2,   -2,    2],
              [ 2,    0,    2,   -1,    1],
              [ 1,    1,    2,    0,    0],
              [ 0,    0,    4,   -1,    2],
              [ 1,    2,    2,    0,    2],
              [-2,    0,    0,    6,    0],
              [ 0,   -1,    0,    4,    1],
              [-2,   -1,    2,    4,    1],
              [ 0,   -2,    2,    2,    1],
              [ 0,   -1,    2,    2,    0],
              [-1,    0,    2,    3,    1],
              [-2,    1,    2,    4,    2],
              [ 2,    0,    0,    2,    2],
              [ 2,   -2,    2,    0,    2],
              [-1,    1,    2,    3,    2],
              [ 3,    0,    2,   -1,    2],
              [ 4,    0,    2,   -2,    1],
              [-1,    0,    0,    6,    0],
              [-1,   -2,    2,    4,    2],
              [-3,    0,    2,    6,    2],
              [-1,    0,    2,    4,    0],
              [ 3,    0,    0,    2,    1],
              [ 3,   -1,    2,    0,    1],
              [ 3,    0,    2,    0,    0],
              [ 1,    0,    4,    0,    2],
              [ 5,    0,    2,   -2,    2],
              [ 0,   -1,    2,    4,    1],
              [ 2,   -1,    2,    2,    1],
              [ 0,    1,    2,    4,    2],
              [ 1,   -1,    2,    4,    2],
              [ 3,   -1,    2,    2,    2],
              [ 3,    0,    2,    2,    1],
              [ 5,    0,    2,    0,    2],
              [ 0,    0,    2,    6,    2],
              [ 4,    0,    2,    2,    2],
              [ 0,   -1,    1,   -1,    1],
              [-1,    0,    1,    0,    3],
              [ 0,   -2,    2,   -2,    3],
              [ 1,    0,   -1,    0,    1],
              [ 2,   -2,    0,   -2,    1],
              [-1,    0,    1,    0,    2],
              [-1,    0,    1,    0,    1],
              [-1,   -1,    2,   -1,    2],
              [-2,    2,    0,    2,    2],
              [-1,    0,    1,    0,    0],
              [-4,    1,    2,    2,    2],
              [-3,    0,    2,    1,    1],
              [-2,   -1,    2,    0,    2],
              [ 1,    0,   -2,    1,    1],
              [ 2,   -1,   -2,    0,    1],
              [-4,    0,    2,    2,    0],
              [-3,    1,    0,    3,    0],
              [-1,    0,   -1,    2,    0],
              [ 0,   -2,    0,    0,    2],
              [ 0,   -2,    0,    0,    2],
              [-3,    0,    0,    3,    0],
              [-2,   -1,    0,    2,    2],
              [-1,    0,   -2,    3,    0],
              [-4,    0,    0,    4,    0],
              [ 2,    1,   -2,    0,    1],
              [ 2,   -1,    0,   -2,    2],
              [ 0,    0,    1,   -1,    0],
              [-1,    2,    0,    1,    0],
              [-2,    1,    2,    0,    2],
              [ 1,    1,    0,   -1,    1],
              [ 1,    0,    1,   -2,    1],
              [ 0,    2,    0,    0,    2],
              [ 1,   -1,    2,   -3,    1],
              [-1,    1,    2,   -1,    1],
              [-2,    0,    4,   -2,    2],
              [-2,    0,    4,   -2,    1],
              [-2,   -2,    0,    2,    1],
              [-2,    0,   -2,    4,    0],
              [ 1,    2,    2,   -4,    1],
              [ 1,    1,    2,   -4,    2],
              [-1,    2,    2,   -2,    1],
              [ 2,    0,    0,   -3,    1],
              [-1,    2,    0,    0,    1],
              [ 0,    0,    0,   -2,    0],
              [-1,   -1,    2,   -2,    2],
              [-1,    1,    0,    0,    2],
              [ 0,    0,    0,   -1,    2],
              [-2,    1,    0,    1,    0],
              [ 1,   -2,    0,   -2,    1],
              [ 1,    0,   -2,    0,    2],
              [-3,    1,    0,    2,    0],
              [-1,    1,   -2,    2,    0],
              [-1,   -1,    0,    0,    2],
              [-3,    0,    0,    2,    0],
              [-3,   -1,    0,    2,    0],
              [ 2,    0,    2,   -6,    1],
              [ 0,    1,    2,   -4,    2],
              [ 2,    0,    0,   -4,    2],
              [-2,    1,    2,   -2,    1],
              [ 0,   -1,    2,   -4,    1],
              [ 0,    1,    0,   -2,    2],
              [-1,    0,    0,   -2,    0],
              [ 2,    0,   -2,   -2,    1],
              [-4,    0,    2,    0,    1],
              [-1,   -1,    0,   -1,    1],
              [ 0,    0,   -2,    0,    2],
              [-3,    0,    0,    1,    0],
              [-1,    0,   -2,    1,    0],
              [-2,    0,   -2,    2,    1],
              [ 0,    0,   -4,    2,    0],
              [-2,   -1,   -2,    2,    0],
              [ 1,    0,    2,   -6,    1],
              [-1,    0,    2,   -4,    2],
              [ 1,    0,    0,   -4,    2],
              [ 2,    1,    2,   -4,    2],
              [ 2,    1,    2,   -4,    1],
              [ 0,    1,    4,   -4,    4],
              [ 0,    1,    4,   -4,    2],
              [-1,   -1,   -2,    4,    0],
              [-1,   -3,    0,    2,    0],
              [-1,    0,   -2,    4,    1],
              [-2,   -1,    0,    3,    0],
              [ 0,    0,   -2,    3,    0],
              [-2,    0,    0,    3,    1],
              [ 0,   -1,    0,    1,    0],
              [-3,    0,    2,    2,    0],
              [ 1,    1,   -2,    2,    0],
              [-1,    1,    0,    2,    2],
              [ 1,   -2,    2,   -2,    1],
              [ 0,    0,    1,    0,    2],
              [ 0,    0,    1,    0,    1],
              [ 0,    0,    1,    0,    0],
              [-1,    2,    0,    2,    1],
              [ 0,    0,    2,    0,    2],
              [-2,    0,    2,    0,    2],
              [ 2,    0,    0,   -1,    1],
              [ 3,    0,    0,   -2,    1],
              [ 1,    0,    2,   -2,    3],
              [ 1,    2,    0,    0,    1],
              [ 2,    0,    2,   -3,    2],
              [-1,    1,    4,   -2,    2],
              [-2,   -2,    0,    4,    0],
              [ 0,   -3,    0,    2,    0],
              [ 0,    0,   -2,    4,    0],
              [-1,   -1,    0,    3,    0],
              [-2,    0,    0,    4,    2],
              [-1,    0,    0,    3,    1],
              [ 2,   -2,    0,    0,    0],
              [ 1,   -1,    0,    1,    0],
              [-1,    0,    0,    2,    0],
              [ 0,   -2,    2,    0,    1],
              [-1,    0,    1,    2,    1],
              [-1,    1,    0,    3,    0],
              [-1,   -1,    2,    1,    2],
              [ 0,   -1,    2,    0,    0],
              [-2,    1,    2,    2,    1],
              [ 2,   -2,    2,   -2,    2],
              [ 1,    1,    0,    1,    1],
              [ 1,    0,    1,    0,    1],
              [ 1,    0,    1,    0,    0],
              [ 0,    2,    0,    2,    0],
              [ 2,   -1,    2,   -2,    1],
              [ 0,   -1,    4,   -2,    1],
              [ 0,    0,    4,   -2,    3],
              [ 0,    1,    4,   -2,    1],
              [ 4,    0,    2,   -4,    2],
              [ 2,    2,    2,   -2,    2],
              [ 2,    0,    4,   -4,    2],
              [-1,   -2,    0,    4,    0],
              [-1,   -3,    2,    2,    2],
              [-3,    0,    2,    4,    2],
              [-3,    0,    2,   -2,    1],
              [-1,   -1,    0,   -2,    1],
              [-3,    0,    0,    0,    2],
              [-3,    0,   -2,    2,    0],
              [ 0,    1,    0,   -4,    1],
              [-2,    1,    0,   -2,    1],
              [-4,    0,    0,    0,    1],
              [-1,    0,    0,   -4,    1],
              [-3,    0,    0,   -2,    1],
              [ 0,    0,    0,    3,    2],
              [-1,    1,    0,    4,    1],
              [ 1,   -2,    2,    0,    1],
              [ 0,    1,    0,    3,    0],
              [-1,    0,    2,    2,    3],
              [ 0,    0,    2,    2,    2],
              [-2,    0,    2,    2,    2],
              [-1,    1,    2,    2,    0],
              [ 3,    0,    0,    0,    2],
              [ 2,    1,    0,    1,    0],
              [ 2,   -1,    2,   -1,    2],
              [ 0,    0,    2,    0,    1],
              [ 0,    0,    3,    0,    3],
              [ 0,    0,    3,    0,    2],
              [-1,    2,    2,    2,    1],
              [-1,    0,    4,    0,    0],
              [ 1,    2,    2,    0,    1],
              [ 3,    1,    2,   -2,    1],
              [ 1,    1,    4,   -2,    2],
              [-2,   -1,    0,    6,    0],
              [ 0,   -2,    0,    4,    0],
              [-2,    0,    0,    6,    1],
              [-2,   -2,    2,    4,    2],
              [ 0,   -3,    2,    2,    2],
              [ 0,    0,    0,    4,    2],
              [-1,   -1,    2,    3,    2],
              [-2,    0,    2,    4,    0],
              [ 2,   -1,    0,    2,    1],
              [ 1,    0,    0,    3,    0],
              [ 0,    1,    0,    4,    1],
              [ 0,    1,    0,    4,    0],
              [ 1,   -1,    2,    1,    2],
              [ 0,    0,    2,    2,    3],
              [ 1,    0,    2,    2,    2],
              [-1,    0,    2,    2,    2],
              [-2,    0,    4,    2,    1],
              [ 2,    1,    0,    2,    1],
              [ 2,    1,    0,    2,    0],
              [ 2,   -1,    2,    0,    0],
              [ 1,    0,    2,    1,    0],
              [ 0,    1,    2,    2,    0],
              [ 2,    0,    2,    0,    3],
              [ 3,    0,    2,    0,    2],
              [ 1,    0,    2,    0,    2],
              [ 1,    0,    3,    0,    3],
              [ 1,    1,    2,    1,    1],
              [ 0,    2,    2,    2,    2],
              [ 2,    1,    2,    0,    0],
              [ 2,    0,    4,   -2,    1],
              [ 4,    1,    2,   -2,    2],
              [-1,   -1,    0,    6,    0],
              [-3,   -1,    2,    6,    2],
              [-1,    0,    0,    6,    1],
              [-3,    0,    2,    6,    1],
              [ 1,   -1,    0,    4,    1],
              [ 1,   -1,    0,    4,    0],
              [-2,    0,    2,    5,    2],
              [ 1,   -2,    2,    2,    1],
              [ 3,   -1,    0,    2,    0],
              [ 1,   -1,    2,    2,    0],
              [ 0,    0,    2,    3,    1],
              [-1,    1,    2,    4,    1],
              [ 0,    1,    2,    3,    2],
              [-1,    0,    4,    2,    1],
              [ 2,    0,    2,    1,    1],
              [ 5,    0,    0,    0,    0],
              [ 2,    1,    2,    1,    2],
              [ 1,    0,    4,    0,    1],
              [ 3,    1,    2,    0,    1],
              [ 3,    0,    4,   -2,    2],
              [-2,   -1,    2,    6,    2],
              [ 0,    0,    0,    6,    0],
              [ 0,   -2,    2,    4,    2],
              [-2,    0,    2,    6,    1],
              [ 2,    0,    0,    4,    1],
              [ 2,    0,    0,    4,    0],
              [ 2,   -2,    2,    2,    2],
              [ 0,    0,    2,    4,    0],
              [ 1,    0,    2,    3,    2],
              [ 4,    0,    0,    2,    0],
              [ 2,    0,    2,    2,    0],
              [ 0,    0,    4,    2,    2],
              [ 4,   -1,    2,    0,    2],
              [ 3,    0,    2,    1,    2],
              [ 2,    1,    2,    2,    1],
              [ 4,    1,    2,    0,    2],
              [-1,   -1,    2,    6,    2],
              [-1,    0,    2,    6,    1],
              [ 1,   -1,    2,    4,    1],
              [ 1,    1,    2,    4,    2],
              [ 3,    1,    2,    2,    2],
              [ 5,    0,    2,    0,    1],
              [ 2,   -1,    2,    4,    2],
              [ 2,    0,    2,    4,    1]])

    # *
    # *  Luni-Solar nutation coefficients, unit 1e-7 arcsec
    # *  longitude (sin, t*sin, cos), obliquity (cos, t*cos, sin)
    # *

    CLS = np.asarray([
      [-172064161.0, -174666.0,  33386.0, 92052331.0,  9086.0, 15377.0],
       [-13170906.0,   -1675.0, -13696.0,  5730336.0, -3015.0, -4587.0],
        [-2276413.0,    -234.0,   2796.0,   978459.0,  -485.0,  1374.0],
         [2074554.0,     207.0,   -698.0,  -897492.0,   470.0,  -291.0],
         [1475877.0,   -3633.0,  11817.0,    73871.0,  -184.0, -1924.0],
         [-516821.0,    1226.0,   -524.0,   224386.0,  -677.0,  -174.0],
         [ 711159.0,      73.0,   -872.0,    -6750.0,     0.0,   358.0],
         [-387298.0,    -367.0,    380.0,   200728.0,    18.0,   318.0],
         [-301461.0,     -36.0,    816.0,   129025.0,   -63.0,   367.0],
         [ 215829.0,    -494.0,    111.0,   -95929.0,   299.0,   132.0],
         [ 128227.0,     137.0,    181.0,   -68982.0,    -9.0,    39.0],
         [ 123457.0,      11.0,     19.0,   -53311.0,    32.0,    -4.0],
         [ 156994.0,      10.0,   -168.0,    -1235.0,     0.0,    82.0],
         [  63110.0,      63.0,     27.0,   -33228.0,     0.0,    -9.0],
         [ -57976.0,     -63.0,   -189.0,    31429.0,     0.0,   -75.0],
         [ -59641.0,     -11.0,    149.0,    25543.0,   -11.0,    66.0],
         [ -51613.0,     -42.0,    129.0,    26366.0,     0.0,    78.0],
         [  45893.0,      50.0,     31.0,   -24236.0,   -10.0,    20.0],
         [  63384.0,      11.0,   -150.0,    -1220.0,     0.0,    29.0],
         [ -38571.0,      -1.0,    158.0,    16452.0,   -11.0,    68.0],
         [  32481.0,       0.0,      0.0,   -13870.0,     0.0,     0.0],
         [ -47722.0,       0.0,    -18.0,      477.0,     0.0,   -25.0],
         [ -31046.0,      -1.0,    131.0,    13238.0,   -11.0,    59.0],
         [  28593.0,       0.0,     -1.0,   -12338.0,    10.0,    -3.0],
         [  20441.0,      21.0,     10.0,   -10758.0,     0.0,    -3.0],
         [  29243.0,       0.0,    -74.0,     -609.0,     0.0,    13.0],
         [  25887.0,       0.0,    -66.0,     -550.0,     0.0,    11.0],
         [ -14053.0,     -25.0,     79.0,     8551.0,    -2.0,   -45.0],
         [  15164.0,      10.0,     11.0,    -8001.0,     0.0,    -1.0],
         [ -15794.0,      72.0,    -16.0,     6850.0,   -42.0,    -5.0],
         [  21783.0,       0.0,     13.0,     -167.0,     0.0,    13.0],
         [ -12873.0,     -10.0,    -37.0,     6953.0,     0.0,   -14.0],
         [ -12654.0,      11.0,     63.0,     6415.0,     0.0,    26.0],
         [ -10204.0,       0.0,     25.0,     5222.0,     0.0,    15.0],
         [  16707.0,     -85.0,    -10.0,      168.0,    -1.0,    10.0],
         [  -7691.0,       0.0,     44.0,     3268.0,     0.0,    19.0],
         [ -11024.0,       0.0,    -14.0,      104.0,     0.0,     2.0],
         [   7566.0,     -21.0,    -11.0,    -3250.0,     0.0,    -5.0],
         [  -6637.0,     -11.0,     25.0,     3353.0,     0.0,    14.0],
         [  -7141.0,      21.0,      8.0,     3070.0,     0.0,     4.0],
         [  -6302.0,     -11.0,      2.0,     3272.0,     0.0,     4.0],
         [   5800.0,      10.0,      2.0,    -3045.0,     0.0,    -1.0],
         [   6443.0,       0.0,     -7.0,    -2768.0,     0.0,    -4.0],
         [  -5774.0,     -11.0,    -15.0,     3041.0,     0.0,    -5.0],
         [  -5350.0,       0.0,     21.0,     2695.0,     0.0,    12.0],
         [  -4752.0,     -11.0,     -3.0,     2719.0,     0.0,    -3.0],
         [  -4940.0,     -11.0,    -21.0,     2720.0,     0.0,    -9.0],
         [   7350.0,       0.0,     -8.0,      -51.0,     0.0,     4.0],
         [   4065.0,       0.0,      6.0,    -2206.0,     0.0,     1.0],
         [   6579.0,       0.0,    -24.0,     -199.0,     0.0,     2.0],
         [   3579.0,       0.0,      5.0,    -1900.0,     0.0,     1.0],
         [   4725.0,       0.0,     -6.0,      -41.0,     0.0,     3.0],
         [  -3075.0,       0.0,     -2.0,     1313.0,     0.0,    -1.0],
         [  -2904.0,       0.0,     15.0,     1233.0,     0.0,     7.0],
         [   4348.0,       0.0,    -10.0,      -81.0,     0.0,     2.0],
         [  -2878.0,       0.0,      8.0,     1232.0,     0.0,     4.0],
         [  -4230.0,       0.0,      5.0,      -20.0,     0.0,    -2.0],
         [  -2819.0,       0.0,      7.0,     1207.0,     0.0,     3.0],
         [  -4056.0,       0.0,      5.0,       40.0,     0.0,    -2.0],
         [  -2647.0,       0.0,     11.0,     1129.0,     0.0,     5.0],
         [  -2294.0,       0.0,    -10.0,     1266.0,     0.0,    -4.0],
         [   2481.0,       0.0,     -7.0,    -1062.0,     0.0,    -3.0],
         [   2179.0,       0.0,     -2.0,    -1129.0,     0.0,    -2.0],
         [   3276.0,       0.0,      1.0,       -9.0,     0.0,     0.0],
         [  -3389.0,       0.0,      5.0,       35.0,     0.0,    -2.0],
         [   3339.0,       0.0,    -13.0,     -107.0,     0.0,     1.0],
         [  -1987.0,       0.0,     -6.0,     1073.0,     0.0,    -2.0],
         [  -1981.0,       0.0,      0.0,      854.0,     0.0,     0.0],
         [   4026.0,       0.0,   -353.0,     -553.0,     0.0,  -139.0],
         [   1660.0,       0.0,     -5.0,     -710.0,     0.0,    -2.0],
         [  -1521.0,       0.0,      9.0,      647.0,     0.0,     4.0],
         [   1314.0,       0.0,      0.0,     -700.0,     0.0,     0.0],
         [  -1283.0,       0.0,      0.0,      672.0,     0.0,     0.0],
         [  -1331.0,       0.0,      8.0,      663.0,     0.0,     4.0],
         [   1383.0,       0.0,     -2.0,     -594.0,     0.0,    -2.0],
         [   1405.0,       0.0,      4.0,     -610.0,     0.0,     2.0],
         [   1290.0,       0.0,      0.0,     -556.0,     0.0,     0.0],
         [  -1214.0,       0.0,      5.0,      518.0,     0.0,     2.0],
         [   1146.0,       0.0,     -3.0,     -490.0,     0.0,    -1.0],
         [   1019.0,       0.0,     -1.0,     -527.0,     0.0,    -1.0],
         [  -1100.0,       0.0,      9.0,      465.0,     0.0,     4.0],
         [   -970.0,       0.0,      2.0,      496.0,     0.0,     1.0],
         [   1575.0,       0.0,     -6.0,      -50.0,     0.0,     0.0],
         [    934.0,       0.0,     -3.0,     -399.0,     0.0,    -1.0],
         [    922.0,       0.0,     -1.0,     -395.0,     0.0,    -1.0],
         [    815.0,       0.0,     -1.0,     -422.0,     0.0,    -1.0],
         [    834.0,       0.0,      2.0,     -440.0,     0.0,     1.0],
         [   1248.0,       0.0,      0.0,     -170.0,     0.0,     1.0],
         [   1338.0,       0.0,     -5.0,      -39.0,     0.0,     0.0],
         [    716.0,       0.0,     -2.0,     -389.0,     0.0,    -1.0],
         [   1282.0,       0.0,     -3.0,      -23.0,     0.0,     1.0],
         [    742.0,       0.0,      1.0,     -391.0,     0.0,     0.0],
         [   1020.0,       0.0,    -25.0,     -495.0,     0.0,   -10.0],
         [    715.0,       0.0,     -4.0,     -326.0,     0.0,     2.0],
         [   -666.0,       0.0,     -3.0,      369.0,     0.0,    -1.0],
         [   -667.0,       0.0,      1.0,      346.0,     0.0,     1.0],
         [   -704.0,       0.0,      0.0,      304.0,     0.0,     0.0],
         [   -694.0,       0.0,      5.0,      294.0,     0.0,     2.0],
         [  -1014.0,       0.0,     -1.0,        4.0,     0.0,    -1.0],
         [   -585.0,       0.0,     -2.0,      316.0,     0.0,    -1.0],
         [   -949.0,       0.0,      1.0,        8.0,     0.0,    -1.0],
         [   -595.0,       0.0,      0.0,      258.0,     0.0,     0.0],
         [    528.0,       0.0,      0.0,     -279.0,     0.0,     0.0],
         [   -590.0,       0.0,      4.0,      252.0,     0.0,     2.0],
         [    570.0,       0.0,     -2.0,     -244.0,     0.0,    -1.0],
         [   -502.0,       0.0,      3.0,      250.0,     0.0,     2.0],
         [   -875.0,       0.0,      1.0,       29.0,     0.0,     0.0],
         [   -492.0,       0.0,     -3.0,      275.0,     0.0,    -1.0],
         [    535.0,       0.0,     -2.0,     -228.0,     0.0,    -1.0],
         [   -467.0,       0.0,      1.0,      240.0,     0.0,     1.0],
         [    591.0,       0.0,      0.0,     -253.0,     0.0,     0.0],
         [   -453.0,       0.0,     -1.0,      244.0,     0.0,    -1.0],
         [    766.0,       0.0,      1.0,        9.0,     0.0,     0.0],
         [   -446.0,       0.0,      2.0,      225.0,     0.0,     1.0],
         [   -488.0,       0.0,      2.0,      207.0,     0.0,     1.0],
         [   -468.0,       0.0,      0.0,      201.0,     0.0,     0.0],
         [   -421.0,       0.0,      1.0,      216.0,     0.0,     1.0],
         [    463.0,       0.0,      0.0,     -200.0,     0.0,     0.0],
         [   -673.0,       0.0,      2.0,       14.0,     0.0,     0.0],
         [    658.0,       0.0,      0.0,       -2.0,     0.0,     0.0],
         [   -438.0,       0.0,      0.0,      188.0,     0.0,     0.0],
         [   -390.0,       0.0,      0.0,      205.0,     0.0,     0.0],
         [    639.0,     -11.0,     -2.0,      -19.0,     0.0,     0.0],
         [    412.0,       0.0,     -2.0,     -176.0,     0.0,    -1.0],
         [   -361.0,       0.0,      0.0,      189.0,     0.0,     0.0],
         [    360.0,       0.0,     -1.0,     -185.0,     0.0,    -1.0],
         [    588.0,       0.0,     -3.0,      -24.0,     0.0,     0.0],
         [   -578.0,       0.0,      1.0,        5.0,     0.0,     0.0],
         [   -396.0,       0.0,      0.0,      171.0,     0.0,     0.0],
         [    565.0,       0.0,     -1.0,       -6.0,     0.0,     0.0],
         [   -335.0,       0.0,     -1.0,      184.0,     0.0,    -1.0],
         [    357.0,       0.0,      1.0,     -154.0,     0.0,     0.0],
         [    321.0,       0.0,      1.0,     -174.0,     0.0,     0.0],
         [   -301.0,       0.0,     -1.0,      162.0,     0.0,     0.0],
         [   -334.0,       0.0,      0.0,      144.0,     0.0,     0.0],
         [    493.0,       0.0,     -2.0,      -15.0,     0.0,     0.0],
         [    494.0,       0.0,     -2.0,      -19.0,     0.0,     0.0],
         [    337.0,       0.0,     -1.0,     -143.0,     0.0,    -1.0],
         [    280.0,       0.0,     -1.0,     -144.0,     0.0,     0.0],
         [    309.0,       0.0,      1.0,     -134.0,     0.0,     0.0],
         [   -263.0,       0.0,      2.0,      131.0,     0.0,     1.0],
         [    253.0,       0.0,      1.0,     -138.0,     0.0,     0.0],
         [    245.0,       0.0,      0.0,     -128.0,     0.0,     0.0],
         [    416.0,       0.0,     -2.0,      -17.0,     0.0,     0.0],
         [   -229.0,       0.0,      0.0,      128.0,     0.0,     0.0],
         [    231.0,       0.0,      0.0,     -120.0,     0.0,     0.0],
         [   -259.0,       0.0,      2.0,      109.0,     0.0,     1.0],
         [    375.0,       0.0,     -1.0,       -8.0,     0.0,     0.0],
         [    252.0,       0.0,      0.0,     -108.0,     0.0,     0.0],
         [   -245.0,       0.0,      1.0,      104.0,     0.0,     0.0],
         [    243.0,       0.0,     -1.0,     -104.0,     0.0,     0.0],
         [    208.0,       0.0,      1.0,     -112.0,     0.0,     0.0],
         [    199.0,       0.0,      0.0,     -102.0,     0.0,     0.0],
         [   -208.0,       0.0,      1.0,      105.0,     0.0,     0.0],
         [    335.0,       0.0,     -2.0,      -14.0,     0.0,     0.0],
         [   -325.0,       0.0,      1.0,        7.0,     0.0,     0.0],
         [   -187.0,       0.0,      0.0,       96.0,     0.0,     0.0],
         [    197.0,       0.0,     -1.0,     -100.0,     0.0,     0.0],
         [   -192.0,       0.0,      2.0,       94.0,     0.0,     1.0],
         [   -188.0,       0.0,      0.0,       83.0,     0.0,     0.0],
         [    276.0,       0.0,      0.0,       -2.0,     0.0,     0.0],
         [   -286.0,       0.0,      1.0,        6.0,     0.0,     0.0],
         [    186.0,       0.0,     -1.0,      -79.0,     0.0,     0.0],
         [   -219.0,       0.0,      0.0,       43.0,     0.0,     0.0],
         [    276.0,       0.0,      0.0,        2.0,     0.0,     0.0],
         [   -153.0,       0.0,     -1.0,       84.0,     0.0,     0.0],
         [   -156.0,       0.0,      0.0,       81.0,     0.0,     0.0],
         [   -154.0,       0.0,      1.0,       78.0,     0.0,     0.0],
         [   -174.0,       0.0,      1.0,       75.0,     0.0,     0.0],
         [   -163.0,       0.0,      2.0,       69.0,     0.0,     1.0],
         [   -228.0,       0.0,      0.0,        1.0,     0.0,     0.0],
         [     91.0,       0.0,     -4.0,      -54.0,     0.0,    -2.0],
         [    175.0,       0.0,      0.0,      -75.0,     0.0,     0.0],
         [   -159.0,       0.0,      0.0,       69.0,     0.0,     0.0],
         [    141.0,       0.0,      0.0,      -72.0,     0.0,     0.0],
         [    147.0,       0.0,      0.0,      -75.0,     0.0,     0.0],
         [   -132.0,       0.0,      0.0,       69.0,     0.0,     0.0],
         [    159.0,       0.0,    -28.0,      -54.0,     0.0,    11.0],
         [    213.0,       0.0,      0.0,       -4.0,     0.0,     0.0],
         [    123.0,       0.0,      0.0,      -64.0,     0.0,     0.0],
         [   -118.0,       0.0,     -1.0,       66.0,     0.0,     0.0],
         [    144.0,       0.0,     -1.0,      -61.0,     0.0,     0.0],
         [   -121.0,       0.0,      1.0,       60.0,     0.0,     0.0],
         [   -134.0,       0.0,      1.0,       56.0,     0.0,     1.0],
         [   -105.0,       0.0,      0.0,       57.0,     0.0,     0.0],
         [   -102.0,       0.0,      0.0,       56.0,     0.0,     0.0],
         [    120.0,       0.0,      0.0,      -52.0,     0.0,     0.0],
         [    101.0,       0.0,      0.0,      -54.0,     0.0,     0.0],
         [   -113.0,       0.0,      0.0,       59.0,     0.0,     0.0],
         [   -106.0,       0.0,      0.0,       61.0,     0.0,     0.0],
         [   -129.0,       0.0,      1.0,       55.0,     0.0,     0.0],
         [   -114.0,       0.0,      0.0,       57.0,     0.0,     0.0],
         [    113.0,       0.0,     -1.0,      -49.0,     0.0,     0.0],
         [   -102.0,       0.0,      0.0,       44.0,     0.0,     0.0],
         [    -94.0,       0.0,      0.0,       51.0,     0.0,     0.0],
         [   -100.0,       0.0,     -1.0,       56.0,     0.0,     0.0],
         [     87.0,       0.0,      0.0,      -47.0,     0.0,     0.0],
         [    161.0,       0.0,      0.0,       -1.0,     0.0,     0.0],
         [     96.0,       0.0,      0.0,      -50.0,     0.0,     0.0],
         [    151.0,       0.0,     -1.0,       -5.0,     0.0,     0.0],
         [   -104.0,       0.0,      0.0,       44.0,     0.0,     0.0],
         [   -110.0,       0.0,      0.0,       48.0,     0.0,     0.0],
         [   -100.0,       0.0,      1.0,       50.0,     0.0,     0.0],
         [     92.0,       0.0,     -5.0,       12.0,     0.0,    -2.0],
         [     82.0,       0.0,      0.0,      -45.0,     0.0,     0.0],
         [     82.0,       0.0,      0.0,      -45.0,     0.0,     0.0],
         [    -78.0,       0.0,      0.0,       41.0,     0.0,     0.0],
         [    -77.0,       0.0,      0.0,       43.0,     0.0,     0.0],
         [      2.0,       0.0,      0.0,       54.0,     0.0,     0.0],
         [     94.0,       0.0,      0.0,      -40.0,     0.0,     0.0],
         [    -93.0,       0.0,      0.0,       40.0,     0.0,     0.0],
         [    -83.0,       0.0,     10.0,       40.0,     0.0,    -2.0],
         [     83.0,       0.0,      0.0,      -36.0,     0.0,     0.0],
         [    -91.0,       0.0,      0.0,       39.0,     0.0,     0.0],
         [    128.0,       0.0,      0.0,       -1.0,     0.0,     0.0],
         [    -79.0,       0.0,      0.0,       34.0,     0.0,     0.0],
         [    -83.0,       0.0,      0.0,       47.0,     0.0,     0.0],
         [     84.0,       0.0,      0.0,      -44.0,     0.0,     0.0],
         [     83.0,       0.0,      0.0,      -43.0,     0.0,     0.0],
         [     91.0,       0.0,      0.0,      -39.0,     0.0,     0.0],
         [    -77.0,       0.0,      0.0,       39.0,     0.0,     0.0],
         [     84.0,       0.0,      0.0,      -43.0,     0.0,     0.0],
         [    -92.0,       0.0,      1.0,       39.0,     0.0,     0.0],
         [    -92.0,       0.0,      1.0,       39.0,     0.0,     0.0],
         [    -94.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [     68.0,       0.0,      0.0,      -36.0,     0.0,     0.0],
         [    -61.0,       0.0,      0.0,       32.0,     0.0,     0.0],
         [     71.0,       0.0,      0.0,      -31.0,     0.0,     0.0],
         [     62.0,       0.0,      0.0,      -34.0,     0.0,     0.0],
         [    -63.0,       0.0,      0.0,       33.0,     0.0,     0.0],
         [    -73.0,       0.0,      0.0,       32.0,     0.0,     0.0],
         [    115.0,       0.0,      0.0,       -2.0,     0.0,     0.0],
         [   -103.0,       0.0,      0.0,        2.0,     0.0,     0.0],
         [     63.0,       0.0,      0.0,      -28.0,     0.0,     0.0],
         [     74.0,       0.0,      0.0,      -32.0,     0.0,     0.0],
         [   -103.0,       0.0,     -3.0,        3.0,     0.0,    -1.0],
         [    -69.0,       0.0,      0.0,       30.0,     0.0,     0.0],
         [     57.0,       0.0,      0.0,      -29.0,     0.0,     0.0],
         [     94.0,       0.0,      0.0,       -4.0,     0.0,     0.0],
         [     64.0,       0.0,      0.0,      -33.0,     0.0,     0.0],
         [    -63.0,       0.0,      0.0,       26.0,     0.0,     0.0],
         [    -38.0,       0.0,      0.0,       20.0,     0.0,     0.0],
         [    -43.0,       0.0,      0.0,       24.0,     0.0,     0.0],
         [    -45.0,       0.0,      0.0,       23.0,     0.0,     0.0],
         [     47.0,       0.0,      0.0,      -24.0,     0.0,     0.0],
         [    -48.0,       0.0,      0.0,       25.0,     0.0,     0.0],
         [     45.0,       0.0,      0.0,      -26.0,     0.0,     0.0],
         [     56.0,       0.0,      0.0,      -25.0,     0.0,     0.0],
         [     88.0,       0.0,      0.0,        2.0,     0.0,     0.0],
         [    -75.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [     85.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [     49.0,       0.0,      0.0,      -26.0,     0.0,     0.0],
         [    -74.0,       0.0,     -3.0,       -1.0,     0.0,    -1.0],
         [    -39.0,       0.0,      0.0,       21.0,     0.0,     0.0],
         [     45.0,       0.0,      0.0,      -20.0,     0.0,     0.0],
         [     51.0,       0.0,      0.0,      -22.0,     0.0,     0.0],
         [    -40.0,       0.0,      0.0,       21.0,     0.0,     0.0],
         [     41.0,       0.0,      0.0,      -21.0,     0.0,     0.0],
         [    -42.0,       0.0,      0.0,       24.0,     0.0,     0.0],
         [    -51.0,       0.0,      0.0,       22.0,     0.0,     0.0],
         [    -42.0,       0.0,      0.0,       22.0,     0.0,     0.0],
         [     39.0,       0.0,      0.0,      -21.0,     0.0,     0.0],
         [     46.0,       0.0,      0.0,      -18.0,     0.0,     0.0],
         [    -53.0,       0.0,      0.0,       22.0,     0.0,     0.0],
         [     82.0,       0.0,      0.0,       -4.0,     0.0,     0.0],
         [     81.0,       0.0,     -1.0,       -4.0,     0.0,     0.0],
         [     47.0,       0.0,      0.0,      -19.0,     0.0,     0.0],
         [     53.0,       0.0,      0.0,      -23.0,     0.0,     0.0],
         [    -45.0,       0.0,      0.0,       22.0,     0.0,     0.0],
         [    -44.0,       0.0,      0.0,       -2.0,     0.0,     0.0],
         [    -33.0,       0.0,      0.0,       16.0,     0.0,     0.0],
         [    -61.0,       0.0,      0.0,        1.0,     0.0,     0.0],
         [     28.0,       0.0,      0.0,      -15.0,     0.0,     0.0],
         [    -38.0,       0.0,      0.0,       19.0,     0.0,     0.0],
         [    -33.0,       0.0,      0.0,       21.0,     0.0,     0.0],
         [    -60.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [     48.0,       0.0,      0.0,      -10.0,     0.0,     0.0],
         [     27.0,       0.0,      0.0,      -14.0,     0.0,     0.0],
         [     38.0,       0.0,      0.0,      -20.0,     0.0,     0.0],
         [     31.0,       0.0,      0.0,      -13.0,     0.0,     0.0],
         [    -29.0,       0.0,      0.0,       15.0,     0.0,     0.0],
         [     28.0,       0.0,      0.0,      -15.0,     0.0,     0.0],
         [    -32.0,       0.0,      0.0,       15.0,     0.0,     0.0],
         [     45.0,       0.0,      0.0,       -8.0,     0.0,     0.0],
         [    -44.0,       0.0,      0.0,       19.0,     0.0,     0.0],
         [     28.0,       0.0,      0.0,      -15.0,     0.0,     0.0],
         [    -51.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [    -36.0,       0.0,      0.0,       20.0,     0.0,     0.0],
         [     44.0,       0.0,      0.0,      -19.0,     0.0,     0.0],
         [     26.0,       0.0,      0.0,      -14.0,     0.0,     0.0],
         [    -60.0,       0.0,      0.0,        2.0,     0.0,     0.0],
         [     35.0,       0.0,      0.0,      -18.0,     0.0,     0.0],
         [    -27.0,       0.0,      0.0,       11.0,     0.0,     0.0],
         [     47.0,       0.0,      0.0,       -1.0,     0.0,     0.0],
         [     36.0,       0.0,      0.0,      -15.0,     0.0,     0.0],
         [    -36.0,       0.0,      0.0,       20.0,     0.0,     0.0],
         [    -35.0,       0.0,      0.0,       19.0,     0.0,     0.0],
         [    -37.0,       0.0,      0.0,       19.0,     0.0,     0.0],
         [     32.0,       0.0,      0.0,      -16.0,     0.0,     0.0],
         [     35.0,       0.0,      0.0,      -14.0,     0.0,     0.0],
         [     32.0,       0.0,      0.0,      -13.0,     0.0,     0.0],
         [     65.0,       0.0,      0.0,       -2.0,     0.0,     0.0],
         [     47.0,       0.0,      0.0,       -1.0,     0.0,     0.0],
         [     32.0,       0.0,      0.0,      -16.0,     0.0,     0.0],
         [     37.0,       0.0,      0.0,      -16.0,     0.0,     0.0],
         [    -30.0,       0.0,      0.0,       15.0,     0.0,     0.0],
         [    -32.0,       0.0,      0.0,       16.0,     0.0,     0.0],
         [    -31.0,       0.0,      0.0,       13.0,     0.0,     0.0],
         [     37.0,       0.0,      0.0,      -16.0,     0.0,     0.0],
         [     31.0,       0.0,      0.0,      -13.0,     0.0,     0.0],
         [     49.0,       0.0,      0.0,       -2.0,     0.0,     0.0],
         [     32.0,       0.0,      0.0,      -13.0,     0.0,     0.0],
         [     23.0,       0.0,      0.0,      -12.0,     0.0,     0.0],
         [    -43.0,       0.0,      0.0,       18.0,     0.0,     0.0],
         [     26.0,       0.0,      0.0,      -11.0,     0.0,     0.0],
         [    -32.0,       0.0,      0.0,       14.0,     0.0,     0.0],
         [    -29.0,       0.0,      0.0,       14.0,     0.0,     0.0],
         [    -27.0,       0.0,      0.0,       12.0,     0.0,     0.0],
         [     30.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [    -11.0,       0.0,      0.0,        5.0,     0.0,     0.0],
         [    -21.0,       0.0,      0.0,       10.0,     0.0,     0.0],
         [    -34.0,       0.0,      0.0,       15.0,     0.0,     0.0],
         [    -10.0,       0.0,      0.0,        6.0,     0.0,     0.0],
         [    -36.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [     -9.0,       0.0,      0.0,        4.0,     0.0,     0.0],
         [    -12.0,       0.0,      0.0,        5.0,     0.0,     0.0],
         [    -21.0,       0.0,      0.0,        5.0,     0.0,     0.0],
         [    -29.0,       0.0,      0.0,       -1.0,     0.0,     0.0],
         [    -15.0,       0.0,      0.0,        3.0,     0.0,     0.0],
         [    -20.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [     28.0,       0.0,      0.0,        0.0,     0.0,    -2.0],
         [     17.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [    -22.0,       0.0,      0.0,       12.0,     0.0,     0.0],
         [    -14.0,       0.0,      0.0,        7.0,     0.0,     0.0],
         [     24.0,       0.0,      0.0,      -11.0,     0.0,     0.0],
         [     11.0,       0.0,      0.0,       -6.0,     0.0,     0.0],
         [     14.0,       0.0,      0.0,       -6.0,     0.0,     0.0],
         [     24.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [     18.0,       0.0,      0.0,       -8.0,     0.0,     0.0],
         [    -38.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [    -31.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [    -16.0,       0.0,      0.0,        8.0,     0.0,     0.0],
         [     29.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [    -18.0,       0.0,      0.0,       10.0,     0.0,     0.0],
         [    -10.0,       0.0,      0.0,        5.0,     0.0,     0.0],
         [    -17.0,       0.0,      0.0,       10.0,     0.0,     0.0],
         [      9.0,       0.0,      0.0,       -4.0,     0.0,     0.0],
         [     16.0,       0.0,      0.0,       -6.0,     0.0,     0.0],
         [     22.0,       0.0,      0.0,      -12.0,     0.0,     0.0],
         [     20.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [    -13.0,       0.0,      0.0,        6.0,     0.0,     0.0],
         [    -17.0,       0.0,      0.0,        9.0,     0.0,     0.0],
         [    -14.0,       0.0,      0.0,        8.0,     0.0,     0.0],
         [      0.0,       0.0,      0.0,       -7.0,     0.0,     0.0],
         [     14.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [     19.0,       0.0,      0.0,      -10.0,     0.0,     0.0],
         [    -34.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [    -20.0,       0.0,      0.0,        8.0,     0.0,     0.0],
         [      9.0,       0.0,      0.0,       -5.0,     0.0,     0.0],
         [    -18.0,       0.0,      0.0,        7.0,     0.0,     0.0],
         [     13.0,       0.0,      0.0,       -6.0,     0.0,     0.0],
         [     17.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [    -12.0,       0.0,      0.0,        5.0,     0.0,     0.0],
         [     15.0,       0.0,      0.0,       -8.0,     0.0,     0.0],
         [    -11.0,       0.0,      0.0,        3.0,     0.0,     0.0],
         [     13.0,       0.0,      0.0,       -5.0,     0.0,     0.0],
         [    -18.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [    -35.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [      9.0,       0.0,      0.0,       -4.0,     0.0,     0.0],
         [    -19.0,       0.0,      0.0,       10.0,     0.0,     0.0],
         [    -26.0,       0.0,      0.0,       11.0,     0.0,     0.0],
         [      8.0,       0.0,      0.0,       -4.0,     0.0,     0.0],
         [    -10.0,       0.0,      0.0,        4.0,     0.0,     0.0],
         [     10.0,       0.0,      0.0,       -6.0,     0.0,     0.0],
         [    -21.0,       0.0,      0.0,        9.0,     0.0,     0.0],
         [    -15.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [      9.0,       0.0,      0.0,       -5.0,     0.0,     0.0],
         [    -29.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [    -19.0,       0.0,      0.0,       10.0,     0.0,     0.0],
         [     12.0,       0.0,      0.0,       -5.0,     0.0,     0.0],
         [     22.0,       0.0,      0.0,       -9.0,     0.0,     0.0],
         [    -10.0,       0.0,      0.0,        5.0,     0.0,     0.0],
         [    -20.0,       0.0,      0.0,       11.0,     0.0,     0.0],
         [    -20.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [    -17.0,       0.0,      0.0,        7.0,     0.0,     0.0],
         [     15.0,       0.0,      0.0,       -3.0,     0.0,     0.0],
         [      8.0,       0.0,      0.0,       -4.0,     0.0,     0.0],
         [     14.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [    -12.0,       0.0,      0.0,        6.0,     0.0,     0.0],
         [     25.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [    -13.0,       0.0,      0.0,        6.0,     0.0,     0.0],
         [    -14.0,       0.0,      0.0,        8.0,     0.0,     0.0],
         [     13.0,       0.0,      0.0,       -5.0,     0.0,     0.0],
         [    -17.0,       0.0,      0.0,        9.0,     0.0,     0.0],
         [    -12.0,       0.0,      0.0,        6.0,     0.0,     0.0],
         [    -10.0,       0.0,      0.0,        5.0,     0.0,     0.0],
         [     10.0,       0.0,      0.0,       -6.0,     0.0,     0.0],
         [    -15.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [    -22.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [     28.0,       0.0,      0.0,       -1.0,     0.0,     0.0],
         [     15.0,       0.0,      0.0,       -7.0,     0.0,     0.0],
         [     23.0,       0.0,      0.0,      -10.0,     0.0,     0.0],
         [     12.0,       0.0,      0.0,       -5.0,     0.0,     0.0],
         [     29.0,       0.0,      0.0,       -1.0,     0.0,     0.0],
         [    -25.0,       0.0,      0.0,        1.0,     0.0,     0.0],
         [     22.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [    -18.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [     15.0,       0.0,      0.0,        3.0,     0.0,     0.0],
         [    -23.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [     12.0,       0.0,      0.0,       -5.0,     0.0,     0.0],
         [     -8.0,       0.0,      0.0,        4.0,     0.0,     0.0],
         [    -19.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [    -10.0,       0.0,      0.0,        4.0,     0.0,     0.0],
         [     21.0,       0.0,      0.0,       -9.0,     0.0,     0.0],
         [     23.0,       0.0,      0.0,       -1.0,     0.0,     0.0],
         [    -16.0,       0.0,      0.0,        8.0,     0.0,     0.0],
         [    -19.0,       0.0,      0.0,        9.0,     0.0,     0.0],
         [    -22.0,       0.0,      0.0,       10.0,     0.0,     0.0],
         [     27.0,       0.0,      0.0,       -1.0,     0.0,     0.0],
         [     16.0,       0.0,      0.0,       -8.0,     0.0,     0.0],
         [     19.0,       0.0,      0.0,       -8.0,     0.0,     0.0],
         [      9.0,       0.0,      0.0,       -4.0,     0.0,     0.0],
         [     -9.0,       0.0,      0.0,        4.0,     0.0,     0.0],
         [     -9.0,       0.0,      0.0,        4.0,     0.0,     0.0],
         [     -8.0,       0.0,      0.0,        4.0,     0.0,     0.0],
         [     18.0,       0.0,      0.0,       -9.0,     0.0,     0.0],
         [     16.0,       0.0,      0.0,       -1.0,     0.0,     0.0],
         [    -10.0,       0.0,      0.0,        4.0,     0.0,     0.0],
         [    -23.0,       0.0,      0.0,        9.0,     0.0,     0.0],
         [     16.0,       0.0,      0.0,       -1.0,     0.0,     0.0],
         [    -12.0,       0.0,      0.0,        6.0,     0.0,     0.0],
         [     -8.0,       0.0,      0.0,        4.0,     0.0,     0.0],
         [     30.0,       0.0,      0.0,       -2.0,     0.0,     0.0],
         [     24.0,       0.0,      0.0,      -10.0,     0.0,     0.0],
         [     10.0,       0.0,      0.0,       -4.0,     0.0,     0.0],
         [    -16.0,       0.0,      0.0,        7.0,     0.0,     0.0],
         [    -16.0,       0.0,      0.0,        7.0,     0.0,     0.0],
         [     17.0,       0.0,      0.0,       -7.0,     0.0,     0.0],
         [    -24.0,       0.0,      0.0,       10.0,     0.0,     0.0],
         [    -12.0,       0.0,      0.0,        5.0,     0.0,     0.0],
         [    -24.0,       0.0,      0.0,       11.0,     0.0,     0.0],
         [    -23.0,       0.0,      0.0,        9.0,     0.0,     0.0],
         [    -13.0,       0.0,      0.0,        5.0,     0.0,     0.0],
         [    -15.0,       0.0,      0.0,        7.0,     0.0,     0.0],
         [      0.0,       0.0,  -1988.0,        0.0,     0.0, -1679.0],
         [      0.0,       0.0,    -63.0,        0.0,     0.0,   -27.0],
         [     -4.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [      0.0,       0.0,      5.0,        0.0,     0.0,     4.0],
         [      5.0,       0.0,      0.0,       -3.0,     0.0,     0.0],
         [      0.0,       0.0,    364.0,        0.0,     0.0,   176.0],
         [      0.0,       0.0,  -1044.0,        0.0,     0.0,  -891.0],
         [     -3.0,       0.0,      0.0,        1.0,     0.0,     0.0],
         [      4.0,       0.0,      0.0,       -2.0,     0.0,     0.0],
         [      0.0,       0.0,    330.0,        0.0,     0.0,     0.0],
         [      5.0,       0.0,      0.0,       -2.0,     0.0,     0.0],
         [      3.0,       0.0,      0.0,       -2.0,     0.0,     0.0],
         [     -3.0,       0.0,      0.0,        1.0,     0.0,     0.0],
         [     -5.0,       0.0,      0.0,        2.0,     0.0,     0.0],
         [      3.0,       0.0,      0.0,       -1.0,     0.0,     0.0],
         [      3.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [      3.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [      0.0,       0.0,      5.0,        0.0,     0.0,     0.0],
         [      0.0,       0.0,      0.0,        1.0,     0.0,     0.0],
         [      4.0,       0.0,      0.0,       -2.0,     0.0,     0.0],
         [      6.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [      5.0,       0.0,      0.0,       -2.0,     0.0,     0.0],
         [     -7.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [    -12.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [      5.0,       0.0,      0.0,       -3.0,     0.0,     0.0],
         [      3.0,       0.0,      0.0,       -1.0,     0.0,     0.0],
         [     -5.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [      3.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [     -7.0,       0.0,      0.0,        3.0,     0.0,     0.0],
         [      7.0,       0.0,      0.0,       -4.0,     0.0,     0.0],
         [      0.0,       0.0,    -12.0,        0.0,     0.0,   -10.0],
         [      4.0,       0.0,      0.0,       -2.0,     0.0,     0.0],
         [      3.0,       0.0,      0.0,       -2.0,     0.0,     0.0],
         [     -3.0,       0.0,      0.0,        2.0,     0.0,     0.0],
         [     -7.0,       0.0,      0.0,        3.0,     0.0,     0.0],
         [     -4.0,       0.0,      0.0,        2.0,     0.0,     0.0],
         [     -3.0,       0.0,      0.0,        1.0,     0.0,     0.0],
         [      0.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [     -3.0,       0.0,      0.0,        1.0,     0.0,     0.0],
         [      7.0,       0.0,      0.0,       -3.0,     0.0,     0.0],
         [     -4.0,       0.0,      0.0,        2.0,     0.0,     0.0],
         [      4.0,       0.0,      0.0,       -2.0,     0.0,     0.0],
         [     -5.0,       0.0,      0.0,        3.0,     0.0,     0.0],
         [      5.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [     -5.0,       0.0,      0.0,        2.0,     0.0,     0.0],
         [      5.0,       0.0,      0.0,       -2.0,     0.0,     0.0],
         [     -8.0,       0.0,      0.0,        3.0,     0.0,     0.0],
         [      9.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [      6.0,       0.0,      0.0,       -3.0,     0.0,     0.0],
         [     -5.0,       0.0,      0.0,        2.0,     0.0,     0.0],
         [      3.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [     -7.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [     -3.0,       0.0,      0.0,        1.0,     0.0,     0.0],
         [      5.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [      3.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [     -3.0,       0.0,      0.0,        2.0,     0.0,     0.0],
         [      4.0,       0.0,      0.0,       -2.0,     0.0,     0.0],
         [      3.0,       0.0,      0.0,       -1.0,     0.0,     0.0],
         [     -5.0,       0.0,      0.0,        2.0,     0.0,     0.0],
         [      4.0,       0.0,      0.0,       -2.0,     0.0,     0.0],
         [      9.0,       0.0,      0.0,       -3.0,     0.0,     0.0],
         [      4.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [      4.0,       0.0,      0.0,       -2.0,     0.0,     0.0],
         [     -3.0,       0.0,      0.0,        2.0,     0.0,     0.0],
         [     -4.0,       0.0,      0.0,        2.0,     0.0,     0.0],
         [      9.0,       0.0,      0.0,       -3.0,     0.0,     0.0],
         [     -4.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [     -4.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [      3.0,       0.0,      0.0,       -2.0,     0.0,     0.0],
         [      8.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [      3.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [     -3.0,       0.0,      0.0,        2.0,     0.0,     0.0],
         [      3.0,       0.0,      0.0,       -1.0,     0.0,     0.0],
         [      3.0,       0.0,      0.0,       -1.0,     0.0,     0.0],
         [     -3.0,       0.0,      0.0,        1.0,     0.0,     0.0],
         [      6.0,       0.0,      0.0,       -3.0,     0.0,     0.0],
         [      3.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [     -3.0,       0.0,      0.0,        1.0,     0.0,     0.0],
         [     -7.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [      9.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [     -3.0,       0.0,      0.0,        2.0,     0.0,     0.0],
         [     -3.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [     -4.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [     -5.0,       0.0,      0.0,        3.0,     0.0,     0.0],
         [    -13.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [     -7.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [     10.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [      3.0,       0.0,      0.0,       -1.0,     0.0,     0.0],
         [     10.0,       0.0,     13.0,        6.0,     0.0,    -5.0],
         [      0.0,       0.0,     30.0,        0.0,     0.0,    14.0],
         [      0.0,       0.0,   -162.0,        0.0,     0.0,  -138.0],
         [      0.0,       0.0,     75.0,        0.0,     0.0,     0.0],
         [     -7.0,       0.0,      0.0,        4.0,     0.0,     0.0],
         [     -4.0,       0.0,      0.0,        2.0,     0.0,     0.0],
         [      4.0,       0.0,      0.0,       -2.0,     0.0,     0.0],
         [      5.0,       0.0,      0.0,       -2.0,     0.0,     0.0],
         [      5.0,       0.0,      0.0,       -3.0,     0.0,     0.0],
         [     -3.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [     -3.0,       0.0,      0.0,        2.0,     0.0,     0.0],
         [     -4.0,       0.0,      0.0,        2.0,     0.0,     0.0],
         [     -5.0,       0.0,      0.0,        2.0,     0.0,     0.0],
         [      6.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [      9.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [      5.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [     -7.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [     -3.0,       0.0,      0.0,        1.0,     0.0,     0.0],
         [     -4.0,       0.0,      0.0,        2.0,     0.0,     0.0],
         [      7.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [     -4.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [      4.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [     -6.0,       0.0,     -3.0,        3.0,     0.0,     1.0],
         [      0.0,       0.0,     -3.0,        0.0,     0.0,    -2.0],
         [     11.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [      3.0,       0.0,      0.0,       -1.0,     0.0,     0.0],
         [     11.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [     -3.0,       0.0,      0.0,        2.0,     0.0,     0.0],
         [     -1.0,       0.0,      3.0,        3.0,     0.0,    -1.0],
         [      4.0,       0.0,      0.0,       -2.0,     0.0,     0.0],
         [      0.0,       0.0,    -13.0,        0.0,     0.0,   -11.0],
         [      3.0,       0.0,      6.0,        0.0,     0.0,     0.0],
         [     -7.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [      5.0,       0.0,      0.0,       -3.0,     0.0,     0.0],
         [     -3.0,       0.0,      0.0,        1.0,     0.0,     0.0],
         [      3.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [      5.0,       0.0,      0.0,       -3.0,     0.0,     0.0],
         [     -7.0,       0.0,      0.0,        3.0,     0.0,     0.0],
         [      8.0,       0.0,      0.0,       -3.0,     0.0,     0.0],
         [     -4.0,       0.0,      0.0,        2.0,     0.0,     0.0],
         [     11.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [     -3.0,       0.0,      0.0,        1.0,     0.0,     0.0],
         [      3.0,       0.0,      0.0,       -1.0,     0.0,     0.0],
         [     -4.0,       0.0,      0.0,        2.0,     0.0,     0.0],
         [      8.0,       0.0,      0.0,       -4.0,     0.0,     0.0],
         [      3.0,       0.0,      0.0,       -1.0,     0.0,     0.0],
         [     11.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [     -6.0,       0.0,      0.0,        3.0,     0.0,     0.0],
         [     -4.0,       0.0,      0.0,        2.0,     0.0,     0.0],
         [     -8.0,       0.0,      0.0,        4.0,     0.0,     0.0],
         [     -7.0,       0.0,      0.0,        3.0,     0.0,     0.0],
         [     -4.0,       0.0,      0.0,        2.0,     0.0,     0.0],
         [      3.0,       0.0,      0.0,       -1.0,     0.0,     0.0],
         [      6.0,       0.0,      0.0,       -3.0,     0.0,     0.0],
         [     -6.0,       0.0,      0.0,        3.0,     0.0,     0.0],
         [      6.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [      6.0,       0.0,      0.0,       -1.0,     0.0,     0.0],
         [      5.0,       0.0,      0.0,       -2.0,     0.0,     0.0],
         [     -5.0,       0.0,      0.0,        2.0,     0.0,     0.0],
         [     -4.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [     -4.0,       0.0,      0.0,        2.0,     0.0,     0.0],
         [      4.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [      6.0,       0.0,      0.0,       -3.0,     0.0,     0.0],
         [     -4.0,       0.0,      0.0,        2.0,     0.0,     0.0],
         [      0.0,       0.0,    -26.0,        0.0,     0.0,   -11.0],
         [      0.0,       0.0,    -10.0,        0.0,     0.0,    -5.0],
         [      5.0,       0.0,      0.0,       -3.0,     0.0,     0.0],
         [    -13.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [      3.0,       0.0,      0.0,       -2.0,     0.0,     0.0],
         [      4.0,       0.0,      0.0,       -2.0,     0.0,     0.0],
         [      7.0,       0.0,      0.0,       -3.0,     0.0,     0.0],
         [      4.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [      5.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [     -3.0,       0.0,      0.0,        2.0,     0.0,     0.0],
         [     -6.0,       0.0,      0.0,        2.0,     0.0,     0.0],
         [     -5.0,       0.0,      0.0,        2.0,     0.0,     0.0],
         [     -7.0,       0.0,      0.0,        3.0,     0.0,     0.0],
         [      5.0,       0.0,      0.0,       -2.0,     0.0,     0.0],
         [     13.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [     -4.0,       0.0,      0.0,        2.0,     0.0,     0.0],
         [     -3.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [      5.0,       0.0,      0.0,       -2.0,     0.0,     0.0],
         [    -11.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [      5.0,       0.0,      0.0,       -2.0,     0.0,     0.0],
         [      4.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [      4.0,       0.0,      0.0,       -2.0,     0.0,     0.0],
         [     -4.0,       0.0,      0.0,        2.0,     0.0,     0.0],
         [      6.0,       0.0,      0.0,       -3.0,     0.0,     0.0],
         [      3.0,       0.0,      0.0,       -2.0,     0.0,     0.0],
         [    -12.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [      4.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [     -3.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [     -4.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [      3.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [      3.0,       0.0,      0.0,       -1.0,     0.0,     0.0],
         [     -3.0,       0.0,      0.0,        1.0,     0.0,     0.0],
         [      0.0,       0.0,     -5.0,        0.0,     0.0,    -2.0],
         [     -7.0,       0.0,      0.0,        4.0,     0.0,     0.0],
         [      6.0,       0.0,      0.0,       -3.0,     0.0,     0.0],
         [     -3.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [      5.0,       0.0,      0.0,       -3.0,     0.0,     0.0],
         [      3.0,       0.0,      0.0,       -1.0,     0.0,     0.0],
         [      3.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [     -3.0,       0.0,      0.0,        1.0,     0.0,     0.0],
         [     -5.0,       0.0,      0.0,        3.0,     0.0,     0.0],
         [     -3.0,       0.0,      0.0,        2.0,     0.0,     0.0],
         [     -3.0,       0.0,      0.0,        2.0,     0.0,     0.0],
         [     12.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [      3.0,       0.0,      0.0,       -1.0,     0.0,     0.0],
         [     -4.0,       0.0,      0.0,        2.0,     0.0,     0.0],
         [      4.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [      6.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [      5.0,       0.0,      0.0,       -3.0,     0.0,     0.0],
         [      4.0,       0.0,      0.0,       -2.0,     0.0,     0.0],
         [     -6.0,       0.0,      0.0,        3.0,     0.0,     0.0],
         [      4.0,       0.0,      0.0,       -2.0,     0.0,     0.0],
         [      6.0,       0.0,      0.0,       -3.0,     0.0,     0.0],
         [      6.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [     -6.0,       0.0,      0.0,        3.0,     0.0,     0.0],
         [      3.0,       0.0,      0.0,       -2.0,     0.0,     0.0],
         [      7.0,       0.0,      0.0,       -4.0,     0.0,     0.0],
         [      4.0,       0.0,      0.0,       -2.0,     0.0,     0.0],
         [     -5.0,       0.0,      0.0,        2.0,     0.0,     0.0],
         [      5.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [     -6.0,       0.0,      0.0,        3.0,     0.0,     0.0],
         [     -6.0,       0.0,      0.0,        3.0,     0.0,     0.0],
         [     -4.0,       0.0,      0.0,        2.0,     0.0,     0.0],
         [     10.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [     -4.0,       0.0,      0.0,        2.0,     0.0,     0.0],
         [      7.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [      7.0,       0.0,      0.0,       -3.0,     0.0,     0.0],
         [      4.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [     11.0,       0.0,      0.0,        0.0,     0.0,     0.0],
         [      5.0,       0.0,      0.0,       -2.0,     0.0,     0.0],
         [     -6.0,       0.0,      0.0,        2.0,     0.0,     0.0],
         [      4.0,       0.0,      0.0,       -2.0,     0.0,     0.0],
         [      3.0,       0.0,      0.0,       -2.0,     0.0,     0.0],
         [      5.0,       0.0,      0.0,       -2.0,     0.0,     0.0],
         [     -4.0,       0.0,      0.0,        2.0,     0.0,     0.0],
         [     -4.0,       0.0,      0.0,        2.0,     0.0,     0.0],
         [     -3.0,       0.0,      0.0,        2.0,     0.0,     0.0],
         [      4.0,       0.0,      0.0,       -2.0,     0.0,     0.0],
         [      3.0,       0.0,      0.0,       -1.0,     0.0,     0.0],
         [     -3.0,       0.0,      0.0,        1.0,     0.0,     0.0],
         [     -3.0,       0.0,      0.0,        1.0,     0.0,     0.0],
         [     -3.0,       0.0,      0.0,        2.0,     0.0,     0.0] ])

# 
#  Planetary argument multipliers
#              L   L'  F   D   Om  Me  Ve  E  Ma  Ju  Sa  Ur  Ne  pre

    NAPL = np.asarray([
             [ 0,  0,  0,  0,  0,  0,  0,  8,-16,  4,  5,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0, -8, 16, -4, -5,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  8,-16,  4,  5,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  2,  2],
             [ 0,  0,  0,  0,  0,  0,  0, -4,  8, -1, -5,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  4, -8,  3,  0,  0,  0,  1],
             [ 0,  0,  1, -1,  1,  0,  0,  3, -8,  3,  0,  0,  0,  0],
             [-1,  0,  0,  0,  0,  0, 10, -3,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  0,  0, -2,  6, -3,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  4, -8,  3,  0,  0,  0,  0],
             [ 0,  0,  1, -1,  1,  0,  0, -5,  8, -3,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0, -4,  8, -3,  0,  0,  0,  1],
             [ 0,  0,  0,  0,  0,  0,  0,  4, -8,  1,  5,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0, -5,  6,  4,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  2, -5,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  2, -5,  0,  0,  1],
             [ 0,  0,  1, -1,  1,  0,  0, -1,  0,  2, -5,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  2, -5,  0,  0,  0],
             [ 0,  0,  1, -1,  1,  0,  0, -1,  0, -2,  5,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  0,  0, -2,  5,  0,  0,  1],
             [ 0,  0,  0,  0,  0,  0,  0,  0,  0, -2,  5,  0,  0,  2],
             [ 2,  0, -1, -1,  0,  0,  0,  3, -7,  0,  0,  0,  0,  0],
             [ 1,  0,  0, -2,  0,  0, 19,-21,  3,  0,  0,  0,  0,  0],
             [ 0,  0,  1, -1,  1,  0,  2, -4,  0, -3,  0,  0,  0,  0],
             [ 1,  0,  0, -1,  1,  0,  0, -1,  0,  2,  0,  0,  0,  0],
             [ 0,  0,  1, -1,  1,  0,  0, -1,  0, -4, 10,  0,  0,  0],
             [-2,  0,  0,  2,  1,  0,  0,  2,  0,  0, -5,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  3, -7,  4,  0,  0,  0,  0,  0],
             [ 0,  0, -1,  1,  0,  0,  0,  1,  0,  1, -1,  0,  0,  0],
             [-2,  0,  0,  2,  1,  0,  0,  2,  0, -2,  0,  0,  0,  0],
             [-1,  0,  0,  0,  0,  0, 18,-16,  0,  0,  0,  0,  0,  0],
             [-2,  0,  1,  1,  2,  0,  0,  1,  0, -2,  0,  0,  0,  0],
             [-1,  0,  1, -1,  1,  0, 18,-17,  0,  0,  0,  0,  0,  0],
             [-1,  0,  0,  1,  1,  0,  0,  2, -2,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0, -8, 13,  0,  0,  0,  0,  0,  2],
             [ 0,  0,  2, -2,  2,  0, -8, 11,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0, -8, 13,  0,  0,  0,  0,  0,  1],
             [ 0,  0,  1, -1,  1,  0, -8, 12,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  8,-13,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  1, -1,  1,  0,  8,-14,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  8,-13,  0,  0,  0,  0,  0,  1],
             [-2,  0,  0,  2,  1,  0,  0,  2,  0, -4,  5,  0,  0,  0],
             [-2,  0,  0,  2,  2,  0,  3, -3,  0,  0,  0,  0,  0,  0],
             [-2,  0,  0,  2,  0,  0,  0,  2,  0, -3,  1,  0,  0,  0],
             [ 0,  0,  0,  0,  1,  0,  3, -5,  0,  2,  0,  0,  0,  0],
             [-2,  0,  0,  2,  0,  0,  0,  2,  0, -4,  3,  0,  0,  0],
             [ 0,  0, -1,  1,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  1,  0,  0, -1,  2,  0,  0,  0,  0,  0],
             [ 0,  0,  1, -1,  2,  0,  0, -2,  2,  0,  0,  0,  0,  0],
             [-1,  0,  1,  0,  1,  0,  3, -5,  0,  0,  0,  0,  0,  0],
             [-1,  0,  0,  1,  0,  0,  3, -4,  0,  0,  0,  0,  0,  0],
             [-2,  0,  0,  2,  0,  0,  0,  2,  0, -2, -2,  0,  0,  0],
             [-2,  0,  2,  0,  2,  0,  0, -5,  9,  0,  0,  0,  0,  0],
             [ 0,  0,  1, -1,  1,  0,  0, -1,  0,  0,  0, -1,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0],
             [ 0,  0,  1, -1,  1,  0,  0, -1,  0,  0,  0,  0,  2,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  1],
             [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  2],
             [-1,  0,  0,  1,  0,  0,  0,  3, -4,  0,  0,  0,  0,  0],
             [ 0,  0, -1,  1,  0,  0,  0,  1,  0,  0,  2,  0,  0,  0],
             [ 0,  0,  1, -1,  2,  0,  0, -1,  0,  0,  2,  0,  0,  0],
             [ 0,  0,  0,  0,  1,  0,  0, -9, 17,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  2,  0, -3,  5,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  1, -1,  1,  0,  0, -1,  0, -1,  2,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  1, -2,  0,  0,  0],
             [ 1,  0,  0, -2,  0,  0, 17,-16,  0, -2,  0,  0,  0,  0],
             [ 0,  0,  1, -1,  1,  0,  0, -1,  0,  1, -3,  0,  0,  0],
             [-2,  0,  0,  2,  1,  0,  0,  5, -6,  0,  0,  0,  0,  0],
             [ 0,  0, -2,  2,  0,  0,  0,  9,-13,  0,  0,  0,  0,  0],
             [ 0,  0,  1, -1,  2,  0,  0, -1,  0,  0,  1,  0,  0,  0],
             [ 0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  1,  0,  0,  0],
             [ 0,  0, -1,  1,  0,  0,  0,  1,  0,  0,  1,  0,  0,  0],
             [ 0,  0, -2,  2,  0,  0,  5, -6,  0,  0,  0,  0,  0,  0],
             [ 0,  0, -1,  1,  1,  0,  5, -7,  0,  0,  0,  0,  0,  0],
             [-2,  0,  0,  2,  0,  0,  6, -8,  0,  0,  0,  0,  0,  0],
             [ 2,  0,  1, -3,  1,  0, -6,  7,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  2,  0,  0,  0,  0,  1,  0,  0,  0,  0],
             [ 0,  0, -1,  1,  1,  0,  0,  1,  0,  1,  0,  0,  0,  0],
             [ 0,  0,  1, -1,  1,  0,  0, -1,  0,  0,  0,  2,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  1],
             [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0, -8, 15,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0, -8, 15,  0,  0,  0,  0,  1],
             [ 0,  0,  1, -1,  1,  0,  0, -9, 15,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  8,-15,  0,  0,  0,  0,  0],
             [ 1,  0, -1, -1,  0,  0,  0,  8,-15,  0,  0,  0,  0,  0],
             [ 2,  0,  0, -2,  0,  0,  2, -5,  0,  0,  0,  0,  0,  0],
             [-2,  0,  0,  2,  0,  0,  0,  2,  0, -5,  5,  0,  0,  0],
             [ 2,  0,  0, -2,  1,  0,  0, -6,  8,  0,  0,  0,  0,  0],
             [ 2,  0,  0, -2,  1,  0,  0, -2,  0,  3,  0,  0,  0,  0],
             [-2,  0,  1,  1,  0,  0,  0,  1,  0, -3,  0,  0,  0,  0],
             [-2,  0,  1,  1,  1,  0,  0,  1,  0, -3,  0,  0,  0,  0],
             [-2,  0,  0,  2,  0,  0,  0,  2,  0, -3,  0,  0,  0,  0],
             [-2,  0,  0,  2,  0,  0,  0,  6, -8,  0,  0,  0,  0,  0],
             [-2,  0,  0,  2,  0,  0,  0,  2,  0, -1, -5,  0,  0,  0],
             [-1,  0,  0,  1,  0,  0,  0,  1,  0, -1,  0,  0,  0,  0],
             [-1,  0,  1,  1,  1,  0,-20, 20,  0,  0,  0,  0,  0,  0],
             [ 1,  0,  0, -2,  0,  0, 20,-21,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  1,  0,  0,  8,-15,  0,  0,  0,  0,  0],
             [ 0,  0,  2, -2,  1,  0,  0,-10, 15,  0,  0,  0,  0,  0],
             [ 0,  0, -1,  1,  0,  0,  0,  1,  0,  1,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  1,  0,  0,  0,  0,  1,  0,  0,  0,  0],
             [ 0,  0,  1, -1,  2,  0,  0, -1,  0,  1,  0,  0,  0,  0],
             [ 0,  0,  1, -1,  1,  0,  0, -1,  0, -2,  4,  0,  0,  0],
             [ 2,  0,  0, -2,  1,  0, -6,  8,  0,  0,  0,  0,  0,  0],
             [ 0,  0, -2,  2,  1,  0,  5, -6,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  1],
             [ 0,  0,  1, -1,  1,  0,  0, -1,  0,  0, -1,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0],
             [ 0,  0,  1, -1,  1,  0,  0, -1,  0,  0,  1,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  1],
             [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  2],
             [ 0,  0,  2, -2,  1,  0,  0, -9, 13,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  1,  0,  0,  7,-13,  0,  0,  0,  0,  0],
             [-2,  0,  0,  2,  0,  0,  0,  5, -6,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  9,-17,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0, -9, 17,  0,  0,  0,  0,  2],
             [ 1,  0,  0, -1,  1,  0,  0, -3,  4,  0,  0,  0,  0,  0],
             [ 1,  0,  0, -1,  1,  0, -3,  4,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  2,  0,  0, -1,  2,  0,  0,  0,  0,  0],
             [ 0,  0, -1,  1,  1,  0,  0,  0,  2,  0,  0,  0,  0,  0],
             [ 0,  0, -2,  2,  0,  1,  0, -2,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  3, -5,  0,  2,  0,  0,  0,  0],
             [-2,  0,  0,  2,  1,  0,  0,  2,  0, -3,  1,  0,  0,  0],
             [-2,  0,  0,  2,  1,  0,  3, -3,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  1,  0,  8,-13,  0,  0,  0,  0,  0,  0],
             [ 0,  0, -1,  1,  0,  0,  8,-12,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  2, -2,  1,  0, -8, 11,  0,  0,  0,  0,  0,  0],
             [-1,  0,  0,  1,  0,  0,  0,  2, -2,  0,  0,  0,  0,  0],
             [-1,  0,  0,  0,  1,  0, 18,-16,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  1, -1,  1,  0,  0, -1,  0, -1,  1,  0,  0,  0],
             [ 0,  0,  0,  0,  1,  0,  3, -7,  4,  0,  0,  0,  0,  0],
             [-2,  0,  1,  1,  1,  0,  0, -3,  7,  0,  0,  0,  0,  0],
             [ 0,  0,  1, -1,  2,  0,  0, -1,  0, -2,  5,  0,  0,  0],
             [ 0,  0,  0,  0,  1,  0,  0,  0,  0, -2,  5,  0,  0,  0],
             [ 0,  0,  0,  0,  1,  0,  0, -4,  8, -3,  0,  0,  0,  0],
             [ 1,  0,  0,  0,  1,  0,-10,  3,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  2, -2,  1,  0,  0, -2,  0,  0,  0,  0,  0,  0],
             [-1,  0,  0,  0,  1,  0, 10, -3,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  1,  0,  0,  4, -8,  3,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  1,  0,  0,  0,  0,  2, -5,  0,  0,  0],
             [ 0,  0, -1,  1,  0,  0,  0,  1,  0,  2, -5,  0,  0,  0],
             [ 2,  0, -1, -1,  1,  0,  0,  3, -7,  0,  0,  0,  0,  0],
             [-2,  0,  0,  2,  0,  0,  0,  2,  0,  0, -5,  0,  0,  0],
             [ 0,  0,  0,  0,  1,  0, -3,  7, -4,  0,  0,  0,  0,  0],
             [-2,  0,  0,  2,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0],
             [ 1,  0,  0,  0,  1,  0,-18, 16,  0,  0,  0,  0,  0,  0],
             [-2,  0,  1,  1,  1,  0,  0,  1,  0, -2,  0,  0,  0,  0],
             [ 0,  0,  1, -1,  2,  0, -8, 12,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  1,  0, -8, 13,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  1, -2,  0,  0,  0,  0,  1],
             [ 0,  0,  1, -1,  1,  0,  0,  0, -2,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  1, -2,  0,  0,  0,  0,  0],
             [ 0,  0,  1, -1,  1,  0,  0, -2,  2,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0, -1,  2,  0,  0,  0,  0,  1],
             [-1,  0,  0,  1,  1,  0,  3, -4,  0,  0,  0,  0,  0,  0],
             [-1,  0,  0,  1,  1,  0,  0,  3, -4,  0,  0,  0,  0,  0],
             [ 0,  0,  1, -1,  1,  0,  0, -1,  0,  0, -2,  0,  0,  0],
             [ 0,  0,  1, -1,  1,  0,  0, -1,  0,  0,  2,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  1],
             [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  2],
             [ 0,  0,  1, -1,  0,  0,  3, -6,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  1,  0, -3,  5,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  1, -1,  2,  0, -3,  4,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  1,  0,  0, -2,  4,  0,  0,  0,  0,  0],
             [ 0,  0,  2, -2,  1,  0, -5,  6,  0,  0,  0,  0,  0,  0],
             [ 0,  0, -1,  1,  0,  0,  5, -7,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  1,  0,  5, -8,  0,  0,  0,  0,  0,  0],
             [-2,  0,  0,  2,  1,  0,  6, -8,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  1,  0,  0, -8, 15,  0,  0,  0,  0,  0],
             [-2,  0,  0,  2,  1,  0,  0,  2,  0, -3,  0,  0,  0,  0],
             [-2,  0,  0,  2,  1,  0,  0,  6, -8,  0,  0,  0,  0,  0],
             [ 1,  0,  0, -1,  1,  0,  0, -1,  0,  1,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  3, -5,  0,  0,  0],
             [ 0,  0,  1, -1,  1,  0,  0, -1,  0, -1,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  0,  1],
             [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  1],
             [ 0,  0,  1, -1,  1,  0,  0, -1,  0,  1,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  1],
             [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  2],
             [ 0,  0,  1, -1,  2,  0,  0, -1,  0,  0, -1,  0,  0,  0],
             [ 0,  0,  0,  0,  1,  0,  0,  0,  0,  0, -1,  0,  0,  0],
             [ 0,  0, -1,  1,  0,  0,  0,  1,  0,  0, -1,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0, -7, 13,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  7,-13,  0,  0,  0,  0,  0],
             [ 2,  0,  0, -2,  1,  0,  0, -5,  6,  0,  0,  0,  0,  0],
             [ 0,  0,  2, -2,  1,  0,  0, -8, 11,  0,  0,  0,  0,  0],
             [ 0,  0,  2, -2,  1, -1,  0,  2,  0,  0,  0,  0,  0,  0],
             [-2,  0,  0,  2,  0,  0,  0,  4, -4,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  0],
             [ 0,  0,  1, -1,  1,  0,  0, -1,  0,  0,  3,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  1],
             [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  2],
             [-2,  0,  0,  2,  0,  0,  3, -3,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  2,  0,  0, -4,  8, -3,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  2,  0,  0,  4, -8,  3,  0,  0,  0,  0],
             [ 2,  0,  0, -2,  1,  0,  0, -2,  0,  2,  0,  0,  0,  0],
             [ 0,  0,  1, -1,  2,  0,  0, -1,  0,  2,  0,  0,  0,  0],
             [ 0,  0,  1, -1,  2,  0,  0,  0, -2,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  1,  0,  0,  1, -2,  0,  0,  0,  0,  0],
             [ 0,  0, -1,  1,  0,  0,  0,  2, -2,  0,  0,  0,  0,  0],
             [ 0,  0, -1,  1,  0,  0,  0,  1,  0,  0, -2,  0,  0,  0],
             [ 0,  0,  2, -2,  1,  0,  0, -2,  0,  0,  2,  0,  0,  0],
             [ 0,  0,  1, -1,  1,  0,  3, -6,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  3, -5,  0,  0,  0,  0,  0,  1],
             [ 0,  0,  0,  0,  0,  0,  3, -5,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  1, -1,  1,  0, -3,  4,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0, -3,  5,  0,  0,  0,  0,  0,  1],
             [ 0,  0,  0,  0,  0,  0, -3,  5,  0,  0,  0,  0,  0,  2],
             [ 0,  0,  2, -2,  2,  0, -3,  3,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0, -3,  5,  0,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  2, -4,  0,  0,  0,  0,  1],
             [ 0,  0,  1, -1,  1,  0,  0,  1, -4,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  2, -4,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0, -2,  4,  0,  0,  0,  0,  1],
             [ 0,  0,  1, -1,  1,  0,  0, -3,  4,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0, -2,  4,  0,  0,  0,  0,  1],
             [ 0,  0,  0,  0,  0,  0,  0, -2,  4,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0, -5,  8,  0,  0,  0,  0,  0,  2],
             [ 0,  0,  2, -2,  2,  0, -5,  6,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0, -5,  8,  0,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0, -5,  8,  0,  0,  0,  0,  0,  1],
             [ 0,  0,  1, -1,  1,  0, -5,  7,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0, -5,  8,  0,  0,  0,  0,  0,  1],
             [ 0,  0,  0,  0,  0,  0,  5, -8,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  1, -1,  2,  0,  0, -1,  0, -1,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  1,  0,  0,  0,  0, -1,  0,  0,  0,  0],
             [ 0,  0, -1,  1,  0,  0,  0,  1,  0, -1,  0,  0,  0,  0],
             [ 0,  0,  2, -2,  1,  0,  0, -2,  0,  1,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0, -6, 11,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  6,-11,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0, -1,  0,  4,  0,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  1,  0, -4,  0,  0,  0,  0,  0,  0],
             [ 2,  0,  0, -2,  1,  0, -3,  3,  0,  0,  0,  0,  0,  0],
             [-2,  0,  0,  2,  0,  0,  0,  2,  0,  0, -2,  0,  0,  0],
             [ 0,  0,  2, -2,  1,  0,  0, -7,  9,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  4, -5,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  1],
             [ 0,  0,  1, -1,  1,  0,  0, -1,  0,  2,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  1],
             [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  2],
             [ 0,  0,  2, -2,  2,  0,  0, -2,  0,  2,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  5,  0,  0,  2],
             [ 0,  0,  0,  0,  1,  0,  3, -5,  0,  0,  0,  0,  0,  0],
             [ 0,  0, -1,  1,  0,  0,  3, -4,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  2, -2,  1,  0, -3,  3,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  1,  0,  0,  2, -4,  0,  0,  0,  0,  0],
             [ 0,  0,  2, -2,  1,  0,  0, -4,  4,  0,  0,  0,  0,  0],
             [ 0,  0,  1, -1,  2,  0, -5,  7,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  3, -6,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0, -3,  6,  0,  0,  0,  0,  1],
             [ 0,  0,  1, -1,  1,  0,  0, -4,  6,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0, -3,  6,  0,  0,  0,  0,  1],
             [ 0,  0,  0,  0,  0,  0,  0, -3,  6,  0,  0,  0,  0,  2],
             [ 0,  0, -1,  1,  0,  0,  2, -2,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  1,  0,  2, -3,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0, -5,  9,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0, -5,  9,  0,  0,  0,  0,  1],
             [ 0,  0,  0,  0,  0,  0,  0,  5, -9,  0,  0,  0,  0,  0],
             [ 0,  0, -1,  1,  0,  0,  0,  1,  0, -2,  0,  0,  0,  0],
             [ 0,  0,  2, -2,  1,  0,  0, -2,  0,  2,  0,  0,  0,  0],
             [-2,  0,  1,  1,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0],
             [ 0,  0, -2,  2,  0,  0,  3, -3,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0, -6, 10,  0,  0,  0,  0,  0,  1],
             [ 0,  0,  0,  0,  0,  0, -6, 10,  0,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0, -2,  3,  0,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0, -2,  3,  0,  0,  0,  0,  0,  1],
             [ 0,  0,  1, -1,  1,  0, -2,  2,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  2, -3,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  2, -3,  0,  0,  0,  0,  0,  1],
             [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  0,  1],
             [ 0,  0,  1, -1,  1,  0,  0, -1,  0,  3,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  0,  1],
             [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  4, -8,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0, -4,  8,  0,  0,  0,  0,  2],
             [ 0,  0, -2,  2,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0, -4,  7,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0, -4,  7,  0,  0,  0,  0,  1],
             [ 0,  0,  0,  0,  0,  0,  0,  4, -7,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  1,  0, -2,  3,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  2, -2,  1,  0,  0, -2,  0,  3,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0, -5, 10,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  1,  0, -1,  2,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  4,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0, -3,  5,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0, -3,  5,  0,  0,  0,  0,  1],
             [ 0,  0,  0,  0,  0,  0,  0,  3, -5,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  1, -2,  0,  0,  0,  0,  0,  1],
             [ 0,  0,  1, -1,  1,  0,  1, -3,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  1, -2,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0, -1,  2,  0,  0,  0,  0,  0,  1],
             [ 0,  0,  0,  0,  0,  0, -1,  2,  0,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0, -7, 11,  0,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0, -7, 11,  0,  0,  0,  0,  0,  1],
             [ 0,  0, -2,  2,  0,  0,  4, -4,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  2, -3,  0,  0,  0,  0,  0],
             [ 0,  0,  2, -2,  1,  0, -4,  4,  0,  0,  0,  0,  0,  0],
             [ 0,  0, -1,  1,  0,  0,  4, -5,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0, -4,  7,  0,  0,  0,  0,  0,  1],
             [ 0,  0,  1, -1,  1,  0, -4,  6,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0, -4,  7,  0,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0, -4,  6,  0,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0, -4,  6,  0,  0,  0,  0,  0,  1],
             [ 0,  0,  1, -1,  1,  0, -4,  5,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0, -4,  6,  0,  0,  0,  0,  0,  1],
             [ 0,  0,  0,  0,  0,  0,  4, -6,  0,  0,  0,  0,  0,  0],
             [-2,  0,  0,  2,  0,  0,  2, -2,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0],
             [ 0,  0, -1,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  1,  0,  1, -1,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0, -1,  0,  5,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  1, -3,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0, -1,  3,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0, -7, 12,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0, -1,  1,  0,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0, -1,  1,  0,  0,  0,  0,  0,  1],
             [ 0,  0,  1, -1,  1,  0, -1,  0,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0,  1],
             [ 0,  0,  1, -1,  1,  0,  1, -2,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0, -2,  5,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0, -1,  0,  4,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  1,  0, -4,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  1,  0, -1,  1,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0, -6, 10,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0, -6, 10,  0,  0,  0,  0,  0],
             [ 0,  0,  2, -2,  1,  0,  0, -3,  0,  3,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0, -3,  7,  0,  0,  0,  0,  2],
             [-2,  0,  0,  2,  0,  0,  4, -4,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0, -5,  8,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  5, -8,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0, -1,  0,  3,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0, -1,  0,  3,  0,  0,  0,  1],
             [ 0,  0,  0,  0,  0,  0,  0,  1,  0, -3,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  2, -4,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0, -2,  4,  0,  0,  0,  0,  0,  1],
             [ 0,  0,  1, -1,  1,  0, -2,  3,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0, -2,  4,  0,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0, -6,  9,  0,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0, -6,  9,  0,  0,  0,  0,  0,  1],
             [ 0,  0,  0,  0,  0,  0,  6, -9,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  1,  0,  0,  1,  0, -2,  0,  0,  0,  0],
             [ 0,  0,  2, -2,  1,  0, -2,  2,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0, -4,  6,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  4, -6,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  1,  0,  3, -4,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0, -1,  0,  2,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  1,  0, -2,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  1,  0,  0,  1,  0, -1,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0, -5,  9,  0,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  3, -4,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0, -3,  4,  0,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0, -3,  4,  0,  0,  0,  0,  0,  1],
             [ 0,  0,  0,  0,  0,  0,  3, -4,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  3, -4,  0,  0,  0,  0,  0,  1],
             [ 0,  0,  0,  0,  1,  0,  0,  2, -2,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  1,  0,  0, -1,  0,  2,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  1,  0,  0, -3,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  1,  0,  1, -5,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0, -1,  0,  1,  0,  0,  0,  1],
             [ 0,  0,  0,  0,  0,  0,  0,  1,  0, -1,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  1,  0, -1,  0,  0,  0,  1],
             [ 0,  0,  0,  0,  0,  0,  0,  1,  0, -3,  5,  0,  0,  0],
             [ 0,  0,  0,  0,  1,  0, -3,  4,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  1,  0,  0, -2,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  1,  0,  0, -1,  0,  0,  0],
             [ 0,  0,  0,  0,  1,  0,  0, -1,  0,  1,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  1,  0,  0, -2,  2,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0, -8, 14,  0,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  1,  0,  2, -5,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  5, -8,  3,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  5, -8,  3,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  0,  0,  0,  1],
             [ 0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  3, -8,  3,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0, -3,  8, -3,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  1,  0, -2,  5,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0, -8, 12,  0,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0, -8, 12,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  1,  0,  1, -2,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  1,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  2,  0,  0,  2],
             [ 0,  0,  2, -2,  1,  0, -5,  5,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  0,  0,  1],
             [ 0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  3, -6,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0, -3,  6,  0,  0,  0,  0,  0,  1],
             [ 0,  0,  0,  0,  0,  0, -3,  6,  0,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0, -1,  4,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0, -5,  7,  0,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0, -5,  7,  0,  0,  0,  0,  0,  1],
             [ 0,  0,  1, -1,  1,  0, -5,  6,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  5, -7,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  2, -2,  1,  0,  0, -1,  0,  1,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0, -1,  0,  1,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0, -1,  0,  3,  0,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  1,  0,  2,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0, -2,  6,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  1,  0,  2, -2,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0, -6,  9,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  6, -9,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0, -2,  2,  0,  0,  0,  0,  0,  1],
             [ 0,  0,  1, -1,  1,  0, -2,  1,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  0,  0,  0,  1],
             [ 0,  0,  0,  0,  0,  0,  0,  1,  0,  3,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0, -5,  7,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  5, -7,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  1,  0, -2,  2,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  4, -5,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  1, -3,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0, -1,  3,  0,  0,  0,  0,  0,  1],
             [ 0,  0,  1, -1,  1,  0, -1,  2,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0, -1,  3,  0,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0, -7, 10,  0,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0, -7, 10,  0,  0,  0,  0,  0,  1],
             [ 0,  0,  0,  0,  0,  0,  0,  3, -3,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0, -4,  8,  0,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0, -4,  5,  0,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0, -4,  5,  0,  0,  0,  0,  0,  1],
             [ 0,  0,  0,  0,  0,  0,  4, -5,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0, -2,  0,  5,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0, -9, 13,  0,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0, -1,  5,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0, -2,  0,  4,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  2,  0, -4,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0, -2,  7,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  2,  0, -3,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0, -2,  5,  0,  0,  0,  0,  0,  1],
             [ 0,  0,  0,  0,  0,  0, -2,  5,  0,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0, -6,  8,  0,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0, -6,  8,  0,  0,  0,  0,  0,  1],
             [ 0,  0,  0,  0,  0,  0,  6, -8,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  1,  0,  0,  2,  0, -2,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0, -3,  9,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  5, -6,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  5, -6,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  2,  0, -2,  0,  0,  0,  1],
             [ 0,  0,  0,  0,  0,  0,  0,  2,  0, -2,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0, -5, 10,  0,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  4, -4,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  4, -4,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0, -3,  3,  0,  0,  0,  0,  0,  1],
             [ 0,  0,  0,  0,  0,  0,  3, -3,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  3, -3,  0,  0,  0,  0,  0,  1],
             [ 0,  0,  0,  0,  0,  0,  3, -3,  0,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  2,  0,  0, -3,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0, -5, 13,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  2,  0, -1,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  2,  0, -1,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  2,  0,  0, -2,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  2,  0,  0, -2,  0,  0,  1],
             [ 0,  0,  0,  0,  0,  0,  0,  3, -2,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  3, -2,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  2,  0,  0, -1,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0, -6, 15,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0, -8, 15,  0,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0, -3,  9, -4,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  2,  0,  2, -5,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0, -2,  8, -1, -5,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  6, -8,  3,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  1],
             [ 0,  0,  1, -1,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  1],
             [ 0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0, -6, 16, -4, -5,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0, -2,  8, -3,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0, -2,  8, -3,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  6, -8,  1,  5,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  2,  0, -2,  5,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  3, -5,  4,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0, -8, 11,  0,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0, -8, 11,  0,  0,  0,  0,  0,  1],
             [ 0,  0,  0,  0,  0,  0, -8, 11,  0,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0, 11,  0,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  1,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  3, -3,  0,  2,  0,  0,  0,  2],
             [ 0,  0,  2, -2,  1,  0,  0,  4, -8,  3,  0,  0,  0,  0],
             [ 0,  0,  1, -1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  2, -2,  1,  0,  0, -4,  8, -3,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  1,  2,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  2,  0,  1,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0, -3,  7,  0,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  0,  4,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0, -5,  6,  0,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0, -5,  6,  0,  0,  0,  0,  0,  1],
             [ 0,  0,  0,  0,  0,  0,  5, -6,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  5, -6,  0,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  2,  0,  2,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0, -1,  6,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  7, -9,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  2, -1,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  2, -1,  0,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  6, -7,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  5, -5,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0, -1,  4,  0,  0,  0,  0,  0,  1],
             [ 0,  0,  0,  0,  0,  0, -1,  4,  0,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0, -7,  9,  0,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0, -7,  9,  0,  0,  0,  0,  0,  1],
             [ 0,  0,  0,  0,  0,  0,  0,  4, -3,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  3, -1,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0, -4,  4,  0,  0,  0,  0,  0,  1],
             [ 0,  0,  0,  0,  0,  0,  4, -4,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  4, -4,  0,  0,  0,  0,  0,  1],
             [ 0,  0,  0,  0,  0,  0,  4, -4,  0,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  2,  1,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0, -3,  0,  5,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0,  0,  1],
             [ 0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0, -9, 12,  0,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  3,  0, -4,  0,  0,  0,  0],
             [ 0,  0,  2, -2,  1,  0,  1, -1,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  7, -8,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  3,  0, -3,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  3,  0, -3,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0, -2,  6,  0,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0, -6,  7,  0,  0,  0,  0,  0,  1],
             [ 0,  0,  0,  0,  0,  0,  6, -7,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  6, -6,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  3,  0, -2,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  3,  0, -2,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  5, -4,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  3, -2,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  3, -2,  0,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  3,  0, -1,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  3,  0, -1,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  3,  0,  0, -2,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  4, -2,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  3,  0,  0, -1,  0,  0,  2],
             [ 0,  0,  2, -2,  1,  0,  0,  1,  0, -1,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0, -8, 16,  0,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  3,  0,  2, -5,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  7, -8,  3,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0, -5, 16, -4, -5,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0, -1,  8, -3,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0, -8, 10,  0,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0, -8, 10,  0,  0,  0,  0,  0,  1],
             [ 0,  0,  0,  0,  0,  0, -8, 10,  0,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  2,  2,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  3,  0,  1,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0, -3,  8,  0,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0, -5,  5,  0,  0,  0,  0,  0,  1],
             [ 0,  0,  0,  0,  0,  0,  5, -5,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  5, -5,  0,  0,  0,  0,  0,  1],
             [ 0,  0,  0,  0,  0,  0,  5, -5,  0,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  1],
             [ 0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  7, -7,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  7, -7,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  6, -5,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  7, -8,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  5, -3,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  4, -3,  0,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  1,  2,  0,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0, -9, 11,  0,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0, -9, 11,  0,  0,  0,  0,  0,  1],
             [ 0,  0,  0,  0,  0,  0,  0,  4,  0, -4,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  4,  0, -3,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0, -6,  6,  0,  0,  0,  0,  0,  1],
             [ 0,  0,  0,  0,  0,  0,  6, -6,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  6, -6,  0,  0,  0,  0,  0,  1],
             [ 0,  0,  0,  0,  0,  0,  0,  4,  0, -2,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  6, -4,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  3, -1,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  3, -1,  0,  0,  0,  0,  0,  1],
             [ 0,  0,  0,  0,  0,  0,  3, -1,  0,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  4,  0, -1,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  4,  0,  0, -2,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  5, -2,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  4,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  8, -9,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  5, -4,  0,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  2,  1,  0,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  2,  1,  0,  0,  0,  0,  0,  1],
             [ 0,  0,  0,  0,  0,  0,  2,  1,  0,  0,  0,  0,  0,  1],
             [ 0,  0,  0,  0,  0,  0, -7,  7,  0,  0,  0,  0,  0,  1],
             [ 0,  0,  0,  0,  0,  0,  7, -7,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  4, -2,  0,  0,  0,  0,  0,  1],
             [ 0,  0,  0,  0,  0,  0,  4, -2,  0,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  4, -2,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  4, -2,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  5,  0, -4,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  5,  0, -3,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  5,  0, -2,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  3,  0,  0,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0, -8,  8,  0,  0,  0,  0,  0,  1],
             [ 0,  0,  0,  0,  0,  0,  8, -8,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  5, -3,  0,  0,  0,  0,  0,  1],
             [ 0,  0,  0,  0,  0,  0,  5, -3,  0,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0, -9,  9,  0,  0,  0,  0,  0,  1],
             [ 0,  0,  0,  0,  0,  0, -9,  9,  0,  0,  0,  0,  0,  1],
             [ 0,  0,  0,  0,  0,  0, -9,  9,  0,  0,  0,  0,  0,  1],
             [ 0,  0,  0,  0,  0,  0,  9, -9,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  6, -4,  0,  0,  0,  0,  0,  1],
             [ 0,  0,  0,  0,  0,  0,  0,  6,  0,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  6,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  6,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  6,  0,  0,  0,  0,  0,  1],
             [ 0,  0,  0,  0,  0,  0,  0,  6,  0,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  6,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  6,  0,  0,  0,  0,  0,  1],
             [ 0,  0,  0,  0,  0,  0,  0,  6,  0,  0,  0,  0,  0,  2],
             [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2],
             [ 1,  0,  0, -2,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0],
             [ 1,  0,  0, -2,  0,  0,  2, -2,  0,  0,  0,  0,  0,  0],
             [ 1,  0,  0, -2,  0,  0,  0,  1,  0, -1,  0,  0,  0,  0],
             [ 1,  0,  0, -2,  0,  0,  1, -1,  0,  0,  0,  0,  0,  0],
             [-1,  0,  0,  0,  0,  0,  3, -3,  0,  0,  0,  0,  0,  0],
             [-1,  0,  0,  0,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0],
             [-1,  0,  0,  2,  0,  0,  0,  4, -8,  3,  0,  0,  0,  0],
             [ 1,  0,  0, -2,  0,  0,  0,  4, -8,  3,  0,  0,  0,  0],
             [-2,  0,  0,  2,  0,  0,  0,  4, -8,  3,  0,  0,  0,  0],
             [-1,  0,  0,  0,  0,  0,  0,  2,  0, -3,  0,  0,  0,  0],
             [-1,  0,  0,  0,  0,  0,  0,  1,  0, -1,  0,  0,  0,  0],
             [-1,  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0,  0],
             [-1,  0,  0,  2,  0,  0,  2, -2,  0,  0,  0,  0,  0,  0],
             [ 1,  0, -1,  1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0],
             [-1,  0,  0,  2,  0,  0,  0,  2,  0, -3,  0,  0,  0,  0],
             [-2,  0,  0,  0,  0,  0,  0,  2,  0, -3,  0,  0,  0,  0],
             [ 1,  0,  0,  0,  0,  0,  0,  4, -8,  3,  0,  0,  0,  0],
             [-1,  0,  1, -1,  1,  0,  0, -1,  0,  0,  0,  0,  0,  0],
             [ 1,  0,  1, -1,  1,  0,  0, -1,  0,  0,  0,  0,  0,  0],
             [-1,  0,  0,  0,  0,  0,  0,  4, -8,  3,  0,  0,  0,  0],
             [-1,  0,  0,  2,  1,  0,  0,  2,  0, -2,  0,  0,  0,  0],
             [ 0,  0,  0,  0,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0],
             [-1,  0,  0,  2,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0],
             [-1,  0,  0,  2,  0,  0,  3, -3,  0,  0,  0,  0,  0,  0],
             [ 1,  0,  0, -2,  1,  0,  0, -2,  0,  2,  0,  0,  0,  0],
             [ 1,  0,  2, -2,  2,  0, -3,  3,  0,  0,  0,  0,  0,  0],
             [ 1,  0,  2, -2,  2,  0,  0, -2,  0,  2,  0,  0,  0,  0],
             [ 1,  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0,  0],
             [ 1,  0,  0,  0,  0,  0,  0,  1,  0, -1,  0,  0,  0,  0],
             [ 0,  0,  0, -2,  0,  0,  2, -2,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0, -2,  0,  0,  0,  1,  0, -1,  0,  0,  0,  0],
             [ 0,  0,  2,  0,  2,  0, -2,  2,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  2,  0,  2,  0,  0, -1,  0,  1,  0,  0,  0,  0],
             [ 0,  0,  2,  0,  2,  0, -1,  1,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  2,  0,  2,  0, -2,  3,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  0,  2,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0],
             [ 0,  0,  1,  1,  2,  0,  0,  1,  0,  0,  0,  0,  0,  0],
             [ 1,  0,  2,  0,  2,  0,  0,  1,  0,  0,  0,  0,  0,  0],
             [-1,  0,  2,  0,  2,  0, 10, -3,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  1,  1,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0],
             [ 1,  0,  2,  0,  2,  0,  0,  1,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  2,  0,  2,  0,  0,  4, -8,  3,  0,  0,  0,  0],
             [ 0,  0,  2,  0,  2,  0,  0, -4,  8, -3,  0,  0,  0,  0],
             [-1,  0,  2,  0,  2,  0,  0, -4,  8, -3,  0,  0,  0,  0],
             [ 2,  0,  2, -2,  2,  0,  0, -2,  0,  3,  0,  0,  0,  0],
             [ 1,  0,  2,  0,  1,  0,  0, -2,  0,  3,  0,  0,  0,  0],
             [ 0,  0,  1,  1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0],
             [-1,  0,  2,  0,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0],
             [-2,  0,  2,  2,  2,  0,  0,  2,  0, -2,  0,  0,  0,  0],
             [ 0,  0,  2,  0,  2,  0,  2, -3,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  2,  0,  2,  0,  1, -1,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  2,  0,  2,  0,  0,  1,  0, -1,  0,  0,  0,  0],
             [ 0,  0,  2,  0,  2,  0,  2, -2,  0,  0,  0,  0,  0,  0],
             [-1,  0,  2,  2,  2,  0,  0, -1,  0,  1,  0,  0,  0,  0],
             [ 1,  0,  2,  0,  2,  0, -1,  1,  0,  0,  0,  0,  0,  0],
             [-1,  0,  2,  2,  2,  0,  0,  2,  0, -3,  0,  0,  0,  0],
             [ 2,  0,  2,  0,  2,  0,  0,  2,  0, -3,  0,  0,  0,  0],
             [ 1,  0,  2,  0,  2,  0,  0, -4,  8, -3,  0,  0,  0,  0],
             [ 1,  0,  2,  0,  2,  0,  0,  4, -8,  3,  0,  0,  0,  0],
             [ 1,  0,  1,  1,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  2,  0,  2,  0,  0,  1,  0,  0,  0,  0,  0,  0],
             [ 2,  0,  2,  0,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0],
             [-1,  0,  2,  2,  2,  0,  0,  2,  0, -2,  0,  0,  0,  0],
             [-1,  0,  2,  2,  2,  0,  3, -3,  0,  0,  0,  0,  0,  0],
             [ 1,  0,  2,  0,  2,  0,  1, -1,  0,  0,  0,  0,  0,  0],
             [ 0,  0,  2,  2,  2,  0,  0,  2,  0, -2,  0,  0,  0,  0] ])

    # *
    # *  Planetary nutation coefficients, unit 1e-7 arcsec
    # *  longitude (sin, cos), obliquity (sin, cos)
    # *
    ICPL = np.asarray([
          [  1440,          0,          0,          0],
          [    56,       -117,        -42,        -40],
          [   125,        -43,          0,        -54],
          [     0,          5,          0,          0],
          [     3,         -7,         -3,          0],
          [     3,          0,          0,         -2],
          [  -114,          0,          0,         61],
          [  -219,         89,          0,          0],
          [    -3,          0,          0,          0],
          [  -462,       1604,          0,          0],
          [    99,          0,          0,        -53],
          [    -3,          0,          0,          2],
          [     0,          6,          2,          0],
          [     3,          0,          0,          0],
          [   -12,          0,          0,          0],
          [    14,       -218,        117,          8],
          [    31,       -481,       -257,        -17],
          [  -491,        128,          0,          0],
          [ -3084,       5123,       2735,       1647],
          [ -1444,       2409,      -1286,       -771],
          [    11,        -24,        -11,         -9],
          [    26,         -9,          0,          0],
          [   103,        -60,          0,          0],
          [     0,        -13,         -7,          0],
          [   -26,        -29,        -16,         14],
          [     9,        -27,        -14,         -5],
          [    12,          0,          0,         -6],
          [    -7,          0,          0,          0],
          [     0,         24,          0,          0],
          [   284,          0,          0,       -151],
          [   226,        101,          0,          0],
          [     0,         -8,         -2,          0],
          [     0,         -6,         -3,          0],
          [     5,          0,          0,         -3],
          [   -41,        175,         76,         17],
          [     0,         15,          6,          0],
          [   425,        212,       -133,        269],
          [  1200,        598,        319,       -641],
          [   235,        334,          0,          0],
          [    11,        -12,         -7,         -6],
          [     5,         -6,          3,          3],
          [    -5,          0,          0,          3],
          [     6,          0,          0,         -3],
          [    15,          0,          0,          0],
          [    13,          0,          0,         -7],
          [    -6,         -9,          0,          0],
          [   266,        -78,          0,          0],
          [  -460,       -435,       -232,        246],
          [     0,         15,          7,          0],
          [    -3,          0,          0,          2],
          [     0,        131,          0,          0],
          [     4,          0,          0,          0],
          [     0,          3,          0,          0],
          [     0,          4,          2,          0],
          [     0,          3,          0,          0],
          [   -17,        -19,        -10,          9],
          [    -9,        -11,          6,         -5],
          [    -6,          0,          0,          3],
          [   -16,          8,          0,          0],
          [     0,          3,          0,          0],
          [    11,         24,         11,         -5],
          [    -3,         -4,         -2,          1],
          [     3,          0,          0,         -1],
          [     0,         -8,         -4,          0],
          [     0,          3,          0,          0],
          [     0,          5,          0,          0],
          [     0,          3,          2,          0],
          [    -6,          4,          2,          3],
          [    -3,         -5,          0,          0],
          [    -5,          0,          0,          2],
          [     4,         24,         13,         -2],
          [   -42,         20,          0,          0],
          [   -10,        233,          0,          0],
          [    -3,          0,          0,          1],
          [    78,        -18,          0,          0],
          [     0,          3,          1,          0],
          [     0,         -3,         -1,          0],
          [     0,         -4,         -2,          1],
          [     0,         -8,         -4,         -1],
          [     0,         -5,          3,          0],
          [    -7,          0,          0,          3],
          [   -14,          8,          3,          6],
          [     0,          8,         -4,          0],
          [     0,         19,         10,          0],
          [    45,        -22,          0,          0],
          [    -3,          0,          0,          0],
          [     0,         -3,          0,          0],
          [     0,          3,          0,          0],
          [     3,          5,          3,         -2],
          [    89,        -16,         -9,        -48],
          [     0,          3,          0,          0],
          [    -3,          7,          4,          2],
          [  -349,        -62,          0,          0],
          [   -15,         22,          0,          0],
          [    -3,          0,          0,          0],
          [   -53,          0,          0,          0],
          [     5,          0,          0,         -3],
          [     0,         -8,          0,          0],
          [    15,         -7,         -4,         -8],
          [    -3,          0,          0,          1],
          [   -21,        -78,          0,          0],
          [    20,        -70,        -37,        -11],
          [     0,          6,          3,          0],
          [     5,          3,          2,         -2],
          [   -17,         -4,         -2,          9],
          [     0,          6,          3,          0],
          [    32,         15,         -8,         17],
          [   174,         84,         45,        -93],
          [    11,         56,          0,          0],
          [   -66,        -12,         -6,         35],
          [    47,          8,          4,        -25],
          [     0,          8,          4,          0],
          [    10,        -22,        -12,         -5],
          [    -3,          0,          0,          2],
          [   -24,         12,          0,          0],
          [     5,         -6,          0,          0],
          [     3,          0,          0,         -2],
          [     4,          3,          1,         -2],
          [     0,         29,         15,          0],
          [    -5,         -4,         -2,          2],
          [     8,         -3,         -1,         -5],
          [     0,         -3,          0,          0],
          [    10,          0,          0,          0],
          [     3,          0,          0,         -2],
          [    -5,          0,          0,          3],
          [    46,         66,         35,        -25],
          [   -14,          7,          0,          0],
          [     0,          3,          2,          0],
          [    -5,          0,          0,          0],
          [   -68,        -34,        -18,         36],
          [     0,         14,          7,          0],
          [    10,         -6,         -3,         -5],
          [    -5,         -4,         -2,          3],
          [    -3,          5,          2,          1],
          [    76,         17,          9,        -41],
          [    84,        298,        159,        -45],
          [     3,          0,          0,         -1],
          [    -3,          0,          0,          2],
          [    -3,          0,          0,          1],
          [   -82,        292,        156,         44],
          [   -73,         17,          9,         39],
          [    -9,        -16,          0,          0],
          [     3,          0,         -1,         -2],
          [    -3,          0,          0,          0],
          [    -9,         -5,         -3,          5],
          [  -439,          0,          0,          0],
          [    57,        -28,        -15,        -30],
          [     0,         -6,         -3,          0],
          [    -4,          0,          0,          2],
          [   -40,         57,         30,         21],
          [    23,          7,          3,        -13],
          [   273,         80,         43,       -146],
          [  -449,        430,          0,          0],
          [    -8,        -47,        -25,          4],
          [     6,         47,         25,         -3],
          [     0,         23,         13,          0],
          [    -3,          0,          0,          2],
          [     3,         -4,         -2,         -2],
          [   -48,       -110,        -59,         26],
          [    51,        114,         61,        -27],
          [  -133,          0,          0,         57],
          [     0,          4,          0,          0],
          [   -21,         -6,         -3,         11],
          [     0,         -3,         -1,          0],
          [   -11,        -21,        -11,          6],
          [   -18,       -436,       -233,          9],
          [    35,         -7,          0,          0],
          [     0,          5,          3,          0],
          [    11,         -3,         -1,         -6],
          [    -5,         -3,         -1,          3],
          [   -53,         -9,         -5,         28],
          [     0,          3,          2,          1],
          [     4,          0,          0,         -2],
          [     0,         -4,          0,          0],
          [   -50,        194,        103,         27],
          [   -13,         52,         28,          7],
          [   -91,        248,          0,          0],
          [     6,         49,         26,         -3],
          [    -6,        -47,        -25,          3],
          [     0,          5,          3,          0],
          [    52,         23,         10,        -23],
          [    -3,          0,          0,          1],
          [     0,          5,          3,          0],
          [    -4,          0,          0,          0],
          [    -4,          8,          3,          2],
          [    10,          0,          0,          0],
          [     3,          0,          0,         -2],
          [     0,          8,          4,          0],
          [     0,          8,          4,          1],
          [    -4,          0,          0,          0],
          [    -4,          0,          0,          0],
          [    -8,          4,          2,          4],
          [     8,         -4,         -2,         -4],
          [     0,         15,          7,          0],
          [  -138,          0,          0,          0],
          [     0,         -7,         -3,          0],
          [     0,         -7,         -3,          0],
          [    54,          0,          0,        -29],
          [     0,         10,          4,          0],
          [    -7,          0,          0,          3],
          [   -37,         35,         19,         20],
          [     0,          4,          0,          0],
          [    -4,          9,          0,          0],
          [     8,          0,          0,         -4],
          [    -9,        -14,         -8,          5],
          [    -3,         -9,         -5,          3],
          [  -145,         47,          0,          0],
          [   -10,         40,         21,          5],
          [    11,        -49,        -26,         -7],
          [ -2150,          0,          0,        932],
          [   -12,          0,          0,          5],
          [    85,          0,          0,        -37],
          [     4,          0,          0,         -2],
          [     3,          0,          0,         -2],
          [   -86,        153,          0,          0],
          [    -6,          9,          5,          3],
          [     9,        -13,         -7,         -5],
          [    -8,         12,          6,          4],
          [   -51,          0,          0,         22],
          [   -11,       -268,       -116,          5],
          [     0,         12,          5,          0],
          [     0,          7,          3,          0],
          [    31,          6,          3,        -17],
          [   140,         27,         14,        -75],
          [    57,         11,          6,        -30],
          [   -14,        -39,          0,          0],
          [     0,         -6,         -2,          0],
          [     4,         15,          8,         -2],
          [     0,          4,          0,          0],
          [    -3,          0,          0,          1],
          [     0,         11,          5,          0],
          [     9,          6,          0,          0],
          [    -4,         10,          4,          2],
          [     5,          3,          0,          0],
          [    16,          0,          0,         -9],
          [    -3,          0,          0,          0],
          [     0,          3,          2,         -1],
          [     7,          0,          0,         -3],
          [   -25,         22,          0,          0],
          [    42,        223,        119,        -22],
          [   -27,       -143,        -77,         14],
          [     9,         49,         26,         -5],
          [ -1166,          0,          0,        505],
          [    -5,          0,          0,          2],
          [    -6,          0,          0,          3],
          [    -8,          0,          1,          4],
          [     0,         -4,          0,          0],
          [   117,          0,          0,        -63],
          [    -4,          8,          4,          2],
          [     3,          0,          0,         -2],
          [    -5,          0,          0,          2],
          [     0,         31,          0,          0],
          [    -5,          0,          1,          3],
          [     4,          0,          0,         -2],
          [    -4,          0,          0,          2],
          [   -24,        -13,         -6,         10],
          [     3,          0,          0,          0],
          [     0,        -32,        -17,          0],
          [     8,         12,          5,         -3],
          [     3,          0,          0,         -1],
          [     7,         13,          0,          0],
          [    -3,         16,          0,          0],
          [    50,          0,          0,        -27],
          [     0,         -5,         -3,          0],
          [    13,          0,          0,          0],
          [     0,          5,          3,          1],
          [    24,          5,          2,        -11],
          [     5,        -11,         -5,         -2],
          [    30,         -3,         -2,        -16],
          [    18,          0,          0,         -9],
          [     8,        614,          0,          0],
          [     3,         -3,         -1,         -2],
          [     6,         17,          9,         -3],
          [    -3,         -9,         -5,          2],
          [     0,          6,          3,         -1],
          [  -127,         21,          9,         55],
          [     3,          5,          0,          0],
          [    -6,        -10,         -4,          3],
          [     5,          0,          0,          0],
          [    16,          9,          4,         -7],
          [     3,          0,          0,         -2],
          [     0,         22,          0,          0],
          [     0,         19,         10,          0],
          [     7,          0,          0,         -4],
          [     0,         -5,         -2,          0],
          [     0,          3,          1,          0],
          [    -9,          3,          1,          4],
          [    17,          0,          0,         -7],
          [     0,         -3,         -2,         -1],
          [   -20,         34,          0,          0],
          [   -10,          0,          1,          5],
          [    -4,          0,          0,          2],
          [    22,        -87,          0,          0],
          [    -4,          0,          0,          2],
          [    -3,         -6,         -2,          1],
          [   -16,         -3,         -1,          7],
          [     0,         -3,         -2,          0],
          [     4,          0,          0,          0],
          [   -68,         39,          0,          0],
          [    27,          0,          0,        -14],
          [     0,         -4,          0,          0],
          [   -25,          0,          0,          0],
          [   -12,         -3,         -2,          6],
          [     3,          0,          0,         -1],
          [     3,         66,         29,         -1],
          [   490,          0,          0,       -213],
          [   -22,         93,         49,         12],
          [    -7,         28,         15,          4],
          [    -3,         13,          7,          2],
          [   -46,         14,          0,          0],
          [    -5,          0,          0,          0],
          [     2,          1,          0,          0],
          [     0,         -3,          0,          0],
          [   -28,          0,          0,         15],
          [     5,          0,          0,         -2],
          [     0,          3,          0,          0],
          [   -11,          0,          0,          5],
          [     0,          3,          1,          0],
          [    -3,          0,          0,          1],
          [    25,        106,         57,        -13],
          [     5,         21,         11,         -3],
          [  1485,          0,          0,          0],
          [    -7,        -32,        -17,          4],
          [     0,          5,          3,          0],
          [    -6,         -3,         -2,          3],
          [    30,         -6,         -2,        -13],
          [    -4,          4,          0,          0],
          [   -19,          0,          0,         10],
          [     0,          4,          2,         -1],
          [     0,          3,          0,          0],
          [     4,          0,          0,         -2],
          [     0,         -3,         -1,          0],
          [    -3,          0,          0,          0],
          [     5,          3,          1,         -2],
          [     0,         11,          0,          0],
          [   118,          0,          0,        -52],
          [     0,         -5,         -3,          0],
          [   -28,         36,          0,          0],
          [     5,         -5,          0,          0],
          [    14,        -59,        -31,         -8],
          [     0,          9,          5,          1],
          [  -458,          0,          0,        198],
          [     0,        -45,        -20,          0],
          [     9,          0,          0,         -5],
          [     0,         -3,          0,          0],
          [     0,         -4,         -2,         -1],
          [    11,          0,          0,         -6],
          [     6,          0,          0,         -2],
          [   -16,         23,          0,          0],
          [     0,         -4,         -2,          0],
          [    -5,          0,          0,          2],
          [  -166,        269,          0,          0],
          [    15,          0,          0,         -8],
          [    10,          0,          0,         -4],
          [   -78,         45,          0,          0],
          [     0,         -5,         -2,          0],
          [     7,          0,          0,         -4],
          [    -5,        328,          0,          0],
          [     3,          0,          0,         -2],
          [     5,          0,          0,         -2],
          [     0,          3,          1,          0],
          [    -3,          0,          0,          0],
          [    -3,          0,          0,          0],
          [     0,         -4,         -2,          0],
          [ -1223,        -26,          0,          0],
          [     0,          7,          3,          0],
          [     3,          0,          0,          0],
          [     0,          3,          2,          0],
          [    -6,         20,          0,          0],
          [  -368,          0,          0,          0],
          [   -75,          0,          0,          0],
          [    11,          0,          0,         -6],
          [     3,          0,          0,         -2],
          [    -3,          0,          0,          1],
          [   -13,        -30,          0,          0],
          [    21,          3,          0,          0],
          [    -3,          0,          0,          1],
          [    -4,          0,          0,          2],
          [     8,        -27,          0,          0],
          [   -19,        -11,          0,          0],
          [    -4,          0,          0,          2],
          [     0,          5,          2,          0],
          [    -6,          0,          0,          2],
          [    -8,          0,          0,          0],
          [    -1,          0,          0,          0],
          [   -14,          0,          0,          6],
          [     6,          0,          0,          0],
          [   -74,          0,          0,         32],
          [     0,         -3,         -1,          0],
          [     4,          0,          0,         -2],
          [     8,         11,          0,          0],
          [     0,          3,          2,          0],
          [  -262,          0,          0,        114],
          [     0,         -4,          0,          0],
          [    -7,          0,          0,          4],
          [     0,        -27,        -12,          0],
          [   -19,         -8,         -4,          8],
          [   202,          0,          0,        -87],
          [    -8,         35,         19,          5],
          [     0,          4,          2,          0],
          [    16,         -5,          0,          0],
          [     5,          0,          0,         -3],
          [     0,         -3,          0,          0],
          [     1,          0,          0,          0],
          [   -35,        -48,        -21,         15],
          [    -3,         -5,         -2,          1],
          [     6,          0,          0,         -3],
          [     3,          0,          0,         -1],
          [     0,         -5,          0,          0],
          [    12,         55,         29,         -6],
          [     0,          5,          3,          0],
          [  -598,          0,          0,          0],
          [    -3,        -13,         -7,          1],
          [    -5,         -7,         -3,          2],
          [     3,          0,          0,         -1],
          [     5,         -7,          0,          0],
          [     4,          0,          0,         -2],
          [    16,         -6,          0,          0],
          [     8,         -3,          0,          0],
          [     8,        -31,        -16,         -4],
          [     0,          3,          1,          0],
          [   113,          0,          0,        -49],
          [     0,        -24,        -10,          0],
          [     4,          0,          0,         -2],
          [    27,          0,          0,          0],
          [    -3,          0,          0,          1],
          [     0,         -4,         -2,          0],
          [     5,          0,          0,         -2],
          [     0,         -3,          0,          0],
          [   -13,          0,          0,          6],
          [     5,          0,          0,         -2],
          [   -18,        -10,         -4,          8],
          [    -4,        -28,          0,          0],
          [    -5,          6,          3,          2],
          [    -3,          0,          0,          1],
          [    -5,         -9,         -4,          2],
          [    17,          0,          0,         -7],
          [    11,          4,          0,          0],
          [     0,         -6,         -2,          0],
          [    83,         15,          0,          0],
          [    -4,          0,          0,          2],
          [     0,       -114,        -49,          0],
          [   117,          0,          0,        -51],
          [    -5,         19,         10,          2],
          [    -3,          0,          0,          0],
          [    -3,          0,          0,          2],
          [     0,         -3,         -1,          0],
          [     3,          0,          0,          0],
          [     0,         -6,         -2,          0],
          [   393,          3,          0,          0],
          [    -4,         21,         11,          2],
          [    -6,          0,         -1,          3],
          [    -3,          8,          4,          1],
          [     8,          0,          0,          0],
          [    18,        -29,        -13,         -8],
          [     8,         34,         18,         -4],
          [    89,          0,          0,          0],
          [     3,         12,          6,         -1],
          [    54,        -15,         -7,        -24],
          [     0,          3,          0,          0],
          [     3,          0,          0,         -1],
          [     0,         35,          0,          0],
          [  -154,        -30,        -13,         67],
          [    15,          0,          0,          0],
          [     0,          4,          2,          0],
          [     0,          9,          0,          0],
          [    80,        -71,        -31,        -35],
          [     0,        -20,         -9,          0],
          [    11,          5,          2,         -5],
          [    61,        -96,        -42,        -27],
          [    14,          9,          4,         -6],
          [   -11,         -6,         -3,          5],
          [     0,         -3,         -1,          0],
          [   123,       -415,       -180,        -53],
          [     0,          0,          0,        -35],
          [    -5,          0,          0,          0],
          [     7,        -32,        -17,         -4],
          [     0,         -9,         -5,          0],
          [     0,         -4,          2,          0],
          [   -89,          0,          0,         38],
          [     0,        -86,        -19,         -6],
          [     0,          0,        -19,          6],
          [  -123,       -416,       -180,         53],
          [     0,         -3,         -1,          0],
          [    12,         -6,         -3,         -5],
          [   -13,          9,          4,          6],
          [     0,        -15,         -7,          0],
          [     3,          0,          0,         -1],
          [   -62,        -97,        -42,         27],
          [   -11,          5,          2,          5],
          [     0,        -19,         -8,          0],
          [    -3,          0,          0,          1],
          [     0,          4,          2,          0],
          [     0,          3,          0,          0],
          [     0,          4,          2,          0],
          [   -85,        -70,        -31,         37],
          [   163,        -12,         -5,        -72],
          [   -63,        -16,         -7,         28],
          [   -21,        -32,        -14,          9],
          [     0,         -3,         -1,          0],
          [     3,          0,          0,         -2],
          [     0,          8,          0,          0],
          [     3,         10,          4,         -1],
          [     3,          0,          0,         -1],
          [     0,         -7,         -3,          0],
          [     0,         -4,         -2,          0],
          [     6,         19,          0,          0],
          [     5,       -173,        -75,         -2],
          [     0,         -7,         -3,          0],
          [     7,        -12,         -5,         -3],
          [    -3,          0,          0,          2],
          [     3,         -4,         -2,         -1],
          [    74,          0,          0,        -32],
          [    -3,         12,          6,          2],
          [    26,        -14,         -6,        -11],
          [    19,          0,          0,         -8],
          [     6,         24,         13,         -3],
          [    83,          0,          0,          0],
          [     0,        -10,         -5,          0],
          [    11,         -3,         -1,         -5],
          [     3,          0,          1,         -1],
          [     3,          0,          0,         -1],
          [    -4,          0,          0,          0],
          [     5,        -23,        -12,         -3],
          [  -339,          0,          0,        147],
          [     0,        -10,         -5,          0],
          [     5,          0,          0,          0],
          [     3,          0,          0,         -1],
          [     0,         -4,         -2,          0],
          [    18,         -3,          0,          0],
          [     9,        -11,         -5,         -4],
          [    -8,          0,          0,          4],
          [     3,          0,          0,         -1],
          [     0,          9,          0,          0],
          [     6,         -9,         -4,         -2],
          [    -4,        -12,          0,          0],
          [    67,        -91,        -39,        -29],
          [    30,        -18,         -8,        -13],
          [     0,          0,          0,          0],
          [     0,       -114,        -50,          0],
          [     0,          0,          0,         23],
          [   517,         16,          7,       -224],
          [     0,         -7,         -3,          0],
          [   143,         -3,         -1,        -62],
          [    29,          0,          0,        -13],
          [    -4,          0,          0,          2],
          [    -6,          0,          0,          3],
          [     5,         12,          5,         -2],
          [   -25,          0,          0,         11],
          [    -3,          0,          0,          1],
          [     0,          4,          2,          0],
          [   -22,         12,          5,         10],
          [    50,          0,          0,        -22],
          [     0,          7,          4,          0],
          [     0,          3,          1,          0],
          [    -4,          4,          2,          2],
          [    -5,        -11,         -5,          2],
          [     0,          4,          2,          0],
          [     4,         17,          9,         -2],
          [    59,          0,          0,          0],
          [     0,         -4,         -2,          0],
          [    -8,          0,          0,          4],
          [    -3,          0,          0,          0],
          [     4,        -15,         -8,         -2],
          [   370,         -8,          0,       -160],
          [     0,          0,         -3,          0],
          [     0,          3,          1,          0],
          [    -6,          3,          1,          3],
          [     0,          6,          0,          0],
          [   -10,          0,          0,          4],
          [     0,          9,          4,          0],
          [     4,         17,          7,         -2],
          [    34,          0,          0,        -15],
          [     0,          5,          3,          0],
          [    -5,          0,          0,          2],
          [   -37,         -7,         -3,         16],
          [     3,         13,          7,         -2],
          [    40,          0,          0,          0],
          [     0,         -3,         -2,          0],
          [  -184,         -3,         -1,         80],
          [    -3,          0,          0,          1],
          [    -3,          0,          0,          0],
          [     0,        -10,         -6,         -1],
          [    31,         -6,          0,        -13],
          [    -3,        -32,        -14,          1],
          [    -7,          0,          0,          3],
          [     0,         -8,         -4,          0],
          [     3,         -4,          0,          0],
          [     0,          4,          0,          0],
          [     0,          3,          1,          0],
          [    19,        -23,        -10,          2],
          [     0,          0,          0,        -10],
          [     0,          3,          2,          0],
          [     0,          9,          5,         -1],
          [    28,          0,          0,          0],
          [     0,         -7,         -4,          0],
          [     8,         -4,          0,         -4],
          [     0,          0,         -2,          0],
          [     0,          3,          0,          0],
          [    -3,          0,          0,          1],
          [    -9,          0,          1,          4],
          [     3,         12,          5,         -1],
          [    17,         -3,         -1,          0],
          [     0,          7,          4,          0],
          [    19,          0,          0,          0],
          [     0,         -5,         -3,          0],
          [    14,         -3,          0,         -1],
          [     0,          0,         -1,          0],
          [     0,          0,          0,         -5],
          [     0,          5,          3,          0],
          [    13,          0,          0,          0],
          [     0,         -3,         -2,          0],
          [     2,          9,          4,          3],
          [     0,          0,          0,         -4],
          [     8,          0,          0,          0],
          [     0,          4,          2,          0],
          [     6,          0,          0,         -3],
          [     6,          0,          0,          0],
          [     0,          3,          1,          0],
          [     5,          0,          0,         -2],
          [     3,          0,          0,         -1],
          [    -3,          0,          0,          0],
          [     6,          0,          0,          0],
          [     7,          0,          0,          0],
          [    -4,          0,          0,          0],
          [     4,          0,          0,          0],
          [     6,          0,          0,          0],
          [     0,         -4,          0,          0],
          [     0,         -4,          0,          0],
          [     5,          0,          0,          0],
          [    -3,          0,          0,          0],
          [     4,          0,          0,          0],
          [    -5,          0,          0,          0],
          [     4,          0,          0,          0],
          [     0,          3,          0,          0],
          [    13,          0,          0,          0],
          [    21,         11,          0,          0],
          [     0,         -5,          0,          0],
          [     0,         -5,         -2,          0],
          [     0,          5,          3,          0],
          [     0,         -5,          0,          0],
          [    -3,          0,          0,          2],
          [    20,         10,          0,          0],
          [   -34,          0,          0,          0],
          [   -19,          0,          0,          0],
          [     3,          0,          0,         -2],
          [    -3,          0,          0,          1],
          [    -6,          0,          0,          3],
          [    -4,          0,          0,          0],
          [     3,          0,          0,          0],
          [     3,          0,          0,          0],
          [     4,          0,          0,          0],
          [     3,          0,          0,         -1],
          [     6,          0,          0,         -3],
          [    -8,          0,          0,          3],
          [     0,          3,          1,          0],
          [    -3,          0,          0,          0],
          [     0,         -3,         -2,          0],
          [   126,        -63,        -27,        -55],
          [    -5,          0,          1,          2],
          [    -3,         28,         15,          2],
          [     5,          0,          1,         -2],
          [     0,          9,          4,          1],
          [     0,          9,          4,         -1],
          [  -126,        -63,        -27,         55],
          [     3,          0,          0,         -1],
          [    21,        -11,         -6,        -11],
          [     0,         -4,          0,          0],
          [   -21,        -11,         -6,         11],
          [    -3,          0,          0,          1],
          [     0,          3,          1,          0],
          [     8,          0,          0,         -4],
          [    -6,          0,          0,          3],
          [    -3,          0,          0,          1],
          [     3,          0,          0,         -1],
          [    -3,          0,          0,          1],
          [    -5,          0,          0,          2],
          [    24,        -12,         -5,        -11],
          [     0,          3,          1,          0],
          [     0,          3,          1,          0],
          [     0,          3,          2,          0],
          [   -24,        -12,         -5,         10],
          [     4,          0,         -1,         -2],
          [    13,          0,          0,         -6],
          [     7,          0,          0,         -3],
          [     3,          0,          0,         -1],
          [     3,          0,          0,         -1] ])

    # * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # *  Interval between fundamental date J2000.0 and given date (JC).
    T = ( ( DATE1-DJ00 ) + DATE2) / DJC

    # *  -------------------
    # *  LUNI-SOLAR NUTATION
    # *  -------------------
    
    # *
    # *  Fundamental (Delaunay) arguments
    # *
    
    # *  Mean anomaly of the Moon (IERS 2003).
    EL = iau_FAL03 ( T )

    # *  Mean anomaly of the Sun (MHB2000).
    ELP = fmod(      1287104.79305 + 
                  T*( 129596581.0481 + 
                  T*(       - 0.5532 + 
                  T*(         0.000136 + 
                  T*(       - 0.00001149 )))), TURNAS ) * DAS2R

    # *  Mean longitude of the Moon minus that of the ascending node
    # *  (IERS 2003).
    F = iau_FAF03 ( T )

    # *  Mean elongation of the Moon from the Sun (MHB2000).
    D =    fmod(    1072260.70369 +
                T*( 1602961601.2090 +
                T*(        - 6.3706 +
                T*(          0.006593 +
                T*(        - 0.00003169 )))), TURNAS ) * DAS2R

    # *  Mean longitude of the ascending node of the Moon (IERS 2003).
    OM = iau_FAOM03 ( T )

    # *  Initialize the nutation values.
    DP = 0.0
    DE = 0.0

    # *  Summation of luni-solar nutation series (in reverse order).
    for I in range( NLS-1,-1,-1):

        # *     Argument and functions.
        ARG = fmod( float( NALS[I,0] ) * EL  + 
                       float( NALS[I,1] ) * ELP + 
                       float( NALS[I,2] ) * F   + 
                       float( NALS[I,3] ) * D   + 
                       float( NALS[I,4] ) * OM, D2PI )
        SARG = sin(ARG)
        CARG = cos(ARG)

        # *     Term.
        DP = DP + ( CLS[I,0] + CLS[I,1] * T ) * SARG +  CLS[I,2] * CARG
        DE = DE + ( CLS[I,3] + CLS[I,4] * T ) * CARG +  CLS[I,5] * SARG


    # *  Convert from 0.1 microarcsec units to radians.
    DPSILS = DP * U2R
    DEPSLS = DE * U2R

    # *  ------------------
    # *  PLANETARY NUTATION
    # *  ------------------
    
    # *  n.b.  The MHB2000 code computes the luni-solar and planetary nutation
    # *        in different routines, using slightly different Delaunay
    # *        arguments in the two cases.  This behaviour is faithfully
    # *        reproduced here.  Use of the IERS 2003 expressions for both
    # *        cases leads to negligible changes, well below
    # *        0.1 microarcsecond.
    
    # *  Mean anomaly of the Moon (MHB2000).
    AL = fmod( 2.35555598 + 8328.6914269554 * T, D2PI )

    # *  Mean anomaly of the Sun (MHB2000).
    ALSU = fmod( 6.24006013 + 628.301955 * T, D2PI )

    # *  Mean longitude of the Moon minus that of the ascending node
    # * (MHB2000).
    AF = fmod( 1.627905234 + 8433.466158131 * T, D2PI )

    # *  Mean elongation of the Moon from the Sun (MHB2000).
    AD = fmod( 5.198466741 + 7771.3771468121 * T, D2PI )

    # *  Mean longitude of the ascending node of the Moon (MHB2000).
    AOM = fmod( 2.18243920 - 33.757045 * T, D2PI )

    # *  General accumulated precession in longitude (IERS 2003).
    APA = iau_FAPA03 ( T )

    # *  Planetary longitudes, Mercury through Uranus (IERS 2003).
    ALME = iau_FAME03 ( T )
    ALVE = iau_FAVE03 ( T )
    ALEA = iau_FAE03 ( T )
    ALMA = iau_FAMA03 ( T )
    ALJU = iau_FAJU03 ( T )
    ALSA = iau_FASA03 ( T )
    ALUR = iau_FAUR03 ( T )

    # *  Neptune longitude (MHB2000).
    ALNE = fmod( 5.321159000 + 3.8127774000 * T, D2PI )
    #ALNE = iau_FANE03 ( T )],

    # *  Initialize the nutation values.
    DP = 0.0
    DE = 0.0

    # *  Summation of planetary nutation series (in reverse order).
    for I in range( NPL-1,-1,-1):

        # *     Argument and functions.
        ARG =    fmod( float( NAPL[I, 0] ) * AL   +
                       float( NAPL[I, 1] ) * ALSU +
                       float( NAPL[I, 2] ) * AF   +
                       float( NAPL[I, 3] ) * AD   +
                       float( NAPL[I, 4] ) * AOM  +
                       float( NAPL[I, 5] ) * ALME +
                       float( NAPL[I, 6] ) * ALVE +
                       float( NAPL[I, 7] ) * ALEA +
                       float( NAPL[I, 8] ) * ALMA +
                       float( NAPL[I, 9] ) * ALJU +
                       float( NAPL[I,10] ) * ALSA +
                       float( NAPL[I,11] ) * ALUR +
                       float( NAPL[I,12] ) * ALNE +
                       float( NAPL[I,13] ) * APA, D2PI )
        SARG = sin(ARG)
        CARG = cos(ARG)
        
        # *     Term.
        DP = DP + float( ICPL[I,0]) * SARG + float( ICPL[I,1]) * CARG
        DE = DE + float( ICPL[I,2]) * SARG + float( ICPL[I,3]) * CARG

    # *  Convert from 0.1 microarcsec units to radians.
    DPSIPL = DP * U2R
    DEPSPL = DE * U2R

    ## *  -------
    ## *  RESULTS
    ## *  -------

    # *  Add luni-solar and planetary components.
    DPSI = DPSILS + DPSIPL
    DEPS = DEPSLS + DEPSPL

    return DPSI, DEPS
    # *  Finished.



def iau_PN00A ( DATE1, DATE2 ):
    '''
    # *+
    # *  - - - - - - - - - -
    # *   i a u _ P N 0 0 A
    # *  - - - - - - - - - -
    # *
    # *  Precession-nutation, IAU 2000A remel:  a multi-purpose routine,
    # *  supporting classical (equinox-based) use directly and CIO-based
    # *  use indirectly.
    # *
    # *  This routine is part of the International Astronomical Union's
    # *  SOFA (Standards of Fundamental Astronomy) software collection.
    # *
    # *  Status:  support routine.
    # *
    # *  Given:
    # *     DATE1,DATE2   d       TT as a 2-part Julian Date (Note 1)
    # *
    # *  Returned:
    # *     DPSI,DEPS     d       nutation (Note 2)
    # *     EPSA          d       mean obliquity (Note 3)
    # *     RB          d(3,3)    frame bias matrix (Note 4)
    # *     RP          d(3,3)    precession matrix (Note 5)
    # *     RBP         d(3,3)    bias-precession matrix (Note 6)
    # *     RN          d(3,3)    nutation matrix (Note 7)	
    # *     RBPN        d(3,3)    GCRS-to-true matrix (Notes 8,9)
    # *
    # *  Notes:
    # *
    # *  1) The TT date DATE1+DATE2 is a Julian Date, apportioned in any
    # *     convenient way between the two arguments.  For example,
    # *     JD(TT)=2450123.7 could be expressed in any of these ways,
    # *     among others:
    # *
    # *            DATE1          DATE2
    # *
    # *         2450123.7.0        0.0        (JD method)
    # *          2451545.0      -1421.3.0     (J2000 method)
    # *         2400000.5.0     50123.2.0     (MJD method)
    # *         2450123.5.0       0.2.0       (date & time method)
    # *
    # *     The JD method is the most natural and convenient to use in
    # *     cases where the loss of several decimal digits of resolution
    # *     is acceptable.  The J2000 method is best matched to the way
    # *     the argument is handled internally and will deliver the
    # *     optimum resolution.  The MJD method and the date & time methods
    # *     are both good compromises between resolution and convenience.
    # *
    # *  2) The nutation components (luni-solar + planetary, IAU 2000A) in
    # *     longitude and obliquity are in radians and with respect to the
    # *     equinox and ecliptic of date.  Free core nutation is omitted],  for
    # *     the utmost accuracy, use the iau_PN00 routine, where the nutation
    # *     components are caller-specified.  For faster but slightly less
    # *     accurate results, use the iau_PN00B routine.
    # *
    # *  3) The mean obliquity is consistent with the IAU 2000 precession.
    # *
    # *  4) The matrix RB transforms vectors from GCRS to J2000.0 mean equator
    # *     and equinox by applying frame bias.
    # *
    # *  5) The matrix RP transforms vectors from J2000.0 mean equator and
    # *     equinox to mean equator and equinox of date by applying
    # *     precession.
    # *
    # *  6) The matrix RBP transforms vectors from GCRS to mean equator and
    # *     equinox of date by applying frame bias then precession.  It is the
    # *     product RP x RB.
    # *
    # *  7) The matrix RN transforms vectors from mean equator and equinox of
    # *     date to true equator and equinox of date by applying the nutation
    # *     (luni-solar + planetary).
    # *
    # *  8) The matrix RBPN transforms vectors from GCRS to true equator and
    # *     equinox of date.  It is the product RN x RBP, applying frame bias,
    # *     precession and nutation in that order.
    # *
    # *  9) The X,Y,Z coordinates of the IAU 2000A Celestial Intermediate Pole
    # *     are elements (3,1-3) of the matrix RBPN.
    # *
    # *  Called:
    # *     iau_NUT00A   nutation, IAU 2000A
    # *     iau_PN00     bias/precession/nutation results, IAU 2000
    # *
    # *  Reference:
    # *
    # *     Capitaine, N., Chapront, J., Lambert, S. and Wallace, P.,
    # *     "Expressions for the Celestial Intermediate Pole and Celestial
    # *     Ephemeris Origin consistent with the IAU 2000A precession-nutation
    # *     remel", Astron.Astrophys. 400, 1145-1154 (2003).
    # *
    # *     n.b. The celestial ephemeris origin (CEO) was renamed "celestial
    # *          intermediate origin" (CIO) by IAU 2006 Resolution 2.
    # *
    # *  This revision:  2010 January 18
    # *
    # *  SOFA release 2010-12-01
    # *
    # *  Copyright (C) 2010 IAU SOFA Board.  See notes at end.
    # *
    # * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    '''
    # *  Nutation.
    DPSI, DEPS = iau_NUT00A ( DATE1, DATE2 )

    # *  Remaining results.
    EPSA, RB, RP, RBP, RN, RBPN = iau_PN00 ( DATE1, DATE2, DPSI, DEPS )

    return DPSI, DEPS, EPSA, RB, RP, RBP, RN, RBPN
    # *  Finished.



def iau_PN00 ( DATE1, DATE2, DPSI, DEPS ):
    '''
    # *+
    # *  - - - - - - - - -
    # *   i a u _ P N 0 0
    # *  - - - - - - - - -
    # *
    # *  Precession-nutation, IAU 2000 remel:  a multi-purpose routine,
    # *  supporting classical (equinox-based) use directly and CIO-based
    # *  use indirectly.
    # *
    # *  This routine is part of the International Astronomical Union's
    # *  SOFA (Standards of Fundamental Astronomy) software collection.
    # *
    # *  Status:  support routine.
    # *
    # *  Given:
    # *     DATE1,DATE2   d       TT as a 2-part Julian Date (Note 1)
    # *     DPSI,DEPS     d       nutation (Note 2)
    # *
    # *  Returned:
    # *     EPSA          d       mean obliquity (Note 3)
    # *     RB          d(3,3)    frame bias matrix (Note 4)
    # *     RP          d(3,3)    precession matrix (Note 5)
    # *     RBP         d(3,3)    bias-precession matrix (Note 6)
    # *     RN          d(3,3)    nutation matrix (Note 7)
    # *     RBPN        d(3,3)    GCRS-to-true matrix (Note 8)
    # *
    # *  Notes:
    # *
    # *  1) The TT date DATE1+DATE2 is a Julian Date, apportioned in any
    # *     convenient way between the two arguments.  For example,
    # *     JD(TT)=2450123.7 could be expressed in any of these ways,
    # *     among others:
    # *
    # *            DATE1          DATE2
    # *
    # *         2450123.7.0        0.0        (JD method)
    # *          2451545.0      -1421.3.0     (J2000 method)
    # *         2400000.5.0     50123.2.0     (MJD method)
    # *         2450123.5.0       0.2.0       (date & time method)
    # *
    # *     The JD method is the most natural and convenient to use in
    # *     cases where the loss of several decimal digits of resolution
    # *     is acceptable.  The J2000 method is best matched to the way
    # *     the argument is handled internally and will deliver the
    # *     optimum resolution.  The MJD method and the date & time methods
    # *     are both good compromises between resolution and convenience.
    # *
    # *  2) The caller is responsible for providing the nutation components],
    # *     they are in longitude and obliquity, in radians and are with
    # *     respect to the equinox and ecliptic of date.  For high-accuracy
    # *     applications, free core nutation should be included as well as
    # *     any other relevant corrections to the position of the CIP.
    # *
    # *  3) The returned mean obliquity is consistent with the IAU 2000
    # *     precession-nutation remels.
    # *
    # *  4) The matrix RB transforms vectors from GCRS to J2000.0 mean equator
    # *     and equinox by applying frame bias.
    # *
    # *  5) The matrix RP transforms vectors from J2000.0 mean equator and
    # *     equinox to mean equator and equinox of date by applying
    # *     precession.
    # *
    # *  6) The matrix RBP transforms vectors from GCRS to mean equator and
    # *     equinox of date by applying frame bias then precession.  It is the
    # *     product RP x RB.
    # *
    # *  7) The matrix RN transforms vectors from mean equator and equinox of
    # *     date to true equator and equinox of date by applying the nutation
    # *     (luni-solar + planetary).
    # *
    # *  8) The matrix RBPN transforms vectors from GCRS to true equator and
    # *     equinox of date.  It is the product RN x RBP, applying frame bias,
    # *     precession and nutation in that order.
    # *
    # *  Called:
    # *     iau_PR00     IAU 2000 precession adjustments
    # *     iau_OBL80    mean obliquity, IAU 1980
    # *     iau_BP00     frame bias and precession matrices, IAU 2000
    # *     iau_NUMAT    form nutation matrix
    # *     iau_RXR      product of two r-matrices
    # *
    # *  Reference:
    # *
    # *     Capitaine, N., Chapront, J., Lambert, S. and Wallace, P.,
    # *     "Expressions for the Celestial Intermediate Pole and Celestial
    # *     Ephemeris Origin consistent with the IAU 2000A precession-nutation
    # *     remel", Astron.Astrophys. 400, 1145-1154 (2003)
    # *
    # *     n.b. The celestial ephemeris origin (CEO) was renamed "celestial
    # *          intermediate origin" (CIO) by IAU 2006 Resolution 2.
    # *
    # *  This revision:  2010 January 18
    # *
    # *  SOFA release 2010-12-01
    # *
    # *  Copyright (C) 2010 IAU SOFA Board.  See notes at end.
    # *
    # * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    '''
    # *  IAU 2000 precession-rate adjustments.
    _, DEPSPR = iau_PR00 ( DATE1, DATE2 )
    
    # *  Mean obliquity, consistent with IAU 2000 precession-nutation.
    EPSA = iau_OBL80 ( DATE1, DATE2 ) + DEPSPR
    
    # *  Frame bias and precession matrices and their product.
    RB, RP, RBP = iau_BP00 ( DATE1, DATE2 )

    # *  Nutation matrix.
    RN = iau_NUMAT ( EPSA, DPSI, DEPS )

    # *  Bias-precession-nutation matrix (classical).
    RBPN = dot(RN, RBP)
    
    return EPSA, RB, RP, RBP, RN, RBPN
    # *  Finished.



def iau_PNM00A ( DATE1, DATE2 ):
    '''
    # *+
    # *  - - - - - - - - - - -
    # *   i a u _ P N M 0 0 A
    # *  - - - - - - - - - - -
    # *
    # *  Form the matrix of precession-nutation for a given date (including
    # *  frame bias), equinox-based, IAU 2000A remel.
    # *
    # *  This routine is part of the International Astronomical Union's
    # *  SOFA (Standards of Fundamental Astronomy) software collection.
    # *
    # *  Status:  support routine.
    # *
    # *  Given:
    # *     DATE1,DATE2    d       TT as a 2-part Julian Date (Note 1)
    # *
    # *  Returned:
    # *     RBPN         d(3,3)    classical NPB matrix (Note 2)
    # *
    # *  Notes:
    # *
    # *  1) The TT date DATE1+DATE2 is a Julian Date, apportioned in any
    # *     convenient way between the two arguments.  For example,
    # *     JD(TT)=2450123.7 could be expressed in any of these ways,
    # *     among others:
    # *
    # *            DATE1          DATE2
    # *
    # *         2450123.7        0.0        (JD method)
    # *          2451545.0      -1421.3     (J2000 method)
    # *         2400000.5     50123.2     (MJD method)
    # *         2450123.5       0.2       (date & time method)
    # *
    # *     The JD method is the most natural and convenient to use in
    # *     cases where the loss of several decimal digits of resolution
    # *     is acceptable.  The J2000 method is best matched to the way
    # *     the argument is handled internally and will deliver the
    # *     optimum resolution.  The MJD method and the date & time methods
    # *     are both good compromises between resolution and convenience.
    # *
    # *  2) The matrix operates in the sense V(date) = RBPN * V(GCRS), where
    # *     the p-vector V(date) is with respect to the true equatorial triad
    # *     of date DATE1+DATE2 and the p-vector V(GCRS) is with respect to
    # *     the Geocentric Celestial Reference System (IAU, 2000).
    # *
    # *  3) A faster, but slightly less accurate result (about 1 mas), can be
    # *     obtained by using instead the iau_PNM00B routine.
    # *
    # *  Called:
    # *     iau_PN00A    bias/precession/nutation, IAU 2000A
    # *
    # *  Reference:
    # *
    # *     IAU: Trans. International Astronomical Union, Vol. XXIVB],  Proc.
    # *     24th General Assembly, Manchester, UK.  Resolutions B1.3, B1.6.
    # *     (2000)
    # *
    # *  This revision:  2009 December 21
    # *
    # *  SOFA release 2010-12-01
    # *
    # *  Copyright (C) 2010 IAU SOFA Board.  See notes at end.
    # *
    # *-----------------------------------------------------------------------
    '''
    # * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    
    # *  Obtain the required matrix (discarding other results).
    _, _, _, _, _, _, _, RBPN = iau_PN00A ( DATE1, DATE2 )
    return RBPN
    # *  Finished.


def iau_POM00 ( XP, YP, SP ):
    '''
    # *+
    # *  - - - - - - - - - - -
    # *   i a u _ P O M 0 0
    # *  - - - - - - - - - - -
    # *
    # *  Form the matrix of polar motion for a given date, IAU 2000.
    # *
    # *  This routine is part of the International Astronomical Union's
    # *  SOFA (Standards of Fundamental Astronomy) software collection.
    # *
    # *  Status:  support routine.
    # *
    # *  Given:
    # *     XP,YP      d      coordinates of the pole (radians, Note 1)
    # *     SP         d      the TIO locator s' (radians, Note 2)
    # *
    # *  Returned:
    # *     RPOM     d(3,3)   polar-motion matrix (Note 3)
    # *
    # *  Notes:
    # *
    # *  1) XP and YP are the coordinates (in radians) of the Celestial
    # *     Intermediate Pole with respect to the International Terrestrial
    # *     Reference System (see IERS Conventions 2003), measured along the
    # *     meridians to 0 and 90 deg west respectively.
    # *
    # *  2) SP is the TIO locator s', in radians, which positions the
    # *     Terrestrial Intermediate Origin on the equator.  It is obtained
    # *     from polar motion observations by numerical integration, and so is
    # *     in essence unpredictable.  However, it is dominated by a secular
    # *     drift of about 47 microarcseconds per century, and so can be taken
    # *     into account by using s' = -47*t, where t is centuries since
    # *     J2000.0.  The routine iau_SP00 implements this approximation.
    # *
    # *  3) The matrix operates in the sense V(TRS) = RPOM * V(CIP), meaning
    # *     that it is the final rotation when computing the pointing
    # *     direction to a celestial source.
    # *
    # *  Called:
    # *     iau_IR       initialize r-matrix to identity
    # *     iau_RZ       rotate around Z-axis
    # *     iau_RY       rotate around Y-axis
    # *     iau_RX       rotate around X-axis
    # *
    # *  Reference:
    # *
    # *     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    # *     IERS Technical Note No. 32, BKG (2004)
    # *
    # *  This revision:  2009 December 15
    # *
    # *  SOFA release 2010-12-01
    # *
    # *  Copyright (C) 2010 IAU SOFA Board.  See notes at end.
    # *
    # *-----------------------------------------------------------------------
    '''
    # *  Construct the matrix.
    RPOM = np.identity(3)
    RPOM = iau_RZ ( SP, RPOM )
    RPOM = iau_RY ( -XP, RPOM )
    RPOM = iau_RX ( -YP, RPOM )

    return RPOM
    # *  Finished.


#@jit(double[:,:](double, double[:,:]))
def rx (phi, R):
    '''
    # *+
    # *  - - - - - - -
    # *   i a u _ R X
    # *  - - - - - - -
    # *
    # *  Rotate an r-matrix about the x-axis.
    # *
    # *  This routine is part of the International Astronomical Union's
    # *  SOFA (Standards of Fundamental Astronomy) software collection.
    # *
    # *  Status:  vector/matrix support routine.
    # *
    # *  Given:
    # *     PHI      d         angle (radians)
    # *
    # *  Given and returned:
    # *     R        d(3,3)    r-matrix
    # *
    # *  Sign convention:  The matrix can be used to rotate the
    # *  reference frame of a vector.  Calling this routine with
    # *  positive PHI incorporates in the matrix an additional
    # *  rotation, about the x-axis, anticlockwise as seen looking
    # *  towards the origin from positive x.
    # *
    # *  This revision:  2006 November 13
    # *
    # *  SOFA release 2010-12-01
    # *
    # *  Copyright (C) 2010 IAU SOFA Board.  See notes at end.
    # *
    # *-----------------------------------------------------------------------
    '''
    # *  Matrix representing new rotation.
    S = sin(phi)
    C = cos(phi)
    A = np.identity(3)
    A[1,1] = C
    A[2,1] = -S
    A[1,2] = S
    A[2,2] = C

    # *  Rotate.
    return dot(A,R)

# *  Finished.

#@jit(double[:,:](double, double[:,:]))
def ry (phi, R):
    '''
    # *+
    # *  - - - - - - -
    # *   i a u _ R Y
    # *  - - - - - - -
    # *
    # *  Rotate an r-matrix about the x-axis.
    # *
    # *  This routine is part of the International Astronomical Union's
    # *  SOFA (Standards of Fundamental Astronomy) software collection.
    # *
    # *  Status:  vector/matrix support routine.
    # *
    # *  Given:
    # *     PHI      d         angle (radians)
    # *
    # *  Given and returned:
    # *     R        d(3,3)    r-matrix
    # *
    # *  Sign convention:  The matrix can be used to rotate the
    # *  reference frame of a vector.  Calling this routine with
    # *  positive PHI incorporates in the matrix an additional
    # *  rotation, about the x-axis, anticlockwise as seen looking
    # *  towards the origin from positive x.
    # *
    # *  This revision:  2006 November 13
    # *
    # *  SOFA release 2010-12-01
    # *
    # *  Copyright (C) 2010 IAU SOFA Board.  See notes at end.
    # *
    # *-----------------------------------------------------------------------
    '''
    # *  Matrix representing new rotation.
    S = sin(phi)
    C = cos(phi)
    A  = np.identity(3)
    A[0,0] = C
    A[2,0] = S
    A[0,2] = -S
    A[2,2] = C

    # *  Rotate.
    return dot(A,R)

    # *  Finished.

def iau_RX(phi, R):
    Rx = np.array([[ 1, 0, 0], \
                   [0, cos(phi), sin(phi)], \
                   [0, -sin(phi), cos(phi)]])
    return np.dot(Rx, R)
def iau_RY(phi, R):
    Ry = np.array([[cos(phi), 0, -sin(phi)], \
                   [0, 1, 0], \
                   [sin(phi), 0, cos(phi)]])
    return np.dot(Ry, R)
def iau_RZ(phi, R):
    Rz = np.array([[ cos(phi), sin(phi), 0], \
                   [-sin(phi), cos(phi), 0], \
                   [0,0,1]])
    return np.dot(Rz, R)
    
def rz ( phi, R ):
    '''
    # *+
    # *  - - - - - - -
    # *   i a u _ R Z
    # *  - - - - - - -
    # *
    # *  Rotate an r-matrix about the z-axis.
    # *
    # *  This routine is part of the International Astronomical Union's
    # *  SOFA (Standards of Fundamental Astronomy) software collection.
    # *
    # *  Status:  vector/matrix support routine.
    # *
    # *  Given:
    # *     PSI      d         angle (radians)
    # *
    # *  Given and returned:
    # *     R        d(3,3)    r-matrix, rotated
    # *
    # *  Sign convention:  The matrix can be used to rotate the
    # *  reference frame of a vector.  Calling this routine with
    # *  positive PSI incorporates in the matrix an additional
    # *  rotation, about the z-axis, anticlockwise as seen looking
    # *  towards the origin from positive z.
    # *
    # *  Called:
    # *     iau_IR       initialize r-matrix to identity
    # *     iau_RXR      product of two r-matrices
    # *     iau_CR       copy r-matrix
    # *
    # *  This revision:  2006 November 13
    # *
    # *  SOFA release 2010-12-01
    # *
    # *  Copyright (C) 2010 IAU SOFA Board.  See notes at end.
    # *
    # *-----------------------------------------------------------------------
    '''
    # *  Matrix representing new rotation.
    S = sin(phi)
    C = cos(phi)
    A = np.identity(3)
    A[0,0] = C
    A[1,0] = -S
    A[0,1] = S
    A[1,1] = C
    
    # *  Rotate.
    return dot(A,R)

# *  Finished.


def iau_SP00 ( DATE1, DATE2 ):
    '''
    # *+
    # *  - - - - - - - - -
    # *   i a u _ S P 0 0
    # *  - - - - - - - - -
    # *
    # *  The TIO locator s', positioning the Terrestrial Intermediate Origin
    # *  on the equator of the Celestial Intermediate Pole.
    # *
    # *  This routine is part of the International Astronomical Union's
    # *  SOFA (Standards of Fundamental Astronomy) software collection.
    # *
    # *  Status:  canonical remel.
    # *
    # *  Given:
    # *     DATE1,DATE2    d      TT as a 2-part Julian Date (Note 1)
    # *
    # *  Returned:
    # *     iau_SP00       d      the TIO locator s' in radians (Note 2)
    # *
    # *  Notes:
    # *
    # *  1) The TT date DATE1+DATE2 is a Julian Date, apportioned in any
    # *     convenient way between the two arguments.  For example,
    # *     JD(TT)=2450123.7 could be expressed in any of these ways,
    # *     among others:
    # *
    # *            DATE1          DATE2
    # *
    # *         2450123.7.0        0.0        (JD method)
    # *          2451545.0      -1421.3.0     (J2000 method)
    # *         2400000.5.0     50123.2.0     (MJD method)
    # *         2450123.5.0       0.2.0       (date & time method)
    # *
    # *     The JD method is the most natural and convenient to use in
    # *     cases where the loss of several decimal digits of resolution
    # *     is acceptable.  The J2000 method is best matched to the way
    # *     the argument is handled internally and will deliver the
    # *     optimum resolution.  The MJD method and the date & time methods
    # *     are both good compromises between resolution and convenience.
    # *
    # *  2) The TIO locator s' is obtained from polar motion observations by
    # *     numerical integration, and so is in essence unpredictable.
    # *     However, it is dominated by a secular drift of about
    # *     47 microarcseconds per century, which is the approximation
    # *     evaluated by the present routine.
    # *
    # *  Reference:
    # *
    # *     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    # *     IERS Technical Note No. 32, BKG (2004)
    # *
    # *  This revision:  2009 December 15
    # *
    # *  SOFA release 2010-12-01
    # *
    # *  Copyright (C) 2010 IAU SOFA Board.  See notes at end.
    # *
    # *-----------------------------------------------------------------------
    '''
    # *  Arcseconds to radians
    DAS2R = 4.848136811095359935899141e-6

    # *  Reference epoch (J2000.0), JD
    DJ00 = 2451545.0

    # *  Days per Julian century
    DJC = 36525.0

    # *  Time since J2000.0, in Julian centuries T

    # * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # *  Interval between fundamental epoch J2000.0 and current date (JC).
    T = ( ( DATE1-DJ00 ) + DATE2 ) / DJC

    # *  Approximate s'.
    return -47e-6 * T * DAS2R

    # *  Finished.


def iau_XYS00A ( DATE1, DATE2 ):
    '''
    # *+
    # *  - - - - - - - - - - -
    # *   i a u _ X Y S 0 0 A
    # *  - - - - - - - - - - -
    # *
    # *  For a given TT date, compute the X,Y coordinates of the Celestial
    # *  Intermediate Pole and the CIO locator s, using the IAU 2000A
    # *  precession-nutation remel.
    # *
    # *  This routine is part of the International Astronomical Union's
    # *  SOFA (Standards of Fundamental Astronomy) software collection.
    # *
    # *  Status:  support routine.
    # *
    # *  Given:
    # *     DATE1,DATE2   d    TT as a 2-part Julian Date (Note 1)
    # *
    # *  Returned:
    # *     X,Y           d    Celestial Intermediate Pole (Note 2)
    # *     S             d    the CIO locator s (Note 2)
    # *
    # *  Notes:
    # *
    # *  1) The TT date DATE1+DATE2 is a Julian Date, apportioned in any
    # *     convenient way between the two arguments.  For example,
    # *     JD(TT)=2450123.7 could be expressed in any of these ways,
    # *     among others:
    # *
    # *            DATE1          DATE2
    # *
    # *         2450123.7.0        0.0        (JD method)
    # *          2451545.0      -1421.3.0     (J2000 method)
    # *         2400000.5.0     50123.2.0     (MJD method)
    # *         2450123.5.0       0.2.0       (date & time method)
    # *
    # *     The JD method is the most natural and convenient to use in
    # *     cases where the loss of several decimal digits of resolution
    # *     is acceptable.  The J2000 method is best matched to the way
    # *     the argument is handled internally and will deliver the
    # *     optimum resolution.  The MJD method and the date & time methods
    # *     are both good compromises between resolution and convenience.
    # *
    # *  2) The Celestial Intermediate Pole coordinates are the x,y components
    # *     of the unit vector in the Geocentric Celestial Reference System.
    # *
    # *  3) The CIO locator s (in radians) positions the Celestial
    # *     Intermediate Origin on the equator of the CIP.
    # *
    # *  4) A faster, but slightly less accurate result (about 1 mas for X,Y),
    # *     can be obtained by using instead the iau_XYS00B routine.
    # *
    # *  Called:
    # *     iau_PNM00A   classical NPB matrix, IAU 2000A
    # *     iau_BPN2XY   extract CIP X,Y coordinates from NPB matrix
    # *     iau_S00      the CIO locator s, given X,Y, IAU 2000A
    # *
    # *  Reference:
    # *
    # *     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    # *     IERS Technical Note No. 32, BKG (2004)
    # *
    # *  This revision:  2006 November 13
    # *
    # *  SOFA release 2010-12-01
    # *
    # *  Copyright (C) 2010 IAU SOFA Board.  See notes at end.
    # *
    # *-----------------------------------------------------------------------
    '''
    # *  Form the bias-precession-nutation matrix, IAU 2000A.
    RBPN = iau_PNM00A ( DATE1, DATE2 )

    # *  Extract X,Y.
    X = RBPN[2,0]
    Y = RBPN[2,1]

    # *  Obtain s.
    S = iau_S00 ( DATE1, DATE2, X, Y )
    
    return X, Y, S
    # *  Finished.


def iau_PR00 ( DATE1, DATE2 ):
    '''
    # *+
    # *  - - - - - - - - -
    # *   i a u _ P R 0 0
    # *  - - - - - - - - -
    # *
    # *  Precession-rate part of the IAU 2000 precession-nutation remels
    # *  (part of MHB2000).
    # *
    # *  This routine is part of the International Astronomical Union's
    # *  SOFA (Standards of Fundamental Astronomy) software collection.
    # *
    # *  Status:  canonical remel.
    # *
    # *  Given:
    # *     DATE1,DATE2    d   TT as a 2-part Julian Date (Note 1)
    # *
    # *  Returned:
    # *     DPSIPR,DEPSPR  d   precession corrections (Notes 2,3)
    # *
    # *  Notes:
    # *
    # *  1) The TT date DATE1+DATE2 is a Julian Date, apportioned in any
    # *     convenient way between the two arguments.  For example,
    # *     JD(TT)=2450123.7 could be expressed in any of these ways,
    # *     among others
    # *
    # *            DATE1          DATE2
    # *
    # *         2450123.7.0        0.0        (JD method)
    # *          2451545.0      -1421.3.0     (J2000 method)
    # *         2400000.5.0     50123.2.0     (MJD method)
    # *         2450123.5.0       0.2.0       (date & time method)
    # *
    # *     The JD method is the most natural and convenient to use in
    # *     cases where the loss of several decimal digits of resolution
    # *     is acceptable.  The J2000 method is best matched to the way
    # *     the argument is handled internally and will deliver the
    # *     optimum resolution.  The MJD method and the date & time methods
    # *     are both good compromises between resolution and convenience.
    # *
    # *  2) The precession adjustments are expressed as "nutation components",
    # *     corrections in longitude and obliquity with respect to the J2000.0
    # *     equinox and ecliptic.
    # *
    # *  3) Although the precession adjustments are stated to be with respect
    # *     to Lieske et al. (1977), the MHB2000 remel does not specify which
    # *     set of Euler angles are to be used and how the adjustments are to
    # *     be applied.  The most literal and straightforward procedure is to
    # *     adopt the 4-rotation epsilon_0, psi_A, omega_A, xi_A option, and
    # *     to add DPSIPR to psi_A and DEPSPR to both omega_A and eps_A.
    # *
    # *  4) This is an implementation of one aspect of the IAU 2000A nutation
    # *     remel, formally adopted by the IAU General Assembly in 2000,
    # *     namely MHB2000 (Mathews et al. 2002).
    # *
    # *  References:
    # *
    # *     Lieske, J.H., Lederle, T., Fricke, W. & Morando, B., "Expressions
    # *     for the precession quantities based upon the IAU (1976) System of
    # *     Astronomical Constants", Astron.Astrophys., 58, 1-16 (1977)
    # *
    # *     Mathews, P.M., Herring, T.A., Buffet, B.A., "remeling of nutation
    # *     and precession   New nutation series for nonrigid Earth and
    # *     insights into the Earth's interior", J.Geophys.Res., 107, B4,
    # *     2002.  The MHB2000 code itself was obtained on 9th September 2002
    # *     from ftp://maia.usno.navy.mil/conv2000/chapter5/IAU2000A.
    # *
    # *     Wallace, P.T., "Software for Implementing the IAU 2000
    # *     Resolutions", in IERS Workshop 5.1 (2002).
    # *
    # *  This revision:  2009 December 15
    # *
    # *  SOFA release 2010-12-01
    # *
    # *  Copyright (C) 2010 IAU SOFA Board.  See notes at end.
    # *
    # *-----------------------------------------------------------------------
    '''
    # *  Arcseconds to radians
    DAS2R = 4.848136811095359935899141e-6

    # *  Reference epoch (J2000.0), JD
    DJ00 = 2451545.0

    # *  Days per Julian century
    DJC = 36525.0

    # *  ------------------------------------
    # *  Precession and obliquity corrections (radians per century)
    # *  ------------------------------------
    
    PRECOR = -0.29965 * DAS2R
    OBLCOR = -0.02524 * DAS2R
    
    # *  Interval between fundamental epoch J2000.0 and given date (JC).
    T = ( ( DATE1-DJ00 ) + DATE2 ) / DJC

    # *  Precession rate contributions with respect to IAU 1976/80.
    DPSIPR = PRECOR * T
    DEPSPR = OBLCOR * T
    
    return DPSIPR, DEPSPR
    # *  Finished.



def iau_OBL80 ( DATE1, DATE2 ):
    '''
    # *+
    # *  - - - - - - - - - -
    # *   i a u _ O B L 8 0
    # *  - - - - - - - - - -
    # *
    # *  Mean obliquity of the ecliptic, IAU 1980 remel.
    # *
    # *  This routine is part of the International Astronomical Union's
    # *  SOFA (Standards of Fundamental Astronomy) software collection.
    # *
    # *  Status:  canonical remel.
    # *
    # *  Given:
    # *     DATE1,DATE2     d      TT as a 2-part Julian Date (Note 1)
    # *
    # *  Returned:
    # *     iau_OBL80       d      obliquity of the ecliptic (radians, Note 2)
    # *
    # *  Notes:
    # *
    # *  1) The date DATE1+DATE2 is a Julian Date, apportioned in any
    # *     convenient way between the two arguments.  For example,
    # *     JD(TDB)=2450123.7 could be expressed in any of these ways,
    # *     among others:
    # *
    # *             DATE1         DATE2
    # *
    # *         2450123.7.0        0.0        (JD method)
    # *          2451545.0      -1421.3.0     (J2000 method)
    # *         2400000.5.0     50123.2.0     (MJD method)
    # *         2450123.5.0       0.2.0       (date & time method)
    # *
    # *     The JD method is the most natural and convenient to use in
    # *     cases where the loss of several decimal digits of resolution
    # *     is acceptable.  The J2000 method is best matched to the way
    # *     the argument is handled internally and will deliver the
    # *     optimum resolution.  The MJD method and the date & time methods
    # *     are both good compromises between resolution and convenience.
    # *
    # *  2) The result is the angle between the ecliptic and mean equator of
    # *     date DATE1+DATE2.
    # *
    # *  Reference:
    # *
    # *     Explanatory Supplement to the Astronomical Almanac,
    # *     P. Kenneth Seidelmann (ed), University Science Books (1992),
    # *     Expression 3.222-1 (p114).
    # *
    # *  This revision:  2009 December 15
    # *
    # *  SOFA release 2010-12-01
    # *
    # *  Copyright (C) 2010 IAU SOFA Board.  See notes at end.
    # *
    # *-----------------------------------------------------------------------
    '''
    # *  Arcseconds to radians
    DAS2R = 4.848136811095359935899141e-6

    # *  Reference epoch (J2000.0), JD
    DJ00 = 2451545.0

    # *  Days per Julian century
    DJC = 36525.0

    # *  Interval between fundamental epoch J2000.0 and given date (JC).
    T = ( ( DATE1-DJ00 ) + DATE2 ) / DJC

    # *  Mean obliquity of date.
    return DAS2R * ( 84381.448 + ( -46.8150 + ( -0.00059 + 
                               0.001813 * T ) * T ) * T )

    # *  Finished.



def iau_NUMAT ( EPSA, DPSI, DEPS ):
    '''
    # *+
    # *  - - - - - - - - - -
    # *   i a u _ N U M A T
    # *  - - - - - - - - - -
    # *
    # *  Form the matrix of nutation.
    # *
    # *  This routine is part of the International Astronomical Union's
    # *  SOFA (Standards of Fundamental Astronomy) software collection.
    # *
    # *  Status:  support routine.
    # *
    # *  Given:
    # *     EPSA          d       mean obliquity of date (Note 1)
    # *     DPSI,DEPS     d       nutation (Note 2)
    # *
    # *  Returned:
    # *     RMATN       d(3,3)    nutation matrix (Note 3)
    # *
    # *  Notes:
    # *
    # *
    # *  1) The supplied mean obliquity EPSA, must be consistent with the
    # *     precession-nutation remels from which DPSI and DEPS were obtained.
    # *
    # *  2) The caller is responsible for providing the nutation components],
    # *     they are in longitude and obliquity, in radians and are with
    # *     respect to the equinox and ecliptic of date.
    # *
    # *  3) The matrix operates in the sense V(true) = RMATN * V(mean),
    # *     where the p-vector V(true) is with respect to the true
    # *     equatorial triad of date and the p-vector V(mean) is with
    # *     respect to the mean equatorial triad of date.
    # *
    # *  Called:
    # *     iau_IR       initialize r-matrix to identity
    # *     iau_RX       rotate around X-axis
    # *     iau_RZ       rotate around Z-axis
    # *
    # *  Reference:
    # *
    # *     Explanatory Supplement to the Astronomical Almanac,
    # *     P. Kenneth Seidelmann (ed), University Science Books (1992),
    # *     Section 3.222-3 (p114).
    # *
    # *  This revision:  2006 November 13
    # *
    # *  SOFA release 2010-12-01
    # *
    # *  Copyright (C) 2010 IAU SOFA Board.  See notes at end.
    # *
    # *-----------------------------------------------------------------------
    '''
    # *  Build the rotation matrix.
    RMATN = np.identity(3)
    RMATN = iau_RX ( EPSA, RMATN )
    RMATN = iau_RZ ( -DPSI, RMATN )
    RMATN = iau_RX ( -(EPSA+DEPS), RMATN )

    return RMATN
    # *  Finished.


def iau_BP00 ( DATE1, DATE2 ):
    '''
    # *+
    # *  - - - - - - - - -
    # *   i a u _ B P 0 0
    # *  - - - - - - - - -
    # *
    # *  Frame bias and precession, IAU 2000.
    # *
    # *  This routine is part of the International Astronomical Union's
    # *  SOFA (Standards of Fundamental Astronomy) software collection.
    # *
    # *  Status:  canonical remel.
    # *
    # *  Given:
    # *     DATE1,DATE2    d       TT as a 2-part Julian Date (Note 1)
    # *
    # *  Returned:
    # *     RB           d(3,3)    frame bias matrix (Note 2)
    # *     RP           d(3,3)    precession matrix (Note 3)
    # *     RBP          d(3,3)    bias-precession matrix (Note 4)
    # *
    # *  Notes:
    # *
    # *  1) The TT date DATE1+DATE2 is a Julian Date, apportioned in any
    # *     convenient way between the two arguments.  For example,
    # *     JD(TT)=2450123.7 could be expressed in any of these ways,
    # *     among others:
    # *
    # *            DATE1          DATE2
    # *
    # *         2450123.7.0        0.0        (JD method)
    # *          2451545.0      -1421.3.0     (J2000 method)
    # *         2400000.5.0     50123.2.0     (MJD method)
    # *         2450123.5.0       0.2.0       (date & time method)
    # *
    # *     The JD method is the most natural and convenient to use in
    # *     cases where the loss of several decimal digits of resolution
    # *     is acceptable.  The J2000 method is best matched to the way
    # *     the argument is handled internally and will deliver the
    # *     optimum resolution.  The MJD method and the date & time methods
    # *     are both good compromises between resolution and convenience.
    # *
    # *  2) The matrix RB transforms vectors from GCRS to mean J2000.0 by
    # *     applying frame bias.
    # *
    # *  3) The matrix RP transforms vectors from J2000.0 mean equator and
    # *     equinox to mean equator and equinox of date by applying
    # *     precession.
    # *
    # *  4) The matrix RBP transforms vectors from GCRS to mean equator and
    # *     equinox of date by applying frame bias then precession.  It is the
    # *     product RP x RB.
    # *
    # *  Called:
    # *     iau_BI00     frame bias components, IAU 2000
    # *     iau_PR00     IAU 2000 precession adjustments
    # *     iau_IR       initialize r-matrix to identity
    # *     iau_RX       rotate around X-axis
    # *     iau_RY       rotate around Y-axis
    # *     iau_RZ       rotate around Z-axis
    # *     iau_RXR      product of two r-matrices
    # *
    # *  Reference:
    # *
    # *     Capitaine, N., Chapront, J., Lambert, S. and Wallace, P.,
    # *     "Expressions for the Celestial Intermediate Pole and Celestial
    # *     Ephemeris Origin consistent with the IAU 2000A precession-nutation
    # *     remel", Astron.Astrophys. 400, 1145-1154 (2003)
    # *
    # *     n.b. The celestial ephemeris origin (CEO) was renamed "celestial
    # *          intermediate origin" (CIO) by IAU 2006 Resolution 2.
    # *
    # *  This revision:  2010 January 18
    # *
    # *  SOFA release 2010-12-01
    # *
    # *  Copyright (C) 2010 IAU SOFA Board.  See notes at end.
    # *
    # *-----------------------------------------------------------------------
    '''
    # *  Arcseconds to radians
    DAS2R = 4.848136811095359935899141e-6

    # *  Reference epoch (J2000.0), JD
    DJ00 = 2451545.0

    # *  Days per Julian century
    DJC = 36525.0

    # *  J2000.0 obliquity (Lieske et al. 1977)
    EPS0 = 84381.448 * DAS2R

    # *  Interval between fundamental epoch J2000.0 and current date (JC).
    T = ( ( DATE1-DJ00 ) + DATE2 ) / DJC

    # *  Frame bias.
    DPSIBI, DEPSBI, DRA0 = iau_BI00()

    # *  Precession angles (Lieske et al. 1977)
    PSIA77 =          ( 5038.7784 + 
                      (   -1.07259 + 
                      (   -0.001147 ) * T ) * T ) * T * DAS2R
    OMA77  =   EPS0 + ( 
                      (    0.05127 + 
                      (   -0.007726 ) * T ) * T ) * T * DAS2R
    CHIA   =          (   10.5526 + 
                      (   -2.38064 + 
                      (   -0.001125 ) * T ) * T ) * T * DAS2R

    # *  Apply IAU 2000 precession corrections.
    DPSIPR, DEPSPR = iau_PR00 ( DATE1, DATE2 )
    PSIA = PSIA77 + DPSIPR
    OMA  = OMA77  + DEPSPR

    # *  Frame bias matrix: GCRS to J2000.0.
    RB = np.identity(3)
    RB = iau_RZ ( DRA0, RB )
    RB = iau_RY ( DPSIBI*sin(EPS0), RB )
    RB = iau_RX ( -DEPSBI, RB )

    # *  Precession matrix: J2000.0 to mean of date.
    RP = np.identity(3)
    RP = iau_RX ( EPS0, RP )
    RP = iau_RZ ( -PSIA, RP )
    RP = iau_RX ( -OMA, RP )
    RP = iau_RZ ( CHIA, RP )

    # *  Bias-precession matrix: GCRS to mean of date.
    RBP = dot(RP, RB)

    return RB, RP, RBP
    # *  Finished.


def iau_BI00():
    '''
    # *+
    # *  - - - - - - - - -
    # *   i a u _ B I 0 0
    # *  - - - - - - - - -
    # *
    # *  Frame bias components of IAU 2000 precession-nutation remels (part of
    # *  MHB2000 with additions).
    # *
    # *  This routine is part of the International Astronomical Union's
    # *  SOFA (Standards of Fundamental Astronomy) software collection.
    # *
    # *  Status:  canonical remel.
    # *
    # *  Returned:
    # *     DPSIBI,DEPSBI  d   longitude and obliquity corrections
    # *     DRA            d   the ICRS RA of the J2000.0 mean equinox
    # *
    # *  Notes:
    # *
    # *  1) The frame bias corrections in longitude and obliquity (radians)
    # *     are required in order to correct for the offset between the GCRS
    # *     pole and the J2000.0 mean pole.  They define, with respect to the
    # *     GCRS frame, a J2000.0 mean pole that is consistent with the rest
    # *     of the IAU 2000A precession-nutation remel.
    # *
    # *  2) In addition to the displacement of the pole, the complete
    # *     description of the frame bias requires also an offset in right
    # *     ascension.  This is not part of the IAU 2000A remel, and is from
    # *     Chapront et al. (2002).  It is returned in radians.
    # *
    # *  3) This is a supplemented implementation of one aspect of the IAU
    # *     2000A nutation remel, formally adopted by the IAU General Assembly
    # *     in 2000, namely MHB2000 (Mathews et al. 2002).
    # *
    # *  References:
    # *
    # *     Chapront, J., Chapront-Touze, M. & Francou, G., Astron.Astrophys.,
    # *     387, 700, 2002.
    # *
    # *     Mathews, P.M., Herring, T.A., Buffet, B.A., "remeling of nutation
    # *     and precession   New nutation series for nonrigid Earth and
    # *     insights into the Earth's interior", J.Geophys.Res., 107, B4,
    # *     2002.  The MHB2000 code itself was obtained on 9th September 2002
    # *     from ftp://maia.usno.navy.mil/conv2000/chapter5/IAU2000A.
    # *
    # *  This revision:  2009 December 15
    # *
    # *  SOFA release 2010-12-01
    # *
    # *  Copyright (C) 2010 IAU SOFA Board.  See notes at end.
    # *
    # *-----------------------------------------------------------------------
    '''
    # *  Arcseconds to radians
    DAS2R = 4.848136811095359935899141e-6

    # *  The frame bias corrections in longitude and obliquity
    DPBIAS = -0.041775 * DAS2R
    DEBIAS = -0.0068192 * DAS2R

    # *  The ICRS RA of the J2000.0 equinox (Chapront et al., 2002)
    DRA0 = -0.0146 * DAS2R

    # *  Return the results (which are fixed).
    DPSIBI = DPBIAS
    DEPSBI = DEBIAS
    DRA = DRA0
    
    return DPSIBI, DEPSBI, DRA
    # *  Finished.

#@jit
def iau_S00 ( DATE1, DATE2, X, Y ):
    '''
    # *+
    # *  - - - - - - - -
    # *   i a u _ S 0 0
    # *  - - - - - - - -
    # *
    # *  The CIO locator s, positioning the Celestial Intermediate Origin on
    # *  the equator of the Celestial Intermediate Pole, given the CIP's X,Y
    # *  coordinates.  Compatible with IAU 2000A precession-nutation.
    # *
    # *  This routine is part of the International Astronomical Union's
    # *  SOFA (Standards of Fundamental Astronomy) software collection.
    # *
    # *  Status:  canonical model.
    # *
    # *  Given:
    # *     DATE1,DATE2    d      TT as a 2-part Julian Date (Note 1)
    # *     X,Y            d      CIP coordinates (Note 3)
    # *
    # *  Returned:
    # *     iau_S00        d      the CIO locator s in radians (Note 2)
    # *
    # *  Notes:
    # *
    # *  1) The TT date DATE1+DATE2 is a Julian Date, apportioned in any
    # *     convenient way between the two arguments.  For example,
    # *     JD(TT)=2450123.7 could be expressed in any of these ways,
    # *     among others:
    # *
    # *            DATE1          DATE2
    # *
    # *         2450123.7.0        0.0        (JD method)
    # *          2451545.0      -1421.3.0     (J2000 method)
    # *         2400000.5.0     50123.2.0     (MJD method)
    # *         2450123.5.0       0.2.0       (date & time method)
    # *
    # *     The JD method is the most natural and convenient to use in
    # *     cases where the loss of several decimal digits of resolution
    # *     is acceptable.  The J2000 method is best matched to the way
    # *     the argument is handled internally and will deliver the
    # *     optimum resolution.  The MJD method and the date & time methods
    # *     are both good compromises between resolution and convenience.
    # *
    # *  2) The CIO locator s is the difference between the right ascensions
    # *     of the same point in two systems:  the two systems are the GCRS
    # *     and the CIP,CIO, and the point is the ascending node of the
    # *     CIP equator.  The quantity s remains below 0.1 arcsecond
    # *     throughout 1900-2100.
    # *
    # *  3) The series used to compute s is in fact for s+XY/2, where X and Y
    # *     are the x and y components of the CIP unit vector],  this series is
    # *     more compact than a direct series for s would be.  This routine
    # *     requires X,Y to be supplied by the caller, who is responsible for
    # *     providing values that are consistent with the supplied date.
    # *
    # *  4) The model is consistent with the IAU 2000A precession-nutation.
    # *
    # *  Called:
    # *     iau_FAL03    mean anomaly of the Moon
    # *     iau_FALP03   mean anomaly of the Sun
    # *     iau_FAF03    mean argument of the latitude of the Moon
    # *     iau_FA.03    mean elongation of the Moon from the Sun
    # *     iau_FAOM03   mean longitude of the Moon's ascending node
    # *     iau_FAVE03   mean longitude of Venus
    # *     iau_FAE03    mean longitude of Earth
    # *     iau_FAPA03   general accumulated precession in longitude
    # *
    # *  References:
    # *
    # *     Capitaine, N., Chapront, J., Lambert, S. and Wallace, P.,
    # *     "Expressions for the Celestial Intermediate Pole and Celestial
    # *     Ephemeris Origin consistent with the IAU 2000A precession-nutation
    # *     model", Astron.Astrophys. 400, 1145-1154 (2003)
    # *
    # *     n.b. The celestial ephemeris origin (CEO) was renamed "celestial
    # *          intermediate origin" (CIO) by IAU 2006 Resolution 2.
    # *
    # *     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    # *     IERS Technical Note No. 32, BKG (2004)
    # *
    # *  This revision:  2010 January 18
    # *
    # *  SOFA release 2010-12-01
    # *
    # *  Copyright (C) 2010 IAU SOFA Board.  See notes at end.
    # *
    # *-----------------------------------------------------------------------
    '''
    # *  Arcseconds to radians
    DAS2R = 4.848136811095359935899141e-6

    # *  Reference epoch (J2000.0), JD
    DJ00 = 2451545.0

    # *  Days per Julian century
    DJC = 36525.0

    # *  ---------------------
    # *  The series for s+XY/2
    # *  ---------------------

    # *  Number of terms in the series
    NS0=33
    NS1=3
    NS2=25
    NS3=4

    # *  Polynomial coefficients
    SP = [94e-6, 3808.35e-6, -119.94e-6, -72574.09e-6, 27.70e-6, 15.61e-6 ]

    # *  Argument coefficients for t^0
    KS0 = np.asarray([
      [ 0,  0,  0,  0,  1,  0,  0,  0],
      [ 0,  0,  0,  0,  2,  0,  0,  0],
      [ 0,  0,  2, -2,  3,  0,  0,  0],
      [ 0,  0,  2, -2,  1,  0,  0,  0],
      [ 0,  0,  2, -2,  2,  0,  0,  0],
      [ 0,  0,  2,  0,  3,  0,  0,  0],
      [ 0,  0,  2,  0,  1,  0,  0,  0],
      [ 0,  0,  0,  0,  3,  0,  0,  0],
      [ 0,  1,  0,  0,  1,  0,  0,  0],
      [ 0,  1,  0,  0, -1,  0,  0,  0],
      [ 1,  0,  0,  0, -1,  0,  0,  0],
      [ 1,  0,  0,  0,  1,  0,  0,  0],
      [ 0,  1,  2, -2,  3,  0,  0,  0],
      [ 0,  1,  2, -2,  1,  0,  0,  0],
      [ 0,  0,  4, -4,  4,  0,  0,  0],
      [ 0,  0,  1, -1,  1, -8, 12,  0],
      [ 0,  0,  2,  0,  0,  0,  0,  0],
      [ 0,  0,  2,  0,  2,  0,  0,  0],
      [ 1,  0,  2,  0,  3,  0,  0,  0],
      [ 1,  0,  2,  0,  1,  0,  0,  0],
      [ 0,  0,  2, -2,  0,  0,  0,  0],
      [ 0,  1, -2,  2, -3,  0,  0,  0],
      [ 0,  1, -2,  2, -1,  0,  0,  0],
      [ 0,  0,  0,  0,  0,  8,-13, -1],
      [ 0,  0,  0,  2,  0,  0,  0,  0],
      [ 2,  0, -2,  0, -1,  0,  0,  0],
      [ 0,  1,  2, -2,  2,  0,  0,  0],
      [ 1,  0,  0, -2,  1,  0,  0,  0],
      [ 1,  0,  0, -2, -1,  0,  0,  0],
      [ 0,  0,  4, -2,  4,  0,  0,  0],
      [ 0,  0,  2, -2,  4,  0,  0,  0],
      [ 1,  0, -2,  0, -3,  0,  0,  0],
      [ 1,  0, -2,  0, -1,  0,  0,  0] ])

    # *  Argument coefficients for t^1
    KS1 = np.asarray([
      [ 0,  0,  0,  0,  2,  0,  0,  0],
      [ 0,  0,  0,  0,  1,  0,  0,  0],
      [ 0,  0,  2, -2,  3,  0,  0,  0] ])

    # *  Argument coefficients for t^2
    KS2 = np.asarray([
      [ 0,  0,  0,  0,  1,  0,  0,  0],
      [ 0,  0,  2, -2,  2,  0,  0,  0],
      [ 0,  0,  2,  0,  2,  0,  0,  0],
      [ 0,  0,  0,  0,  2,  0,  0,  0],
      [ 0,  1,  0,  0,  0,  0,  0,  0],
      [ 1,  0,  0,  0,  0,  0,  0,  0],
      [ 0,  1,  2, -2,  2,  0,  0,  0],
      [ 0,  0,  2,  0,  1,  0,  0,  0],
      [ 1,  0,  2,  0,  2,  0,  0,  0],
      [ 0,  1, -2,  2, -2,  0,  0,  0],
      [ 1,  0,  0, -2,  0,  0,  0,  0],
      [ 0,  0,  2, -2,  1,  0,  0,  0],
      [ 1,  0, -2,  0, -2,  0,  0,  0],
      [ 0,  0,  0,  2,  0,  0,  0,  0],
      [ 1,  0,  0,  0,  1,  0,  0,  0],
      [ 1,  0, -2, -2, -2,  0,  0,  0],
      [ 1,  0,  0,  0, -1,  0,  0,  0],
      [ 1,  0,  2,  0,  1,  0,  0,  0],
      [ 2,  0,  0, -2,  0,  0,  0,  0],
      [ 2,  0, -2,  0, -1,  0,  0,  0],
      [ 0,  0,  2,  2,  2,  0,  0,  0],
      [ 2,  0,  2,  0,  2,  0,  0,  0],
      [ 2,  0,  0,  0,  0,  0,  0,  0],
      [ 1,  0,  2, -2,  2,  0,  0,  0],
      [ 0,  0,  2,  0,  0,  0,  0,  0] ])

    # *  Argument coefficients for t^3
    KS3 = np.asarray([
      [ 0,  0,  0,  0,  1,  0,  0,  0],
      [ 0,  0,  2, -2,  2,  0,  0,  0],
      [ 0,  0,  2,  0,  2,  0,  0,  0],
      [ 0,  0,  0,  0,  2,  0,  0,  0] ])

    # *  Argument coefficients for t^4
    KS4 = np.asarray([ 0,  0,  0,  0,  1,  0,  0,  0])

# *  sine and cosine coefficients for t^0
    SS0 = np.asarray([
                 [-2640.73e-6,          +0.39e-6],
                 [  -63.53e-6,          +0.02e-6],
                 [  -11.75e-6,          -0.01e-6],
                 [  -11.21e-6,          -0.01e-6],
                 [   +4.57e-6,           0.00e-6],
                 [   -2.02e-6,           0.00e-6],
                 [   -1.98e-6,           0.00e-6],
                 [   +1.72e-6,           0.00e-6],
                 [   +1.41e-6,          +0.01e-6],
                 [   +1.26e-6,          +0.01e-6],
                 [   +0.63e-6,           0.00e-6],
                 [   +0.63e-6,           0.00e-6],
                 [   -0.46e-6,           0.00e-6],
                 [   -0.45e-6,           0.00e-6],
                 [   -0.36e-6,           0.00e-6],
                 [   +0.24e-6,          +0.12e-6],
                 [   -0.32e-6,           0.00e-6],
                 [   -0.28e-6,           0.00e-6],
                 [   -0.27e-6,           0.00e-6],
                 [   -0.26e-6,           0.00e-6],
                 [   +0.21e-6,           0.00e-6],
                 [   -0.19e-6,           0.00e-6],
                 [   -0.18e-6,           0.00e-6],
                 [   +0.10e-6,          -0.05e-6],
                 [   -0.15e-6,           0.00e-6],
                 [   +0.14e-6,           0.00e-6],
                 [   +0.14e-6,           0.00e-6],
                 [   -0.14e-6,           0.00e-6],
                 [   -0.14e-6,           0.00e-6],
                 [   -0.13e-6,           0.00e-6],
                 [   +0.11e-6,           0.00e-6],
                 [   -0.11e-6,           0.00e-6],
                 [   -0.11e-6,           0.00e-6] ])

    # *  sine and cosine coefficients for t^1
    SS1 = np.asarray([
                 [   -0.07e-6,          +3.57e-6],
                 [   +1.71e-6,          -0.03e-6],
                 [    0.00e-6,          +0.48e-6] ])

    # *  sine and cosine coefficients for t^2
    SS2 = np.asarray([
                 [ +743.53e-6,          -0.17e-6],
                 [  +56.91e-6,          +0.06e-6],
                 [   +9.84e-6,          -0.01e-6],
                 [   -8.85e-6,          +0.01e-6],
                 [   -6.38e-6,          -0.05e-6],
                 [   -3.07e-6,           0.00e-6],
                 [   +2.23e-6,           0.00e-6],
                 [   +1.67e-6,           0.00e-6],
                 [   +1.30e-6,           0.00e-6],
                 [   +0.93e-6,           0.00e-6],
                 [   +0.68e-6,           0.00e-6],
                 [   -0.55e-6,           0.00e-6],
                 [   +0.53e-6,           0.00e-6],
                 [   -0.27e-6,           0.00e-6],
                 [   -0.27e-6,           0.00e-6],
                 [   -0.26e-6,           0.00e-6],
                 [   -0.25e-6,           0.00e-6],
                 [   +0.22e-6,           0.00e-6],
                 [   -0.21e-6,           0.00e-6],
                 [   +0.20e-6,           0.00e-6],
                 [   +0.17e-6,           0.00e-6],
                 [   +0.13e-6,           0.00e-6],
                 [   -0.13e-6,           0.00e-6],
                 [   -0.12e-6,           0.00e-6],
                 [   -0.11e-6,           0.00e-6] ])

    # *  sine and cosine coefficients for t^3
    SS3 =np.asarray([
                 [   +0.30e-6,         -23.51e-6],
                 [   -0.03e-6,          -1.39e-6],
                 [   -0.01e-6,          -0.24e-6],
                 [    0.00e-6,          +0.22e-6] ])

    # *  sine and cosine coefficients for t^4
    SS4 = np.asarray([ -0.26e-6, -0.01e-6 ])

    # * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # *  Interval between fundamental epoch J2000.0 and current date (JC).
    T = ( ( DATE1-DJ00 ) + DATE2 ) / DJC

    # *  Fundamental Arguments (from IERS Conventions 2003)
    
    FA = np.zeros(8)
    
    # *  Mean anomaly of the Moon.
    FA[0] = iau_FAL03 ( T )

    # *  Mean anomaly of the Sun.
    FA[1] = iau_FALP03 ( T )

    # *  Mean longitude of the Moon minus that of the ascending node.
    FA[2] = iau_FAF03 ( T )

    # *  Mean elongation of the Moon from the Sun.
    FA[3] = iau_FAD03 ( T )

    # *  Mean longitude of the ascending node of the Moon.
    FA[4] = iau_FAOM03 ( T )

    # *  Mean longitude of Venus.
    FA[5] = iau_FAVE03 ( T )

    # *  Mean longitude of Earth.
    FA[6] = iau_FAE03 ( T )

    # *  General precession in longitude.
    FA[7] = iau_FAPA03 ( T )

    # *  Evaluate s.
    S0 = SP[0]
    S1 = SP[1]
    S2 = SP[2]
    S3 = SP[3]
    S4 = SP[4]
    S5 = SP[5]

    for I in range(NS0-1,-1,-1):
        A = 0.0
        for J in range(0,8):
            A = A + float(KS0[I,J])*FA[J]
        S0 = S0 + ( SS0[I,0]*sin(A) + SS0[I,1]*cos(A) )
        
    for I in range(NS1-1,-1,-1):
        A = 0.0
        for J in range(0,8):
            A = A + float(KS1[I,J])*FA[J]
        S1 = S1 + ( SS1[I,0]*sin(A) + SS1[I,1]*cos(A) )

    for I in range(NS2-1,-1,-1):
        A = 0.0
        for J in range(0,8):
            A = A + float(KS2[I,J])*FA[J]
        S2 = S2 + ( SS2[I,0]*sin(A) + SS2[I,1]*cos(A) )

    for I in range(NS3-1,-1,-1):
        A = 0.0
        for J in range(0,8):
            A = A + float(KS3[I,J])*FA[J]
        S3 = S3 + ( SS3[I,0]*sin(A) + SS3[I,1]*cos(A) )

    A = 0.0
    for J in range(0,8):
        A = A + float(KS4[J])*FA[J]
    S4 = S4 + ( SS4[0]*sin(A) + SS4[1]*cos(A) )

    return ( S0 + ( S1 + ( S2 + ( S3 + ( S4 + 
                  S5 * T ) * T ) * T ) * T ) * T ) * DAS2R - X*Y/2.0

    # *  Finished.

'''
#==============================================================================
# Download/update EOP, meteo and iono data
#==============================================================================
'''
def doup(do_trp_calc, do_ion_calc, cat_eop, meteo_cat, ion_cat, \
             date_start, date_stop, iono_model='igs'):
    ''' Download/update EOP, meteo and iono data '''
    ''' EOPs '''
    eop_update(cat_eop, 3)

    ''' Don't do anything if the dates are in the future '''
    if (datetime.datetime.now() - date_start).total_seconds() < 0:
        print 'Experiment start date is in the future. Thus no tropo/iono data.'
        return
    
    ''' load data from first day to last+1 '''
    # set dates
    day_start = datetime.datetime(date_start.year,date_start.month,\
                                  date_start.day,0,0,0)
    day_stop = datetime.datetime(date_stop.year,date_stop.month,\
                      date_stop.day,0,0,0) + datetime.timedelta(days=1)
    dd = (day_stop - day_start).days

    # make list with datetime objects ranging from 1st day to last+1
    dates = [day_start]
    for d in range(1,dd+1):
        dates.append(day_start+ datetime.timedelta(days=d))

    ''' meteo data '''
    if do_trp_calc:
        # for each day
        for day in dates:
            year = day.year
            doy = day.timetuple().tm_yday # day of year
            
            ''' vmf1 files with precomputed site-specific data: '''
            # check file existance and download if necessary:
            vmf_file = '{:4d}{:03d}.vmf1_r'.format(year, doy)
            if not os.path.isfile(meteo_cat+'/'+vmf_file):
                try:
                    print 'vmf1 meteo file '+vmf_file+\
                          ' not found, downloading...'
                    met_url = 'http://ggosatm.hg.tuwien.ac.at/DELAY/SITE/VLBI/'+\
                              str(year)+'/'+vmf_file
                    response = urllib2.urlopen(met_url)
                    met_file = response.read()
                    # print it to file:
                    with open(meteo_cat+'/'+vmf_file,'w') as out:
                        for line in met_file:
                            out.write(line)
                except Exception, err:
                    print str(err)
                    print 'no troposphere available this time.'
                    
            
            ''' tropospheric gradient files: '''
            # check file existance and download if necessary:
            lhg_file = '{:4d}{:03d}.lhg_r'.format(year, doy)
            if not os.path.isfile(meteo_cat+'/'+lhg_file):
                try:
                    print 'tropo gradient file '+lhg_file+\
                        ' not found, downloading...'
                    met_url = 'http://ggosatm.hg.tuwien.ac.at/DELAY/ETC/LHG/VLBI/'+\
                                str(year)+'/'+lhg_file
                    response = urllib2.urlopen(met_url)
                    met_file = response.read()
                    # print it to file:
                    with open(meteo_cat+'/'+lhg_file,'w') as out:
                        for line in met_file:
                            out.write(line)
                except Exception, err:
                    print str(err)
                    print 'no tropo gradients available this time.'
                        
            ''' vmf1 grid files: '''
            vmf_grid_file = 'VMFG_{:4d}{:02d}{:02d}.H'\
                                .format(year, day.month, day.day)
            #(check the last file):
            if not os.path.isfile(meteo_cat+'/'+vmf_grid_file+'18'): 
                print 'vmf1 grid files '+vmf_grid_file+\
                      '00-18 not found, downloading...'
                for hh in ('00','06','12','18'):
                    try:
                        met_url = 'http://ggosatm.hg.tuwien.ac.at/DELAY/GRID/VMFG/'+\
                                  str(year) + '/' + vmf_grid_file + hh
                        response = urllib2.urlopen(met_url)
                        met_file = response.read()
                        # print it to file:
                        with open(meteo_cat+'/'+vmf_grid_file+hh,'w') as out:
                            for line in met_file:
                                out.write(line)
                    except Exception, err:
                        print str(err)
                        print 'TU Wien server is down or does not have needed'+\
                              ' files. Trying Canadian stuff..'
                        try:
                            met_url = 'http://unb-vmf1.gge.unb.ca/pub/unbvmfG/'+\
                                      str(year) + '/UNB' + vmf_grid_file + hh
                            response = urllib2.urlopen(met_url)
                            met_file = response.read()
                            # print it to file:
                            with open(meteo_cat+'/'+vmf_grid_file+hh,'w') as out:
                                for line in met_file:
                                    out.write(line)
                        except Exception, err:
                            print str(err)
                            print 'no vmf1 grid data this time'
            
    ''' ionospheric data '''
    if do_ion_calc:
        # for each day
        for day in dates:
            year = day.year
            doy = day.timetuple().tm_yday # day of year
            # download and uncompress the ionex files, if necessary
            ionex_zip = iono_model+'g{:03d}'.format(doy)+'0.'+str(year)[2:]+'i.Z'
            if not os.path.isfile(ion_cat+'/'+ionex_zip[:-2]):
                print 'iono TEC file: '+ionex_zip+' not found, fetching...'
                try:
                    ftp = FTP('cddis.nasa.gov')
                    ftp.login() # user anonymous, passwd anonymous
                    ftp.cwd('pub/gps/products/ionex/{:4d}/{:03d}'.format(year,doy))
                    if ionex_zip in ftp.nlst(): # check that the file is there
                        ftp.retrbinary('RETR {:s}'.format(ionex_zip), \
                                   open(ion_cat+'/'+ionex_zip, 'wb').write)
                        # uncompress:
                        print 'uncompressing: ' + ionex_zip + '...'
                        os.system('uncompress -f {:s}'.format(ion_cat+'/'+ionex_zip))
                    else:
                        print 'file {:s} not found on the server. no iono this time'.\
                                format(ionex_zip)
                    ftp.quit()
                    
                except Exception, err:
                    print str(err)
                    print 'no ionosphere this time'

'''
#==============================================================================
# Update eops, if they are out of date (older, than n days)
#==============================================================================
'''
def eop_update(cat_eop, n=0):
    do_update = False
    if os.path.isfile(cat_eop):
        age = datetime.datetime.now() - \
              datetime.datetime.utcfromtimestamp(os.path.getmtime(cat_eop))
        if age.days > n:
            do_update = True
            print 'EOP file: '+cat_eop+' is out of date, updating...'
    else:
        do_update = True
        print 'EOP file: '+cat_eop+' is missing, fetching...'
    # if the file is older than n days:
    if do_update:
        try:
            ftp = FTP('hpiers.obspm.fr')
            ftp.login() # user anonymous, passwd anonymous
            ftp.cwd('iers/series/opa')
            ftp.retrbinary('RETR eopc04_IAU2000', open(cat_eop, 'wb').write)
            ftp.quit()
        except Exception, err:
            print str(err)
            pass
        
        
#==============================================================================
# >>>>>>>>>  VLBI In the Near-field Toolkit -- VINT  <<<<<<<<<<<
#==============================================================================
#@jit(locals=dict(mydouble=double))
def vint_s(ob):
    '''
    Main engine used to calculate dudes (Delays, Uvws, Doppler, Etc.)
    
    v0.1 24/10/2013
    v0.2 23/11/2013
    v0.3 19/12/2013
    
    coded by Dmitry Duev (JIVE),
    as everything else here, pretty much
    '''
    # exp_name is used to create a folder in the _out directory. if empty,
    # the output goes straight into the _out dir.
    
    # tell the user what you're doing
    print ob
    
    ''' load input sittings: '''
#    inp = inp_set('inp.cfg')
    inp = ob.inp
    
    ''' init a dude object for output '''    
    dud = dude()
    
    # Deal with 'Zero-baseline' in delay mode: produce zeros
    # This is used when a station-centric mode of SFXC is wanted
    # in other words, when the phase center is at one of the stations
    if inp['delay_calc'] and ob.sta[0]==ob.sta[1]:
        if ob.sta[0]!='RA':
            dud.delay = np.zeros((len(ob.tstamps), 5))
        else:
            # in case some freak wants a RadioAstron-centric delays
            dud.delay = np.zeros((len(ob.tstamps), 7))
        if inp['uvw_calc']:
            dud.uvw = np.zeros((len(ob.tstamps), 3))
        # exit the function
        return dud
    
    ''' mkdir '_out/exp_name' if non existend '''
    exp_name = ob.exp_name

    if not os.path.isdir(os.path.join(inp['out_path'], exp_name)) and \
            (exp_name is not None or len(exp_name)>0):
        os.makedirs(os.path.join(inp['out_path'], exp_name))
    
    ''' init constants '''
    const = constants(inp['jpl_eph'])
    
    ''' load s/c ephemerides '''
    if ob.sou_type!='C':
        eph = load_sc_eph(ob.sou_type, ob.source, \
                          ob.tstamps[0], ob.tstamps[-1], inp, \
                          inp['uvw_calc'], \
                          inp['sc_rhophitheta'], \
                          inp['sc_xyz'])
    elif ('RA' in ob.sta): # RA was observing
        eph = load_sc_eph(ob.sou_type, 'RA', \
                          ob.tstamps[0], ob.tstamps[-1], inp, \
                          inp['uvw_calc'], \
                          inp['sc_rhophitheta'], \
                          inp['sc_xyz'])
    else:
        eph = None
    
    ''' load S/C freq ramping parameters if 2(3)-way Doppler '''
    if ob.sou_type!='C' and inp['doppler_calc'] and \
            inp['dop_model'] == 'bary3way':
        freq_ramp = np.array(freqRamp(cat_dir=inp['f_ramp'], \
                                       sc=ob.source.lower(), tx_type='3way'))
        freq_type = 'proper'
        # cut the time range of ob.tstamps:
        obs_start = Time(str(ob.tstamps[0]), format='iso', scale='utc')
        obs_stop = Time(str(ob.tstamps[-1]), format='iso', scale='utc')
        _, _, lt_sec = eph.RaDec_bc_sec(obs_start.tdb.jd1, \
                                    obs_start.tdb.jd2*86400.0, inp['jpl_eph'])
        numlt = datetime.timedelta(seconds=int(3.5*lt_sec))
#        print numlt
        obs_ind_start = np.logical_and(freq_ramp[:,0]<=obs_start.datetime-numlt,\
                                 freq_ramp[:,1]>=obs_start.datetime-numlt)
        ind_start = np.argmax(obs_ind_start)
        obs_ind_stop = np.logical_and(freq_ramp[:,0]<=obs_stop.datetime,\
                                     freq_ramp[:,1]>=obs_stop.datetime)
        ind_stop = np.argmax(obs_ind_stop)
        freq_ramp = freq_ramp[ind_start:ind_stop+1,:]
#        print freq_ramp
        # add uplink stations to self.sta list:
        ob.sta.append(list(set(freq_ramp[:,4])))
        # flatten the list
        ob.sta = list(flatten(ob.sta))
    
    ''' load S/C freq if 1-way Doppler '''
    if ob.sou_type!='C' and inp['doppler_calc'] and \
            inp['dop_model'] == 'bary1way':
        try:
            # is it a 'ramped' s/c? let's find out:
            freq_ramp = np.array(freqRamp(cat_dir=inp['f_ramp1w'], \
                                       sc=ob.source.lower(), tx_type='1way'))
            freq_type = 'proper'
        except Exception, err:
            print str(err)
            print 'No ramp table found for 1-way Tx of {:s}, trying sc.freq...'\
                    .format(ob.source.lower())
            try:
                freq, freq_type = freqConst(cat_file=inp['f_gc'], \
                                            sc=ob.source.lower())
                freq_ramp = None
            except:
                raise Exception('Can\'t find s/c Tx freq in cats/sc.freq')
        

    ''' read all kinds of catalogs '''
    # update/(down)load eops, meteo and iono data
    if internet_on():
        try:
            doup(inp['do_trp_calc'], inp['do_ion_calc'], \
                 inp['cat_eop'], inp['meteo_cat'], inp['ion_cat'], \
                 ob.tstamps[0], ob.tstamps[-1], inp['iono_model'])
        except Exception, err:
            print str(err)
            print 'catalogue update failed'
    else:
#        print 'no internet connection. can\'t update catalogues'
        pass
    # update eop file, if it's out of date, namely, older than (n) days:
#    eop_update(inp['cat_eop'], 3)
    
    # if RA was observing, add Pu or Gt to the station list, it'll be used for
    # calculating formatter time offset for the RA downlink stream
#    if 'RA' in ob.sta:
#        ob.sta.append('PUSHCHIN')
    # now do read everything:
    sou, sta, eops = load_cats(inp, ob.source, ob.sou_type, \
                               ob.sta, ob.tstamps[0], ob.sou_radec)
    ''' calculate site positions in geodetic coordinate frame
    + transformation matrix VW from VEN to the Earth-fixed coordinate frame '''
#    sta = geodetic(sta, const)
    for st in sta:
        st.geodetic(const)
    
    ''' calculate [Az, El] UTC series for the S/C case '''
    # lt-corrected ITRF ephemeris of the S/C is used here for each station:
#    if ob.sou_type!='C' or ('RA' in ob.sta):
#        t_0 = (ob.tstamps[0].hour + ob.tstamps[0].minute/60.0 + \
#                ob.tstamps[0].second/3600.0)/24.0
#        for st in sta:
#            t_obs = [(x - ob.tstamps[0]).total_seconds()/86400.0 + t_0 \
#                     for x in ob.tstamps]
#            st.AzEl(eph.gtrs, JD, eph.UT, t_obs, inp['jpl_eph'], interpolants=True)
        
    ''' Load MF & Meteo data into the site-objects '''
    if inp['do_trp_calc']:
        for st in sta:
            try:
                st.addMet(ob.tstamps[0], ob.tstamps[-1], inp)
            except Exception, err:
                print str(err)
                print 'No troposphere for ' + st.name + \
                      ' this time. Sorry about that..'
                # check met-params len when calculating tropo delay later!
                continue

    ''' load ionospheric TEC data into an ion object '''
    if inp['do_ion_calc']:
        try:
            iono = ion(ob.tstamps[0], ob.tstamps[-1], inp)
        except Exception, err:
            print str(err)
            print 'Unable to load vTEC data. Ionospheric delay will not be computed.'
            inp['do_ion_calc'] = False
    
    
    ''' ########  main loop ######## '''
    mjd_start = mjuliandate(ob.tstamps[0].year,\
                               ob.tstamps[0].month,ob.tstamps[0].day)
    
#    from astropy.utils.iers import IERS_A, IERS_A_URL
#    from astropy.utils.iers import IERS_B, IERS_B_URL
#    from astropy.utils.data import download_file
#    iers_b_file = download_file(IERS_B_URL, cache=True)  
#    iers_b = IERS_B.open(iers_b_file)                     

    uv_const = const.C  / (inp['mas_step']*(pi/648000000.0))

#    tic = _time()
    for ti, tstamp in enumerate(ob.tstamps):
#        tic = _time()
#        print tstamp
        ''' obs freq if not Doppler '''
        if not inp['doppler_calc']:
#            print ob.freqs
            freq = [f[2] for f in ob.freqs if f[0]<=tstamp<=f[1]][0]
    
        ''' set dates: '''
        mjd = mjuliandate(tstamp.year, tstamp.month, tstamp.day)
        dd = mjd - mjd_start
        UTC = (tstamp.hour + tstamp.minute/60.0 + \
                (tstamp.second + tstamp.microsecond*1e-6)/3600.0)/24.0
        JD = mjd + 2400000.5
        
#        toc = _time()
#        print toc-tic
#        tic = _time()

        ''' compute tai & tt '''
        TAI, TT = taitime(mjd, UTC)
#        print 'TT = {:.18f}'.format(TT)
#        toc = _time()
#        print toc-tic
#        tic = _time()

        
        ''' interpolate eops to tstamp '''
        UT1, eop_int = eop_iers(mjd, UTC, eops)

#        toc = _time()
#        print toc-tic
#        tic = _time()


        ''' compute coordinate time fraction of CT day at 1st observing site '''
        if sta[0].name == 'RA':
#            ra_gtrs = [eph.fGtrs[0](UTC+dd),\
#                       eph.fGtrs[1](UTC+dd),\
#                       eph.fGtrs[2](UTC+dd)]
            ra_gtrs = [np.interp(UTC+dd, eph.UT, eph.gtrs[:,6]),\
                       np.interp(UTC+dd, eph.UT, eph.gtrs[:,7]),\
                       np.interp(UTC+dd, eph.UT, eph.gtrs[:,8])]
            lon_gcen = atan2(ra_gtrs[1], ra_gtrs[0])
            if (lon_gcen < 0.0):
                lon_gcen += 2.0*pi
            sta[0].lon_gcen = lon_gcen
            sta[0].u = sqrt(ra_gtrs[0]**2 + ra_gtrs[1]**2 ) *1e-3
            sta[0].v = ra_gtrs[2]*1e-3
        # subj
        CT, dTAIdCT = t_eph(JD, UT1, TT, sta[0].lon_gcen, sta[0].u, sta[0].v)
#        tt_tdb = dtdb(JD, TT, UT1, sta[0].lon_gcen, sta[0].u, sta[0].v)
#        print 'dtdb = {:.18f}'.format(TT+tt_tdb/86400)
#        print 'UT1 = {:.18f}'.format(UT1)
#        print 'CT_mine = {:.18f}'.format(CT)
        
        astro_tstamp = Time(str(tstamp), format='iso', scale='utc', precision=9,
                 location=EarthLocation.from_geocentric(*sta[0].r_GTRS, unit=units.m))
#        t = Time(str(tstamp), format='iso', scale='utc', precision=9, \
#                 location=EarthLocation.from_geocentric(*sta[0].r_GTRS, unit=units.m))
#        t.delta_ut1_utc = t.get_delta_ut1_utc(iers_b)
#        print 'delta_ut1_utc = {:.18f}'.format(float(t.delta_ut1_utc))
#        print t._get_delta_tdb_tt
#        print t.tdb.iso
#        print t.tdb.jd
#        CT = (t.tdb.jd-0.5) - np.floor(t.tdb.jd-0.5)
#        print 'CT = {:.18f}'.format(CT)
        
#        toc = _time()
#        print toc-tic
#        tic = _time()

        ''' cut s/c ephem for faster computations '''
        if ob.sou_type!='C' or 'RA' in ob.sta:
            # create ephem-object to keep cut ephem for faster computations:
            eph_cut = ephem(eph.sc_name)
#            numlt = eph.fLt_gc(UTC+dd)/86400.0
            _, _, numlt = eph.RaDec_bc_sec(JD, CT*86400.0, inp['jpl_eph'])
            numlt /= 86400.0
            # should cover uplink sta for 3-way Doppler (2lt back) - another lt
            numlt = max(15.0/86400.0, 3.5*numlt) 
            # bcrs
            mini = max(np.searchsorted(eph.CT, CT+dd-numlt)-1, 0)
            maxi = np.searchsorted(eph.CT, CT+dd+numlt)
#            print mini, maxi
            eph_cut.CT = eph.CT[mini:maxi]
            eph_cut.CT_sec = eph.CT_sec[mini:maxi]
            for i, v in enumerate(eph.bcrs):
                eph_cut.bcrs.append(v[mini:maxi,:])
#            print len(eph_cut.bcrs)
            # gcrs
            mini = max(np.searchsorted(eph.UT, UTC+dd-numlt)-1, 0)
            maxi = np.searchsorted(eph.UT, UTC+dd+numlt)
            eph_cut.UT = eph.UT[mini:maxi]
            eph_cut.gcrs = eph.gcrs[mini:maxi,:]
            eph_cut.gtrs = eph.gtrs[mini:maxi,:]

        ''' BCRS state vectors of celestial bodies at JD+CT, [m, m/s]: '''
        ## Earth:
#        print JD
        rrd = pleph(JD+CT, 3, 12, inp['jpl_eph'])
        earth = np.reshape(np.asarray(rrd), (3, 2), 'F') * 1e3
        # Earth's acceleration in m/s**2:
        v_plus = np.array(pleph(JD+CT+1.0/86400.0, 3, 12, inp['jpl_eph'])[3:])
        v_minus = np.array(pleph(JD+CT-1.0/86400.0, 3, 12, inp['jpl_eph'])[3:])
        a = (v_plus - v_minus)*1e3 / 2.0
        a = np.array(np.matrix(a).T)
        earth = np.hstack((earth, a))
#        print earth[:,0]
        ## Sun:
        rrd = pleph(JD+CT, 11, 12, inp['jpl_eph'])
        sun = np.reshape(np.asarray(rrd), (3,2), 'F') * 1e3
#        print sun[:,0]
        ## Moon:
        rrd = pleph(JD+CT, 10, 12, inp['jpl_eph'])
        moon = np.reshape(np.asarray(rrd), (3,2), 'F') * 1e3
        ## Planets:
#        planets = []
#        for jj in (1,2,4,5,6,7,8):
#            rrd = pleph(JD+CT, jj, 12, inp['jpl_eph'])
#            planets.append(np.reshape(np.asarray(rrd), (3,2), 'F') * 1e3)
        ## All SSBC vectors in one list for near-field:
        if ob.sou_type!='C':
            state_ss = []
            for jj in (1,2,4,5,6,7,8,9):
                rrd = pleph(JD+CT, jj, 12, inp['jpl_eph'])
                state_ss.append(np.reshape(np.asarray(rrd), (3,2), 'F') * 1e3)
            state_ss.insert(2,earth)
            state_ss.append(moon)
            state_ss.append(sun)
#        toc = _time()
#        print toc-tic
#        tic = _time()

        ''' rotation matrix IERS '''
        r2000 = ter2cel(tstamp, eop_int, dTAIdCT, 'iau2000')
#        print r2000[:,:,0]
#        print r2000[:,:,1]
#        toc = _time()
#        print toc-tic
#        tic = _time()

        ''' displacements due to geophysical effects '''
        for st in sta:
            if st.name == 'GEOCENTR' or st.name == 'RA':
                continue
            # displacement due to solid Earth tides:
#            tic = _time()
            st = dehanttideinel(st, tstamp, earth, sun, moon, r2000)
#            print _time()-tic
            # displacement due to ocean loading:
#            tic = _time()
            st = hardisp(st, tstamp, r2000)
#            print _time()-tic
            # rotational deformation due to pole tide:
#            tic = _time()
            st = poletide(st, tstamp, eop_int, r2000)
#            print _time()-tic
            # displacement due to atmospheric loading:
            # Should be the one by Ray and Ponte (2003)
            # http://geophy.uni.lu/ggfc-atmosphere/tide-loading-calculator.html
        
#        toc = _time()
#        print toc-tic
#        tic = _time()

        
        ''' add up geophysical corrections and convert sta state to J2000 '''
        for st in sta:
#            if 'eph' in locals():
            if eph != None:
                st.j2000gp(r2000, eph_cut.gcrs, eph_cut.UT, UTC+dd)
            else:
                st.j2000gp(r2000)
#            print st.r_GCRS
#            print st.v_GCRS
#            print st.a_GCRS

#        toc = _time()
#        print toc-tic
#        tic = _time()

        ''' uplink sta and f for 2(3)-way Doppler '''
        if ob.sou_type!='C' and inp['doppler_calc'] and \
            inp['dop_model']=='bary3way':
            _, _, lt_sec = eph.RaDec_bc_sec(JD, CT*86400.0, inp['jpl_eph'])
#            lt = datetime.timedelta(seconds=int(eph.fLt_gc(UTC+dd))) # lt [sec]
            lt = datetime.timedelta(seconds=int(lt_sec)) # lt [sec]
#            print lt, tstamp-2*lt
            try:
                # tstart, tstop in UTC; f_0, df, up_sta
#                rampTslot = [x for x in freq_ramp if \
#                                x[0] <= tstamp-2*lt <= x[1]][-1]
                rampTslot = freq_ramp[np.logical_and(freq_ramp[:,0]<=tstamp-2*lt,\
                                        freq_ramp[:,1]>=tstamp-2*lt),:][0]
#                print lt, rampTslot
                upSta = rampTslot[4]
            except Exception, err:
                print str(err)
                raise Exception('Not found ramp params for '+\
                                ob.source+' at ' + str(tstamp))
                                
        ''' coarse! f ramp table for 1-way deep space Doppler (for iono) '''
        if ob.sou_type!='C' and inp['doppler_calc'] and \
            inp['dop_model']=='bary1way' and freq_ramp is not None:
            _, _, lt_sec = eph.RaDec_bc_sec(JD, CT*86400.0, inp['jpl_eph'])
#            lt = datetime.timedelta(seconds=int(eph.fLt_gc(UTC+dd))) # lt [sec]
            lt = datetime.timedelta(seconds=int(lt_sec)) # lt [sec]
            try:
                # tstart, tstop in TDB; f_0, df
#                print astro_tstamp.tdb.datetime-lt
#                tic = _time()
#                rampTslot = [x for x in freq_ramp if \
#                             x[0] <= astro_tstamp.tdb.datetime-lt <= x[1]][-1]
#                print rampTslot, _time()-tic
#                tic = _time()
#                print np.logical_and(\
#                                freq_ramp[:,0]<astro_tstamp.tdb.datetime-lt,\
#                                freq_ramp[:,1]>=astro_tstamp.tdb.datetime-lt)
                rampTslot = freq_ramp[np.logical_and(\
                                freq_ramp[:,0]<=astro_tstamp.tdb.datetime-lt,\
                                freq_ramp[:,1]>=astro_tstamp.tdb.datetime-lt),:][0]
#                print rampTslot#, _time()-tic
            except Exception, err:
                print str(err)
                raise Exception('Not found ramp params for '+\
                                ob.source+' at ' + str(tstamp))
        
        ''' delays due to instrumental and propagation effects '''
        for si, st in enumerate(sta):
            if st.name == 'GEOCENTR' or st.name == 'RA':
                continue
            ## Source's [current] [el,az] and [ra,dec] for current station:
            # s/c observation
            if ob.sou_type!='C':
                # for 3-way Doppler: az, el 2lt's ago (transmitter):
                if inp['doppler_calc'] and inp['dop_model']=='bary3way' and \
                                si>0:
                    continue # this is all done once precise transmission t is found

                else:
                    # 'normal' case ('receiving station')
                    az, el, _ = st.AzEl2(eph_cut.gcrs, eph_cut.UT, JD, UTC,\
                                     r2000[:,:,0], inp['jpl_eph'])
#                    print tstamp, st.name, az*180/np.pi, el*180/np.pi
                    ra, dec, _ = eph_cut.RaDec_bc_sec(JD, CT*86400.0, inp['jpl_eph'])
#                    print ra*180/np.pi, dec*180/np.pi
#                    radec = np.hstack((Angle(ra, unit=u.rad).hms, \
#                                        Angle(dec, unit=u.rad).dms))
#                    print 'ra = {:02.0f}h{:02.0f}m{:010.7f}s  dec = {:-3.0f}d{:02.0f}\'{:010.7f}\"\n'\
#                            .format(*radec)
#                    raw_input()
            # recording station for RA downlink:
            elif ob.sou_type=='C' and ('RA' in ob.sta) \
                    and ((st.name == 'PUSHCHIN') or (st.name == 'GREENBNK')):
#                az = st.azel_interp[0](UTC+dd)
#                el = st.azel_interp[1](UTC+dd)
#                az, el = st.AzEl(eph_cut.gtrs, JD, eph_cut.UT, UTC+dd, inp['jpl_eph'])
                az, el, _ = st.AzEl2(eph_cut.gcrs, eph_cut.UT, JD, UTC,\
                                     r2000[:,:,0], inp['jpl_eph'])
                ra = sou.ra
                dec = sou.dec
            # far-field source:
            else:
                az, el = aber_source(st.v_GCRS, st.vw, sou.K_s, r2000, earth)
                ra = sou.ra
                dec = sou.dec
                
            # current met params:
            if len(st.met['TC'])>0:
                T_site = st.fMet['fT'](mjd+UTC)
                P_site = st.fMet['fP'](mjd+UTC)
                H_site = st.fMet['fH'](mjd+UTC)
            else:
                T_site = 0.0
                P_site = 1000.0
                H_site = 100.0
                
            ## thermal deformation of telescopes:
            st = thermal_def(st, ra, dec, el, T_site, const.C)
            ## axis offset:
            st = mount_tel(st, r2000[:,:,0], el, az, T_site, P_site, H_site, const)
            ## tropospheric delay:
            if inp['do_trp_calc'] and len(st.met['ahz'])>0:
                st = tropo_wien(st, el, az, mjd+UTC, const, inp['do_trp_grad_calc'])
            ## ionospheric delay:
            if inp['do_ion_calc'] and len(iono.fVTEC)>0:
                if ob.sou_type=='C' and ((st.name == 'PUSHCHIN') \
                                      or (st.name == 'GREENBNK')):
                    f_0 = 15e9 # RadioAstron data stream downlink freq
                elif inp['doppler_calc'] and inp['dop_model']=='bary3way' and \
                       si>0 and st.name==upSta and len(iono.fVTEC)>0:
                    continue # this is done once precise transmission t is found
#                    f_0 = rampTslot[2] # use uplink freq for 3-way Doppler, upsta
                elif inp['doppler_calc'] and inp['dop_model']=='bary3way' and \
                       si==0 and len(iono.fVTEC)>0:
                    f_0 = rampTslot[2]*inp['tr'] # downlink freq for 3-Dop, downsta
                elif inp['doppler_calc'] and inp['dop_model']=='bary1way' and \
                       freq_ramp is not None and len(iono.fVTEC)>0:
                    f_0 = rampTslot[2] # downlink freq for ramped 1-Dop
                else:
                    f_0 = freq

                st = ion_igs(st, iono, el, az, JD, UTC, f_0)
        
#        toc = _time()
#        print toc-tic
#        tic = _time()
        
        ''' geometric delay '''
        if inp['delay_calc']:
            ## Concensus model IERS Conventions 2010:
            if ob.sou_type=='C' and ('RA' not in ob.sta):
                dtau = delay_iers(JD, CT, sta[0].r_GCRS, sta[1].r_GCRS, \
                                  sta[1].v_GCRS, earth, sun, sou.K_s, inp['jpl_eph'],\
                                  const.GM, const.TDB_TCB, const.L_C, const.C)
            # this correction is rather small, thus neglected: [see tn36 11.10]
    #        print sta[1].dtau_tropo*dot(sou.K_s,sta[1].v_GCRS-sta[0].v_GCRS)/const.C
            
            ## Vlasov/Zharov/Sazhin far-field for RadioAdtron:
            if ob.sou_type=='C' and ('RA' in ob.sta):
                dtau, dtau_dt_min1, lt_downlink = delay_ra(JD, CT, UTC, \
                                  sta[0].r_GCRS, sta[1].r_GCRS, \
                                  sta[1].v_GCRS, sta[1].a_GCRS, sta[2].r_GCRS, \
                                  earth, sun, \
                                  sou.K_s, inp['jpl_eph'], const.GM, const.TDB_TCB, \
                                  const.L_C, const.C, const.AE)
            
            ## Near-field Fukushima:
            if ob.sou_type!='C' and inp['nf_model']=='Fukushima':
                dtau = delay_nf_fukushima(JD, CT, dd, \
                                    state_ss, eph_cut.CT_sec, eph_cut.bcrs[0], \
                                    const, sta[0], sta[1])
#            dtau_ = dtau
#            print 'dtau_ = {:.18f}'.format(dtau)

            # Near-field Moyer/Duev aka LTea42:
            if ob.sou_type!='C' and inp['nf_model']=='Moyer':
#                dtau = delay_moyer(JD, CT*86400.0, dd*86400.0, sta[0].r_GCRS, \
#                                  sta[1].r_GCRS, sta[1].v_GCRS, sta[1].a_GCRS,\
#                                  state_ss, eph_cut.CT_sec, eph_cut.bcrs[0], \
#                                  const.GM, const.TDB_TCB, \
#                                  const.L_C, const.C, inp)
#                print 'dtau__ = {:.18f}'.format(dtau)
#                dtau__ = dtau
                dtau = delay_nf_moyer(JD, CT, dd, \
                                      state_ss, eph_cut.CT_sec, eph_cut.bcrs[0], \
                                      const, sta[0], sta[1], inp, UTC)
#                print 'dtau = {:.18f}'.format(dtau)
#                print tstamp, dtau_ - dtau, dtau_ - dtau__
#            print dtau_, dtau
#            raw_input()
            
            ''' stack delays due to all effects together: '''
            # RadioAstron far-field:
            # [(geometry, geophysics); instrumental(thermal, axis offset),
            #  propagation(troposphere, ionosphere),
            #  dtau/dt - 1, LT to Ra for tracking station]
            # instrumental effects are given for receiving station
            if ob.sou_type=='C' and ('RA' in ob.sta):
                if ob.sta[1]=='RA':
                    delay = [dtau, \
                         (sta[2].dtau_therm - sta[0].dtau_therm),\
                         (sta[2].dtau_ao - sta[0].dtau_ao),\
                         (sta[2].dtau_tropo - sta[0].dtau_tropo),\
                         (sta[2].dtau_iono - sta[0].dtau_iono),\
                         dtau_dt_min1, lt_downlink]
                # perverted case of 'Ra-centric'-delays
                if ob.sta[0]=='RA':
                    # indices of sta[xx] are not set at random,
                    # I actually did some thinking about the issue!
                    delay = [dtau, \
                         (sta[1].dtau_therm - sta[2].dtau_therm),\
                         (sta[1].dtau_ao - sta[2].dtau_ao),\
                         (sta[1].dtau_tropo - sta[2].dtau_tropo),\
                         (sta[1].dtau_iono - sta[2].dtau_iono),\
                         dtau_dt_min1, lt_downlink]

            else:
            # usual far-field or near-field:
            # [(geometry, geophysics), instrumental(thermal, axis offset),
            #  propagation(troposphere, ionosphere)]
                delay = [dtau,\
                         (sta[1].dtau_therm - sta[0].dtau_therm),\
                         (sta[1].dtau_ao - sta[0].dtau_ao),\
                         (sta[1].dtau_tropo - sta[0].dtau_tropo),\
                         (sta[1].dtau_iono - sta[0].dtau_iono)]
            
                     
            # the rug really tied the room together, did it not?:
            dud.delay.append(delay)
        
#        toc = _time()
#        print toc-tic
#        tic = _time()

        
        ''' far-field uvw '''
        # compute numerically as c*(cos(dec)*dtau/dra,dtau/ddec,tau)
        if inp['uvw_calc'] and ob.sou_type=='C':
            # u, v
            K_plus_ra = sph2cart([1.0, dec, ra + 1e-5]) # u
            K_plus_dec = sph2cart([1.0, dec + 1e-4, ra]) # v
            if 'RA' not in ob.sta:
                dtau_plus_ra = delay_iers(JD, CT, sta[0].r_GCRS, sta[1].r_GCRS, \
                              sta[1].v_GCRS, earth, sun, K_plus_ra, inp['jpl_eph'],\
                              const.GM, const.TDB_TCB, const.L_C, const.C)
                dtau_plus_dec = delay_iers(JD, CT, sta[0].r_GCRS, sta[1].r_GCRS, \
                              sta[1].v_GCRS, earth, sun, K_plus_dec, inp['jpl_eph'],\
                              const.GM, const.TDB_TCB, const.L_C, const.C)
            else:
                dtau_plus_ra  = delay_ra(JD, CT, UTC, \
                              sta[0].r_GCRS, sta[1].r_GCRS, \
                              sta[1].v_GCRS, sta[1].a_GCRS, sta[2].r_GCRS, \
                              earth, sun, \
                              K_plus_ra, inp['jpl_eph'], const.GM, const.TDB_TCB, \
                              const.L_C, const.C, const.AE, uv=True)
                dtau_plus_dec  = delay_ra(JD, CT, UTC, \
                              sta[0].r_GCRS, sta[1].r_GCRS, \
                              sta[1].v_GCRS, sta[1].a_GCRS, sta[2].r_GCRS, \
                              earth, sun, \
                              K_plus_dec, inp['jpl_eph'], const.GM, const.TDB_TCB, \
                              const.L_C, const.C, const.AE, uv=True)
            
            u = const.C * (dtau_plus_ra - dtau)/(1e-5) / cos(dec)
            v = const.C * (dtau_plus_dec - dtau)/(1e-4)
            # w
            w = dtau*const.C
            
            # Donny, please:
            dud.uvw.append([u, v, w])
            
        ''' near-field uvw '''
        # compute numerically
#        tic1=_time()
        if inp['uvw_calc'] and inp['sc_rhophitheta'] and ob.sou_type!='C':
            dtaus = []
            for bcrs in eph_cut.bcrs[1:]:
#                dtaus.append(delay_moyer(JD, CT*86400.0, dd*86400.0, sta[0].r_GCRS, \
#                                  sta[1].r_GCRS, sta[1].v_GCRS, sta[1].a_GCRS,\
#                                  state_ss, eph_cut.CT_sec, bcrs, \
#                                  const.GM, const.TDB_TCB, \
#                                  const.L_C, const.C, inp))
                if inp['nf_model']=='Moyer':
                    dtaus.append(delay_nf_moyer(JD, CT, dd, \
                                    state_ss, eph_cut.CT_sec, bcrs, \
                                    const, sta[0], sta[1],\
                                    inp, UTC))
                elif inp['nf_model']=='Fukushima':
                    dtaus.append(delay_nf_fukushima(JD, CT, dd, \
                                    state_ss, eph_cut.CT_sec, bcrs, \
                                    const, sta[0], sta[1]))
#            dtaus.append(dtau)
#            dud.delay.append(dtaus)

            # this is u (c*dtau/dtheta):
            u = uv_const * (dtaus[2] - dtau) / cos(dec)
            # this is v (c*dtau/dphi):
            v = uv_const * (dtaus[1] - dtau)
            # w (dtau/drho)
            w = (dtaus[0] - dtau)/inp['m_step']
                        
#            u = uv_const * (dtaus[4] - dtaus[5]) / (2.0*cos(dec))
#            v = uv_const * (dtaus[2] - dtaus[3]) / 2.0
#            w = (dtaus[0] - dtaus[1])/(2.0*inp['m_step'])

            # Donny, please:
            dud.uvw.append([u, v, w])
#        print _time()-tic1

        ''' Doppler '''
        if inp['doppler_calc']:
            # 1-way geocentric model
            if inp['dop_model'] == 'geo1way':
                raise NotImplemented

            # 1-way barycentric model
            if inp['dop_model'] == 'bary1way':
                if freq_type=='gc': # gnss
                    fonafe, _ = doppler_bc(JD, CT*86400.0, dd*86400.0, \
                                        state_ss, eph_cut.CT_sec, eph_cut.bcrs[0], \
                                        const.GM, const.TDB_TCB, \
                                        const.L_C, const.C, 'one', freq_type, UT1, \
                                        sta[0], None, const.L_G, const.AE, const.J_2, \
                                        eph_cut.UT*86400.0, eph_cut.gcrs, UTC*86400.0,\
                                        inp=inp)
                else:
                    fonafe, t_tr = doppler_bc(JD, CT*86400.0, dd*86400.0, \
                                        state_ss, eph_cut.CT_sec, eph_cut.bcrs[0], \
                                        const.GM, const.TDB_TCB, \
                                        const.L_C, const.C, 'one', freq_type, UT1, \
                                        sta[0], inp=inp)
#                                    sta[0].r_GCRS, sta[0].v_GCRS)
                ''' stack all effects together: '''
                # no ramping:
                if freq_ramp is None:
                    doppler = [freq*fonafe, \
                               UTC+dd, \
                               sta[0].dtau_therm, sta[0].dtau_ao,\
                               sta[0].dtau_tropo, sta[0].dtau_iono]
                else:
                    # ramped frequency at the TDB moment of signal transmission:
#                    print rampTslot, t_tr.datetime, tstamp - t_tr.datetime
                    # reload rampTslot, just in case:
                    rampTslot = freq_ramp[np.logical_and(\
                                            freq_ramp[:,0]<t_tr.datetime,\
                                            freq_ramp[:,1]>=t_tr.datetime),:][0]
#                    rampTslot = [x for x in freq_ramp if \
#                                x[0] <= t_tr.datetime <= x[1]][-1]
#                    print rampTslot, t_tr.datetime, tstamp, '\n'
#                    raw_input()
                    freq = rampTslot[2] + rampTslot[3]*\
                               ((t_tr.datetime-rampTslot[0]).total_seconds())
                    doppler = [freq*fonafe, \
                               UTC+dd, \
                               sta[0].dtau_therm, sta[0].dtau_ao,\
                               sta[0].dtau_tropo, sta[0].dtau_iono]
                           
            # 3-way barycentric model               
            if inp['dop_model'] == 'bary3way':
                # current uplink station index in sta list:
                cus = [x.name for x in sta[1:]].index(upSta) + 1
                fonafe, t_tr = doppler_bc(JD, CT*86400.0, dd*86400.0, \
                                    state_ss, eph_cut.CT_sec, eph_cut.bcrs[0], \
                                    const.GM, const.TDB_TCB, \
                                    const.L_C, const.C, 'three', freq_type, \
                                    UT1, sta[0], sta[cus], \
                                    AE=const.AE, J_2=const.J_2, L_G=const.L_G, \
                                    eops=eops, inp=inp)
#                                    sta[0].r_GCRS, sta[0].v_GCRS, \
#                                    sta[cus].r_GCRS, sta[cus].v_GCRS, \
#                                    sta[cus].a_GCRS)
                #reload rampTslot, just in case:
#                print t_tr.datetime
                rampTslot = freq_ramp[np.logical_and(\
                                        freq_ramp[:,0]<t_tr.datetime,\
                                        freq_ramp[:,1]>=t_tr.datetime),:][0]
#                print rampTslot
                # ramped frequency at the moment of uplink signal transmission:
                freq = rampTslot[2] + rampTslot[3]*\
                               ((t_tr.datetime-rampTslot[0]).total_seconds())
#                print tstamp
#                print tstamp, 'f = {:.9f}\n'.format(freq*fonafe*inp['tr'])
#                raw_input()
                ''' calculate corrections for the uplink station: '''
#                print 'bylo: ', sta[cus].dtau_therm, sta[cus].dtau_ao, \
#                               sta[cus].dtau_tropo, sta[cus].dtau_iono
#                az, el = sta[cus].AzEl(eph_cut.gtrs, JD, eph_cut.UT, \
#                                     t_tr.jd2 + dd, inp['jpl_eph'])
#                print t_tr.datetime
                JD_UTC_tr = np.floor(t_tr.mjd) + 2400000.5
                UTC_tr = (t_tr.datetime.hour + t_tr.datetime.minute/60.0 + \
                            (t_tr.datetime.second + \
                                t_tr.datetime.microsecond*1e-6)/3600.0)/24.0
                lt_tr = (tstamp-t_tr.datetime).total_seconds()
#                print eph_cut.UT, JD, UTC+dd, JD_UTC_tr, UTC_tr+dd, lt_tr
                az, el, _ = sta[cus].AzEl2(eph_cut.gcrs, eph_cut.UT, JD_UTC_tr, UTC_tr,\
                                     r2000[:,:,0]-lt_tr*r2000[:,:,1], inp['jpl_eph'])
#                print tstamp, sta[cus].name, az*180/np.pi, el*180/np.pi
                JD_tr = np.floor(t_tr.tdb.mjd) + 2400000.5
                CT_tr = (t_tr.tdb.datetime.hour + t_tr.tdb.datetime.minute/60.0 + \
                            (t_tr.tdb.datetime.second + \
                                t_tr.tdb.datetime.microsecond*1e-6)/3600.0)/24.0
#                print JD, CT, JD_tr, CT_tr, dd
#                if CT_tr > CT: # if tr a day before the first day of obs:
#                    CT_tr -= 1
                                
#                print tstamp, JD, CT, JD_tr, CT_tr, dd
                ra, dec, _ = eph_cut.RaDec_bc_sec(JD_tr, CT_tr*86400.0, inp['jpl_eph'])
#                print ra, dec
                # current met params:
                if len(sta[cus].met['TC'])>0:
                    T_site = sta[cus].fMet['fT'](t_tr.mjd)
                    P_site = sta[cus].fMet['fP'](t_tr.mjd)
                    H_site = sta[cus].fMet['fH'](t_tr.mjd)
                else:
                    T_site = 0.0
                    P_site = 1000.0
                    H_site = 100.0
                
                # thermal deformation:
                sta[cus] = thermal_def(sta[cus], ra, dec, el, T_site, const.C)
                ## axis offset:
                sta[cus] = mount_tel(sta[cus], r2000[:,:,0]-lt_tr*r2000[:,:,1], \
                                     el, az, T_site, P_site, H_site, const)
                ## tropospheric delay:
                if inp['do_trp_calc'] and len(sta[cus].met['ahz'])>0:
                    sta[cus] = tropo_wien(sta[cus], el, az, t_tr.mjd, const,\
                                          inp['do_trp_grad_calc'])
                ## ionospheric delay:
                if inp['do_ion_calc'] and len(iono.fVTEC)>0:
                    f_0 = rampTslot[2] # use uplink freq for 3-way Doppler, upsta
                    sta[cus] = ion_igs(sta[cus], iono, el, az, JD_UTC_tr, UTC_tr, f_0)
                
#                print 'stalo: ', sta[cus].dtau_therm, sta[cus].dtau_ao, \
#                               sta[cus].dtau_tropo, sta[cus].dtau_iono
#                print tstamp, sta[0].dtau_tropo, sta[cus].dtau_tropo
#                print tstamp, sta[0].dtau_iono, sta[cus].dtau_iono
                ''' stack all effects together: '''
                doppler = [freq*fonafe*inp['tr'], \
                           UTC+dd, \
                           (sta[0].dtau_therm + sta[cus].dtau_therm),\
                           (sta[0].dtau_ao + sta[cus].dtau_ao),\
                           (sta[0].dtau_tropo + sta[cus].dtau_tropo),\
                           (sta[0].dtau_iono + sta[cus].dtau_iono)]
                
            ''' append to dude-obj '''
            dud.doppler.append(doppler)
#        toc = _time()
#        print toc-tic

#    toc = _time()
#    print toc-tic

    ''' output '''
    if inp['delay_calc']:
        dud.delay = np.array(dud.delay)
    if inp['uvw_calc']:
        dud.uvw = np.array(dud.uvw)
    if inp['doppler_calc']:
        # list -> np.array
        doppler_raw = np.array(dud.doppler)
#        ''' correction of geom prediction epoch due to instr/prop delay '''
#        t_corr = doppler_raw[:,1]*86400.0 - \
#                  np.sum(doppler_raw[:,2:], axis=1)
#        
#        doppler = np.interp(t_corr, doppler_raw[:,1]*86400.0, doppler_raw[:,0])
#    
#        doppler = np.vstack((doppler_raw[:,0], doppler - doppler_raw[:,0]))

        ''' turn delays due to instr. and propagation effects
            into delay rates in Hz '''
        tp = doppler_raw[:,1]*86400.0
        doppler = doppler_raw[:,0]
        
        for ii in range(2,6):
            # make cheb polyfit to delays, then differentiate it to get rates
#            cheb_order = min(len(doppler_raw[:,ii]), 10)
#            p = cheb.chebfit(tp, doppler_raw[:,ii], cheb_order)
#            drate = cheb.chebval(tp, cheb.chebder(p))
            # optimal high-order Cheb-fit:
#            p = optimalFit(tp, doppler_raw[:,ii], \
#                               min_order=10, max_order=15, fit_type='cheb')
#            drate = cheb.chebval(tp, cheb.chebder(p.best_estimator_.coef_))
            # local polyfit:
#            fig = plt.figure()
#            ax = fig.add_subplot(211)
#            ax.plot(tp, doppler_raw[:,ii], '.')
            drate = derivative(tp, doppler_raw[:,ii], points=7, poly=2)
#            ax = fig.add_subplot(212)
#            ax.plot(tp, drate, '.')
            doppler = np.vstack((doppler, -drate*doppler_raw[:,0]))
        dud.doppler = doppler.T
    return dud
    
    
#==============================================================================
#     
#==============================================================================
def load_cats(inp, sou_name, sou_type, sta_names, date_t_start, sou_radec=None):
    '''    
     dyear - difference btw epoch of observation and epoch of station catalog (EPOCH)
       Arrays:
           /ocean/
            1) amp_ocean(11,3,N) - ocean loading amplitudes for 11 tides in
               next order: M2  S2  N2  K2  K1  O1  P1  Q1  MF  MM SSA
               at each VLBI site.(M)
               J = 1,...,11- number of tide
               K = 1,2,3   - Vertical, East-West, North-South components !!!!
               L = 1,...N  - number of sites
            2) phs_ocean(11,3,N) - ocean loading phases. (rad)
       /eops/ 
            1) eops(7,7) -  on dates closest to date of observation
               eops(1:7,1) = MJD 
               eops(1:7,2) = UT1-UTC (Sec)
               eops(1:7,3) = UT1-TAI (Sec)
               eops(1:7,4) = x       (Arcsec)
               eops(1:7,5) = y       (Arcsec)
               eops(1:7,6) = dX      (Arcsec)
               eops(1:7,7) = dY      (Arcsec)
    '''
    ## Source info
    sou = source(sou_name, sou_type) # init source object
    # check source name
#    if sou.name.lower()[0] == 'j' or sou.name.lower()[0] == 'w':
    if sou.sou_type=='C':
        with open(inp['source_nam'], 'r') as f: #open cat with alternative sou names
            f_lines = f.readlines()
        # remove comments
        f_lines = [l for l in f_lines if l[0]!='#']
        # look up alternative names, pick the IVS one
        try:
            sou.ivsname = [x[0:8].strip() for x in f_lines \
                                if sou.name.lower() in x.lower()][0]
        except:
            pass
#        for f_line in f_lines:
#            if f_line[0]!='#' and (sou.name.lower() in f_line.lower()):
#                sou.ivsname = f_line[0:8]
#                break
    # far-field source position:
    if sou.sou_type=='C':
        with open(inp['source_cat'], 'r') as f: #open cat with alternative sou names
            f_lines = f.readlines()
        if not any(sou.ivsname.lower() in s.lower() for s in f_lines):
            if sou_radec==None:
                raise Exception('Source position for '+sou.name+\
                                ' not found in '+inp['source_cat']+'.')
            else:
                print 'Source position for '+sou.name+\
                      ' not found in '+inp['source_cat']+'. '+\
                      'Using the one from the vex-file.'
#                sou.radec.append(sou_radec[0])
#                sou.radec.append(sou_radec[1])
                sou.ra = sou_radec[0]
                sou.dec = sou_radec[1]
        else:
            for f_line in f_lines:
                if f_line[0]!='#' and (sou.ivsname.lower() in f_line.lower()):
                    #ra:
                    sou.radec.append([float(i) for i in f_line[14:29].strip().split()]) 
                    #dec:
                    sou.radec.append([float(i) for i in f_line[34:49].strip().split()])
            # ra/dec in rad:
            sou.ra = pi*(sou.radec[0][0] + sou.radec[0][1]/60.0 +
                         sou.radec[0][2]/3600.0)/12.0
            sou.dec = pi*(abs(sou.radec[1][0]) + sou.radec[1][1]/60.0 +
                          sou.radec[1][2]/3600.0)/180.0
            if sou.radec[1][0]<0: sou.dec *= -1
        # J2000.0 source unit vector:        
        sou.K_s = np.array([cos(sou.dec)*cos(sou.ra),
                            cos(sou.dec)*sin(sou.ra),
                            sin(sou.dec)])
    
    ## Station info
    sta = [site(sta_name) for sta_name in sta_names] # init site objects
    # first, get the cat epoch:
    with open(inp['sta_xyz'], 'r') as f:
        f_lines = f.readlines()
    EPOCH = None
    for f_line in f_lines:
        if 'EPOCH' in f_line:
            ep_vec = datetime.datetime.strptime(f_line[f_line.index('EPOCH')+6:-1], '%Y.%d.%m')
            ep_0 = datetime.datetime.strptime(f_line[f_line.index('EPOCH')+6:\
                                        f_line.index('EPOCH')+10], '%Y')
            ep_0p1 = datetime.datetime(ep_0.year+1,1,1)
            decim = float((ep_vec - ep_0).days) / float((ep_0p1 - ep_0).days)
            EPOCH = ep_0.year + decim
#    if 'EPOCH' not in locals(): # var EPOCH doesn't exist:
    if EPOCH is None:
        raise Exception('Station positions catalogue EPOCH not found. Check '\
               +inp['sta_xyz']+' file')
    
    # load cats for station data:
    with open(inp['sta_xyz'], 'r') as f: xyz_lines = f.readlines()   
    with open(inp['sta_vxvyvz'], 'r') as f: vxvyvz_lines = f.readlines()
    with open(inp['sta_thermdef'], 'r') as f: thermdef_lines = f.readlines()
    with open(inp['oc_load'], 'r') as f: oclo_lines = f.readlines()

    for st in sta:
    # xyz
    # should include a proper handling of the date columns, even though it's
    # not yet needed
        if st.name.upper() != 'GEOCENTR' and st.name.upper() != 'RA':
            if not any(st.name in s for s in xyz_lines):
                raise Exception('Station position for '+st.name+' not found.')
            # st.name might be in the cat multiple times:
            matching = [s for s in xyz_lines if st.name in s and len(s)>0 and s[0]!='$'][0]
            st.r_GTRS = np.array([float(i) for i in matching[15:61].split()])
    # vxvyvz
            if not any(st.name in s for s in vxvyvz_lines):
                raise Exception('Station velocity for '+st.name+' not found.')
            # st.name might be in the cat multiple times:
            matching = [s for s in vxvyvz_lines if st.name in s and len(s)>0 and s[0]!='$'][0]
            st.v_GTRS = np.array([float(i) for i in matching[15:61].split()])
            st.v_GTRS = st.v_GTRS*1e-3 # mm/year -> m/year
    # antennae Thermal deformations coefficients + some auxiliary info, e.g. mount type
            if not any(st.name in s for s in thermdef_lines):
                #raise Exception('Antenna therm_def data for '+st.name+' not found.')
                print 'Antenna therm_def and auxiliary data for '+st.name+' not found.'
                print 'Set mount type for {:s} to AltAz.'.format(st.name)
                st.mount_type = 'MO_AZEL'
            else:
                for f_line in thermdef_lines:
                    if len(f_line)>0 and f_line[0]!='#' and \
                                ('ANTENNA_INFO' in f_line) and (st.name in f_line):
                        matching = f_line
                        break
                # alternative way:
                # matching = [s for s in thermdef_lines if st.name in s and len(s)>0 and s[0]!='#'][-1]
                strz = matching[0:55].split()         
                fltz = [float(i) for i in matching[55:].split()]
                st.ivs_name, st.focus_type, \
                             st.mount_type, st.radome, st.meas_type = strz[1:]
                st.T0, st.sin_T, st.cos_T, st.h0, st.ant_diam, st.hf, st.df, \
                       st.gamma_hf, st.hp, st.gamma_hp, st.AO, st.gamma_AO, \
                       st.hv, st.gamma_hv, st.hs, st.gamma_hs = fltz    
    # load ocean loading parameters:
            for jj in range(len(oclo_lines)):
                f_line = oclo_lines[jj]
                if  len(f_line)>0 and f_line[0]!='$' and (st.name in f_line):
                    jj+=1
                    f_line = oclo_lines[jj] 
                    while f_line[0]=='$':
                        jj+=1
                        f_line = oclo_lines[jj]
                    for ii in range(3):
                        st.amp_ocean[:,ii] = np.array([float(i) \
                                             for i in oclo_lines[jj+ii].split()])
                        st.phs_ocean[:,ii] = np.array([float(i) \
                                             for i in oclo_lines[jj+ii+3].split()])
                    break
    # load atmospheric loading parameters: (to be added..)

    # Compute dyear - difference btw the epoch of the first observation 
    # and the epoch of the station catalog (necessary for calculation of 
    # tectonic motion of the sites)
    UTC_start = (float(date_t_start.hour) + date_t_start.minute/60.0 + \
                 date_t_start.second/3600.0) / 24.0
    dyear = ((date_t_start-ep_0).days + UTC_start) / 365.25
    
    # sta crds in TRF to date of obs:
    for st in sta: 
        st.r_GTRS_date = st.r_GTRS + dyear*st.v_GTRS
    
    ## load eops
    # get the relevant eop entries from the catalogue:
    with open(inp['cat_eop'], 'r') as fc:
        fc_lines = fc.readlines()
    mjd_start = mjuliandate(date_t_start.year, \
                            date_t_start.month, date_t_start.day)
    eops = np.zeros((7,7)) # +/- 3 days
    for jj in range(len(fc_lines)):
        if fc_lines[jj][0]!=' ' and fc_lines[jj][0]!='*':
            entry = [float(x) for x in fc_lines[jj].split()]
            if len(entry) > 0 and entry[3] == np.floor(mjd_start) - 3:
                for kk in range(7):
                    eops[kk,0] = entry[3] # mjd
                    eops[kk,1] = entry[6] # UT1-UTC
                    eops[kk,2] = entry[6] - nsec(entry[3]) # UT1 - TAI
                    eops[kk,3] = entry[4] # x_p
                    eops[kk,4] = entry[5] # y_p
                    eops[kk,5] = entry[8] # dX
                    eops[kk,6] = entry[9] # dY
                    entry = [float(x) for x in fc_lines[jj+kk+1].split()]
                break #exit loop
                    
    # done
    return sou, sta, eops
    
    
#==============================================================================
# 
#==============================================================================
def geodetic(sta, const):
    
    '''
    function geodetic calculates the site position in geodetic coordinate
    systems for the stations participating in the current observation.
    Transformation VW from local geodetic coordinate system (Vertical,East,North)
    to the Earth-fixed coordinate system is calculated for each cite.
    % 
    For transformation to geodetic coordinate system the 
    IAG1999 Reference Ellipsoid (Table 1.1, IERS 2000 Conventions) is used
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
    
    for st in sta:
        if st.name == 'RA' or st.name == 'GEOCENTR':
            continue
        # Compute the site spherical radius
        st.sph_rad = sqrt( sum(x*x for x in st.r_GTRS_date) )
        # Compute geocentric latitude
        st.lat_gcen = asin( st.r_GTRS_date[2] / st.sph_rad )
        # Compute geocentric longitude       
        st.lon_gcen = atan2( st.r_GTRS_date[1], st.r_GTRS_date[0] )
        if st.lon_gcen < 0.0:
            st.lon_gcen = st.lon_gcen + 2.0*pi
        # Compute the stations distance from the Earth spin axis and
        # from the equatorial plane (in KM)
        req = sqrt( sum(x*x for x in st.r_GTRS_date[:-1]) )
        st.u = req*1e-3
        st.v = st.r_GTRS_date[2]*1e-3
        # Compute geodetic latitude and height.
        # The geodetic longitude is equal to the geocentric longitude
        st.lat_geod, st.h_geod = geoid( req, st.r_GTRS_date[2], \
                                        const.AE, const.F )
        # Compute the local VEN-to-crust-fixed rotation matrices by rotating
        # about the geodetic latitude and the longitude. 
        # w - rotation matrix by an angle lat_geod around the y axis
        w = R_123(2, st.lat_geod)
        # v - rotation matrix by an angle -lon_gcen around the z axis
        v = R_123(3, -st.lon_gcen)
        # product of the two matrices:        
        st.vw = dot(v, w)
        
    return sta


def geoid(r, z, a, fr):
    '''
    Program to transform Cartesian to geodetic coordinates
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
    
#==============================================================================
# 
#==============================================================================
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
    
    
#==============================================================================
# 
#==============================================================================
def load_sc_eph(sou_type, source, t_start, t_end, inp,
                uvw_calc=False, sc_rhophitheta=False, sc_xyz=False,
                load=True, forcePaddLeft=None, forcePaddRight=None):
    '''
    (Down)Load spacecraft ephemerides into an ephem-class object
    Might be worth recoding some parts of this into class ephem methods
    '''
    # load constants:
#    const = constants()

    # split if observations go overnight
    # take padding into account (mind the LT!)
    if forcePaddLeft is None:
        if sou_type=='S' and (source.lower()!='gaia' and source.lower()!='ce3'):        
            if source.lower()=='mex': paddLeft=50
            if source.lower()=='vex': paddLeft=30
            if source.lower()=='her': paddLeft=10
            if source.lower()=='rosetta': paddLeft=90
        else:
            paddLeft = 1 # much more than enough for GNSS/RA/Gaia
    else:
        paddLeft = forcePaddLeft
    t_start_padd = t_start - datetime.timedelta(minutes=paddLeft)
    
    # tdb is more than 1 minute ahead of utc. padd the right-hand side therefore!
    if forcePaddRight is None:
        paddRight = 2 # should do the job
    else:
        paddRight = forcePaddRight
    t_end_padd = t_end + datetime.timedelta(minutes=paddRight)
    
    dd = ( datetime.datetime(t_end_padd.year, t_end_padd.month, t_end_padd.day) - \
           datetime.datetime(t_start_padd.year, t_start_padd.month, \
                             t_start_padd.day) ).days
    # start/stop times broken into days
    days = []
    for di in range(dd+1):
        days.append([ datetime.datetime(t_start_padd.year, t_start_padd.month, \
                       t_start_padd.day) + datetime.timedelta(days=di), 
                      datetime.datetime(t_start_padd.year, t_start_padd.month, \
                       t_start_padd.day, 23, 59, 59) + datetime.timedelta(days=di) ])
    days = np.array(days)
#    print days
    # fix for t_start and t_end for Gaia and RadioAstron
    if source.lower()=='gaia' or source.lower()=='ra':
        days[0,0] = t_start_padd
        days[-1,1] = t_end_padd

    # eph_file_names for each (full in VEX/MEX case) day
    eph_files = []
    for t_s, t_e in days:
        # ESA's spacecraft (VEX, MEX, HERSHEL):
        if sou_type=='S' and (source.lower()!='gaia' and source.lower()!='ce3'):        
            eph_file_names = esa_sc_eph_make(source, t_s, t_e, inp, paddLeft=0, paddRight=0)
        # RadioAstron:
        elif (sou_type=='C' and source.lower()=='ra') or sou_type=='R':
            eph_file_names = ra_eph_down(source, t_s, t_e, inp)
        # Gaia, Chang'E-3:
        elif sou_type=='S' and (source.lower()=='gaia' or source.lower()=='ce3'):
            eph_file_names = ra_eph_down(source, t_s, t_e, inp)
        # GNSS:
        elif sou_type=='G':
#            raise NotImplemented
            eph_file_names = ra_eph_down(source, t_s, t_e, inp)
        # append to eph_files name list
        eph_files.append(eph_file_names)

#    # RadioAstron:
#    if (sou_type=='C' and source.lower()=='ra') or sou_type=='R':
#        eph_file_names = ra_eph_down(source, t_start, t_end, inp)
#    # Gaia, Chang'E-3:
#    elif sou_type=='S' and (source.lower()=='gaia' or source.lower()=='ce3'):
#        eph_file_names = ra_eph_down(source, t_start, t_end, inp)
#    # ESA's spacecraft (VEX, MEX, HERSHEL):
#    elif sou_type=='S' and (source.lower()!='gaia' and source.lower()!='ce3'):
##        eph_file_names = esa_sc_eph_down(source, t_start, t_end, inp)
#        if source.lower()=='mex': padd=50
#        if source.lower()=='vex': padd=30
#        eph_file_names = esa_sc_eph_make(source, t_start, t_end, inp, \
#                                         paddLeft=padd)
#    # GNSS:
#    elif sou_type=='G':
#        raise NotImplemented
    
    # return if user doesn't want to load that shit
    if not load:
        return
    
    ''' now load the ephs into an ephem-object '''
    if sou_type=='C': source = 'RA'
    eph = ephem(source) # initialise ephem-object
    # load each day separately, stack the result together
    for f_bcrs, f_gtrs, f_gcrs in eph_files:
        # bcrs:
        with open(os.path.join(inp['sc_eph_cat'], f_bcrs),'r') as f:
            f_float = []
            for f_line in f.readlines():
                f_float.append([float(x) for x in f_line[:-1].split()])
            try:
                eph.bcrs[0] = np.vstack((eph.bcrs[0], np.array(f_float)))
            except:
                eph.bcrs.append(np.array(f_float))
        # gtrs:
        with open(os.path.join(inp['sc_eph_cat'], f_gtrs),'r') as f:
            f_float = []
            for f_line in f.readlines():
                f_float.append([float(x) for x in f_line[:-1].split()])
            try:
                eph.gtrs = np.vstack((eph.gtrs, np.array(f_float)))
            except:
                eph.gtrs = np.array(f_float)
        # gcrs:
        with open(os.path.join(inp['sc_eph_cat'], f_gcrs),'r') as f:
            f_float = []
            for f_line in f.readlines():
                f_float.append([float(x) for x in f_line[:-1].split()])
            try:
                eph.gcrs = np.vstack((eph.gcrs, np.array(f_float)))
            except:
                eph.gcrs = np.array(f_float)

    # convert to metres
    eph.bcrs[0][:,6:] *= 1e3
    eph.gcrs[:,6:] *= 1e3
    eph.gtrs[:,6:] *= 1e3

    # build up time stamps:
    # decimal days:
    eph.UT = (eph.gtrs[:,3] + eph.gtrs[:,4]/60.0 + eph.gtrs[:,5]/3600.0)/24.0
    # allow for 'overnighters':
    t0 = datetime.datetime(t_start.year, t_start.month, t_start.day)
    for nn, _ in enumerate(eph.UT):
        # should be wrt day of t_start (t0), i.e. allow for negative values
        # if padding goes one day earlier
        eph.UT[nn] += (datetime.datetime(*map(int,eph.gtrs[nn,0:3])) - t0).days

    eph.CT = (eph.bcrs[0][:,3] + eph.bcrs[0][:,4]/60.0 + \
              eph.bcrs[0][:,5]/3600.0)/24.0
    eph.CT_sec = eph.bcrs[0][:,3]*3600.0 + eph.bcrs[0][:,4]*60.0 + \
                 eph.bcrs[0][:,5]

    # allow for 'overnighters':
    for nn in range(len(eph.CT)):
        # should be wrt t_start, i.e. allow for negative values
        # if padding goes one day earlier
        eph.CT[nn] += (datetime.datetime(*map(int,eph.bcrs[0][nn,0:3])) - t0).days
        eph.CT_sec[nn] += 86400.0*\
                   (datetime.datetime(*map(int,eph.bcrs[0][nn,0:3])) - t0).days
#    eph.CT_sec = 86400.0*eph.CT

    # cut a proper piece
    start = (t_start_padd - t0).total_seconds()/86400.0
    stop  = (t_end_padd - t0).total_seconds()/86400.0
    # create 1D mask
    # TOOD: check correctness
    mask = np.ma.masked_inside(eph.CT, start, stop).mask
    mask = np.array(mask, dtype=bool)
    # remove masked points
#    print eph.CT.shape
    eph.CT = eph.CT[mask]
    eph.CT_sec = eph.CT_sec[mask]
    eph.UT = eph.UT[mask]
#    print eph.CT.shape    
    # remove rows in ephs according to 1D CT mask
#    print eph.gtrs.shape
#    # (first, turn 1D mask into a 2D matrix)
#    mask2d = np.expand_dims(mask, 0).T # turn a row into a column
#    mask2d = np.repeat(mask2d, eph.gtrs.shape[1])
#    eph.gtrs = np.ma.array(eph.gtrs, mask = mask2d)
#    eph.gtrs = np.ma.compress_rows(eph.gtrs)
    eph.gtrs = eph.gtrs[mask,:]
    eph.gcrs = eph.gcrs[mask,:]
    eph.bcrs[0] = eph.bcrs[0][mask,:]
#    print eph.gtrs.shape

    # if perturbed orbits for the s/c are needed, calculate them:
    if sou_type!='C' and uvw_calc:
        if sc_rhophitheta:
            eph = sc_uvw_eph(eph, inp)
        elif sc_xyz:
            eph = sc_xyz_eph(eph, inp)

    ## make spline interpolants
    # gtrs:
#    for jj in range(6,12):
#        f = sp.interpolate.splrep(eph.UT, eph.gtrs[:,jj], s=0)
#        eph.fGtrs.append(f)
    # gcrs:
#    for jj in range(6,15):
#        f = sp.interpolate.splrep(eph.UT, eph.gcrs[:,jj], s=0)
#        eph.fGcrs.append(f)
    # bcrs:
#    for bce in eph.bcrs:
#        fBce = []
#        for jj in range(6,12):
#            f = sp.interpolate.splrep(eph.CT, bce[:,jj], s=0)
#            fBce.append(f)
#        eph.fBcrs.append(fBce) # list of lists
    # evaluate: sp.interpolate.splev(xnew,eph.fGtrs[jj],der=0)
        
    # return the eph object
    return eph
    
    
#==============================================================================
# 
#==============================================================================
def sc_uvw_eph(eph, inp):
    '''
    Make perturbed orbits for the s/c uwvs calculation
    '''
    jpl_eph = inp['jpl_eph']
    mas_step = inp['mas_step']
    m_step = inp['m_step']
    JD = 2400000.5 + mjuliandate(*eph.bcrs[0][0,0:3])
    # calculate state vector of the Earth at JD+CT, in [m, m/s]
    earth_ssb = [pleph(JD+CT, 3, 12, jpl_eph)[0:3] * 1e3 for CT in eph.CT]
    sc_gc = eph.bcrs[0][:,6:9] - earth_ssb

    # convert sc_gc to spherical coordinates (transform to theta,phi,rho):
    sc_gc_rpt = cart2sph(sc_gc)

    # allocate space for additional orbits:
#    for jj in range(6):
    for jj in range(3):
        eph.bcrs.append(np.copy(eph.bcrs[0]))

    # rho + m_step m
    tmp = np.copy(sc_gc_rpt)
    tmp[:,0] += m_step
    eph.bcrs[1][:,6:9] = sph2cart(tmp) + earth_ssb
    
    # phi + mas_step
    tmp = np.copy(sc_gc_rpt)
    tmp[:,1] += mas_step*(pi/648000000.0)
    eph.bcrs[2][:,6:9] = sph2cart(tmp) + earth_ssb
    
    # theta + mas_step
    tmp = np.copy(sc_gc_rpt)
    tmp[:,2] += mas_step*(pi/648000000.0)
    eph.bcrs[3][:,6:9] = sph2cart(tmp) + earth_ssb
    
#    # rho +- m_step m
#    tmp = np.copy(sc_gc_rpt)
#    tmp[:,0] += m_step
#    #sph2cart(np.hstack((sc_gc_rpt[0:1],sc_gc_rpt[2]+m_step)))
#    eph.bcrs[1][:,6:9] = sph2cart(tmp) + earth_ssb
#    
#    tmp = np.copy(sc_gc_rpt)
#    tmp[:,0] -= m_step
#    eph.bcrs[2][:,6:9] = sph2cart(tmp) + earth_ssb
#    
#    # phi +- mas_step
#    tmp = np.copy(sc_gc_rpt)
#    tmp[:,1] += mas_step*(pi/648000000.0)
#    eph.bcrs[3][:,6:9] = sph2cart(tmp) + earth_ssb
#    
#    tmp = np.copy(sc_gc_rpt)
#    tmp[:,1] -= mas_step*(pi/648000000.0)
#    eph.bcrs[4][:,6:9] = sph2cart(tmp) + earth_ssb
#    
#    # theta +- mas_step
#    tmp = np.copy(sc_gc_rpt)
#    tmp[:,2] += mas_step*(pi/648000000.0)
#    eph.bcrs[5][:,6:9] = sph2cart(tmp) + earth_ssb
#    
#    tmp = np.copy(sc_gc_rpt)
#    tmp[:,2] -= mas_step*(pi/648000000.0)
#    eph.bcrs[6][:,6:9] = sph2cart(tmp) + earth_ssb
    
#    plt.subplot(2,1,1)
#    plt.plot(eph.CT, eph.bcrs[0][:,6]-eph.bcrs[5][:,6],\
#            eph.CT, eph.bcrs[0][:,7]-eph.bcrs[5][:,7],\
#            eph.CT, eph.bcrs[0][:,8]-eph.bcrs[5][:,8])
#    plt.subplot(2,1,2)
#    plt.plot(eph.CT, eph.bcrs[0][:,6]-eph.bcrs[6][:,6],\
#            eph.CT, eph.bcrs[0][:,7]-eph.bcrs[6][:,7],\
#            eph.CT, eph.bcrs[0][:,8]-eph.bcrs[6][:,8])
#    plt.show()
    
    return eph


#==============================================================================
# 
#==============================================================================
def cart2sph(xyz):
    '''
    Cartesian to spherical crd transformation
    Input - an N (rows) by 3 (columns) array
    '''
    if type(xyz) != type(np.zeros(3)):
        xyz = np.array(xyz)
    rpt = np.zeros(xyz.shape)
    try:
        xy = xyz[:,0]**2 + xyz[:,1]**2
        rpt[:,0] = np.sqrt(xy + xyz[:,2]**2) # rho
        # for elevation angle defined from Z-axis down:
        #rpt[:,1] = np.arctan2(np.sqrt(xy), xyz[:,2])
        # for elevation angle defined from XY-plane up:
        #rpt[:,1] = np.arctan2(xyz[:,2], np.sqrt(xy)) # phi
        rpt[:,1] = np.arctan2(xyz[:,2], np.sqrt(xy)) # phi ('elevation')
        rpt[:,2] = np.arctan2(xyz[:,1], xyz[:,0]) # theta ('azimuth')
    except IndexError: # 1d-array case
        xy = xyz[0]**2 + xyz[1]**2
        rpt[0] = sqrt(xy + xyz[2]**2) # rho
        rpt[1] = atan2(xyz[2], sqrt(xy)) # phi ('elevation')
        rpt[2] = atan2(xyz[1], xyz[0]) # theta ('azimuth')
    return rpt
    
def sph2cart(rpt):
    '''
    Spherical to cartesian crd transformation
    Input - an N (rows) by 3 (columns) array
    '''
    if type(rpt) != type(np.zeros(3)):
        rpt = np.array(rpt)
    xyz = np.zeros(rpt.shape)
    try:
        xyz[:,0] = rpt[:,0]*np.cos(rpt[:,1])*np.cos(rpt[:,2])
        xyz[:,1] = rpt[:,0]*np.cos(rpt[:,1])*np.sin(rpt[:,2])
        xyz[:,2] = rpt[:,0]*np.sin(rpt[:,1])
    except IndexError: # 1d-array case
        xyz[0] = rpt[0]*cos(rpt[1])*cos(rpt[2])
        xyz[1] = rpt[0]*cos(rpt[1])*sin(rpt[2])
        xyz[2] = rpt[0]*sin(rpt[1])
    return xyz


#==============================================================================
# 
#==============================================================================
def sc_xyz_eph(eph, inp):
    '''
    Make perturbed s/c orbits wrt xyz
    '''
    m_step = inp['m_step']

#    for jj in range(6):
    for jj in range(3):
        eph.bcrs.append(np.copy(eph.bcrs[0]))
    # x + m_step m
    eph.bcrs[1][:,6] += m_step
    # y + m_step m
    eph.bcrs[2][:,7] += m_step
    # z + m_step m
    eph.bcrs[3][:,8] += m_step
    
#    # x +- m_step m
#    eph.bcrs[1][:,6] += m_step
#    eph.bcrs[2][:,6] -= m_step
#    # y +- m_step m
#    eph.bcrs[3][:,7] += m_step
#    eph.bcrs[4][:,7] -= m_step
#    # z +- m_step m
#    eph.bcrs[5][:,8] += m_step
#    eph.bcrs[6][:,8] -= m_step
    
    return eph
    
    
#==============================================================================
# 
#==============================================================================
def dehanttideinel(sta, t, earth, sun, moon, r2000, calc_vel=True):
    '''
      - - - - - - - - - - - - - - - 
       D E H A N T T I D E I N E L
      - - - - - - - - - - - - - - - 
    
      This routine is part of the International Earth Rotation and
      Reference Systems Service (IERS) Conventions software collection.
    
      This subroutine computes the tidal corrections of station displacements
      caused by lunar and solar gravitational attraction (see References). 
      The computations are calculated by the following steps:
    
      Step 1): General degree 2 and degree 3 corrections + CALL ST1IDIU 
      + CALL ST1ISEM + CALL ST1L1.
      
      Step 2): CALL STEP2DIU + CALL STEP2LON
    
      It has been decided that the Step 3 non-correction for permanent tide
      would not be applied in order to avoid a jump in the reference frame.
      This Step 3 must be added in order to get the non-tidal station position
      and to conform with the IAG Resolution.
    
      In general, Class 1, 2, and 3 models represent physical effects that
      act on geodetic parameters while canonical models provide lower-level
      representations or basic computations that are used by Class 1, 2, or
      3 models.
     
      Status: Class 1
    
         Class 1 models are those recommended to be used a priori in the
         reduction of raw space geodetic data in order to determine
         geodetic parameter estimates.
         Class 2 models are those that eliminate an observational
         singularity and are purely conventional in nature.
         Class 3 models are those that are not required as either Class
         1 or 2.
         Canonical models are accepted as is and cannot be classified as a
         Class 1, 2, or 3 model.
    
      Given:
         sta - site object
         t - Datetime object with t_obs
         earth, sun, moon, - BCRS state vecs of Earth, Sun and Moon
         r2000 - rotation matrix ITRS -> ICRS
         calc_vel - numerically calculate velocity or not
         
       which is converted to:  
         XSTA          d(3)   Geocentric position of the station (Note 1)
         XSUN          d(3)   Geocentric position of the Sun (Note 2)
         XMON          d(3)   Geocentric position of the Moon (Note 2)
         YR            i      Year (Note 3)
         MONTH         i      Month (Note 3)
         DAY           i      Day of Month (Note 3)
         FHR           d      Hour in the day (Note 4)
    
      Returned:
         dx_tide        d(3)   Displacement vector (Note 5)
    
      Notes:
    
      1) The IGS station is in ITRF co-rotating frame.  All coordinates,
         X, Y, and Z, are expressed in meters. 
      
      2) The position is in Earth Centered Earth Fixed (ECEF) frame.  All
         coordinates are expressed in meters.
 !!!     This means that the earth-centered positions of the Sun and the Moon
         should be rotated by r2000T!!! Note that in the example below the
         radius-vectors are given in a wrong frame (~ICRS).
    
      3) The values are expressed in Coordinated Universal Time (UTC).
    
      4) The fractional hours in the day is computed as the hour + minutes/60.0
         + sec/3600.0.  The unit is expressed in Universal Time (UT).
    
      5) The displacement vector is in the geocentric ITRF.  All components are
         expressed in meters.
    
      Called:
         SPROD             Finds the scalar product and unit vector of two vectors 
         ZERO_VEC8         Returns the zero vector
         ST1IDIU           Corrects for the out-of-phase part of Love numbers
                           for the diurnal band
         ST1ISEM           Same as above for the semi-diurnal band
         ST1L1             Corrects for the latitude dependence of Love numbers
         CAL2JD            Computes Julian Date from Gregorian calendar date
         DAT               Computes the difference TAI-UTC
         STEP2DIU          Computes in-phase and out-of-phase corrections in
                           the diurnal band
         STEP2LON          Same as above for the long period band
    
      Test case:
         given input: XSTA(1) = 4075578.385 meters
                      XSTA(2) =  931852.890 meters
                      XSTA(3) = 4801570.154 meters   
                      XSUN(1) = 137859926952.015 meters
                      XSUN(2) = 54228127881.4350 meters
                      XSUN(3) = 23509422341.6960 meters
                      XMON(1) = -179996231.920342 meters
                      XMON(2) = -312468450.131567 meters
                      XMON(3) = -169288918.592160 meters
                      YR      = 2009
                      MONTH   = 4
                      DAY     = 13
                      FHR     = 0.00 hours 
                      
         expected output:  dx_tide(1) = 0.7700420357108125891e-01 meters
                           dx_tide(2) = 0.6304056321824967613e-01 meters
                           dx_tide(3) = 0.5516568152597246810e-01 meters
    
    % XSTA=[4075578.385 931852.890 4801570.154];
    % XSUN=[137859926952.015 54228127881.4350 23509422341.6960];
    % XMON=[-179996231.920342 -312468450.131567 -169288918.592160];
    % YR=2009; MONTH=4; DAY=13; FHR=0.0;
    % dehanttideinel(XSTA,YR,MONTH,DAY,FHR,XSUN,XMON,1)
      References:
    
         Groten, E., 2000, Geodesists Handbook 2000, Part 4,
         http://www.gfy.ku.dk/~iag/HB2000/part4/groten.htm. See also
         ''Parameters of Common Relevance of Astronomy, Geodesy, and
         Geodynamics," J. Geod., 74, pp. 134-140
    
         Mathews, P. M., Dehant, V., and Gipson, J. M., 1997, ''Tidal station
         displacements," J. Geophys. Res., 102(B9), pp. 20,469-20,477
    
         Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
         IERS Technical Note No. 36, BKG (2010)
         
         Pitjeva, E. and Standish, E. M., 2009, ''Proposals for the masses
         of the three largest asteroids, the Moon-Earth mass ratio and the
         Astronomical Unit," Celest. Mech. Dyn. Astr., 103, pp. 365-372
    
         Ries, J. C., Eanes, R. J., Shum, C. K. and Watkins, M. M., 1992,
         ''Progress in the Determination of the Gravitational Coefficient
         of the Earth," Geophys. Res. Lett., 19(6), pp. 529-531
    
      Revisions:
      1996 March    23 V. Dehant      Original code
                       P. M. Mathews
                       J. Gipson
      2000 May      17 V. Dehant      Last modifications
                       P. M. Mathews
      2006 February 06 J. Ray         Header comments modified to clarify
                                      input/output units and systems
      2006 February 06 J. Ray         Subroutine DUTC modified for leap
                                      second on 2006.0 and to correct 
                                      do 5 i=1,87 from 84 to 87
      2006 August   31 G. Petit       Correct DUTC for dates after 2007
      2007 June     20 H. Manche      Modified DUTC to correct past mistake
                                      and corrected DE line in subroutine
                                      STEP2DIU
      2007 October  23 H. Manche      Replace subroutines DUTC and FJLDY with
                       G. Petit       SOFA subroutines iau_CAL2JD and iau_DAT
                                      and correct time arguments of subroutine
                                      STEP2DIU
      2009 February 19 G. Petit       Update routine iau_DAT for 2009.0 leap
                                      second
      2009 August   06 B.E. Stetzler  Initial standardization of code 
      2009 August   07 B.E. Stetzler  Updated MASS_RATIO_SUN, 
                                      MASS_RATIO_MOON and RE to CBEs and
                                      provided a test case
      2009 August  07  B.E. Stetzler  Capitalized all variables for Fortran
                                      77 compatibility
      2009 September 01 B.E. Stetzler Removed 'iau_' from redistributed SOFA
                                      subroutines
    -----------------------------------------------------------------------
    '''
    XSTA = sta.r_GTRS_date
    YR = t.year
    MONTH = t.month
    DAY = t.day
    FHR = (t.hour + t.minute/60.0 + t.second/3600.0)
    XSUN = dot(r2000[:,:,0].T, sun[:,0] - earth[:,0]) # in GTRS
    XMON = dot(r2000[:,:,0].T, moon[:,0] - earth[:,0]) # in GTRS
    
# test values
#    XSTA=np.array([4075578.385, 931852.890, 4801570.154])
#    XSUN=np.array([137859926952.015, 54228127881.4350, 23509422341.6960])
#    XMON=np.array([-179996231.920342, -312468450.131567, -169288918.592160])
#    YR=2009; MONTH=4; DAY=13; FHR=0.0;
    
    #----------------------------------------------------------------------  
    # NOMINAL SECOND DEGREE AND THIRD DEGREE LOVE NUMBERS AND SHIDA NUMBERS  
    #----------------------------------------------------------------------  
    H20 = 0.6078
    L20 = 0.0847
    H3 = 0.292
    L3 = 0.015
    #----------------------------------------------------------------------  
    # SCALAR PRODUCT OF STATION VECTOR WITH SUN/MOON VECTOR  
    #----------------------------------------------------------------------  
    SCS = dot(XSTA,XSUN)
    RSUN = norm(XSUN)
    SCM = dot(XSTA,XMON)
    RSTA = norm(XSTA)
    RMON = norm(XMON)
    SCSUN=SCS/RSTA/RSUN
    SCMON=SCM/RSTA/RMON
    #----------------------------------------------------------------------   
    # COMPUTATION OF NEW H2 AND L2  
    #----------------------------------------------------------------------  
    COSPHI=sqrt(XSTA[0]**2+XSTA[1]**2)/RSTA
    H2=H20-0.0006*(1.0-3.0/2.0*COSPHI**2)
    L2=L20+0.0002*(1.0-3.0/2.0*COSPHI**2)

    # P2 term  
    P2SUN=3.0*(H2/2.0-L2)*SCSUN**2-H2/2.0
    P2MON=3.0*(H2/2.0-L2)*SCMON**2-H2/2.0

    # P3 term  
    P3SUN=5.0/2.0*(H3-3.0*L3)*SCSUN**3+3.0/2.0*(L3-H3)*SCSUN
    P3MON=5.0/2.0*(H3-3.0*L3)*SCMON**3+3.0/2.0*(L3-H3)*SCMON

    #----------------------------------------------------------------------  
    # TERM IN DIRECTION OF SUN/MOON VECTOR  
    #----------------------------------------------------------------------  
    X2SUN=3.0*L2*SCSUN
    X2MON=3.0*L2*SCMON
    X3SUN=3.0*L3/2.0*(5.0*SCSUN**2-1.0)
    X3MON=3.0*L3/2.0*(5.0*SCMON**2-1.0)
    #----------------------------------------------------------------------  
    # FACTORS FOR SUN/MOON USING IAU CURRENT BEST ESTIMATES (SEE REFERENCES) 
    #----------------------------------------------------------------------  
    MASS_RATIO_SUN=332946.0482
    MASS_RATIO_MOON=0.0123000371
    RE=6378136.6
    FAC2SUN=MASS_RATIO_SUN*RE*(RE/RSUN)**3
    FAC2MON=MASS_RATIO_MOON*RE*(RE/RMON)**3
    FAC3SUN=FAC2SUN*(RE/RSUN)
    FAC3MON=FAC2MON*(RE/RMON)
  
    # TOTAL DISPLACEMENT
    dx_tide = np.zeros(3)
    for I in range(3):
        dx_tide[I]=FAC2SUN*( X2SUN*XSUN[I]/RSUN + P2SUN*XSTA[I]/RSTA ) +\
                   FAC2MON*( X2MON*XMON[I]/RMON + P2MON*XSTA[I]/RSTA ) +\
                   FAC3SUN*( X3SUN*XSUN[I]/RSUN + P3SUN*XSTA[I]/RSTA ) +\
                   FAC3MON*( X3MON*XMON[I]/RMON + P3MON*XSTA[I]/RSTA )

    #+---------------------------------------------------------------------  
    # CORRECTIONS FOR THE OUT-OF-PHASE PART OF LOVE NUMBERS (PART H_2^(0)I  
    # AND L_2^(0)I )  
    #----------------------------------------------------------------------

    # FIRST, FOR THE DIURNAL BAND       
    XCORSTA = ST1IDIU(XSTA,XSUN,XMON,FAC2SUN,FAC2MON)
    dx_tide = dx_tide + XCORSTA

    # SECOND, FOR THE SEMI-DIURNAL BAND       
    XCORSTA = ST1ISEM(XSTA,XSUN,XMON,FAC2SUN,FAC2MON)
    dx_tide = dx_tide + XCORSTA

    #+---------------------------------------------------------------------
    # CORRECTIONS FOR THE LATITUDE DEPENDENCE OF LOVE NUMBERS (PART L^(1) )  
    #----------------------------------------------------------------------   
    XCORSTA = ST1L1(XSTA,XSUN,XMON,FAC2SUN,FAC2MON)
    dx_tide = dx_tide + XCORSTA
    
    if calc_vel: 
        dx_tide0 = dx_tide
    
    # CONSIDER CORRECTIONS FOR STEP 2  

    #+---------------------------------------------------------------------  
    # CORRECTIONS FOR THE DIURNAL BAND:  
    # 
    #  FIRST, WE NEED TO KNOW THE DATE CONVERTED IN JULIAN CENTURIES 
    #        
    #   1) CALL THE SUBROUTINE COMPUTING THE JULIAN DATE 
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
    JJM0 = 2400000.5
    JJM1 = mjuliandate( YR, MONTH, DAY )
    T=((JJM0-2451545.0)+JJM1+FHR/24.0)/36525.0

    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
    #   2) CALL THE SUBROUTINE COMPUTING THE CORRECTION OF UTC TIME  
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
    DTT = nsec(JJM1)
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    DTT = DTT + 32.184
    #     CONVERSION OF T IN TT TIME
    T=T+DTT/(3600.0*24.0*36525.0)

    #  SECOND, WE CAN CALL THE SUBROUTINE STEP2DIU, FOR THE DIURNAL BAND
    #  CORRECTIONS, (in-phase and out-of-phase frequency dependence):
    XCORSTA = STEP2DIU(XSTA,FHR,T)
    dx_tide=dx_tide+XCORSTA
    
#    np.set_printoptions(precision=16)
#    print dx_tide
    #  CORRECTIONS FOR THE LONG-PERIOD BAND,
    #  (in-phase and out-of-phase frequency dependence):  
    XCORSTA = STEP2LON(XSTA,T)
    dx_tide=dx_tide+XCORSTA
    
    # CONSIDER CORRECTIONS FOR STEP 3
    
    #----------------------------------------------------------------------
    # UNCORRECT FOR THE PERMANENT TIDE  
    #  
    #      SINPHI=XSTA(3)/RSTA  
    #      COSPHI=sqrt(XSTA(1)**2+XSTA(2)**2)/RSTA
    #      COSLA=XSTA(1)/COSPHI/RSTA  
    #      SINLA=XSTA(2)/COSPHI/RSTA  
    #      DR=-sqrt(5.0/4.0/PI)*H2*0.31460*(3.0/2.0*SINPHI**2-0.5)
    #      DN=-sqrt(5.0/4.0/PI)*L2*0.31460*3.0*COSPHI*SINPHI
    #      dx_tide(1)=dx_tide(1)-DR*COSLA*COSPHI+DN*COSLA*SINPHI
    #      dx_tide(2)=dx_tide(2)-DR*SINLA*COSPHI+DN*SINLA*SINPHI  
    #      dx_tide(3)=dx_tide(3)-DR*SINPHI      -DN*COSPHI
    #-----------------------------------------------------------------------

# check test values
#    print dx_tide
    
    #calculate dv_tide, if requested
    dv_tide = np.zeros(3)
    dx_td = []
    if (calc_vel):
        sec = datetime.timedelta(seconds=1)
        # +/- 1 second:
        for jj in (-1, 1):
            ts = t + jj*sec
            YR = ts.year
            MONTH = ts.month
            DAY = ts.day
            FHR = (ts.hour + ts.minute/60.0 + ts.second/3600.0)
            XSUN = dot((r2000[:,:,0] + jj*r2000[:,:,1]).T, \
                          sun[:,0] - earth[:,0] + \
                          jj*(sun[:,1] - earth[:,1])) # in GTRS
            XMON = dot((r2000[:,:,0] + jj*r2000[:,:,1]).T, \
                          moon[:,0] - earth[:,0] + \
                          jj*(moon[:,1] - earth[:,1])) # in GTRS
            
            JJM0 = 2400000.5
            JJM1 = mjuliandate( YR, MONTH, DAY )
            T=((JJM0-2451545.0)+JJM1+FHR/24.0)/36525.0
            DTT = nsec(JJM1)
            DTT = DTT + 32.184
            T=T+DTT/(3600.0*24.0*36525.0)
            XCORSTA = STEP2DIU(XSTA, FHR, T)
            dx_ti = dx_tide0 + XCORSTA
            XCORSTA = STEP2LON(XSTA, T)
            dx_ti = dx_ti + XCORSTA
            dx_td.append(dx_ti)
        dv_tide = (dx_td[1] - dx_td[0])/2.0 # m/s in GTRS
    
    # dx_tide and dv_tide are in GTRS. We want them in GCRS:
    sta.dr_tide = dot(r2000[:,:,0], dx_tide)
    sta.dv_tide = dot(r2000[:,:,1], dx_tide) + dot(r2000[:,:,0], dv_tide)
    
    return sta

#@jit(double[:](double[:], double[:], double[:], double, double))
def ST1IDIU (XSTA,XSUN,XMON,FAC2SUN,FAC2MON):
    '''
    +
      - - - - - - - - - - -
       S T 1 I D I U
      - - - - - - - - - - -
    
      This routine is part of the International Earth Rotation and
      Reference Systems Service (IERS) Conventions software collection.
    
      This subroutine gives the out-of-phase corrections induced by
      mantle anelasticity in the diurnal band. 
    
      In general, Class 1, 2, and 3 models represent physical effects that
      act on geodetic parameters while canonical models provide lower-level
      representations or basic computations that are used by Class 1, 2, or
      3 models.
     
      Status: Class 1
    
         Class 1 models are those recommended to be used a priori in the
         reduction of raw space geodetic data in order to determine
         geodetic parameter estimates.
         Class 2 models are those that eliminate an observational
         singularity and are purely conventional in nature.
         Class 3 models are those that are not required as either Class
         1 or 2.
         Canonical models are accepted as is and cannot be classified as a
         Class 1, 2, or 3 model.
    
      Given:
         XSTA          d(3)   Geocentric position of the IGS station (Note 1)
         XSUN          d(3)   Geocentric position of the Sun (Note 2)
         XMON          d(3)   Geocentric position of the Moon (Note 2)
         FAC2SUN       d      Degree 2 TGP factor for the Sun (Note 3)      
         FAC2MON       d      Degree 2 TGP factor for the Moon (Note 3) 
    
      Returned:
         XCORSTA       d(3)   Out of phase station corrections for diurnal band
    
      Notes:
    
      1) The IGS station is in ITRF co-rotating frame.  All coordinates are
         expressed in meters. 
      
      2) The position is in Earth Centered Earth Fixed (ECEF) frame.  All
         coordinates are expressed in meters.
    
      3) The expressions are computed in the main program.  TGP is the tide
         generated potential.  The units are inverse meters. 
    
      Test case:
         given input: XSTA(1) = 4075578.385D0 meters
                      XSTA(2) =  931852.890D0 meters
                      XSTA(3) = 4801570.154D0 meters   
                      XSUN(1) = 137859926952.015D0 meters
                      XSUN(2) = 54228127881.4350D0 meters
                      XSUN(3) = 23509422341.6960D0 meters
                      XMON(1) = -179996231.920342D0 meters
                      XMON(2) = -312468450.131567D0 meters
                      XMON(3) = -169288918.592160D0 meters
                      FAC2SUN =  0.163271964478954D0 1/meters     
                      FAC2MON =  0.321989090026845D0 1/meters    
                      
         expected output:  XCORSTA(1) = -0.2836337012840008001D-03 meters
                           XCORSTA(2) =  0.1125342324347507444D-03 meters
                           XCORSTA(3) = -0.2471186224343683169D-03 meters
    
      References:
    
         Mathews, P. M., Dehant, V., and Gipson, J. M., 1997, ''Tidal station
         displacements," J. Geophys. Res., 102(B9), pp. 20,469-20,477
    
         Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
         IERS Technical Note No. 36, BKG (2010)
    
      Revisions:
      1996 March    23 V. Dehant      Original code
      2009 July     30 B.E. Stetzler  Initial standardization of code 
      2009 July     31 B.E. Stetzler  Provided a test case
    -----------------------------------------------------------------------
    '''
    DHI = -0.0025
    DLI = -0.0007

    # Compute the normalized position vector of the IGS station.
    RSTA = norm(XSTA)
    SINPHI = XSTA[2]/RSTA
    COSPHI = sqrt(XSTA[0]*XSTA[0]+XSTA[1]*XSTA[1])/RSTA
    COS2PHI = COSPHI*COSPHI - SINPHI*SINPHI
    SINLA = XSTA[1]/COSPHI/RSTA
    COSLA = XSTA[0]/COSPHI/RSTA
    # Compute the normalized position vector of the Moon.
    RMON = norm(XMON)
    # Compute the normalized position vector of the Sun.
    RSUN = norm(XSUN)

    DRSUN=-3.0*DHI*SINPHI*COSPHI*FAC2SUN*XSUN[2]*(XSUN[0]*\
                 SINLA-XSUN[1]*COSLA)/RSUN**2

    DRMON=-3.0*DHI*SINPHI*COSPHI*FAC2MON*XMON[2]*(XMON[0]*\
                 SINLA-XMON[1]*COSLA)/RMON**2

    DNSUN=-3.0*DLI*COS2PHI*FAC2SUN*XSUN[2]*(XSUN[0]*SINLA-\
                 XSUN[1]*COSLA)/RSUN**2

    DNMON=-3.0*DLI*COS2PHI*FAC2MON*XMON[2]*(XMON[0]*SINLA-\
                 XMON[1]*COSLA)/RMON**2

    DESUN=-3.0*DLI*SINPHI*FAC2SUN*XSUN[2]*\
          (XSUN[0]*COSLA+XSUN[1]*SINLA)/RSUN**2

    DEMON=-3.0*DLI*SINPHI*FAC2MON*XMON[2]*\
        (XMON[0]*COSLA+XMON[1]*SINLA)/RMON**2

    DR = DRSUN + DRMON
    DN = DNSUN + DNMON
    DE = DESUN + DEMON

    #  Compute the corrections for the station.
    XCORSTA = np.zeros(3)
    XCORSTA[0] = DR*COSLA*COSPHI-DE*SINLA-DN*SINPHI*COSLA
    XCORSTA[1] = DR*SINLA*COSPHI+DE*COSLA-DN*SINPHI*SINLA
    XCORSTA[2] = DR*SINPHI+DN*COSPHI

    return XCORSTA


def ST1ISEM (XSTA,XSUN,XMON,FAC2SUN,FAC2MON):
    '''
    +
      - - - - - - - - - - -
       S T 1 I S E M
      - - - - - - - - - - -
    
      This routine is part of the International Earth Rotation and
      Reference Systems Service (IERS) Conventions software collection.
    
      This subroutine gives the out-of-phase corrections induced by
      mantle anelasticity in the semi-diurnal band. 
    
      In general, Class 1, 2, and 3 models represent physical effects that
      act on geodetic parameters while canonical models provide lower-level
      representations or basic computations that are used by Class 1, 2, or
      3 models.
     
      Status: Class 1
    
         Class 1 models are those recommended to be used a priori in the
         reduction of raw space geodetic data in order to determine
         geodetic parameter estimates.
         Class 2 models are those that eliminate an observational
         singularity and are purely conventional in nature.
         Class 3 models are those that are not required as either Class
         1 or 2.
         Canonical models are accepted as is and cannot be classified as a
         Class 1, 2, or 3 model.
    
      Given:
         XSTA          d(3)   Geocentric position of the IGS station (Note 1)
         XSUN          d(3)   Geocentric position of the Sun (Note 2)
         XMON          d(3)   Geocentric position of the Moon (Note 2)
         FAC2SUN       d      Degree 2 TGP factor for the Sun (Note 3)      
         FAC2MON       d      Degree 2 TGP factor for the Moon (Note 3) 
    
      Returned:
         XCORSTA       d(3)   Out of phase station corrections for
                              semi-diurnal band
    
      Notes:
    
      1) The IGS station is in ITRF co-rotating frame.  All coordinates are
         expressed in meters. 
      
      2) The position is in Earth Centered Earth Fixed (ECEF) frame.  All
         coordinates are expressed in meters.
    
      3) The expressions are computed in the main program.  TGP is the tide
         generated potential.  The units are inverse meters. 
    
      Test case:
         given input: XSTA(1) = 4075578.385D0 meters
                      XSTA(2) =  931852.890D0 meters
                      XSTA(3) = 4801570.154D0 meters   
                      XSUN(1) = 137859926952.015D0 meters
                      XSUN(2) = 54228127881.4350D0 meters
                      XSUN(3) = 23509422341.6960D0 meters
                      XMON(1) = -179996231.920342D0 meters
                      XMON(2) = -312468450.131567D0 meters
                      XMON(3) = -169288918.592160D0 meters
                      FAC2SUN =  0.163271964478954D0 1/meters     
                      FAC2MON =  0.321989090026845D0 1/meters    
                      
         expected output:  XCORSTA(1) = -0.2801334805106874015D-03 meters
                           XCORSTA(2) =  0.2939522229284325029D-04 meters
                           XCORSTA(3) = -0.6051677912316721561D-04 meters
    
      References:
    
         Mathews, P. M., Dehant, V., and Gipson, J. M., 1997, ''Tidal station
         displacements," J. Geophys. Res., 102(B9), pp. 20,469-20,477
    
         Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
         IERS Technical Note No. 36, BKG (2010)
    
      Revisions:
      1996 March    23 V. Dehant      Original code
      2009 July     31 B.E. Stetzler  Initial standardization of code 
      2009 July     31 B.E. Stetzler  Provided a test case
    -----------------------------------------------------------------------
    '''

    DHI = -0.0022; DLI = -0.0007;

    # Compute the normalized position vector of the IGS station.
    RSTA = norm(XSTA)
    SINPHI = XSTA[2]/RSTA
    COSPHI = sqrt(XSTA[0]*XSTA[0]+XSTA[1]*XSTA[1])/RSTA
    SINLA=XSTA[1]/COSPHI/RSTA
    COSLA=XSTA[0]/COSPHI/RSTA
    COSTWOLA=COSLA*COSLA-SINLA*SINLA
    SINTWOLA=2.0*COSLA*SINLA
    # Compute the normalized position vector of the Moon.
    RMON = norm(XMON)
    # Compute the normalized position vector of the Sun.
    RSUN = norm(XSUN)

    DRSUN=-3.0/4.0*DHI*COSPHI**2*FAC2SUN*((XSUN[0]**2-XSUN[1]**2)*\
        SINTWOLA-2.0*XSUN[0]*XSUN[1]*COSTWOLA)/RSUN**2

    DRMON=-3.0/4.0*DHI*COSPHI**2*FAC2MON*((XMON[0]**2-XMON[1]**2)*\
       SINTWOLA-2.0*XMON[0]*XMON[1]*COSTWOLA)/RMON**2

    DNSUN=3.0/2.0*DLI*SINPHI*COSPHI*FAC2SUN*((XSUN[0]**2-XSUN[1]**2)*\
       SINTWOLA-2.0*XSUN[0]*XSUN[1]*COSTWOLA)/RSUN**2

    DNMON=3.0/2.0*DLI*SINPHI*COSPHI*FAC2MON*((XMON[0]**2-XMON[1]**2)*\
       SINTWOLA-2.0*XMON[0]*XMON[1]*COSTWOLA)/RMON**2

    DESUN=-3.0/2.0*DLI*COSPHI*FAC2SUN*((XSUN[0]**2-XSUN[1]**2)*\
       COSTWOLA+2.0*XSUN[0]*XSUN[1]*SINTWOLA)/RSUN**2

    DEMON=-3.0/2.0*DLI*COSPHI*FAC2MON*((XMON[0]**2-XMON[1]**2)*\
       COSTWOLA+2.0*XMON[0]*XMON[1]*SINTWOLA)/RMON**2

    DR = DRSUN + DRMON
    DN = DNSUN + DNMON
    DE = DESUN + DEMON

    XCORSTA = np.zeros(3)
    XCORSTA[0] = DR*COSLA*COSPHI-DE*SINLA-DN*SINPHI*COSLA
    XCORSTA[1] = DR*SINLA*COSPHI+DE*COSLA-DN*SINPHI*SINLA
    XCORSTA[2] = DR*SINPHI+DN*COSPHI

    return XCORSTA

def ST1L1 (XSTA,XSUN,XMON,FAC2SUN,FAC2MON):
    '''
    +
      - - - - - - - - - - -
       S T 1 L 1
      - - - - - - - - - - -
    
      This routine is part of the International Earth Rotation and
      Reference Systems Service (IERS) Conventions software collection.
    
      This subroutine gives the corrections induced by the latitude 
      dependence given by L^1 in Mathews et al. 1991 (See References).
    
      In general, Class 1, 2, and 3 models represent physical effects that
      act on geodetic parameters while canonical models provide lower-level
      representations or basic computations that are used by Class 1, 2, or
      3 models.
     
      Status: Class 1
    
         Class 1 models are those recommended to be used a priori in the
         reduction of raw space geodetic data in order to determine
         geodetic parameter estimates.
         Class 2 models are those that eliminate an observational
         singularity and are purely conventional in nature.
         Class 3 models are those that are not required as either Class
         1 or 2.
         Canonical models are accepted as is and cannot be classified as a
         Class 1, 2, or 3 model.
    
      Given:
         XSTA          d(3)   Geocentric position of the IGS station (Note 1)
         XSUN          d(3)   Geocentric position of the Sun (Note 2)
         XMON          d(3)   Geocentric position of the Moon (Note 2)
         FAC2SUN       d      Degree 2 TGP factor for the Sun (Note 3)      
         FAC2MON       d      Degree 2 TGP factor for the Moon (Note 3) 
    
      Returned:
         XCORSTA       d(3)   Out of phase station corrections for
                              semi-diurnal band
    
      Notes:
    
      1) The IGS station is in ITRF co-rotating frame.  All coordinates are
         expressed in meters. 
      
      2) The position is in Earth Centered Earth Fixed (ECEF) frame.  All
         coordinates are expressed in meters.
    
      3) The expressions are computed in the main program. TGP is the tide
         generated potential.  The units are inverse meters. 
    
      Test case:
         given input: XSTA(1) = 4075578.385D0 meters
                      XSTA(2) =  931852.890D0 meters
                      XSTA(3) = 4801570.154D0 meters   
                      XSUN(1) = 137859926952.015D0 meters
                      XSUN(2) = 54228127881.4350D0 meters
                      XSUN(3) = 23509422341.6960D0 meters
                      XMON(1) = -179996231.920342D0 meters
                      XMON(2) = -312468450.131567D0 meters
                      XMON(3) = -169288918.592160D0 meters
                      FAC2SUN =  0.163271964478954D0 1/meters     
                      FAC2MON =  0.321989090026845D0 1/meters    
                      
         expected output:  XCORSTA(1) = 0.2367189532359759044D-03 meters
                           XCORSTA(2) = 0.5181609907284959182D-03 meters
                           XCORSTA(3) = -0.3014881422940427977D-03 meters
    
      References:
    
         Mathews, P. M., Buffett, B. A., Herring, T. A., Shapiro, I. I.,
         1991b, Forced nutations of the Earth: Influence of inner core
         Dynamics 2. Numerical results and comparisons, J. Geophys. Res.,
         96, 8243-8257
    
         Mathews, P. M., Dehant, V., and Gipson, J. M., 1997, ''Tidal station
         displacements," J. Geophys. Res., 102(B9), pp. 20,469-20,477
    
         Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
         IERS Technical Note No. 36, BKG (2010)
    
      Revisions:
      1996 March    23 V. Dehant      Original code
      2009 July     31 B.E. Stetzler  Initial standardization of code 
      2009 July     31 B.E. Stetzler  Provided a test case and Mathews
                                      reference
    -----------------------------------------------------------------------
    '''

    L1D = 0.0012
    L1SD = 0.0024

    # Compute the normalized position vector of the IGS station.
    RSTA = norm(XSTA)
    SINPHI = XSTA[2]/RSTA
    COSPHI = sqrt(XSTA[0]**2+XSTA[1]**2)/RSTA
    SINLA = XSTA[1]/COSPHI/RSTA
    COSLA = XSTA[0]/COSPHI/RSTA

    # Compute the normalized position vector of the Moon.
    RMON = norm(XMON)

    # Compute the normalized position vector of the Sun.
    RSUN = norm(XSUN)

    # Compute the station corrections for the diurnal band.
    L1=L1D
    DNSUN=-L1*SINPHI**2*FAC2SUN*XSUN[2]*(XSUN[0]*COSLA+XSUN[1]*SINLA)/RSUN**2
    DNMON=-L1*SINPHI**2*FAC2MON*XMON[2]*(XMON[0]*COSLA+XMON[1]*SINLA)/RMON**2
    DESUN=L1*SINPHI*(COSPHI**2-SINPHI**2)*FAC2SUN*XSUN[2]*\
       (XSUN[0]*SINLA-XSUN[1]*COSLA)/RSUN**2
    DEMON=L1*SINPHI*(COSPHI**2-SINPHI**2)*FAC2MON*XMON[2]*\
       (XMON[0]*SINLA-XMON[1]*COSLA)/RMON**2

    DE = 3.0*(DESUN+DEMON)
    DN = 3.0*(DNSUN+DNMON)

    XCORSTA = np.zeros(3)
    XCORSTA[0] = -DE*SINLA-DN*SINPHI*COSLA
    XCORSTA[1] = DE*COSLA-DN*SINPHI*SINLA
    XCORSTA[2] = DN*COSPHI
   
   # Compute the station corrections for the semi-diurnal band.
  
    L1=L1SD;
    COSTWOLA=COSLA**2-SINLA**2
    SINTWOLA=2.*COSLA*SINLA

    DNSUN=-L1/2.0*SINPHI*COSPHI*FAC2SUN*((XSUN[0]**2-XSUN[1]**2)*\
       COSTWOLA+2.0*XSUN[0]*XSUN[1]*SINTWOLA)/RSUN**2

    DNMON=-L1/2.0*SINPHI*COSPHI*FAC2MON*((XMON[0]**2-XMON[1]**2)*\
       COSTWOLA+2.0*XMON[0]*XMON[1]*SINTWOLA)/RMON**2

    DESUN=-L1/2.0*SINPHI**2*COSPHI*FAC2SUN*((XSUN[0]**2-XSUN[1]**2)*\
       SINTWOLA-2.0*XSUN[0]*XSUN[1]*COSTWOLA)/RSUN**2

    DEMON=-L1/2.0*SINPHI**2*COSPHI*FAC2MON*((XMON[0]**2-XMON[1]**2)*\
       SINTWOLA-2.0*XMON[0]*XMON[1]*COSTWOLA)/RMON**2

    DE = 3.0*(DESUN+DEMON)
    DN = 3.0*(DNSUN+DNMON) 

    XCORSTA[0]=XCORSTA[0]-DE*SINLA-DN*SINPHI*COSLA
    XCORSTA[1]=XCORSTA[1]+DE*COSLA-DN*SINPHI*SINLA
    XCORSTA[2]=XCORSTA[2]+DN*COSPHI

    return XCORSTA
    

def STEP2DIU (XSTA,FHR,T):
    '''
    +
      - - - - - - - - - - -
       S T E P 2 D I U
      - - - - - - - - - - -
    
      This routine is part of the International Earth Rotation and
      Reference Systems Service (IERS) Conventions software collection.
    
      This subroutine gives the in-phase and out-of-phase corrections
      induced by mantle anelasticity in the diurnal band. 
    
      In general, Class 1, 2, and 3 models represent physical effects that
      act on geodetic parameters while canonical models provide lower-level
      representations or basic computations that are used by Class 1, 2, or
      3 models.
     
      Status: Class 1
    
         Class 1 models are those recommended to be used a priori in the
         reduction of raw space geodetic data in order to determine
         geodetic parameter estimates.
         Class 2 models are those that eliminate an observational
         singularity and are purely conventional in nature.
         Class 3 models are those that are not required as either Class
         1 or 2.
         Canonical models are accepted as is and cannot be classified as a
         Class 1, 2, or 3 model.
    
      Given:
         XSTA          d(3)   Geocentric position of the IGS station (Note 1)
         FHR           d      Fractional hours in the day (Note 2)
         T             d      Centuries since J2000
    
      Returned:
         XCORSTA       d(3)   In phase and out of phase station corrections
                              for diurnal band (Note 4)
    
      Notes:
    
      1) The IGS station is in ITRF co-rotating frame.  All coordinates are
         expressed in meters. 
      
      2) The fractional hours in the day is computed as the hour + minutes/60.0
         + sec/3600.0.  The unit is expressed in Universal Time (UT).
    
      4) All coordinates are expressed in meters.
    
      Test case:
         given input: XSTA(1) = 4075578.385D0 meters
                      XSTA(2) =  931852.890D0 meters
                      XSTA(3) = 4801570.154D0 meters 
                      FHR     = 0.00D0 hours
                      T       = 0.1059411362080767D0 Julian centuries
                      
         expected output:  XCORSTA(1) = 0.4193085327321284701D-02 meters
                           XCORSTA(2) = 0.1456681241014607395D-02 meters
                           XCORSTA(3) = 0.5123366597450316508D-02 meters
    
      References:
    
         Mathews, P. M., Dehant, V., and Gipson, J. M., 1997, ''Tidal station
         displacements," J. Geophys. Res., 102(B9), pp. 20,469-20,477
    
         Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
         IERS Technical Note No. 36, BKG (2010)
    
      Revisions:
      1996 March    23 V. Dehant      Original code
      2009 July     31 B.E. Stetzler  Initial standardization of code 
      2009 August   06 B.E. Stetzler  Provided a test case
      2009 August   06 B.E. Stetzler  Capitalized all variables for 
                                      Fortran 77 compatibility
      2010 October  20 B.E. Stetzler  Input T corrected to be number of
                                      centuries since J2000
    -----------------------------------------------------------------------
    '''
    D2PI = 6.283185307179586476925287

    DATDI  = np.array([\
      [-3.0, 0.0, 2.0, 0.0, 0.0,-0.01, 0.0, 0.0, 0.0],\
      [-3.0, 2.0, 0.0, 0.0, 0.0,-0.01, 0.0, 0.0, 0.0],\
      [-2.0, 0.0, 1.0,-1.0, 0.0,-0.02, 0.0, 0.0, 0.0],\
      [-2.0, 0.0, 1.0, 0.0, 0.0,-0.08, 0.0,-0.01, 0.01],\
      [-2.0, 2.0,-1.0, 0.0, 0.0,-0.02, 0.0, 0.0, 0.0],\
      [-1.0, 0.0, 0.0,-1.0, 0.0,-0.10, 0.0, 0.0, 0.0],\
      [-1.0, 0.0, 0.0, 0.0, 0.0,-0.51, 0.0,-0.02, 0.03],\
      [-1.0, 2.0, 0.0, 0.0, 0.0, 0.01, 0.0, 0.0, 0.0],\
      [ 0.0,-2.0, 1.0, 0.0, 0.0, 0.01, 0.0, 0.0, 0.0],\
      [ 0.0, 0.0,-1.0, 0.0, 0.0, 0.02, 0.0, 0.0, 0.0],\
      [ 0.0, 0.0, 1.0, 0.0, 0.0, 0.06, 0.0, 0.0, 0.0],\
      [ 0.0, 0.0, 1.0, 1.0, 0.0, 0.01, 0.0, 0.0, 0.0],\
      [ 0.0, 2.0,-1.0, 0.0, 0.0, 0.01, 0.0, 0.0, 0.0],\
      [ 1.0,-3.0, 0.0, 0.0, 1.0,-0.06, 0.0, 0.0, 0.0],\
      [ 1.0,-2.0, 0.0,-1.0, 0.0, 0.01, 0.0, 0.0, 0.0],\
      [ 1.0,-2.0, 0.0, 0.0, 0.0,-1.23,-0.07, 0.06, 0.01],\
      [ 1.0,-1.0, 0.0, 0.0,-1.0, 0.02, 0.0, 0.0, 0.0],\
      [ 1.0,-1.0, 0.0, 0.0, 1.0, 0.04, 0.0, 0.0, 0.0],\
      [ 1.0, 0.0, 0.0,-1.0, 0.0,-0.22, 0.01, 0.01, 0.0],\
      [ 1.0, 0.0, 0.0, 0.0, 0.0,12.00,-0.80,-0.67,-0.03],\
      [ 1.0, 0.0, 0.0, 1.0, 0.0, 1.73,-0.12,-0.10, 0.0],\
      [ 1.0, 0.0, 0.0, 2.0, 0.0,-0.04, 0.0, 0.0, 0.0],\
      [ 1.0, 1.0, 0.0, 0.0,-1.0,-0.50,-0.01, 0.03, 0.0],\
      [ 1.0, 1.0, 0.0, 0.0, 1.0, 0.01, 0.0, 0.0, 0.0],\
      [ 0.0, 1.0, 0.0, 1.0,-1.0,-0.01, 0.0, 0.0, 0.0],\
      [ 1.0, 2.0,-2.0, 0.0, 0.0,-0.01, 0.0, 0.0, 0.0],\
      [ 1.0, 2.0, 0.0, 0.0, 0.0,-0.11, 0.01, 0.01, 0.0],\
      [ 2.0,-2.0, 1.0, 0.0, 0.0,-0.01, 0.0, 0.0, 0.0],\
      [ 2.0, 0.0,-1.0, 0.0, 0.0,-0.02, 0.0, 0.0, 0.0],\
      [ 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],\
      [ 3.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0 ]])

    DEG2RAD = D2PI/360.0

    #  Compute the phase angles in degrees.
    S = 218.31664563 + (481267.88194 + (-0.0014663889 + (0.00000185139)*T)*T)*T

    TAU = FHR*15.0 + 280.4606184 + (36000.7700536 +\
          (0.00038793 + (-0.0000000258)*T)*T)*T + (-S)

    PR = (1.396971278 + (0.000308889 + (0.000000021 + (0.000000007)*T)*T)*T)*T

    S = S + PR

    H = 280.46645 + (36000.7697489 + \
          (0.00030322222 + (0.000000020 + (-0.00000000654)*T)*T)*T)*T

    P = 83.35324312 + (4069.01363525 + \
          (-0.01032172222 + (-0.0000124991 + (0.00000005263)*T)*T)*T)*T

    ZNS = 234.95544499 + (1934.13626197 + \
          (-0.00207561111 + (-0.00000213944 + (0.00000001650)*T)*T)*T)*T

    PS = 282.93734098 + (1.71945766667 + \
          (0.00045688889 + (-0.00000001778 + (-0.00000000334)*T)*T)*T)*T

    # Reduce angles to between the range 0 and 360.
    S =  fmod(S,360.0)
    TAU = fmod(TAU,360.0)
    H =  fmod(H,360.0)
    P =  fmod(P,360.0)
    ZNS = fmod(ZNS,360.0)
    PS = fmod(PS,360.0)

    RSTA = sqrt(XSTA[0]**2+XSTA[1]**2+XSTA[2]**2)
    SINPHI = XSTA[2]/RSTA
    COSPHI = sqrt(XSTA[0]**2+XSTA[1]**2)/RSTA

    COSLA = XSTA[0]/COSPHI/RSTA
    SINLA = XSTA[1]/COSPHI/RSTA
    ZLA = atan2(XSTA[1],XSTA[0])

    XCORSTA = np.zeros(3)
    
    THETAF=(TAU+DATDI[:,0]*S+DATDI[:,1]*H+DATDI[:,2]*P+\
                DATDI[:,3]*ZNS+DATDI[:,4]*PS)*DEG2RAD
        
    DR=DATDI[:,5]*2.0*SINPHI*COSPHI*np.sin(THETAF+ZLA)+\
        DATDI[:,6]*2.0*SINPHI*COSPHI*np.cos(THETAF+ZLA)

    DN=DATDI[:,7]*(COSPHI**2-SINPHI**2)*np.sin(THETAF+ZLA)+\
        DATDI[:,8]*(COSPHI**2-SINPHI**2)*np.cos(THETAF+ZLA)

    DE=DATDI[:,7]*SINPHI*np.cos(THETAF+ZLA)-\
        DATDI[:,8]*SINPHI*np.sin(THETAF+ZLA)
        
    XCORSTA[0] = np.sum(DR*COSLA*COSPHI-DE*SINLA-DN*SINPHI*COSLA)
    XCORSTA[1] = np.sum(DR*SINLA*COSPHI+DE*COSLA-DN*SINPHI*SINLA)
    XCORSTA[2] = np.sum(DR*SINPHI+DN*COSPHI)
    
    XCORSTA = XCORSTA/1000.0

    return XCORSTA
    

def STEP2LON (XSTA,T):
    '''
    +
      - - - - - - - - - - -
       S T E P 2 L O N
      - - - - - - - - - - -
    
      This routine is part of the International Earth Rotation and
      Reference Systems Service (IERS) Conventions software collection.
    
      This subroutine gives the in-phase and out-of-phase corrections
      induced by mantle anelasticity in the long period band. 
    
      In general, Class 1, 2, and 3 models represent physical effects that
      act on geodetic parameters while canonical models provide lower-level
      representations or basic computations that are used by Class 1, 2, or
      3 models.
     
      Status: Class 1
    
         Class 1 models are those recommended to be used a priori in the
         reduction of raw space geodetic data in order to determine
         geodetic parameter estimates.
         Class 2 models are those that eliminate an observational
         singularity and are purely conventional in nature.
         Class 3 models are those that are not required as either Class
         1 or 2.
         Canonical models are accepted as is and cannot be classified as a
         Class 1, 2, or 3 model.
    
      Given:
         XSTA          d(3)   Geocentric position of the IGS station (Note 1)
         T             d      Centuries since J2000
    
      Returned:
         XCORSTA       d(3)   In phase and out of phase station corrections
                              for diurnal band (Note 2)
    
      Notes:
    
      1) The IGS station is in ITRF co-rotating frame.  All coordinates are
         expressed in meters. 
      
      2) All coordinates are expressed in meters.
    
      Test case:
         given input: XSTA(1) = 4075578.385D0 meters
                      XSTA(2) =  931852.890D0 meters
                      XSTA(3) = 4801570.154D0 meters 
                      T       = 0.1059411362080767D0 Julian centuries
                      
         expected output:  XCORSTA(1) = -0.9780962849562107762D-04 meters
                           XCORSTA(2) = -0.2236349699932734273D-04 meters
                           XCORSTA(3) =  0.3561945821351565926D-03 meters
    
      References:
    
         Mathews, P. M., Dehant, V., and Gipson, J. M., 1997, ''Tidal station
         displacements," J. Geophys. Res., 102(B9), pp. 20,469-20,477
    
         Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
         IERS Technical Note No. 36, BKG (2010)
    
      Revisions:
      1996 March    23 V. Dehant      Original code
      2009 August   07 B.E. Stetzler  Initial standardization of code
                                      and found unnecessary variables tau
                                      and fhr 
      2009 August   07 B.E. Stetzler  Provided a test case
      2009 August   07 B.E. Stetzler  Capitalized all variables for 
                                      Fortran 77 compatibility
      2010 October  20 B.E. Stetzler  Input T corrected to be number of 
                                      centuries since J2000
    -----------------------------------------------------------------------
    '''

    D2PI = 6.283185307179586476925287

    DATDI  = np.array([\
        [0, 0, 0, 1, 0,   0.47, 0.23, 0.16, 0.07],\
        [0, 2, 0, 0, 0,  -0.20,-0.12,-0.11,-0.05],\
        [1, 0,-1, 0, 0,  -0.11,-0.08,-0.09,-0.04],\
        [2, 0, 0, 0, 0,  -0.13,-0.11,-0.15,-0.07],\
        [2, 0, 0, 1, 0,  -0.05,-0.05,-0.06,-0.03]])


    DEG2RAD = D2PI/360.0

    #  Compute the phase angles in degrees.
    S = 218.31664563 + (481267.88194 + (-0.0014663889 + (0.00000185139)*T)*T)*T

    PR = (1.396971278 + (0.000308889 + (0.000000021 + (0.000000007)*T)*T)*T)*T

    S = S + PR

    H = 280.46645 + (36000.7697489 + (0.00030322222 + (0.000000020\
                         + (-0.00000000654)*T)*T)*T)*T

    P = 83.35324312 + (4069.01363525 + (-0.01032172222 + (-0.0000124991\
                         + (0.00000005263)*T)*T)*T)*T

    ZNS = 234.95544499 + (1934.13626197 + (-0.00207561111\
             + (-0.00000213944 + (0.00000001650)*T)*T)*T)*T

    PS = 282.93734098 + (1.71945766667 + (0.00045688889\
             + (-0.00000001778 + (-0.00000000334)*T)*T)*T)*T

    RSTA = sqrt(XSTA[0]**2+XSTA[1]**2+XSTA[2]**2)
    SINPHI = XSTA[2]/RSTA
    COSPHI = sqrt(XSTA[0]**2+XSTA[1]**2)/RSTA

    COSLA=XSTA[0]/COSPHI/RSTA
    SINLA=XSTA[1]/COSPHI/RSTA

    # Reduce angles to between the range 0 and 360.
    S =  fmod(S,360.0)
#    TAU = fmod(TAU,360.0)
    H =  fmod(H,360.0)
    P =  fmod(P,360.0)
    ZNS = fmod(ZNS,360.0)
    PS = fmod(PS,360.0)

    XCORSTA = np.zeros(3)
    
    THETAF=(DATDI[:,0]*S+DATDI[:,1]*H+DATDI[:,2]*P+\
                DATDI[:,3]*ZNS+DATDI[:,4]*PS)*DEG2RAD
        
    DR=DATDI[:,5]*(3.0*SINPHI**2-1.0)/2.0*np.cos(THETAF)+\
           DATDI[:,7]*(3.0*SINPHI**2-1.0)/2.0*np.sin(THETAF)

    DN=DATDI[:,6]*(COSPHI*SINPHI*2.0)*np.cos(THETAF)+\
           DATDI[:,8]*(COSPHI*SINPHI*2.0)*np.sin(THETAF)

    XCORSTA[0] = np.sum(DR*COSLA*COSPHI-DN*SINPHI*COSLA)
    XCORSTA[1] = np.sum(DR*SINLA*COSPHI-DN*SINPHI*SINLA)
    XCORSTA[2] = np.sum(DR*SINPHI+DN*COSPHI)
    
    XCORSTA /= 1000.0

    return XCORSTA
    
    
#==============================================================================
# 
#==============================================================================
def hardisp(sta, t, r2000, calc_vel=True):
    '''
     input:
     sta               site object with:
         sta.amp_ocean(11,3)   Amplitudes
         sta.phs_ocean(11,3)   Phases
         sta.vw                Transformation matrix from VEN to the 
                               Earth-fixed coordinate system
     r2000             Transformation matrix from crust-fixed to J2000 system
     calc_vel          Swich whether to calculate velocity numerically or not 
                       (as it's very small..)
    
     output:
      sta object with dx_octide, dv_octide saved in dr_oclo, dv_oclo
     
    '''    
    amp_ocean = sta.amp_ocean
    phs_ocean = sta.phs_ocean
    vw = sta.vw
    
    dx_octide = hardisp_calc(amp_ocean, phs_ocean, vw, t, r2000)
    
    # calculate velocity, if requested:
    dv_octide = np.zeros(3)
    
    # calculate velocity numerically:
    if (calc_vel):
        dx_ol = []
        sec = datetime.timedelta(seconds=1)
        # +/- 1 sec:
        for jj in (-1, 1):
            dx_ol.append(hardisp_calc(amp_ocean, phs_ocean, vw,\
                                      t + jj*sec, r2000))
        dv_octide = (dx_ol[1] - dx_ol[0])/2.0 # m/s in GCRS    
    
    sta.dr_oclo = np.array(dx_octide)
    sta.dv_oclo = np.array(dv_octide)
    
    return sta


def hardisp_calc(amp_ocean, phs_ocean, vw, t, r2000):
    '''
     input:
     amp_ocean(11,3)   Amplitudes
     phs_ocean(11,3)   Phases
     vw                Transformation matrix from VEN to the 
                       Earth-fixed coordinate system
     r2000             Transformation matrix from crust-fixed to J2000 system
     calc_vel          Swich whether to calculate velocity numerically or not 
                       (as it's very small..)
    
     output:
      dx_octide - displacement due to ocean loading
     
    '''
    # parameters
    NT = 342
    ntin = 11
          
    IDT = np.array([ [2, 0, 0, 0, 0, 0],  [2, 2,-2, 0, 0, 0],   [2,-1, 0, 1, 0, 0],\
            [2, 2, 0, 0, 0, 0],  [1, 1, 0, 0, 0, 0],   [1,-1, 0, 0, 0, 0],\
            [1, 1,-2, 0, 0, 0],  [1,-2, 0, 1, 0, 0],   [0, 2, 0, 0, 0, 0],\
            [0, 1, 0,-1, 0, 0],  [0, 0, 2, 0, 0, 0] ], order='F')
    IDT = np.transpose(IDT)
    
    year = t.year
    hh = t.hour
    mm = t.minute
    ss = t.second
    
    doy = (t - datetime.datetime(year,1,1)).days + 1
    it = np.array([year, doy, hh, mm, ss])

    # Find amplitudes and phases for all constituents, for each of the three 
    # displacements. Note that the same frequencies are returned each time.
    # 
    # BLQ format order is vertical, horizontal EW, horizontal NS
    AZ, _, PZ, _ = admint2(amp_ocean[:,0],IDT,-phs_ocean[:,0],ntin,it)
    AW, _, PW, _ = admint2(amp_ocean[:,1],IDT,-phs_ocean[:,1],ntin,it)
    AS, F, PS, ntout = admint2(amp_ocean[:,2],IDT,-phs_ocean[:,2],ntin,it)

    # set up for recursion, by normalizing frequencies, and converting
    # phases to radians
    
    DR = 0.01745329252
    
    PZ = DR*PZ
    PS = DR*PS
    PW = DR*PW
    F = pi*F/43200.0
    
    # Set up harmonic coefficients, and compute tide
    HCZ = np.zeros(2*NT)
    HCS = np.zeros(2*NT)
    HCW = np.zeros(2*NT)
    
    HCZ[::2]    = AZ*np.cos(PZ)
    HCZ[1::2]   = -AZ*np.sin(PZ)
    HCS[::2]    = AS*np.cos(PS)
    HCS[1::2]   = -AS*np.sin(PS)
    HCW[::2]    = AW*np.cos(PW)
    HCW[1::2]   = -AW*np.sin(PW)
    
    dz, _ = recurs(HCZ, ntout, F)
    ds, _ = recurs(HCS, ntout, F)
    dw, _ = recurs(HCW, ntout, F)

    # Change the signs on the horizontal displacements to convert the 
    # Vertical-West-South system to the Vertical-East-North system:
    dzne = np.array([dz, -ds, -dw])
    
    # Transformation of the displacements from local VEN to the 
    # Earth-fixed coordinate system
    # Then rotate vectors from crust-fixed to J2000 system.
    dx_octide = dot(r2000[:,:,0], dot(vw, dzne))
    
    return dx_octide
    
#@jit
#@numba.jit((double[:], double, double[:]))
@numba.jit('(f8[:], f8, f8[:])')
def recurs(HC,NF,OM):
    '''
   +
     - - - - - - - - -
      R E C U R S
     - - - - - - - - -
   
     This routine is part of the International Earth Rotation and
     Reference Systems Service (IERS) Conventions software collection.
   
     The purpose of the subroutine is to perform sine and cosine recursion
     to fill in data x, of length n, for nf sines and cosines with frequencies
     om. 
   
     In general, Class 1, 2, and 3 models represent physical effects that
     act on geodetic parameters while canonical models provide lower-level
     representations or basic computations that are used by Class 1, 2, or
     3 models.
    
     Status: Canonical model	
    
        Class 1 models are those recommended to be used a priori in the
        reduction of raw space geodetic data in order to determine
        geodetic parameter estimates.
        Class 2 models are those that eliminate an observational
        singularity and are purely conventional in nature.
        Class 3 models are those that are not required as either Class
        1 or 2.
        Canonical models are accepted as is and cannot be classified as a
        Class 1, 2, or 3 model.
   
     Given: This is a support routine of the main program HARDISP.F.
        x              d      data provided from a file given as standard
                              input from the MAIN program HARDISP.F (Note 1)
        n  = 1 !!!     i      length of the data file x
        hc             d      array containing alternating cosine and sine
                              coefficients
        nf             i      number of sine and cosine terms
        om             d      sine and cosine frequencies (Note 2)  
   
     Returned:
        scr            d      scratch array of length 3 times nf which is
                              returned as the recursion cr
     Notes:
   
     1) See the MAIN program HARDISP.F header comments for detailed information.
    
     2) The frequencies are normalized so that the Nyquist frequency is pi.
   
     Called:
        None
   
     Test case:
        Not provided for this subroutine.
   
     References:
   
        Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
        IERS Technical Note No. 36, BKG (2010)
   
     Revisions:
     2009 June 05 B.E. Stetzler    Added header and copyright, used Dcos
                                   and Dsin exclusively, and replaced END 
                                   DO statements with CONTINUE statements
     2009 August 19 B.E. Stetzler  Capitalized all variables for FORTRAN
                                   77 compatibility
   -----------------------------------------------------------------------
    '''
    #  Set up for start of recursion by computing harmonic values
    #  at starting point and just before it
    SCR = np.zeros(3*NF)
    SCR[::3] = HC[::2]
    SCR[1::3] = HC[::2]*np.cos(OM) - HC[1::2]*np.sin(OM)
    SCR[2::3] = 2.0*np.cos(OM)

    #  Do recursion over data
    X = 0.0
    #  Then do recursive computation for each harmonic
    for J in range(NF):
        X = X + SCR[3*J]
        SC = SCR[3*J]
        SCR[3*J] = SCR[3*J+2]*SC-SCR[3*J+1]
        SCR[3*J+1] = SC        
        
    return X, SCR
    
    
#==============================================================================
# 
#==============================================================================
def poletide(sta, t, eops, r2000, calc_vel=True):
    '''
    Rotarional deformation due to polar motion
    IERS 2010 Conventions, see Chapter 7

    input:
        date, lat_geod, lon_gcen, x_p, y_p, r2000
    output:
        dx_poltide, dv_poltide [m] in GCRS
    '''
    lat = sta.lat_geod
    lon = sta.lon_gcen
    
    dx_poltide = poletide_calc(lat, lon, t, eops, r2000)
    
    # calculate velocity, if requested:
    dv_poltide = np.zeros(3)
    
    # calculate velocity numerically:
    if (calc_vel):
        dx_pt = []
        sec = datetime.timedelta(seconds=1)
        # +/- 1 sec:
        for jj in (-1, 1):
            dx_pt.append(poletide_calc(lat, lon, t + jj*sec, eops, r2000))
        dv_poltide = (dx_pt[1] - dx_pt[0])/2.0 # m/s in GCRS

    sta.dr_poltide = dx_poltide
    sta.dv_poltide = dv_poltide

    return sta

#@jit
def poletide_calc(lat, lon, t, eops, r2000):
    
    x_p = eops[1]
    y_p = eops[2]
    
    mns = np.array([ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
    ds = 365.0
    if ( fmod(t.year, 4) == 0 ):
        mns[1] = 29
        ds = 366.0
    t = t.year + ( np.sum(mns[0:t.month-1]) + \
        t.day + (t.hour + t.minute/60.0 + t.second/3600.0)/24.0 ) / ds
    t_0 = 2000.0
    
    # until 2010:
    # stupid IERS' cubic model:
    if t<2010.0:
        x_p_mean_i = np.array([55.974, 1.8243, 0.18413, 0.007024])
        y_p_mean_i = np.array([346.346, 1.7896, -0.10729, -0.000908])
    # linear model after 2010.0:
    else:
        x_p_mean_i = np.array([23.513, 7.6141])
        y_p_mean_i = np.array([358.891, -0.6287])
    
    x_p_mean = 0.0
    y_p_mean = 0.0 # initialise
    for ii, vv in enumerate(x_p_mean_i):
        x_p_mean += vv * (t-t_0)**ii
    for ii, vv in enumerate(y_p_mean_i):
        y_p_mean += vv * (t-t_0)**ii

    x_p_mean *= 1e-3 # mas to as
    y_p_mean *= 1e-3 # mas to as

    m_1 = x_p - x_p_mean
    m_2 = y_p - y_p_mean
    # in Vertical-East-South system, [m]
    S_theta  =  -9.0 * cos(2.0*lat) * (m_1*cos(lon) + m_2*sin(lon)) * 1e-3
    S_lambda =   9.0 * cos(lat) * (m_1*sin(lon) - m_2*cos(lon)) * 1e-3
    S_r      = -33.0 * sin(2.0*lat) * (m_1*cos(lon)+m_2*sin(lon)) * 1e-3

    S = np.array([S_theta, S_lambda, S_r])
    # transposed transformation matrix from VES to equatorial system (ITRF)
    R_T = np.array([[cos(lat)*cos(lon), cos(lat)*sin(lon), -sin(lat)],\
                     [-sin(lon), cos(lon), 0.0],\
                     [sin(lat)*cos(lon), sin(lat)*sin(lon), cos(lat)]]).T
    # rotate to GCRS:
    dx_poltide = dot(r2000[:,:,0],dot(R_T,S))
    
    return dx_poltide
    
    
#==============================================================================
# 
#==============================================================================
def aber_source(v_GCRS, vw, K_s, r2000, earth):
    '''
    Compute source position corrected for annual and diurnal aberration.
    
    Returns apparent geodetic elevation angle and azimuth 
    measured East from North
    '''
    const = constants()
    #In order to calculate the elevation angle E it is necessary to rotate
    #the J2000.0 source unit vector to the topocentric system.
    
    # Calculation of diurnal and annual aberration 
    V_total = earth[:,1] + v_GCRS
    VdotS = dot(K_s, V_total)
    
    # Use Eq.6-66 from Zharov's book Spherical Astronomy
    abs_vec_V = norm(V_total)
    # normalise velocity vector:
    unit_V = V_total / norm(V_total)
    ndotS = dot(K_s, unit_V)
    
    rel_p = 1.0/sqrt( 1.0 - (abs_vec_V/const.C)**2 )
    w1 = ( rel_p - 1.0 ) / rel_p
    w2 = 1.0 /( 1.0 + VdotS/const.C )
    
    K_star_aber = w2 * ( K_s/rel_p + V_total/const.C + w1*ndotS*unit_V )
    # normalise aberrated vector:
    K_unit_aber = K_star_aber / norm(K_star_aber)
    
    # compute the rotation matrix which rotates from the geocentric 
    # crust fixed system to the VEN system
    
    # Rotate the aberrated vector to the crust fixed system:
    Crust_star = dot(r2000[:,:,0].T, K_unit_aber)
    ven_star = dot(vw.T, Crust_star)
    
    el = asin( ven_star[0] )
    az = atan2(ven_star[1],ven_star[2])
    if az < 0.0:
        az += 2.0*pi
    
    return az, el
    
    
#==============================================================================
# 
#==============================================================================
def thermal_def(sta, alpha, delta, elv, T, C):
    '''
    thermal_def computes delta_tau due to the 
    thermal deformation effect of the telescope

    input:
       sta - site object
       alpha, delta, elv, T - right ascention, declination, elevation, 
                               air temperature in C
       const
    output delay:
        dt_thermal
    '''
    dl = -0.12*np.pi/180.0
    phi0 = 39.06*np.pi/180.0
    
    # Antenna focus factor
    if sta.focus_type == 'FO_PRIM':
        Fa = 0.9
    else:
        Fa = 1.8

    dt_thermal = 0.0
    # Alt-azimuth
    if sta.mount_type == 'MO_AZEL':
        dt_thermal = ( sta.gamma_hf * (T - sta.T0) * (sta.hf * sin(elv)) + \
                       sta.gamma_hp * (T - sta.T0) * (sta.hp * sin(elv) + \
		             sta.AO * cos(elv) + sta.hv - Fa * sta.hs) ) / C
    # Equatorial
    elif sta.mount_type == 'MO_EQUA':
        dt_thermal = ( sta.gamma_hf * (T - sta.T0) * (sta.hf * sin(elv)) + \
                       sta.gamma_hp * (T - sta.T0) * (sta.hp * sin(elv) + \
                       sta.AO * cos(delta) + sta.hv - Fa * sta.hs) ) / C
    # XY north
    elif sta.mount_type == 'MO_XYNO':
        dt_thermal = ( sta.gamma_hf * (T - sta.T0) * (sta.hf * sin(elv)) + \
                       sta.gamma_hp * (T - sta.T0) * (sta.hp * sin(elv) + \
		 sta.AO * sqrt( 1.0 - cos(elv)*cos(elv)*cos(alpha)*cos(alpha) ) + \
		 sta.hv - Fa * sta.hs) ) / C
    # XY east
    elif sta.mount_type == 'MO_XYEA':
        dt_thermal = ( sta.gamma_hf * (T - sta.T0) * (sta.hf * sin(elv)) + \
                       sta.gamma_hp * (T - sta.T0) * (sta.hp * sin(elv) + \
		 sta.AO * sqrt( 1.0 - cos(elv)*cos(elv)*cos(alpha)*cos(alpha) ) + \
		 sta.hv - Fa * sta.hs) ) / C
    # misplaced equatorial RICHMOND
    elif sta.mount_type == 'MO_RICH':
        dt_thermal = ( sta.gamma_hf * (T - sta.T0) * (sta.hf * sin(elv)) + \
                       sta.gamma_hp * (T - sta.T0) * (sta.hp * sin(elv) + \
                       sta.AO * sqrt( 1.0 - ( sin(elv)*sin(phi0) + \
		 cos(elv)*cos(phi0)*(cos(alpha)*cos(dl) + sin(alpha)*sin(dl)) )**2 ) + \
		 sta.hv - Fa * sta.hs) ) / C    

    sta.dtau_therm = dt_thermal
    return sta
    
    
#==============================================================================
# 
#==============================================================================
def mount_tel(sta, r2000, el, az, T, P, H, const):
    '''
    MOUNT_TEL computes instrumental delay that caused by
    the axis offset. The offset changes the site position
    due to axis offset orientation.
    In subroutine the axis offset vector and its time
    derivative for different telescope mounting is calculated.

    Usially, the unit source vector corrected for aberration and
    atmospheric refraction is parallel to the symmetry axis of antenna.
    One of the ends of this axis is point of rotation of telescope;
    and axis of rotation is perpendicular to the symmetry axis.
    The axis of rotation is offset by some distance from a second
    rotation axis, all points of which are fixed relative to the Earth.
    Direction of fixed axis depends on the telescope mounting.
    '''
    RICHM = [39.06, 0.12]
    
    # Determine the unit vector representing the antenna fixed axis in a
    # topocentric VEN system    
    # Alt-azimuth
    if sta.mount_type == 'MO_AZEL':
        unit_I = np.array([1.0, 0.0, 0.0]) # xyz (VEN)
    # Equatorial
    elif sta.mount_type == 'MO_EQUA':
        unit_I = np.array([sin(sta.lat_geod), 0.0, cos(sta.lat_geod)])
    # XY north
    elif sta.mount_type == 'MO_XYNO':
        unit_I = np.array([0.0, 0.0, 1.0])
    # XY east
    elif sta.mount_type == 'MO_XYEA':
        unit_I = np.array([0.0, 1.0, 0.0])
    # misplaced equatorial RICHMOND
    elif sta.mount_type == 'MO_RICH':
        w1 = RICHM[0]*const.CDEGRAD
        w2 = RICHM[1]*const.CDEGRAD
        unit_I = np.array([sin(w1), -cos(w1)*sin(w2), cos(w1)*cos(w2)])
    else:
        print 'Unknown mount type for {:s}, guessing AltAz'.format(sta.name)
        unit_I = np.array([1.0, 0.0, 0.0]) # xyz (VEN)
    
    # Correct the aberrated topocentric source unit vector for 
    # atmospheric radio refraction. 
    
    # The zenith angle of the aberrated source:
    Z =  pi/2.0 - el
    
    Temp_K = T + 273.16
    Humid_F = H / 100.0
    Press_Hg = P * 760.0/1013.25

    # index of refraction in air
    N_air = 77.6e-6*P/Temp_K + 1.0
    # Compute atmospheric bending:
    rho = sbend(el, Temp_K, Humid_F, Press_Hg)
#    print rho
    # The apparent (aberrated+refracted) topocentric star unit
    # vector: Z_app = Z_true - rho.
    app = np.array([cos(Z - rho), \
                    sin(Z - rho)*sin(az),\
                    sin(Z - rho)*cos(az)])
    
    # Normalise aberrated+refracted vector:
    star_unit_app =  app/norm(app)    

    ## Delay computation
    # Compute topocentric axis offset vector:
    work1 = np.cross(star_unit_app, unit_I)
    vec_L = np.cross(unit_I, work1)
#    abs_vec_L = norm(vec_L)

    # Normalise topocentric axis offset vector
    unit_vec_L = vec_L/norm(vec_L)
    
    #Rotate axis offset vector to the crust fixed frame
    unit_cff = dot(sta.vw, unit_vec_L)
    # Then rotate axis offset vector to the J2000 frame:
    unit_ax2000 = dot(r2000, unit_cff)

    work2 = dot(sta.vw, star_unit_app)
    star_ab2000 = dot(r2000, work2)
    
    doff_dl = dot(star_ab2000, unit_ax2000)
    d_dax = -doff_dl/const.C*N_air
    # why minus? because the telescope 'gets closer' to the source if dl>0
    dtau_off = d_dax*sta.AO
    sta.dtau_ao = dtau_off

    return sta
    
#@numba.jit('f8(f8, f8, f8, f8)')
def sbend(El_Rad, Temp_K, Humid_F, Press_Hg):
    '''
    input:
        El_rad   -- elevation angle in radians 
        Press_Hg -- Pressure in mm of Mercury (Hg)
        Temp_K   -- Temperature in Kelvins
        Humid_F  -- relative humidity (percent)

    output   --
        bend_ang -- bending angle in radians.
    '''
    CDEGRAD = 0.017453292519943295
    CARCRAD = 4.84813681109536e-06
    
    a1 = 0.40816
    a2 = 112.30
    b1 = 0.12820
    b2 = 142.88
    c1 = 0.80000
    c2 = 99.344
    e = [46.625, 45.375, 4.1572, 1.4468,\
         0.25391, 2.2716, -1.3465, -4.3877,\
         3.1484, 4.5201, -1.8982, 0.89000 ]
    p1 = 760.0
    t1 = 273.0
    w = [ 22000.0, 17.149, 4684.1, 38.450 ]
    z1 = 91.870

    #Zenith angle in degrees
    z2 = 90.0 - El_Rad/CDEGRAD
    # Temperature in Kelvins
    t2 = Temp_K
    # Fractional humidity (0.0 -> 1.0)
    r = Humid_F
    # Pressure in mm of Hg
    p2 = Press_Hg

    # CALCULATE CORRECTIONS FOR PRES, TEMP, AND WETNESS 
    d3 = 1.0 + (z2-z1)*exp(c1*(z2-c2))
    fp = (p2/p1)*(1.0-(p2-p1)*exp(a1*(z2-a2))/d3)
    ft = (t1/t2)*(1.0-  (t2-t1)*exp(b1*(z2-b2)))
    fw = 1.0+(w[0]*r*exp((w[1]*t2-w[2])/(t2-w[3]))/(t2*p2))

    # CALCULATE OPTICAL REFRACTION 
    u = (z2-e[0])/e[1]
    x = e[10]
    for i in range(8): 
        x = e[9-i] + u*x

    # COMBINE FACTORS AND FINISH OPTICAL FACTOR
    bend_ang = ft*fp*fw*(exp(x/d3)-e[11]) 

    # BACK TO RADIANS FROM ARC SECONDS
    bend_ang = bend_ang*CARCRAD
      
    return bend_ang
    
    
#==============================================================================
# 
#==============================================================================
def tropo_wien(st, el, az, dmjd, const, do_trp_grad_calc=False):
    '''
    Tropospheric delay caclulation using VMF1 mapping functions
    '''
    # kostyl':
    # FIXME: could be out-of-bounds for uplink stations in 3-way Doppler
    # case if near midnight on the first day
    if dmjd < st.met['mjd'][0]:
        dmjd = st.met['mjd'][0] 
    # if vmf1 site data are present:
    ah = st.fMet['fAhz'](dmjd)
    aw = st.fMet['fAwz'](dmjd)

    dzh = st.fMet['fZhz'](dmjd)
    dzw = st.fMet['fZwz'](dmjd)

    lat_gcen = st.lat_gcen
    h_geod = st.h_geod

    if len(st.met['gnh'])>0:
        vmf1h, vmf1w = vmf1(ah, aw, dmjd, lat_gcen, np.pi/2.0 - el)
    else:
        vmf1h, vmf1w = vmf1_ht(ah, aw, dmjd, lat_gcen, h_geod, np.pi/2.0 - el)
    
    dtau_tropo = (dzh*vmf1h + dzw*vmf1w) / const.C
    
    # tropospheric gradients:
    if len(st.met['gnh'])>0 and do_trp_grad_calc:
        nh = st.fMet['fGnh'](dmjd)
        eh = st.fMet['fGeh'](dmjd)
        nw = st.fMet['fGnw'](dmjd)
        ew = st.fMet['fGew'](dmjd)

        # horizaontal gradient:
        grad = 1e-3 * ( 0.53*vmf1h*(nh*cos(az)+eh*sin(az))/tan(el) + \
                0.71*vmf1w*(nw*cos(az)+ew*sin(az))/tan(el) ) / const.C
        # add horizaontal gradient to the path delay:
        dtau_tropo += grad
    
    st.dtau_tropo = dtau_tropo
    
    return st
    
#@jit
def vmf1(ah, aw, dmjd, dlat, zd):
    '''
    This subroutine determines the VMF1 (Vienna Mapping Functions 1)
    for specific sites.
    Reference: Boehm, J., B. Werl, H. Schuh (2006), 
    Troposphere mapping functions for GPS and very long baseline interferometry 
    from European Centre for Medium-Range Weather Forecasts operational analysis data,
    J. Geoph. Res., Vol. 111, B02406, doi:10.1029/2005JB003629.

    Please mind that the coefficients in this paper are wrong. The corrected version of
    the paper can be found at: 
    http://ggosatm.hg.tuwien.ac.at/DOCS/PAPERS/2006Boehm_etal_VMF1.pdf


    input data
    ----------
    ah:   hydrostatic coefficient a (http://ggosatm.hg.tuwien.ac.at/DELAY/SITE/)
    aw:   wet coefficient a         (http://ggosatm.hg.tuwien.ac.at/DELAY/SITE/)  
    dmjd: modified julian date
    dlat: ellipsoidal latitude in radians
    zd:   zenith distance in radians

    output data
    -----------
    vmf1h: hydrostatic mapping function
    vmf1w: wet mapping function

    Johannes Boehm, 2005 October 2
    rev 2011 July 21: latitude -> ellipsoidal latitude

    reference day is 28 January
    this is taken from Niell (1996) to be consistent
    '''
    doy = dmjd  - 44239.0 + 1 - 28

    bh = 0.0029
    c0h = 0.062
    if (dlat<0.0):      #   ! southern hemisphere
        phh  = pi
        c11h = 0.007
        c10h = 0.002
    else:             #   ! northern hemisphere
        phh  = 0.0
        c11h = 0.005
        c10h = 0.001
          
    ch = c0h + ((cos(doy/365.25*2.0*pi + phh)+1.0)*c11h/2.0 \
             + c10h)*(1.0-cos(dlat))

    sine   = sin(pi/2.0 - zd)
    beta   = bh/( sine + ch  )
    gamma  = ah/( sine + beta)
    topcon = (1.0 + ah/(1.0 + bh/(1.0 + ch)))
    vmf1h   = topcon/(sine+gamma)

    bw = 0.00146
    cw = 0.04391
    beta   = bw/( sine + cw )
    gamma  = aw/( sine + beta)
    topcon = (1.0 + aw/(1.0 + bw/(1.0 + cw)))
    vmf1w   = topcon/(sine+gamma)

    return vmf1h, vmf1w
    
#@jit
def vmf1_ht(ah, aw, dmjd, dlat, ht, zd):
    '''
    !!! This is the version with height correction !!!  
    !!! It has to be used with the grid !!!    

    This subroutine determines the VMF1 (Vienna Mapping Functions 1)
    Reference: Boehm, J., B. Werl, H. Schuh (2006), 
    Troposphere mapping functions for GPS and very long baseline interferometry 
    from European Centre for Medium-Range Weather Forecasts operational analysis data,
    J. Geoph. Res., Vol. 111, B02406, doi:10.1029/2005JB003629.

    Please mind that the coefficients in this paper are wrong. The corrected version of
    the paper can be found at: 
    http://ggosatm.hg.tuwien.ac.at/DOCS/PAPERS/2006Boehm_etal_VMF1.pdf

    input data
    ----------
    ah:   hydrostatic coefficient a (http://ggosatm.hg.tuwien.ac.at/DELAY/GRID/)
    aw:   wet coefficient a         (http://ggosatm.hg.tuwien.ac.at/DELAY/GRID/)  
    dmjd: modified julian date
    dlat: ellipsoidal latitude in radians
    ht:   ellipsoidal height in meter
    zd:   zenith distance in radians

    output data
    -----------
    vmf1h: hydrostatic mapping function
    vmf1w: wet mapping function

    Johannes Boehm, 2005 October 2
    Rev 2011 July 21: latitude -> ellipsoidal latitude

    reference day is 28 January
    this is taken from Niell (1996) to be consistent
    '''
    doy = dmjd  - 44239.0 + 1 - 28
      
    bh = 0.0029
    c0h = 0.062
    if (dlat<0.0):      #   ! southern hemisphere
        phh  = pi
        c11h = 0.007
        c10h = 0.002
    else:               #   ! northern hemisphere
        phh  = 0.0
        c11h = 0.005
        c10h = 0.001
          
    ch = c0h + ((cos(doy/365.25*2.0*pi + phh)+1.0)*c11h/2.0 \
             + c10h)*(1.0-cos(dlat))

    sine   = sin(pi/2.0 - zd)
    beta   = bh/( sine + ch  )
    gamma  = ah/( sine + beta)
    topcon = (1.0 + ah/(1.0 + bh/(1.0 + ch)))
    vmf1h   = topcon/(sine+gamma)

    # height correction for hydrotatic part [Niell, 1996]     
    a_ht = 2.53e-5
    b_ht = 5.49e-3
    c_ht = 1.14e-3
    hs_km        = ht/1000.0
    beta         = b_ht/( sine + c_ht)
    gamma        = a_ht/( sine + beta)
    topcon       = (1.0 + a_ht/(1.0 + b_ht/(1.0 + c_ht)))
    ht_corr_coef = 1.0/sine - topcon/(sine + gamma)
    ht_corr      = ht_corr_coef * hs_km
    vmf1h        = vmf1h + ht_corr
    
    bw = 0.00146
    cw = 0.04391
    beta   = bw/( sine + cw )
    gamma  = aw/( sine + beta)
    topcon = (1.0 + aw/(1.0 + bw/(1.0 + cw)))
    vmf1w   = topcon/(sine+gamma)

    return vmf1h, vmf1w
    
    
#==============================================================================
# 
#==============================================================================
def ion_igs(sta, iono, elv, azi, jd, UT, f_0):
    '''
    Calculate ionospheric delay using IGS vertical TEC maps
    '''
    # FIXME: ionosphere must be thoroughly checked!
    UT_tec = iono.UT_tec
    fVTEC = iono.fVTEC
    
    date_tec_0 = Time(str(iono.date_tec[0]), format='iso', scale='utc')
    t_UT = Time(jd+UT, format='jd', scale='utc')
    dd = (t_UT.datetime - date_tec_0.datetime).days
#    print 'ion_igs:', t_UT.datetime, date_tec_0.datetime, dd
    
    # calculate distance to point J of ray pierce into ionosphere from site
#    H = 438.1*1e3 # m - mean height of the ionosphere above R_E
    H = 450*1e3 # m - height of ionosphere above R_E as stated in IONEX files
    R_E = 6371.0*1e3 # m - Earth's radius from the TEC map 

    alpha = 0.9782 # Schaer JPL, but this is for H = 506.7 km
    
    if 1==0:
    #    lat_geod = sta.lat_geod
        # geocentric lat is used! 
        # for reference see https://igscb.jpl.nasa.gov/igscb/data/format/ionex1.pdf
        lat_gcen = sta.lat_gcen
    #    lon_gcen = sta.lon_gcen
        h_geod = sta.h_geod
        
        #WGS84 Ellipsoid
        a = 6378137.0 # m
        f = 1.0/298.2572235630
        b = a*(1.0-f) # m
        ec = sqrt((a**2-b**2)/(a**2))
        
        R_oscul = a*sqrt(1.0-ec**2)/(1.0-(ec*sin(lat_gcen))**2) # m
        
        source_vec = np.array([sin(pi/2.0-elv)*cos(azi),\
                               sin(pi/2.0-elv)*sin(azi),\
                               cos(pi/2.0-elv) ])
        
        # slanted distance btw the ground and the iono layer
        ds = (R_oscul+h_geod)*sin(-elv) + \
             0.5 * sqrt( (2.0*(R_oscul+h_geod)*sin(-elv))**2 - \
             4.0*((R_oscul+h_geod)**2 - (R_E+H)**2) )
        # cart crds of the starting point
    #    rpt = [R_oscul+h_geod, lat_gcen, lon_gcen]
    #    r0 = sph2cart(rpt)
        r0 = sta.r_GTRS
        # cart crds of the ionospheric pierce point
        r1 = r0 + ds*source_vec

    # follow Petrov! vtd manual-2014 p.32 (error!), thus vtd_iono_delay.f
    beta = np.arcsin(cos(elv)/(1 + H/R_E))
    D = R_E*sqrt(1 + (1+H/R_E)**2 - 2*(1+H/R_E)*sin(elv+beta))
    r0 = sta.r_GTRS
    source_vec = np.array([sin(pi/2.0-elv)*cos(azi),\
                           sin(pi/2.0-elv)*sin(azi),\
                           cos(pi/2.0-elv) ])
    r1 = r0 + D*source_vec
    
    # lat/long of the pierce point
    rlalo = cart2sph(r1)
    lon = 180.0*rlalo[2]/pi
    lat = 180.0*rlalo[1]/pi
#    print lat, lon
    
    # yet another way of computing lon/lat
    #http://www.tudelft.nl/fileadmin/Faculteit/CiTG/
    #Over_de_faculteit/Afdelingen/Afdeling_Geoscience_and_Remote_Sensing/
    #MSc_theses/wienia08_msc.pdf p.49:
    a = 6378137.0 # m
    f = 1.0/298.2572235630
    b = a*(1.0-f) # m
    ec = sqrt((a**2-b**2)/(a**2))
    
    R_oscul = a*sqrt(1.0-ec**2)/(1.0-(ec*sin(sta.lat_geod))**2) # m
    z = pi/2 - elv
    z_prime = np.arcsin(R_oscul/(R_E+H)*sin(z))
    dz = z - z_prime
#    print 'dz', dz
    lon = (sta.lon_gcen + np.arcsin(sin(dz)*sin(azi)/sin(sta.lat_gcen)))*180/pi
#    print sin(sta.lat_gcen)*cos(dz)+cos(sta.lat_gcen)*sin(dz)*cos(azi)
    lat = (np.arcsin(sin(sta.lat_gcen)*cos(dz)+cos(sta.lat_gcen)*sin(dz)*cos(azi)))*180/pi
#    print lat, lon
#    raw_input()
#    print UT
    
#    print sta.lat_gcen*180/pi, sta.lon_gcen*180/pi
#    print lat, lon
#    print azi*180/pi, elv*180/pi
    
    # find closest epoch in TEC data:
    # easy case - UT is in UT_tec
    n0 = np.searchsorted(UT_tec, UT+dd)
    if UT+dd in UT_tec:
        # only need to interpolate TECz to the Lat/Lon of the pierce point
        TEC_z = float(fVTEC[n0](lon, lat))
    else:
        # else take 4 consecutive epochs and interpolate:
        N_epoch = len(fVTEC)
        if n0==1 or n0==0:
            nl = 0; nr = 4
        elif n0==N_epoch-2 or n0==N_epoch-1:
            nl = N_epoch-4; nr = N_epoch
        else:
            nl = n0-2; nr = n0+2;
        TEC_z = []
        for nn in range(nl,nr):
            TEC_z.append(float(fVTEC[nn](lon, lat)))
        # interpolate zenith TECz to the epoch of observation:
#        fTEC = sp.interpolate.interp1d(UT_tec[nl:nr], TEC_z, kind='cubic')
#        TEC_z = fTEC(UT)
#        print UT, UT_tec[nl:nr]
        TEC_z = np.interp(UT+dd, UT_tec[nl:nr], TEC_z)
        # TODO: fix at some point
        # note: when calculating tau_iono for an uplink station when t_obs 
        # starts at 00h, this will yield TEC_z[0], as it doesn't 
        # do extrapolation. No big deal...
    
    # calculate slanted TEC
    if 1==0:
        TEC = TEC_z / cos(asin( (R_oscul+h_geod)*sin(alpha*(pi/2.0-elv))/(R_E+H) ))
    # follow Petrov:
    TEC = TEC_z / cos(beta)
    TEC_tecu = 0.1*TEC # in TEC units

#    print TEC_tecu
#    raw_input()
    
    # calculate ionspheric delay for the source
    delay_ion = 5.308018e10*TEC_tecu/(4.0*pi**2*f_0*f_0)

    sta.dtau_iono = delay_ion
    
    return sta
    

#==============================================================================
# 
#==============================================================================
def iau_tdbtcb( TDB1, TDB2 ):
    '''
    +
      - - - - - - - - - - -
       i a u _ T D B T C B
      - - - - - - - - - - -
    
      Time scale transformation:  Barycentric Dynamical Time, TDB, to
      Barycentric Coordinate Time, TCB.
    
      This routine is part of the International Astronomical Union's
      SOFA (Standards of Fundamental Astronomy) software collection.
    
      Status:  canonical.
    
      Given:
         TDB1,TDB2    d      TDB as a 2-part Julian Date
    
      Returned:
         TCB1,TCB2    d      TCB as a 2-part Julian Date
    
      Notes:
    
      1  TDB1+TDB2 is Julian Date, apportioned in any convenient way
         between the two arguments, for example where TDB1 is the Julian
         Day Number and TDB2 is the fraction of a day.  The returned
         TCB1,TCB2 follow suit.
    
      2  The 2006 IAU General Assembly introduced a conventional linear
         transformation between TDB and TCB.  This transformation
         compensates for the drift between TCB and terrestrial time TT,
         and keeps TDB approximately centered on TT.  Because the
         relationship between TT and TCB depends on the adopted solar
         system ephemeris, the degree of alignment between TDB and TT over
         long intervals will vary according to which ephemeris is used.
         Former definitions of TDB attempted to avoid this problem by
         stipulating that TDB and TT should differ only by periodic
         effects.  This is a good description of the nature of the
         relationship but eluded precise mathematical formulation.  The
         conventional linear relationship adopted in 2006 sidestepped
         these difficulties whilst delivering a TDB that in practice was
         consistent with values before that date.
    
      3  TDB is essentially the same as Teph, the time argument for the
         JPL solar system ephemerides.
    
      Reference:
    
         IAU 2006 Resolution B3
    
      This revision:  2010 September 10
    
      SOFA release 2010-12-01
    
      Copyright (C) 2010 IAU SOFA Board.  See notes at end.
    
    -----------------------------------------------------------------------
    '''
    # 1977 Jan 1.0 TAI = 1977/1/1 00:00:32.184 TCB, as two-part JD
    T77TD = 2443144.0
    T77TF = 0.5003725

    # L_B, and TDB0 (d)
    ELB = 1.550519768e-8
    TDB0 = -6.55e-5/86400.0

    # TDB to TCB rate
    ELBB = ELB/(1.0-ELB)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # Result, preserving date format but safeguarding precision.
    if ( TDB1>TDB2 ):
        D = T77TD - TDB1
        F = TDB2 - TDB0
        TCB1 = TDB1
        TCB2 = F - ( D - ( F - T77TF ) ) * ELBB
    else:
        D = T77TD - TDB2
        F = TDB1 - TDB0
        TCB1 = F - (D-(F-T77TF))*ELBB
        TCB2 = TDB2
    
    return TCB1, TCB2

@memoize
def iau_tcbtdb(TCB1, TCB2):
    '''
    +
      - - - - - - - - - - -
       i a u _ T C B T D B
      - - - - - - - - - - -
    
      Time scale transformation:  Barycentric Coordinate Time, TCB, to
      Barycentric Dynamical Time, TDB.
    
      This routine is part of the International Astronomical Union's
      SOFA (Standards of Fundamental Astronomy) software collection.
    
      Status:  canonical.
    
      Given:
         TCB1,TCB2    d      TCB as a 2-part Julian Date
    
      Returned:
         TDB1,TDB2    d      TDB as a 2-part Julian Date
         J            i      status:  0 = OK
    
      Notes:
    
      1  TCB1+TCB2 is Julian Date, apportioned in any convenient way
         between the two arguments, for example where TCB1 is the Julian
         Day Number and TCB2 is the fraction of a day.  The returned
         TDB1,TDB2 follow suit.
    
      2  The 2006 IAU General Assembly introduced a conventional linear
         transformation between TDB and TCB.  This transformation
         compensates for the drift between TCB and terrestrial time TT,
         and keeps TDB approximately centered on TT.  Because the
         relationship between TT and TCB depends on the adopted solar
         system ephemeris, the degree of alignment between TDB and TT over
         long intervals will vary according to which ephemeris is used.
         Former definitions of TDB attempted to avoid this problem by
         stipulating that TDB and TT should differ only by periodic
         effects.  This is a good description of the nature of the
         relationship but eluded precise mathematical formulation.  The
         conventional linear relationship adopted in 2006 sidestepped
         these difficulties whilst delivering a TDB that in practice was
         consistent with values before that date.
    
      3  TDB is essentially the same as Teph, the time argument for the
         JPL solar system ephemerides.
    
      Reference:
    
         IAU 2006 Resolution B3
    
      This revision:  2010 May 13
    
      SOFA release 2010-12-01
    
      Copyright (C) 2010 IAU SOFA Board.  See notes at end.
    
    -----------------------------------------------------------------------
    '''
    T77TD = 2443144.0
    T77TF = 0.5003725

    ELB = 1.550519768e-8
    TDB0 = -6.55e-5/86400.0

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # Result, safeguarding precision.
    if ( TCB1>TCB2 ):
        D = TCB1 - T77TD
        TDB1 = TCB1
        TDB2 = TCB2 + TDB0 - ( D + ( TCB2-T77TF ) ) * ELB
    else:
        D = TCB2 - T77TD
        TDB1 = TCB1 + TDB0 - ( D + ( TCB1-T77TF ) ) * ELB
        TDB2 = TCB2

    return TDB1, TDB2
    
    
#==============================================================================
# 
#==============================================================================
def delay_iers(tjd, CT, r_1, r_2, v_2, earth, sun, K_s, jpl_eph, \
               GM, TDB_TCB, L_C, C):
    '''
    VLBI delay calculation following IERS Conventions 2010
    '''    
    # G*masses in TCB-frame! (see constants class)

    gamma_PPN = 1
    
#    b = sta[1].r_GCRS - sta[0].r_GCRS
    b = r_2 - r_1
    
    # follow Kopeikin's notations
#    _, t = iau_tdbtcb(tjd, CT) # reception time in TCB
    
    # Find potential U:
    U = GM[10]/norm(sun[:,0]-earth[:,0])
    
    # BCRS radius vectors of the reception sites at t:
    x_1 = earth[:,0] + (1.0 - U/(C**2) - L_C)*r_1 - \
          dot(earth[:,1],r_1)*earth[:,1] / (2.0*C**2)
    x_2 = earth[:,0] + (1.0 - U/(C**2) - L_C)*r_2 - \
          dot(earth[:,1],r_2)*earth[:,1] / (2.0*C**2)

    ''' This is more precise, but it's irrelevant at a ps level
    # Kopeikin's algorithm for finding gravitational delay
    precision = 1e-16 # precision of the Newton-Raphson iterations
    n_max = 3 # maximum number of the Newton-Raphson iterations
    ## derive retarded times s
    # >> Gravitating bodies and gravitational delay
    # (0=Mercury,1=Venus,2=Earth,3=Mars,4=Jupiter,5=Saturn,
    #  6=Uranus,7=Neptune,8=Pluto,9=Moon,10=Sun)
#    tic = _time()
    s_A_tmp = np.zeros(11)
    s_A = np.zeros(11)
    x_A = np.zeros((3,11))
    v_A = np.zeros((3,11))
    r_A = np.zeros((3,11))
    r_2_A = np.zeros((3,11))

    T_g = 0.0
    for kk in range(11):
        if ( (norm(r_1)==0.0 or norm(r_2)==0.0) and kk==2):
            continue
        nn = 1
        s_A[kk] = t
        while (abs(s_A[kk] - s_A_tmp[kk]) > precision) and (nn < n_max):
            s_A_tmp[kk] = s_A[kk]
            _, s_A_TDB = iau_tcbtdb(tjd, s_A[kk])
            rrd = pleph(tjd+s_A_TDB, kk+1, 12, jpl_eph)
            x_A[:,kk] = rrd[0:3]*1e3 * TDB_TCB
            v_A[:,kk] = rrd[3:]*1e3
            r_A[:,kk] = x_1 - x_A[:,kk]
            s_A[kk] = (s_A[kk]*86400.0 - \
                      (s_A[kk]*86400.0 - t*86400.0 + \
                        norm(r_A[:,kk])/C) / \
                      (1.0 - dot(r_A[:,kk], v_A[:,kk]) / \
                       (norm(r_A[:,kk])*C))) / 86400.0;
            nn += 1
        # calc position and velocity vectors at the found s
        _, s_A_TDB = iau_tcbtdb(tjd, s_A[kk])
        rrd = pleph(tjd+s_A_TDB, kk+1, 12, jpl_eph)
        x_A[:,kk] = rrd[0:3]*1e3 * TDB_TCB
        v_A[:,kk] = rrd[3:]*1e3
        r_A[:,kk] = x_1 - x_A[:,kk]
        r_2_A[:,kk] = x_2 - dot(K_s,b)*earth[:,1]/C - x_A[:,kk]
        # gravitational delay
        T_g = T_g + (2.0*GM[kk]/(C**3))*\
              log( (norm(r_A[:,kk])+dot(K_s,r_A[:,kk])) / \
                   (norm(r_2_A[:,kk])+dot(K_s,r_2_A[:,kk])) )
        # the difference between T_g and T_g_IERS is negligible at ps level
#    print 'aa', _time()-tic
#    print 'T_g = ', T_g
    '''
    
    # IERS algorithm for finding gravitational delay
#    tic = _time()
    T_g = 0.0
    for kk in range(11):
        if ( (norm(r_1)==0.0 or norm(r_2)==0.0) and kk==2):
            continue
        # X_J at t_1:
        rrd = pleph(tjd+CT, kk+1, 12, jpl_eph)
        x_J = rrd[0:3]*1e3
#        v_J = rrd[3:]*1e3
        t_J = min(CT, CT - dot(K_s, x_J-x_1)/C/86400.0)
        rrd = pleph(tjd+t_J, kk+1, 12, jpl_eph)
        x_J = rrd[0:3]*1e3
        r_1_J = x_1 - x_J
        r_2_J = x_2 - dot(K_s, b)*earth[:,1]/C - x_J
        T_g = T_g + (2.0*GM[kk]/(C**3))*\
              log( (norm(r_1_J)+dot(K_s, r_1_J)) / \
                   (norm(r_2_J)+dot(K_s, r_2_J)) )
#    print 'bb', _time()-tic
#    print 'T_g = ', T_g

    # Lorenz-transform from BCRS to GCRS
    # full geometric delay in the TT-frame
    dtau = ( T_g - (dot(K_s, b)/C) * \
                    ( 1.0 - (1.0+gamma_PPN)*U/C**2 - \
                      0.5*(norm(earth[:,1])/C)**2 -\
                      dot(earth[:,1], v_2)/C**2 ) - \
                   (dot(earth[:,1], b)/C**2) * \
                    (1.0 + 0.5*dot(K_s, earth[:,1])/C) \
            ) / (1.0 + dot(K_s, earth[:,1] + v_2)/C)
    
    return dtau
    
    
#==============================================================================
# 
#==============================================================================
def delay_ra(tjd, CT, UTC, r_1, r_2, v_2, a_2, r_3, earth, sun,\
             K_s, jpl_eph, GM, TDB_TCB, L_C, C, AE, uv=False):
    '''
    VLBI delay calculation following Vlasov, Zharov, Sazhin (2012)
    '''    
    # G*masses in TCB-frame! (see constants class)
    
    b = r_2 - r_1
    
    # Find potential U:
    U = GM[10]/norm(sun[:,0]-earth[:,0])

#    earth = earth*TDB_TCB
#    sun = sun*TDB_TCB
    
    # BCRS radius vectors of the reception sites at t:
    x_1 = earth[:,0] + (1.0 - U/(C**2) - L_C)*r_1 - \
          dot(earth[:,1],r_1)*earth[:,1] / (2.0*C**2)
    x_2 = earth[:,0] + (1.0 - U/(C**2) - L_C)*r_2 - \
          dot(earth[:,1],r_2)*earth[:,1] / (2.0*C**2)

    
    # IERS algorithm for finding gravitational delay
    T_g = 0.0
    for kk in range(11):
        if ( (norm(r_1)==0.0 or norm(r_2)==0.0) and kk==2):
            continue
        # X_J at t_1:
        rrd = pleph(tjd+CT, kk+1, 12, jpl_eph)
        x_J = rrd[0:3]*1e3
#        v_J = rrd[3:]*1e3
        t_J = min(CT, CT - dot(K_s, x_J-x_1)/C/86400.0)
        rrd = pleph(tjd+t_J, kk+1, 12, jpl_eph)
        x_J = rrd[0:3]*1e3
        r_1_J = x_1 - x_J
        r_2_J = x_2 - dot(K_s, b)*earth[:,1]/C - x_J
        T_g = T_g + (2.0*GM[kk]/(C**3))*\
              log( (norm(r_1_J)+dot(K_s, r_1_J)) / \
                   (norm(r_2_J)+dot(K_s, r_2_J)) )
        if kk==9:
            r_moon = x_J - earth[:,0]
        if kk==10:
            r_sun = x_J - earth[:,0]
    
    # Lorenz-transform from BCRS to GCRS
    # full geometric delay in the TT-frame
    dtau = ( T_g - (dot(K_s, b)/C) * \
             ( 1.0 - 0.5*(norm(earth[:,1])/C)**2  - \
             dot(v_2, earth[:,1])/C**2 - 2.0*U/C**2 - \
             dot(earth[:,2], r_2)/C**2 + \
             dot(K_s, b)*(dot(K_s,a_2)+dot(K_s, earth[:,2]))/(2.0*C**2) ) - \
             dot(b, earth[:,1])*(1.0+dot(K_s, earth[:,1])/(2.0*C))/(C**2) + \
             dot(K_s, r_2)*dot(r_2, earth[:,2])/C**3 - \
             dot(K_s, earth[:,2])*(norm(r_2)**2)/(2.0*C**3)  - \
             dot(K_s, r_1)*dot(r_1, earth[:,2])/C**3 + \
             dot(K_s, earth[:,2])*(norm(r_1)**2)/(2.0*C**3) ) \
            / (1.0 + dot(K_s, earth[:,1] + v_2)/C)
    
    if uv:
        return dtau
    
    ''' calculate dtau_spacecraft/dt (proper time by geocentric coordinate time)
        to precision roughly 10^-15 '''
    # dtau/dTCB:
    J_2 = 1.0826e-3
    # -(1+L_G)
    dtau_dt_min1 = - ( 0.5*norm(v_2)**2 + GM[2]/norm(r_2) + \
                      GM[2]*AE**2*J_2 * \
                         (1.0-3.0*(r_2[2]/norm(r_2))**2)/(2.0*norm(r_2)**3) + \
                      GM[10]*(3.0*(dot(r_2,r_sun))**2/(norm(r_sun)**2) - \
                         norm(r_2)**2)/(2.0*norm(r_sun)**3) + \
                      GM[9]/norm(r_moon-r_2) - GM[9]/norm(r_moon) + \
                         GM[9]*dot(r_2,r_moon)/(norm(r_moon)**3) ) / C**2
    
    ''' calculate downlink light-time from RA to Pushchino at CT '''
    # Get Solar system bodies r, v at CT wrt Earth
    r_B = []
    v_B = []
    for kk in range(11):
        rrd = pleph(tjd+CT, kk+1, 3, jpl_eph)
        r_B.append(rrd[0:3]*1e3)
        v_B.append(rrd[3:]*1e3)
    
    precision = 1e-13
    n_max = 2
    
    # initial approximation:
    nn = 0
    lt_tmp = 0.0

    lt = UTC - (norm(r_3 - r_2)/C)/86400.0

    while (abs(lt - lt_tmp) > precision) and (nn < n_max):
        lt_tmp = lt
        dlt = (UTC - lt)*86400.0
        r = r_2 - dlt*v_2 + (dlt**2)*a_2/2.0
        v = v_2 - dlt*a_2
        # vectors needed for RLT calculation
        r_23 = r_3 - r
        # >> SS bodies
        RLT = 0.0
        for ii, (rb, vb) in enumerate(zip(r_B, v_B)):
            if ii==3: continue
            r_2_B  = r - (rb - dlt*vb)
            r_3_B  = r_3 - rb
            r_23_B = r_3_B - r_2_B
            RLT = RLT + (2.0*GM[ii]/C**3) * \
                  log(  ( norm(r_2_B) + norm(r_3_B) + norm(r_23_B) + \
                          2.0*GM[ii]/C**2 ) / \
                        ( norm(r_2_B) + norm(r_3_B) - norm(r_23_B) + \
                          2.0*GM[ii]/C**2 ) )
        p_dot = np.dot(r_23, v) / norm(r_23)
        dlt = (dlt - norm(r_23)/C - RLT) / (1.0 - p_dot/C)
        lt = lt + dlt/86400.0
        nn += 1
    
    lt_downlink = (UTC - lt)*86400.0

    return dtau, dtau_dt_min1, lt_downlink
    
    

#==============================================================================
# 
#==============================================================================
#@numba.jit
def delay_moyer(tjd, t_1, dd, r_1, r_2, v_2, a_2, state_ss, tdb, bcrs,
                GM, TDB_TCB, L_C, C, inp):
    '''
    VLBI delay calculation following Moyer/Duev
    For reference see Duev et al. 2012 A&A 43
    
    tjd - Julian Date of observation
    t_1 - epoch in TDB of observation, decimal days [0,1] * *86400
    dd - number of days since the start epoch of the ephemeris *86400
    r_1 - radius-vector of station 1 in GCRS (usually, Geocenter)
    r_2, v_2, a_2 - radius-vector, velocity, acceleration of station 2 in GCRS
    state_ss  - Solar system bodies r, v (and a for Earth) at t_1 wrt SSBC
    
    tdb - "time coordinate" of the S/C ephemeris, decimal days *86400
    bcrs - state of the S/C
    
    GM - G*M for Solar System bodies
    TDB_TCB - obvious
    L_C = 1 -d(TCG)/d(TCB)
    C - the speed of light
    '''
    earth = state_ss[2]
    sun = state_ss[-1]
    
    # Find potential U:
    U = GM[10]/norm(sun[:,0]-earth[:,0])

    # BCRS radius vectors of the reception sites at t:
    R_1 = earth[:,0] + (1.0 - U/(C**2) - L_C)*r_1 - \
          dot(earth[:,1],r_1)*earth[:,1] / (2.0*C**2)
    R_2 = earth[:,0] + (1.0 - U/(C**2) - L_C)*r_2 - \
          dot(earth[:,1],r_2)*earth[:,1] / (2.0*C**2)
    V_2 = earth[:,1] + \
            (1.0 - 2.0*U/C**2 - 0.5*(norm(earth[:,1])/C)**2 - \
             dot(earth[:,1], v_2)/C**2) * v_2 - \
            0.5*dot(earth[:,1], v_2)*earth[:,1]/C**2
    A_2 = earth[:,2] + \
            (1.0 - 3.0*U/C**2 - (norm(earth[:,1])/C)**2 + L_C - \
             2.0*dot(earth[:,1], v_2)/C**2) * a_2 - \
            0.5*dot(earth[:,1], a_2)*(earth[:,1] + 2.0*v_2)/C**2

    ''' calculate 1st downleg light-time from S/C to Receiver 1
        to find signal transmission t_0 time given the reception time t_1
    '''    
    precision = 1e-15
    n_max = 3
    lag_order = 5
    
    # initial approximation:
    nn = 0
    lt_01_tmp = 0.0
    
    # s/c:

    x, _ = lagint(lag_order, tdb, bcrs[:,6], dd+t_1)
    y, _ = lagint(lag_order, tdb, bcrs[:,7], dd+t_1)
    z, _ = lagint(lag_order, tdb, bcrs[:,8], dd+t_1)
    vx, _ = lagint(lag_order, tdb, bcrs[:,9], dd+t_1)
    vy, _ = lagint(lag_order, tdb, bcrs[:,10], dd+t_1)
    vz, _ = lagint(lag_order, tdb, bcrs[:,11], dd+t_1)
    R_0 = np.hstack((x,y,z))
    V_0 = np.hstack((vx,vy,vz))

#    print 't_1 = {:.18f}'.format(t_1)
    lt_01 = (norm(R_1 - R_0)/C)
#    print '{:.18f}'.format(lt_01)

    while (abs(lt_01 - lt_01_tmp) > precision) and (nn < n_max):
        lt_01_tmp = lt_01
        t_0 = t_1 - lt_01
        
        x, _ = lagint(lag_order, tdb, bcrs[:,6], dd+t_0)
        y, _ = lagint(lag_order, tdb, bcrs[:,7], dd+t_0)
        z, _ = lagint(lag_order, tdb, bcrs[:,8], dd+t_0)
        vx, _ = lagint(lag_order, tdb, bcrs[:,9], dd+t_0)
        vy, _ = lagint(lag_order, tdb, bcrs[:,10], dd+t_0)
        vz, _ = lagint(lag_order, tdb, bcrs[:,11], dd+t_0)
        R_0 = np.hstack((x,y,z))
        V_0 = np.hstack((vx,vy,vz))
        
        # vector needed for RLT calculation
        R_01 = R_1 - R_0

        # >> SS bodies
        RLT = 0.0
        
        ''' BCRS state vectors of celestial bodies at t_0, [m, m/s]: '''
        ## Earth:
        JD = tjd
        rrd = pleph(JD+t_0/86400.0, 3, 12, inp['jpl_eph'])
        earth = np.reshape(np.asarray(rrd), (3,2), 'F') * 1e3
        ## Sun:
        rrd = pleph(JD+t_0/86400.0, 11, 12, inp['jpl_eph'])
        sun = np.reshape(np.asarray(rrd), (3,2), 'F') * 1e3
        ## Moon:
        rrd = pleph(JD+t_0/86400.0, 10, 12, inp['jpl_eph'])
        moon = np.reshape(np.asarray(rrd), (3,2), 'F') * 1e3
    
        state_ss_t0 = []
        for jj in (1,2,4,5,6,7,8,9):
            rrd = pleph(JD+t_0/86400.0, jj, 12, inp['jpl_eph'])
            state_ss_t0.append(np.reshape(np.asarray(rrd), (3,2), 'F') * 1e3)
        state_ss_t0.insert(2,earth)
        state_ss_t0.append(moon)
        state_ss_t0.append(sun)
        
#        for ii, state in enumerate(state_ss):
        for ii, (state_t1, state_t0) in enumerate(zip(state_ss, state_ss_t0)):
            if not (ii==2 and norm(r_1)<1e-3):
                rb_t1 = state_t1[:,0]
                rb_t0 = state_t0[:,0]
    #            vb = state[:,1]
    #            R_0_B  = R_0 - (rb - lt_01*vb)
                R_0_B  = R_0 - rb_t0
                R_1_B  = R_1 - rb_t1
                R_01_B = R_1_B - R_0_B
    
                RLT += (2.0*GM[ii]/C**3) * \
                      log(  ( norm(R_0_B) + norm(R_1_B) + norm(R_01_B) + \
                              2.0*GM[ii]/C**2 ) / \
                            ( norm(R_0_B) + norm(R_1_B) - norm(R_01_B) + \
                              2.0*GM[ii]/C**2 ) )

        lt_01 = lt_01 - (lt_01 - norm(R_01)/C - RLT) / \
                    ( 1.0 - dot(R_01, V_0)/(C*norm(R_01)) )

        t_0 = t_1 - lt_01
#        print 't_0 = {:.18f}'.format(t_0)
        nn += 1
    
    x, _ = lagint(lag_order, tdb, bcrs[:,6], dd+t_0)
    y, _ = lagint(lag_order, tdb, bcrs[:,7], dd+t_0)
    z, _ = lagint(lag_order, tdb, bcrs[:,8], dd+t_0)
    vx, _ = lagint(lag_order, tdb, bcrs[:,9], dd+t_0)
    vy, _ = lagint(lag_order, tdb, bcrs[:,10], dd+t_0)
    vz, _ = lagint(lag_order, tdb, bcrs[:,11], dd+t_0)
    R_0 = np.hstack((x,y,z))
    V_0 = np.hstack((vx,vy,vz))
#    print R_0
#    raw_input()
    ''' BCRS state vectors of celestial bodies at t_0, [m, m/s]: '''
    ## Earth:
    JD = tjd
    rrd = pleph(JD+t_0/86400.0, 3, 12, inp['jpl_eph'])
    earth = np.reshape(np.asarray(rrd), (3,2), 'F') * 1e3
    ## Sun:
    rrd = pleph(JD+t_0/86400.0, 11, 12, inp['jpl_eph'])
    sun = np.reshape(np.asarray(rrd), (3,2), 'F') * 1e3
    ## Moon:
    rrd = pleph(JD+t_0/86400.0, 10, 12, inp['jpl_eph'])
    moon = np.reshape(np.asarray(rrd), (3,2), 'F') * 1e3

    state_ss_t0 = []
    for jj in (1,2,4,5,6,7,8,9):
        rrd = pleph(JD+t_0/86400.0, jj, 12, inp['jpl_eph'])
        state_ss_t0.append(np.reshape(np.asarray(rrd), (3,2), 'F') * 1e3)
    state_ss_t0.insert(2,earth)
    state_ss_t0.append(moon)
    state_ss_t0.append(sun)
    
    ''' calculate 2st downleg light-time from S/C to Receiver 2
        to find signal reception time t_2 given the transmission time t_0
    '''    
    # initial approximation:
    nn = 0
    lt_02_tmp = 0.0
    # receiver #2:
    R_2_T1 = R_2
    V_2_T1 = V_2
    A_2_T1 = A_2

    lt_02 = (norm(R_2_T1 - R_0)/C)

    while (abs(lt_02 - lt_02_tmp) > precision) and (nn < n_max):
        lt_02_tmp = lt_02
        
        # vector needed for RLT calculation
        t_2 = t_0 + lt_02
        dt_12 = t_2 - t_1
#        print t_2, t_1
        R_2 = R_2_T1 + V_2_T1*dt_12 + A_2_T1*dt_12**2/2
        V_2 = V_2_T1 + A_2_T1*dt_12
        R_02 = R_2 - R_0

        # >> SS bodies
        RLT = 0.0
        for ii, (state_t1, state_t0) in enumerate(zip(state_ss, state_ss_t0)):
            if not (ii==2 and norm(r_2)<1e-3):
                rb_t1 = state_t1[:,0]
                vb_t1 = state_t1[:,1]
                rb_t0 = state_t0[:,0]
                R_0_B  = R_0 - rb_t0
                R_2_B  = R_2 - (rb_t1 + dt_12*vb_t1)
                R_02_B = R_2_B - R_0_B
                RLT = RLT + (2.0*GM[ii]/C**3) * \
                      log(  ( norm(R_0_B) + norm(R_2_B) + norm(R_02_B) + \
                              2.0*GM[ii]/C**2 ) / \
                            ( norm(R_0_B) + norm(R_2_B) - norm(R_02_B) + \
                              2.0*GM[ii]/C**2 ) )
#        print 'RLT = {:.18f}'.format(RLT)
        lt_02 = lt_02 - (lt_02 - norm(R_02)/C - RLT) / \
                    ( 1.0 - dot(R_02, V_2)/(C*norm(R_02)) )
        
        nn += 1
#        print 't_2 = {:.18f}'.format(t_2)
#    print 'f(x)=0?: {:.18f}'.format(t_2 - t_0 - norm(R_02)/C - RLT)
#    raw_input()
    ''' delay TDB -> TT '''
#    T2_T1 = t_2 - t_1
    T2_T1 = lt_02 - lt_01
    # Fukushima
    t2_t1 = T2_T1 * \
            ( 1 - (norm(earth[:,1])**2/2.0 + \
                   dot(earth[:,1],v_2) + U)/C**2) / (1.0 - L_C) - \
                 dot(earth[:,1],r_2-r_1)/C**2
    print 't2_t1 = {:.18f}'.format(t2_t1)
    # Myself
    t2_t1 = ( T2_T1 * \
             ( 1 - ((norm(earth[:,1])**2)/2.0 + U)/C**2) / (1.0 - L_C) - \
             dot(earth[:,1],r_2-r_1)/C**2 ) / \
             ( 1 + dot(earth[:,1],v_2) / C**2 )
#    print T2_T1
    print 't2-t1 = {:.18f}'.format(t2_t1)
#    raw_input()
    
    return t2_t1
    
    
#==============================================================================
# Station short names
#==============================================================================
def shname(staz, shnames_cat, shnames_cat_igs=None):
    '''
    Station short names
    '''
    staz_short = []
    
    with open(shnames_cat,'r') as cat:
        cat_lines = cat.readlines()
    cat_lines = [l for l in cat_lines if l[0]!='#']
    # IGS catalogue:
    if shnames_cat_igs is not None:
        with open(shnames_cat_igs,'r') as cat_igs:
            cat_igs_lines = cat_igs.readlines()
            cat_igs_lines = [l for l in cat_igs_lines if l[0]!='*']
    
    # if single string was input, convert it to a list
    if type(staz)==type('A'):
        staz = [staz]
        
    # find short name for each station:
    for sta in staz:
        matching = [line for line in cat_lines \
                    if line[:8].strip() == sta]
        if len(matching)==0:
            if shnames_cat_igs is None:
                raise Exception('Short name for station '+sta+' not found.\n'+\
                                'Check catalogue '+shnames_cat+'.')
            print 'Short name for station '+sta+' not found in '+shnames_cat
            print 'Trying IGS catalogue '+shnames_cat_igs
#            print [line[4:12].strip() for line in cat_igs_lines]
            matching = [line for line in cat_igs_lines \
                    if line[4:12].strip() == sta]            
            matching = matching[-1].split()[0].strip()
            print matching
            if len(matching)==0:
                raise Exception('Short name for station '+sta+' not found.\n'+\
                  'Check catalogues '+shnames_cat+', '+shnames_cat_igs+'.')
            staz_short.append(matching)
        else:
            matching = matching[-1].strip() # skip comments, which will show up first here
            staz_short.append(matching[-2:])
        
    return staz_short
    
#==============================================================================
# Station long names
#==============================================================================
def loname(staz, shnames_cat, shnames_cat_igs=None):
    '''
    Station short names
    '''
    staz_long = []
    
    with open(shnames_cat,'r') as cat:
        cat_lines = cat.readlines()
    cat_lines = [l for l in cat_lines if l[0]!='#']
    # IGS catalogue:
    if shnames_cat_igs is not None:
        with open(shnames_cat_igs,'r') as cat_igs:
            cat_igs_lines = cat_igs.readlines()
        cat_igs_lines = [l for l in cat_igs_lines if l[0]!='*']
    
    # if single string was input, convert it to a list
    if type(staz)==type('A'):
        staz = [staz]
        
    # find long name for each station:
    for sta in staz:
        matching = [line for line in cat_lines \
                    if line[11:].strip() == sta.lower()]
        if len(matching)==0:
            if shnames_cat_igs is None:
                raise Exception('Long name for station '+sta+' not found.\n'+\
                                'Check catalogue '+shnames_cat+'.')
            print 'Long name for station '+sta+' not found in '+shnames_cat
            print 'Trying IGS catalogue '+shnames_cat_igs
#            print [line[4:12].strip() for line in cat_igs_lines if line[0]!='*']
            matching = [line for line in cat_igs_lines \
                    if line[:4].strip().lower() == sta.lower()]
            matching = matching[-1].split()[1].strip()
            print matching
            if len(matching)==0:
                raise Exception('Long name for station '+sta+' not found.\n'+\
                  'Check catalogues '+shnames_cat+', '+shnames_cat_igs+'.')
            staz_long.append(matching)
        else:
            matching = matching[-1].strip() # skip comments, which will show up first here
            staz_long.append(matching.split()[0])
            
    return staz_long
    
#==============================================================================
# 
#==============================================================================
def delay_nf_moyer(jd, t_1_days, dd_days, state_ss_t1, tdb, bcrs, const,\
                    sta1, sta2, inp, t_1_UTC):
    '''
    NF delay calculation
    '''
    debug = False
    
    GM = const.GM
    C = const.C
    L_C = const.L_C
    r_1 = sta1.r_GCRS
    
    earth = state_ss_t1[2]
    sun = state_ss_t1[-1]
    
    # Find potential U:
    U = GM[10]/norm(sun[:,0]-earth[:,0])
#    if debug: print 'U = {:.18f}'.format(U)
    
    # BCRS radius vectors of the first reception site at t_1:
    R_1 = earth[:,0] + (1.0 - U/(C**2) - L_C)*r_1 - \
          dot(earth[:,1], r_1)*earth[:,1] / (2.0*C**2)
          
    ''' calculate downleg light-time from S/C to Receiver
        to find signal transmission t_0 time given the reception time t_1
    '''    
    precision = 1e-16
    n_max = 3
    lag_order = 9
    
    # initial approximation:
    nn = 0
    lt_01_tmp = 0.0
    
    # s/c:    
    t_1, dd = t_1_days*86400, dd_days*86400
    x, _ = lagint(lag_order, tdb, bcrs[:,6], dd+t_1)
    y, _ = lagint(lag_order, tdb, bcrs[:,7], dd+t_1)
    z, _ = lagint(lag_order, tdb, bcrs[:,8], dd+t_1)
    vx, _ = lagint(lag_order, tdb, bcrs[:,9], dd+t_1)
    vy, _ = lagint(lag_order, tdb, bcrs[:,10], dd+t_1)
    vz, _ = lagint(lag_order, tdb, bcrs[:,11], dd+t_1)
    R_0 = np.hstack((x,y,z))
    V_0 = np.hstack((vx,vy,vz))

    if debug: print 't_1 = {:.18f}'.format(t_1)
    lt_01 = (norm(R_1 - R_0)/C)
    if debug: print 'lt_01 = {:.18f}'.format(lt_01)
    t_0 = t_1 - lt_01
    if debug: print 't_0_0 = {:.18f}'.format(t_0)

    while (abs(lt_01 - lt_01_tmp) > precision) and (nn < n_max):
        lt_01_tmp = lt_01
        t_0 = t_1 - lt_01
        
        x, _ = lagint(lag_order, tdb, bcrs[:,6], dd+t_0)
#        x, _ = lagint(lag_order, tdb/86400.0, bcrs[:,6], (dd+t_0)/86400.0)
        y, _ = lagint(lag_order, tdb, bcrs[:,7], dd+t_0)
        z, _ = lagint(lag_order, tdb, bcrs[:,8], dd+t_0)
        vx, _ = lagint(lag_order, tdb, bcrs[:,9], dd+t_0)
        vy, _ = lagint(lag_order, tdb, bcrs[:,10], dd+t_0)
        vz, _ = lagint(lag_order, tdb, bcrs[:,11], dd+t_0)
        R_0 = np.hstack((x,y,z))
        V_0 = np.hstack((vx,vy,vz))

        # vector needed for RLT calculation
        R_01 = R_1 - R_0

        # >> SS bodies
        RLT = 0.0
        for ii, state in enumerate(state_ss_t1):
            if not (ii==2 and norm(r_1)==0.0):
                rb = state[:,0]
                vb = state[:,1]
                R_0_B  = R_0 - (rb - lt_01*vb)
                R_1_B  = R_1 - rb
                R_01_B = R_1_B - R_0_B
    #            print R_0_B, R_1_B, R_01_B
                if debug: 
                    print 'rlt = {:.18f}'.format((2.0*GM[ii]/C**3) * \
                      log(  ( norm(R_0_B) + norm(R_1_B) + norm(R_01_B) + \
                              2.0*GM[ii]/C**2 ) / \
                            ( norm(R_0_B) + norm(R_1_B) - norm(R_01_B) + \
                              2.0*GM[ii]/C**2 ) ))
                RLT += (2.0*GM[ii]/C**3) * \
                      log(  ( norm(R_0_B) + norm(R_1_B) + norm(R_01_B) + \
                              2.0*GM[ii]/C**2 ) / \
                            ( norm(R_0_B) + norm(R_1_B) - norm(R_01_B) + \
                              2.0*GM[ii]/C**2 ) )
        
        lt_01 = lt_01 - (lt_01 - norm(R_01)/C - RLT) / \
                    ( 1.0 - dot(R_01, V_0)/(C*norm(R_01)) )

        t_0 = t_1 - lt_01
        if debug: print 't_0 = {:.18f}'.format(t_0)
        nn += 1
        if debug: print 'delta = {:.18f}'.format(abs(lt_01 - lt_01_tmp))
    
    x, _ = lagint(lag_order, tdb, bcrs[:,6], dd+t_0)
    y, _ = lagint(lag_order, tdb, bcrs[:,7], dd+t_0)
    z, _ = lagint(lag_order, tdb, bcrs[:,8], dd+t_0)
    vx, _ = lagint(lag_order, tdb, bcrs[:,9], dd+t_0)
    vy, _ = lagint(lag_order, tdb, bcrs[:,10], dd+t_0)
    vz, _ = lagint(lag_order, tdb, bcrs[:,11], dd+t_0)
    R_0 = np.hstack((x,y,z))
    V_0 = np.hstack((vx,vy,vz))
    
    if debug: 
        print 't_0 = {:.18f}'.format(t_0)
        print 'RLT = {:.18f}'.format(RLT)
        print 't_1 = {:.18f}'.format(t_1)
        print 'dd = {:.18f}'.format(dd)
    
    ''' B/GCRS state of second reception station at t_1'''
    # BCRS radius vectors of the transmitting site at t_1:
    r_2 = sta2.r_GCRS
    v_2 = sta2.v_GCRS
    a_2 = sta2.a_GCRS
    R_2_t_1 = earth[:,0] + (1.0 - U/(C**2) - L_C)*r_2 - \
          dot(earth[:,1], r_2)*earth[:,1] / (2.0*C**2)
    V_2_t_1 = earth[:,1] + \
            (1.0 - 2.0*U/C**2 - 0.5*(norm(earth[:,1])/C)**2 - \
             dot(earth[:,1], v_2)/C**2) * v_2 - \
            0.5*dot(earth[:,1], v_2)*earth[:,1]/C**2
    A_2_t_1 = earth[:,2] + \
            (1.0 - 3.0*U/C**2 - (norm(earth[:,1])/C)**2 + L_C - \
             2.0*dot(earth[:,1], v_2)/C**2) * a_2 - \
            0.5*dot(earth[:,1], a_2)*(earth[:,1] + 2.0*v_2)/C**2
    
#    print 'R2_t1 = [[{:.18f}],[{:.18f}],[{:.18f}]]'.format(*R_2_t_1)
    
    if debug:
        print R_2_t_1
        print V_2_t_1
        print A_2_t_1
        print 'lalala\n'
    ##############
        
    ''' BCRS state vectors of celestial bodies at t_0, [m, m/s]: '''
    ## Earth:
    rrd = pleph(jd+t_0/86400.0, 3, 12, inp['jpl_eph'])
    earth = np.reshape(np.asarray(rrd), (3,2), 'F') * 1e3
    # Earth's acceleration in m/s**2:
    v_plus = np.array(pleph(jd+t_0/86400.0+1.0/86400.0, 3, 12, inp['jpl_eph'])[3:])
    v_minus = np.array(pleph(jd+t_0/86400.0-1.0/86400.0, 3, 12, inp['jpl_eph'])[3:])
    a = (v_plus - v_minus)*1e3 / 2.0
    a = np.array(np.matrix(a).T)
    earth = np.hstack((earth, a))
    ## Sun:
    rrd = pleph(jd+t_0/86400.0, 11, 12, inp['jpl_eph'])
    sun = np.reshape(np.asarray(rrd), (3,2), 'F') * 1e3
    ## Moon:
    rrd = pleph(jd+t_0/86400.0, 10, 12, inp['jpl_eph'])
    moon = np.reshape(np.asarray(rrd), (3,2), 'F') * 1e3

    state_ss_t0 = []
    for jj in (1,2,4,5,6,7,8,9):
        rrd = pleph(jd+t_0/86400.0, jj, 12, inp['jpl_eph'])
        state_ss_t0.append(np.reshape(np.asarray(rrd), (3,2), 'F') * 1e3)
    state_ss_t0.insert(2,earth)
    state_ss_t0.append(moon)
    state_ss_t0.append(sun)
        
    # initial approximation
    nn = 0
    lt_02 = norm(R_2_t_1 - R_0)/C

    t_2 = t_0 + lt_02
    if debug: print 'initial t_2 = {:.18f}'.format(t_2)
#    lt_02_tmp = 0.0
    t_2_tmp = 0.0
    
#    while (abs(lt_02 - lt_02_tmp) > precision) and (nn < n_max):
    while (abs(t_2 - t_2_tmp) > precision) and (nn < n_max):
#        if debug: print 'lt_02 - lt_02_tmp = {:.18f}'.format(lt_02 - lt_02_tmp)
#        lt_02_tmp = deepcopy(lt_02)      
        if debug: print 't_2 - t_2_tmp = {:.18f}'.format(t_2 - t_2_tmp)
        t_2_tmp = deepcopy(t_2)
            
        R_2 = R_2_t_1 + V_2_t_1*(t_2-t_1) + 0.5*A_2_t_1*(t_2-t_1)**2
        V_2 = V_2_t_1 + A_2_t_1*(t_2-t_1)
        
        R_02 = R_2 - R_0
        
        # >> SS bodies
        ## Earth:
        rrd = pleph(jd+t_2/86400.0, 3, 12, inp['jpl_eph'])
        earth = np.reshape(np.asarray(rrd), (3,2), 'F') * 1e3
        ## Sun:
        rrd = pleph(jd+t_2/86400.0, 11, 12, inp['jpl_eph'])
        sun = np.reshape(np.asarray(rrd), (3,2), 'F') * 1e3
        ## Moon:
        rrd = pleph(jd+t_2/86400.0, 10, 12, inp['jpl_eph'])
        moon = np.reshape(np.asarray(rrd), (3,2), 'F') * 1e3
    
        state_ss_t2 = []
        for jj in (1,2,4,5,6,7,8,9):
            rrd = pleph(jd+t_2/86400.0, jj, 12, inp['jpl_eph'])
            state_ss_t2.append(np.reshape(np.asarray(rrd), (3,2), 'F') * 1e3)
        state_ss_t2.insert(2,earth)
        state_ss_t2.append(moon)
        state_ss_t2.append(sun)

        RLT = 0.0
        for ii, (state_t0, state_t2) in enumerate(zip(state_ss_t0,state_ss_t2)):
            if not (ii==2 and norm(r_2)<1e-3):
                rb_t0 = state_t0[:,0]
                rb_t2 = state_t2[:,0]
                R_0_B  = R_0 - rb_t0
                R_2_B  = R_2 - rb_t2
                R_02_B = R_2_B - R_0_B
                if debug:
                    print 'rlt = {:.18f}'.format((2.0*GM[ii]/C**3) * \
                      log(  ( norm(R_2_B) + norm(R_0_B) + norm(R_02_B) + \
                              2.0*GM[ii]/C**2 ) / \
                            ( norm(R_2_B) + norm(R_0_B) - norm(R_02_B) + \
                              2.0*GM[ii]/C**2 ) ))
                RLT += (2.0*GM[ii]/C**3) * \
                      log(  ( norm(R_2_B) + norm(R_0_B) + norm(R_02_B) + \
                              2.0*GM[ii]/C**2 ) / \
                            ( norm(R_2_B) + norm(R_0_B) - norm(R_02_B) + \
                              2.0*GM[ii]/C**2 ) )
                          
        t_2 = t_2 - (t_2 - t_0 - norm(R_02)/C - RLT) / \
                    ( 1.0 + dot(R_02, V_2)/(C*norm(R_02)) )

#        lt_02 = lt_02 - (lt_02 - norm(R_02)/C - RLT) / \
#                    ( 1.0 - dot(R_02, V_2)/(C*norm(R_02)) )
#        t_2 = t_0 + lt_02
        
        nn += 1
        if debug: print 't_2 = {:.18f}'.format(t_2)
        
    if debug: print 'final t_2 = {:.18f}'.format(t_2)
    
    if debug: print 'f(x)=0?: {:.18f}'.format(t_2 - t_0 - norm(R_02)/C - RLT)
    
    # gravitational effect on the ray path from Fukushima
    R_0_1 = R_0 - R_1
    R_0_2 = R_0 - R_2
    T_g_21 = 0.0
    for ii, (state_t0, state_t1, state_t2) in \
                enumerate(zip(state_ss_t0, state_ss_t1, state_ss_t2)):
        if ii!=2:
            rb_t0 = state_t0[:,0]
            rb_t1 = state_t1[:,0]
            rb_t2 = state_t2[:,0]
            R_0_B = R_0 - rb_t0
            R_1_B = R_1 - rb_t1
            R_2_B = R_2 - rb_t2
            T_g_21 += (2.0*GM[ii]/C**3) * \
                  log( (norm(R_2_B) + norm(R_0_B) + norm(R_0_2)) *\
                       (norm(R_1_B) + norm(R_0_B) - norm(R_0_1)) /\
                       (norm(R_2_B) + norm(R_0_B) - norm(R_0_2)) /
                       (norm(R_1_B) + norm(R_0_B) + norm(R_0_1)) )
            if debug: 
                print 'i={:d}, T_g_21_i = {:.18f}'.format(ii, (2.0*GM[ii]/C**3) * \
                  log( (norm(R_2_B) + norm(R_0_B) + norm(R_0_2)) *\
                       (norm(R_1_B) + norm(R_0_B) - norm(R_0_1)) /\
                       (norm(R_2_B) + norm(R_0_B) - norm(R_0_2)) /
                       (norm(R_1_B) + norm(R_0_B) + norm(R_0_1)) ))    
    
#    raw_input()
#    astropy_t_2 = Time(jd - 2400000.5, t_2/86400.0, \
#                format='mjd', scale='tdb', precision=9, \
#                 location=EarthLocation.from_geocentric(*sta2.r_GTRS, \
#                                                         unit=units.m))
    astropy_t_2 = Time(jd - 2400000.5, (t_2+T_g_21)/86400.0, \
                format='mjd', scale='tdb', precision=9, \
                 location=EarthLocation.from_geocentric(*sta2.r_GTRS, \
                                                         unit=units.m))
#    print 't_1 = {:.18f}'.format(t_1)
#    astropy_t_1 = Time(jd - 2400000.5, t_1/86400.0, \
#                format='mjd', scale='tdb', precision=9, \
#                 location=EarthLocation.from_geocentric(*sta1.r_GTRS, \
#                                                         unit=units.m))
#    astropy_t_1_a = Time(jd - 2400000.5, t_1_UTC, \
#                format='mjd', scale='utc', precision=9, \
#                 location=EarthLocation.from_geocentric(*sta1.r_GTRS, \
#                                                         unit=units.m))
#    print '{:.18f}'.format(t_1_UTC), '{:.18f}'.format(astropy_t_1.utc.jd2), \
#                                       '{:.18f}'.format(astropy_t_1_a.utc.jd2)
#    print '{:.18f}'.format(t_1_UTC*86400.0), \
#                                    astropy_t_1.utc.iso, astropy_t_1_a.utc.iso
#    print '{:.18f}'.format(t_1_days), '{:.18f}'.format(astropy_t_1_a.tdb.jd2)
#    print (t_1_UTC - astropy_t_1.utc.jd2)*86400.0

    t_2_UTC = astropy_t_2.utc.jd2
#    t_2_UTC = astropy_t_2.utc.mjd - np.floor(astropy_t_2.utc.mjd)
#    if t_2_UTC<0: t_2_UTC += 1
#    print t_2_UTC, astropy_t_2.utc.jd1, jd, t_1_UTC, dd_days
    # fixes on a day change
    if astropy_t_2.utc.jd1 < jd and t_2_UTC>0: t_2_UTC -= 1
    if astropy_t_2.utc.jd1 > jd and t_2_UTC<0: t_2_UTC += 1
#    if debug: print 't_2_UTC = {:.18f}'.format(t_2_UTC*86400.0)
    return (t_2_UTC-t_1_UTC)*86400.0

#==============================================================================
# 
#==============================================================================
def delay_nf_fukushima(jd, t_1_days, dd_days, state_ss_t1, tdb, bcrs, const,\
                    sta1, sta2):
    '''
    NF delay calculation
    '''
    debug = False
    
    GM = const.GM
    C = const.C
    L_C = const.L_C
    r_1 = sta1.r_GCRS
    r_2 = sta2.r_GCRS
    v_2 = sta2.v_GCRS
    
    earth = state_ss_t1[2]
    sun = state_ss_t1[-1]
    
    # Find potential U:
    U = GM[10]/norm(sun[:,0]-earth[:,0])
#    if debug: print 'U = {:.18f}'.format(U)
    
    # BCRS radius vectors of the reception sites at t_1:
    R_1 = earth[:,0] + (1.0 - U/(C**2) - L_C)*r_1 - \
          dot(earth[:,1], r_1)*earth[:,1] / (2.0*C**2)
    R_2 = earth[:,0] + (1.0 - U/(C**2) - L_C)*r_2 - \
          dot(earth[:,1], r_2)*earth[:,1] / (2.0*C**2)
    V_2 = earth[:,1] + \
            (1.0 - 2.0*U/C**2 - 0.5*(norm(earth[:,1])/C)**2 - \
             dot(earth[:,1], v_2)/C**2) * v_2 - \
            0.5*dot(earth[:,1], v_2)*earth[:,1]/C**2
          
    ''' calculate downleg light-time from S/C to Receiver
        to find signal transmission t_0 time given the reception time t_1
    '''    
    precision = 1e-16
    n_max = 3
    lag_order = 9
    
    # initial approximation:
    nn = 0
    lt_01_tmp = 0.0
    
    # s/c:    
    t_1, dd = t_1_days*86400, dd_days*86400
    x, _ = lagint(lag_order, tdb, bcrs[:,6], dd+t_1)
    y, _ = lagint(lag_order, tdb, bcrs[:,7], dd+t_1)
    z, _ = lagint(lag_order, tdb, bcrs[:,8], dd+t_1)
    vx, _ = lagint(lag_order, tdb, bcrs[:,9], dd+t_1)
    vy, _ = lagint(lag_order, tdb, bcrs[:,10], dd+t_1)
    vz, _ = lagint(lag_order, tdb, bcrs[:,11], dd+t_1)
    R_0 = np.hstack((x,y,z))
    V_0 = np.hstack((vx,vy,vz))

    if debug: print 't_1 = {:.18f}'.format(t_1)
    lt_01 = (norm(R_1 - R_0)/C)
    if debug: print 'lt_01 = {:.18f}'.format(lt_01)
    t_0 = t_1 - lt_01
    if debug: print 't_0_0 = {:.18f}'.format(t_0)

    while (abs(lt_01 - lt_01_tmp) > precision) and (nn < n_max):
        lt_01_tmp = lt_01
        t_0 = t_1 - lt_01
        
        x, _ = lagint(lag_order, tdb, bcrs[:,6], dd+t_0)
#        x, _ = lagint(lag_order, tdb/86400.0, bcrs[:,6], (dd+t_0)/86400.0)
        y, _ = lagint(lag_order, tdb, bcrs[:,7], dd+t_0)
        z, _ = lagint(lag_order, tdb, bcrs[:,8], dd+t_0)
        vx, _ = lagint(lag_order, tdb, bcrs[:,9], dd+t_0)
        vy, _ = lagint(lag_order, tdb, bcrs[:,10], dd+t_0)
        vz, _ = lagint(lag_order, tdb, bcrs[:,11], dd+t_0)
        R_0 = np.hstack((x,y,z))
        V_0 = np.hstack((vx,vy,vz))

        # vector needed for RLT calculation
        R_01 = R_1 - R_0

        # >> SS bodies
        RLT = 0.0
        for ii, state in enumerate(state_ss_t1):
            if not (ii==2 and norm(r_1)==0.0):
                rb = state[:,0]
                vb = state[:,1]
                R_0_B  = R_0 - (rb - lt_01*vb)
                R_1_B  = R_1 - rb
                R_01_B = R_1_B - R_0_B
    #            print R_0_B, R_1_B, R_01_B
                if debug: 
                    print 'rlt = {:.18f}'.format((2.0*GM[ii]/C**3) * \
                      log(  ( norm(R_0_B) + norm(R_1_B) + norm(R_01_B) + \
                              2.0*GM[ii]/C**2 ) / \
                            ( norm(R_0_B) + norm(R_1_B) - norm(R_01_B) + \
                              2.0*GM[ii]/C**2 ) ))
                RLT += (2.0*GM[ii]/C**3) * \
                      log(  ( norm(R_0_B) + norm(R_1_B) + norm(R_01_B) + \
                              2.0*GM[ii]/C**2 ) / \
                            ( norm(R_0_B) + norm(R_1_B) - norm(R_01_B) + \
                              2.0*GM[ii]/C**2 ) )
        
        lt_01 = lt_01 - (lt_01 - norm(R_01)/C - RLT) / \
                    ( 1.0 - dot(R_01, V_0)/(C*norm(R_01)) )

        t_0 = t_1 - lt_01
        if debug: print 't_0 = {:.18f}'.format(t_0)
        nn += 1
        if debug: print 'delta = {:.18f}'.format(abs(lt_01 - lt_01_tmp))
    
    x, _ = lagint(lag_order, tdb, bcrs[:,6], dd+t_0)
    y, _ = lagint(lag_order, tdb, bcrs[:,7], dd+t_0)
    z, _ = lagint(lag_order, tdb, bcrs[:,8], dd+t_0)
    vx, _ = lagint(lag_order, tdb, bcrs[:,9], dd+t_0)
    vy, _ = lagint(lag_order, tdb, bcrs[:,10], dd+t_0)
    vz, _ = lagint(lag_order, tdb, bcrs[:,11], dd+t_0)
    R_0 = np.hstack((x,y,z))
    V_0 = np.hstack((vx,vy,vz))
    
    if debug: 
        print 't_0 = {:.18f}'.format(t_0)
        print 'RLT = {:.18f}'.format(RLT)
        print 't_1 = {:.18f}'.format(t_1)
        print 'dd = {:.18f}'.format(dd)
    
    # That's the confusing way Fukushima defines it. suck it off!
    R_0_1 = R_0 - R_1
    R_0_2 = R_0 - R_2
    
    # psuedo source-vector
    K = (R_0_1 + R_0_2)/(norm(R_0_1) + norm(R_0_2))
    
    R_2_hat = R_0_2/norm(R_0_2)
    b = r_2 - r_1
    H = dot(K, b)/(2.0*norm(R_0_2)) * norm(np.cross(V_2/C, R_2_hat))**2
    
    if debug:
        print 'R_2_hat = {:s}'.format(R_2_hat)
        print 'b = {:s}'.format(b)
        print 'H = {:.18f}'.format(H)
        print 'K = {:s}'.format(K)
    
    # gravitational effect on the ray path
    T_g_21 = 0.0
    for ii, state in enumerate(state_ss_t1):
        if ii!=2:
            rb = state[:,0]
            vb = state[:,1]
            R_0_B = R_0 - (rb - lt_01*vb)
            R_1_B = R_1 - rb
            R_2_B = R_2 - rb
            T_g_21 += (2.0*GM[ii]/C**3) * \
                  log( (norm(R_2_B) + norm(R_0_B) + norm(R_0_2)) *\
                       (norm(R_1_B) + norm(R_0_B) - norm(R_0_1)) /\
                       (norm(R_2_B) + norm(R_0_B) - norm(R_0_2)) /
                       (norm(R_1_B) + norm(R_0_B) + norm(R_0_1)) )
            if debug: 
                print 'i={:d}, T_g_21_i = {:.18f}'.format(ii, (2.0*GM[ii]/C**3) * \
                  log( (norm(R_2_B) + norm(R_0_B) + norm(R_0_2)) *\
                       (norm(R_1_B) + norm(R_0_B) - norm(R_0_1)) /\
                       (norm(R_2_B) + norm(R_0_B) - norm(R_0_2)) /
                       (norm(R_1_B) + norm(R_0_B) + norm(R_0_1)) ))

    if debug: 
        print 'T_g_21 = {:.18f}'.format(T_g_21)

    # delay in the TT-frame
    V_E = earth[:,1]
    dtau = (-(1 - 2*U/C**2 - \
            (norm(V_E)**2 + 2*dot(V_E,v_2))/(2*C**2))*dot(K,b)/C  \
            - dot(V_E, b)/C**2 * (1 + dot(R_2_hat, V_2)/C - \
            dot(V_E+2*v_2, K)/(2*C)) + T_g_21) \
            / ((1 + dot(R_2_hat, V_2)/C)*(1 + H))
    if debug: print 'dtau = {:.18f}'.format(dtau)

    return dtau


#==============================================================================
# 
#==============================================================================
def doppler_bc(tjd, t_1, dd, state_ss_t1, tdb, bcrs,
                GM, TDB_TCB, L_C, C, x_way, freq_type,
                ut1, sta1, sta2=None, L_G=None, AE=None, J_2=None,
                utc=None, gcrs=None, t_utc=None, eops=None, inp=None):
#                r_1, v_1, r_2=None, v_2=None, a_2=None):
    '''
    Doppler calculation following Moyer/Duev
    For reference see Duev PhD thesis, MSU 2012.
    
    dd - number of days since the start epoch of the ephemeris * 86400
    state_ss_t1  - Solar system bodies r, v (and a for Earth) at t_1 wrt SSBC
    '''
    debug = False
    
    r_1 = sta1.r_GCRS
    v_1 = sta1.v_GCRS
    if sta2!=None:
        r_2 = sta2.r_GCRS
        v_2 = sta2.v_GCRS
        a_2 = sta2.a_GCRS
    
    earth = state_ss_t1[2]
    sun = state_ss_t1[-1]
    
    # Find potential U:
    U = GM[10]/norm(sun[:,0]-earth[:,0])
    if debug: print 't_1 = {:.18f}'.format(t_1)

    # BCRS radius vectors of the reception site at t_1:
    R_1 = earth[:,0] + (1.0 - U/(C**2) - L_C)*r_1 - \
          dot(earth[:,1],r_1)*earth[:,1] / (2.0*C**2)
#    V_1 = earth[:,1] + \
#            (1.0 - 2.0*U/C**2 - 0.5*(norm(earth[:,1])/C)**2 - \
#             dot(earth[:,1], v_1)/C**2) * v_1 - \
#            0.5*dot(earth[:,1], v_1)*earth[:,1]/C**2
    V_1 = earth[:,1] + ( (1.0 - U/(C**2) - L_C)*v_1 - \
          dot(earth[:,1],v_1)*earth[:,1] / (2.0*C**2) )* \
          (1.0 - (U + v_1**2/2.0 - L_C)/C**2)
    

    ''' calculate downleg light-time from S/C to Receiver
        to find signal transmission t_0 time given the reception time t_1
    '''    
    precision = 1e-16
    n_max = 3
    lag_order = 9
    
    # initial approximation:
    nn = 0
    lt_01_tmp = 0.0
    
    # s/c:    
    x, _ = lagint(lag_order, tdb, bcrs[:,6], dd+t_1)
    y, _ = lagint(lag_order, tdb, bcrs[:,7], dd+t_1)
    z, _ = lagint(lag_order, tdb, bcrs[:,8], dd+t_1)
    vx, _ = lagint(lag_order, tdb, bcrs[:,9], dd+t_1)
    vy, _ = lagint(lag_order, tdb, bcrs[:,10], dd+t_1)
    vz, _ = lagint(lag_order, tdb, bcrs[:,11], dd+t_1)
    R_0 = np.hstack((x,y,z))
    V_0 = np.hstack((vx,vy,vz))

    lt_01 = (norm(R_1 - R_0)/C)
    if debug: print 'lt_01 = {:.18f}'.format(lt_01)
    t_0 = t_1 - lt_01
    if debug: print 't_0_0 = {:.18f}'.format(t_0)

    while (abs(lt_01 - lt_01_tmp) > precision) and (nn < n_max):
        lt_01_tmp = lt_01
        t_0 = t_1 - lt_01
        
        x, _ = lagint(lag_order, tdb, bcrs[:,6], dd+t_0)
        y, _ = lagint(lag_order, tdb, bcrs[:,7], dd+t_0)
        z, _ = lagint(lag_order, tdb, bcrs[:,8], dd+t_0)
        vx, _ = lagint(lag_order, tdb, bcrs[:,9], dd+t_0)
        vy, _ = lagint(lag_order, tdb, bcrs[:,10], dd+t_0)
        vz, _ = lagint(lag_order, tdb, bcrs[:,11], dd+t_0)
        R_0 = np.hstack((x,y,z))
        V_0 = np.hstack((vx,vy,vz))

        # vector needed for RLT calculation
        R_01 = R_1 - R_0

        # >> SS bodies
        RLT = 0.0
        for ii, state in enumerate(state_ss_t1):
#            if debug: 
#                print 'keeping Sun only'
#                if ii != 10: continue
            if ii==2 and norm(r_1)==0.0: continue
            rb = state[:,0]
            vb = state[:,1]
            R_0_B  = R_0 - (rb - lt_01*vb)
            R_1_B  = R_1 - rb
            R_01_B = R_1_B - R_0_B
#            print R_0_B, R_1_B, R_01_B
            if debug: 
                print 'rlt = {:.18f}'.format((2.0*GM[ii]/C**3) * \
                  log(  ( norm(R_0_B) + norm(R_1_B) + norm(R_01_B) + \
                          2.0*GM[ii]/C**2 ) / \
                        ( norm(R_0_B) + norm(R_1_B) - norm(R_01_B) + \
                          2.0*GM[ii]/C**2 ) ))
            RLT += (2.0*GM[ii]/C**3) * \
                  log(  ( norm(R_0_B) + norm(R_1_B) + norm(R_01_B) + \
                          2.0*GM[ii]/C**2 ) / \
                        ( norm(R_0_B) + norm(R_1_B) - norm(R_01_B) + \
                          2.0*GM[ii]/C**2 ) )
        
        lt_01 = lt_01 - (lt_01 - norm(R_01)/C - RLT) / \
                    ( 1.0 - dot(R_01, V_0)/(C*norm(R_01)) )

        t_0 = t_1 - lt_01
        if debug: print 't_0 = {:.18f}'.format(t_0)
        nn += 1
        if debug: print 'delta = {:.18f}'.format(abs(lt_01 - lt_01_tmp))
    
    x, _ = lagint(lag_order, tdb, bcrs[:,6], dd+t_0)
    y, _ = lagint(lag_order, tdb, bcrs[:,7], dd+t_0)
    z, _ = lagint(lag_order, tdb, bcrs[:,8], dd+t_0)
    vx, _ = lagint(lag_order, tdb, bcrs[:,9], dd+t_0)
    vy, _ = lagint(lag_order, tdb, bcrs[:,10], dd+t_0)
    vz, _ = lagint(lag_order, tdb, bcrs[:,11], dd+t_0)
    R_0 = np.hstack((x,y,z))
    V_0 = np.hstack((vx,vy,vz))
    
    if debug: 
        print 't_0 = {:.18f}'.format(t_0)
        print 'RLT = {:.18f}'.format(RLT)
        print 't_1 = {:.18f}'.format(t_1)
        print 'dd = {:.18f}'.format(dd)
    
    if debug:
        rrd = pleph(2456266.5+t_0/86400.0, 2, 12, 'jpl_eph/JPLEPH.421')*1e3
        print 'xsat = [[{:.18f}],[{:.18f}],[{:.18f}]]'.format(*(R_0-rrd[:3]))
        print 'obs_sat = [[{:.18f}],[{:.18f}],[{:.18f}]]'.format(*(R_1 - R_0))
    
    ''' BCRS state vectors of celestial bodies at t_0, [m, m/s]: '''
    ## Earth:
    JD = tjd
    rrd = pleph(JD+t_0/86400.0, 3, 12, inp['jpl_eph'])
    earth = np.reshape(np.asarray(rrd), (3,2), 'F') * 1e3
    # Earth's acceleration in m/s**2:
#    v_plus = np.array(pleph(JD+t_0/86400.0+1.0/86400.0, 3, 12, inp['jpl_eph'])[3:])
#    v_minus = np.array(pleph(JD+t_0/86400.0-1.0/86400.0, 3, 12, inp['jpl_eph'])[3:])
#    a = (v_plus - v_minus)*1e3 / 2.0
#    a = np.array(np.matrix(a).T)
#    earth = np.hstack((earth, a))
    ## Sun:
    rrd = pleph(JD+t_0/86400.0, 11, 12, inp['jpl_eph'])
    sun = np.reshape(np.asarray(rrd), (3,2), 'F') * 1e3
    ## Moon:
    rrd = pleph(JD+t_0/86400.0, 10, 12, inp['jpl_eph'])
    moon = np.reshape(np.asarray(rrd), (3,2), 'F') * 1e3

    state_ss_t0 = []
    for jj in (1,2,4,5,6,7,8,9):
        rrd = pleph(JD+t_0/86400.0, jj, 12, inp['jpl_eph'])
        state_ss_t0.append(np.reshape(np.asarray(rrd), (3,2), 'F') * 1e3)
    state_ss_t0.insert(2,earth)
    state_ss_t0.append(moon)
    state_ss_t0.append(sun)
    
    ''' My algorithm from PhD thesis '''
    # direction vector
#    print R_0, R_1
    n_b = (R_0 - R_1)/norm(R_0 - R_1)
#    print n_b
    # >> calculate dtau_s(t_s)/dTCB, dtau_o(t_s)/dTCB and f_b(t_o)/f_b(t_s)
    GMnaR_s = 0.0
    GMnaR_o = 0.0
    z_sh = 0.0
    for ii, (state_t0, state_t1, gm) in enumerate(zip(state_ss_t0, state_ss_t1, GM)):
        if ii==2 and norm(r_1)==0.0: continue
        rb_t0 = state_t0[:,0]
        vb_t0 = state_t0[:,1]
        rb_t1 = state_t1[:,0]
        vb_t1 = state_t1[:,1]
        R_0_B  = R_0 - rb_t0
        R_1_B  = R_1 - rb_t1
        R_01_B = R_1_B - R_0_B
        V_0_B  = V_0 - vb_t0
        V_1_B  = V_1 - vb_t1
        GMnaR_s += gm/norm(R_0_B)
        GMnaR_o += gm/norm(R_1_B)
        z_sh += (4.0*gm/C**3) * \
            ( (norm(R_1_B) + norm(R_0_B))*dot(R_01_B, V_1_B-V_0_B)/norm(R_01_B) - \
            norm(R_01_B)*(dot(R_1_B,V_1_B)/norm(R_1_B) + \
            dot(R_0_B,V_0_B)/norm(R_0_B) + 2.0*gm/C**2) ) / \
            ( (norm(R_1_B)+norm(R_0_B) + 2.0*gm/C**2)**2 - norm(R_01_B)**2 )

    z_sh = -z_sh
#    print 'z_sh = {:.18e}'.format(z_sh)
    
    dtau_s_po_dTCB = 1.0 - (GMnaR_s + (norm(V_0)**2)/2.0)/C**2
#    print 'dtau_s_po_dTCB = {:.18f}'.format(dtau_s_po_dTCB)    
    dtau_o_po_dTCB = 1.0 - (GMnaR_o + (norm(V_1)**2)/2.0)/C**2
#    print 'dtau_o_po_dTCB = {:.18f}'.format(dtau_o_po_dTCB)

    # >> put everything together
    fonafe = (1.0+z_sh) * dtau_s_po_dTCB * (1.0 + dot(n_b,V_1)/C) / \
         ( dtau_o_po_dTCB * (1.0 + dot(n_b,V_0)/C) )
#    print 'fonafe = {:.18f}'.format(fonafe)
    
    # correction:
#    U_1 = GM[10]/norm(sun[:,0]-R_1)
#    U_E_1 = GM[2]/norm(r_1) + \
#        GM[2]*AE**2*J_2*(1-3*(r_1[2]/norm(r_1))**2)/(2.0*norm(r_1)**3)
#    
#    UTC_to_proper = 1 - norm(V_1)**2/(2*C**2) - ((norm(v_1)**2)/2.0 + U_E_1)/C**2
    
#    fonafe /= UTC_to_proper
    
    ''' correct fonafe if f is given in GC (i.e. it's not proper): '''
    if x_way=='one' and freq_type=='gc':
#        raise NotImplemented
        x, _ = lagint(lag_order, utc, gcrs[:,6], dd+t_utc-lt_01)
        y, _ = lagint(lag_order, utc, gcrs[:,7], dd+t_utc-lt_01)
        z, _ = lagint(lag_order, utc, gcrs[:,8], dd+t_utc-lt_01)
        vx, _ = lagint(lag_order, utc, gcrs[:,9], dd+t_utc-lt_01)
        vy, _ = lagint(lag_order, utc, gcrs[:,10], dd+t_utc-lt_01)
        vz, _ = lagint(lag_order, utc, gcrs[:,11], dd+t_utc-lt_01)
        r_sc = np.hstack((x,y,z))
        v_sc = np.hstack((vx,vy,vz))
        
        U_E_sc = GM[2]/norm(r_sc) + \
            GM[2]*AE**2*J_2*(1-3*(r_sc[2]/norm(r_sc))**2)/(2.0*norm(r_sc)**3)

        fonafe /= (1 + L_G - ((norm(v_sc)**2)/2.0 + U_E_sc)/C**2)
        return fonafe, None
    
    if x_way=='one':
        # transmission event in TDB (used if frequency is ramped):
        astropy_t_0 = Time(tjd - 2400000.5, t_0/86400.0, \
                            format='mjd', scale='tdb', precision=9)
        return fonafe, astropy_t_0 # else simply pass
        
    ''' 2(3)-way Doppler:
        calculate upleg light-time from Transmitter to S/C
        to find signal transmission time t_2 given the reception time t_0
    '''
    # BCRS radius vectors of the transmitting site at t_1:
    R_2_t_1 = earth[:,0] + (1.0 - U/(C**2) - L_C)*r_2 - \
          dot(earth[:,1],r_2)*earth[:,1] / (2.0*C**2)
#    V_2_t_1 = earth[:,1] + \
#            (1.0 - 2.0*U/C**2 - 0.5*(norm(earth[:,1])/C)**2 - \
#             dot(earth[:,1], v_2)/C**2) * v_2 - \
#            0.5*dot(earth[:,1], v_2)*earth[:,1]/C**2
#    A_2_t_1 = earth[:,2] + \
#            (1.0 - 3.0*U/C**2 - (norm(earth[:,1])/C)**2 + L_C - \
#             2.0*dot(earth[:,1], v_2)/C**2) * a_2 - \
#            0.5*dot(earth[:,1], a_2)*(earth[:,1] + 2.0*v_2)/C**2
    
#    print 'R2_t1 = [[{:.18f}],[{:.18f}],[{:.18f}]]'.format(*R_2_t_1)
    
    ''' calculate upleg light-time from Transmitter to S/C
        to find signal transmission t_2 time given the reception time t_0
    '''
    # initial approximation:
    nn = 0
    lt_20_tmp = 0.0
    
    # initial approximation using R_2_t_1:
    lt_20 = (norm(R_0 - R_2_t_1)/C)
    if debug: print 'lt_20_0 = {:.18f}'.format(lt_20)
    t_2 = t_0 - lt_20
    if debug: print 't_2_0 = {:.18f}'.format(t_2)


    ###############
    JD = tjd
    mjd = tjd - 2400000.5
#    print JD, mjd
#    astropy_t_2 = Time(mjd, t_2/86400.0, \
#                format='mjd', scale='tdb', precision=9, \
#                 location=(sta2.lon_gcen*180/np.pi*units.deg, \
#                           sta2.lat_geod*180/np.pi*units.deg,\
#                           sta2.h_geod))
    astropy_t_2 = Time(mjd, t_2/86400.0, \
                format='mjd', scale='tdb', precision=9, \
                 location=EarthLocation.from_geocentric(*sta2.r_GTRS, \
                                                         unit=units.m))
    # t_2 might be negative. redefine JD and mjd therefore? wozu?
#    if t_2 < 0:
#        mjd = np.floor(astropy_t_2.utc.mjd)
#        JD = mjd + 2400000.5
#    print JD, mjd, astropy_t_2.utc.jd1, astropy_t_2.utc.jd2
#    raw_input()
#    UTC = astropy_t_2.utc.mjd - np.floor(astropy_t_2.utc.mjd)
    UTC = astropy_t_2.utc.jd2
    if debug: print 'UTC = {:.18f}'.format(UTC)
    mjd_keep = mjd
    if UTC<0:
#        UTC, JD, mjd = UTC+1, JD-1, mjd-1
        UTC, JD, mjd = UTC+1, astropy_t_2.utc.jd1-1, astropy_t_2.utc.mjd-1
    
    ''' compute tai & tt '''
    TAI, TT = taitime(mjd, UTC)
    if debug: 
        print 'UTC = {:.18f}'.format(UTC)
        print 'TAI = {:.18f}'.format(TAI)
        print 'TT = {:.18f}'.format(TT)
    ''' interpolate eops to tstamp '''
    UT1, eop_int = eop_iers(mjd, UTC, eops)
    if debug: print 'UT1 = {:.18f}'.format(UT1)
    
    ''' compute coordinate time fraction of CT day at 2nd observing site '''
    CT, dTAIdCT = t_eph(JD, UT1, TT, sta2.lon_gcen, sta2.u, sta2.v)

    ''' BCRS state vectors of celestial bodies at JD+CT, [m, m/s]: '''
    ## Earth:
    rrd = pleph(JD+CT, 3, 12, inp['jpl_eph'])
    earth = np.reshape(np.asarray(rrd), (3,2), 'F') * 1e3
    # Earth's acceleration in m/s**2:
    v_plus = np.array(pleph(JD+CT+1.0/86400.0, 3, 12, inp['jpl_eph'])[3:])
    v_minus = np.array(pleph(JD+CT-1.0/86400.0, 3, 12, inp['jpl_eph'])[3:])
    a = (v_plus - v_minus)*1e3 / 2.0
    a = np.array(np.matrix(a).T)
    earth = np.hstack((earth, a))
    ## Sun:
    rrd = pleph(JD+CT, 11, 12, inp['jpl_eph'])
    sun = np.reshape(np.asarray(rrd), (3,2), 'F') * 1e3
    ## Moon:
    rrd = pleph(JD+CT, 10, 12, inp['jpl_eph'])
    moon = np.reshape(np.asarray(rrd), (3,2), 'F') * 1e3

    state_ss_t2 = []
    for jj in (1,2,4,5,6,7,8,9):
        rrd = pleph(JD+CT, jj, 12, inp['jpl_eph'])
        state_ss_t2.append(np.reshape(np.asarray(rrd), (3,2), 'F') * 1e3)
    state_ss_t2.insert(2,earth)
    state_ss_t2.append(moon)
    state_ss_t2.append(sun)

    ''' rotation matrix IERS '''
    tstamp = astropy_t_2.utc.datetime
    r2000 = ter2cel(tstamp, eop_int, dTAIdCT, 'iau2000')

    ''' displacements due to geophysical effects '''
    if sta2.name == 'GEOCENTR':
        pass
    else:
    # displacement due to solid Earth tides:
        sta2 = dehanttideinel(sta2, tstamp, earth, sun, moon, r2000)
    # displacement due to ocean loading:
        sta2 = hardisp(sta2, tstamp, r2000)
    # rotational deformation due to pole tide:
        sta2 = poletide(sta2, tstamp, eop_int, r2000)

    
    ''' add up geophysical corrections and convert sta state to J2000 '''
    sta2.j2000gp(r2000)
    
    r_2 = sta2.r_GCRS
    v_2 = sta2.v_GCRS
    a_2 = sta2.a_GCRS
    
    # BCRS radius vectors of the transmitting site at t_2_0:
    R_2_t_2 = earth[:,0] + (1.0 - U/(C**2) - L_C)*r_2 - \
          dot(earth[:,1],r_2)*earth[:,1] / (2.0*C**2)
    V_2_t_2 = earth[:,1] + \
            (1.0 - 2.0*U/C**2 - 0.5*(norm(earth[:,1])/C)**2 - \
             dot(earth[:,1], v_2)/C**2) * v_2 - \
            0.5*dot(earth[:,1], v_2)*earth[:,1]/C**2
    A_2_t_2 = earth[:,2] + \
            (1.0 - 3.0*U/C**2 - (norm(earth[:,1])/C)**2 + L_C - \
             2.0*dot(earth[:,1], v_2)/C**2) * a_2 - \
            0.5*dot(earth[:,1], a_2)*(earth[:,1] + 2.0*v_2)/C**2
    
    t_2_0 = deepcopy(t_2)
    if debug:
        print R_2_t_2
        print V_2_t_2
        print A_2_t_2
        print 'lalala'
    ##############

    while (abs(lt_20 - lt_20_tmp) > precision) and (nn < n_max):
        lt_20_tmp = deepcopy(lt_20)
        t_2 = t_0 - lt_20
#        lt_201 = lt_20 + lt_01
        
#        R_2 = R_2_t_1 - V_2_t_1*lt_201 + 0.5*A_2_t_1*lt_201**2
#        V_2 = V_2_t_1 - A_2_t_1*lt_201
        
        R_2 = R_2_t_2 + V_2_t_2*(t_2-t_2_0) + 0.5*A_2_t_2*(t_2-t_2_0)**2
        V_2 = V_2_t_2 + A_2_t_2*(t_2-t_2_0)
        
        if debug:
            print R_2
            print V_2
        # vector needed for RLT calculation
        R_20 = R_0 - R_2

        # >> SS bodies
        RLT = 0.0
        for ii, (state_t0, state_t2) in enumerate(zip(state_ss_t0,state_ss_t2)):
#            if debug: 
#                print 'keeping Sun and Earth only'
#                if ii != 10: continue
            if ii==2 and norm(r_2)==0.0: continue
            rb_t0 = state_t0[:,0]
#            vb_t0 = state_t0[:,1]
            rb_t2 = state_t2[:,0]
            vb_t2 = state_t2[:,1]
            R_0_B  = R_0 - rb_t0
            R_2_B  = R_2 - (rb_t2 + (t_2-t_2_0)*vb_t2)
            R_20_B = R_0_B - R_2_B
            if debug:
                print 'rlt = {:.18f}'.format((2.0*GM[ii]/C**3) * \
                  log(  ( norm(R_2_B) + norm(R_0_B) + norm(R_20_B) + \
                          2.0*GM[ii]/C**2 ) / \
                        ( norm(R_2_B) + norm(R_0_B) - norm(R_20_B) + \
                          2.0*GM[ii]/C**2 ) ))
            RLT += (2.0*GM[ii]/C**3) * \
                  log(  ( norm(R_2_B) + norm(R_0_B) + norm(R_20_B) + \
                          2.0*GM[ii]/C**2 ) / \
                        ( norm(R_2_B) + norm(R_0_B) - norm(R_20_B) + \
                          2.0*GM[ii]/C**2 ) )
        
        lt_20 = lt_20 - (lt_20 - norm(R_20)/C - RLT) / \
                    ( 1.0 - dot(R_20, V_2)/(C*norm(R_20)) )

        t_2 = t_0 - lt_20
        if debug: print 't_2 = {:.18f}'.format(t_2)
        nn += 1
        if debug: print 'delta = {:.18f}'.format(abs(lt_20 - lt_20_tmp))
    
    # transmitter state at found t_2
#    lt_201 = lt_20 + lt_01
#    R_2 = R_2_t_1 - V_2_t_1*lt_201 + 0.5*A_2_t_1*lt_201**2
#    V_2 = V_2_t_1 - A_2_t_1*lt_201
    R_2 = R_2_t_2 + V_2_t_2*(t_2-t_2_0) + 0.5*A_2_t_2*(t_2-t_2_0)**2
    V_2 = V_2_t_2 + A_2_t_2*(t_2-t_2_0)
    
    if debug:
        print 't_0 = {:.18f}'.format(t_0)
        print 'RLT = {:.18f}'.format(RLT)
        print 't_2 = {:.18f}'.format(t_2)
        print 'dd = {:.18f}'.format(dd)
    
        print 'obs_sat = [[{:.18f}],[{:.18f}],[{:.18f}]]'.format(*(R_2 - R_0))
        print 'R2_t2 = [[{:.18f}],[{:.18f}],[{:.18f}]]'.format(*R_2)

    # unfix mjd if necessary:
    mjd = mjd_keep
#    astropy_t_2 = Time(mjd, t_2/86400.0, \
#                format='mjd', scale='tdb', precision=9, \
#                 location=(sta2.lon_gcen*180/np.pi, \
#                           sta2.lat_geod*180/np.pi,\
#                           sta2.h_geod))
    if debug: print 'mjd = {:.18f}, t_2 = {:.18f}'.format(mjd, t_2/86400.0)
    astropy_t_2 = Time(mjd, t_2/86400.0, \
                format='mjd', scale='tdb', precision=9, \
                 location=EarthLocation.from_geocentric(*sta2.r_GTRS, \
                                                         unit=units.m))
    if debug: print astropy_t_2.datetime
#    t_2_UTC = astropy_t_2.utc.mjd - np.floor(astropy_t_2.utc.mjd)
    t_2_UTC = astropy_t_2.utc.jd2
    if debug:
        print 't_2_UTC = {:.18f}'.format(t_2_UTC*86400.0)
#    astropy_t_2 = Time(mjd, t_2/86400.0, \
#                format='mjd', scale='tdb', precision=9, \
#                 location=(sta1.lon_gcen*180/np.pi, \
#                           sta1.lat_geod*180/np.pi,\
#                           sta1.h_geod))
#    t_2_UTC = astropy_t_2.utc.jd2
#    if debug:
#        print 't_2_UTC = {:.18f}'.format(t_2_UTC*86400.0)
#    astropy_t_2 = Time(mjd, t_2/86400.0, \
#                format='mjd', scale='tdb', precision=9)
#    t_2_UTC = astropy_t_2.utc.jd2
#    if debug:
#        print 't_2_UTC = {:.18f}'.format(t_2_UTC*86400.0)
#        raw_input()
    
    ''' My algorithm from PhD thesis '''
    # direction vector
#    print R_2, R_0
    n_b = (R_2 - R_0)/norm(R_2 - R_0)
#    print n_b
    # >> calculate dtau_s(t_s)/dTCB, dtau_o(t_s)/dTCB and f_b(t_o)/f_b(t_s)
    GMnaR_s = 0.0
    GMnaR_o = 0.0
    z_sh = 0.0
    for ii, (state_t0, state_t2, gm) in enumerate(zip(state_ss_t0, state_ss_t2, GM)):
        if ii==2 and norm(r_2)==0.0: continue
        rb_t0 = state_t0[:,0]
        vb_t0 = state_t0[:,1]
        rb_t2 = state_t2[:,0]
        vb_t2 = state_t2[:,1]
        R_0_B  = R_0 - rb_t0
        R_2_B  = R_2 - (rb_t2 + (t_2-t_2_0)*vb_t2)
        R_20_B = R_0_B - R_2_B
        V_0_B  = V_0 - vb_t0
        V_2_B  = V_2 - vb_t2
        GMnaR_s += gm/norm(R_2_B)
        GMnaR_o += gm/norm(R_0_B)
        z_sh += (4.0*gm/C**3) * \
            ( (norm(R_0_B) + norm(R_2_B))*dot(R_20_B, V_0_B-V_2_B)/norm(R_20_B) - \
            norm(R_20_B)*(dot(R_0_B,V_0_B)/norm(R_0_B) + \
            dot(R_2_B,V_2_B)/norm(R_2_B) + 2.0*gm/C**2) ) / \
            ( (norm(R_0_B)+norm(R_2_B) + 2.0*gm/C**2)**2 - norm(R_20_B)**2 )

    z_sh = -z_sh
#    print 'z_sh = {:.18e}'.format(z_sh)
    
    dtau_s_po_dTCB = 1.0 - (GMnaR_s + (norm(V_2)**2)/2.0)/C**2
#    print 'dtau_s_po_dTCB = {:.18f}'.format(dtau_s_po_dTCB)    
    dtau_o_po_dTCB = 1.0 - (GMnaR_o + (norm(V_0)**2)/2.0)/C**2
#    print 'dtau_o_po_dTCB = {:.18f}'.format(dtau_o_po_dTCB)

    # >> put everything together
    fonafe *= (1.0+z_sh) * dtau_s_po_dTCB * (1.0 + dot(n_b,V_0)/C) / \
         ( dtau_o_po_dTCB * (1.0 + dot(n_b,V_2)/C) )
#    print 'fonafe = {:.18f}'.format(fonafe)
#    raw_input()
    
    # correction if station 1 is GC:
    if sta1.name == 'GEOCENTR':
        fonafe *= (1 + L_G)
    
    # correction:
#    U_2 = GM[10]/norm((sun[:,0]+(t_2-t_2_0)*sun[:,1])-R_2)
#    r_ = r_2 + (t_2-t_2_0)*v_2 + (t_2-t_2_0)**2*a_2/2.0
#    v_ = v_2 + (t_2-t_2_0)*a_2
#    U_E_2 = GM[2]/norm(r_) + \
#        GM[2]*AE**2*J_2*(1-3*(r_[2]/norm(r_))**2)/(2.0*norm(r_)**3)
#    if debug:        
#        print 'U_2 = {:.18f}'.format(U_2)
#        print 'U_2/C**2 = {:.18f}'.format(U_2/C**2)
#        print 'V_2**2/(2*C**2) = {:.18f}'.format(norm(V_2)**2/(2*C**2))
#        print 'L_G, U_E_2/C**2 = {:.18f} {:.18f}'.format(L_G, U_E_2/C**2)
#        print '{:.18f} {:.18f}'.format(L_G, ((norm(v_)**2)/2.0 + U_E_2)/C**2)
#    
#    UTC_to_proper = 1 - norm(V_2)**2/(2*C**2) - ((norm(v_)**2)/2.0 + U_E_2)/C**2
#    
#    if debug:
#        print 't_2_proper = {:.18f}'.format(t_2_UTC*86400.0*UTC_to_proper)
#        raw_input()

#    fonafe *= UTC_to_proper

#    U_2 = GM[10]/norm(sun[:,0]-R_2)
#    fonafe /= (1 + L_G - ((norm(V_2)**2)/2.0 + U_2)/C**2)
#    fonafe /= (1 + L_G)
#    U_E_1 = GM[2]/norm(r_1) + \
#            GM[2]*AE**2*J_2*(1-3*(r_1[2]/norm(r_1))**2)/(2.0*norm(r_1)**3)
#    fonafe *= (1 + ((norm(v_1)**2)/2.0 + U_E_1)/C**2)
    
#    U_E_2 = GM[2]/norm(r_2) + \
#            GM[2]*AE**2*J_2*(1-3*(r_2[2]/norm(r_2))**2)/(2.0*norm(r_2)**3)
#    fonafe /= (1 + ((norm(v_2)**2)/2.0 + U_E_2)/C**2)
    
    # convert t1 and t2 to TT
#    t1_dtdb = dtdb(tjd, t_1/86400., ut1/86400., sta1.lon_gcen, sta1.u, sta1.v)
#    t1_tt = t_1 + t1_dtdb
#    t_2 = t_0 - lt_20
#    t2_dtdb = dtdb(tjd, t_2/86400., (ut1-lt_201)/86400., \
#                      sta2.lon_gcen, sta2.u, sta2.v)
#    t2_tt = t_2 + t2_dtdb
    # lt_201 in TT. (t1_tt-t2_tt) = t_1-t_2 + (t1_dtdb-t2_dtdb)
#    lt_201_tt = lt_201 + (t1_dtdb - t2_dtdb)

    return fonafe, astropy_t_2.utc
    
'''
#==============================================================================
# 
#==============================================================================
'''
def pointings(source, stations, date_t_start, date_t_stop, t_step, cfg,
              output=False):
    '''
    Compute pointings on spacecraft for a list of stations
    
    date_t_start and date_t_stop - datetime.datetime objects
                                   defining obs time slot
    t_step - time step in seconds
    
    cfg - input config file location

    output - produce a txt-damp or not
    
    output path is set in cfg['out_path']
    
    Function outputs an array of datetime objects with obs epochs,
    RA/Decs for J2000 and date in radians and Az/El in degrees

    '''
    
    ''' load input sittings: '''
    inp = inp_set(cfg)
    inp = inp.get_section('all')
    
    ''' mkdir '_out' if non existend '''
    if inp['out_path'][0] == '/': # abs path?
        if not os.path.isdir(os.path.join(inp['out_path'])):
            os.makedirs(os.path.join(inp['out_path']))
    else:
        abs_path = os.path.dirname(inspect.getfile(inspect.currentframe()))
        if not os.path.isdir(os.path.join(abs_path, inp['out_path'])):
            os.makedirs(os.path.join(abs_path, inp['out_path']))
    
    const = constants()
    
    # known spacecraft:
    spacecraft = {'MEX':'S', 'VEX':'S', 'ROSETTA':'S', 'GAIA':'S', 'HER':'S',\
                  'RA':'R'}
    # GNSS (GLONASS)
    spacecraft.update({'PR{:02d}'.format(i):'G' for i in range(40)})
    spacecraft.update({'PG{:02d}'.format(i):'G' for i in range(40)})

    source = source.upper()
    sou_type = spacecraft[source]
    
    ''' stations '''
    if type(stations) != type([]):        
        stations = [stations]
    
    ''' fake obs '''
    ob = obs(['DUMMY'],'DUMMY','C')
    ob.addScan(date_t_start, t_step, stop=date_t_stop)
    print ob
    ''' update/(down)load eops, meteo and iono data '''
    if internet_on(): # check internet connection
        try:
            doup(inp['do_trp_calc'], inp['do_ion_calc'], \
                 inp['cat_eop'], inp['meteo_cat'], inp['ion_cat'],\
                 date_t_start, date_t_stop, inp['iono_model'])
        except Exception, err:
            print str(err)
            print 'catalogue updates failed'
    
    ''' load cats '''
    _, sta, eops = load_cats(inp, source, sou_type, stations, date_t_start)
    
    ''' calculate site positions in geodetic coordinate frame
    + transformation matrix VW from VEN to the Earth-fixed coordinate frame '''
    for st in sta:
        st.geodetic(const) 
    
    # load ephemeris
    eph = load_sc_eph(sou_type, source, date_t_start, date_t_stop, inp)
        
    #%% 
    ''' actual pointings '''
    
    pointingsJ2000 = [] # RA/Decs at J2000
    pointingsDate = []  # apparent RA/Decs (precessed and nutated to date)
    azels = [] # azimuth/elevations
    
    mjd_start = mjuliandate(ob.tstamps[0].year,\
                                   ob.tstamps[0].month,ob.tstamps[0].day)
    
    for tstamp in ob.tstamps:
        ''' set dates: '''
        mjd = mjuliandate(tstamp.year, tstamp.month, tstamp.day)
        dd = mjd - mjd_start
        UTC = (tstamp.hour + tstamp.minute/60.0 + tstamp.second/3600.0)/24.0
        JD = mjd + 2400000.5
    
        ''' compute tai & tt '''
        TAI, TT = taitime(mjd, UTC)
    
        ''' interpolate eops to tstamp '''
        UT1, eop_int = eop_iers(mjd, UTC, eops)
    
        ''' compute coordinate time fraction of CT day at 1st observing site '''
        CT, dTAIdCT = t_eph(JD, UT1, TT, sta[0].lon_gcen, sta[0].u, sta[0].v)
    
        ''' BCRS state vectors of celestial bodies at JD+CT, [m, m/s]: '''
        ## Earth:
        rrd = pleph(JD+CT, 3, 12, inp['jpl_eph'])
        earth = np.reshape(np.asarray(rrd), (3,2), 'F') * 1e3
        # Earth's acceleration in m/s**2:
        v_plus = np.array(pleph(JD+CT+1.0/86400.0, 3, 12, inp['jpl_eph'])[3:])
        v_minus = np.array(pleph(JD+CT-1.0/86400.0, 3, 12, inp['jpl_eph'])[3:])
        a = (v_plus - v_minus)*1e3 / 2.0
        a = np.array(np.matrix(a).T)
        earth = np.hstack((earth, a))
        ## Sun:
        rrd = pleph(JD+CT, 11, 12, inp['jpl_eph'])
        sun = np.reshape(np.asarray(rrd), (3,2), 'F') * 1e3
        ## Moon:
        rrd = pleph(JD+CT, 10, 12, inp['jpl_eph'])
        moon = np.reshape(np.asarray(rrd), (3,2), 'F') * 1e3
    
        ''' rotation matrix IERS '''
        r2000 = ter2cel(tstamp, eop_int, dTAIdCT, 'iau2000')
    
        ''' displacements due to geophysical effects '''
        for st in sta:
            if st.name == 'GEOCENTR' or st.name == 'RA':
                continue
            st = dehanttideinel(st, tstamp, earth, sun, moon, r2000)
            st = hardisp(st, tstamp, r2000)
            st = poletide(st, tstamp, eop_int, r2000)
        
        ''' add up geophysical corrections and convert sta state to J2000 '''
        for st in sta:
            if st.name == 'GEOCENTR' or st.name == 'RA':
                continue
            st.j2000gp(r2000)
        
        ''' calculate [Ra, Dec] for cross-eyed scheduling: '''
        pnt_J2000_sta = []
        pnt_Date_sta = []
        azel_sta = []
        for st in sta:
            # J2000 RA/Dec:
#            print eph.CT[0]<= CT + dd <= eph.CT[-1]
            _, ra, dec = st.LT_radec_bc(eph.bcrs[0], eph.CT, JD, CT+dd, \
                                        inp['jpl_eph'])
            pnt_J2000_sta.append([ra, dec])
    #        print st.name, ra, dec
            # RA/Dec to date:
            xyz2000 = sph2cart(np.array([1.0, dec, ra]))
            rDate = iau_PNM00A(JD, TT)
            xyzDate = np.dot(rDate, xyz2000)
            dec, ra = cart2sph(xyzDate)[1:]
            if ra < 0: ra += 2.0*np.pi
    #        print st.name, ra, dec
            pnt_Date_sta.append([ra, dec])
            if st.name == 'GEOCENTR' or st.name == 'RA':
                azel_sta.append([0, 0])
            else:
                az, el, _ = st.AzEl2(eph.gcrs, eph.UT, JD, UTC,\
                                         r2000[:,:,0], inp['jpl_eph'])
                azel_sta.append([az*180/np.pi, el*180/np.pi])
    
        pointingsJ2000.append(pnt_J2000_sta)
        pointingsDate.append(pnt_Date_sta)
        azels.append(azel_sta)
        
    
    pointingsJ2000 = np.array(pointingsJ2000)
    pointingsDate = np.array(pointingsDate)
    azels = np.array(azels)
    
    #%% 
    ''' save to files '''
    if output:
        stations_short = shname(stations, inp['shnames_cat'], \
                                            inp['shnames_cat_igs'])
        
        date_string = date_t_start.strftime("%y%m%d")
        
        for jj, stash in enumerate(stations_short):
            print 'Outputting {:s}'.format(stash)
            # save AzEls to human-readable text-files
            with open(inp['out_path']+'/azel.'+source.lower()+'.'+\
                      date_string+'.'+stash.lower(),'w') as f:
                for ii, tstamp in enumerate(ob.tstamps):
                    line = str(tstamp)+ '  '
                    az = azels[ii,jj,0]
                    el = azels[ii,jj,1]
                    azel = np.hstack((az, el))
                    azel_str = \
                       'az = {:010.6f}  el = {:10.6f}\n'\
                       .format(*azel) # '{:06.2f} {:6.2f}\n'\
                    line += azel_str
                    f.write(line)
            # save AzEls to human-readable text-files in GreenBank format
            if stash.lower()=='gb':
                with open(inp['out_path']+'/azel.gbt.'+source.lower()+'.'+\
                          date_string+'.'+stash.lower(),'w') as f:
                    f.write('# {:s} tracking table (angles in degrees)\n'.\
                            format(source.upper()))
                    f.write('format=ephemeris\n')
                    f.write('name={:s}\n'.format(source.upper()))
                    f.write('coordmode=azel\n')
                    f.write('head=date    utc        az          el\n')
                    for ii, tstamp in enumerate(ob.tstamps):
                        line = tstamp.strftime('%Y-%m-%d %H:%M:%S')+ '   '
                        az = azels[ii,jj,0]
                        el = azels[ii,jj,1]
                        azel = np.hstack((az, el))
                        azel_str = '{:10.6f}  {:10.6f}\n'.format(*azel)
                        line += azel_str
                        f.write(line)
            # save pointingsJ2000 (J2000 RA/Decs) to human-readable text-files
            with open(inp['out_path']+'/pointing.J2000.'+source.lower()+'.'+\
                      date_string+'.'+stash.lower(),'w') as f:
                for ii, tstamp in enumerate(ob.tstamps):
                    line = tstamp.strftime('%Y-%m-%dT%H:%M:%S') + '    '
                    ra = Angle(pointingsJ2000[ii,jj,0], unit=units.rad)
                    dec = Angle(pointingsJ2000[ii,jj,1], unit=units.rad)
                    radec = np.hstack((ra.hms, dec.dms))
                    radec[4:] = abs(radec[4:]) # minus doesn't belong to everyone..
                    radec_str = \
              'ra = {:02.0f}h{:02.0f}m{:010.7f}s  dec = {:-3.0f}d{:02.0f}\'{:010.7f}\"\n'\
                       .format(*radec)
                    line += radec_str
                    f.write(line)
            # save pointingsJ2000 to Giuseppe-friendly text-files
            fout = os.path.join(inp['out_path'], ''.join((source.lower(), '_', \
                                           stash.lower(), '_', date_string, '.txt')) )
            with open(fout, 'w') as f:
                for ii, tstamp in enumerate(ob.tstamps):
                    line = ''.join((' source=\'', tstamp.strftime("%H%M%S"), '\' \t '))
                    ra = Angle(pointingsJ2000[ii,jj,0], unit=units.rad)
                    dec = Angle(pointingsJ2000[ii,jj,1], unit=units.rad)
                    radec = np.hstack((ra.hms, dec.dms))
                    radec[4:] = abs(radec[4:]) # minus doesn't belong to everyone..
                    radec_str = 'ra={:02.0f}:{:02.0f}:{:05.2f} \t '.format(*radec[0:3]) 
                    radec_str += 'dec={:+.0f}:{:02.0f}:{:05.2f} \t '.format(*radec[3:])
                    radec_str +=  'equinox=\'j2000\' /\n'
                    line += radec_str
                    f.write(line)
            # save pointingsDate (apparent RA/Decs) to human-readable text-files
            with open(inp['out_path']+'/pointing.apparent.'+source.lower()+'.'+\
                      date_string+'.'+stash.lower(),'w') as f:
                for ii, tstamp in enumerate(ob.tstamps):
                    line = tstamp.strftime('%Y-%m-%dT%H:%M:%S') + '    '
                    ra = Angle(pointingsDate[ii,jj,0], unit=units.rad)
                    dec = Angle(pointingsDate[ii,jj,1], unit=units.rad)
                    radec = np.hstack((ra.hms, dec.dms))
        
                    radec[4:] = abs(radec[4:]) # minus doesn't belong to everyone..
                    radec_str = \
              'ra = {:02.0f}h{:02.0f}m{:010.7f}s  dec = {:-3.0f}d{:02.0f}\'{:010.7f}\"\n'\
                       .format(*radec)
                    line += radec_str
                    f.write(line)
                    
    return ob.tstamps, pointingsJ2000, pointingsDate, azels
