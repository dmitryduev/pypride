# -*- coding: utf-8 -*-
"""
Produce pointings on a spacecraft for a set of telescopes 
for a given time period.

For a list of supported spacecraft, see function 
pointings definition in vintlib

Function pointings outputs an array of datetime objects with obs epochs,
RA/Decs for J2000 and date in radians and Az/El in degrees

Created on Fri Jan 17 22:52:40 2014

@author: oasis
"""

from pypride.vintlib import pointings
import datetime

#%% 
'''
#==============================================================================
# example set-up
#==============================================================================
'''
if __name__=='__main__':
    cfg = 'inp.cfg'
    
    source = 'MEX'
    
    stations = ['BADARY', 'CEDUNA']
    
    ''' time slot '''
    date_t_start = datetime.datetime(2013,12,28,23,57,0)
    date_t_stop  = datetime.datetime(2013,12,29,0,1,0)
    t_step = 10 # seconds
    
    # output txt dumps?
    output = False
    
    ''' run calculation '''
    # output is placed in cfg['out_path']
    times, RaDecJ2000, RaDecDate, AzEl = \
        pointings(source, stations, date_t_start, date_t_stop, t_step, \
                  cfg, output)