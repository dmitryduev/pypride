# -*- coding: utf-8 -*-
"""
Produce pointings on a spacecraft for a set of telescopes 
for a given time period.

For a list of supported spacecraft, see function 
pointings definition in vintlib

Function pointings outputs an array of datetime objects with obs epochs,
RA/Decs for J2000 and date in radians and Az/El in degrees

Created on Fri Jan 17 22:52:40 2014

@author: Dmitry A. Duev
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
    
    source = 'RA'
    
#    stations = ['GBT-VLBA', 'HN-VLBA', 'NL-VLBA']
#    stations = ['GBT-VLBA']
    stations = ['ONSALA60']
#    stations = ['EFLSBERG', 'SVETLOE', 'ZELENCHK', 'WETTZELL']
    
    ''' time slot '''
#    date_t_start = datetime.datetime(2013,12,28,23,57,0)
#    date_t_stop  = datetime.datetime(2013,12,29,0,1,0)
#    t_step = 10 # seconds
    
    date_t_start = datetime.datetime(2015,10,8,11,0,0)
    date_t_stop  = datetime.datetime(2015,10,8,14,0,0)
#    date_t_start = datetime.datetime(2015,10,10,7,0,0)
#    date_t_stop  = datetime.datetime(2015,10,10,9,0,0)
#    date_t_start = datetime.datetime(2015,10,24,14,0,0)
#    date_t_stop  = datetime.datetime(2015,10,24,15,0,0)
#    date_t_start = datetime.datetime(2015,10,25,14,0,0)
#    date_t_stop  = datetime.datetime(2015,10,25,15,0,0)
#    date_t_start = datetime.datetime(2015,10,31,22,0,0)
#    date_t_stop  = datetime.datetime(2015,10,31,23,0,0)
#    date_t_start = datetime.datetime(2015,11,2,18,0,0)
#    date_t_stop  = datetime.datetime(2015,11,2,19,0,0)
    t_step = 60 # seconds
    
    # output txt dumps?
    output = True
    
    ''' run calculation '''
    # output is placed in cfg['out_path']
    times, RaDecJ2000, RaDecDate, AzEl = \
        pointings(source, stations, date_t_start, date_t_stop, t_step, \
                  cfg, output)