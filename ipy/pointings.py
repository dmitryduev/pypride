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

'''
#==============================================================================
# example set-up
#==============================================================================
'''
if __name__ == '__main__':
    cfg = '/Users/dmitryduev/anaconda/lib/python2.7/site-packages/pypride/inp.cfg'
    
    source = 'RA'
    
    # stations = ['GBT-VLBA', 'HN-VLBA', 'NL-VLBA']
    # stations = ['WETTZELL', 'WETTZ13N']
    # stations = ['EFLSBERG', 'WETTZELL', 'ZELENCHK', 'SVETLOE']
    # stations = ['ONSALA60']
    stations = ['MATERA']
    # stations = ['HARTRAO']
    # stations = ['EFLSBERG', 'HARTRAO', 'WETTZELL', 'ZELENCHK', 'SVETLOE']
    # stations = ['YARRA12M']
    
    ''' time slot '''
    # date_t_start = datetime.datetime(2016, 6, 14, 18, 5, 0)
    # date_t_stop = datetime.datetime(2016, 6, 14, 19, 55, 0)
    # date_t_start = datetime.datetime(2016, 6, 15, 10, 5, 0)
    # date_t_stop = datetime.datetime(2016, 6, 15, 16, 58, 0)
    date_t_start = datetime.datetime(2016, 6, 24, 0, 5, 0)
    date_t_stop = datetime.datetime(2016, 6, 25, 20, 55, 0)
    # date_t_start = datetime.datetime(2016, 6, 24, 4, 0, 0)
    # date_t_stop = datetime.datetime(2016, 6, 24, 5, 0, 0)
    t_step_in = 60  # seconds

    # interpolate to a denser frig if t_step_out < t_step_in:
    t_step_out = 1  # seconds
    scan_len = 600  # split time range into 'scans'
    
    # output txt dumps?
    output = True
    
    ''' run calculation '''
    # output is placed in cfg['out_path']
    times, RaDecJ2000, RaDecDate, AzEl = \
        pointings(source, stations, date_t_start, date_t_stop,
                  t_step_in, t_step_out, scan_len, cfg, output)
