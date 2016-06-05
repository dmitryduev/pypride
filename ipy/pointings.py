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
    # stations = ['GBT-VLBA']
    # stations = ['ONSALA60']
    # stations = ['ZELENCHK', 'SVETLOE']
    # stations = ['EFLSBERG', 'WETTZELL']
    # stations = ['WETTZELL', 'WETTZ13N']
    # stations = ['EFLSBERG', 'WETTZELL', 'ZELENCHK', 'SVETLOE']
    # stations = ['ONSALA60']
    stations = ['HARTRAO']
    
    ''' time slot '''
    # date_t_start = datetime.datetime(2016, 5, 28, 10, 0, 5)
    # date_t_stop = datetime.datetime(2016, 5, 28, 12, 10, 5)
    # date_t_start = datetime.datetime(2016, 5, 29, 11, 0, 5)
    # date_t_stop = datetime.datetime(2016, 5, 29, 13, 10, 5)
    # date_t_start = datetime.datetime(2016, 5, 29, 20, 0, 0)
    # date_t_stop = datetime.datetime(2016, 5, 29, 22, 10, 0)
    date_t_start = datetime.datetime(2016, 6, 4, 0, 5, 0)
    # date_t_start = datetime.datetime(2016, 6, 4, 23, 58, 0)
    date_t_stop = datetime.datetime(2016, 6, 7, 23, 55, 0)
    t_step = 10  # seconds
    
    # output txt dumps?
    output = True
    
    ''' run calculation '''
    # output is placed in cfg['out_path']
    times, RaDecJ2000, RaDecDate, AzEl = \
        pointings(source, stations, date_t_start, date_t_stop, t_step, cfg, output)
