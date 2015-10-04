#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 10 13:48:36 2015

@author: oasis

Convert RadioAstron's blt-file to scf-file(s).
"""

import argparse
import os
import datetime

#%%
def make_scf_header(t_start, t_stop):
    scf_header = '''CCSDS_OEM_VERS = 2.0
CREATION DATE = {:s}
ORIGINATOR    = KIAM
META_START
OBJECT_NAME   = RADIOASTRON
OBJECT_ID     = 2011-37A
CENTER_NAME   = Earth
REF_FRAME     = EME2000
TIME_SYSTEM   = UTC
START_TIME    = {:s}
STOP_TIME     = {:s}
META_STOP
'''.format(datetime.datetime.strftime(datetime.datetime.now(), \
                                      '%Y-%m-%dT%H:%M:%S'),\
            datetime.datetime.strftime(t_start,'%Y-%m-%dT%H:%M:%S.%f')[:-3],\
            datetime.datetime.strftime(t_stop, '%Y-%m-%dT%H:%M:%S.%f')[:-3])
            
    return scf_header
    

if __name__ == '__main__':
    '''
        If not imported, but run from the command line
    '''
    # create parser
    parser = argparse.ArgumentParser(prog='python blt2scf.py', \
                formatter_class=argparse.RawDescriptionHelpFormatter,\
                description='Convert blt-file to scf-file(s).')
    # optional arguments
    parser.add_argument('-s', '--split', action='store_true',
                        help='split to multiple files if it goes overnight')
    parser.add_argument('-t', '--dmtutc', default=3,
                        help='Moscow Time - UTC in h. defaults to 3 h.')
    # positional argument
    parser.add_argument('blt', type=str, help='input blt-file')    
    args = parser.parse_args()
    
    with open(args.blt, 'r') as f:
        f_lines = f.readlines()
    path = os.path.dirname(args.blt)
        
    header = f_lines.pop(0).split()
#    day_start = datetime.datetime.strptime(header[0], '%Y%m%d')
    t_start = datetime.datetime.strptime(header[0]+header[1], '%Y%m%d%H%M%S.%f')
    t_start = t_start - datetime.timedelta(hours=int(args.dmtutc)) # Moscow time to UTC
    day_start = datetime.datetime(t_start.year, t_start.month, t_start.day)
    t_step = datetime.timedelta(seconds=float(header[-1]))
    t_stop = t_start+t_step*(len(f_lines)-1)
#    print t_start, t_step
    
    if not args.split: # make one big file
        scf_header = make_scf_header(t_start, t_stop)
        
        f_out = 'RA{:s}j.scf'.format(datetime.datetime.strftime(t_start,'%y%m%d'))
        with open(os.path.join(path, f_out), 'w') as f:
            f.write('{:s}\n'.format(scf_header))
            for ii, line in enumerate(f_lines):
                rv = map(float, line.split()[0:6])
                t = t_start + t_step*ii
                f.write('{:s} '.format(datetime.datetime.strftime(t,\
                                        '%Y-%m-%dT%H:%M:%S.%f')[:-3]))
                f.write('{:14.6f} {:14.6f} {:14.6f} {:14.9f} {:14.9f} {:14.9f}\n'.\
                        format(*rv))
    
    else: # split on a daily basis
        dd = (datetime.datetime(t_stop.year, t_stop.month, t_stop.day) - \
             datetime.datetime(t_start.year, t_start.month, t_start.day)).days
        times = [t_start + t_step*ii for ii in range(len(f_lines))]
#        print dd
        
        for d in range(dd+1):
            # current day indices:
            day = [ii for ii,tt in enumerate(times) if (datetime.datetime(tt.year,\
                    tt.month, tt.day)-day_start).days==d]
            f_out = 'RA{:s}j.scf'.format(datetime.datetime.strftime(t_start+\
                                        datetime.timedelta(days=d),'%y%m%d'))
            scf_header = make_scf_header(times[day[0]], times[day[-1]])
            with open(os.path.join(path, f_out), 'w') as f:
                f.write('{:s}\n'.format(scf_header))
                for ii in day:
                    rv = map(float, f_lines[ii].split()[0:6])
                    t = t_start + t_step*ii
                    f.write('{:s} '.format(datetime.datetime.strftime(t,\
                                            '%Y-%m-%dT%H:%M:%S.%f')[:-3]))
                    f.write('{:14.6f} {:14.6f} {:14.6f} {:14.9f} {:14.9f} {:14.9f}\n'.\
                            format(*rv))
                    