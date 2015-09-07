#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 13:41:50 2015

@author: oasis

Convert scf-file to bsp-file for creating a SPICE-kernel.
"""

import numpy as np
import argparse
import os
import datetime

#%%
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

if __name__ == '__main__':
    '''
        If not imported, but run from the command line
    '''
    # create parser
    parser = argparse.ArgumentParser(prog='python scf2bsp.py', \
                formatter_class=argparse.RawDescriptionHelpFormatter,\
                description='Convert scf-file to bsp-file for creating a SPICE-kernel.')
    # positional argument
    parser.add_argument('scf', type=str, help='input scf-file')    
    args = parser.parse_args()
    

    path = os.path.dirname(args.scf)
        
    eph = load_scf(args.scf)
    
    t_start = datetime.datetime(*map(int, eph[0,0:3]))
    
    f_out = 'RA{:s}j.bsp'.format(datetime.datetime.strftime(t_start,'%y%m%d'))
    with open(os.path.join(path, f_out), 'w') as f:
        for entry in eph:
            line = '{:4.0f}-{:02.0f}-{:02.0f} {:02.0f}:{:02.0f}:{:06.3f}'.\
                    format(*entry[0:6])
            line += ',{:.6f},{:.6f},{:.6f},{:.10f},{:.10f},{:.10f},\n'\
                    .format(*entry[6:12])
            f.write(line)
    