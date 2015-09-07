# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 14:46:37 2015

Convert ESA deep space s/c 1-way ramp table to pyp format
Output is placedin the same folder as input

@author: oasis
"""
import argparse
import os
from datetime import datetime, timedelta

#esa_in = '/Users/oasis/_jive/Mex/gr035/tati-1way-Doppler/PSA_test_ingress.txt'

if __name__ == '__main__':
    '''
        If not imported, but run from the command line
    '''
    # create parser
    parser = argparse.ArgumentParser(prog='python ramp1w_esa2pyp.py', \
                formatter_class=argparse.RawDescriptionHelpFormatter,\
                description='Convert ESA deep space s/c 1-way ramp table to pyp format.')
    # optional arguments
    parser.add_argument('-o', '--outfile', type=str, default='ramp1w.sc',\
                        help='output file name (location is the same as for esa_in)')
    # positional argument
    parser.add_argument('esa_in', type=str, help='input ESA file')    
    args = parser.parse_args()
    
    out_dir = os.path.dirname(args.esa_in)
    
    if args.outfile:
        pyp_out = args.outfile
    
    ramp = []
    with open(args.esa_in, 'r') as fin, \
            open(os.path.join(out_dir, pyp_out), 'w') as fout:
        f_lines = fin.readlines()
        f_lines = [l.split() for l in f_lines if len(l)>10]
        
        for i, l in enumerate(f_lines[:-1]):
            t_start = datetime.strptime(l[5], "%Y-%m-%dT%H:%M:%S.%f")
            t_stop = datetime.strptime(f_lines[i+1][5], "%Y-%m-%dT%H:%M:%S.%f")
            ramp.append([t_start, t_stop, l[6], l[7]])

            fout.write('{:s}  {:s}  {:22.16e}  {:23.16e}\n'.format(\
                        datetime.strftime(t_start, "%Y-%m-%d %H:%M:%S.%f"),\
                        datetime.strftime(t_stop, "%Y-%m-%d %H:%M:%S.%f"),\
                        float(l[6]), float(l[7]) ))

        # last interval:
        t_start = datetime.strptime(f_lines[-1][5], "%Y-%m-%dT%H:%M:%S.%f")
        fout.write('{:s}  {:s}  {:22.16e}  {:23.16e}\n'.format(\
                        datetime.strftime(t_start, "%Y-%m-%d %H:%M:%S.%f"),\
                        datetime.strftime(t_start + timedelta(seconds=0.1), \
                                            "%Y-%m-%d %H:%M:%S.%f"),\
                        float(l[6]), float(l[7]) ))
        ramp.append([t_start, t_start + timedelta(seconds=0.1), l[6], l[7]])
        