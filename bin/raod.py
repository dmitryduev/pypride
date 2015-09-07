#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 18:52:05 2015

@author: oasis

Create and maintain a mirror of RA OD data FTP-server
"""

from ftplib import FTP
import datetime
import time

import argparse
import os
import numpy as np

import logging

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
    

#%%
if __name__=='__main__':
    # create parser
    parser = argparse.ArgumentParser(prog='python raod.py',\
                formatter_class=argparse.RawDescriptionHelpFormatter,\
                description='Create and maintain a mirror of RA OD data FTP-server.')
    # optional arguments
    parser.add_argument('-o', '--out', type=str,\
                        default='/jop75_0/duev/radioastron_oddata',\
                        help='local output folder. '+\
                          'defaults to /jop75_0/duev/radioastron_oddata')
    parser.add_argument('-n', '--nodown', action='store_true',\
                        help='do not download updates, only check database')
    # positional argument
#    parser.add_argument('vexfile', type=str, help='input vex-file')
    args = parser.parse_args()
    
    # set-up log
#    path = os.path.dirname(os.path.abspath(__file__))
    path = args.out
    logging.basicConfig(filename=os.path.join(path, 'raod.log'),\
                        format='%(asctime)s %(message)s',\
                        datefmt='%d/%m/%Y %H:%M:%S',\
                        level=logging.DEBUG)
    
    # updates wanted?
    if not args.nodown:
        try:
            ftp = FTP('webinet.asc.rssi.ru')
            ftp.login(user='radik', passwd='Zap069_qq')
#            ftp.cwd('radioastron/oddata/reconstr/new_test')
            ftp.cwd('radioastron/oddata/reconstr/')
    #        print ftp.nlst()
            for f in ftp.nlst():
                resp = ftp.sendcmd('MDTM {:s}'.format(f))
                if resp[0] == '2':
                    timestamp = time.mktime(datetime.datetime.strptime(resp[4:18],\
                                                    "%Y%m%d%H%M%S").timetuple())
                else: 
                    timestamp = None
                if f[-4:] == '.scf' and (f not in os.listdir(args.out) or \
                        (timestamp is not None and \
                            os.path.getmtime(os.path.join(args.out,f)) < timestamp)):
#                    print 'fetching {:s}...'.format(f)
                    logging.info('fetching {:s}...'.format(f))
                    ftp.retrbinary('RETR {:s}'.format(f), \
                                        open(os.path.join(args.out,f), 'wb').write)
                    if resp[0] == '2':
                        os.utime(os.path.join(args.out,f), (timestamp, timestamp))
    #                print datetime.datetime.utcfromtimestamp(\
    #                        os.path.getmtime(os.path.join(args.out,f)))
            ftp.quit()
            
        except Exception, err:
#            print str(err)
            logging.info('Error: {:s}'.format(str(err)))
        
    ''' maintain database with boundaries '''
    f_b = os.path.join(args.out, 'boundaries.raw')
   
    ls = os.listdir(args.out)
    ls = [l for l in ls if l[-4:]=='.scf']

    skip = ['RA150502_1700_v02.scf']
    
    # does it exist already?
    if not os.path.isfile(f_b):
#        print 'Boundaries database file {:s} does not exist, creating'.\
#                    format(f_b)
        logging.info('Boundaries database file {:s} does not exist, creating'.\
                    format(f_b))
        entries = {}
        for f_ls in ls:
            if f_ls not in skip:
                try:
#                    print 'Adding {:s} to the database'.format(f_ls)
                    logging.info('Adding {:s} to the database'.format(f_ls))
                    orb_file = load_scf(os.path.join(args.out, f_ls))
                    t_start = datetime.datetime.strftime(\
                                datetime.datetime(*map(int,orb_file[0,0:6])),\
                                '%Y-%m-%dT%H:%M:%S')
                    t_stop  = datetime.datetime.strftime(\
                                datetime.datetime(*map(int,orb_file[-1,0:6])),\
                                '%Y-%m-%dT%H:%M:%S')
                    entries[f_ls] = [t_start, t_stop]
                except Exception, err:
#                    print '{:s}:'.format(f_ls), str(err)
                    logging.info('{:s}:'.format(f_ls), str(err))
                    pass
        if len(ls)>0:
            with open(f_b, 'w') as f:
                for key in sorted(entries.keys()):
                    f.write('{:35s}  {:s}  {:s}\n'.\
                                format(key, entries[key][0],entries[key][1]))
                    
    else:
        with open(f_b, 'r') as f:
            f_lines = f.readlines()
        f_lines = [l.split() for l in f_lines]
        entries = {l[0]:[l[1], l[2]] for l in f_lines}
        ls_new = [l for l in ls if l not in entries.keys()]
        for f_ls in ls_new:
            if f_ls not in skip:
                try:
#                    print 'Adding {:s} to the database'.format(f_ls)
                    logging.info('Adding {:s} to the database'.format(f_ls))
                    orb_file = load_scf(os.path.join(args.out, f_ls))
                    t_start = datetime.datetime.strftime(\
                                datetime.datetime(*map(int,orb_file[0,0:6])),\
                                '%Y-%m-%dT%H:%M:%S')
                    t_stop  = datetime.datetime.strftime(\
                                datetime.datetime(*map(int,orb_file[-1,0:6])),\
                                '%Y-%m-%dT%H:%M:%S')
                    entries[f_ls] = [t_start, t_stop]
                except Exception, err:
#                    print '{:s}:'.format(f_ls), str(err)
                    logging.info('{:s}:'.format(f_ls), str(err))
                    pass
        # got something?
        if len(ls_new)>0:
            with open(f_b, 'w') as f:
                for key in sorted(entries.keys()):
                    f.write('{:35s}  {:s}  {:s}\n'.\
                                format(key, entries[key][0],entries[key][1]))