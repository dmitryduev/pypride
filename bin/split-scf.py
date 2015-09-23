#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 10 13:48:36 2015

@author: oasis
"""

import argparse
import os
import datetime
import numpy as np
import paramiko

'''
#==============================================================================
#  Load eph in IPM-style scf format       
#==============================================================================
'''
def load_scf(f_lines):
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
    
def strictly_increasing(L):
    return all(x<y for x, y in zip(L, L[1:]))

#%%
if __name__ == '__main__':
    '''
        If not imported, but run from the command line
    '''
    # create parser
    parser = argparse.ArgumentParser(prog='split-scf.py', \
                formatter_class=argparse.RawDescriptionHelpFormatter,\
                description='Split scf-file to multiple files if it goes overnight.')
    # optional arguments
    parser.add_argument('-o', '--out', type=str,\
                        default='sc_eph/raw_radioastron',\
                        help='output folder. '+\
                          'defaults to <pypride_path>/sc_eph/raw_radioastron')
    
    # positional argument
    parser.add_argument('scf', type=str, help='input scf-file.'+\
                        'could be accessed over ssh.'+\
                      'provide -u and -p parameters if there\'s no key.')
    parser.add_argument('-u', '--user', type=str, default=None, help='ssh username')
    parser.add_argument('-p', '--pwd', type=str, default=None, help = 'ssh password')
    args = parser.parse_args()
    
    
    # if it's on the network:
    if ':' in args.scf:
        if args.user is None:
            raise Exception('username must be provided for ssh access.')
        client = paramiko.SSHClient()
        client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        host = args.scf[:args.scf.index(':')]
        if args.pwd is None:
            client.connect(hostname=host, username=args.user)
        else:
            client.connect(hostname=host, username=args.user, password=args.pwd)
        # read the database over ssh:
        sftp_client = client.open_sftp()
        path = args.scf[args.scf.index(':')+1:]
        fname = os.path.basename(path)
        with sftp_client.open(path) as f:
            print 'fetching {:s} from {:s}'.format(fname, host)
            f_lines = f.readlines()          
        # close connection
        client.close()
    # stored locally?
    else:
        with open(args.scf, 'r') as f:
            f_lines = f.readlines()
    path = args.out
    import pypride
    path = os.path.join(os.path.dirname(pypride.__file__), path)

    
    # load input scf-file:
    scf_raw = load_scf(f_lines)
    
    # check for bad entries:
    xx = [float(entry[5]) for entry in scf_raw]
    for i,x in enumerate(xx):
        if x>=60.0:
            print 'Bad entry \'{:f}\' on line {:d}'.format(x,i)

    t_raw = np.array([datetime.datetime(*map(int, entry[:5]), \
                      second=int(np.floor(float(entry[5]))), \
                      microsecond=int(np.round((entry[5]-np.fix(entry[5]))*1e6,6))) \
                                                         for entry in scf_raw])
    
    # check for *missed* overlaps:
    t_sec = [(ti-t_raw[0]).total_seconds() for ti in t_raw]
#    print strictly_increasing(t_sec)
    if not strictly_increasing(t_sec):
        raise Exception('Check input file - epochs are not monotonically increasing')
    
    # starting day
    t_0 = datetime.datetime(*map(int, scf_raw[0,:3]))
    
    dd = np.array(sorted(set([(t-t_0).days for t in t_raw])))
    
    # iterate over days:
    for d in dd:
        day_start = t_0 + datetime.timedelta(days=d)
        day_stop  = t_0 + datetime.timedelta(days=d+1)
        # 'cut' the current day:
        t_d_ind = np.logical_and(t_raw>=day_start, t_raw<day_stop)
        t_d = t_raw[t_d_ind]
        scf_raw_d = scf_raw[t_d_ind, :]

        t_start, t_stop = t_d[0], t_d[-1]
        f_out = 'RA{:s}j.scf'.format(datetime.datetime.strftime(day_start,'%y%m%d'))
#        print t_start, t_stop, f_out
        
        # make header:
        scf_header = make_scf_header(t_start, t_stop)
        with open(os.path.join(path, f_out), 'w') as f:
            print 'outputting {:s}'.format(f_out)
            f.write('{:s}\n'.format(scf_header))
            for ti, entry in zip (t_d, scf_raw_d):
                f.write('{:s} '.format(datetime.datetime.strftime(ti,\
                                            '%Y-%m-%dT%H:%M:%S.%f')[:-3]))
                f.write('{:14.6f} {:14.6f} {:14.6f} {:14.9f} {:14.9f} {:14.9f}\n'.\
                            format(*entry[6:]))
    