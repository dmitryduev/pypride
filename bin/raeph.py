#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 14:40:58 2015

@author: oasis

Create raw RA eph files for use with PYP

KIAM usually provides more than one eph for a given date,
thus some manual 'picking' must be done

"""

import datetime
import argparse
import os
import numpy as np
from pypride.vex import Vex
import paramiko

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

def flatten(x):
    result = []
    for el in x:
        if hasattr(el, "__iter__") and not isinstance(el, basestring):
            result.extend(flatten(el))
        else:
            result.append(el)
    return result


#%%
if __name__=='__main__':
    # create parser
    parser = argparse.ArgumentParser(prog='raeph.py',\
                formatter_class=argparse.RawDescriptionHelpFormatter,\
                description='Create raw RA eph files for use with PYPRIDE.')
    # optional arguments
    parser.add_argument('-s', '--scf', type=str,\
                        default='jop75:/jop75_0/duev/radioastron_oddata',\
                        help='RA OD database folder. could be accessed over ssh.'+\
                          'provide -u and -p parameters if there\'s no key.'+\
                          'defaults to jop75:/jop75_0/duev/radioastron_oddata')
    parser.add_argument('-u', '--user', type=str, default=None, help='ssh username')
    parser.add_argument('-p', '--pwd', type=str, default=None, help = 'ssh password')
#    parser.add_argument('-n', '--nodown', action='store_true',\
#                        help='do not download updates, only check database')
    # positional argument
    parser.add_argument('vexfile', type=str, help='input vex-file')
    args = parser.parse_args()
    
    ''' vex-file '''
    vex = Vex(args.vexfile)

    # stations
    stations = {s:vex['STATION'][s]['ANTENNA'] for s in vex['STATION']}
    # sources
    sources_obs = list(set(flatten([vex['SCHED'][s].getall('source') \
                                        for s in vex['SCHED']])))
    sou = ['RA']

    for sousou in sou:
        # RA was observed:
        if sousou in sources_obs: # S/C was actually observed?
            t_begin, t_end = None, None
            for s in vex['SCHED']:
                sou_scan = vex['SCHED'][s].getall('source')
                if sousou in [ss.upper() for ss in sou_scan]:
                    t_start = datetime.datetime.strptime(\
                                        vex['SCHED'][s]['start'], \
                                                  '%Yy%jd%Hh%Mm%Ss')
                    if t_begin is None:
                        t_begin = t_start
                    N_sec = int(vex['SCHED'][s].\
                                getall('station')[0][2].split()[0])
                    t_end = t_start + datetime.timedelta(seconds=N_sec)
            # make the ephem
            if t_begin is not None:
                print 'experiment boundaries:', t_begin, t_end

        # RadioAstron observed itself:
        if sousou in stations.values(): # RA observing scheduled?
            t_begin, t_end = None, None
            for s in vex['SCHED']:
                sta = [st for st in vex['SCHED'][s].getall('station') \
                        if st[0].upper()==sousou]
                if len(sta)>0:
                    sta = sta[0]
                    t_start = datetime.datetime.strptime(\
                                        vex['SCHED'][s]['start'], \
                                                  '%Yy%jd%Hh%Mm%Ss')
                    if t_begin is None:
                        t_begin = t_start
                    
                    N_sec = int(sta[2].split()[0])
                    t_end = t_start + datetime.timedelta(seconds=N_sec)
            # make the ephem
            if t_begin is not None:
                print 'experiment boundaries:', t_begin, t_end

    ''' database + pick the files'''
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
        with sftp_client.open(os.path.join(path, 'boundaries.raw')) as f:
            f_lines = f.readlines()          
        # close connection
        client.close()
    # stored locally?
    else:
        f_b = os.path.join(args.scf, 'boundaries.raw')
        with open(f_b, 'r') as f:
            f_lines = f.readlines()
    
    f_lines = [l.split() for l in f_lines]
    db = {l[0]:[datetime.datetime.strptime(l[1], '%Y-%m-%dT%H:%M:%S'), \
                datetime.datetime.strptime(l[2], '%Y-%m-%dT%H:%M:%S')] \
                for l in f_lines}
    
    print 'Choose from the following files:'
    for key in db.keys():
        if db[key][0] <= t_begin <= db[key][1] or \
            db[key][0] <= t_end <= db[key][1]:
            print '{:s}'.format(key)