#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
"""
Copyright (c) 2015 Joint Institute for VLBI -- ERIC
All rights reserved.

@author: Dmitry Duev <duev@jive.nl>, 2015

Created on Wed Feb 25 17:59:52 2015

"""

from pypride.vex import Vex

import argparse
import json
from collections import OrderedDict
import subprocess
import re
import os
import datetime

if __name__ == '__main__':
    '''
        If not imported, but run from the command line
    '''
    # create parser
    parser = argparse.ArgumentParser(prog='python ftptest.py', \
                formatter_class=argparse.RawDescriptionHelpFormatter,\
                description='Generate stuff for an ftp fringe test.')
    # optional arguments
    parser.add_argument('-a', '--ant', type=str, default='all',
                        help='antenna 2-letter codes')
    parser.add_argument('-r', '--refant', type=str, default=None,
                        help='Reference antenna 2-letter codes')
    parser.add_argument('-c', '--channels', type=str, default='all',
                        help='channels to use, counting from 1')
    parser.add_argument('-f', '--fft', type=int, default=1024,
                        help='number of fft points to use')
    parser.add_argument('-s', '--scans', type=str, default='all',
                        help='ftp-scans to generate stuff for')
    parser.add_argument('-d', '--dummy', action='store_true',
                        help='make a dummy list of file locations')
    parser.add_argument('-g', '--gabri', action='store_true',
                        help='create folders according to Gabriele\'s taste')

    # positional arguments
    parser.add_argument('vex', type=str, help='vex-file name')

    args = parser.parse_args()
    
    vex_file = args.vex
    
    vex = Vex(vex_file)
        
    exp_vex = vex['GLOBAL']['EXPER']
    exp = exp_vex.lower()
    start = vex['EXPER'][exp_vex]['exper_nominal_start']
    stop = vex['EXPER'][exp_vex]['exper_nominal_stop']
    
#    sta = [st for st in vex['STATION']]

    mode = vex['MODE']
    mod = [m for m in mode][0] # pick the first one
    freq = vex['MODE'][mod].getall('FREQ')[0][0] # pick the first one
    cha = vex['FREQ'][freq]
    chan = [ch[4] for ch in cha.getall('chan_def')]

    # get scans:
    scans = OrderedDict()
    for s in vex['SCHED']:
        disk2net = vex['SCHED'][s].getall('data_transfer')
        if len(disk2net)>0:
            start = vex['SCHED'][s]['start']
            start = datetime.datetime.strptime(start, '%Yy%jd%Hh%Mm%Ss')
            beg = datetime.timedelta(seconds=int(disk2net[0][3].split()[0]))
            end = datetime.timedelta(seconds=int(disk2net[0][4].split()[0]))
            scans[s] = OrderedDict()
            scans[s]['start'] = datetime.datetime.strftime(start + beg,\
                                                            '%Yy%jd%Hh%Mm%Ss')
            scans[s]['stop'] = datetime.datetime.strftime(start + end,\
                                                            '%Yy%jd%Hh%Mm%Ss')
            scans[s]['tint'] = int((end-beg).total_seconds())

            scans[s]['sta'] = [st[0] for st in vex['SCHED'][s].getall('station')]
            
            scans[s]['sta_ftp'] = []
#    print scans
    
    # mkdirs if they don't exist!
    if not os.path.isdir('html'):
        os.makedirs('html')
    if not os.path.isdir('output'):
        os.makedirs('output')
    if args.gabri:
        if not os.path.isdir('ctrl'):
            os.makedirs('ctrl')
    
    for s in scans.keys():

        no = int(s[2:])
        
        if args.scans!='all': # skip scan if unwanted
            if no not in map(int, args.scans.split(',')):
                continue
        
        # mkdirs if necessary:
        if not os.path.isdir(os.path.join('html','scan{:02d}'.format(no))):
            os.makedirs(os.path.join('html','scan{:02d}'.format(no)))
        if args.gabri:
            if not os.path.isdir(os.path.join('ctrl','scan{:02d}'.format(no))):
                os.makedirs(os.path.join('ctrl','scan{:02d}'.format(no)))
            # copy vex-file if it's not there already:
            vex_in_ctrl = os.path.join('ctrl','scan{:02d}'.format(no),\
                                                os.path.basename(vex_file))
            if not os.path.isfile(vex_in_ctrl):
                subprocess.Popen(['cp {:s} {:s}'.format(vex_file, vex_in_ctrl)],\
                                  stdout=subprocess.PIPE, shell=True)
        
        if not args.gabri:
            ctrl = '{:s}.scan{:02d}.ctrl'.format(exp, no)
        else:
            ctrl = os.path.join('ctrl','scan{:02d}'.format(no),\
                                '{:s}.scan{:02d}.ctrl'.format(exp, no))
        
        data = OrderedDict()   
        
    #    data['channels'] = [ "CH01", "CH02", "CH03", "CH04", "CH05", "CH06", "CH07",\
    #                        "CH08", "CH09", "CH10", "CH11", "CH12", "CH13", "CH14",\
    #                        "CH15", "CH16" ]
        if args.channels!='all':
            ch = map(int, args.channels.split(','))
            ch = [u'CH{:02d}'.format(c) for c in ch]
            data['channels'] = ch
        else:
            data['channels'] = chan
        
        data['cross_polarize'] = 'true'
        
        data['data_sources'] = OrderedDict()
        if args.dummy:
    #        data['data_sources'] = {st:['file:///scratch/ftp/ftpdata/'+\
    #                                    '{:s}_{:s}_no{:04d}.m5a'.format(exp,st,2)]\
    #                                for st in scans[s]['sta']}
            for st in scans[s]['sta']:
                if st in ('Bd','Zc','Sv'):
                    data['data_sources'][st] = ['file:///scratch/ftp/ftpdata/'+\
                        '{:s}_{:s}_no{:04d}.m5b'.format(exp,st.lower(),no)]
                else:
                    data['data_sources'][st] = ['file:///scratch/ftp/ftpdata/'+\
                        '{:s}_{:s}_no{:04d}.m5a'.format(exp,st.lower(),no)]
        else:
#            proc = subprocess.Popen(['ssh sfxc-h1 ls /scratch/ftp/ftpdata | grep \'{:s}.*{:04d}\''.\
#                        format(exp, no)], stdout=subprocess.PIPE, shell=True)
            proc = subprocess.Popen(['ls log | grep \'{:s}.*{:04d}\''.\
                        format(exp, no)], stdout=subprocess.PIPE, shell=True)
            (ls, err) = proc.communicate()
            
            ls = ls.split()
            
            if len(ls)>0:
                for l in ls:
                    st = re.search(r'_(.*)_', l) # station name
                    st = st.group(1).title()
                    if 'vdif' in l:
                        print 'Found a vdif-file for station {:s}, scan {:d}.'\
                                .format(st, no)+\
                              ' Pay attentioin and handle this manually!'
                        continue
                    else:
                        # add station to actual scan sta list
                        scans[s]['sta_ftp'].append(st)
                        data['data_sources'][st] = ['file:///scratch/ftp/ftpdata/'+\
                        '{:s}'.format(l)]
        
        data['exper_name'] = exp_vex
        data['integr_time'] = scans[s]['tint']
        data['message_level'] = 1
        data['number_channels'] = args.fft
        if not args.gabri:
            data['output_file'] = 'file://output/scan{:02d}.cor'.format(no)
            data['html_output'] = 'file://html/scan{:02d}/'.format(no)
        else:
            data['output_file'] = 'file://../../output/scan{:02d}.cor'.format(no)
            data['html_output'] = 'file://../../html/scan{:02d}/'.format(no)
        data['start'] = scans[s]['start']
        if args.dummy:
            data['stations'] = scans[s]['sta']
            if args.refant is None and len(ls)>0:
                data['reference_station'] = scans[s]['sta'][0]
            else:
                data['reference_station'] = args.refant
        else:
            if args.ant=='all':
                data['stations'] = scans[s]['sta_ftp']
            else:
                data['stations'] = args.ant.split(',')
            if args.refant is None and len(ls)>0:
                data['reference_station'] = scans[s]['sta_ftp'][0]
            else:
                data['reference_station'] = args.refant
        data['stop'] = scans[s]['stop']
        
        with open(ctrl, 'w') as outfile:
            json.dump(data, outfile, sort_keys=False, \
                      indent=4, separators=(',', ': '))