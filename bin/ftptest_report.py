#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""
Copyright (c) 2015 Joint Institute for VLBI in Europe (The Netherlands)
All rights reserved.

@author: Dmitry Duev <duev@jive.nl>, 2015

Created on Thu Mar  5 18:32:39 2015

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
    parser = argparse.ArgumentParser(prog='python ftptest_report.py', \
                formatter_class=argparse.RawDescriptionHelpFormatter,\
                description='Generate output stuff for an ftp fringe test.')
    # optional arguments
#    parser.add_argument('-a', '--ant', type=str, default='all',
#                        help='antenna 2-letter codes')
#    parser.add_argument('-r', '--refant', type=str, default=None,
#                        help='Reference antenna 2-letter codes')
    parser.add_argument('-d', '--dummy', action='store_true',
                        help='make a dummy list of file locations')

    # positional arguments
    parser.add_argument('vex', type=str, help='vex-file name')

    args = parser.parse_args()
    
    vex_file = args.vex
    
    vex = Vex(vex_file)
        
    exp_vex = vex['GLOBAL']['EXPER']
    exp = exp_vex.lower()
    comment = vex['EXPER'][exp_vex]['exper_description']
    start = vex['EXPER'][exp_vex]['exper_nominal_start']
    stop = vex['EXPER'][exp_vex]['exper_nominal_stop']
    
#    sta = [st for st in vex['STATION']]

    mode = vex['MODE']
    mod = [m for m in mode][0] # pick the first one
    freq = vex['MODE'][mod].getall('FREQ')[0][0] # pick the first one
    cha = vex['FREQ'][freq]
    chan = [ch[4] for ch in cha.getall('chan_def')]
    bw = float([ch[3] for ch in cha.getall('chan_def')][0].split()[0]) # MHz
    ncha = len(chan)
    samplesPerSec = float(vex['FREQ'][freq]['sample_rate'].split()[0]) # Msamples/s
    bitPerSample = samplesPerSec/bw

    # get ftp-scans:
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
            
            scans[s]['source'] = vex['SCHED'][s]['source']
#    print scans

    
    ftp_html = '{:s}.ftp.html'.format(exp)
    report_html = '{:s}.report.html'.format(exp)
    
    with open(ftp_html, 'w') as outfile:

        date = datetime.datetime.strptime(\
                        vex['EXPER'][exp_vex]['exper_nominal_start'], \
                        '%Yy%jd%Hh%Mm%Ss')        
        
        outfile.write('''<html>
<head>
<title>{:s}</title>
</head>
<body bgcolor='white' style='color:black'>
<h1 style='color:red'>{:s}</h1>

<p style='font-weight:bold'>{:s} - {:s} - 
selected results from the JIVE SFXC software correlator</p>
<p style='font-weight:bold'>{:s}</p>
'''.format(exp_vex, exp_vex, exp_vex, comment, \
            datetime.datetime.strftime(date, '%b %d, %Y')))
    
        # mode:
        outfile.write('''
<p style='font-weight:bold; color:green'>
<br>Mode: {:.0f} Mbps, {:d} BBC channels, {:.0f} Mhz filters, 2-bit sampling
</p><br>
    '''.format(ncha*samplesPerSec*bitPerSample, ncha, bw))
    
        # the scans
        for s in scans.keys():
    
            no = int(s[2:])

            if not args.dummy:
        #        proc = subprocess.Popen(['ssh sfxc-h1 ls /scratch/ftp/ftpdata | grep \'{:s}.*{:04d}\''.\
        #                    format(exp, no)], stdout=subprocess.PIPE, shell=True)
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
                    staz = scans[s]['sta_ftp']
                else:
                    staz = []
            else:
                staz = scans[s]['sta']
            
            outfile.write('''
<p style='font-weight:bold'>
<a href='http://www.evlbi.org/tog/ftp_fringes/{:s}/scan{:02d}/index.html'>scan{:02d}</a> 
{:s}: {:d} sec integration time, source = {:s}
<br>
</p>
'''.format(exp_vex, no, no, ', '.join(staz), scans[s]['tint'], scans[s]['source']))


        outfile.write('''
</body>
</html>
'''.format(exp_vex, exp_vex))

        
        