#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 16:29:02 2015

@author: Dmitry A. Duev

Calculate RA/Dec of a spacecraft given a vex-file with the schedule

"""

from pypride.classes import ephem
from pypride.vintlib import inp_set, load_sc_eph
from pypride.vex import Vex
import datetime
from astropy.time import Time
from astropy.coordinates import Angle#, SkyCoord, ICRS
from astropy import units as u
import numpy as np
import os
import argparse
    

if __name__ == '__main__':
    # create parser
    parser = argparse.ArgumentParser(prog='python scradec.py',
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description='Calculate RA/Dec of a spacecraft.')
    # optional arguments
    parser.add_argument('-e', '--epochs', type=str, default='middle',
                        help='epochs at which to ouput the RA/Dec\'s. ' +
                             'could be \'middle\' (of each source scan) ' +
                             'or \'all\'. defaults to \'middle\'')

    # positional argument
    parser.add_argument('vexfile', type=str, help='input vex-file')
    parser.add_argument('cfgfile', type=str, help='input config-file')
    parser.add_argument('source', type=str, help='source name')
    args = parser.parse_args()

    source = args.source.upper()
    
    vex = Vex(args.vexfile)

    scans = []
    times = []
    for scan in vex['SCHED']:
        if vex['SCHED'][scan]['source'] == source:
            start = datetime.datetime.strptime(vex['SCHED'][scan]['start'], '%Yy%jd%Hh%Mm%Ss')
            dur = int(vex['SCHED'][scan]['station'][2].split()[0])#/2
    #        print scan, start + datetime.timedelta(seconds=dur)
            scans.append(scan)
            times.append([start, dur])
    
    times = np.array(times)
    
    #%%
    eph = ephem(source)
    
    inp = inp_set(args.cfgfile)
    inp = inp.get_section('all')
    
    t_start = times[0,0]
    t_end   = times[-1,0] + datetime.timedelta(seconds=times[-1,1])
                       
    eph = load_sc_eph('S', source, t_start, t_end, inp,
                      uvw_calc=False, sc_rhophitheta=False, sc_xyz=False, load=True)
    
    epochs = args.epochs # 'middle' or 'all'
#    ra, dec, _ = eph.RaDec_bc_sec(2456655.5, 0.571571838630275593E+05, inp['jpl_eph'])
#    print 'ra={:22.18f} dec={:22.18f}'.format(ra, dec)
#    raise Exception('oops')
    if epochs == 'all':
        with open(os.path.join(inp['out_path'], '{:s}.{:s}.crd'.\
                        format(source.lower(), vex['GLOBAL']['EXPER'].lower())), 'w') as f:
            for scan, (start, dur) in zip(scans, times):
                for ii in range(dur+1):
                    
                    astrot = Time(str(start + datetime.timedelta(seconds=ii)),
                                    format='iso', scale='utc', precision=9)
                    
                    jd = astrot.tdb.jd1
                    ti = astrot.tdb.jd2
                    if ti<0:
                        jd -= 1
                        ti += 1
                    elif ti>=1:
                        jd += 1
                        ti -= 1
    #                print jd, ti
                    
                    ra, dec, _ = eph.RaDec_bc_sec(jd, ti*86400.0, inp['jpl_eph'])
                    
                    ra = Angle(ra, unit=u.rad)
                    dec = Angle(dec, unit=u.rad)
                    
                    radec = np.hstack((Angle(ra, unit=u.rad).hms,
                                        Angle(dec, unit=u.rad).dms))
                    print scan, 'ra = {:02.0f}h{:02.0f}m{:010.7f}s  dec = {:-3.0f}d{:02.0f}\'{:010.7f}\"'\
                            .format(*radec)
                    
        #            f.write('{:s}  {:s}  {:d} {:d} {:.6f}   {:d} {:d} {:.6f}\n'.format(scan,\
        #                    str(start + datetime.timedelta(seconds=ii)),\
        #                    int(ra.hms[0]), int(ra.hms[1]), ra.hms[2], \
        #                    int(dec.dms[0]), np.abs(int(dec.dms[1])), np.abs(dec.dms[2])))
        
                    f.write('{:s}   {:02.0f}h{:02.0f}m{:011.8f}s   {:-3.0f}d{:02.0f}\'{:09.6f}\"\n'.format(
                            datetime.datetime.strftime(start + \
                                datetime.timedelta(seconds=ii),'%Y-%m-%dT%H:%M:%S'),
                            ra.hms[0], ra.hms[1], ra.hms[2],
                            dec.dms[0], np.abs(dec.dms[1]), np.abs(dec.dms[2])))
    
    elif epochs == 'middle':
        with open(os.path.join(inp['out_path'], '{:s}.{:s}.crd'.\
                        format(source.lower(), vex['GLOBAL']['EXPER'].lower())), 'w') as f:
            for scan, (start, dur) in zip(scans, times):
                astrot = Time(str(start + datetime.timedelta(seconds=dur/2)),
                                format='iso', scale='utc', precision=9)
    
                # use barycentric approach - it agrees well with HORIZONS and old TASC.ESA.INT
                astrot = Time(str(start + datetime.timedelta(seconds=dur/2)),
                                format='iso', scale='utc', precision=9)
                jd = astrot.tdb.jd1
                ti = astrot.tdb.jd2
                if ti<0:
                    jd -= 1
                    ti += 1
                elif ti>=1:
                    jd += 1
                    ti -= 1
                ra, dec, _ = eph.RaDec_bc_sec(jd, ti*86400.0, inp['jpl_eph'])
    
                
                ra = Angle(ra, unit=u.rad)
                dec = Angle(dec, unit=u.rad)
                
                radec = np.hstack((Angle(ra, unit=u.rad).hms,
                                    Angle(dec, unit=u.rad).dms))
                print scan, astrot.datetime, \
                        'ra = {:02.0f}h{:02.0f}m{:010.7f}s  dec = {:-3.0f}d{:02.0f}\'{:010.7f}\"'\
                        .format(*radec)
    #            print scan, astrot.datetime, 'ra = {:16.10f}  dec = {:16.10f}'\
    #                    .format(Angle(ra, unit=u.rad).deg, Angle(dec, unit=u.rad).deg)
                    
                f.write('{:s}   {:02.0f}h{:02.0f}m{:011.8f}s   {:-3.0f}d{:02.0f}\'{:09.6f}\"\n'.format(
                        datetime.datetime.strftime(start +
                            datetime.timedelta(seconds=dur/2),'%Y-%m-%dT%H:%M:%S'),
                        ra.hms[0], ra.hms[1], ra.hms[2],
                        dec.dms[0], np.abs(dec.dms[1]), np.abs(dec.dms[2])))
