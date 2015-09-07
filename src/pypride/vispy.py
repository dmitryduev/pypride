#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  4 20:05:45 2015

@author: oasis
"""

import struct # pack/read binary data
from pypride.vintlib import *
from astropy.coordinates import SkyCoord
from time import time as _time
import multiprocessing
import argparse
from collections import Counter
#import os

import sys
##sys.path.append('/Users/oasis/_jive/python/vex')
#sys.path.append('../vex')

from pypride.vex import Vex

'''
#==============================================================================
# Class for parsing vex-files and running delay generation for SFXC
#==============================================================================
'''
class vispy(object):
    '''
        Class for parsing vex-files and running delay generation for SFXC
    '''
    
    def __init__(self, vex_file, inp_file, \
                    parallel=True, showTiming=True):
        '''
        Parse main info contained in the vex-file
        '''
        if showTiming: 
            tic = _time()
        
        # load input sittings:
        try:
            self.inp = inp_set(inp_file)
        except:
            raise Exception('Could not load input settings file.')
        
        # path to vex-file
        self.vex_file = vex_file
        
        # to show or not to show timing:
        self.showTiming = showTiming
        
        # run in parallel what could be run in parallel?
        self.parallel = parallel

        ''' read vex-file '''
        try:
            with open(vex_file, 'r') as f:
                self.f_lines = f.readlines()
        except:
            raise Exception('Failed to read the vex-file')

        
        ''' parse vex-file '''
        vex = Vex(vex_file)
        self.vex = vex
        
        ''' experiment name: '''
        self.exp_name = vex['GLOBAL']['EXPER']

        ''' get the exp date:'''
        start = vex['EXPER'][self.exp_name]['exper_nominal_start']
        self.date_start = datetime.datetime.strptime(start, '%Yy%jd%Hh%Mm%Ss')
        stop = vex['EXPER'][self.exp_name]['exper_nominal_stop']
        self.date_stop = datetime.datetime.strptime(stop, '%Yy%jd%Hh%Mm%Ss')
        # it is desired in a 'yymmdd' format:
        self.date_str = self.date_start.strftime('%y%m%d')
    
        ''' get mode freqs: ''' # first chan; convert from MHz to Hz
        self.mods = {m:1e6*float(vex['FREQ'][vex['MODE'][m].getall('FREQ')[0][0]].\
                     getall('chan_def')[0][1].split()[0])\
                     for m in vex['MODE']} # dictionary

        ''' get station names and short names: '''
        # list of stations with possible multiple names:
        sta_multi = {}
        with open(self.inp.sta_nam) as f:
            sta_multi_lines = f.readlines()
        sta_multi_lines = [l for l in sta_multi_lines if l[0]!='#']
        for f_line in sta_multi_lines:
                sta_multi[f_line[0:8].strip()] = f_line[8:].split()
        
        self.stations = {s:vex['STATION'][s]['ANTENNA'] for s in vex['STATION']}
        # check if the right (for vispy) name is used
        for sk, sv in self.stations.iteritems():
            for k, v in sta_multi.iteritems():
                if sv in v:
                    self.stations[sk] = k

        # station duplicates (might be e.g. the same station with a DBBC):
        dup = [item for item, count in \
                     Counter(self.stations.values()).iteritems() if count > 1]
        self.duplicate_stations = {d:[k for k,v in self.stations.iteritems() if v==d] \
                                        for d in dup}
        # remove duplicates from self.stations
        for _, vdup in self.duplicate_stations.iteritems():
            # keep the first short name only:
            for shn in vdup[1:]:
                del self.stations[shn]
        
        ''' get sources:'''
        self.sources = {}
        for s in vex['SOURCE']:
            c = SkyCoord(vex['SOURCE'][s]['ra'], \
                         vex['SOURCE'][s]['dec'], frame='icrs')
#            self.sources[s] = [c.ra.hms, c.dec.dms]
            self.sources[s] = [c.ra.rad, c.dec.rad]
            
        if self.showTiming: 
            toc = _time()
            print 'Initialising vispy took {:.1f} seconds.'.format(toc-tic)
 
       
    def makeObs(self, t_step=1, \
            delays=True, \
            dopplerPhaseCor=False, \
            dopplers=False, \
            staz='all'):
        ''' 
            Parse vex-file and create obs-objects for vint 
            staz - list of station short names to process or 'all'
        '''

        if self.showTiming: 
            tic = _time()
        
        # not 'all'? check that it's a list then:
        if staz!='all' and type(staz)!=type([]):
            raise Exception('staz should be \'all\' or a list')
        if staz!='all':
            staz = [st.title() for st in staz] # make 'em all St
        
        # initialize obs-objects:
        obsz = []
        
        # set input switches for delays
#        inp_swchs = self.inp.get_section('Switches')
        # get all sections!
        inp_swchs = self.inp.get_section('all')
        inp_swchs['delay_calc'] = True
        inp_swchs['uvw_calc'] = True
        inp_swchs['sc_rhophitheta'] = True
        # turn off forced eph recalculation - it's done one when run self.updates()
        inp_swchs['sc_eph_force_update'] = False
        
        # RadioAstron-Ground session? Get tracking station then:
        # [could be 'Pu' or 'Gt']:
        if 'Ra' in self.stations.keys():
            if 'Pu' in self.stations.keys():
                sta_ra_ts = self.stations['Pu']
            elif 'Gt' in self.stations.keys():
                sta_ra_ts = self.stations['Gt']
            else:
                print 'could not guess RadioAstron tracking station. set to Pu.'
                sta_ra_ts = 'PUSHCHIN'
        
        # obs-objects for delay calculation:
        for staSh, sta in self.stations.iteritems():
            # check whether this station is actually wanted:
            if staz=='all' or staSh in staz:
                for sou, sou_radec in self.sources.iteritems():
                    # source type:
                    if sou.lower()=='ra':
                        obs_type = 'R' # radioastron observations
                    elif sou.lower()[0:2]=='pr' or sou.lower()[0:2]=='pg':
                        obs_type = 'G' # gnss observations
                    elif sou.lower()=='vex' or sou.lower()=='mex' or \
                         sou.lower()=='her' or sou.lower()=='gaia' or \
                         sou.lower()=='ce3':
                        obs_type = 'S' # vex/mex/herschel/gaia/ce3 observations
                    else:
                        obs_type = 'C' # calibrator observations
                    # append
                    if staSh != 'Ra':
                        obsz.append( obs([self.inp.phase_center, sta], \
                                         sou, obs_type, \
                                         self.exp_name, sou_radec, inp=inp_swchs) )
                    else:
                        # if RA was observing, add Pu or Gt to the station list, 
                        # it'll be used for calculating formatter time offset for 
                        # the RA downlink stream
                        obsz.append( obs([self.inp.phase_center, sta, sta_ra_ts], \
                                         sou, obs_type, \
                                         self.exp_name, sou_radec, inp=inp_swchs) )
                                 
        # set input switches for Doppler correction
        if delays and dopplerPhaseCor:
#            inp_swchs = self.inp.get_section('Switches') # reset to default Falses
            inp_swchs = self.inp.get_section('all') # reset to default Falses
            inp_swchs['doppler_calc'] = True
            inp_swchs['sc_eph_force_update'] = False
            
            # obs-objects for GC Doppler correction calculation:
            for sou, sou_radec in self.sources.iteritems():
                # source type:
                if sou.lower()=='ra':
                    obs_type = 'R' # radioastron observations
                elif sou.lower()[0:2]=='pr' or sou.lower()[0:2]=='pg':
                    obs_type = 'G' # gnss observations
                elif sou.lower()=='vex' or sou.lower()=='mex' or \
                     sou.lower()=='her' or sou.lower()=='gaia' or \
                     sou.lower()=='ce3':
                    obs_type = 'S' # vex/mex/herschel/gaia/ce3 observations
                else:
                    obs_type = 'C' # calibrator observations
                    continue # no correction applies
                # append
                # phase center = 'GEOCENTR'
                obsz.append( obs([self.inp.phase_center], sou, obs_type, \
                                     self.exp_name, sou_radec, inp=inp_swchs) )
    
        ''' parse scans '''
        for s in self.vex['SCHED']:
            scan = self.vex['SCHED'][s].getall('station')
            t_scan_start = \
                datetime.datetime.strptime(self.vex['SCHED'][s]['start'], \
                                                      '%Yy%jd%Hh%Mm%Ss')
            mode = self.mods[self.vex['SCHED'][s]['mode']]
            # multiple sources per scan are possible (mult.ph.cen, off-beamGNSS)
            sources = self.vex['SCHED'][s].getall('source')
            for sta in scan:
                st = sta[0] # sta short name
                # offset from t_scan_start:
#                beg = int(sta[1].split()[0])
                # nominal scan start time is wanted instead
                beg = datetime.timedelta(seconds=0)
                t_start = t_scan_start + beg
                scanLength = int(sta[2].split()[0])
                end = datetime.timedelta(seconds=scanLength)
                t_stop = t_scan_start + end
                
                # add the scan to the proper obs object:
                # [delay calculation]
                N_sec = (end - beg).total_seconds()
#                scanLength = max(scanLength, N_sec)
                nobs = N_sec/t_step + 1 # Number of observations
                # check if this station is desired (I am Venused, I am fired):
                if st in self.stations.keys():
                    # if scan is shorter than ~1.5 minute, use t_step=1s
                    step = 1 if nobs<10 and t_step!=1 else t_step
                    for sou in sources:
                        [ob.addScan(t_start,step=step,stop=t_stop,freq=mode) \
                             for ob in obsz \
                                 if not ob.inp['doppler_calc'] and\
                                    ob.sta[1]==self.stations[st] and \
                                    ob.source==sou]
    
            # [doppler correction calculation]
            # should be done for one 'station' only (Geocentre)
            nobs = scanLength/t_step + 1 # use max N_sec
            # if 1-way
            step = 1 if nobs<10 and t_step!=1 else t_step
            for sou in sources:
                [ob.addScan(t_start, step=step, stop=t_stop) \
                     for ob in obsz \
                         if ob.source==sou and ob.sou_type!='C' and \
                            ob.inp['doppler_calc'] and \
                            (self.inp.dop_model=='bary1way' or \
                             self.inp.dop_model=='geo1way')]
            # 2(3)-way
            step = 1 if nobs<10 and t_step!=1 else t_step
            for sou in sources:
                [ob.addScan(t_start, step=step, stop=t_stop) \
                     for ob in obsz \
                         if ob.source==sou and ob.sou_type!='C' and \
                            ob.inp['doppler_calc'] and \
                            self.inp.dop_model=='bary3way']
                        
        ''' remove empty obs-objects '''
        obsz = [ob for ob in obsz if len(ob.tstamps)>0]
        
        self.obsz = obsz
        
        if self.showTiming: 
            toc = _time()
            print 'Creating obs-objects took {:.1f} seconds.'.format(toc-tic)
 
           
    def eop2vex(self, cat_eop):
        ''' 
            Add EOP section to a vex-file, if it's not there 
        '''
        try:
            # is it there already?
            self.vex['EOP']
            return
        except:
            print 'vex-file ' + self.vex_file + \
                  ' lacks EOP section needed for SFXC, adding...'
            eops = np.zeros((3,16))

            ''' get the exp date:'''
            date = self.date_start
            datem1 = self.date_start - datetime.timedelta(days=1)
            dddm1 = datem1.timetuple().tm_yday
            # exp date in mjd:
            mjd = mjuliandate(date.year,date.month,date.day)
    
            ''' get the relevant eop entries from the catalogue: '''
            with open(cat_eop, 'r') as fc:
                fc_lines = fc.readlines()
    
            for jj in range(len(fc_lines)):
                if fc_lines[jj][0]!=' ' and fc_lines[jj][0]!='*':
                    entry = [float(x) for x in fc_lines[jj].split()]
                    if len(entry) > 0 and entry[3] == mjd - 1:
                        eops[0] = entry
                        eops[1] = [float(x) for x in fc_lines[jj+1].split()]
                        eops[2] = [float(x) for x in fc_lines[jj+2].split()]
            
            with open(self.vex_file, 'r') as fv:
                fv_lines = fv.readlines()
            
            pos = fv_lines.index('$GLOBAL;\n')
            fv_lines.insert(pos+2,'     ref $EOP = EOP{:03d};\n'.format(dddm1))
            #finally, add eop section to the end of the file
            if fv_lines[-1]!='\n': fv_lines.append('\n')
            fv_lines.append('*----------------------------------------------------'+\
                            '---------------------------\n')
            fv_lines.append('$EOP;  * here delta_psi and delta_eps = dX and dY\n')
            fv_lines.append('*----------------------------------------------------'+\
                            '---------------------------\n')
            fv_lines.append('*\n')
    #        fv_lines.append('  def EOP'+str(ddd-1).zfill(3)+';\n')
            fv_lines.append('  def EOP{:03d};\n'.format(dddm1))
            fv_lines.append('    TAI-UTC = {:3.1f} sec;\n'.format(nsec(mjd)))
            fv_lines.append('    eop_ref_epoch = {:4d}y{:03d}d00h00m00s;\n'\
                            .format(datem1.year, dddm1))
            fv_lines.append('    num_eop_points = 3;\n')
            fv_lines.append('    eop_interval = 24 hr;\n')
            fv_lines.append('    ut1-utc = ' + '%-10.7f'%eops[0][6]+' sec : '+
                                               '%-10.7f'%eops[1][6]+' sec : '+
                                               '%-10.7f'%eops[2][6]+' sec;\n')
            fv_lines.append('    x_wobble = ' + '%-9.6f'%eops[0][4]+' asec : '+
                                                '%-9.6f'%eops[1][4]+' asec : '+
                                                '%-9.6f'%eops[2][4]+' asec;\n')
            fv_lines.append('    y_wobble = ' + '%-9.6f'%eops[0][5]+' asec : '+
                                                '%-9.6f'%eops[1][5]+' asec : '+
                                                '%-9.6f'%eops[2][5]+' asec;\n')
            fv_lines.append('    delta_psi = ' + '%-9.6f'%eops[0][8]+' asec : '+
                                                 '%-9.6f'%eops[1][8]+' asec : '+
                                                 '%-9.6f'%eops[2][8]+' asec;\n')
            fv_lines.append('    delta_eps = ' + '%-9.6f'%eops[0][9]+' asec : '+
                                                 '%-9.6f'%eops[1][9]+' asec : '+
                                                 '%-9.6f'%eops[2][9]+' asec;\n')
            fv_lines.append('  enddef;\n')
            fv_lines.append('*\n')
            fv_lines.append('*----------------------------------------------------'+\
                            '---------------------------\n')
            with open(self.vex_file,'w') as vexOutFile:
                for line in fv_lines:
                    vexOutFile.write(line)
                    
    def clocks(self):
        '''
            Extract linear clocks from vex-file
        '''
        try:
            self.vex['CLOCK']
            clock = {}
            for sta in self.vex['CLOCK']:
                epoch = datetime.datetime.strptime(\
                                    self.vex['CLOCK'][sta]['clock_early'][2],\
                                    '%Yy%jd%Hh%Mm%Ss')
                offset = float(self.vex['CLOCK'][sta]['clock_early'][1].split()[0])
                rate = float(self.vex['CLOCK'][sta]['clock_early'][3])
                clock[sta.title()] = [epoch, offset*1e-6, rate*1e-6]
            return clock
        except:
            return {}

    def updates(self):
        '''
            Do all kinds of stuff checks and update/download if necessary
        '''
        if self.showTiming: 
            tic = _time()
        
        # add the eop section to the vex-file if it is not there:
        self.eop2vex(self.inp.cat_eop)
        
        # sources that are stated in the vex "header" need not necessarily be in
        # the scan list:
        try:
            sources_obs = [ob.source for ob in self.obsz]
        except:
            raise Exception('Can\'t run update, run makeObs first')
        
        # get GNSS sp3 files (otherwise it could go mad if run in parallel..)
        for source in sources_obs:
            if source[0:2].lower()=='pr' or source[0:2].lower()=='pg':
                beg = datetime.datetime(self.date_start.year,
                                        self.date_start.month,
                                        self.date_start.day)
                end = datetime.datetime(self.date_stop.year,
                                        self.date_stop.month,
                                        self.date_stop.day)
                dd = (end-beg).days
                for d in range(dd+1):
                    load_sp3(self.inp.sc_eph_cat, source, \
                             beg+datetime.timedelta(days=d), load=False)
        
        # spacecraft ephs, if they were observed:
        sou = ['VEX', 'MEX', 'HER', 'RA', 'GAIA', 'CE3']
        soutyp = ['S', 'S', 'S', 'R', 'S', 'S']
        # gnss
        # make a list with GLONASS/GPS satellite names (with a margin, up to 40)
        for gg in range(1,41):
            sou.append('PR'+'{:02d}'.format(gg))
            soutyp.append('G')
            sou.append('PG'+'{:02d}'.format(gg))
            soutyp.append('G')
        
        inps = []
        for sousou, stst in zip(sou, soutyp):
            if sousou in sources_obs: # S/C was actually observed?
                t_begin, t_end = None, None
                for s in self.vex['SCHED']:
                    sou_scan = self.vex['SCHED'][s].getall('source')
                    if sousou in [ss.upper() for ss in sou_scan]:
                        t_start = datetime.datetime.strptime(\
                                            self.vex['SCHED'][s]['start'], \
                                                      '%Yy%jd%Hh%Mm%Ss')
                        if t_begin is None:
                            t_begin = t_start
                        N_sec = int(self.vex['SCHED'][s].\
                                    getall('station')[0][2].split()[0])
                        t_end = t_start + datetime.timedelta(seconds=N_sec)
                # make the ephem
                if t_begin is not None:
                    inp_swchs = self.inp.get_section('all')
                    inp_swchs['sc_eph_force_update'] = False
                    inps.append([stst, sousou, t_begin, t_end, \
                                  inp_swchs, False])
        
        n_sc = len(inps) # number of sc for which to make ephs
        if n_sc>0:
            if self.parallel: # Parallel way
                n_cpu = multiprocessing.cpu_count()
                # create pool
                pool = multiprocessing.Pool(np.min((n_cpu, n_sc)))
                # asyncronously apply sc_eph to each of inps
                pool.map_async(sc_eph, inps)
                # close bassejn
                pool.close() # we are not adding any more processes
                # tell it to wait until all threads are done before going on
                pool.join() 
            else: # Serial way
                for inp in inps:
                    sc_eph(inp)

        
        # check if RadioAstron observed, download/update its eph:
        sou = ['RA']
        
        for sousou in sou:
            if sousou in self.stations.values(): # RA observing scheduled?
                t_begin, t_end = None, None
                for s in self.vex['SCHED']:
                    sta = [st for st in self.vex['SCHED'][s].getall('station') \
                            if st[0].upper()==sousou]
                    if len(sta)>0:
                        sta = sta[0]
                        t_start = datetime.datetime.strptime(\
                                            self.vex['SCHED'][s]['start'], \
                                                      '%Yy%jd%Hh%Mm%Ss')
                        if t_begin is None:
                            t_begin = t_start
                        
                        N_sec = int(sta[2].split()[0])
                        t_end = t_start + datetime.timedelta(seconds=N_sec)
                # make the ephem
                if t_begin is not None:
                    inp_swchs = self.inp.get_section('all')
                    inp_swchs['sc_eph_force_update'] = False
                    load_sc_eph('C', sousou, t_begin, t_end, \
                                inp_swchs, load=False)
                

        ''' update/(down)load eops, meteo and iono data '''
        # check internet connection
        if internet_on():
            try:
                doup(self.inp.do_trp_calc, self.inp.do_ion_calc, \
                     self.inp.cat_eop, self.inp.meteo_cat, self.inp.ion_cat,\
                     self.date_start, self.date_stop, self.inp.iono_model)
            except Exception, err:
                print str(err)
                print 'catalogue updates failed'
        else:
            print 'no internet connection. can\'t update catalogues'
            
        if self.showTiming: 
            toc = _time()
            print 'Updating stuff took {:.1f} seconds.'.format(toc-tic)
    
    
    def calcDelays(self):
        '''
            Run vint for each of the created obs-objects
        '''
        # Delays, Uvws, Doppler, Etc. = DUDE
        if not self.parallel:
            ''' Serial way '''
            if self.showTiming: 
                tic = _time()
            for ob in self.obsz:
                # check that there's actually what to process..
                if len(ob.tstamps)>0:
                    # work 'in place':
                    ob.dude = vint_s(ob)
            if self.showTiming: 
                toc = _time()
                print 'Calculation took {:.1f} seconds.'.format(toc-tic)

        else:
            ''' Parallel way '''
            if self.showTiming:
                tic = _time()
            n_cpu = multiprocessing.cpu_count()
            # create pool
            pool = multiprocessing.Pool(n_cpu)
            # asyncronously apply vint_s to each of obsz
            result = pool.map_async(vint_s, self.obsz)
            # close bassejn
            pool.close() # we are not adding any more processes
            # tell it to wait until all threads are done before going on
            pool.join() 
            # get the ordered results
            dudez = result.get()
            for ob, res in zip(self.obsz, dudez):
                ob.dude = res
            if self.showTiming: 
                toc = _time()
                print 'Calculation took {:.1f} seconds.'.format(toc-tic)

    
    def packDelays(self, delay_type = 'group', \
              smoothing=True, ionPhaseCor=False, dopplerPhaseCor=False,\
              noIonDelay=False):     
        ''' 
            Pack dudes into binary SFXC del-files and human-readable txt-files
        '''
        # load input sittings:
        inp = self.inp
        
        # experiment name:
        exp_name = self.exp_name
        
        # output directory:
        out_path = inp.out_path
#        for ob in self.obsz:
#            print ob.dude.delay[:,0]
        # get station names:
        stations = []
        # 'uniqueify'
        [stations.append(ob.sta[1]) for ob in self.obsz \
                      if not ob.inp['doppler_calc'] and ob.sta[1] not in stations]

        # get station short names:
        stations_short = []
        for st in stations:
            stations_short.append([k for k,v in self.stations.iteritems() \
                                    if v==st][0])
        
        # append duplicates:
        for kdup, vdup in self.duplicate_stations.iteritems():
            for key in vdup[1:]:
                stations.append(kdup)
                stations_short.append(key)
        
        # get short name for the phase center
        phase_center_short = shname(inp.phase_center, inp.shnames_cat,\
                                    inp.shnames_cat_igs)[0]
        
        # uplink station info (for text ouput?)
    #    if dopplerPhaseCor:
    #        uplink_sta = [ob.sta[1] for ob in obsz if ob.inp['doppler_calc']][0]
    #        uplink_sta_short = shname(uplink_sta, inp.shnames_cat)
        
        # load clock offsets/rates at stations from vex-file:
        clock = self.clocks()
        
        # smooth dudes
        if smoothing:
            if self.showTiming: 
                tic = _time()
            if self.parallel:
                n_cpu = multiprocessing.cpu_count()
                # create pool
                pool = multiprocessing.Pool(n_cpu)
                # asyncronously apply smoothie to each of obsz
                result = pool.map_async(smoothie, self.obsz)
                # close bassejn
                pool.close() # we are not adding any more processes
                pool.join() # wait until all threads are done before going on
                # get the ordered results back into obsz
                self.obsz = result.get()
            else:
                for ob in self.obsz:
                    ob.smoothDude(tstep=1)
            if self.showTiming: 
                toc = _time()
                print 'Smoothing took {:.1f} seconds.'.format(toc-tic)
        
        if self.showTiming: 
            tic = _time()
            
        # integrate freq predictions into phase correction on a scan basis
        # 2nd condition: durak-protection -- check that a s/c was actually observed
        sou_scans = {} # init here for parallelisation
        if dopplerPhaseCor and len([ob for ob in self.obsz if ob.sou_type!='C'])>0:
            # dictionary with [tstamps, gc doppler freq prediction] for each s/c
            dpc = {ob.source: np.array(zip(ob.tstamps, np.sum(ob.dude.doppler, 1))) \
                          for ob in self.obsz \
                           if len(ob.tstamps)>0 and ob.inp['doppler_calc']}
            # add int sec t stamps for interpolation
            for k in dpc.keys():
                t_0 = datetime.datetime(dpc[k][0,0].year, \
                                        dpc[k][0,0].month, dpc[k][0,0].day)
                t = np.array([t.hour*3600.0 + t.minute*60.0 + t.second + \
                                (t-t_0).days*86400.0 for t in dpc[k][:,0]])
                dpc[k] = np.column_stack(( t, dpc[k] ))
            # result: [int sec tstamps, datetime tstamps, gc doppler freq prediction]
    
    #        print dpc
            
            # list of s/c - dpc.keys()
            
            # scan start/stop times + source names
            sstz_dpc = np.array([ [[scanStartTime, scanStopTime, ob.source] \
                    for scanStartTime, scanStopTime in \
                            zip(ob.scanStartTimes,ob.scanStopTimes)] \
                    for ob in self.obsz if ob.inp['doppler_calc'] and len(ob.tstamps)>0 ])
            # flatten list, because it can be irregular
            sstz_dpc = flatten(sstz_dpc)
            
            if len(sstz_dpc)>0: # there actually exist non-empty scans
                # reshape
                sstz_dpc = np.reshape(sstz_dpc, (-1, 3))
                # sort according to 0th column, which is time stamp
                sstz_dpc = sstz_dpc[np.argsort(sstz_dpc[:,0])]
    
            
            # make a dictionary out of sstz_dpc
            # {sou_name: [ [t_start, t_stop, [poly]], 
            #               ...
            #              [t_start, t_stop, [poly]] ]}
#            sou_scans = {}
            for sc in dpc.keys():
                sou_scans[sc] = []
            for sstz in sstz_dpc:
                sc = sstz[2]
                # allocate space for phase polies and mu_t for rescaling
                sou_scans[sc].append([sstz[0], sstz[1], [], []])
    #        print sou_scans
            
            # make optimal poly interpolants for freq predictions
            # iterate over s/c
            for sc in sou_scans.keys():
                # iterate over each scan:
                for ii, scan in enumerate(sou_scans[sc]):
                    # create mask using using current scan t_start, t_stop
                    maska = np.ma.masked_outside(dpc[sc][:,1], scan[0], scan[1]).mask
                    # mask t [sec] -> masked np.ma.array
                    t = np.ma.array(dpc[sc][:,0], mask=maska, dtype=float)
                    # -> np.array
                    t = np.ma.compressed(t)
                    # scale t for a more robust fit:
                    mu_t = np.array([np.mean(t), np.std(t)])
                    t = (t-mu_t[0])/mu_t[1]
                    # mask f [Hz] -> masked np.ma.array
                    f = np.ma.array(dpc[sc][:,2], mask=maska, dtype=float)
                    # -> np.array
                    f = np.ma.compressed(f)
                    # subtract f at t[0] (to decrease dynamic range):
                    f -= f[0]
                    # make an optimal polyfit:
                    ffit = optimalFit(t, f, min_order=3, max_order=7, fit_type='poly')
                    # integrate it to a phase poly
                    pint = np.polyint(ffit.best_estimator_.coef_)
                    # save fit to phase by integrating p.best_estimator_.coef_
                    sou_scans[sc][ii][2] = pint
                    sou_scans[sc][ii][3] = mu_t

#            print sou_scans

        # obs-obj w Doppler predictions for phase correction not needed further:
    #    obsz_dop = [ob for ob in self.obsz if ob.inp['doppler_calc']]
        self.obsz = [ob for ob in self.obsz if not ob.inp['doppler_calc']]    
            
        if self.parallel:
            n_cpu = multiprocessing.cpu_count()
            # create pool
            pool = multiprocessing.Pool(n_cpu)
            # asyncronously apply packie to each of inps
            inps = [(st, st_sh, self.obsz, phase_center_short, out_path, \
                     exp_name, sou_scans, ionPhaseCor, dopplerPhaseCor, \
                     noIonDelay, delay_type, clock, self.vex) \
                     for (st, st_sh) in zip(stations, stations_short)]
            pool.map_async(packie, inps)
            # close bassejn
            pool.close() # we are not adding any more processes
            pool.join() # wait until all threads are done before going on

        else: # serial way
            # for each station stick all the data together, then sort it by datetime
            inps = [(st, st_sh, self.obsz, phase_center_short, out_path, \
                     exp_name, sou_scans, ionPhaseCor, dopplerPhaseCor, \
                     noIonDelay, delay_type, clock, self.vex) \
                     for (st, st_sh) in zip(stations, stations_short)]
            for inpu in inps:
                packie(inpu)
                

        if self.showTiming: 
            toc = _time()
            print 'Output took {:.1f} seconds.'.format(toc-tic)

def sc_eph(inps):
    # helper func to make s/c ephs (in parallel)
    stst, sousou, t_begin, t_end, inp, load = inps
    load_sc_eph(stst, sousou, t_begin, t_end, inp, load=load)

def packie(inps):
    # unpack inps:
    sta, station_short, obsz, phase_center_short, out_path, exp_name,\
        sou_scans, ionPhaseCor, dopplerPhaseCor, noIonDelay,\
        delay_type, clock, vex = inps

    # check that there are data to pack:
    if len([np.array(ob.tstamps) for ob in obsz \
            if ob.sta[1]==sta and len(ob.tstamps)>0])==0:
        return # skip 'empty' station

    # L_G for Ra observations
    if station_short.title() == 'Ra':
        const = constants('421')
        L_G = const.L_G
    
    # dump everything to files:
    out_txt = 'vdm.'+exp_name+'.'+\
              phase_center_short + station_short.lower() + '.txt'
    out_del = exp_name+'_'+station_short.title()+'.del'
    print 'output to ' + out_del + ' and ' + out_txt
    
    # make headers for txt output:
    header = '# A priori delay table for baseline '\
             + phase_center_short.title() + '-' \
             + station_short.title() + ', experiment ' \
             + exp_name + '.\n#\n' \
             + '# Produced by vispy software package on ' \
             + str(datetime.datetime.utcnow()) + ' UTC.\n#\n' \
             + '# Column descriptions:\n' \
             + '#   1-6:   UTC time stamp.\n' \
             + '#     7:   Source name.\n' \
             + '#  8-10:   Far-field sources: [u v w],\n'\
             + '#          computed numerically as c*[d(tau)/d(RA)/cos(D),\n' \
             + '#                                     d(tau)/d(D),\n' \
             + '#                                     tau],\n' \
             + '#          c - speed of light in vacuum.\n' \
             + '#          Near-field sources: [c*d(tau)/d(theta)/cos(phi),\n' \
             + '#                               c*d(tau)/d(phi),\n' \
             + '#                               d(tau)/d(rho)].\n' \
             + '#    11:   Geometric delay.\n' \
             + '#    12:   Delay due to thermal expantion of telescope(s).\n' \
             + '#    13:   Delay due to axis offset.\n' \
             + '#    14:   Tropospheric delay (model used: VMF1).\n' \
             + '#    15:   Ionospheric delay (for base freq of first BBC).\n' \
             + '# 16-17:   Phase [rad] and amplitude corrections applied.\n'
    if len(clock)>0:
        header += '#    18:   LO clock offset.\n'
    if station_short.title()=='Ra':
        if len(clock)==0:
            header +=  '#    18:   d(proper time)/d(TT) at Ra position.\n' \
                 + '#    19:   LT from Ra to receiving station \n' \
                 + '#          corrected for instrumental and propagation effects.\n'\
                 + '#    20:   Total delay, which goes into correlator.\n'
        else:
            header +=  '#    19:   d(proper time)/d(TT) at Ra position.\n' \
                 + '#    20:   LT from Ra to receiving station \n' \
                 + '#          corrected for instrumental and propagation effects.\n'\
                 + '#    21:   Total delay, which goes into correlator.\n'\
                 + '#    Note: for Ra, columns 12-15 refer to receiving station.'
    header += '#\n' \
           + '# For details about software and models used therein, see\n' \
           + '# Duev et al. A&A 541, A43 (2012), ' \
           + 'http://dx.doi.org/10.1051/0004-6361/201218885\n#\n'
        
    # build skeleton and use it
    with open(os.path.join(out_path, exp_name, out_txt),'w') as f_txt, \
         open(os.path.join(out_path, exp_name, out_del),'wb') as f_del:

        # txt header:
        f_txt.write(header)

        # binary header (header_size, sta_name):
        line = struct.pack('<i2sx', 3, station_short.title())
        f_del.write(line)

        sched = sorted([s for s in vex['SCHED']])
        for s in sched:
            scan = vex['SCHED'][s].getall('station')
            sta_scan = [sn[0] for sn in scan] # station names for this scan
            # station scheduled for this scan?
            if station_short.title() in sta_scan:
                t_start = datetime.datetime.strptime(vex['SCHED'][s]['start'], \
                                                          '%Yy%jd%Hh%Mm%Ss')
                scanLength = int(scan[0][2].split()[0])
                end = datetime.timedelta(seconds=scanLength)
                t_stop = t_start + end
                # multiple sources per scan are possible (mult.ph.cen, off-beamGNSS)
                sources = vex['SCHED'][s].getall('source')
                for sou in sources:
                    # get this source for this station
                    tz = np.array([ob.tstamps for ob in obsz \
                          if ob.sta[1]==sta and len(ob.tstamps)>0 \
                          and ob.source==sou])
                    tz = np.squeeze(tz)
                    
                    if len(tz)==0:
                        continue

                    # Delays
                    dz = np.array([ob.dude.delay for ob in obsz \
                          if ob.sta[1]==sta and len(ob.tstamps)>0 \
                          and ob.source==sou])
                    dz = np.squeeze(dz)

                    # uvw-'Projections'
                    pz = np.array([ob.dude.uvw for ob in obsz \
                          if ob.sta[1]==sta and len(ob.tstamps)>0 \
                          and ob.source==sou])
                    pz = np.squeeze(pz)
                    
                    # freqs for iono phase correction:
                    fz = np.array([ob.freqs for ob in obsz \
                          if ob.sta[1]==sta and len(ob.tstamps)>0 \
                          and ob.source==sou])
                                        
                    # mask the current scan:
                    maska = np.ma.masked_inside(tz, t_start, t_stop).mask
                    tz = tz[maska]
                    dz = dz[maska,:]
                    if noIonDelay: # don' wan' ion'?
                        dz[:,4] = np.zeros_like(dz[:,4])
                    pz = pz[maska,:]
                    fz = np.array([f for f in fz[0] \
                                    if f[0]>=t_start and f[1]<=t_stop])
                    
                    ## go-go-go
                    # this is for trapz integration step:
                    t_step = (tz[1]-tz[0]).total_seconds()
                    # this is to account for overnighters within a scan
                    # (if tstamp goes over 86400, but only within a scan)
                    t_start = datetime.datetime(tz[0].year, tz[0].month, tz[0].day)
                    # this is to calculate phase correction at scan start
                    # (if first scan is on s/c) and not load poly for each t_stamp
                    if dopplerPhaseCor and sou in sou_scans.keys():
                        t_start_sec = tz[0].hour*3600.0 + \
                                      tz[0].minute*60.0 + tz[0].second
                        scans = sou_scans[sou]
                        scan = [x for x in scans if x[0] <= tz[0] <= x[1]][0]
                        mu_t = scan[3]
                        pint = scan[2]
                        ph_0 = np.polyval(pint, (t_start_sec - mu_t[0]) / mu_t[1])

                    ## binary output:
                    # scan start - write sou_name and mjd then
                    mjd = int(mjuliandate(t_start.year,t_start.month,t_start.day))
                    line = struct.pack('<80sxi', sou.ljust(80), mjd)
                    f_del.write(line)
                    
                    # if not RA:
                    if station_short.title()=='Ra':
                        # set at scan start
                        LT_downlink_sst = dz[0,6] + sum(dz[0,1:5])
                    
                    for jj, _t in enumerate(tz):
                        # take possible overnighters into account
                        tstamp = _t.hour*3600.0 + _t.minute*60.0 + _t.second + \
                                     (_t-t_start).days*86400.0
                        # if not RA:
                        if station_short.title()!='Ra':
                            if delay_type=='group': # group delay +iono
                                delay = sum(dz[jj,0:5])
                            elif delay_type=='phase': # phase delay -iono
                                delay = sum(dz[jj,0:4]) - dz[jj,4]
                        # if RA, integrate dtau/dTT and add LTs to receiving 
                        # station at scan start times
                        else:
                            # geometric part
                            delay = dz[jj,0]
                            # d(proper time @ Ra)/d(TT) [which = d/d(UTC)]:
                            # intgrate dtau/dTT from beginnig of scan to current t.
                            # L_G is to convert TCG to TT, as dtau_dtt was actually
                            # calcuted in coordinate time scale
                            proper_vs_utc = L_G*(_t-tz[0]).total_seconds() + \
                                            np.trapz(dz[0:jj+1,5], dx=t_step)
                            # add to full delay
                            delay += LT_downlink_sst + proper_vs_utc
                        
                        # phase correction
                        phs = 0.0
                        amp = 1.0
                        # ionosphere:
                        if ionPhaseCor and delay_type=='group':
                            # get proper freq value:
                            _f = [x[2] for x in fz if x[0]<=_t<=x[1]][0]
                            # -2 * 2pi * tau_iono * f (subtract what's already there 
                            # and change sign)
                            phs += -2.0 * 2.0*pi*dz[jj,4]*_f
                            # FIXME: WTF?? shouldn't it be simply?:
                            # (coz nothing's there!!)
#                            phs -= 2.0*pi*_f * dz[jj,4]
                            # no, iono delay does 'spoil' phases, thus *2, right??
                        # Doppler smearing:
                        if dopplerPhaseCor and sou in sou_scans.keys():
                            # find correct phase poly for correct source])
                            # also convert it to radians
                            to = (tstamp - mu_t[0])/mu_t[1]
                            phs += 2.0*np.pi*(np.polyval(pint, to) - ph_0)*mu_t[1]
                        
                        # write t u v w delay
                        uvw = pz[jj]
                        line = struct.pack('<5d', *flatten([tstamp, list(uvw), delay]))
                        f_del.write(line)
                        # write phase/amp correction
                        line = struct.pack('<2d', phs, amp)
                        f_del.write(line)
                        
                        ## stack line together for txt output:
                        line = '{:4d} {:02d} {:02d} {:02d} {:02d} {:06.3f} '.\
                                format(_t.year,_t.month,_t.day,\
                                       _t.hour,_t.minute,_t.second) + \
                               '{:>20s}   '.\
                                format(sou) + \
                                '{:23.16e} {:23.16e} {:23.16e}   '.\
                                format(*pz[jj,:]) + \
                                '{:23.16e} {:14.7e} {:14.7e} {:14.7e} {:14.7e}   '.\
                                format(*dz[jj,:5]) + \
                                '{:23.16e} {:23.16e}'.\
                                format(phs, amp)
                        # delay offset+rate at epoch
                        # check that vex-file had $CLOCK section and there's
                        # entry for the current station
                        if len(clock)>0 and (station_short in clock.keys()):
                            t_lo = clock[station_short][0]
                            offset = clock[station_short][1]
                            rate = clock[station_short][2]
                            dt = (_t - t_lo).total_seconds()
                            line += ' {:23.16e}'.format(offset + dt*rate)
                        # additional info for RadioAstron
                        if station_short.title()=='Ra':
                            # add dtau/dTT and LT to receiving station for Ra
                            line += '   {:23.16e} {:23.16e}'.format(*dz[jj,5:])
                            # total delay, which goes into correlator
                            line += '   {:23.16e}'.format(delay)
                        # new line
                        line += '\n'
                        # write it
                        f_txt.write(line)
                        
                    # scan stop - write trailing zeros
                    line = struct.pack('<7d', *list(np.zeros(7)))
                    f_del.write(line)
    
# helper function to perform parallel smoothing
def smoothie(ob):
    # .. to apply class method on list of objects in parallel
    ob.smoothDude(tstep=1)
    return ob
    
#%%
def main():
    # create parser
    parser = argparse.ArgumentParser(prog='vispy.py', \
                formatter_class=argparse.RawDescriptionHelpFormatter,\
                description='Computation of VLBI delays.')
    # optional arguments
    parser.add_argument('-p', '--parallel', action='store_true',
                        help='run computation in parallel mode')
    parser.add_argument('-t', '--showTiming', action='store_true',
                        help='display timing info')
    parser.add_argument('-s', '--stations', type=str, default='all',
                        help='2-letter station code(s) to process. '+\
                             'defaults to \'all\'')
    parser.add_argument('-ts', '--timeStep', metavar='SECONDS',\
                        type=int, default=10,
                        help='time step [sec] for delay calculation. '+\
                             'defaults to 10 sec')
    parser.add_argument('-nid', '--noIonDelay', action='store_true',
                        help='do not include ionospheric delay')                         
    parser.add_argument('-ipc', '--ionPhaseCor', action='store_true',
                        help='compute ionospheric phase correction')
    parser.add_argument('-dpc', '--dopplerPhaseCor', action='store_true',
                        help='compute Doppler phase correction for spacecraft')
    parser.add_argument('-i', '--info', action='store_true',
                        help='display info about the vex-file')
    # positional argument
    parser.add_argument('vexfile', type=str, help='input vex-file')
    parser.add_argument('cfgfile', type=str, help='input config-file')
    args = parser.parse_args()
    
    
    if args.info:
        ''' Display info and exit '''
        # parse vex-file
        info = vispy(args.vexfile, args.cfgfile, parallel=False, showTiming=True)
        print 'Summary of the vex-file {:s}'.format(args.vexfile)
        print 'Station list with codes and names:'
        print info.stations
        print 'Source list:'
        print info.sources
        print 'Clock information:'
        print info.clocks()
        print 'Time in seconds spent on a source on a particular baseline:'
        info.makeObs(t_step=1, dopplerPhaseCor=False, staz='all')
        for ob in info.obsz:
            print ob
        sys.exit(1)
    
    else:
        ''' Run computation '''
        # parse input
        parallel = args.parallel
        showTiming = args.showTiming
        ionPhaseCor = args.ionPhaseCor
        dopplerPhaseCor = args.dopplerPhaseCor
        noIonDelay = args.noIonDelay
        
        # only some specific stations wanted?
        staz = 'all'
        if args.stations=='all':
            staz = args.stations
        else:
            staz = args.stations.split(',')
            # capitalise:
            staz = [st.title() for st in staz]

        # non-standard t_step?
        t_step = args.timeStep        
        
        # run computation:
        v = vispy(args.vexfile, args.cfgfile, parallel=parallel, showTiming=showTiming)
        v.makeObs(t_step=t_step, dopplerPhaseCor=dopplerPhaseCor, staz=staz)
        v.updates()
        v.calcDelays()
        v.packDelays(delay_type = 'group', smoothing=True, \
                     ionPhaseCor=ionPhaseCor, dopplerPhaseCor=dopplerPhaseCor,\
                     noIonDelay=noIonDelay)
        # a happy end
        sys.exit(1)
   
   
#%%    
    
if __name__ == '__main__':
    '''
        If not imported, but run from the command line
    '''
    main()