# -*- coding: utf-8 -*-
"""
Created on Sun Feb  1 13:56:59 2015

@author: oasis
"""

import sys
sys.path.append('/Users/oasis/_jive/python/vex')
#sys.path.append('../vex')

import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('whitegrid') # plot em niice!

from pypride.vex import Vex

from pypride.vintlib import *

from astropy.coordinates import SkyCoord

class gnss(object):
    
    def __init__(self, date_t_start, date_t_end, gps=True, glonass=True):
        self.date_t_start = date_t_start
        self.date_t_end = date_t_end
        
        days = [datetime.datetime(date_t_start.year,date_t_start.month,\
                                date_t_start.day)]
        dd = (days[0]-datetime.datetime(date_t_end.year,date_t_end.month,\
                                    date_t_end.day)).days
    #    if dd>0:
        for d in range(1,dd+1):
            days.append(days[0]+datetime.timedelta(days=d))

        self.days = days
        
        self.sources = {}

        if not gps and not glonass:
            raise Exception('Choose at least one GNSS')
        
        self.gps = gps
        self.glonass = glonass
        
        
        for day in self.days:
            self.sp3(day)
    
    def sp3(self, day):
        inp = inp_set('inp.cfg')
    
        sc_eph_cat = inp.sc_eph_cat
        cat_eop = inp.cat_eop
#        jpl_eph = inp.jpl_eph
        
#        const = constants()

        # set up dates and names:
        mjd = mjuliandate(day.year, day.month, day.day)
        jd = mjd + 2400000.5
        dow = np.fmod(jd + 0.5, 7) + 1
        gnss_week = np.floor( (mjd - mjuliandate(1980,1,7))/7.0 )
        if dow==7:
            dow=0
            gnss_week += 1
        
        gnss_week = int(gnss_week)
        dow = int(dow)
        
        gnss_sp3s = []
        sp3_names = []
        # GLONASS
        if self.glonass:
            gnss_sp3s.append(os.path.join(sc_eph_cat,'raw_gnss/igl{:04d}{:01d}.sp3'.\
                                        format(gnss_week, dow)))
            sp3_names.append('igl{:04d}{:01d}.sp3'.format(gnss_week, dow))
        # GPS
        if self.gps:
            gnss_sp3s.append(os.path.join(sc_eph_cat,'raw_gnss/igs{:04d}{:01d}.sp3'.\
                                        format(gnss_week, dow)))
            sp3_names.append('igs{:04d}{:01d}.sp3'.format(gnss_week, dow))

        for gnss_sp3, sp3_name in zip(gnss_sp3s, sp3_names):
            # doesn't exist? download first then:
            if not os.path.isfile(gnss_sp3):
                try:
                    ftp = FTP('cddis.nasa.gov')
                    ftp.login() # user anonymous, passwd anonymous
                    if 'igl' in gnss_sp3:
                        folder = 'glonass'
                    elif 'igs' in gnss_sp3:
                        folder = 'gps'
                    ftp.cwd(folder+'/products/{:04d}'.format(gnss_week))
                    if '{:s}.Z'.format(sp3_name) in ftp.nlst(): # check that the file is there
                        ftp.retrbinary('RETR {:s}.Z'.format(sp3_name), \
                                   open('{:s}.Z'.format(gnss_sp3), 'wb').write)
                        # uncompress:
                        print 'uncompressing: ' + '{:s}.Z'.format(sp3_name) + '...'
                        os.system('uncompress -f {:s}'.format(gnss_sp3))
                    else:
                        raise Exception('file {:s} not found on the server. fail!'.\
                                format(sp3_name+'.Z'))
                    ftp.quit()
                    
                except Exception, err:
                    print str(err)
                    print 'Failed to download {:s} from cddis.nasa.gov'.format(sp3_name)
            
#            try:
            # load:
            with open(gnss_sp3) as f:
                f_lines = f.readlines()
            
            sources = {}
            
            for line in f_lines:
#                ln = re.findall(r'[\w.\s]+:[\s]+[\d\se.,-]+', line)
                if len(line)>0 and line[0:2]=='+ ':
                    ln = re.findall(r'[G,R]\d+', line)
                    for s in ln:
                        try:
                            # are we loading the next day?
                            sources['P'+s]
                        except:
                            sources['P'+s] = ephem('P'+s)
#            print sources

            # remove comments:
            f_lines = [l for l in f_lines if l[0] not in ['#','+','%','/']]
            # extract time in GPS time scale and convert to UTC:
            utc_gps = nsec(mjd) - nsec(mjuliandate(1980,1,6))
            time = [map(float, t[1:].split()) for t in f_lines if t[0]=='*']
            time_dt = np.array([datetime.datetime(*map(int, t)) - \
                    datetime.timedelta(seconds=utc_gps) for t in time])
            time = [[t.year, t.month, t.day, t.hour, t.minute, t.second] \
                    for t in time_dt]
            
            time = np.array(time)
            
            # done with loading gtrs eph. now, do gcrs
            N_obs = len(time)
            ob = obs(['DUMMY'],'DUMMY','C')
            t_step = (time_dt[1]-time_dt[0]).total_seconds()
            ob.addScan(time_dt[0], t_step, nobs=N_obs)
            
            #load eops
            ''' get the relevant eop entries from the catalogue: '''
            mjd_start = mjuliandate(day.year, day.month, day.day)
            with open(cat_eop, 'r') as fc:
                fc_lines = fc.readlines()
            eops = np.zeros((7,7)) # +/- 3 days
            for jj in range(len(fc_lines)):
                if fc_lines[jj][0]!=' ' and fc_lines[jj][0]!='*':
                    entry = [float(x) for x in fc_lines[jj].split()]
                    if len(entry) > 0 and entry[3] == np.floor(mjd_start) - 3:
                        for kk in range(7):
                            eops[kk,0] = entry[3] # mjd
                            eops[kk,1] = entry[6] # UT1-UTC
                            eops[kk,2] = entry[6] - nsec(entry[3]) # UT1 - TAI
                            eops[kk,3] = entry[4] # x_p
                            eops[kk,4] = entry[5] # y_p
                            eops[kk,5] = entry[8] # dX
                            eops[kk,6] = entry[9] # dY
                            entry = [float(x) for x in fc_lines[jj+kk+1].split()]
                        break #exit loop
    
            ''' go-go-go: '''
            for s in sources.keys():
                #allocate memory:
                sc_gcrs = np.zeros((N_obs,9)) # allocate GCRS pos
                sc_gtrs = np.zeros((N_obs,9))
                sc_gtrs[:, 0:6] = time

                xyz = np.array([map(float, x.split()[1:4]) for x in f_lines \
                        if s.lower() in x.lower()])
                sc_gtrs[:, 6:] = xyz

                if len(sources[s].gtrs)==0:
                    UT = np.array([(t - day).total_seconds()/86400. \
                                    for t in time_dt])
                    sources[s].UT = UT
                    sources[s].gtrs = sc_gtrs
                else: # append to the existing stuff
                    t0 = datatime.datetime(map(int, sources[s].gtrs[0,0:6]))
                    UT = np.array([(t - t0).total_seconds()/86400. \
                                    for t in time_dt])
                    sources[s].UT = np.append(sources[s].UT, UT)
                    sources[s].gtrs = np.append(sources[s].gtrs, \
                                                     sc_gtrs, axis=0)

#                radec = []
                for ii, tstamp in enumerate(ob.tstamps):
                    ''' set dates: '''
                    mjd = mjuliandate(tstamp.year, tstamp.month, tstamp.day)
                    UTC = (tstamp.hour*3600.0 + tstamp.minute*60.0 + \
                            tstamp.second)/86400.0
        
                    ''' compute tai & tt '''
                    TAI, TT = taitime(mjd, UTC)
                    
                    ''' interpolate eops to t_obs '''
                    UT1, eop_int = eop_iers(mjd, UTC, eops)
        
                    ''' rotation matrix IERS '''
                    r2000 = ter2cel(tstamp, eop_int, mode='der')
        
                    ''' GCRS position of GNSS '''
                    sc_gcrs[ii,0:6]  = sc_gtrs[ii,0:6]
                    sc_gcrs[ii,6:9]  = dot(r2000[:,:,0], sc_gtrs[ii,6:9].T)
                    
                    # this is rough
#                    r = sc_gcrs[ii,6:9]
#                    ra = np.arctan2(r[1],r[0]) # right ascention
#                    dec = np.arctan(r[2]/np.sqrt(r[0]**2+r[1]**2)) # declination
#                    radec.append([ra, dec])

                if len(sources[s].gcrs)==0:
#                    sources[s].radec = np.array(radec)
                    sources[s].gcrs = sc_gcrs
                else: # append to the existing stuff
#                    sources[s].radec = np.append(sources[s].radec,\
#                                                      radec, axis=0)
                    sources[s].gcrs = np.append(sources[s].gcrs,\
                                                     sc_gcrs, axis=0)
    #                print s, self.sources[s].radec
            for s in sources.keys():
                self.sources[s] = sources[s]

#            except Exception, err:
#                pass
            
    def minDistance(self, ra, dec, time):
        '''
            Get angular distances from a given point to the gnss satellites
            at a given time. Return the minimum
            http://en.wikipedia.org/wiki/Great-circle_distance
        '''
        pass
        phi1 = dec # Dec
        lam1 = ra # RA
        
        dist = []
        
        for sou in sorted(self.sources.keys()):
            # 1 - because 0 is yesterday due to gps-utc diff!
            t0 = datetime.datetime(*map(int,self.sources[sou].gtrs[1,0:3]))
            ut = (time - t0).total_seconds()/86400.
            
            t = self.sources[sou].UT
            x, _ = lagint(11, t, self.sources[sou].gcrs[:,6], ut)
            y, _ = lagint(11, t, self.sources[sou].gcrs[:,7], ut)
            z, _ = lagint(11, t, self.sources[sou].gcrs[:,8], ut)
            
            ra_ut = np.arctan2(y, x) # right ascention
            dec_ut = np.arctan(z/np.sqrt(x**2+y**2)) # declination
            if ra_ut < 0: ra_ut += 2.0*np.pi
            
            phi2 = dec_ut # Dec
            lam2 = ra_ut # RA
                
            d = np.arccos(np.sin(phi1)*np.sin(phi2) + \
                          np.cos(phi1)*np.cos(phi2)*np.cos(lam2-lam1))

#            d = np.arctan2(np.sqrt((np.cos(phi2)*np.sin(np.abs(lam2-lam1)))**2 + \
#                       (np.cos(phi1)*np.sin(phi2) - \
#                       np.sin(phi1)*np.cos(phi2)*np.cos(np.abs(lam2-lam1)))**2 ), \
#                       np.sin(phi1)*np.sin(phi2) + \
#                       np.cos(phi1)*np.cos(phi2)*np.cos(np.abs(lam2-lam1)))
#
#            d = np.arctan2(np.sqrt((np.cos(phi2)*np.sin(lam2-lam1))**2 + \
#                       (np.cos(phi1)*np.sin(phi2) - \
#                       np.sin(phi1)*np.cos(phi2)*np.cos(lam2-lam1))**2 ), \
#                       np.sin(phi1)*np.sin(phi2) + \
#                       np.cos(phi1)*np.cos(phi2)*np.cos(lam2-lam1))

            dist.append([sou,  float(d*180.0/np.pi)])
        
#        print dist
        
        dist = np.array(dist)
        im = np.argmin(np.array(dist[:,1], dtype=np.float))
            
#        print dist[im,:]
        return dist[im,:]


if __name__ == '__main__':
    vexfile = '_obs/ts026.vex'
#    vexfile = '_obs/g0402.vex'
    
    vex = Vex(vexfile)
    scans = [s for s in vex['SCHED']]
    
    sources = {}
    
    for s in vex['SOURCE']:
        c = SkyCoord(vex['SOURCE'][s]['ra'], \
                     vex['SOURCE'][s]['dec'], frame='icrs')
        sources[s] = [c.ra.rad, c.dec.rad]
    
    start_times = []
    for scan in scans:
        s = vex['SCHED'][scan]['source']
        t = datetime.datetime.strptime(vex['SCHED'][scan]['start'],\
                                                  '%Yy%jd%Hh%Mm%Ss')
        start_times.append([s, t])

    # load gnss ephemerides
    gnss = gnss(start_times[0][1], start_times[-1][1], glonass=False)
    
    # get the satellite for each scan
    ii=1
    src = []
    for (s,t) in start_times:
        ra, dec = sources[s]
        print ii, s,t
        gs, d = gnss.minDistance(ra, dec, t)
        print gc, d
        src.append(gs)
        ii += 1
    print sorted(set(src))