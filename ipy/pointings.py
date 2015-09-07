# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 22:52:40 2014

@author: oasis
"""

#bd, ur, sh, t6, km, ww. 12:30 - 14:30 UTC, interval 1 min

from pypride.vintlib import *
import datetime
#from time import time as _time

from astropy.coordinates import Angle#, SkyCoord, ICRS
from astropy import units as u
#from astropy.time import Time

''' load input sittings: '''
inp = inp_set('inp.cfg')
inp = inp.get_section('all')

const = constants()

#source = 'Gaia'
#sou_type = 'S'
#source = 'RA'
#sou_type = 'R'
source = 'MEX'
sou_type = 'S'
#orb_file = 'sc_eph/raw_gaia/ORB1_20140117_000001.txt'

# load full-length ephemeris
#eph_full = load_scf(orb_file)

''' stations '''
#stations = ['GEOCENTR', \
#            'BADARY', 'HART15M', 'URUMQI', 'SESHAN25', 'TIANMA65', \
#            'MEDICINA', 'KUNMING', 'WARK12M', 'ZELENCHK', 'KVNUS', \
#            'EFLSBERG', 'YEBES40M', 'METSAHOV']
#stations = ['EFLSBERG', 'ZELENCHK', 'BADARY', 'METSAHOV', 'SESHAN25', \
#            'KUNMING']
#stations = ['ZELENCHK', 'BADARY', 'YEBES40M']
stations = ['GBT-VLBA']

''' time slot '''
#date_t_start = datetime.datetime(2015,4,11,21,0,0)
#date_t_stop  = datetime.datetime(2015,4,11,22,0,0)
#t_step = 15 # seconds
date_t_start = datetime.datetime(2013,12,28,23,57,0)
date_t_stop  = datetime.datetime(2013,12,29,0,1,0)
t_step = 10 # seconds

''' fake obs '''
ob = obs(['DUMMY'],'DUMMY','C')
ob.addScan(date_t_start, t_step, stop=date_t_stop)
print ob
''' update/(down)load eops, meteo and iono data '''
if internet_on(): # check internet connection
    try:
        doup(inp['do_trp_calc'], inp['do_ion_calc'], \
             inp['cat_eop'], inp['meteo_cat'], inp['ion_cat'],\
             date_t_start, date_t_stop, inp['iono_model'])
    except Exception, err:
        print str(err)
        print 'catalogue updates failed'

''' load cats '''
_, sta, eops = load_cats(inp, source, sou_type, stations, date_t_start)

''' calculate site positions in geodetic coordinate frame
+ transformation matrix VW from VEN to the Earth-fixed coordinate frame '''
for st in sta:
    st.geodetic(const) 

# load ephemeris
eph = load_sc_eph(sou_type, source, date_t_start, date_t_stop, inp)

''' calculate [Az, El] UTC series for the S/C case '''
## lt-corrected ITRF ephemeris of the S/C is used here for each station:
#for st in sta:
#    st.AzEl(eph.gtrs, eph.UT)
    
#%% 
''' actual pointings '''

pointingsJ2000 = [] # RA/Decs at J2000
pointingsDate = []  # apparent RA/Decs (precessed and nutated to date)
azels = []

mjd_start = mjuliandate(ob.tstamps[0].year,\
                               ob.tstamps[0].month,ob.tstamps[0].day)

for tstamp in ob.tstamps:
    
    ''' set dates: '''
    mjd = mjuliandate(tstamp.year, tstamp.month, tstamp.day)
    dd = mjd - mjd_start
    UTC = (tstamp.hour + tstamp.minute/60.0 + tstamp.second/3600.0)/24.0
    JD = mjd + 2400000.5

    ''' compute tai & tt '''
    TAI, TT = taitime(mjd, UTC)

    ''' interpolate eops to tstamp '''
    UT1, eop_int = eop_iers(mjd, UTC, eops)

    ''' compute coordinate time fraction of CT day at 1st observing site '''
    CT, dTAIdCT = t_eph(JD, UT1, TT, sta[0].lon_gcen, sta[0].u, sta[0].v)

    ''' BCRS state vectors of celestial bodies at JD+CT, [m, m/s]: '''
    ## Earth:
    rrd = pleph(JD+CT, 3, 12, inp['jpl_eph'])
    earth = np.reshape(np.asarray(rrd), (3,2), 'F') * 1e3
    # Earth's acceleration in m/s**2:
    v_plus = np.array(pleph(JD+CT+1.0/86400.0, 3, 12, inp['jpl_eph'])[3:])
    v_minus = np.array(pleph(JD+CT-1.0/86400.0, 3, 12, inp['jpl_eph'])[3:])
    a = (v_plus - v_minus)*1e3 / 2.0
    a = np.array(np.matrix(a).T)
    earth = np.hstack((earth, a))
    ## Sun:
    rrd = pleph(JD+CT, 11, 12, inp['jpl_eph'])
    sun = np.reshape(np.asarray(rrd), (3,2), 'F') * 1e3
    ## Moon:
    rrd = pleph(JD+CT, 10, 12, inp['jpl_eph'])
    moon = np.reshape(np.asarray(rrd), (3,2), 'F') * 1e3

    ''' rotation matrix IERS '''
    r2000 = ter2cel(tstamp, eop_int, dTAIdCT, 'iau2000')

    ''' displacements due to geophysical effects '''
    for st in sta:
        if st.name == 'GEOCENTR' or st.name == 'RA':
            continue
        st = dehanttideinel(st, tstamp, earth, sun, moon, r2000)
        st = hardisp(st, tstamp, r2000)
        st = poletide(st, tstamp, eop_int, r2000)
    
    ''' add up geophysical corrections and convert sta state to J2000 '''
    for st in sta:
        if st.name == 'GEOCENTR' or st.name == 'RA':
            continue
        st.j2000gp(r2000)
    
    ''' calculate [Ra, Dec] for cross-eyed scheduling: '''
    pnt_J2000_sta = []
    pnt_Date_sta = []
    azel_sta = []
    for st in sta:
        # J2000 RA/Dec:
#        _, ra, dec = st.LT_radec(eph.gcrs, eph.UT, JD, UTC+dd, inp['jpl_eph'])
        _, ra, dec = st.LT_radec_bc(eph.bcrs[0], eph.CT, JD, CT+dd, inp['jpl_eph'])
        pnt_J2000_sta.append([ra, dec])
#        print st.name, ra, dec
#        st.LT_radec_obsolete(eph.fGcrs, UTC)
        # RA/Dec to date:
        xyz2000 = sph2cart(np.array([1.0,dec,ra]))
        rDate = iau_PNM00A(JD, TT)
        xyzDate = np.dot(rDate, xyz2000)
        dec, ra = cart2sph(xyzDate)[1:]
#        print st.name, ra, dec
        pnt_Date_sta.append([ra, dec])
        if st.name == 'GEOCENTR' or st.name == 'RA':
            azel_sta.append([0, 0])
        else:
#            az, el = st.AzEl(eph.gtrs, JD, eph.UT, UTC+dd, inp['jpl_eph'])
            az, el, _ = st.AzEl2(eph.gcrs, eph.UT, JD, UTC,\
                                     r2000[:,:,0], inp['jpl_eph'])
            azel_sta.append([az*180/np.pi, el*180/np.pi])

    pointingsJ2000.append(pnt_J2000_sta)
    pointingsDate.append(pnt_Date_sta)
    azels.append(azel_sta)
    

pointingsJ2000 = np.array(pointingsJ2000)
pointingsDate = np.array(pointingsDate)
azels = np.array(azels)

#%% 
''' save to files '''

# gc is the first station, so skip it
#stations_short = shname(stations[1:], inp['shnames_cat'])
stations_short = shname(stations, inp['shnames_cat'], inp['shnames_cat_igs'])

date_string = date_t_start.strftime("%y%m%d")

for jj, stash in enumerate(stations_short):
    print 'Outputting {:s}'.format(stash)
    # save AzEls to human-readable text-files
    with open(inp['out_path']+'/azel.'+source.lower()+'.'+\
              date_string+'.'+stash.lower(),'w') as f:
        for ii, tstamp in enumerate(ob.tstamps):
            line = str(tstamp)+ '  '
            az = azels[ii,jj,0]
            el = azels[ii,jj,1]
            azel = np.hstack((az, el))
            azel_str = \
               'az = {:010.6f}  el = {:10.6f}\n'\
               .format(*azel) # '{:06.2f} {:6.2f}\n'\
            line += azel_str
            f.write(line)
    # save AzEls to human-readable text-files in GreenBank format
    if stash.lower()=='gb':
        with open(inp['out_path']+'/azel.gbt.'+source.lower()+'.'+\
                  date_string+'.'+stash.lower(),'w') as f:
            f.write('# {:s} tracking table (angles in degrees)\n'.\
                    format(source.upper()))
            f.write('format=ephemeris\n')
            f.write('name={:s}\n'.format(source.upper()))
            f.write('coordmode=azel\n')
            f.write('head=date    utc        az          el\n')
            for ii, tstamp in enumerate(ob.tstamps):
                line = tstamp.strftime('%Y-%m-%d %H:%M:%S')+ '   '
                az = azels[ii,jj,0]
                el = azels[ii,jj,1]
                azel = np.hstack((az, el))
                azel_str = '{:10.6f}  {:10.6f}\n'.format(*azel)
                line += azel_str
                f.write(line)
    # save pointingsJ2000 (J2000 RA/Decs) to human-readable text-files
    with open(inp['out_path']+'/pointing.J2000.'+source.lower()+'.'+\
              date_string+'.'+stash.lower(),'w') as f:
        for ii, tstamp in enumerate(ob.tstamps):
            # mjd:
#            astrotime = Time(str(tstamp), format='iso', scale='utc', precision=9)
#            line = '{:.9f}'.format(astrotime.mjd) + '    '
            # utc
#            line = str(tstamp)+ '  '
            line = tstamp.strftime('%Y-%m-%dT%H:%M:%S') + '    '
#            c = ICRS(ra=pointingsJ2000[ii,jj,0], dec=pointingsJ2000[ii,jj,1], \
#                     unit=(u.rad, u.rad))
#            sc = SkyCoord(ra=pointingsJ2000[ii,jj,0], dec=pointingsJ2000[ii,jj,1], \
#                     frame="icrs", unit=(u.rad, u.rad))
            ra = Angle(pointingsJ2000[ii,jj,0], unit=u.rad)
            dec = Angle(pointingsJ2000[ii,jj,1], unit=u.rad)
#            radec = np.hstack((c.ra.hms, c.dec.dms))
            radec = np.hstack((ra.hms, dec.dms))
            radec[4:] = abs(radec[4:]) # minus doesn't belong to everyone..
#            radec_str = \
#               '{:02.0f} {:02.0f} {:011.8f} {:-3.0f} {:02.0f} {:011.8f}\n'\
#               .format(*radec)
            radec_str = \
      'ra = {:02.0f}h{:02.0f}m{:010.7f}s  dec = {:-3.0f}d{:02.0f}\'{:010.7f}\"\n'\
               .format(*radec)
            line += radec_str
            f.write(line)
    # save pointingsJ2000 to Giuseppe-friendly text-files
    fout = os.path.join(inp['out_path'], ''.join((source.lower(), '_', \
                                   stash.lower(), '_', date_string, '.txt')) )
    with open(fout, 'w') as f:
        for ii, tstamp in enumerate(ob.tstamps):
            line = ''.join((' source=\'', tstamp.strftime("%H%M%S"), '\' \t '))
#            c = ICRS(ra=pointingsJ2000[ii,jj,0], dec=pointingsJ2000[ii,jj,1], \
#                     unit=(u.rad, u.rad))
#            radec = np.hstack((c.ra.hms, c.dec.dms))
            ra = Angle(pointingsJ2000[ii,jj,0], unit=u.rad)
            dec = Angle(pointingsJ2000[ii,jj,1], unit=u.rad)
            radec = np.hstack((ra.hms, dec.dms))
            radec[4:] = abs(radec[4:]) # minus doesn't belong to everyone..
            radec_str = 'ra={:02.0f}:{:02.0f}:{:05.2f} \t '.format(*radec[0:3]) 
            radec_str += 'dec={:+.0f}:{:02.0f}:{:05.2f} \t '.format(*radec[3:])
            radec_str +=  'equinox=\'j2000\' /\n'
            line += radec_str
            f.write(line)
    # save pointingsDate (apparent RA/Decs) to human-readable text-files
    with open(inp['out_path']+'/pointing.apparent.'+source.lower()+'.'+\
              date_string+'.'+stash.lower(),'w') as f:
        for ii, tstamp in enumerate(ob.tstamps):
            # mjd:
#            astrotime = Time(str(tstamp), format='iso', scale='utc', precision=9)
#            line = '{:.9f}'.format(astrotime.mjd) + '    '
            # utc
#            line = str(tstamp)+ '  '
            line = tstamp.strftime('%Y-%m-%dT%H:%M:%S') + '    '
#            c = ICRS(ra=pointingsDate[ii,jj,0], dec=pointingsDate[ii,jj,1], \
#                     unit=(u.rad, u.rad))
#            radec = np.hstack((c.ra.hms, c.dec.dms))
            ra = Angle(pointingsDate[ii,jj,0], unit=u.rad)
            dec = Angle(pointingsDate[ii,jj,1], unit=u.rad)
            radec = np.hstack((ra.hms, dec.dms))

            radec[4:] = abs(radec[4:]) # minus doesn't belong to everyone..
#            radec_str = \
#               '{:02.0f} {:02.0f} {:011.8f} {:-3.0f} {:02.0f} {:011.8f}\n'\
#               .format(*radec)
            radec_str = \
      'ra = {:02.0f}h{:02.0f}m{:010.7f}s  dec = {:-3.0f}d{:02.0f}\'{:010.7f}\"\n'\
               .format(*radec)
            line += radec_str
            f.write(line)
    