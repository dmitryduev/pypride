#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 15:18:45 2014

@author: oasis
"""

from pypride.vintlib import *
import os
import multiprocessing as mp
#import matplotlib.pyplot as plt
#import prettyplotlib as ppl
#from time import time as _time
import argparse

'''
#==============================================================================
# Function returning slant tecs (modified ion_igs)
#==============================================================================
'''
def ion_tec(sta, iono, elv, azi, UT, f_0=None):
    '''
    Calculate ionospheric delay using IGS vertical TEC maps
    '''
    # calculate distance to point J of ray pierce into ionosphere from site
#    H = 438.1*1e3 # m - mean height of the ionosphere above R_E
    H = 450*1e3 # m - height of ionosphere above R_E as stated in IONEX files
    R_E = 6371.0*1e3 # m - Earth's radius from the TEC map 

    alpha = 0.9782
    
    lat_geod = sta.lat_geod
    lon_gcen = sta.lon_gcen
    h_geod = sta.h_geod
    
    UT_tec = iono.UT_tec
    fVTEC = iono.fVTEC
    
    #WGS84 Ellipsoid
    a = 6378137.0 # m
    f = 1.0/298.2572235630
    b = a*(1.0-f) # m
    ec = sqrt((a**2-b**2)/(a**2))
    
    R_oscul = a*sqrt(1.0-ec**2)/(1.0-(ec*sin(lat_geod))**2) # m
    
    source_vec = np.array([sin(pi/2.0-elv)*cos(azi),\
                           sin(pi/2.0-elv)*sin(azi),\
                           cos(pi/2.0-elv) ])
    
    # slanted distance btw the ground and the iono layer
    ds = (R_oscul+h_geod)*sin(-elv) + \
         0.5 * sqrt( (2.0*(R_oscul+h_geod)*sin(-elv))**2 - \
         4.0*((R_oscul+h_geod)**2 - (R_E+H)**2) )
    # cart crds of the starting point
    rpt = [R_oscul+h_geod, lat_geod, lon_gcen]
    r0 = sph2cart(rpt)
    # cart crds of the ionospheric pierce point
    r1 = r0 + ds*source_vec
    
    # lat/long of the pierce point
    rlalo = cart2sph(r1)
    lon = 180.0*rlalo[2]/pi
    lat = 180.0*rlalo[1]/pi
    
    # find closest epoch in TEC data:
    # easy case - UT is in UT_tec
    n0 = np.searchsorted(UT_tec, UT)
    if UT in UT_tec:
        # only need to interpolate TECz to the Lat/Lon of the pierce point
        TEC_z = float(fVTEC[n0](lon, lat))
    else:
        # else take 4 consecutive epochs and interpolate:
        N_epoch = len(fVTEC)
        if n0==1 or n0==0:
            nl = 0; nr = 4
        elif n0==N_epoch-2 or n0==N_epoch-1:
            nl = N_epoch-4; nr = N_epoch
        else:
            nl = n0-2; nr = n0+2;
        TEC_z = []
        for nn in range(nl,nr):
            TEC_z.append(float(fVTEC[nn](lon, lat)))
        # interpolate zenith TECz to the epoch of observation
#        fTEC = sp.interpolate.interp1d(UT_tec[nl:nr], TEC_z, kind='linear')
#        TEC_z = fTEC(UT)
        TEC_z = np.interp(UT, UT_tec[nl:nr], TEC_z)
    
    # calculate slanted TEC
    TEC = TEC_z / cos(asin( (R_oscul+h_geod)*sin(alpha*(pi/2.0-elv))/(R_E+H) ))
    TEC_tecu = 0.1*TEC # in TEC units
    
    # calculate ionspheric delay for the source
#    delay_ion = 5.308018e10*TEC_tecu/(4.0*pi**2*f_0*f_0)
    
    return TEC_tecu


'''
#==============================================================================
# Calculate all sorts of TECs
#==============================================================================
'''
def calctec(ins):
    ''' 
    Calculate all sorts of TECs
    '''
    # parse single-variable input:
    record, sta_r, sta_t, sou_type, source, \
                ip_tecs, const, inp, tec_uplink = ins

    t = datetime.datetime(*map(int, record[4:10])) + \
        datetime.timedelta(seconds=record[10]/2.0)

    # t_obs - middle of obs run
    t_obs = (t.hour*3600.0 + t.minute*60.0 + record[10]/2.0)/86400.0

    ''' download iono data '''
#    doup(False, inp['do_ion_calc'], \
#         inp['cat_eop'], inp['meteo_cat'], inp['ion_cat'], \
#         t, t, 'igs')
    
    ''' load vtec maps '''
    try:
        if record[34]==30.0 or record[35]==30.0 or record[36]==1.0:
            iono = ion(t, t, inp)
        ''' calculate tec for downlink '''
        if record[34]==30.0:
            st = sta_r[int(record[3])-1] # receiving station
            el = record[18]*pi/180
            if record[17]>180:
                az = (record[17]-360.0)*pi/180
            else:
                az = record[17]*pi/180
            tec = ion_tec(st, iono, el, az, t_obs)
    #            print 'TEC down:', record[34], '{:4.1f}'.format(tec)
            # replace corresponding field in the record:
            record[34] = tec
        ''' calculate tec for uplink '''
        if record[35]==30.0:
            st = sta_t[int(record[33])-1] # transmittimg station
            # use planetary ephems
            mjd = mjuliandate(t.year, t.month, t.day)
            UTC = (t.hour + t.minute/60.0 + t.second/3600.0)/24.0
            
            JD = mjd + 2400000.5
            TAI, TT = taitime(mjd, UTC)
            with open(inp['cat_eop'], 'r') as fc:
                fc_lines = fc.readlines()
            eops = np.zeros((7,7)) # +/- 3 days
            for jj in range(len(fc_lines)):
                if fc_lines[jj][0]!=' ' and fc_lines[jj][0]!='*':
                    entry = [float(x) for x in fc_lines[jj].split()]
                    if len(entry) > 0 and entry[3] == np.floor(mjd) - 3:
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
            UT1, eop_int = eop_iers(mjd, UTC, eops)
            CT, dTAIdCT = t_eph(JD, UT1, TT, st.lon_gcen, st.u, st.v)
            
            # Earth:
            rrd = pleph(JD+CT, 3, 12, inp['jpl_eph'])
            earth = np.reshape(np.asarray(rrd), (3,2), 'F') * 1e3
            # Venus/Mars
            if source.lower()=='vex':
                rrd = pleph(JD+CT, 2, 12, inp['jpl_eph'])
            elif source.lower()=='mex':
                rrd = pleph(JD+CT, 4, 12, inp['jpl_eph'])
            planet = np.reshape(np.asarray(rrd), (3,2), 'F') * 1e3
            lt = np.linalg.norm(planet[:,0]-earth[:,0])/const.C
            dt = 3*datetime.timedelta(seconds=lt) # 2-way LT

            if tec_uplink=='planet':
                # 3 LTs ago (2LTs - signal round trip + another LT )
    #            print dt
    #            print CT
                CT -= dt.total_seconds()/86400.0
                _, eop_int = eop_iers(mjd, UTC-2.0*lt/86400.0, eops)
    #            print CT
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
                # Venus/Mars
                if source.lower()=='vex':
                    rrd = pleph(JD+CT, 2, 12, inp['jpl_eph'])
                elif source.lower()=='mex':
                    rrd = pleph(JD+CT, 4, 12, inp['jpl_eph'])
                planet = np.reshape(np.asarray(rrd), (3,2), 'F') * 1e3
                
                r2000 = ter2cel(t-dt, eop_int, dTAIdCT, 'iau2000')
                st = dehanttideinel(st, t-dt, earth, sun, moon, r2000)
                st = hardisp(st, t-dt, r2000)
                st = poletide(st, t-dt, eop_int, r2000)
                st.j2000gp(r2000)
    
                r = planet[:,0] - (st.r_GCRS + earth[:,0])
                ra  = np.arctan2(r[1],r[0]) # right ascention
                dec = np.arctan(r[2]/np.sqrt(r[0]**2+r[1]**2)) # declination
                if ra < 0: ra += 2.0*np.pi
                K_s = np.array([cos(dec)*cos(ra), \
                                cos(dec)*sin(ra), \
                                sin(dec)])
    
                az, el = aber_source(st.v_GCRS, st.vw, K_s, r2000, earth)
    #            print map(lambda x:x*180/pi, (az,el))
    #            print t_obs, 2.0*lt/86400.0
                tec = ion_tec(st, iono, el, az, (t_obs - 2.0*lt/86400.0))
#                print t, tec
            
            elif tec_uplink=='sc':
                # load s/c ephemeris:
    #            tic = _time()
                eph = load_sc_eph(sou_type, source, t, t, inp)
    #            print 'loading eph took', _time()-tic, ' s'
                # lt to station in seconds at t_obs
                lt, _, _ = st.LT_radec_bc(eph.bcrs[0], eph.CT, JD, UTC, inp['jpl_eph'])
                # az/el 2LT ago (another LT is accounted for internally)
                r2000 = ter2cel(t, eop_int, dTAIdCT, 'iau2000')
                az, el = st.AzEl2(eph.gtrs, eph.UT, JD, \
                                    UTC - 2.0*lt/86400.0, inp['jpl_eph'])
    #            print map(lambda x:x*180/pi, (az,el))
    #            print t_obs, 2.0*lt/86400.0
                tec = ion_tec(st, iono, el, az, (t_obs - 2.0*lt/86400.0))
#                print t, tec
    #            print 'TEC up:', record[35], '{:4.1f}'.format(tec)
            # replace corresponding field in the record:
            record[35] = tec
        ''' calculate IP tec (intepolate from table) '''
        if record[36]==1.0:
            orb_phase = record[22] + record[23]
            ip_tec, _ = lagint(4, ip_tecs[:,0], ip_tecs[:,3], orb_phase)
            record[36] = ip_tec*(1+880.0/749.0)
                
    except Exception, err:
        print err
        print 'error occured for {:4d}/{:02d}/{:02d}'\
              .format(t.year, t.month, t.day)
    
    finally:
        return record


'''
#==============================================================================
# Run pipeline
#==============================================================================
'''
if __name__ == '__main__':
    '''
        This script is supposed to be run from the command line
    '''
    # create parser
    parser = argparse.ArgumentParser()
    # optional arguments
    parser.add_argument('-s', '--spacecraft', type=str, choices=['vex', 'mex'],
                        help='spacecraft')
    parser.add_argument('-u', '--uppoint', type=str, choices=['sc', 'planet'],
                        default='planet',
                        help='where to point when calculating TEC on uplink:'+\
                        ' \'planet\' to point at the S/C host planet (default),'+\
                        ' \'sc\' to point at the S/C')
    parser.add_argument('-i', '--ionomodel', type=str, choices=['igs', 'igr'],
                        default='igr',
                        help='IGS\' ionospheric TEC model to use: final or '+\
                        'rapid (default)')
    parser.add_argument('-p', '--parallel', action='store_true',
                        help='run computation in parallel mode')
    # positional argument
    parser.add_argument('inpFile', type=str,
                        help="input ScintObsSummary table")    
    args = parser.parse_args()
    scint_table_file = args.inpFile
    
    # which spacecraft?
    if args.spacecraft=='vex':
        source = 'vex'
    elif args.spacecraft=='mex':
        source = 'mex'
    else:
        # try to guess:
        if 'vex' in scint_table_file.lower():
            source = 'vex'
        elif 'mex' in scint_table_file.lower():
            source = 'mex'
        else:
            raise Exception('Spacecraft not set; failed to guess.')

    # where to point uplink?
    if args.uppoint=='sc':
        tec_uplink = 'sc'
    elif args.uppoint=='planet':
        tec_uplink = 'planet'
    else:
        # let's do it quickly by default
        tec_uplink = 'planet'

    # proceed with pipelining:

    inp = inp_set('inp.cfg')
    inp.iono_model = args.ionomodel
    inp = inp.get_section('all')
    const = constants()
    
    receivers = ['METSAHOV', 'MEDICINA', 'MATERA', 'NOTO', 'WETTZELL', \
                 'YEBES40M', 'PUSHCHIN', 'ONSALA60', 'HARTRAO', \
                 'SESHAN25', 'KUNMING', 'URUMQI', 'HART15M', \
                 'SVETLOE', 'ZELENCHK', 'BADARY', 'TIANMA65', 'WARK12M',\
                 'HOBART26', 'HOBART12', 'YARRA12M', 'KATH12M',\
                 'WARK30M', 'WETTZ13S', 'WETTZ13N']
    
    transmitters = ['NWNORCIA', 'CEBREROS', 'MALARGUE', 'TIDBIN64', 'DSS35',\
                    'DSS45', 'DSS34', 'DSS65', 'DSS63', 'DSS14', 'DSS15']

    recv_short = shname(receivers, inp['shnames_cat'], inp['shnames_cat_igs'])
    tran_short = shname(transmitters, inp['shnames_cat'], inp['shnames_cat_igs'])

    ''' load cats '''
    # last argument is dummy
    sou_type = 'S'
    _, sta_r, _ = load_cats(inp, source, sou_type, receivers, \
                            datetime.datetime.now())
    _, sta_t, _ = load_cats(inp, source, sou_type, transmitters, \
                            datetime.datetime.now())
    
    ''' calculate site positions in geodetic coordinate frame
    + transformation matrix VW from VEN to the Earth-fixed coordinate frame '''
    for st_r in sta_r:
        st_r.geodetic(const)
    for st_t in sta_t:
        st_t.geodetic(const)
    
    ''' load IP TEC table '''
    if source == 'vex':
        f_iptec = 'TecVenus.NFS.txt'
    if source == 'mex':
        f_iptec = 'TecMars.NFS.txt'
    with open(os.path.join(inp['ion_cat'], 'scint', f_iptec)) as fipt:
        ip_tecs = np.array([map(float,x.split()) for x in fipt.readlines() \
                            if x[0]!='#'])


    ''' Parse Scint Table '''
    #    with open(os.path.join(inp['ion_cat'], 'scint', scint_table_file)) as f:
    with open(scint_table_file) as f:
        f_lines = f.readlines()
    
    scintTable = []
    f_lines_clean = [l for l in f_lines if l[0]=='|' and len(l)>20]
    for line in f_lines_clean:
        line_parsed = map(float, line.replace('|',' ').split())
        scintTable.append(line_parsed)
    scintTable = np.array(scintTable)


    ''' precompute ephemerides for faster access '''
    # get unique dates in the scintTable
    dates = np.array([datetime.datetime(*map(int, x[4:7])) for x in scintTable])
    dates = np.unique(dates)
    
    # planet position could be used instead!
    if tec_uplink == 'sc':
        # now get "min" and "max" epochs on that day
        startStopTimes = []
        for date in dates:
            span = [x[4:10] for x in scintTable \
                    if datetime.datetime(*map(int,x[4:7])) == date]
            mjds = []
            for spa in span:
                mjds.append(mjuliandate(*spa))
            mjd_min = np.argmin(mjds) # "min" epoch
            mjd_max = np.argmax(mjds) # "max" epoch
        
            start = datetime.datetime(*map(int, span[mjd_min]))
            stop = datetime.datetime(*map(int, span[mjd_max]))
            startStopTimes.append([start, stop])
        #    print date, start, stop
        # create one ephemeris spanning across min-max
        # (not to recalculate it at each step)
        print 'Creating ephemerides...'
        
        # do it in a parallel way!
        def ephmakerWrapper(args):
            ''' wrapper func to unpack arguments '''
            return load_sc_eph(*args)
            
        # check raw eph files boundaries:
        path = os.path.join(inp['sc_eph_cat'], 'raw_'+source.lower())
        checkbound(source, orb_path=path)
        ## make single-var inputs:
        #ins = []
        #for (start, stop) in startStopTimes:
        #    ins.append((sou_type, source, start, stop, inp))
        ## number of cores available
        #n_cpu = mp.cpu_count()
        ## create pool
        #pool = mp.Pool(n_cpu)
        ## asyncronously apply calctec to each of ins
        #pool.map_async(ephmakerWrapper, ins)
        ## close bassejn
        #pool.close() # we are not adding any more processes
        #pool.join() # wait until all threads are done before going on
        
        ## do it in serial way!
        for ii, (start, stop) in enumerate(startStopTimes):
            print len(startStopTimes)-ii, 'ephemerides to go'
            load_sc_eph(sou_type, source, start, stop, inp, load=False)
    
    ''' download iono data if necessary '''
    print 'Fetching ionospheric data...'
    #tec_model = 'igs' # final 'igs' or rapid 'igr' IGS solution
    tec_model = inp['iono_model']
    for t in dates:
        doup(False, inp['do_ion_calc'], \
             inp['cat_eop'], inp['meteo_cat'], inp['ion_cat'], \
             t, t, tec_model)

    ''' iterate over records in scint table '''
    print 'Computing TECs...'
    # make single-var inputs:
    ins = []
    for record in scintTable:
        ins.append((record, sta_r, sta_t, sou_type, source, \
                    ip_tecs, const, inp, tec_uplink))
    ## parallel way:
    if args.parallel:
        # number of cores available
        n_cpu = mp.cpu_count()
        # create pool
        pool = mp.Pool(n_cpu)
        # asyncronously apply calctec to each of ins
        result = pool.map_async(calctec, ins)
        # close bassejn
        pool.close() # we are not adding any more processes
        pool.join() # wait until all threads are done before going on
        # get the ordered results back into obsz
        outs = result.get()
        scintTableOut = []
        for ou in outs:
            scintTableOut.append(ou)
    ## serial way:
    else:    
        scintTableOut = []
        for ii, record in enumerate(ins):
            print len(scintTable)-ii, 'records to go'
            scintTableOut.append(calctec(record))


    #%%
    '''
    #==========================================================================
    # Create output table
    #==========================================================================
    '''
    print 'Outputting...'
    # date string
    #date_string = datetime.datetime.now().strftime("%y%m%d")
    #out_table = ''.join(('ScintTable.', date_string))

    # put in the vispy output folder
    if '/' in scint_table_file:
        slash = [i for i,x in enumerate(scint_table_file) if x=='/']
        out = scint_table_file[slash[-1]+1:]
    else:
        out = scint_table_file
    
    out_table = ''.join((out[:-4], 'i', out[-4:]))
    
    first_entry = [i for i, x in enumerate(f_lines) if x[0]=='|'][0]
    header = [x for x in f_lines[:first_entry] if x[0]=='/']
    last_entry = [i for i, x in enumerate(f_lines) if x[0]=='|'][-1]
    footer = [x for x in f_lines[last_entry:] if x[0]=='/']
    
    with open(os.path.join(inp['out_path'], out_table), 'w') as f:
        # print header:
        for line in header:
            f.write(''.join((line.strip(), '\n')))
    
        # print table:
        for ii, record in enumerate(scintTableOut):
            line = '|{:5d}'.format(ii+1)
            line += '{:4d}{:4d}{:4d}{:8d} {:02d} {:02d}  {:02d} {:02d} {:02d}'\
                    .format(*map(int, record[1:10]) )
            line += '{:6d}{:4d} {:02d}'\
                    .format(*map(int, record[10:13]))
            line += ' {:04.1f}'\
                    .format(float(record[13]))
            line += '{:4d} {:02d}'\
                    .format(*map(int,record[14:16]))
            line += ' {:04.1f}'\
                    .format(float(record[16]))
            line += '{:6.1f}{:6.1f}{:6.1f}{:6.1f}{:8.4f}{:7.2f}{:7.2f} |'\
                    .format(*map(float, record[17:24]))
            line += '{:8.3f}{:7.3f}{:8.3f}{:6.3f}{:10.2f}{:7.2f}{:8.2f}'\
                    .format(*map(float, record[24:31]))
            line += '{:8d} |{:3d}{:7d}'\
                    .format(*map(int, record[31:34]))
            line += '{:6.1f}{:8.1f}{:9.1f} |\n'\
                    .format(*map(float, record[34:])) 
    
            f.write(line)
    
        # print footer:
        for line in footer:
            f.write(''.join((line.strip(), '\n')))