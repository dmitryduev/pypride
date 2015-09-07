# -*- coding: utf-8 -*-
"""
Create input config file with all the stuff necessary for vispy
Created on Tue Sep 17 18:59:22 2013

@author: oasis
"""

import ConfigParser
from datetime import datetime

config = ConfigParser.RawConfigParser(allow_no_value = True)
config.optionxform = str # preserve capital letters in comments

config.add_section('Catalogues')
config.set('Catalogues', '# absolute path to pyp. gets appended to paths below in case they are relative.')
config.set('Catalogues', 'abs_path', '/Users/oasis/_jive/python/vispy')
config.set('Catalogues', '# site positions [compatible with calc and vtd]:')
config.set('Catalogues', 'sta_xyz', 'cats/glo.sit') # site positions
config.set('Catalogues', '# site velocities [compatible with calc and vtd]:')
config.set('Catalogues', 'sta_vxvyvz', 'cats/glo.vel') # site velocities
config.set('Catalogues', '# telescope axis offsets:')
config.set('Catalogues', 'sta_axo', 'cats/rfc_2011a.axof') # axis offsets
config.set('Catalogues', '# ocean loading:')
config.set('Catalogues', 'oc_load', 'cats/tpxo72.blq') # ocean loading
config.set('Catalogues', '# atmospheric loading:')
config.set('Catalogues', 'atm_load', 'cats/atmload.cat') # atmospheric loading
config.set('Catalogues', '# alternative site names:')
config.set('Catalogues', 'sta_nam', 'cats/station.names') # alternative site names
config.set('Catalogues', '# eops IAU2000 (C04 series):')
config.set('Catalogues', 'cat_eop', 'cats/eopc04.cat') # eops IAU2000
config.set('Catalogues', '# thermal deformation coeffs + axis types:')
config.set('Catalogues', 'sta_thermdef', 'cats/antenna-info.txt') # T def coeffs+axis types
config.set('Catalogues', '# source positions [compatible with calc and vtd]:')
config.set('Catalogues', 'source_cat', 'cats/glo.src') # source positions
config.set('Catalogues', '# source names [compatible with calc and vtd]:')
config.set('Catalogues', 'source_nam', 'cats/source.names') # source names
config.set('Catalogues', '# stations 2-letter codes:')
config.set('Catalogues', 'shnames_cat', 'cats/ant.shn') # stations short names
config.set('Catalogues', 'shnames_cat_igs', 'cats/ns-codes.txt') # stations short names
config.set('Catalogues', '# meteo data folder:')
config.set('Catalogues', 'meteo_cat', 'meteo') # meteo data folder
config.set('Catalogues', '# ionospheric data folder:')
config.set('Catalogues', 'ion_cat', 'ion') # ionospheric data folder
config.set('Catalogues', '# deep space s/c freq 3-way ramping (+ s/c name):')
config.set('Catalogues', 'f_ramp', 'cats/ramp.') # freq 3-way ramping + s/c name
config.set('Catalogues', '# deep space s/c freq 1-way ramping (+ s/c name):')
config.set('Catalogues', 'f_ramp1w', 'cats/ramp1w.') # freq 3-way ramping + s/c name
config.set('Catalogues', '# Earth satellites Tx frequencies:')
config.set('Catalogues', 'f_gc', 'cats/sc.freq') # freq 3-way ramping + s/c name

config.add_section('Ephemerides')
config.set('Ephemerides', '# Solar System ephemerides in JPL format: '+\
            u'[\'JPLEPH.403\', \'JPLEPH.405\', \'JPLEPH.421\', \'inpop13c.tdb\']')
config.set('Ephemerides', 'jpl_eph', 'jpl_eph/JPLEPH.421') # JPL DE/LE EPHEMERIDES
config.set('Ephemerides', '# spacecraft ephemerides folder:')
config.set('Ephemerides', 'sc_eph_cat', 'sc_eph') # S/C ephemerides

config.add_section('Directories')
config.set('Directories', '# input (vex-files) folder:')
config.set('Directories', 'obs_path', '_obs')
config.set('Directories', '# output folder:')
config.set('Directories', 'out_path', '_out')

config.add_section('Models')
# Where to place the phase center (station name)
config.set('Models', '# phase center position (station name):')
config.set('Models', 'phase_center', 'GEOCENTR')
# Near-field VLBI delay model:
config.set('Models', '# near field VLBI delay model: [\'Moyer\', \'Fukushima\']')
config.set('Models', 'nf_model', 'Moyer') # Fukushima, Kopeikin or Moyer
#tropo and iono:
config.set('Models', '# tropospheric model: [\'wien\']')
config.set('Models', 'tropo_model', 'wien') # wien (vmf1), TODO: merra or ecmwf
config.set('Models', '# ionospheric IGS TEC maps (final or rapid): [\'igs\', \'igr\']')
config.set('Models', 'iono_model', 'igs') # igs (final) or igr (rapid)
# Doppler prediction model
config.set('Models', '# Doppler prediction model: '+\
            '[\'none\', \'bary1way\', \'geo1way\', \'bary3way\']')
config.set('Models', 'dop_model', 'bary3way') # none, bary1way, geo1way, bary3way
# Doppler mode parametres for 2(3)-way
#config.set('Models', 'uplink_sta', 'None') # Uplink transmitting station
#config.set('Models', 'freq_type', 'proper') # frequency value type ('proper', 'gc')
#config.set('Models', 'freq', '8400000000.2827155') # Uplink frequency, Hz
# set to 1 for 1-way Doppler
config.set('Models', '# Turnaround ratio for 2(3)-way Doppler: [880/749]')
config.set('Models', 'TR', '1.1748998664886516') # Turnaround ratio 
# Jacobians
config.set('Models', '# step size for Jacobian (UVW in near field) calculation:')
config.set('Models', 'mas_step', '1000.0')  # mas
config.set('Models', 'm_step', '10000.0')  # m
# generate additional ephemerides and calc delays for RA/GNSS XYZ position correction
config.set('Models', 'm_step_xyz', '1000.0')  # m

config.add_section('Switches')
''' These are default values (don't calculate anything) '''
#tropo and iono:
config.set('Switches', '# calculate troposphere or not?:')
config.set('Switches', 'do_trp_calc', 'True') # calculate troposphere or not?
config.set('Switches', '# calculate tropospheric gradients or not?:')
config.set('Switches', 'do_trp_grad_calc', 'False') # calculate tropospheric gradients or not?
config.set('Switches', '# calculate ionosphere or not?:')
config.set('Switches', 'do_ion_calc', 'True') #  calculate ionosphere or not?
# Calculate Doppler prediction AND/OR delay prediction?
config.set('Switches', '# calculate delays?:')
config.set('Switches', 'delay_calc', 'False')
# Calculate uvws?
config.set('Switches', '# calculate UVWs/Jacobians?:')
config.set('Switches', 'uvw_calc', 'False')
# note that 2nd station is the transmitter if 2(3)-way
config.set('Switches', '# calculate Doppler?:')
config.set('Switches', '# note that 2nd station is the transmitter if 2(3)-way')
config.set('Switches', 'doppler_calc', 'False')
# Jacobians
# generate additional ephemerides and calc delays for S/C position correction
config.set('Switches', '# generate additional ephemerides and calc delays for S/C position correction?:')
config.set('Switches', 'sc_rhophitheta', 'False')
# generate additional ephemerides and calc delays for RA/GNSS XYZ position correction
config.set('Switches', '# generate additional ephemerides and calc delays for RA/GNSS XYZ position correction?:')
config.set('Switches', 'sc_xyz', 'False')
config.set('Switches', '# force update s/c ephemeris')
config.set('Switches', 'sc_eph_force_update', 'False')

# Write our configuration file to 'inp.cfg'
now = datetime.strftime(datetime.now(), '%Y%m%d%H%M%S')

with open('inp.{:s}.cfg'.format(now), 'wb') as configfile:
    config.write(configfile)