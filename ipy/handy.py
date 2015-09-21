# -*- coding: utf-8 -*-
"""
Created on Sun Jan 26 23:21:44 2014

@author: oasis
"""

from classes import *
from pypride.vintlib import *
from datetime import datetime

''' load input sittings: '''
inp = inp_set('inp.cfg')

'''
#==============================================================================
# obs setup examples
#==============================================================================
'''

''' 1-way ramped Doopler from a deep space S/C'''
#t_start = datetime(2008,1,5,3,32,0,50000)
#t_stop  = datetime(2008,1,5,3,32,1,50000)
#t_step = 0.1 # seconds
#inp_swchs = inp.get_section('all') #default input switches (all False)
#inp_swchs['doppler_calc'] = True
#ob = obs(['NWNORCIA'], 'VEX', 'S', inp=inp_swchs) # freq in Hz!!
#ob.addScan(t_start, t_step, stop=t_stop)

#t_start = datetime(2013,9,18,18,8,0)
#t_stop  = datetime(2013,9,18,18,8,30)
#t_step = 1 # seconds
#ob = obs(['GEOCENTR','METSAHOV'],'W3IRS5_H2O','C')
#ob.addScan(t_start, t_step, stop=t_stop, freq=8.4e9)
#print ob
#ob.addScan(t_start, t_step, nobs=10)
#print ob

#t_start = datetime(2013,3,9,10,5,0)
#t_stop  = datetime(2013,3,10,1,0,0)
#t_step = 10 # seconds
#inp_swchs = inp.get_section('all') #default input switches (all False)
#inp_swchs['doppler_calc'] = True
#ob = obs(['MEDICINA','GEOCENTR'],'RA','R', inp=inp_swchs) # freq in Hz!!
#ob.addScan(t_start, t_step, stop=t_stop, freq=8.4e9 + 0.2827155)
           
#t_start = datetime(2014,9,21,5,0,0)
#t_stop  = datetime(2014,9,21,6,0,0)
#t_step = 10 # seconds
#inp_swchs = inp.get_section('all') #default input switches (all False)
#inp_swchs['doppler_calc'] = True
#ob = obs(['EFLSBERG','GEOCENTR'],'RA','R', inp=inp_swchs) # freq in Hz!!
#ob.addScan(t_start, t_step, stop=t_stop, freq=15e9)
           
#t_start = datetime(2015,2,15,16,0,0)
#t_stop  = datetime(2015,2,15,18,0,0)
#t_step = 60 # seconds
#inp_swchs = inp.get_section('all') #default input switches (all False)
#inp_swchs['doppler_calc'] = True
#ob = obs(['EFLSBERG','GEOCENTR'],'RA','R', inp=inp_swchs) # freq in Hz!!
#ob.addScan(t_start, t_step, stop=t_stop, freq=15e9)

#t_start = datetime(2013,12,28,23,30,0)
#t_stop  = datetime(2013,12,29,9,40,0)
#t_step = 10 # seconds
#inp_swchs = inp.get_section('all') #default input switches all False
#inp_swchs['doppler_calc'] = True
#ob = obs(['METSAHOV'],'MEX','S', inp=inp_swchs)
#ob.addScan(t_start, t_step, stop=t_stop)

#t_start = datetime(2013,12,28,22,30,0)
#t_stop  = datetime(2013,12,29,11,40,0)
#t_start = datetime(2013,12,28,17,30,0)
#t_stop  = datetime(2013,12,29,18,40,0)
#t_step = 10 # seconds
#inp_swchs = inp.get_section('all') #default input switches all False
#inp_swchs['delay_calc'] = True
#ob = obs(['GEOCENTR','HARTRAO'],'MEX','S', inp=inp_swchs)
#ob.addScan(t_start, t_step, stop=t_stop, freq=8.412e9)

#t_start = datetime(2013,12,3,1,10,0)
#t_stop  = datetime(2013,12,3,2,9,0)
#t_step = 10 # seconds
#inp_swchs = inp.get_section('all')
#inp_swchs['delay_calc'] = True
#ob = obs(['GEOCENTR','ZELENCHK'],'MEX','S', inp=inp_swchs)
#ob.addScan(t_start, t_step, stop=t_stop, freq=8.412e9)

#t_start = datetime(2012,12,5,2,39,1)
#t_stop  = datetime(2012,12,5,2,40,1)
#t_start = datetime(2012,12,5,2,38,0)
#t_stop  = datetime(2012,12,5,5,38,0)
#t_step = 10 # seconds
#inp_swchs = inp.get_section('all')
#inp_swchs['doppler_calc'] = True
#ob = obs(['HART15M'], 'VEX', 'S', inp=inp_swchs) # freq in Hz!!
#ob.addScan(t_start, t_step, stop=t_stop)

#t_start = datetime(2013,12,29,2,41,30)
#t_stop  = datetime(2013,12,29,2,43,30)
#t_step = 1 # seconds
#inp_swchs = inp.get_section('all') #default input switches all False
#inp_swchs['delay_calc'] = True
#ob = obs(['GEOCENTR', 'ONSALA60'], 'MEX', 'S', inp=inp_swchs) # freq in Hz!!
#ob.addScan(t_start, t_step, stop=t_stop, freq=8.42e9)

t_start = datetime(2013,12,28,23,59,30)
t_stop  = datetime(2013,12,29,0,1,00)
t_step = 1 # seconds
#t_start = datetime(2013,12,28,21,30,00)
#t_stop  = datetime(2013,12,28,22,30,00)
#t_step = 10 # seconds
inp_swchs = inp.get_section('all') #default input switches (all False)
inp_swchs['delay_calc'] = True
ob = obs(['GEOCENTR', 'TIANMA65'], 'MEX', 'S', inp=inp_swchs) # freq in Hz!!
ob.addScan(t_start, t_step, stop=t_stop, freq=8.42e9)

#t_start = datetime(2013,12,28,18,20,0)
#t_stop  = datetime(2013,12,29,18,30,0)
#t_step = 10 # seconds
##t_start = datetime(2013,12,28,23,59,0)
##t_stop  = datetime(2013,12,29,19,5,0)
##t_step = 60 # seconds
##t_start = datetime(2013,12,28,23,59,20)
##t_stop  = datetime(2013,12,29,0,42,30)
##t_step = 10 # seconds
#inp_swchs = inp.get_section('all')
#inp_swchs['doppler_calc'] = True
#ob = obs(['GEOCENTR'], 'MEX', 'S', inp=inp_swchs) # freq in Hz!!
#ob.addScan(t_start, t_step, stop=t_stop, freq=8.42e9)

#t_start = datetime(2013,12,29,8,20,0)
#t_stop  = datetime(2013,12,29,18,30,0)
#t_step = 10 # seconds
#inp_swchs = inp.get_section('all')
#inp_swchs['doppler_calc'] = True
#ob = obs(['KP-VLBA'], 'MEX', 'S', inp=inp_swchs) # freq in Hz!!
#ob.addScan(t_start, t_step, stop=t_stop, freq=8.42e9)

#t_start = datetime(2013,12,29,11,30,0)
#t_stop  = datetime(2013,12,29,19,5,0)
#t_start = datetime(2013,12,29,8,20,0)
#t_stop  = datetime(2013,12,29,18,30,0)
#t_step = 1 # seconds
#inp_swchs = inp.get_section('all')
#inp_swchs['delay_calc'] = True
#ob = obs(['GEOCENTR','MK-VLBA'],'MEX','S', inp=inp_swchs) # freq in Hz!!
#ob.addScan(t_start, t_step, stop=t_stop, freq=8.42e9)

#t_start = datetime(2013,12,1,2,0,0)
#t_stop  = datetime(2013,12,1,3,0,0)
#t_step = 1 # seconds
#inp_swchs = inp.get_section('all')
#inp_swchs['doppler_calc'] = True
#ob = obs(['GEOCENTR','WETTZELL'],'VEX','S',inp_swchs)
#ob.addScan(t_start, t_step, stop=t_stop, freq=8.4e9) # freq in Hz!!

ob.dude = vint_s(ob)

#%% 
''' 
#==============================================================================
# text output
#==============================================================================
'''

# date string
date_string = t_start.strftime("%y%m%d")

# output directory:
out_path = inp.out_path

# get station short names:
stations_short = shname(ob.sta, inp.shnames_cat, inp.shnames_cat_igs)

''' output delays: '''
if ob.inp['delay_calc']:
    out_txt = 'delay.'+ob.source+'.'+date_string+'.'+stations_short[0]+\
                stations_short[1]+'.txt'
    with open(out_path + '/' + out_txt,'w') as f_txt:
        dz = ob.dude.delay
        for jj, _t in enumerate(ob.tstamps):
            ## stack line together for txt output:
            line = '{:4d} {:02d} {:02d} {:02d} {:02d} {:06.3f}  '.\
                    format(_t.year,_t.month,_t.day,\
                           _t.hour,_t.minute,_t.second+_t.microsecond*1e-6) + \
                    '{:23.16e} {:14.7e} {:14.7e} {:14.7e} {:14.7e}\n'.\
                    format(*dz[jj,:5])
            f_txt.write(line)
            
''' output doppler: '''
if ob.inp['doppler_calc']:
#    if inp.dop_model=='bary1way' or inp.dop_model=='geo1way':
    out_txt = 'doppler.'+ob.source+'.'+date_string+'.'+stations_short[0]+\
                '.txt'
    with open(out_path + '/' + out_txt,'w') as f_txt:
        dz = ob.dude.doppler
        for jj, _t in enumerate(ob.tstamps):
            ## stack line together for txt output:
            line = '{:4d} {:02d} {:02d} {:02d} {:02d} {:06.3f}  '.\
                    format(_t.year,_t.month,_t.day,\
                           _t.hour,_t.minute,_t.second+_t.microsecond*1e-6) + \
                    '{:23.16e} {:14.7e} {:14.7e} {:14.7e} {:14.7e}\n'.\
                    format(*dz[jj,:5])
            f_txt.write(line)