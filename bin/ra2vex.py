#!/usr/bin/env python2.7
from pypride.vex import Vex 
import argparse
import os

def radification(vex_file, out_path):
    
    if out_path is None:
        out_path = os.path.dirname(vex_file)
    
    # parse vex-file:
    vex = Vex(vex_file)
    
    # find out the obs band:
    with open(vex_file,'r') as vexInFile:
        # read the input vex-file into one string
        vexOneString = vexInFile.read()
    with open(vex_file,'r') as vexInFile:
        # read the input vex-file into a list of lines:
        vexIn = vexInFile.readlines()
    
    k_mode = 'ra1cm2'
    c_mode = 'ra6cm2'
    l_mode = 'ra18cm2'
    
    k_band = k_mode in vexOneString
    c_band = c_mode in vexOneString
    l_band = l_mode in vexOneString #or 'L256' in vexOneString
    
    # if none of these is the case, set mode name
    if not (k_band and c_band and l_band):
        pos = [i for i,v in enumerate(vexIn) if '$MODE' in v][0]
        mode = vexIn[pos+2].split()[1]
        if 'L' in mode:
            l_band = True
            l_mode = mode[:-1]
        elif 'C' in mode:
            c_band = True
            c_mode = mode[:-1]
        elif 'K' in mode:
            k_band = True
            k_mode = mode[:-1]
        
    das_ra = '2NONE<' in vexOneString
    
    #print k_band, c_band, l_band
    
    #############################################################################
    ## add obsfreq-unrelated stuff to the vex-file:
    # station definition:
    pos = vexIn.index('$STATION;\n')
    
    vexIn.insert(pos+1,'enddef;\n')
    vexIn.insert(pos+1,'     ref $CLOCK = Ra;\n')
    vexIn.insert(pos+1,'     ref $DAS = 2NONE<;\n')
    vexIn.insert(pos+1,'     ref $ANTENNA = RA;\n')
    vexIn.insert(pos+1,'     ref $SITE = RA;\n')
    vexIn.insert(pos+1,'def Ra;\n')
    vexIn.insert(pos+1,'*\n')
    
    # site definition:
    pos = vexIn.index('$SITE;\n')
    
    vexIn.insert(pos+1,'enddef;\n')
    vexIn.insert(pos+1,'     site_position_epoch =       0;\n')
    vexIn.insert(pos+1,'     site_velocity =  0.000000   m/yr:  0.000000   m/yr:  0.000000  m/yr;\n')
    vexIn.insert(pos+1,'     site_position = 9999999.00000 m: 9999999.00000 m: 9999999.00000 m;\n')
    vexIn.insert(pos+1,'     site_ID = Ra;\n')
    vexIn.insert(pos+1,'     site_name = RA;\n')
    vexIn.insert(pos+1,'     site_type = earth_orbit;\n')
    vexIn.insert(pos+1,'def RA;\n')
    vexIn.insert(pos+1,'*\n')
    
    # antenna definition:
    pos = vexIn.index('$ANTENNA;\n')
    
    vexIn.insert(pos+1,'enddef;\n')
    vexIn.insert(pos+1,'     axis_offset =    0.00000 m;\n')
    vexIn.insert(pos+1,'*     antenna_motion = az :  15.0 deg/min : 20 sec;  * 1000.000 deg/sec/sec\n')
    vexIn.insert(pos+1,'*     antenna_motion = el :  15.0 deg/min : 20 sec;  * 1000.000 deg/sec/sec\n')
    vexIn.insert(pos+1,'     axis_type = sp : ace;\n')
    vexIn.insert(pos+1,'def RA;\n')
    vexIn.insert(pos+1,'*\n')
    
    # das definition, if not yet in the file:
    if not das_ra:
        pos = vexIn.index('$DAS;\n')
    
        vexIn.insert(pos+1,'enddef;\n')
        vexIn.insert(pos+1,'     tape_motion = adaptive : 0 min: 0 min: 10 sec;\n')
        vexIn.insert(pos+1,'     headstack = 2 :            : 1 ;\n')
        vexIn.insert(pos+1,'     headstack = 1 :            : 0 ;\n')
        vexIn.insert(pos+1,'     number_drives = 2;\n')
        vexIn.insert(pos+1,'     electronics_rack_type = none;\n')
        vexIn.insert(pos+1,'     record_transport_type = Mark5B;\n')
        vexIn.insert(pos+1,'def 2NONE<;\n')
        vexIn.insert(pos+1,'*\n')
    
    # bbc definition:
    pos = vexIn.index('$BBC;\n')
    
    vexIn.insert(pos+1,'*\n')
    vexIn.insert(pos+1,'enddef;\n')
    vexIn.insert(pos+1,'     BBC_assign = &BBC02 :  2 : &IF_B;\n')
    vexIn.insert(pos+1,'     BBC_assign = &BBC01 :  1 : &IF_A;\n')
    vexIn.insert(pos+1,'* mode =  *    stations =Ra\n')
    vexIn.insert(pos+1,'def 4BBCs#RA;\n')
    
    # bitstreams definition:
    try:
        pos = vexIn.index('$BITSTREAMS;\n')
    
        vexIn.insert(pos+1,'enddef;\n')
        vexIn.insert(pos+1,'     stream_def = &CH04 : sign : 3:  3;\n')
        vexIn.insert(pos+1,'     stream_def = &CH03 : sign : 2:  2;\n')
        vexIn.insert(pos+1,'     stream_def = &CH02 : sign : 1:  1;\n')
        vexIn.insert(pos+1,'     stream_def = &CH01 : sign : 0:  0;\n')
        vexIn.insert(pos+1,'* stations = Ra\n')
        vexIn.insert(pos+1,'def Mk5B.Ra;\n')
        vexIn.insert(pos+1,'*\n')
    
    except ValueError:
        pos = vexIn.index('$ROLL;\n')
    
        vexIn.insert(pos,'*------------------------------------------------------------------------------\n')
        vexIn.insert(pos,'enddef;\n')
        vexIn.insert(pos,'     stream_def = &CH04 : sign : 3:  3;\n')
        vexIn.insert(pos,'     stream_def = &CH03 : sign : 2:  2;\n')
        vexIn.insert(pos,'     stream_def = &CH02 : sign : 1:  1;\n')
        vexIn.insert(pos,'     stream_def = &CH01 : sign : 0:  0;\n')
        vexIn.insert(pos,'* stations = Ra\n')
        vexIn.insert(pos,'def Mk5B.Ra;\n')
        vexIn.insert(pos,'*\n')
        vexIn.insert(pos,'$BITSTREAMS;\n')
    
    
    ## add obsfreq-related stuff to the vex-file:
    # K-band:
    if k_band:
        # mode definition:
        mode_string = ''.join(('def ', k_mode, ';\n'))
        pos = vexIn.index(mode_string)
        while vexIn[pos] != 'enddef;\n':
            pos = pos + 1
        vexIn.insert(pos,'     ref $BBC = 4BBCs#RA:Ra;\n')
        vexIn.insert(pos,'     ref $IF = LO@22232MHzDPolTone/1:Ra;\n')
        vexIn.insert(pos,'     ref $FREQ = 22228.00MHz4x16MHz#RA:Ra;\n')
        vexIn.insert(pos,'     ref $BITSTREAMS = Mk5B.Ra:Ra;\n')
    
        # freq definition:
        pos = vexIn.index('$FREQ;\n')
        vexIn.insert(pos+1,'enddef;\n')
        vexIn.insert(pos+1,'     chan_def = : 22228.00 MHz : L :  16.00 MHz : &CH04 : &BBC02 : &L_Cal; *Lcp\n')
        vexIn.insert(pos+1,'     chan_def = : 22228.00 MHz : U :  16.00 MHz : &CH03 : &BBC02 : &L_Cal; *Lcp\n')
        vexIn.insert(pos+1,'     chan_def = : 22228.00 MHz : L :  16.00 MHz : &CH02 : &BBC01 : &L_Cal; *Rcp\n')
        vexIn.insert(pos+1,'     chan_def = : 22228.00 MHz : U :  16.00 MHz : &CH01 : &BBC01 : &L_Cal; *Rcp\n')
        vexIn.insert(pos+1,'     sample_rate =  32.000 Ms/sec;  * (1bits/sample)\n')
        vexIn.insert(pos+1,'* mode =  k    stations =Ra\n')
        vexIn.insert(pos+1,'def 22228.00MHz4x16MHz#RA;\n')
        vexIn.insert(pos+1,'*\n')
    
        # if definition:
        pos = vexIn.index('$IF;\n')
        vexIn.insert(pos+1,'enddef;\n')
        vexIn.insert(pos+1,'     if_def = &IF_B : B : L : 22232.0 MHz : L : 1 MHz ; *\n')
        vexIn.insert(pos+1,'     if_def = &IF_A : A : R : 22232.0 MHz : L : 1 MHz ; *\n')
        vexIn.insert(pos+1,'* mode =  k    stations =Ra\n')
        vexIn.insert(pos+1,'def LO@22232MHzDPolTone/1;\n')
        vexIn.insert(pos+1,'*\n')
    
    # C-band:
    if c_band:
        # mode definition:
        mode_string = ''.join(('def ', c_mode, ';\n'))
        pos = vexIn.index(mode_string)
        while vexIn[pos] != 'enddef;\n':
            pos = pos + 1
        vexIn.insert(pos,'     ref $BBC = 4BBCs#RA:Ra;\n')
        vexIn.insert(pos,'     ref $IF = LO@4320MHzDPolTone/1:Ra;\n')
        vexIn.insert(pos,'     ref $FREQ = 4828.00MHz4x16MHz#RA:Ra;\n')
        vexIn.insert(pos,'     ref $BITSTREAMS = Mk5B.Ra:Ra;\n')
    
        # freq definition:
        pos = vexIn.index('$FREQ;\n')
        vexIn.insert(pos+1,'enddef;\n')
        vexIn.insert(pos+1,'     chan_def = : 4828.00 MHz : L :  16.00 MHz : &CH04 : &BBC01 : &L_Cal; *Rcp\n')
        vexIn.insert(pos+1,'     chan_def = : 4828.00 MHz : U :  16.00 MHz : &CH03 : &BBC01 : &L_Cal; *Rcp\n')
        vexIn.insert(pos+1,'     chan_def = : 4828.00 MHz : L :  16.00 MHz : &CH02 : &BBC02 : &L_Cal; *Lcp\n')
        vexIn.insert(pos+1,'     chan_def = : 4828.00 MHz : U :  16.00 MHz : &CH01 : &BBC02 : &L_Cal; *Lcp\n')
        vexIn.insert(pos+1,'     sample_rate =  32.000 Ms/sec;  * (1bits/sample)\n')
        vexIn.insert(pos+1,'* mode =  c    stations =Ra\n')
        vexIn.insert(pos+1,'def 4828.00MHz4x16MHz#RA;\n')
        vexIn.insert(pos+1,'*\n')
    
        # if definition:
        pos = vexIn.index('$IF;\n')
        vexIn.insert(pos+1,'enddef;\n')
        vexIn.insert(pos+1,'     if_def = &IF_B : B : L :  4320.0 MHz : U : 1 MHz ; *\n')
        vexIn.insert(pos+1,'     if_def = &IF_A : A : R :  4320.0 MHz : U : 1 MHz ; *\n')
        vexIn.insert(pos+1,'* mode =  c    stations =Ra\n')
        vexIn.insert(pos+1,'def LO@4320MHzDPolTone/1;\n')
        vexIn.insert(pos+1,'*\n')
    
    # L-band:
    if l_band:
        # mode definition:
        mode_string = ''.join(('def ', l_mode, ';\n'))
        pos = vexIn.index(mode_string)
        while vexIn[pos] != 'enddef;\n':
            pos = pos + 1
        vexIn.insert(pos,'     ref $BBC = 4BBCs#RA:Ra;\n')
        vexIn.insert(pos,'     ref $IF = LO@1152MHzDPolTone/1:Ra;\n')
        vexIn.insert(pos,'     ref $FREQ = 1660.00MHz4x16MHz#RA:Ra;\n')
        vexIn.insert(pos,'     ref $BITSTREAMS = Mk5B.Ra:Ra;\n')
    
        # freq definition:
        pos = vexIn.index('$FREQ;\n')
        vexIn.insert(pos+1,'enddef;\n')
        vexIn.insert(pos+1,'     chan_def = :  1660.00 MHz : L :  16.00 MHz : &CH04 : &BBC01 : &L_Cal; *Rcp\n')
        vexIn.insert(pos+1,'     chan_def = :  1660.00 MHz : U :  16.00 MHz : &CH03 : &BBC01 : &L_Cal; *Rcp\n')
        vexIn.insert(pos+1,'     chan_def = :  1660.00 MHz : L :  16.00 MHz : &CH02 : &BBC02 : &L_Cal; *Lcp\n')
        vexIn.insert(pos+1,'     chan_def = :  1660.00 MHz : U :  16.00 MHz : &CH01 : &BBC02 : &L_Cal; *Lcp\n')
        vexIn.insert(pos+1,'     sample_rate =  32.000 Ms/sec;  * (1bits/sample)\n')
        vexIn.insert(pos+1,'* mode =  l    stations =Ra\n')
        vexIn.insert(pos+1,'def 1660.00MHz4x16MHz#RA;\n')
        vexIn.insert(pos+1,'*\n')
    
        # if definition:
        pos = vexIn.index('$IF;\n')
        vexIn.insert(pos+1,'enddef;\n')
        vexIn.insert(pos+1,'     if_def = &IF_B : B : L :  1152.0 MHz : U : 1 MHz ; *\n')
        vexIn.insert(pos+1,'     if_def = &IF_A : A : R :  1152.0 MHz : U : 1 MHz ; *\n')
        vexIn.insert(pos+1,'* mode =  l    stations =Ra\n')
        vexIn.insert(pos+1,'def LO@1152MHzDPolTone/1;\n')
        vexIn.insert(pos+1,'*\n')
    
    ## add ra scans:
    indices = [i for i, x in enumerate(vexIn) if x == "endscan;\n"]
    # start from the end to prevent spoiling the line numbers where to add ra
    indices.reverse()
    # for each scan create a line to add and add it:
    sched_inv = [s for s in vex['SCHED']][::-1]
    for si, s in enumerate(sched_inv):
        if len(vex['SCHED'][s].getall('station'))==1:
            # tracking or calibration scan, skip:
            continue
        pos = indices[si]
        ra_scan = vexIn[pos-1]
        ra_scan_col = ra_scan.index(':')
        ra_scan_eqs = ra_scan.index('=')
        ra_scan = ra_scan[0:ra_scan_eqs+1] + 'Ra' + ra_scan[ra_scan_col:]
        vexIn.insert(pos, ra_scan)

    
    # output:
    vex_base = os.path.basename(vex_file)
    indices = [i for i, x in enumerate(vex_base) if x == "."]
    vexOut = vex_base[0:indices[-1]] + '_ra' + vex_base[indices[-1]:]
    with open(os.path.join(out_path, vexOut),'w') as vexOutFile:
        print 'outputting to {:s}'.format(os.path.join(out_path, vexOut))
        for line in vexIn:
            vexOutFile.write(line)


if __name__=='__main__':
    # create parser
    parser = argparse.ArgumentParser(prog='ra2vex.py', \
                formatter_class=argparse.RawDescriptionHelpFormatter,\
                description='Add RadioAstron-specific entries to a vex-file.\n'+\
                            'Note that it\'s kind of fake, so check it!')
    # optional argument
    parser.add_argument('-o', '--outpath', type=str, default=None,
                        help='output path, defaults to the path to vex-file')
    # positional argument
    parser.add_argument('vexfile', type=str, help='input vex-file')
    args = parser.parse_args()
     
    # run radification
    radification(args.vexfile, args.outpath)
