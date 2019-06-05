#!/usr/bin/env python

import numpy as np
from obspy.geodetics import gps2dist_azimuth

'''
This module is all about getting the metadata for each trace in a stream
ready for analysis
'''
def even_streams(sta,stb):
    a = []
    b = []
    for tr in sta:
        if tr.stats.network+tr.stats.station+tr.stats.location in a:
            sta.remove(tr)
        else:
            a.append(tr.stats.network+tr.stats.station+tr.stats.location)
    for tr in stb:
        if tr.stats.network+tr.stats.station+tr.stats.location in b:
            stb.remove(tr)
        else:
            b.append(tr.stats.network+tr.stats.station+tr.stats.location)
    c = set(a).intersection(set(b))
    for tr in sta:
        if tr.stats.network+tr.stats.station+tr.stats.location not in c:
            sta.remove(tr)
    for tr in stb:
        if tr.stats.network+tr.stats.station+tr.stats.location not in c:
            stb.remove(tr)

    return sta,stb

def sac_2_h5(st,name=False):
    for tr in st:
        try:
            tr.stats.evdp  = tr.stats.sac['evdp']
            tr.stats.evla  = tr.stats.sac['evla']
            tr.stats.evlo  = tr.stats.sac['evlo']
            tr.stats.stla  = tr.stats.sac['stla']
            tr.stats.stlo  = tr.stats.sac['stlo']
            tr.stats.o     = tr.stats.sac['o']
            tr.stats.gcarc = tr.stats.sac['gcarc']
            tr.stats.az    = tr.stats.sac['az']
            tr.stats.baz   = tr.stats.sac['baz']
        except KeyError:
            print('Metadata missing')
            continue
    st = set_az_gcarc(st)
    if name:
        st.write(name,format='H5')
    return st

def set_az_gcarc(st,**kwargs):
    f = kwargs.get('f',0.0033528106647474805)

    for tr in st:

        a = gps2dist_azimuth(tr.stats.evla,
                             tr.stats.evlo,
                             tr.stats.stla,
                             tr.stats.stlo,f=f)
        tr.stats.baz = a[-1]
        tr.stats.az = a[-2]
        tr.stats.gcarc = a[0]/111195.

    return st

def gemini_stations(st,name='gemini_STATIONS'):
    '''
    Make stations file for gemini run
    '''
    for tr in st:
        tr.stats.location = tr.stats.station+tr.stats.network
    st.sort(['location'])

    f = open(name,'w')
    with open(name,'w') as f:
        f.write('File last updated on \n')
        f.write(time.strftime("%m/%d/%Y_%H:%M:%S")+'\n')
        f.write('00  s_station s_location s_lat       s_long       s_elev   s_site\n')
        for idx,tr in enumerate(st):
            f.write('{:>5} {:11} {:12} {:12} {:12} {:8} {:30}\n'.format(
            str(st[0].stats.starttime.year)[-2::],
            tr.stats.station,
            tr.stats.network,
            tr.stats.sac['stla'],
            tr.stats.sac['stlo'],
            '0',
            'Ann Arbor, MI'))
    f.close()

def axisem_stations(st,name='STATIONS'):
    '''
    Make STATIONS ascii file for a 1D axisem run. remove redundant stations
    '''
    for tr in st:
        tr.stats.location = tr.stats.station+tr.stats.network
    st.sort(['location'])

    f = open(name,'w')
    for idx,tr in enumerate(st):
        #if tr.stats.station+tr.stats.network == \
        #    st[idx+1].stats.station+st[idx].stats.network:
        #    continue
        f.write('{}   {}   {}   {}   {}   {}   {}   {}\n'.format(
        tr.stats.station,
        tr.stats.network,
        round(tr.stats.sac['stla'],2),
        round(tr.stats.sac['stlo'],2),
        round(tr.stats.sac['evdp'],2),
        round(tr.stats.sac['gcarc'],2),
        round(tr.stats.sac['az'],2),
        round(tr.stats.sac['baz'],2)))
    f.close()

def axisem3d_stations(st,name='STATIONS'):
    '''
    Make STATIONS ascii file for a 1D axisem3d run. remove redundant stations
    '''
    for tr in st:
        tr.stats.name = tr.stats.network+tr.stats.network+tr.stats.location
    st.sort(['location'])

    f = open(name,'w')
    for idx,tr in enumerate(st):
        f.write('{}   {}   {}   {}   {}   {}\n'.format(
        tr.stats.station+'_'+tr.stats.location,
        tr.stats.network,
        tr.stats.stla,
        tr.stats.stlo,
        '0.0',
        '0.0'
        ))
    f.close()

def specfem_stations(st,name='STATIONS'):
    '''
    Make STATIONS ascii file for a 1D axisem run. remove redundant stations
    '''
    #for tr in st:
    #    tr.stats.location = tr.stats.station+tr.stats.network+tr.stats.location
    st.sort(['stla'])

    f = open(name,'w')
    for idx,tr in enumerate(st):
        #if tr.stats.station+tr.stats.network == \
        #    st[idx+1].stats.station+st[idx].stats.network:
        #    continue
        f.write('{}   {}   {}   {}   {}   {}\n'.format(
        tr.stats.station,
        tr.stats.network,
        round(tr.stats.stla,3),
        round(tr.stats.stlo,3),
        0.0,
        0.0))
    f.close()


