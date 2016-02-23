#!/usr/bin/env python

import numpy as np
import obspy


def xh2sac(st):
    '''
    add sac dictionary to stream object produced for xh files
    '''

    for tr in st:
        tr.stats.sac = {}
        tr.stats.sac['evla'] = 0.
        tr.stats.sac['evlo'] = 0.
        tr.stats.sac['stla'] = tr.stats.xh['receiver_latitude']
        tr.stats.sac['stlo'] = 0.
        tr.stats.sac['evdp'] = tr.stats.xh['source_depth_in_km']
        tr.stats.sac['o'] = 0.
        tr.stats.sac['baz'] = 0.
        tr.stats.sac['gcarc'] = tr.stats.xh['receiver_latitude']
    return st

def axisem_stations(st):
    '''
    Make STATIONS ascii file for a 1D axisem run. remove redundant stations
    '''
    for tr in st:
        tr.stats.location = tr.stats.station+tr.stats.network
    st.sort(['location'])

    f = open('STATIONS','w')
    for idx,tr in enumerate(st[:-1]):
        if tr.stats.station+tr.stats.network == \
            st[idx+1].stats.station+st[idx].stats.network:
            continue
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

def axisem_stations_2D(st):
    '''
    make STATIONS ascii file for use with axisem. Rotate the source receiver
    geometry so that the source is at the north pole
    '''
    def cart_coord(lat,lon):
        if lat < 0:
            lat = abs(lat)+90.
        if lat > 0:
            lat = 90.-lat
        if lat == 0:
            lat = 90.
        if lon < 0:
            lon += 360.
        return lat, lon
        lat = np.radians(lat)
        lon = np.radians(lon)
        cart_pole = [np.cos(lon)*np.sin(lat),np.sin(lon)*np.sin(lat),np.cos(lat)]
        #return cart_pole

    def get_euler_pole(st):
        north_pole = [0,0,1]
        lat = st[0].stats.sac['evla']
        lon = st[0].stats.sac['evlo']
        event_pole = cart_coord(lat,lon)
        euler_pole = np.cross(event_pole,north_pole)
        return euler_pole

    print cart_coord(10,0)


