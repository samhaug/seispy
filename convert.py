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
        a = float(tr.stats.station.split('_')[0])
        b = float(tr.stats.station.split('_')[1])
        tr.stats.sac['gcarc'] = a+(b/1000.)
    return st

def axisem_stations(st):
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


