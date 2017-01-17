#!/usr/bin/env python

import numpy as np
import obspy
from geopy.distance import great_circle
from obspy.taup import TauPyModel
model = TauPyModel(model="prem")



'''
This module is all about getting the metadata for each trace in a stream
ready for analysis
'''

def mineos_convert(st):
    st = set_baz(st)
    st = SOD_evdp(st)
    st = set_gcarc(st)
    return st

def set_baz(st,**kwargs):

    f = kwargs.get('f',0.0033528106647474805)

    for tr in st:
        baz = obspy.geodetics.gps2dist_azimuth(tr.stats.sac['evla'],
                                               tr.stats.sac['evlo'],
                                               tr.stats.sac['stla'],
                                               tr.stats.sac['stlo'],f=f)[-1]
        az = obspy.geodetics.gps2dist_azimuth(tr.stats.sac['evla'],
                                               tr.stats.sac['evlo'],
                                               tr.stats.sac['stla'],
                                               tr.stats.sac['stlo'],f=f)[-2]
        tr.stats.sac['baz'] = baz
        tr.stats.sac['az'] = az
    return st

def master_set(st):

    for tr in st:
        tr.stats.location = tr.stats.sac['gcarc']
        tr.stats.sortname = tr.stats.network+tr.stats.station+str(tr.stats.location)
    st.sort(['sortname'])
    return st

def set_origin_time(st,**kwargs):
    '''
    set sac['o'] time for events retrieved from SOD.
    '''
    phase = kwargs.get('phase',['P'])
    event_depth = st[0].stats.sac['evdp']
    for tr in st:
        arrivals = model.get_travel_times(source_depth_in_km=event_depth,
                  distance_in_degree=tr.stats.sac['gcarc'],phase_list=phase)
        time = arrivals[0].time
        tr.stats.sac['o'] = -1*time
    return st

def SOD_evdp(st):
    '''
    divide all event depths by 1000
    '''
    for tr in st:
        tr.stats.sac['evdp'] *= 0.001
    return st

def equalize_start_end(st):
    '''
    make all traces in the stream share the same start and end time. Do this
    by cutting off data from longer streams
    '''
    out_st = st.copy()
    starttime_list = []
    endtime_list = []
    for tr in out_st:
        starttime_list.append(tr.stats.starttime)
        endtime_list.append(tr.stats.endtime)
    start = min(starttime_list)
    end = max(endtime_list)
    print start, end
    for tr in out_st:
        t = tr.stats.starttime
        e = tr.stats.endtime
        samp = tr.stats.sampling_rate
        begin_pad = int(np.abs(t-start)*samp)
        end_pad = int(np.abs(e-end)*samp)
        tr.data = np.hstack((np.zeros(begin_pad),tr.data,np.zeros(end_pad)))
    return out_st

def set_gcarc(st):
    '''
    Set gcarc for each trace in stream given station and event lat lon exist
    '''
    for tr in st:
        source = (tr.stats.sac['evla'],tr.stats.sac['evlo'])
        stat = (tr.stats.sac['stla'],tr.stats.sac['stlo'])
        tr.stats.sac['gcarc'] = np.abs(great_circle(source,stat).km/111.195)
    return st

def xh2sac(st):
    '''
    add sac dictionary to stream object produced for xh files
    '''

    for tr in st:
        tr.stats.sac = {}
        tr.stats.sac['evla'] = tr.stats.xh['source_latitude']
        tr.stats.sac['evlo'] = tr.stats.xh['source_longitude']
        tr.stats.sac['stla'] = tr.stats.xh['receiver_latitude']
        tr.stats.sac['stlo'] = tr.stats.xh['receiver_longitude']
        tr.stats.sac['evdp'] = tr.stats.xh['source_depth_in_km']
        tr.stats.sac['o'] = 0.
        tr.stats.sac['az'] = tr.stats.xh['sensor_azimuth']
        tr.stats.sac['baz'] = tr.stats.xh['sensor_azimuth']-180
        source = (tr.stats.xh['source_latitude'],tr.stats.xh['source_longitude'])
        stat = (tr.stats.xh['receiver_latitude'],tr.stats.xh['receiver_longitude'])
        tr.stats.sac['gcarc'] = np.abs(great_circle(source,stat).km/111.195)
    return st

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


