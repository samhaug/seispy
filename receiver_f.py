#!/usr/bin/env python

import obspy
import numpy as np
import filter

def make_dict(st_n, st_e, **kwargs):
    '''
    Make dictionary of trace objects so distances can be matched
    '''

    e_dict = {}
    n_dict = {}
    for tr in st_n:
        n_dict[tr.stats.network+tr.stats.station] = tr
    for tr in st_e:
        e_dict[tr.stats.network+tr.stats.station] = tr

    return n_dict, e_dict

def rotate_tr(tr_n, tr_e, **kwargs):
    '''
    Rotate trace from NE to RT coordinates
    '''

    tr_r = tr_n.copy()
    tr_t = tr_n.copy()

    r,t = obspy.signal.rotate.rotate_NE_RT(tr_n.data,
           tr_e.data,tr_n.stats.sac['baz'])

    tr_r.data = r
    tr_t.data = t

    tr_r.channel = 'BHR'
    tr_t.channel = 'BHT'

    return tr_r, tr_t

def rotate_st(st_n, st_e, **kwargs):
    '''
    Use sorted streams to rotate components
    '''

    n_dict, e_dict = make_dict(st_n,st_e)
    st_r = obspy.core.Stream()
    st_t = obspy.core.Stream()

    for keys in n_dict:
        if keys in e_dict:
            tr_r, tr_t = rotate_tr(n_dict[keys],e_dict[keys])
            st_r.append(tr_r)
            st_t.append(tr_t)
    return st_r, st_t

def express_rt():
    '''
    get r and t
    '''
    stn = obspy.read('*BHN*filtered')
    ste = obspy.read('*BHE*filtered')
    stz = obspy.read('*BHZ*filtered')
    stn = filter.gimp_filter(stn)
    ste = filter.gimp_filter(ste)
    stz = filter.gimp_filter(stz)
    if len(ste[0]) < len(stn[0]):
        for tr in stn:
            tr.data = tr.data[0:len(ste[0])]
    if len(stn[0]) < len(ste[0]):
        for tr in ste:
            tr.data = tr.data[0:len(stn[0])]
    str, stt = rotate_st(stn,ste)
    str.differentiate()
    stt.differentiate()
    stz.differentiate()
    str = filter.gimp_filter(str)
    stt = filter.gimp_filter(stt)
    stz = filter.gimp_filter(stz)
    for idx, tr in enumerate(str):
        tr.stats.channel = 'BHR'
        stt[idx].stats.channel = 'BHT'
        tr.stats.location = tr.stats.sac['gcarc']
        stt[idx].stats.location = stt[idx].stats.sac['gcarc']
    return str, stt, stz

