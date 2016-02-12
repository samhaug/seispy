#!/usr/bin/env python

import obspy
import numpy as np
import filter
import data
import seispy.data
from obspy.taup import TauPyModel
model = TauPyModel(model="premd")

def make_dict(st_1, st_2, **kwargs):
    '''
    Make dictionary of trace objects so distances can be matched
    '''

    dict_1 = {}
    dict_2 = {}
    for tr in st_1:
        dict_1[tr.stats.network+tr.stats.station] = tr
    for tr in st_2:
        dict_2[tr.stats.network+tr.stats.station] = tr

    return dict_1, dict_2

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

    tr_r.stats.location = tr_n.stats.sac['gcarc']
    tr_t.stats.location = tr_n.stats.sac['gcarc']

    tr_r.stats.azimuth = tr_n.stats.sac['az']
    tr_t.stats.azimuth = tr_t.stats.sac['az']

    tr_r.stats.channel = 'BHR'
    tr_t.stats.channel = 'BHT'

    return tr_r, tr_t

def rotate_st(stn, ste, **kwargs):
    '''
    Use sorted streams to rotate components
    '''
    str = obspy.core.Stream()
    stt = obspy.core.Stream()

    for idx,tr in enumerate(stn):
        trr, trt = rotate_tr(tr,ste[idx])
        str.append(trr)
        stt.append(trt)
    return str, stt

def express_zne():
    '''
    get z,n,e components
    '''
    stz = obspy.read('*BHZ*filtered')
    stz = filter.gimp_filter(stz)
    stn = obspy.read('*BHN*filtered')
    stn = filter.gimp_filter(stn)
    ste = obspy.read('*BHE*filtered')
    ste = filter.gimp_filter(ste)
    stz.normalize()
    stn.normalize()
    ste.normalize()
    stz = data.align_on_phase(stz)
    stn = data.align_on_phase(stn)
    ste = data.align_on_phase(ste)
    return stz, stn, ste

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

    stn_name_list = []
    ste_name_list = []
    stz_name_list = []

    for tr in stn:
        tr.stats.full_name = tr.stats.network+'.'+tr.stats.station
        stn_name_list.append(tr.stats.network+'.'+tr.stats.station)
    for tr in ste:
        tr.stats.full_name = tr.stats.network+'.'+tr.stats.station
        ste_name_list.append(tr.stats.network+'.'+tr.stats.station)
    for tr in stz:
        tr.stats.full_name = tr.stats.network+'.'+tr.stats.station
        stz_name_list.append(tr.stats.network+'.'+tr.stats.station)

    common_name = set(stn_name_list) & set(ste_name_list) & set(stz_name_list)

    for tr in stn:
        if tr.stats.full_name not in common_name:
            stn.remove(tr)
    for tr in ste:
        if tr.stats.full_name not in common_name:
            ste.remove(tr)
    for tr in stz:
        if tr.stats.full_name not in common_name:
            stz.remove(tr)

    stz.sort(['full_name'])
    ste.sort(['full_name'])
    stn.sort(['full_name'])

    str, stt = rotate_st(stn,ste)

    return str, stt, stz

def rz_2_lq(str,stz):
    '''
    convert radial/Z component to max P/ max S component
    find maximum P wave amplitude to determine incidence angle for every
    trace
    '''

    stl = stz.copy()
    stq = str.copy()

    for idx,tr in enumerate(str):

        R = seispy.data.phase_window(tr,['P'],(-10,10)).data.min()
        Z = seispy.data.phase_window(stz[idx],['P'],(-10,10)).data.min()
        if R < 0 and Z < 0:
            R*=-1
            Z*=-1
            deg = np.arctan2(R,Z)
        elif R < 0 and Z > 0:
            deg = np.arctan2(R,Z)
        else:
            deg = -1*np.arctan2(R,Z)

        stq[idx].data = np.cos(deg)*tr.data-np.sin(deg)*stz[idx].data
        stl[idx].data = np.sin(deg)*tr.data+np.cos(deg)*stz[idx].data
        stq[idx].stats.channel = 'BHS'
        stl[idx].stats.channel = 'BHP'

    return stl, stq














