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
    if len(ste[0]) < len(stn[0]):
        for tr in stn:
            tr.data = tr.data[0:len(ste[0])]
    if len(stn[0]) < len(ste[0]):
        for tr in ste:
            tr.data = tr.data[0:len(stn[0])]
    str, stt = rotate_st(stn,ste)

    str = filter.gimp_filter(str)
    stt = filter.gimp_filter(stt)
    stz = filter.gimp_filter(stz)

    for idx, tr in enumerate(str):
        tr.stats.channel = 'BHR'
        stt[idx].stats.channel = 'BHT'
        tr.stats.location = tr.stats.sac['gcarc']
        stt[idx].stats.location = stt[idx].stats.sac['gcarc']
        tr.stats.azimuth = tr.stats.sac['az']
        stt[idx].stats.azimuth = stt[idx].stats.sac['az']
    z_dict,t_dict =make_dict(stz,stt)

    new_stz = obspy.core.Stream()
    for ii in t_dict.keys():
        new_stz += z_dict[ii]
    str.sort(['station'])
    stt.sort(['station'])
    new_stz.sort(['station'])

    return str, stt, new_stz

def rz_2_ps(str,stz):
    '''
    convert radial/Z component to max P/ max S component
    find maximum P wave amplitude to determine incidence angle for every
    trace
    '''

    stp = stz.copy()
    sts = str.copy()

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
        print R, Z

        sts[idx].data = np.cos(deg)*tr.data-np.sin(deg)*stz[idx].data
        stp[idx].data = np.sin(deg)*tr.data+np.cos(deg)*stz[idx].data
        sts[idx].stats.channel = 'BHS'
        stp[idx].stats.channel = 'BHP'

    return sts, stp














