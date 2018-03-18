#!/usr/bin/env python

import obspy
from obspy.signal.rotate import rotate_NE_RT
import numpy as np
import filter
import data
import seispy.data
from obspy.taup import TauPyModel
model = TauPyModel(model="prem")

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

def obspy_rotate_ne_rt(stn,ste):
    str = stn.copy()
    stt = stn.copy()
    for idx,tr in enumerate(stn):
        r,t = obspy.signal.rotate.rotate_ne_rt(stn[idx].data,
                 ste[idx].data,stn[idx].stats.sac['baz'])
        str[idx].data = r
        str[idx].stats.channel = 'R'
        stt[idx].data = t
        stt[idx].stats.channel = 'T'
    return str,stt

def rotate_ne_rt(stn,ste):
    str = stn.copy()
    stt = stn.copy()
    for idx,tr in enumerate(stn):
        n = stn[idx].data
        e = ste[idx].data
        ne = np.vstack((n,e))
        rot = np.array([[np.cos(np.radians(tr.stats.sac['az'])),np.sin(np.radians(tr.stats.sac['az']))],
                        [-1*np.sin(np.radians(tr.stats.sac['az'])),np.cos(np.radians(tr.stats.sac['az']))]])
        rt = np.dot(rot,ne)
        str[idx].data = rt[0,:]
        stt[idx].data = rt[1,:]
    return str,stt

def rotate_tr(tr_n, tr_e, **kwargs):
    '''
    Rotate trace from NE to RT coordinates
    '''
    tr_r = tr_n.copy()
    tr_t = tr_n.copy()

    r,t = rotate_NE_RT(tr_n.data,tr_e.data,tr_n.stats.sac['baz'])

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

def rotate_LQT(stn, ste, stz, **kwargs):
    stn.sort(['station'])
    ste.sort(['station'])
    stz.sort(['station'])
    str, stt = rotate_st(stn,ste)
    stl, stq = rz_2_lq(str,stz)
    return stl, stq, stt

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

def express_all():
    '''
    get r and t z n e
    '''
    stn = obspy.read('*BHN*filtered')
    ste = obspy.read('*BHE*filtered')
    stz = obspy.read('*BHZ*filtered')

    stn_name_list = []
    ste_name_list = []
    stz_name_list = []

    stn = filter.gimp_filter(stn)
    ste = filter.gimp_filter(ste)
    stz = filter.gimp_filter(stz)
    for tr in stn:
        tr.stats.full_name = tr.stats.network+'.'+tr.stats.station+'.'
                             #str(tr.stats.calib)
        if tr.stats.full_name in stn_name_list:
            stn.remove(tr)
        stn_name_list.append(tr.stats.full_name)
    for tr in ste:
        tr.stats.full_name = tr.stats.network+'.'+tr.stats.station+'.'
                             #str(tr.stats.calib)
        if tr.stats.full_name in ste_name_list:
            ste.remove(tr)
        ste_name_list.append(tr.stats.full_name)
    for tr in stz:
        tr.stats.full_name = tr.stats.network+'.'+tr.stats.station+'.'
                             #str(tr.stats.calib)
        if tr.stats.full_name in stz_name_list:
            stz.remove(tr)
        stz_name_list.append(tr.stats.full_name)

    #common_name = set(stn_name_list)|set(ste_name_list)|set(stz_name_list)
    #full_name = set(stz_name_list+ste_name_list+stn_name_list)
    #common_list = min(len(set(stn_name_list)),len(set(stz_name_list)),len(set(ste_name_list)))
    #return common_name,full_name
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


    print len(stn),len(ste),len(stz)
    if (len(stn[0].data) == len(ste[0].data) == len(stz[0].data)) == False:
        new_len = min(len(stn[0].data),len(ste[0].data),len(stz[0].data))
        for idx,tr in enumerate(stn):
            tr.data = tr.data[0:new_len-1]
            ste[idx].data = ste[idx].data[0:new_len-1]
            stz[idx].data = stz[idx].data[0:new_len-1]

    strad, stt = rotate_st(stn,ste)

    #return common_name
    return  strad, stt, stz, ste, stn

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

def find_incidence_angle(trr,trz,**kwargs):
    '''
    Return incidence angle in radians
    '''

    pol = kwargs.get('polarity','min')

    R_min = seispy.data.phase_window(trr,['P'],(-10,10)).data.min()
    Z_min = seispy.data.phase_window(trz,['P'],(-10,10)).data.min()
    R_max = seispy.data.phase_window(trr,['P'],(-10,10)).data.max()
    Z_max = seispy.data.phase_window(trz,['P'],(-10,10)).data.max()

    if pol == 'min':
        R = R_min
        Z = Z_min
    elif pol == 'max':
        R = R_max
        Z = Z_max
    if R < 0 and Z < 0:
        R*=-1
        Z*=-1
        deg = np.arctan2(R,Z)
    elif R < 0 and Z > 0:
        deg = np.arctan2(R,Z)
    else:
        deg = -1*np.arctan2(R,Z)
    return deg

def express_zne():
    '''
    read zne components
    '''
    stz = obspy.read('*BHZ*filtered')
    ste = obspy.read('*BHE*filtered')
    stn = obspy.read('*BHN*filtered')
    stz = seispy.filter.gimp_filter(stz)
    ste = seispy.filter.gimp_filter(ste)
    stn = seispy.filter.gimp_filter(stn)
    return stz, ste, stn

def M3D(i,ba):
    '''
    Return matrix to rotate from ZEN to LQT
    i and baz in radians
    '''
    m3d = np.zeros((3,3))
    m3d[0,0] = np.cos(i)
    m3d[1,1] = np.cos(i)*np.sin(ba)
    m3d[2,2] = np.sin(ba)
    m3d[1,0] = np.sin(i)
    m3d[2,0] = 0
    m3d[0,1] = -np.sin(i)*np.sin(ba)
    m3d[0,2] = -np.sin(i)*np.cos(ba)
    m3d[1,2] =  np.cos(i)*np.cos(ba)
    m3d[2,1] =  -np.cos(ba)

    return m3d

def M2D(ba):
    '''
    make matrix to rotate L and Q components for beamforming
    '''
    ba = np.radians(ba)
    m2d = np.zeros((2,2))
    m2d[0,0] = np.cos(ba)
    m2d[1,1] = np.cos(ba)
    m2d[1,0] = np.sin(ba)
    m2d[0,1] = -np.sin(ba)

    return m2d

def beamform_LQ(trl,trq,baz):
    '''
    rotate azimuth to beamform traces a la Niu Kawakatsu 1997
    '''
    if trl.stats.sac['baz'] != trq.stats.sac['baz']:
        print "Traces do not match stations"
    m2d = M2D(baz)

    rotated = np.dot(m2d,np.vstack((trq.data,trl.data)))
    q_data = rotated[0,:]
    l_data = rotated[1,:]
    return l_data,q_data

def express_LQT():
    '''
    Execute in directory to get stl, stq, stt
    '''

    def main():
        str, stt, stz, ste, stn = express_all()
        stl = obspy.core.Stream()
        stq = obspy.core.Stream()
        for ii in range(len(stz)):
            trl, trq, trt = apply_rotation(stz[ii],ste[ii],stn[ii],
                            str[ii])
            stl.append(trl)
            stq.append(trq)

        return stl, stq, stt


    def apply_rotation(trz,tre,trn,trr):
        trl = trz.copy()
        trl.stats.channel = 'BHL'
        trq = tre.copy()
        trq.stats.channel = 'BHQ'
        trt = trn.copy()
        trt.stats.channel = 'BHT'

        deg = find_incidence_angle(trr,trz)
        baz = np.radians(trz.stats.sac['baz'])

        m3d = M3D(deg,baz)

        zen = np.vstack((([trz.data]),
                         ([tre.data]),
                         ([trn.data])))
        lqt = np.dot(m3d,zen)
        trl.data = lqt[0,:]
        trq.data = lqt[1,:]
        trt.data = lqt[2,:]

        return trl,trq,trt

    return main()




























