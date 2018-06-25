#!/home/samhaug/anaconda2/bin/python

import scipy
import obspy
import numpy as np
from obspy.taup import TauPyModel
from matplotlib import pyplot as plt
from scipy.signal import argrelextrema
model = TauPyModel(model="prem")
import itertools
import seispy.convert
from scipy.signal import tukey
from scipy.signal import cosine


def water_level(a,b,alpha=0.1,plot=False):
    #ross's code
    t = np.linspace(0,len(a)/10.,len(a))

    #Convert to frequency domain-------------------------------
    a_omega = np.fft.fft(a)
    b_omega = np.fft.fft(b)

    #Perform division------------------------------------------
    F_omega = b_omega*b_omega.conjugate()
    Phi_ss  = np.maximum(F_omega,alpha*(np.amax(F_omega)))
    try:
        H_omega = ((a_omega*b_omega.conjugate())/Phi_ss)
    except RuntimeWarning:
        return np.zeros(len(a))

    #Convert back to time domain-------------------------------
    #rf = np.zeros(len(H_omega))
    rf = np.fft.ifft(H_omega)

    #Plots-----------------------------------------------------
    if plot==True:
        fig,ax = plt.subplots(2,1)
        ax[0].plot(t,a)
        ax[0].plot(t,b)
        ax[1].plot(t,rf)
        plt.show()

    #roll = int(10*50-np.argmax(rf))
    #rf = np.roll(rf,roll)
    return t,np.real(rf)

def mask_depth_phase(st,**kwargs):
    offset = kwargs.get('offset',0)
    phase_list = ['pPdiff','sPdiff','pPKP','sPKP','pPKIKP','sPKIKP',
                  'pPKiKP','sPKiKP','pPP','sPP','pSKP','sSKP','pP','sP']

    for tr in st:
        arr = model.get_travel_times(tr.stats.sac['evdp'],
                                     tr.stats.sac['gcarc'],phase_list)
        for ii in arr:
            mask = np.zeros(tr.stats.npts)
            t = int((ii.time+tr.stats.sac['o'])*tr.stats.sampling_rate)
            mask[t-offset:t-offset+int(30*tr.stats.sampling_rate)] = tukey(int(30*tr.stats.sampling_rate))
            tr.data*= 1.-mask
    return st

def rotate_phase(stz,stn,ste,phase,**kwargs):
    '''
    Rotate a three component trace in zne to coordinate system for specific
    phase.
    do not normalize or process the array in any way.
    '''
    model = kwargs.get('model','prem50')
    model = TauPyModel(model=model)

    def make_R(i,baz):
        R = np.matrix([[np.cos(i),-1*np.sin(i)*np.sin(baz),-1*np.sin(i)*np.cos(baz)],
                       [np.sin(i),np.cos(i)*np.sin(baz),np.cos(i)*np.cos(baz)],
                       [0,-1*np.cos(baz),np.sin(baz)]
                      ])
        return R

    stz = seispy.convert.master_set(stz)
    stn = seispy.convert.master_set(stn)
    ste = seispy.convert.master_set(ste)
    stz_list = []
    ste_list = []
    stn_list = []

    for tr in stz:
        stz_list.append(tr.stats.sortname)
    for tr in stn:
        stn_list.append(tr.stats.sortname)
    for tr in ste:
        ste_list.append(tr.stats.sortname)

    d = [stz_list,ste_list,stn_list]
    'Find common stations in array'
    common = list(reduce(set.intersection, [set(item) for item in d ]))

    stz_list = []
    ste_list = []
    stn_list = []
    for tr in stz:
        if tr.stats.sortname not in common:
            stz.remove(tr)
        elif tr.stats.sortname in stz_list:
            stz.remove(tr)
        stz_list.append(tr.stats.sortname)
    for tr in ste:
        if tr.stats.sortname not in common:
            ste.remove(tr)
        elif tr.stats.sortname in ste_list:
            ste.remove(tr)
        ste_list.append(tr.stats.sortname)
    for tr in stn:
        if tr.stats.sortname not in common:
            stn.remove(tr)
        elif tr.stats.sortname in stn_list:
            stn.remove(tr)
        stn_list.append(tr.stats.sortname)

    stz = seispy.convert.master_set(stz)
    stn = seispy.convert.master_set(stn)
    ste = seispy.convert.master_set(ste)

    stl = stz.copy()
    stq = stz.copy()
    stt = stz.copy()

    for idx,tr in enumerate(stz):

        baz = np.radians(obspy.geodetics.base.gps2dist_azimuth(
              tr.stats.evla,
              tr.stats.evlo,
              tr.stats.stla,
              tr.stats.stlo)[-1])
        gcarc = tr.stats.gcarc
        h = tr.stats.evdp
        arrivals = model.get_travel_times(source_depth_in_km=h,
                                          distance_in_degree=gcarc,
                                          phase_list=phase)
        i = np.radians(arrivals[0].incident_angle)
        R = make_R(i,baz)
        zen = np.vstack((stz[idx].data,ste[idx].data,stn[idx].data))
        lqt = np.dot(R,zen)
        #zen = [stz[idx].data,ste[idx].data,stn[idx].data]
        #lqt = R*zen
        stl[idx].data = np.array(lqt[0])[0,:]
        stq[idx].data = np.array(lqt[1])[0,:]
        stt[idx].data = np.array(lqt[2])[0,:]

    return stl,stq,stt

def roll_zero(array,n):
    '''
    Roll and pad with zeros
    '''
    if n < 0:
        array = np.roll(array,n)
        array[n::] = 0
    else:
        array = np.roll(array,n)
        array[0:n] = 0
    return array

def clip_traces(st_in):
    st = st_in.copy()
    st_out = obspy.core.Stream()
    depth = st[0].stats.sac['evdp']
    samp = st[0].stats.sampling_rate
    for tr in st:
        starttime = tr.stats.starttime
        gcarc = tr.stats.sac['gcarc']
        arrival = model.get_travel_times(distance_in_degree=gcarc,
                               source_depth_in_km=depth,
                               phase_list=['P'])
        phasetime = arrival[0].time
        cuttimestart = starttime+phasetime-20+tr.stats.sac['o']
        cuttimeend = starttime+phasetime+300+tr.stats.sac['o']

        st_out.append(tr.slice(cuttimestart,cuttimeend))
        #time_ind = int(time*samp)
        #start_ind = int((time-100)*samp)
        #end_ind = int((time+500)*samp)
        #tr.data = tr.data[start_ind:end_ind]
    return st_out

def roll(st,seconds):
    '''
    Shift all traces in stream by seconds
    '''
    for tr in st:
        shift = int(tr.stats.sampling_rate*seconds)
        tr.data = np.roll(tr.data,shift)
    return st

def phase_window(tr,phase,**kwargs):
    '''
    return window around PKIKP phase
    '''
    tr = tr.copy()
    window_tuple = kwargs.get('window',(-10,10))
    in_model = kwargs.get('model','prem50')

    if type(in_model) == str:
        model = TauPyModel(model=in_model)
    else:
        model = in_model

    tr.stats.distance = tr.stats.sac['gcarc']
    origin_time = tr.stats.sac['o']
    start = tr.stats.starttime

    time = model.get_travel_times(source_depth_in_km = tr.stats.sac['evdp'],
                                       distance_in_degree = tr.stats.sac['gcarc'],
                                       phase_list = phase)
    t = time[0].time+origin_time
    PKPPKP_tr = tr.slice(start+t+window_tuple[0],start+t+window_tuple[1])
    #PKPPKP_tr.stats.sac['o'] += -1*window_tuple[0]
    return PKPPKP_tr

def phase_mask(tr_in,phase,**kwargs):
    '''
    return window around PKIKP phase
    '''
    tr = tr_in.copy()
    window_len = kwargs.get('window_len',20)
    in_model = kwargs.get('model','prem50')

    if type(in_model) == str:
        model = TauPyModel(model=in_model)
    else:
        model = in_model

    tr.stats.distance = tr.stats.sac['gcarc']
    origin_time = tr.stats.sac['o']
    start = tr.stats.starttime

    sr = tr.stats.sampling_rate
    time = model.get_travel_times(source_depth_in_km = tr.stats.sac['evdp'],
                                       distance_in_degree = tr.stats.sac['gcarc'],
                                       phase_list = phase)
    t = time[0].time+origin_time
    window = tukey(int(sr*window_len),0.2)
    mask = np.zeros(len(tr.data))
    mask[int(sr*t):int(sr*t)+int(sr*window_len)] = window
    tr.data *= mask
    return tr

def tr_align_on_phase(tr, **kwargs):
    '''
    Use to precisely align seismogram on phase
    '''
    phase = kwargs.get('phase',['P'])
    a_min = kwargs.get('min',True)
    window_tuple = kwargs.get('window',(-20,20))
    in_model = kwargs.get('model','prem50')

    model = TauPyModel(model=in_model)
    def roll_zero(array,n):
        if n < 0:
            array = np.roll(array,n)
            array[n::] = 0
        else:
            array = np.roll(array,n)
            array[0:n] = 0
        return array

    arrivals = model.get_travel_times(distance_in_degree=tr.stats.sac['gcarc'],
                     source_depth_in_km=tr.stats.sac['evdp'],
                     phase_list = phase)
    P = arrivals[0]
    t = tr.stats.starttime
    o = tr.stats.sac['o']
    window_data = (tr.slice(t+P.time+window_tuple[0]+o,
                    t+P.time+window_tuple[1]+o).data)
    if a_min:
        min_P = window_data.min()
        imin = np.argmin(np.abs(min_P-window_data))
        shift = int(len(window_data)/2.)-imin
        tr.data = np.roll(tr.data,(1*shift))
    else:
        max_P = window_data.max()
        imax = np.argmin(np.abs(max_P-window_data))
        shift = int(len(window_data)/2.)-imax
        tr.data = np.roll(tr.data,(1*shift))

    return tr

def align_on_phase(st, **kwargs):
    '''
    Use to precisely align seismogram on phase
    '''
    phase = kwargs.get('phase',['P'])
    a_min = kwargs.get('min',True)
    window_tuple = kwargs.get('window',(-30,30))
    in_model = kwargs.get('model','prem')

    model = TauPyModel(model=in_model)
    def roll_zero(array,n):
        if n < 0:
            array = np.roll(array,n)
            array[n::] = 0
        else:
            array = np.roll(array,n)
            array[0:n] = 0
        return array

    for tr in st:

        arrivals = model.get_travel_times(distance_in_degree=tr.stats.sac['gcarc'],
                         source_depth_in_km=tr.stats.sac['evdp'],
                         phase_list = phase)
        P = arrivals[0]
        t = tr.stats.starttime
        o = tr.stats.sac['o']
        window_data = (tr.slice(t+P.time+window_tuple[0]+o,
                        t+P.time+window_tuple[1]+o).data)
        if a_min:
            min_P = window_data.min()
            imin = np.argmin(np.abs(min_P-window_data))
            shift = int(len(window_data)/2.)-imin
            tr.data = np.roll(tr.data,(1*shift))
        else:
            max_P = window_data.max()
            imax = np.argmin(np.abs(max_P-window_data))
            shift = int(len(window_data)/2.)-imax
            tr.data = np.roll(tr.data,(1*shift))

    return st

def phase_align(st,**kwargs):
    phase = kwargs.get('phase',['P'])
    pol = kwargs.get('pol','min')
    window = kwargs.get('window',(-30,30))
    in_model = kwargs.get('model',model)

    for tr in st:

        arrivals = in_model.get_travel_times(distance_in_degree=tr.stats.sac['gcarc'],
                         source_depth_in_km=tr.stats.sac['evdp'],
                         phase_list = phase)
        P = arrivals[0]
        t = tr.stats.starttime
        o = tr.stats.sac['o']
        sr = tr.stats.sampling_rate

        if pol == 'min':
            sl = tr.slice(t+P.time+o+window[0],t+P.time+o+window[1]).data
            mi = np.argmin(sl)
        if pol == 'max':
            sl = tr.slice(t+P.time+o+window[0],t+P.time+o+window[1]).data
            mi = np.argmax(sl)

        shift = int(len(sl)/2.)-mi
        s_shift = shift/sr
        o *= -1*shift
    return st

def align_on_correlation(st, **kwargs):
    '''
    Use to precisely align seismogram on phase
    '''

    phase = kwargs.get('phase',['P'])
    a_min = kwargs.get('min',True)
    window_tuple = kwargs.get('window_tuple',(-5,5))
    ref_tuple = kwargs.get('ref_window',(-50,50))

    tr = st[0]
    arrivals = model.get_travel_times(distance_in_degree=tr.stats.sac['gcarc'],
                     source_depth_in_km=tr.stats.sac['evdp'],
                     phase_list = phase)
    P = arrivals[0]
    t = tr.stats.starttime
    o = tr.stats.sac['o']

    ref_data = (tr.slice(t+P.time+ref_tuple[0]+o,
                    t+P.time+ref_tuple[1]+o).data)


    for tr in st[1::]:
        arrivals = model.get_travel_times(distance_in_degree=tr.stats.sac['gcarc'],
                     source_depth_in_km=tr.stats.sac['evdp'],
                     phase_list = phase)
        P = arrivals[0]
        t = tr.stats.starttime
        o = tr.stats.sac['o']
        tr_data = (tr.slice(t+P.time+ref_tuple[0]+o,
                    t+P.time+ref_tuple[1]+o).data)
        cor = scipy.signal.fftconvolve(ref_data,tr_data[::-1])
        midpoint = cor.shape[0]/2
        imax = np.where(cor == cor.max())[0][0]
        roll = -1*(midpoint-imax)
        tr.data = np.roll(tr.data,roll)
    return st

def normalize_on_phase(st,**kwargs):
    '''
    normalize traces in stream based on maximum value in phase window
    '''
    phase = kwargs.get('phase',['P'])
    window_tuple = kwargs.get('window_tuple',(-10,10))
    in_model = kwargs.get('model','prem50')

    for tr in st:
        window = phase_window(tr,phase,window=window_tuple,model=in_model)
        tr.data = tr.data/np.abs(window.data).max()
    return st

def trace_normalize_on_phase(tr,**kwargs):
    '''
    normalize traces in stream based on maximum value in phase window
    '''
    phase = kwargs.get('phase',['P'])
    window_tuple = kwargs.get('window_tuple',(-10,10))

    window = phase_window(tr,phase,window=window_tuple)
    tr.data = tr.data/np.abs(window.data).max()
    return tr

def normalize_on_phase_range(st,**kwargs):
    phase = kwargs.get('phase',['P'])
    window_tuple = kwargs.get('window_tuple',(-100,100))

    for tr in st:
        window = phase_window(tr,phase,window=window_tuple)
        tr.data = tr.data/np.mean(np.abs(window.data))
    return st

def normalize_on_envelope(st,**kwargs):
    phase = kwargs.get('phase',['S'])
    window_tuple = kwargs.get('window_tuple',(-100,100))

    env_st = st.copy()

    for idx,tr in enumerate(st):
       env_st[idx].data = obspy.signal.filter.envelope(tr.data)
       window = phase_window(env_st[idx],phase,window=window_tuple)
       tr.data *= 1./np.max(window.data)

    return st

def periodic_corr(data, deconvolution):
    '''
    Periodic correlation, implemented using the FFT.
    data and deconvolution must be real sequences with the same length.

    Designed to align deconvolved trace with convolved trace.

    Use np.roll to shift deconvolution by value returned.
    '''
    corr = np.fft.ifft(np.fft.fft(data) * np.fft.fft(deconvolution).conj()).real
    shift = np.where(corr == corr.max())[0][0]
    return shift


def RT_ratio(stt,str):
    '''
    find ratio of Radial to Transverse component for every azimuth
    '''
    rtlist = []
    for idx,tr in enumerate(stt):
        if tr.stats.station != str[idx].stats.station:
            print 'Mismatch'
            continue
        trt = phase_window(tr,['S'],(-10,10))
        trr = phase_window(str[idx],['S'],(-10,10))
        T = trt.data.max()
        R = trr.data.max()
        diff = (T-R)**2
        dist = tr.stats.sac['gcarc']
        rtlist.append([diff,dist])
        rtarray = np.array(rtlist)

    fig, ax = plt.subplots()
    ax.scatter(rtarray[:,1],np.log(rtarray[:,0]),alpha=0.8,marker='D')
    #ax.set_ylim(np.log(rtarray[:,0]).mean()-20,np.log(rtarray[:,0].mean()+20))
    ax.set_xlabel('Source reciever azimuthal distance (deg)')
    ax.set_ylabel(u'$log((Amplitude(S_{T})-Amplitude(S_{R}))^{2})$')
    ax.grid()
    coeff = np.polyfit(rtarray[:,1],np.log(rtarray[:,0]),1)
    p = np.poly1d(coeff)
    x = np.linspace(rtarray[:,1].min(),rtarray[:,1].max())
    plt.plot(x,p(x),c='k')
    plt.show()

def slant(st_in,slowness):
    '''
    shift stack by slowness relative to center slowness
    '''
    st = st_in.copy()
    gc_list = []
    for tr in st:
        gc_list.append(tr.stats.gcarc)

    mean_dist = np.mean(gc_list)
    samp_rate = st[0].stats.sampling_rate

    for tr in st:
        dist = tr.stats.gcarc-mean_dist
        tr.data = roll_zero(tr.data,int(slowness*dist*samp_rate))
    return st

def stack(st):
    '''
    shift stack by slowness relative to center slowness
    '''
    s = []
    for tr in st:
        s.append(tr.data)

    return np.mean(s,axis=0)




