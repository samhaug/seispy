#!/usr/bin/env python

import scipy
import obspy
from obspy.taup import TauPyModel
import numpy as np
model = TauPyModel(model="prem")
import scipy.optimize
from matplotlib import pyplot as plt

'''
Samuel Haugland 01/19/16

seis_filter.py includes functions needed to remove unwanted traces from
streams based on various criteria. All functions should take a stream object
and arguments and return the filtered stream object
'''

def kurtosis_filter(st, **kwargs):
    '''
    remove traces from phase based on kurtosis
    '''

    alpha = kwargs.get('alpha', False)
    if alpha is not False:
        alpha = 0.5

    k = []
    for tr in st:
        ki = scipy.stats.kurtosis(tr.data)
        if np.isnan(ki):
            st.remove(tr)
            continue
        else:
            k.append(ki)
    mean_k = sum(k)/len(st)

    for tr in st:
        if scipy.stats.kurtosis(tr.data) < (alpha*mean_k):
            st.remove(tr)
    return st

def dirty_filter(st,**kwargs):
    '''
    Remove trace from stream if noise is too much
    a,b  refer to the time windows before the phase
    c,d  refer to the time windows after the phase
    '''

    a = kwargs.get('a',50)
    b = kwargs.get('b',20)
    c = kwargs.get('c',10)
    d  = kwargs.get('d',30)

    phase = kwargs.get('phase',['P'])

    pre_limit = kwargs.get('pre_limit',0.3)
    post_limit = kwargs.get('post_limit',0.3)

    for tr in st:
        arrivals = model.get_travel_times(
                   distance_in_degree=tr.stats.sac['gcarc'],
                   source_depth_in_km=tr.stats.sac['evdp'],
                   phase_list=phase)

        P = arrivals[0]
        t = tr.stats.starttime
        o = tr.stats.sac['o']
        max_P = abs(tr.slice(t+P.time-3+o, t+P.time+3+o).data).max()
        pre_noise = abs(tr.slice(t+P.time-a+o, t+P.time-b+o).data).max()
        post_noise = abs(tr.slice(t+P.time+c+o, t+P.time+d+o).data).max()
        if (pre_noise > max_P*pre_limit) or (post_noise > max_P*post_limit):
            st.remove(tr)

    return st

def gimp_filter(st,**kwargs):
    len_cutoff = kwargs.get('len_cutoff',0.8)
    '''
    Removes seismograms from trace if they have lengths too short. Makes all
    seismograms the same length and same sampling rate
    '''

    len_list = []
    for tr in st:
        len_list.append(tr.stats.endtime-tr.stats.starttime)
        #if tr.stats.npts < 500:
        #    st.remove(tr)

    mx_len = np.max(len_list)
    mn_len = np.min(len_list)

    for tr in st:
        if tr.stats.npts < mx_len*len_cutoff:
            st.remove(tr)
            continue
        if np.isnan(sum(tr.data)):
            st.remove(tr)
    st.interpolate(20)
    return st

def range_filter(st, range_tuple):
    '''
    Removes seismograms from trace if they fall outside the range limits
    of range_tuple

    range_tuple = (30,50) removes all traces outside of 30 to 50 degrees from
    source
    '''

    for tr in st:
        if not range_tuple[0] <= tr.stats.gcarc <= range_tuple[1]:
            st.remove(tr)

    return st

def az_filter(st, az_tuple):
    '''
    Removes seismograms from trace if they fall outside the azimuth limits
    of az_tuple
    '''

    for tr in st:
        if not az_tuple[0] <= tr.stats.az <= az_tuple[1]:
            st.remove(tr)

    return st

def monochrome(tr,**kwargs):
    '''Fit sin to the input time sequence, and return fitting parameters "amp", "omega", "phase", "offset", "freq", "period" and "fitfunc"'''
    cutoff = kwargs.get('cutoff',30)
    guess = kwargs.get('guess',False)

    tt = np.linspace(0,tr.stats.endtime-tr.stats.starttime,num=tr.stats.npts)
    yy = tr.data
    ff = np.fft.fftfreq(len(tt), (tt[1]-tt[0]))   # assume uniform spacing
    Fyy = abs(np.fft.fft(yy))
    if guess == False:
        guess_freq = abs(ff[np.argmax(Fyy[1:])+1])   # excluding the zero frequency "peak", which is related to offset
    else:
        guess_freq = guess
    guess_amp = np.std(yy) * 2.**0.5
    guess_offset = np.mean(yy)
    guess = np.array([guess_amp, 2.*np.pi*guess_freq, 0., guess_offset])

    def sinfunc(t, A, w, p, c):  return A * np.sin(w*t + p) + c
    try:
        popt, pcov = scipy.optimize.curve_fit(sinfunc, tt, yy, p0=guess)
    except RuntimeError:
        return tr
    A, w, p, c = popt
    f = w/(2.*np.pi)
    fitfunc = lambda t: A * np.sin(w*t + p) + c
    res =  {"amp": A, "omega": w, "phase": p, "offset": c, "freq": f, "period": 1./f,
             "fitfunc": fitfunc, "maxcov": np.max(pcov), "rawres": (guess,popt,pcov)}
    if res['period'] > cutoff:
        tr.data += -1*res['fitfunc'](tt)
        return tr
    else:
        return tr


