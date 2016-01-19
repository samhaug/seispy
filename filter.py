#!/usr/bin/env python

import scipy
import obspy
from obspy.taup import TauPyModel
import numpy as np
model = TauPyModel(model="prem")

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


def dirty_filter(st):
    '''
    Remove trace from stream if h
    '''

    for tr in st:
        arrivals = model.get_travel_times(
                   distance_in_degree=tr.stats.sac['gcarc'],
                   source_depth_in_km=tr.stats.sac['evdp'],
                   phase_list=['P'])
        P = arrivals[0]
        t = tr.stats.starttime
        o = tr.stats.sac['o']
        max_P = abs(tr.slice(t+P.time-5+o, t+P.time+5+o).data).max()
        pre_noise = abs(tr.slice(t+P.time-50+o, t+P.time-20+o).data).max()
        post_noise = abs(tr.slice(t+P.time+10+o, t+P.time+30+o).data).max()
        if (pre_noise > max_P*0.30) or (post_noise > max_P*0.3):
            st.remove(tr)

    return st

def gimp_filter(st):
    '''
    Removes seismograms from trace if they have lengths too short. Makes all
    seismograms the same length and same sampling rate
    '''

    def max_len(st):
        a = []
        for tr in st:
            a.append(tr.data.shape[0])
        return max(a)

    def min_len(st):
        a = []
        for tr in st:
            a.append(tr.data.shape[0])
        return min(a)

    for tr in st:
        if tr.data.shape[0] < 100:
            st.remove(tr)

    st.interpolate(sampling_rate=50.0)

    mx_len = max_len(st)

    for tr in st:
        if tr.data.shape[0] < mx_len-10:
            st.remove(tr)

    mn_len = min_len(st)

    for tr in st:
        tr.data = tr.data[0: mn_len]

    return st

def range_filter(st, range_tuple):
    '''
    Removes seismograms from trace if they fall outside the range limits
    of range_tuple

    range_tuple = (30,50) removes all traces outside of 30 to 50 degrees from
    source
    '''

    for tr in st:
        if (tr.stats.sac['gcarc'] < range_tuple[0]):
            st.remove(tr)
        elif (tr.stats.sac['gcarc'] > range_tuple[1]):
            st.remove(tr)

    return st
