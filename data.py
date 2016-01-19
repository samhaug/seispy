#!/usr/bin/env python

import scipy
import obspy
import numpy as np
from obspy.taup import TauPyModel
model = TauPyModel(model="prem")

###############################################################################
def align_on_phase(st,phase):
###############################################################################
    '''
    Use to precisely align seismogram on phase
    '''
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
        t+P.time+o
        window_data = (tr.slice(t+P.time-10+o,t+P.time+10+o).data)
        max_P = window_data[window_data > 0].max()
        imax = np.argmin(np.abs(max_P-window_data))
        shift = int(len(window_data)/2.)-imax
        tr.data = np.roll(tr.data,(1*shift))
    return st

###############################################################################
def periodic_corr(data, deconvolution):
###############################################################################
    '''
    Periodic correlation, implemented using the FFT.
    data and deconvolution must be real sequences with the same length.

    Designed to align deconvolved trace with convolved trace.

    Use np.roll to shift deconvolution by value returned.
    '''
    corr = np.fft.ifft(np.fft.fft(data) * np.fft.fft(deconvolution).conj()).real
    shift = np.where(corr == corr.max())[0][0]
    return shift

###############################################################################
def waterlevel_deconvolve(tr):
###############################################################################
    '''
    Water_level_deconvolve trace
    Sssume seismic wavelet is the P phase.
    works for both .xh and .sac formats
    '''

    def isolate_source(tr):
        '''
        Find source time function by multiplying P arrival with Tukey window
        '''

        stats_dict = {}

        if tr.stats._format == 'XH':
            stats_dict['depth'] = tr.stats.xh['source_depth_in_km']
            stats_dict['distance'] = float(tr.stats.station.split('_')[0])+float(tr.stats.station.split('_')[1])/100.
            stats_dict['start_time'] = tr.stats.starttime
            stats_dict['end_time'] = tr.stats.endtime
            stats_dict['offset'] = 0.
            stats_dict['sampling_rate'] = tr.stats.sampling_rate
        elif tr.stats._format == 'SAC':
            stats_dict['depth'] = tr.stats.sac['evdp']
            stats_dict['distance'] = tr.stats.sac['gcarc']
            stats_dict['start_time'] = tr.stats.starttime
            stats_dict['end_time'] = tr.stats.endtime
            stats_dict['offset'] = tr.stats.sac['o']
            stats_dict['sampling_rate'] = tr.stats.sampling_rate

        if stats_dict['end_time']-stats_dict['start_time'] <= 0:
            #return 'REMOVE'
            raise ValueError('starttime is larger than endtime')

        taup_time = model.get_travel_times(source_depth_in_km = stats_dict['depth'],
                                      distance_in_degree = stats_dict['distance'],
                                      phase_list = ['P'])
        P_arrival_time = taup_time[0].time+stats_dict['offset']


        begin_pad = tr.slice(stats_dict['start_time'],stats_dict['start_time']+P_arrival_time).data.size
        P_array = tr.slice(stats_dict['start_time']+P_arrival_time,stats_dict['start_time']+P_arrival_time+5).data.size
        end_pad = tr.slice(stats_dict['start_time']+P_arrival_time+5,stats_dict['end_time']).data.size

        tukey = scipy.signal.tukey(P_array,alpha=0.6)
        tukey_pad = np.pad(tukey,(begin_pad-1,end_pad-1),'constant',constant_values=(0,0))
        wavelet = tukey_pad*tr.data
        wavelet = np.roll(wavelet,-1*begin_pad-(P_array/2))

        return wavelet,tr.data,stats_dict

    def apply_waterlevel(tukey_pad,trace_data,alpha):
        '''
        Apply water level to fill spectral holes. see Practical Seismic Data Analysis
        page 182
        '''

        tukey_omega = np.fft.fft(tukey_pad)
        trace_omega = np.fft.fft(trace_data)

        F_omega = tukey_omega*tukey_omega.conjugate()
        F_omega[F_omega < alpha*F_omega.max()] = alpha*F_omega.max()
        out = (trace_omega*tukey_omega.conjugate())/F_omega
        out = np.fft.ifft(out)
        return out.real

    pad, data, stats_dict = isolate_source(tr)
    out = apply_waterlevel(pad,data,0.1)
    out = out-out.mean()
    out = out/out.max()

    if np.isnan(np.sum(out)):
        #return 'REMOVE'
        raise ValueError('NaN in deconvolved trace')

    return [stats_dict,out]

