#!/usr/bin/env python

import h5py
import numpy as np
import obspy
from geopy.distance import great_circle
from obspy.taup import TauPyModel
model = TauPyModel(model="prem_50")
from matplotlib import pyplot as plt

def read_norm():

    f = h5py.File('env.h5','r')
    e = f['env'][...]
    d = f['num'][...]
    return (e,d)

def env_plot(env_array):
    log_env = np.log10(np.flipud(env_array))
    fig, ax = plt.subplots()
    image = ax.imshow(log_env,aspect='auto',interpolation='none',vmin=-1.8,vmax=-0.0,
                     cmap='Spectral_r',extent=[-20,20,125,140])

    cbar = fig.colorbar(ax=ax,mappable=image)
    ax.set_xlabel('Seconds before PKIKP')
    ax.set_ylabel('Range')
    plt.show()

def env_line(env_array):
    for idx,ii in enumerate(env_array):
        plt.plot(ii+idx,color='k')
        plt.axhline(idx,color='k')
    plt.show()

def section(st):
    fig, ax = plt.subplots()
    time = np.linspace(-200,200,num=len(st[0].data))
    for tr in st:
        ax.plot(time,tr.data+tr.stats.sac['gcarc'],color='k',alpha=0.5)
    plt.show()

def sod_express(st):
    st = SOD_evdp(st)
    st = set_origin_time(st)
    return st

def equalize_start(st):
    '''
    equalize_start_time
    '''
    start = st[0].stats.starttime
    for tr in st:
        tr.stats.starttime = start
    return st

def remove_noise(st):
    '''
    remove signals with noise before arrival
    '''
    for tr in st:
        dat_len = len(tr.data)
        signal = tr.data[int(dat_len/2.-200):int(dat_len/2.+200)]
        noise = tr.data[2000:6000]
        noise_l = np.mean(abs(noise))
        sig_l = np.mean(abs(signal))
        if noise_l > 0.2*sig_l:
            st.remove(tr)
    return st

def save_envelope(env_array,bootstrap_array,h5_file):
    f = h5py.File(h5_file,'w')
    f.create_dataset('env',data=env_array)
    f.create_dataset('num',data=bootstrap_array)
    f.close()

def set_envelope(st):
    '''
    set envelope for PKIKP wave
    '''
    st.sort(['location'])
    env_st = obspy.core.Stream()
    rb_min = np.arange(125,140.5,0.5)
    zero_shape = st[0].data.shape
    env_list = []
    bootstrap_env = []
    for ii in rb_min:
        new = st.select(location = str(ii))
        bootstrap_env.append([ii, len(new)])
        print ii, len(new)
        if len(new) == 0:
            env_list.append(np.zeros(zero_shape))
        else:
            avg_env = np.zeros(new[0].data.shape)
            for tr in new:
                tr.data = np.nan_to_num(tr.data)
                env = obspy.signal.envelope(tr.data)
                env = np.nan_to_num(env)
                env = env - np.min(np.abs(env[1000:3000]))
                avg_env += env
            avg_env = avg_env/len(new)
            env_list.append(avg_env)
    env_array = np.array(env_list)
    env_array = env_array[:,8000:12000]

    for idx, ii in enumerate(env_array):
        env_array[idx] = ii-ii.min()
    return env_array,np.array(bootstrap_env)

def PKIKP_range_bin(st):
    '''
    get envelope for data in range bins
    '''

    def find_nearest(array,value):
        idx = (np.abs(array-value)).argmin()
        return array[idx]

    rb_min = np.arange(125,140.5,0.5)

    for tr in st:
        nearest = find_nearest(rb_min,np.round(tr.stats.sac['gcarc'],1))
        tr.stats.sac['gcarc'] =round(nearest,1)
        tr.stats.location = str(nearest)
    return st

def zero_PKIKP(st):
    '''
    For a trace requested from SOD, remove any noise at the end of the trace
    so a normalizatino won't be skewed by later arrivals
    '''
    for tr in st:
        datalen = len(tr.data)
        tr.data[int(datalen*(2./3.))::] = 0
    return st

def set_origin_time(st):
    '''
    set sac['o'] time for events retrieved from SOD.
    '''
    event_depth = st[0].stats.sac['evdp']
    for tr in st:
        arrivals = model.get_travel_times(source_depth_in_km=event_depth,
                  distance_in_degree=tr.stats.sac['gcarc'],phase_list=['PKIKP'])
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
