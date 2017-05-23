#!/home/samhaug/anaconda2/bin/python

'''
==============================================================================

File Name : io.py
Purpose : read_write functions
Creation Date : 23-05-2017
Last Modified : Tue 23 May 2017 02:13:28 PM EDT
Created By : Samuel M. Haugland

==============================================================================
'''

import numpy as np
import h5py
import obspy
import seispy

def h5write(st,name):
    '''write a ghetto h5py file to save obspy stream'''
    f = h5py.File(name,'w')
    for idx,tr in enumerate(st):
        f.create_group(str(idx))
        f[str(idx)].create_dataset('data',data=tr.data)
        f[str(idx)].create_group('sampling_rate/./'+tr.stats.sampling_rate)
        f[str(idx)].create_group('o/./'+tr.stats.sac['o'])
        f[str(idx)].create_group('evdp/./'+tr.stats.sac['evdp'])
        f[str(idx)].create_group('stla/./'+tr.stats.sac['stla'])
        f[str(idx)].create_group('stlo/./'+tr.stats.sac['stlo'])
        f[str(idx)].create_group('evla/./'+tr.stats.sac['evla'])
        f[str(idx)].create_group('evlo/./'+tr.stats.sac['evlo'])
        f[str(idx)].create_group('gcarc/./'+tr.stats.sac['gcarc'])
        f[str(idx)].create_group('az/./'+tr.stats.sac['az'])
        f[str(idx)].create_group('bz/./'+tr.stats.sac['baz'])
        f[str(idx)].create_group('starttime/./'+tr.stats.starttime)
        f[str(idx)].create_group('location/./'+tr.stats.location)
        f[str(idx)].create_group('station/./'+tr.stats.station)
        f[str(idx)].create_group('network/./'+tr.stats.station)



