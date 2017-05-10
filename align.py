#!/home/samhaug/anaconda2/bin/python

'''
==============================================================================

File Name : align.py
Purpose : module for manually aligning data
Creation Date : 09-05-2017
Last Modified : Tue 09 May 2017 04:37:53 PM EDT
Created By : Samuel M. Haugland

==============================================================================
'''

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import pylab as p
from subprocess import call
from os import listdir
import h5py
import obspy
import seispy
from obspy.taup import TauPyModel
model = TauPyModel(model='prem')

def align(st_in,**kwargs):
    sec_type = kwargs.get('section_type','dist')
    a_list = kwargs.get('a_list',True)

    st = st_in.copy()
    drag_dict = {}
    for tr in st:
        tr.stats.identify = tr.stats.station
        drag_dict[tr.stats.identify] = 0

    st.interpolate(50)
    simple_section(st,section_type=sec_type,a_list=a_list)
    dragh = DragHandler(st)
    plt.show()

    shift_list = zip(dragh.station,dragh.shift)
    for ii in shift_list:
        drag_dict[ii[0]] += ii[1]

    for keys in drag_dict:
        for idx,tr in enumerate(st):
            if tr.stats.identify == keys:
                print keys,tr.stats.identify,int(drag_dict[keys]*tr.stats.sampling_rate)
                st[idx].data = np.roll(st[idx].data,int(drag_dict[keys]*tr.stats.sampling_rate))
            else:
                #print 'fuck'
                continue
    return st

def simple_section(st,**kwargs):
    '''
    Simpler section plotter for obspy stream object
    '''
    a_list = kwargs.get('a_list',True)
    fig = kwargs.get('fig',None)
    ax = kwargs.get('ax',None)
    color = kwargs.get('color','k')
    save = kwargs.get('save',False)
    picker = kwargs.get('picker',False)
    sec_type = kwargs.get('section_type','az')

    if sec_type == 'az':
        for tr in st:
            tr.stats.azimuth = tr.stats.sac['az']
        st.sort(['azimuth'])
    st = st[::-1]


    if fig == None and ax == None:
        fig,ax = plt.subplots(figsize=(10,15))
        ax.axvline(0)
    else:
        print('using outside figure')

    def plot(tr,o,ax):
        e = tr.stats.npts/tr.stats.sampling_rate
        t = np.linspace(o,o+e,num=tr.stats.npts)
        if sec_type == 'az':
            XX = tr.stats.sac['az']
        elif sec_type == 'dist':
            XX = tr.stats.sac['gcarc']
        else:
            XX = (tr.stats.sac['gcarc']*np.sin(np.radians(sec_type))+
                 tr.stats.sac['az']*np.cos(np.radians(sec_type)))

        ax.plot(t,(tr.data/0.5)+(XX),alpha=0.5,
                color=color,label=tr.stats.identify,
                picker=10)

    if a_list == True:
        for tr in st:
            plot(tr,0,ax)

    elif type(a_list) == list:
        if len(a_list) != 1:
            print('Must have phase identifier string of len = 1')
            return
        else:
            for tr in st:
                evdp = tr.stats.sac['evdp']
                gcarc = tr.stats.sac['gcarc']
                P = model.get_travel_times(distance_in_degree=gcarc,
                    source_depth_in_km=evdp,
                    phase_list = a_list)
                P_time = P[0].time
                plot(tr,-1*(P_time+tr.stats.sac['o']),ax)
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Epicentral Distance (deg)')
    ax.set_title(st[0].stats.network)

class DragHandler(object):
    """ A simple class to handle Drag n Drop.

    This is a simple example, which works for Text objects only.
    """
    def __init__(self, st, figure=None) :
        if figure is None : figure = p.gcf()
        # simple attibute to store the dragged text object
        self.dragged = None
        self.current_stream = None
        self.station = []
        self.shift = []
        self.srate = st[0].stats.sampling_rate

        # Connect events and callbacks
        figure.canvas.mpl_connect("pick_event", self.on_pick_event)
        figure.canvas.mpl_connect("button_release_event", self.on_release_event)

    def on_pick_event(self, event):
        #print type(event.artist.get_label())

        if type(event.artist.get_label()) == unicode:
            self.dragged = event.artist
            self.xdata = event.artist.get_data()[0]
            self.ydata = event.artist.get_data()[1]
            self.pick_pos = event.mouseevent.xdata
            self.station.append(event.artist.get_label())
        return True

    def on_release_event(self, event):
        newx = event.xdata
        try:
            newy = np.roll(self.ydata,int((newx-self.pick_pos)*(self.srate)))
            self.dragged.set_data(self.xdata,newy)
            self.dragged = None
            p.draw()
            self.shift.append(newx-self.pick_pos)
            return True
        except (AttributeError,TypeError):
            return True



