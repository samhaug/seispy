#!/usr/bin/env python

import numpy as np
from matplotlib import pyplot as plt
import obspy
import seispy
from obspy.taup import TauPyModel
model = TauPyModel(model='prem')

def gui_pick(st):

    coordx_list = []
    coordy_list = []
    label_list = []
    def coord_pick(event):
        x = event.xdata
        y = event.ydata
        print x,y
        coordx_list.append(round(x,3))
        coordy_list.append(round(y,3))
    def label_pick(event):
        artist = event.artist
        trace_label = artist.get_label()
        print trace_label
        label_list.append(str(trace_label))

    for ii in range(0,(len(st)-len(st)%5),5):
        fig,ax = plt.subplots(5,figsize=(10,15))

        for idx, axes in enumerate(ax):
            tlen0 = st[ii+idx].stats.endtime-st[ii+idx].stats.starttime
            start0 = -1*st[ii+idx].stats.sac['o']
            end0 = tlen0 - st[ii+idx].stats.sac['o']
            time0 = np.linspace(start0,end0,num=st[ii+idx].stats.npts)
            label0 = str(st[ii+idx].stats.network)+str(st[ii+idx].stats.station)
            axes.plot(time0,st[ii+idx].data,color='k',label=label0,picker=10)
            arrivals = model.get_travel_times(
                       distance_in_degree = st[ii+idx].stats.sac['gcarc'],
                       source_depth_in_km = st[ii+idx].stats.sac['evdp'],
                       phase_list =['P'])
            if len(arrivals) != 0:
                p = arrivals[0]
                axes.axvline(p.time,color='b')
                axes.set_xlim(p.time-10,p.time+10)
        fig.canvas.mpl_connect('button_press_event', coord_pick)
        fig.canvas.mpl_connect('pick_event', label_pick)
        plt.show()

    coord_list = zip(coordx_list,coordy_list)

    return coord_list, label_list

def stream_remove(st,coord_list,label_list):

    def make_stat_list(st):
        stat_list = []
        for tr in st:
            stat_list.append(str(tr.stats.network)+str(tr.stats.station))
        return stat_list

    def edit_data(tr,arrival_time,selected_time,flip):
        tlen0 = tr.stats.endtime-tr.stats.starttime
        start0 = -1*tr.stats.sac['o']
        end0 = tlen0 - tr.stats.sac['o']
        time0 = np.linspace(start0,end0,num=tr.stats.npts)
        samp = tr.stats.sampling_rate

        select_idx = np.argmin(abs(time0-selected_time))
        arrive_idx = np.argmin(abs(time0-arrival_time))

        window = np.abs(tr.data[select_idx-(int(0.1*samp)):select_idx+(int(0.1*samp))])
        win_max = np.where(window == window.max())[0][0]
        pick_idx = select_idx-int(0.1*samp)+win_max

        tr.data = np.roll(tr.data,(arrive_idx-pick_idx))
        if flip > 0:
            return tr
        elif flip < 0:
            tr.data *= -1
            return tr

    for idx, tr in enumerate(st):
        if str(tr.stats.network)+str(tr.stats.station) not in label_list:
            st.remove(tr)
        else:
            ii = label_list.index(str(tr.stats.network)+str(tr.stats.station))
            selected_time = coord_list[ii][0]
            flip = coord_list[ii][1]
            arrivals = model.get_travel_times(
                       distance_in_degree = tr.stats.sac['gcarc'],
                       source_depth_in_km = tr.stats.sac['evdp'],
                       phase_list =['P'])
            p = arrivals[0]
            arrival_time = p.time

            st[ii] = edit_data(st[ii],arrival_time,selected_time,flip)
    return st












