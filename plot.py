#!/usr/bin/env python

'''
seis_plot.py includes all functions for plotting seismic data
'''

import numpy as np
import obspy
from matplotlib import pyplot as plt
from obspy.taup import TauPyModel
from mpl_toolkits.basemap import Basemap
model = TauPyModel(model="prem")

def mapplot(st,**kwargs):
    gc = kwargs.get('great_circle',False)
    width = 28000000
    latav = []
    lonav = []
    for tr in st:
        latav.append(tr.stats.stla)
        lonav.append(tr.stats.stlo)
    latav = np.mean(latav)
    lonav = np.mean(lonav)

    fig1 = plt.figure(figsize=(10,10))
    lat_s = st[0].stats.evla
    lon_s = st[0].stats.evlo

    m = Basemap(projection='aeqd',lat_0=lat_s,lon_0=lon_s)
    xpt, ypt = m(lon_s,lat_s)
    m.scatter(xpt,ypt,s=99,c='red',marker='o',lw=1)
    m.drawcoastlines(linewidth=0.5)
    for tr in st:
        lat_0 = tr.stats.stla
        lon_0 = tr.stats.stlo
        xpt, ypt = m(lon_0, lat_0)
        if gc == True:
            m.drawgreatcircle(lon_s,lat_s,lon_0,lat_0,color='k',alpha=0.1,lw=0.4)
        m.scatter(xpt,ypt,s=5,c='green',marker='o',lw=0)

    '''
    fig2 = plt.figure(figsize=(10,10))
    m = Basemap(projection='ortho',lat_0=latav,lon_0=lonav)
    xpt, ypt = m(st[0].stats.sac['evlo'], st[0].stats.sac['evla'])
    m.scatter(xpt,ypt,s=99,c='red',marker='o',lw=1)
    # fill background.
    #m.drawmapboundary(fill_color='aqua')
    # draw coasts and fill continents.
    m.drawcoastlines(linewidth=0.5)
    #m.fillcontinents(color='coral',lake_color='aqua')
    # 20 degree graticule.
    m.drawparallels(np.arange(-80,81,20))
    m.drawmeridians(np.arange(-180,180,20))
    # draw a black dot at the center.
    for tr in st:
        lat_0 = tr.stats.sac['stla']
        lon_0 = tr.stats.sac['stlo']
        xpt, ypt = m(lon_0, lat_0)
        m.scatter(xpt,ypt,s=5,c='green',marker='o',lw=0)
    '''

    plt.show()

def plot(tr,**kwargs):
    '''
    plot trace object
    phase = optional list of phases. list of strings corresponding
            to taup names
    '''
    phases = kwargs.get('phase_list',False)
    window = kwargs.get('window',False)
    t_model = kwargs.get('model',model)

    fig, ax = plt.subplots(figsize=(23,9))
    if phases != False:
        arrivals = t_model.get_travel_times(
                   distance_in_degree=tr.stats.gcarc,
                   source_depth_in_km=tr.stats.evdp,
                   phase_list = phases)
        window_list = []
        colors = ['b','g','r','c','m','y','k']
        if len(arrivals) != 0:
            for idx, ii in enumerate(arrivals):
                ax.axvline(ii.time,label=ii.purist_name,c=np.random.rand(3,1))
                window_list.append(ii.time)
    ax.legend()

    time = np.linspace(-1*tr.stats.o,
           (tr.stats.delta*tr.stats.npts)-tr.stats.o,
           num=tr.stats.npts)
    ax.plot(time,tr.data,c='k')
    ax.grid()
    ax.set_title('Network: {}, Station: {}, Channel: {},\
    Dist (deg): {}, Depth (km): {} \nStart: {} \nEnd: {}'.format(
                  tr.stats.network,
                  tr.stats.station,
                  tr.stats.channel,
                  round(tr.stats.gcarc,3),
                  round(tr.stats.evdp,3),
                  tr.stats.starttime,
                  tr.stats.endtime))
    ax.set_xlabel('Time (s), sampling_rate: {}, npts: {}'
                 .format(tr.stats.sampling_rate,tr.stats.npts))

    if window != False:
        ax.set_xlim([min(window_list)+window[0],max(window_list)+window[1]])
        #ax.set_xlim([min(window_list)-300,max(window_list)+300])

    plt.show()

def compare_section(std_in,sts_in,**kwargs):
    '''
    compare two streams, std is data and sts is synthetic
    '''
    print('black is first, red is second')
    a_list = kwargs.get('a_list',True)
    fig = kwargs.get('fig',None)
    ax = kwargs.get('ax',None)
    model = kwargs.get('model','prem')
    mode = kwargs.get('mode','gcarc')
    model = TauPyModel(model=model)
    xlim = kwargs.get('x_lim',None)
    ylim = kwargs.get('y_lim',None)
    label1 = kwargs.get('label_1',None)
    label2 = kwargs.get('label_2',None)
    roll2 = kwargs.get('roll_2',0)
    roll1 = kwargs.get('roll_1',0)
    mult = kwargs.get('mult',1.0)
    norm = kwargs.get('norm',True)

    std = std_in.copy()
    sts = sts_in.copy()
    if norm:
        std.normalize()
        sts.normalize()

    if fig == None and ax == None:
        fig,ax = plt.subplots(figsize=(9,12))
        plt.tight_layout()

    for idx,tr in enumerate(std):
        o = tr.stats.o
        ds = tr.stats.starttime
        ddelt = tr.stats.delta
        dpts = tr.stats.npts
        P_time = 0
        if a_list != True:

            evdp = tr.stats.evdp
            gcarc = tr.stats.gcarc
            P = model.get_travel_times(distance_in_degree=gcarc,
                source_depth_in_km=evdp,
                phase_list = a_list)
            P_time += P[0].time
        t_dat = np.linspace(-o,dpts*ddelt-o,num=dpts)
        if roll1 != 0:
            tr.data = np.roll(tr.data,int(roll1*tr.stats.sampling_rate))
        if idx == 0:
            if mode == 'gcarc':
                ax.plot(t_dat-P_time,
                        tr.data*mult+tr.stats.gcarc,alpha=0.5,color='k',label=label1,
                        lw=0.5)
            elif mode == 'az':
                ax.plot(t_dat-P_time,
                        tr.data*mult+tr.stats.az,alpha=0.5,color='k',label=label1,
                        lw=0.5)
        else:
            if mode == 'gcarc':
                ax.plot(t_dat-P_time,
                        tr.data*mult+tr.stats.gcarc,alpha=0.5,color='k',lw=0.5)
            elif mode == 'az':
                ax.plot(t_dat-P_time,
                        tr.data*mult+tr.stats.az,alpha=0.5,color='k',lw=0.5)

    for idx,tr in enumerate(sts):
        o = tr.stats.o
        ss = tr.stats.starttime
        sdelt = tr.stats.delta
        spts = tr.stats.npts
        P_time = 0
        if a_list != True:
            evdp = tr.stats.evdp
            gcarc = tr.stats.gcarc
            P = model.get_travel_times(distance_in_degree=gcarc,
                                       source_depth_in_km=evdp,
                                       phase_list=a_list)
            P_time += P[0].time
        t_syn = np.linspace(-o,spts*sdelt-o,num=spts)
        if roll2 != 0:
            tr.data = np.roll(tr.data,int(roll2*tr.stats.sampling_rate))
        if idx == 0:
            if mode == 'gcarc':
                ax.plot(t_syn-P_time,tr.data*mult+tr.stats.gcarc,
                        alpha=0.5,color='r',label=label2,lw=0.5)
            elif mode == 'az':
                ax.plot(t_syn-P_time,tr.data*mult+tr.stats.az,
                        alpha=0.5,color='r',label=label2,lw=0.5)
        else:
            if mode == 'gcarc':
                ax.plot(t_syn-P_time,tr.data*mult+tr.stats.gcarc,
                        alpha=0.5,color='r',lw=0.5)
            elif mode == 'az':
                ax.plot(t_syn-P_time,tr.data*mult+tr.stats.az,
                        alpha=0.5,color='r',lw=0.5)

    plt.legend()
    ax.set_ylabel(mode)
    ax.set_xlabel('Time (s)')
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    plt.show()

def section(st,**kwargs):
    '''
    Simpler section plotter for obspy stream object
    '''
    plt.ioff()
    a_list = kwargs.get('a_list',True)
    fig = kwargs.get('fig',None)
    ax = kwargs.get('ax',None)
    color = kwargs.get('color','k')
    save = kwargs.get('save',False)
    picker = kwargs.get('picker',False)
    rainbow = kwargs.get('rainbow',False)
    mode = kwargs.get('mode','gcarc')
    xlim = kwargs.get('x_lim',None)
    ylim = kwargs.get('y_lim',None)
    title = kwargs.get('title','')
    model = kwargs.get('model','prem')
    mult = kwargs.get('mult',1.0)
    model = TauPyModel(model=model)
    figsize = kwargs.get('figsize',(9,12))

    if fig == None and ax == None:
        fig,ax = plt.subplots(figsize=figsize)
        plt.tight_layout()
    else:
        print('using outside figure')
    ax.set_xlabel('Time (s)')
    ax.set_ylabel(mode)
    ax.set_title(title)

    def plot(tr,o,ax,color,**kwargs):
        alpha = kwargs.get('alpha',0.5)
        e = tr.stats.npts/tr.stats.sampling_rate
        t = np.linspace(o,o+e,num=tr.stats.npts)

        ax.plot(t,(mult*tr.data/tr.data.max())+tr.stats.gcarc,alpha=alpha,
                color=color,label=tr.stats.network+'.'+tr.stats.station,
                picker=10,lw=0.8)

    def randcolor():
        c_list = ['#1f77b4','#ff7f0e','#2ca02c',
                   '#d62728','#9467bd','#8c564b',
                   '#e377c2','#7f7f7f','#bcbd22','#17becf']
        return c_list[np.random.randint(len(c_list))]

    if a_list == True:
        for tr in st:
            if rainbow == True:
                plot(tr,0,ax,randcolor(),alpha=0.7)
            else:
                plot(tr,0,ax,color)

    elif type(a_list) == list:
        if len(a_list) != 1:
            print('Must have phase identifier string of len = 1')
            return
        else:
            for tr in st:
                evdp = tr.stats.evdp
                gcarc = tr.stats.gcarc
                P = model.get_travel_times(distance_in_degree=gcarc,
                    source_depth_in_km=evdp,
                    phase_list = a_list)
                P_time = P[0].time
                if rainbow == True:
                    plot(tr,-1*(P_time+tr.stats.o),ax,randcolor())
                else:
                    plot(tr,-1*(P_time+tr.stats.o),ax,color)

    if picker == True:
        remove_list = []
        def on_pick(event):
            artist = event.artist
            artist.set_c('white')
            artist.set_alpha(0.0)
            remove_list.append(artist.get_label())
            fig.canvas.draw()
        fig.canvas.mpl_connect('pick_event', on_pick)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        plt.show()
        for tr in st:
            if tr.stats.network+'.'+tr.stats.station in remove_list:
                st.remove(tr)

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    if save == False:
        plt.show()
    else:
        plt.savefig(save)
    plt.ion()

def az_section(st,**kwargs):
    '''
    Simpler section plotter for obspy stream object
    '''
    plt.ioff()
    a_list = kwargs.get('a_list',True)
    fig = kwargs.get('fig',None)
    ax = kwargs.get('ax',None)
    color = kwargs.get('color','k')
    save = kwargs.get('save',False)
    picker = kwargs.get('picker',False)
    rainbow = kwargs.get('rainbow',False)
    mode = kwargs.get('mode','gcarc')
    xlim = kwargs.get('x_lim',None)
    ylim = kwargs.get('y_lim',None)
    title = kwargs.get('title','')
    model = kwargs.get('model','prem')
    mult = kwargs.get('mult',1.0)
    model = TauPyModel(model=model)

    if fig == None and ax == None:
        fig,ax = plt.subplots(figsize=(9,12))
        plt.tight_layout()
    else:
        print('using outside figure')
    ax.set_xlabel('Time (s)')
    ax.set_ylabel(mode)
    ax.set_title(title)

    def plot(tr,o,ax,color,**kwargs):
        alpha = kwargs.get('alpha',0.5)
        e = tr.stats.npts/tr.stats.sampling_rate
        t = np.linspace(o,o+e,num=tr.stats.npts)

        ax.plot(t,(mult*tr.data/tr.data.max())+tr.stats.az,alpha=alpha,
                color=color,label=tr.stats.network+'.'+tr.stats.station,
                picker=10,lw=0.8)

    def randcolor():
        c_list = ['#1f77b4','#ff7f0e','#2ca02c',
                   '#d62728','#9467bd','#8c564b',
                   '#e377c2','#7f7f7f','#bcbd22','#17becf']
        return c_list[np.random.randint(len(c_list))]

    if a_list == True:
        for tr in st:
            if rainbow == True:
                plot(tr,0,ax,randcolor(),alpha=0.7)
            else:
                plot(tr,0,ax,color)

    elif type(a_list) == list:
        if len(a_list) != 1:
            print('Must have phase identifier string of len = 1')
            return
        else:
            for tr in st:
                evdp = tr.stats.evdp
                gcarc = tr.stats.gcarc
                P = model.get_travel_times(distance_in_degree=gcarc,
                    source_depth_in_km=evdp,
                    phase_list = a_list)
                P_time = P[0].time
                if rainbow == True:
                    plot(tr,-1*(P_time+tr.stats.o),ax,randcolor())
                else:
                    plot(tr,-1*(P_time+tr.stats.o),ax,color)

    if picker == True:
        remove_list = []
        def on_pick(event):
            artist = event.artist
            artist.set_c('white')
            artist.set_alpha(0.0)
            remove_list.append(artist.get_label())
            fig.canvas.draw()
        fig.canvas.mpl_connect('pick_event', on_pick)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        plt.show()
        for tr in st:
            if tr.stats.network+'.'+tr.stats.station in remove_list:
                st.remove(tr)

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    if save == False:
        plt.show()
    else:
        plt.savefig(save)
    plt.ion()

def fft(tr,**kwargs):
    '''
    plot fast fourier transform of trace object
    '''
    freqmin = kwargs.get('freqmin',0.0)
    freqmax = kwargs.get('freqmax',2.0)
    plot = kwargs.get('plot',True)

    Fs = tr.stats.sampling_rate  # sampling rate
    Ts = tr.stats.delta # sampling interval
    time = tr.stats.endtime - tr.stats.starttime
    t = np.arange(0,time,Ts) # time vector

    n = len(tr.data) # length of the signal
    k = np.arange(n)
    T = n/Fs
    frq = k/T # two sides frequency range
    frq = frq[range(n/2)] # one side frequency range

    Y = np.fft.fft(tr.data)/n # fft computing and normalization
    Y = Y[range(n/2)]
    if tr.data.shape[0] - t.shape[0] == 1:
        t = np.hstack((t,0))

    if plot == False:
        return frq, abs(Y)
    else:

        fig, ax = plt.subplots(2, 1,figsize=(15,8))
        ax[0].plot(t,tr.data,'k')
        ax[0].set_xlabel('Time')
        ax[0].set_ylabel('Amplitude')
        ax[0].grid()
        ax[1].plot(frq,abs(Y),'r') # plotting the spectrum
        ax[1].set_xlim([freqmin,freqmax])
        ax[1].set_xticks(np.arange(freqmin,freqmax,0.25))
        ax[1].set_xlabel('Freq (Hz)')
        ax[1].set_ylabel('|Y(freq)|')
        ax[1].grid()
        plt.show()


