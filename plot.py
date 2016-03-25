#!/usr/bin/env python

'''
seis_plot.py includes all functions for plotting seismic data
'''

import matplotlib
import scipy
import obspy
import h5py
import numpy as np

from seispy.data import phase_window
from matplotlib import pyplot as plt
from obspy.taup import TauPyModel
import obspy.signal.filter
import obspy.signal
model = TauPyModel(model="prem_50")
from matplotlib import colors, ticker, cm
from matplotlib.patches import Polygon
from matplotlib.colors import LogNorm
from matplotlib.colors import ListedColormap
import os
from scipy.optimize import curve_fit
import math
from matplotlib.colors import LightSource
from mpl_toolkits.basemap import Basemap
from cycler import cycler
import multiprocessing
import shutil
from seispy import mapplot
from seispy import data

ppt_font =  {'family' : 'sans-serif',
             'style' : 'normal',
             'variant' : 'normal',
             'weight' : 'bold',
             'size' : 'xx-large'}
paper_font =  {'family' : 'serif',
             'style' : 'normal',
             'variant' : 'normal',
             'weight' : 'medium',
             'size' : 'large'}

def new_precursor_PKIKP(st_in,**kwargs):


    second = kwargs.get('second',True)
    time_window = kwargs.get('time',(-35,20))
    range_window = kwargs.get('range',(130,140))
    band_filter = kwargs.get('filter',None)
    plot = kwargs.get('plot',True)
    differentiate = kwargs.get('diff',False)
    align = kwargs.get('align',True)
    shift = kwargs.get('shift',0)
    interp = kwargs.get('interp','lanczos')

    def main():

        st = st_in.copy()
        st = sort_stream(st)

        if range_window != None:
            for tr in st:
                if tr.stats.gcarc < range_window[0] or \
                tr.stats.gcarc > range_window[1]:
                    st.remove(tr)

        st.normalize()
        if align == True:
            st = data.align_on_phase(st,phase=['PKIKP'])

        if second == True:
            for tr in st:
                if tr.stats.gcarc%0.5 != 0:
                    st.remove(tr)

        if band_filter != None:
            st.filter('bandpass',freqmin=band_filter[0],freqmax=band_filter[1])


        start = st[0].stats.starttime
        end = st[0].stats.endtime
        time = np.linspace(0,end-start,num=st[0].stats.npts)
        event_depth = tr.stats.sac['evdp']
        samp = tr.stats.sampling_rate

        envelope_array = []
        waveform_array = []

        for tr in st:
            print tr.stats.gcarc
            hil = scipy.fftpack.hilbert(tr.data)
            data_envelope = pow((hil**2+tr.data**2),0.5)
            arrivals = model.get_travel_times(source_depth_in_km=event_depth,
                       distance_in_degree=tr.stats.sac['gcarc'],phase_list=['PKIKP'])
            PKIKP = arrivals[0].time+tr.stats.sac['o']+shift
            closest_time = np.argmin(abs(PKIKP-time))
            begin =  closest_time+(int(time_window[0]*samp))
            stop = closest_time+(int(time_window[1]*samp))
            envelope_array.append(data_envelope[begin:stop])
            waveform_array.append(tr.data[begin:stop])

        envelope = (np.array(envelope_array))
        envelope = envelope/envelope.max()
        waveform = ((np.array(waveform_array)))

        if plot == False:
            return envelope, waveform

        plot_env(envelope,range_window,time_window)
        plot_wave(waveform,range_window,time_window)
        plt.show()
        return envelope, waveform


    def sort_stream(st):
        for tr in st:
            tr.stats.gcarc = round(tr.stats.sac['gcarc'],1)
        st.sort(['gcarc'])
        return st

    def plot_env(envelope,rw,tw):
        e = np.log10(envelope)
        fig, ax = plt.subplots()
        image = ax.imshow(np.log10(np.flipud(envelope)+1e-8),aspect='auto',
            cmap='Spectral_r',extent=[(tw[0]),(tw[1]),rw[0],rw[1]],
            interpolation=interp,vmin=-3.0,vmax=e.max())
        ax.set_ylabel('Range (degrees)',fontsize=14)
        ax.set_xlabel(r'Seconds before $PKPdf$',fontsize=14)
        #ax.text(121,-19,'CMB',color='w',fontsize=14)
        cbar = fig.colorbar(image,ax=ax)
        cbar.solids.set_rasterized(True)
        cbar.set_label(r'$\log_{10}(Amplitude)$',fontsize=14)
        ax.grid(color='k')

    def plot_wave(waveform,rw,tw):
        w = waveform/waveform.max()+rw[0]
        t = np.linspace(tw[0],tw[1],num=w.shape[1])
        fig1,ax1 = plt.subplots()
        for ii in range(0,w.shape[0]):
            fill = (rw[0]+ii/2.0)*np.ones(w.shape[1])
            ax1.fill_between(t,w[ii,:]+ii/2.0,fill,where=fill<=w[ii,:]+ii/2.0,
                             facecolor='goldenrod',alpha=0.5)
            ax1.fill_between(t,w[ii,:]+ii/2.0,fill,where=fill>=w[ii,:]+ii/2.0,
                             facecolor='blue',alpha=0.5)
            ax1.plot(t,w[ii,:]+ii/2.0,color='k',alpha=0.5)
            #ax1.plot(t,fill,color='k')
            #ax1.fill_between(t,w[ii,:],fill,where=w[ii,:] <=fill,
            #                 facecolor='blue',alpha=0.5)
        ax1.set_ylabel('Range (degrees)',fontsize=14)
        ax1.set_xlabel(r'Seconds before $PKPdf$',fontsize=14)

    return  main()


def precursor_PKIKP(seis_stream_in,precursor_onset,time_offset,name=False,plot=True):
    '''
    Produce colorbar of PKPdf scatterers from Mancinelli & Shearer 2011
    PARAMETERS
    __________

    seis_stream_in: obspy.core.stream.Stream object.

    precursor_onset = h5py file

    time_offset = integer/float amount of seconds to shift the arrival
                 to match PKPdf on zero

    name = default = False. If a string is passed, it will save the
                 image to pdf in the
           directory work4/IMAGES
    '''
    event_depth = seis_stream_in[0].stats.sac['evdp']

    f = h5py.File(precursor_onset,'r')

    #Find reciever latitude and correct for coordinate system
    def receiver_latitude(tr):
        station = tr.stats.station.split('_')
    #   lat = float(station[0]+'.'+station[1])
        lat = tr.stats.sac['gcarc']
        return lat

    seis_stream = seis_stream_in.copy()
    shearer_stream = obspy.core.stream.Stream()
    envelope_dict = {}
    waveform_dict = {}

    t = seis_stream[0].stats.starttime
    samp_rate = seis_stream[0].stats.sampling_rate

    for tr in seis_stream:
        lat = receiver_latitude(tr)
        if lat in np.arange(123.,143.):
        #if 120. <= lat <= 145.:
            tr.differentiate()
            tr.filter('bandpass',freqmin=1.5,freqmax=2.5)
            arrivals = model.get_travel_times(source_depth_in_km=event_depth,
                       distance_in_degree = lat, phase_list = ['PKiKP'])
            PKiKP = arrivals[0].time+time_offset
            #array = tr.slice(t+PKiKP-20,t+PKiKP+10)
            #array_norm = tr.slice(t+PKiKP-18,t+PKiKP-15).data.max()
            array_norm = 1
            #data_envelope = obspy.signal.filter.envelope(array.data)
            #hil = scipy.fftpack.hilbert(array.data)
            #data_envelope = pow((hil**2+array.data**2),0.5)
            hil = scipy.fftpack.hilbert(tr.data)
            data_envelope = pow((hil**2+tr.data**2),0.5)
            tr.data = data_envelope
            array = tr.slice(t+PKiKP-20,t+PKiKP+10)
            data_envelope = array.data 
            #data_envelope = obspy.signal.filter.envelope(array.data)
            data_waveform = (array.data/array_norm).clip(min=0)*2e-2
            data_waveform = (array.data/array_norm)*2e-1
            envelope_dict[lat] = data_envelope
            waveform_dict[lat] = data_waveform

    print envelope_dict.keys()
    rows = envelope_dict[130].shape[0]
    cols = len(envelope_dict.keys())
    waveform_array = np.zeros((rows,cols))
    envelope_array = np.zeros((rows,cols))

    for idx, keys in enumerate(sorted(envelope_dict)):
        length = envelope_dict[keys][0:rows].size
        envelope_array[0:length,idx] = envelope_dict[keys][0:rows]
        waveform_array[0:length,idx] = waveform_dict[keys][0:rows]

    waveform_array = np.fliplr(np.transpose(np.flipud(waveform_array)))
    precursor_onset = int(waveform_array.shape[1]*3./5.)

    # Adjust precursor amplitude
    waveform_array[:,0:precursor_onset] *= 1
    time = np.linspace(-25,25,num=waveform_array.shape[1])

    if plot == False:
        return envelope_array

    fig, (ax,ax1) = plt.subplots(1,2,figsize=plt.figaspect(0.5))

    # Waveform axis
    for ii in range(0,waveform_array.shape[0])[::-5]:
        waveform_array[ii,:] += (120+(ii/10.))
        base = (120+(ii/10.))*np.ones(waveform_array[ii,:].shape)
        #ax.plot(time,waveform_array[ii,:],c='k',lw=0)
        #ax.plot(time,base,c='w',alpha=0.0)
        ax.fill_between(time, waveform_array[ii,:], base,
                where=waveform_array[ii,:] >= base,
                facecolor='#E6A817', interpolate=True,lw=0.6,alpha=0.8)
        ax.fill_between(time, waveform_array[ii,:], base,
                where=waveform_array[ii,:] <= base,
                facecolor='blue', interpolate=True,lw=0.6,alpha=0.8)

    x = np.linspace(120,145,num=f['2891'][...].shape[0])
    for items in f:
       onset = f[items][...]
       onset[onset > 0] = np.nan
       ax.plot(onset,x,color='k',lw=2)
       ax1.plot(x,onset,color='w',lw=2)

    ax.set_ylim([120,148])
    ax.set_xlim([-25,25])
    ax.set_ylabel('Range (degrees)',fontsize=14)
    ax.set_xlabel(r'Seconds after $PKPdf$',fontsize=14)
    #ax.plot(y,x,color='k',lw=2)

    # Colorbar axis
    envelope_array = envelope_array/envelope_array.max()
    image = ax1.imshow(np.log10(np.flipud(envelope_array)+1e-8),aspect='auto',
            cmap='Spectral_r',extent=[120,145,-25,25],
            interpolation='none',vmin=-3.5,vmax=-2.5)
    ax1.set_xlabel('Range (degrees)',fontsize=14)
    ax1.set_ylabel(r'Seconds after $PKPdf$',fontsize=14)
    ax1.set_ylim(-25,25)
    ax1.set_xlim(120,145)
    ax1.text(121,-19,'CMB',color='w',fontsize=14)
    #ax1.plot(x,y,color='w',lw=2)
    cbar = fig.colorbar(image,ax=ax1)
    cbar.solids.set_rasterized(True)
    cbar.set_label(r'$\log_{10}(Amplitude)$',fontsize=14)

    return np.flipud(envelope_array)
    #if name == False:
    #    plt.show()
    #else:
    #    fig.suptitle(str(name),fontsize=16)
    #    plt.savefig('/home/samhaug/Documents/Figures/'+name+'.pdf')

def vespagram(st_in,**kwargs):
    '''
    make vespagram of stream object. Recommended to pass stream object
    through remove_dirty
    '''

    window_tuple = kwargs.get('window_tuple',(-10,230))
    window_phase = kwargs.get('window_phase',['P'])
    window_slowness = kwargs.get('window_slowness',(-1.5,1.5))
    slowness_tick = kwargs.get('slowness_tick',-0.1)
    phase_list = kwargs.get('phase_list',False)
    plot_line = kwargs.get('plot_line',False)
    n_root = kwargs.get('n_root',2.0)
    save = kwargs.get('save',False)
    clim = kwargs.get('clim',(-3,0))
    cmap = kwargs.get('cmap','gnuplot')
    font = kwargs.get('font',paper_font)
    plot = kwargs.get('plot',True)
    title = kwargs.get('title',True)

    st = obspy.core.Stream()
    for idx, tr in enumerate(st_in):
        st += phase_window(tr,window_phase,window_tuple)

    def mean_range(st):
        a = []
        for tr in st:
            a.append(tr.stats.sac['gcarc'])
        mn_r = np.mean(a)
        return mn_r

    def roll_zero(array,n):
        if n < 0:
            array = np.roll(array,n)
            array[n::] = 0
        else:
            array = np.roll(array,n)
            array[0:n] = 0
        return array

    def slant_stack(st,mn_range,slowness,n):
        d = st[0].stats.delta
        R = np.zeros(st[0].data.shape[0])
        for tr in st:
            az = tr.stats.sac['gcarc']-mn_range
            shift_in_sec = slowness*az
            shift_in_bin = int(shift_in_sec/d)
            x = roll_zero(tr.data,shift_in_bin)
            R += np.sign(x)*pow(np.abs(x),1./n)
        R = R/float(len(st))
        yi = R*pow(abs(R),n-1)
        hil = scipy.fftpack.hilbert(yi)
        yi = pow(hil**2+R**2,1/2.)
        return yi,R

    def phase_plot(ax,evdp,degree,phases,text_color):
        P_arrivals = model.get_travel_times(distance_in_degree = degree,
        source_depth_in_km = evdp,
        phase_list = ['P'])

        P_slowness = P_arrivals[0].ray_param_sec_degree
        P_time = P_arrivals[0].time

        arrivals = model.get_travel_times(distance_in_degree=degree,
        source_depth_in_km=evdp,
        phase_list = phases)
        if len(arrivals) != 0:
            colors = ['b','g','r','c','m','y','k']
            for idx, ii in enumerate(arrivals):
                p = ii.ray_param_sec_degree-P_slowness
                time = ii.time-P_time
                ax.scatter(time,p,s=22,marker='D',c='w',zorder=20)
                #ax.text(time,p,ii.name,fontsize=16,color=text_color)

    mn_r = mean_range(st)
    evdp = st[0].stats.sac['evdp']

    #st.normalize()

    vesp_y = np.linspace(0,0,num=st[0].data.shape[0])
    vesp_R = np.linspace(0,0,num=st[0].data.shape[0])
    for ii in np.arange(window_slowness[0],window_slowness[1],0.1):
        yi,R = slant_stack(st,mn_r,ii,1.0)
        vesp_y= np.vstack((vesp_y,yi))
        vesp_R= np.vstack((vesp_R,R))
    vesp_y = vesp_y[1::,:]
    vesp_R = vesp_R[1::,:]
    if plot == False:
        return vesp_R
    vesp_y = vesp_y/ vesp_y.max()


    fig, ax = plt.subplots(2,sharex=True,figsize=(15,10))

    image_0 = ax[0].imshow(np.log10(vesp_y),aspect='auto',interpolation='lanczos',
           extent=[window_tuple[0],window_tuple[1],window_slowness[0],window_slowness[1]],
           cmap=cmap,vmin=clim[0],vmax=clim[1])
    one_stack = vesp_y

    vesp_y = np.linspace(0,0,num=st[0].data.shape[0])
    vesp_R = np.linspace(0,0,num=st[0].data.shape[0])
    for ii in np.arange(window_slowness[0],window_slowness[1],0.1):
        yi,R = slant_stack(st,mn_r,ii,n_root)
        vesp_y= np.vstack((vesp_y,yi))
        vesp_R= np.vstack((vesp_R,R))
    vesp_y = vesp_y[1::,:]
    vesp_R = vesp_R[1::,:]
    vesp_y = vesp_y/ vesp_y.max()

    image_1 = ax[1].imshow(np.log10(vesp_y), aspect='auto',
         interpolation='lanczos', extent=[window_tuple[0],
         window_tuple[1],window_slowness[0],window_slowness[1]],cmap=cmap, vmin=clim[0],vmax=clim[1])

    two_stack = vesp_y

    cbar_0 = fig.colorbar(image_0,ax=ax[0])
    cbar_0.set_label('Log(Seismic energy)',fontdict=font)
    cbar_1 = fig.colorbar(image_1,ax=ax[1])
    cbar_1.set_label('Log(Seismic energy)',fontdict=font)
    ax[0].set_xlim(window_tuple)
    ax[1].set_xlim(window_tuple)
    ax[1].xaxis.set(ticks=range(window_tuple[0],window_tuple[1],10))
    ax[0].xaxis.set(ticks=range(window_tuple[0],window_tuple[1],10))
    ax[0].set_ylim([window_slowness[0],window_slowness[1]])
    ax[1].set_ylim([window_slowness[0],window_slowness[1]])
    ax[0].grid(color='w',lw=2,alpha=0.9)
    ax[1].grid(color='w',lw=2,alpha=0.9)
    ax[1].set_ylabel('Slowness (s/deg)',fontdict=font)
    ax[0].set_ylabel('Slowness (s/deg)',fontdict=font)
    ax[1].set_xlabel('Seconds after {}'.format(window_phase[0]),fontdict=font)
    if title == True:
        ax[0].set_title('Start: {} \n Source Depth: {} km, Ref_dist: {} deg, {} \
                     \n Bottom : N-root = {} Top: N-root = 1'
                  .format(st[0].stats.starttime,
                  round(st[0].stats.sac['evdp'],3),
                  round(mn_r,3), os.getcwd().split('-')[3],
                  str(n_root)))

    if phase_list:
        phase_plot(ax[0],evdp,mn_r,phase_list,text_color='green')
        phase_plot(ax[1],evdp,mn_r,phase_list,text_color='green')
        #phase_plot(axR,evdp,mn_r,phase_list,text_color='black')

    if save != False:
        plt.savefig(save+'/vespagram.pdf',format='pdf')

    time_vec = np.linspace(window_tuple[0], window_tuple[1],
               num=vesp_R.shape[1])

    figR, axR= plt.subplots(1,figsize=(14,7))

    #vesp_R *= 2
    for idx, ii in enumerate(np.arange(window_slowness[1],window_slowness[0],
                             slowness_tick)):
        vesp_R[idx,:] += ii
        axR.fill_between(time_vec,ii,vesp_R[idx,:],where=vesp_R[idx,:] >= ii,
                         facecolor='goldenrod',alpha=0.8,lw=0.5)
        axR.fill_between(time_vec,ii,vesp_R[idx,:],where=vesp_R[idx,:] <= ii,
                         facecolor='blue',alpha=0.8,lw=0.5)
        #axR.plot(time_vec,vesp_R[idx,:])

    axR.set_xlim(window_tuple)
    axR.set_ylim([window_slowness[0],window_slowness[1]])
    axR.set_ylabel('Slowness (s/deg)')
    axR.set_xlabel('Seconds after P')
    if title == True:
        axR.set_title('Start: {} \n Source Depth: {} km, Ref_dist: {} deg, {}'
                  .format(st[0].stats.starttime,
                   round(st[0].stats.sac['evdp'],3),round(mn_r,3),
                   os.getcwd().split('-')[3]))
    axR.set_xticks(np.arange(window_tuple[0],window_tuple[1],10))
    axR.grid()


    if save != False:
        plt.savefig(save+'/wave.pdf',format='pdf')
    else:
        plt.show()

    return vesp_R
def plot(tr,**kwargs):
    '''
    plot trace object
    phase = optional list of phases. list of strings corresponding
            to taup names
    '''
    phases = kwargs.get('phase_list',False)
    window = kwargs.get('window',False)

    fig, ax = plt.subplots(figsize=(23,9))
    if phases != False:
        arrivals = model.get_travel_times(
                   distance_in_degree=tr.stats.sac['gcarc'],
                   source_depth_in_km=tr.stats.sac['evdp'],
                   phase_list = phases)
        window_list = []
        colors = ['b','g','r','c','m','y','k']
        if len(arrivals) != 0:
            for idx, ii in enumerate(arrivals):
                ax.axvline(ii.time,label=ii.name,c=colors[idx])
                window_list.append(ii.time)
    ax.legend()

    time = np.linspace(-1*tr.stats.sac['o'],
           (tr.stats.delta*tr.stats.npts)-tr.stats.sac['o'],
           num=tr.stats.npts)
    ax.plot(time,tr.data,c='k')
    ax.grid()
    ax.set_title('Network: {}, Station: {}, Channel: {},\
 Dist (deg): {}, Depth (km): {} \nStart: {} \nEnd: {}'.format(
                  tr.stats.network,
                  tr.stats.station,
                  tr.stats.channel,
                  round(tr.stats.sac['gcarc'],3),
                  round(tr.stats.sac['evdp'],3),
                  tr.stats.starttime,
                  tr.stats.endtime))
    ax.set_xlabel('Time (s), sampling_rate: {}, npts: {}'
                 .format(tr.stats.sampling_rate,tr.stats.npts))

    if window != False:
        ax.set_xlim([min(window_list)+window[0],max(window_list)+window[1]])
        #ax.set_xlim([min(window_list)-300,max(window_list)+300])

    plt.show()

def component_plot(tr_list,**kwargs):
    '''
    Plot three component section
    '''
    if tr_list[0].stats.station != tr_list[1].stats.station:
        print 'Components from different stations'
    separate = kwargs.get('separate',True)
    phase_list = kwargs.get('phase_list',[])
    window_tuple = kwargs.get('window_tuple',False)

    arrivals = model.get_travel_times(
               distance_in_degree=tr_list[0].stats.sac['gcarc'],
               source_depth_in_km=tr_list[0].stats.sac['evdp'],
               phase_list = phase_list)
    colors = ['b','g','r','c','m','y','k','b','g','r','c','m','y']
    trace_c = ['k','m','goldenrod']

    if separate:
        fig, ax = plt.subplots(len(tr_list), sharey=True,sharex=True,
                               figsize=(23,15))
        for idx, tr in enumerate(tr_list):
            t_len = tr.stats.endtime - tr.stats.starttime
            start = -1*tr.stats.sac['o']
            end = t_len-tr.stats.sac['o']
            time = np.linspace(start,end,num=tr.stats.npts)
            data = tr.data
            trace = ax[idx].plot(time,data,color='k',zorder=99)
            ax[idx].text(1, 1, tr.stats.channel,
                           horizontalalignment='right',
                           verticalalignment='bottom',
                           transform=ax[idx].transAxes,
                           size='x-large')
            ax[idx].grid()
            ax[idx].legend(loc=2)
            ax[idx].set_xticks(np.arange(time[0],time[-1],25))
            #add arrivals
            for jdx, ii in enumerate(arrivals):
                ax[idx].axvline(ii.time,label=ii.name,c=colors[jdx])
                ax[idx].legend(loc=2)

        t_len = tr_list[0].stats.endtime - tr_list[0].stats.starttime
        start = -1*tr_list[0].stats.sac['o']
        end = t_len-tr_list[0].stats.sac['o']
        if window_tuple:
            ax[0].set_xlim((arrivals[0].time+window_tuple[0],
                        arrivals[0].time+window_tuple[1]))
            ax[1].set_xlim((arrivals[0].time+window_tuple[0],
                        arrivals[0].time+window_tuple[1]))
        ax[-1].set_xlabel('Seconds after event')
        ax[0].set_title('{} \n Depth (km): {} Dist (deg): {}, Az (deg): {}, {}, {}'.format(
              tr_list[0].stats.starttime,
              str(round(tr_list[0].stats.sac['evdp'],3)),
              str(round(tr_list[0].stats.sac['gcarc'],3)),
              str(round(tr_list[0].stats.sac['az'],3)),
              tr_list[0].stats.network+tr_list[0].stats.station,
              os.getcwd().split('-')[3]))
        for ii in ax:
            ii.ticklabel_format(style='sci',scilimits=(0,0),axis='y')
        plt.show()



    elif separate != True:
        fig, ax = plt.subplots(figsize=(23,9))
        for idx, tr in enumerate(tr_list):
            t_len = tr.stats.endtime - tr.stats.starttime
            start = -1*tr.stats.sac['o']
            end = t_len-tr.stats.sac['o']
            time = np.linspace(start,end,num=tr.stats.npts)
            data = tr.data
            trace = ax.plot(time,data,label=tr.stats.channel,color=trace_c[idx])
            ax.legend(loc=2)
            #add arrivals
        for jdx, ii in enumerate(arrivals):
            ax.axvline(ii.time,label=ii.name,c=colors[jdx])
            ax.legend(loc=2)

        t_len = tr[0].stats.endtime - tr[0].stats.starttime
        start = -1*tr[0].stats.sac['o']
        end = t_len-tr[0].stats.sac['o']
        if window_tuple:
            ax.set_xlim((start+arrivals[0].time+window_tuple[0],
                        start+arrivals[0].time+window_tuple[1]))
        ax.grid()
        ax.set_xlabel('Seconds after event')
        ax.set_xticks(np.arange(start-50,end,25))
        ax.set_title('{} \n Depth (km): {} Dist (deg): {}, Az (deg): {}, {}'.format(
            tr_list[0].stats.starttime,
            str(round(tr_list[0].stats.sac['evdp'],3)),
            str(round(tr_list[0].stats.sac['gcarc'],3)),
            str(round(tr_list[0].stats.sac['baz'],3)),
            os.getcwd().split('-')[3]))
        plt.show()

def parallel_add_to_axes(trace_tuple):
    '''
    parallel plotting function to be used by section
    '''
    data = trace_tuple[0]
    time = trace_tuple[1]
    dist = trace_tuple[2]
    line = matplotlib.lines.Line2D(time,data+dist,alpha=0.5,c='k',lw=1)
    return line

def section(st,**kwargs):
    '''
    Plot record section of obspy stream object

    labels = kwargs.get('labels',False)
    phases = kwargs.get('phase_list',False)
    fill = kwargs.get('fill',False)
    shift = kwargs.get('shift',False)
    save = kwargs.get('save',False)
    title = kwargs.get('title',True)
    x_lim = kwargs.get('x_lim',(-50,1000))
    color = kwargs.get('color',False)
    picker = kwargs.get('picker',False)
    '''
    labels = kwargs.get('labels',False)
    phases = kwargs.get('phase_list',False)
    fill = kwargs.get('fill',False)
    shift = kwargs.get('shift',False)
    save = kwargs.get('save',False)
    title = kwargs.get('title',True)
    x_lim = kwargs.get('x_lim',(-50,1000))
    color = kwargs.get('color',False)
    picker = kwargs.get('picker',False)
    align_phase = kwargs.get('align_phase',['P'])
    name = kwargs.get('name_plot',False)

    def main():
        p_list,name_list,dist_list = p_list_maker(st)
        lim_tuple = ax_limits(p_list)

        fig, ax = plt.subplots(figsize =(10,15))
        ax.set_xlim((x_lim[0],x_lim[1]))
        for idx,ii in enumerate(p_list):
            add_to_axes(ii,ax)

        if phases != False:
            range_list = []
            for tr in st:
                tr.stats.location = tr.stats.sac['gcarc']
            st.sort(['location'])
            ymin = st[0].stats.location
            ymax = st[-1].stats.location
            t = st[len(st)/2]
            phase_plot(lim_tuple,70.,t.stats.sac['evdp'],phases,ax,
                       t.stats.sac['o'])
            ax.set_ylim((ymin-2,ymax+2))

        ax.set_ylabel('Distance (deg)')
        ax.set_xlabel('Seconds After Event')

        if shift:
            ax.set_xlabel('Seconds')
        if title == True:
            ax.set_title("""Start: {} Coord: ({}$^\circ$,{}$^\circ$)
Depth (km): {}, Channel: {}, Mag: {}""".format(
             st[0].stats.starttime.isoformat(),
             round(st[0].stats.sac['evla'],2),
             round(st[0].stats.sac['evlo'],2),
             round(st[0].stats.sac['evdp'],3),
             st[0].stats.channel,
             os.getcwd().split('-')[3][-2]+'.'+os.getcwd().split('-')[3][-1]))
        if labels:
            y1, y2 = ax.get_ylim()
            ax_n = ax.twinx()
            ax_n.set_yticks(dist_list)
            ax_n.set_yticklabels(name_list)
            ax_n.set_ylim(y1,y2)
            for ii in dist_list:
                ax.axhline(ii,alpha=0.5,c='k')

        if picker == True:
            rmfile = file('./REMOVE_LIST','w')
            remove_list = []
            def on_pick(event):
                artist = event.artist
                artist.set_c('white')
                artist.set_alpha(0.0)
                remove_list.append(artist.get_label())
                fig.canvas.draw()
            fig.canvas.mpl_connect('pick_event', on_pick)
            plt.show()
            for tr in st:
                if tr.stats.network+'.'+tr.stats.station in remove_list:
                    st.remove(tr)
            for item in remove_list:
                rmfile.write("%s\n" % item)
            rmfile.close()

        if save != False:
            plt.savefig(save+'/section.pdf',format='pdf')
        if save == False:
            plt.show()

    def phase_plot(lim_tuple,ref_degree,evdp,phases,ax,o):
        arrivals = model.get_travel_times(distance_in_degree=ref_degree,
                   source_depth_in_km=evdp,
                   phase_list = phases)
        P = model.get_travel_times(distance_in_degree=ref_degree,
                   source_depth_in_km=evdp,
                   phase_list = phases)
        P_slow = P[0].ray_param_sec_degree
        P_time = P[0].time
        ax.axvline(0,c='b',alpha=0.2,lw=2.0)
        if len(arrivals) != 0:
            colors = ['b','g','r','c','m','y','b','g','r','c']
            for idx, ii in enumerate(arrivals):
                p =  ii.ray_param_sec_degree - P_slow
                time = ii.time-P_time
                x = np.linspace(time-500,time+500)
                y = (1/p)*(x-time)+ref_degree
                ax.plot(x,y,alpha=0.2,label=ii.purist_name,c=colors[idx],lw=2.0)
                ax.legend()

    def add_to_axes(trace_tuple,ax):
        data = trace_tuple[0]
        time = trace_tuple[1]
        dist = trace_tuple[2]
        name = trace_tuple[3]
        if color == True and picker == True:
            ax.plot(time,data+dist,alpha=0.5,lw=1,picker=5,label=name)
        if color == True and picker != True:
            ax.plot(time,data+dist,alpha=0.5,lw=1)
        if color != True and picker == True:
            ax.plot(time,data+dist,alpha=0.5,lw=1,picker=5,label=name)
        if color != True and picker != True:
            ax.plot(time,data+dist,alpha=0.5,lw=1,c='k')
        if fill:
            ax.fill_between(time, dist, data+dist, where=data+dist <= dist,
                            facecolor='r', alpha=0.5, interpolate=True)
            ax.fill_between(time, dist, data+dist, where=data+dist >= dist,
                            facecolor='g', alpha=0.5, interpolate=True)

    def p_list_maker(st):
        p_list = []
        name_list = []
        dist_list = []
        for tr in st:
            o = tr.stats.sac['o']
            data = tr.data
            start = tr.stats.starttime+o
            end = tr.stats.endtime+o
            if align_phase == None:
                p_time = o
            else:
                arrivals = model.get_travel_times(
                   distance_in_degree=tr.stats.sac['gcarc'],
                   source_depth_in_km=tr.stats.sac['evdp'],
                   phase_list = align_phase)
                p = arrivals[0]
                p_time = p.time+o
            time = np.linspace(-1*p_time,end-start-p_time,num=tr.stats.npts)
            name = (str(tr.stats.network)+'.'+str(tr.stats.station))
            name_list.append(name)
            dist_list.append(tr.stats.sac['gcarc']+tr.data[0])
            p_list.append((data,time,tr.stats.sac['gcarc'],name))
        return p_list,name_list,dist_list

    def ax_limits(p_list):
        range_list = []
        for ii in p_list:
            range_list.append([ii[2],ii[1][0],ii[1][-1]])
        range_list = np.array(range_list)
        min_range = min(range_list[:,0])
        max_range = max(range_list[:,0])
        min_time = min(range_list[:,1])
        max_time = max(range_list[:,2])
        return ([min_time,max_time],[min_range-3,max_range+3])

    main()

def az_section(st, phase, **kwargs):
    '''
    Plot traces as a function of azimuthal change
    '''
    window_tuple = kwargs.get('window_tuple',(-40,40))
    tr_list = []
    for tr in st:
        d = data.phase_window(tr,['S'],window_tuple).normalize().data
        az = round(tr.stats.sac['az'],3)
        tr_list.append((az,d))


    fig, ax = plt.subplots(figsize=(10,15))
    ax.set_ylabel('Source-reciever azimuth')
    ax.set_xlabel('Seconds after PREM predicted {} phase'.format(phase[0]))
    ax.set_title('{} \n Channel: {}, Depth: {} km'.format(
                 st[10].stats.starttime,
                 st[10].stats.channel,
                 round(st[10].stats.sac['evdp'],3)))
    ax.grid()

    for ii in tr_list:
        time = np.linspace(window_tuple[0],window_tuple[1],
                           num=len(ii[1]))
        ax.plot(time,ii[0]+ii[1],'k',alpha=0.5)
    plt.show()

def fft(tr, freqmin=0.0, freqmax=2.0):
    '''
    plot fast fourier transform of trace object
    '''

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

def express_plot(st, **kwargs):
    phase_list = kwargs.pop('phase_list',['P'])
    name = os.getcwd().split('/')[-1]
    fig_dir = '/home/samhaug/work1/Figures/'
    if os.path.exists(fig_dir+name):
        #os.rmdir(fig_dir+name)
        shutil.rmtree(fig_dir+name)
    os.mkdir(fig_dir+name)

    #vespagram(st,phase_list=phase_list,save=fig_dir+name)
    vespagram(st,phase_list=['P','pP','sP','S660P','S860P','S1260P','S1060P',
                          'S1460P','S1660P'],window_tuple=(-10,200),
                           save=fig_dir+name)
    section(st,shift=True,save=fig_dir+name)
    mapplot.source_reciever_plot(st,save=fig_dir+name)
    st.write(fig_dir+name+'/'+name+'.SAC',format='SAC')














