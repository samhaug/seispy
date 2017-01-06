#!/usr/bin/env python
import numpy as np
from matplotlib import pyplot as plt
import plot
import obspy
import multiprocessing as mp

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

def bootstrap_compute(st,**kwargs):
    repeat = kwargs.get('repeat',100)
    i_resamp = kwargs.get('resamp',len(st))
    p_lim = kwargs.get('p_lim',(-1.5,1.5))
    x_lim = kwargs.get('x_lim',(-10,160))
    p_tick = kwargs.get('p_tick',-0.1)
    vesp_row = kwargs.get('vesp_row',1)

    def random_gen(i_resamp):
        rand_list = []
        for ii in range(i_resamp):
            rand_list.append(np.random.randint(0,i_resamp))
        return rand_list

    def resample_stream(st,rand_list):
        st_resamp = obspy.core.Stream()
        for ii in rand_list:
            st_resamp.append(st[ii])
        return st_resamp

    vesp_list = []
    for ii in range(repeat):
        rand_list = random_gen(len(st))
        st_resamp = resample_stream(st,rand_list)
        vesp = plot.vespagram(st_resamp,p_lim=p_lim,
                              x_lim=x_lim,p_tick=p_tick,plot=False)
        vesp_list.append(vesp)

    std_array = np.std(vesp_list,axis=0)
    return std_array 

def bootstrap_color(std_array,**kwargs):
    window_tuple = kwargs.get('window_tuple',(-10,230))
    window_slowness = kwargs.get('window_slowness',(-1.5,1.5))
    slowness_tick = kwargs.get('slowness_tick',-0.1)
    plot_line = kwargs.get('plot_line',False)
    save = kwargs.get('save',False)
    clim = kwargs.get('clim',(np.log10(std_array).min(),
                              np.log10(std_array).max()))
    cmap = kwargs.get('cmap','gnuplot')
    font = kwargs.get('font',paper_font)
    plot = kwargs.get('plot',True)


    fig, ax = plt.subplots(figsize=(15,5))
    image = ax.imshow(np.log10(std_array),aspect='auto',interpolation='lanczos',
                      extent=[window_tuple[0],window_tuple[1],
                      window_slowness[0],window_slowness[1]],cmap=cmap,
                      vmin=clim[0],vmax=clim[1])
    cbar = fig.colorbar(image,ax=ax)
    cbar.set_label('Log(Standard Deviation)',fontdict=font)
    ax.set_xlim(window_tuple)
    ax.xaxis.set(ticks=range(window_tuple[0],window_tuple[1],10))
    ax.set_ylim(window_slowness)
    ax.grid(color='w',lw=2,alpha=0.6)
    ax.set_xlabel('Seconds after P',fontdict=font)
    ax.set_ylabel('Slowness (s/deg)',fontdict=font)
    plt.show()

def bootstrap_wave(vesp,std_array,**kwargs):
    window_tuple = kwargs.get('window_tuple',(-10,230))
    window_slowness = kwargs.get('window_slowness',(-1.5,1.5))
    slowness_tick = kwargs.get('slowness_tick',-0.1)
    plot_line = kwargs.get('plot_line',False)
    save = kwargs.get('save',False)
    clim = kwargs.get('clim',(np.log10(std_array).min(),
                              np.log10(std_array).max()))
    cmap = kwargs.get('cmap','gnuplot')
    font = kwargs.get('font',paper_font)
    plot = kwargs.get('plot',True)

    vesp_wave = vesp.copy()
    fig, ax = plt.subplots(figsize=(15,5))
    ax.set_xlim(window_tuple)
    ax.set_ylim((window_slowness[0],window_slowness[1]+0.1))
    ax.set_xlabel('Seconds after P',fontdict=font)
    ax.set_ylabel('Slowness (s/deg)',fontdict=font)
    ax.xaxis.set(ticks=range(window_tuple[0],window_tuple[1],10))

    time_vec = np.linspace(window_tuple[0],window_tuple[1],
                           num=vesp_wave.shape[1])
    for idx, ii in enumerate(np.arange(window_slowness[1],window_slowness[0],
                                      slowness_tick)):
        vesp_wave[idx,:] +=ii
        std_upper = vesp_wave[idx,:]+std_array[idx,:]
        std_lower = vesp_wave[idx,:]-std_array[idx,:]
        ax.plot(time_vec,vesp_wave[idx,:],color='k',lw=0.5)
        ax.fill_between(time_vec,std_lower,std_upper,alpha=0.5,color='r')

    fig2, ax2 = plt.subplots(figsize=(15,5))
    ax2.set_xlim(window_tuple)
    ax2.set_ylim((window_slowness[0],window_slowness[1]+0.1))
    ax2.set_xlabel('Seconds after P',fontdict=font)
    ax2.set_ylabel('Stacked Amplitude',fontdict=font)
    ax2.xaxis.set(ticks=range(window_tuple[0],window_tuple[1],10))

    vesp_single = vesp_wave[15,:]/np.abs(vesp_wave[15.:]).max()
    std_upper = vesp_single+std_array[15,:]/np.abs(vesp_wave[15,:]).max()
    std_lower = vesp_single-std_array[15,:]/np.abs(vesp_wave[15,:]).max()

    ax2.plot(time_vec,vesp_single,color='k',lw=0.5)
    ax2.fill_between(time_vec,std_lower,std_upper,alpha=0.5,color='r')

    plt.show()









