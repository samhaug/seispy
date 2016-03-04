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

def bootleg_compute(st,**kwargs):
    repeat = kwargs.get('repeat',100)
    i_resamp = kwargs.get('resamp',len(st))

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
        vesp = plot.vespagram(st,plot=False)
        vesp_list.append(vesp)
    vesp_array = np.zeros((len(vesp_list),vesp_list[0].shape[0],vesp_list[0].shape[1]))
    for ii in range(len(vesp_list)):
        vesp_array[ii,:,:] = vesp_list[ii]
    std_array = np.std(vesp_array,axis=0)
    return std_array

def bootleg_plot(std,**kwargs):
    window_tuple = kwargs.get('window_tuple',(-10,230))
    window_slowness = kwargs.get('window_slowness',(-1.5,1.5))
    slowness_tick = kwargs.get('slowness_tick',-0.1)
    plot_line = kwargs.get('plot_line',False)
    save = kwargs.get('save',False)
    clim = kwargs.get('clim',(-18,-14))
    cmap = kwargs.get('cmap','gnuplot')
    font = kwargs.get('font',paper_font)
    plot = kwargs.get('plot',True)

    std_array = std+1e-18

    fig, ax = plt.subplots(figsize=(19,6))
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











