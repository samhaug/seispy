#!/usr/bin/env python

import numpy as np
from matplotlib import pyplot as plt
from obspy.taup import TauPyModel
model = TauPyModel(model="prem")
import seispy
import scipy


def signal_extract(st,phase,**kwargs):
    window_tuple = kwargs.get('window_tuple',(-10,10))

    data_list = []

    for tr in st:
        a = seispy.data.phase_window(tr,phase,window_tuple)
        data_list.append(a.data)
    return np.array(data_list)

def make_corr_list(data_array):

    def correlate_max(data_1,data_2):
        corr = scipy.signal.correlate(data_1,data_2,mode='same')
        shift = np.argmax(corr)-int(data_array.shape[1]/2.)
        return shift

    a,b = np.triu_indices(data_array.shape[0],1)
    delta_t = np.zeros((data_array.shape[0],data_array.shape[0]))

    for ii in range(0,a.shape[0]):
        #t_cor.append(correlate_max(data_array[a[ii],:],data_array[b[ii],:]))
        delta_t[a[ii],b[ii]] = correlate_max(data_array[a[ii],:],data_array[b[ii],:])
    #delta_t = delta_t[:-1,:]
    #delta_t = np.vstack((np.zeros(delta_t.shape[1]),delta_t))

    return delta_t


def corr_align(delta_array,data_array):

    t_shift = []
    for i in range(1,data_array.shape[0]):
        one_sum = 0
        for j in range(1,i):
            one_sum += delta_array[j,i]
        two_sum = 0
        for j in range(i+1,data_array.shape[0]):
            two_sum += delta_array[i,j]
        t_shift.append(int(1./data_array.shape[0]*(two_sum-one_sum)))
    t_shift.append(0)
    print t_shift

    align_data = np.zeros(data_array.shape)

    for ii in range(0,len(t_shift)):
        align_data[ii,:] = np.roll(data_array[ii,:],-1*t_shift[ii])

    return align_data

def PCA(data_array,**kwargs):
    plot = kwargs.get('plot',False)
    limit = kwargs.get('limit',5)
    window_tuple = kwargs.get('window_tuple',(-10,10))
    phase = kwargs.get('phase','P')

    data_cov = np.cov(data_array)
    val,vec = np.linalg.eig(data_cov)

    fig = plt.figure(figsize=(20,10))
    fig.suptitle('PCA for phase {}'.format(phase),fontsize=14, fontweight='bold')

    ax1 = plt.subplot(321)
    ax1.grid()
    ax1.set_ylabel('First Eigenvector')
    ax2 = plt.subplot(323)
    ax2.grid()
    ax2.set_ylabel('First two Eigenvectors')
    ax3 = plt.subplot(325)
    ax3.grid()
    ax3.set_xlabel('Seconds after {}'.format(phase))
    ax3.set_ylabel('First {} Eigenvectors'.format(str(limit)))
    ax_eig = plt.subplot(122)
    ax_eig.set_xlabel('Eigenvalues')
    ax_eig.set_ylim(0,val.max()+0.001)
    ax_eig.grid()

    ax_eig.scatter(range(len(val)),val,marker='D',color='k')
    ax_eig.set_xlim(-1,len(val))
    #ax_eig.set_ylim(0,len(val))

    time = np.linspace(-10,10,num=data_array.shape[1])
    cut_vec = vec[:,0]
    data_project = np.dot(np.transpose(cut_vec),data_array)
    ax1.plot(time,-1*data_project,color='k')


    cut_vec = vec[:,0:2]
    data_project = np.dot(np.transpose(cut_vec),data_array)
    simple_data = np.sum(data_project,axis=0)
    ax2.plot(time,-1*simple_data,color='k')

    cut_vec = vec[:,0:limit]
    data_project = np.dot(np.transpose(cut_vec),data_array)
    simple_data = np.sum(data_project,axis=0)
    ax3.plot(time,-1*simple_data,color='k')

    plt.show()


