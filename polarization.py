import numpy as np
import numpy.linalg
import obspy
import seispy
import cmocean
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm


def polarization(tr_vert,tr_rad,**kwargs):
    time_len = kwargs.get('window_len',10)
    time_shift = kwargs.get('time_shift',2)
    plot = kwargs.get('plot',False)

    def window_shift(tr_vert,tr_rad,time_len,time_shift):
        s_len = tr_vert.stats.npts
        s_rate = tr_vert.stats.sampling_rate
        window_length = int(time_len*s_rate)
        index_shift = int(time_shift*s_rate)
        parallel_pol_vec = []
        for ii in range(0,s_len,index_shift):
            window = range(ii-window_length,ii+window_length)
            v_win = tr_vert.data.take(window,mode='wrap')
            r_win = tr_rad.data.take(window,mode='wrap')
            parallel_pol_vec.append((r_win,v_win,ii/s_rate))
        return parallel_pol_vec

    def eigen_vector((rad,vert,time)):
        cov = np.cov(rad,vert)
        eigval, eigvec = numpy.linalg.eig(cov)
        eigval_1 = eigval[0]
        eigval_2 = eigval[1]
        eigvec_1 = eigvec[:,0]
        eigvec_2 = eigvec[:,1]
        linearity = 1 - (min(eigval_1,eigval_2)/max(eigval_1,eigval_2))
        if max(eigval_1,eigval_2) == eigval_1:
            theta = np.degrees(np.arctan2(eigvec_1[0],eigvec_1[1]))
        elif max(eigval_1,eigval_2) == eigval_2:
            theta = np.degrees(np.arctan2(eigvec_2[0],eigvec_2[1]))
        return (linearity, theta, time)

    in_list = window_shift(tr_vert,tr_rad,time_len,time_shift)
    out = []
    for idx,ii in enumerate(in_list):
        out.append(eigen_vector(ii))
    out = np.array(out)
    linearity = out[:,0]
    t = out[:,1]
    time = out[:,2]

    t[t<0] += 180
    #t += -180
    #t[t>-270][t[t>-270]<-180] += 180
    #t[t>90][t[t>90]<180] += -180
    #t[t<0] += 360
    #t[t>90][t[t>90]<180] += 180
    #t[t>180][t[t>180]<270] += -180
    #t[t > 270] = 360 - t[t > 270]
    if plot != True:
        return linearity, t, time
    elif plot == True:
        plot_polarization(tr_vert,tr_rad,linearity,t,time)

def plot_polarization(tr_vert,tr_rad,lin,theta,time):
    fig, ax = plt.subplots(3,1,sharex=True,figsize=(23,12))
    time = time+tr_vert.stats.sac['o']*-1

    points = np.array([time, lin]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    lc = LineCollection(segments,cmap=plt.get_cmap('Spectral'))
    #lc = LineCollection(segments,cmap=cmocean.cm.phase)
    lc.set_array(theta)
    lc.set_linewidth(4)
    ax[0].add_collection(lc)
    cax = fig.add_axes([0.91,0.66470,0.01,0.2353])
    axcb = plt.colorbar(lc,cax=cax)
    axcb.set_label('Angle from vertical')
    ax[0].set_xlim(time[0],time[-1])
    ax[0].set_ylim(-0.1,1.1)
    ax[0].set_ylabel('Linearity')
    timevec = np.linspace(-1*tr_vert.stats.sac['o'],
           (tr_vert.stats.delta*tr_vert.stats.npts)-tr_vert.stats.sac['o'],
           num=tr_vert.stats.npts)
    ax[1].plot(timevec,tr_vert.data,color='k')
    ax[1].text(1, 1, tr_vert.stats.channel,
                  horizontalalignment='right',
                  verticalalignment='bottom',
                  transform=ax[1].transAxes,
                  size='x-large')
    ax[2].plot(timevec,tr_rad.data,color='k')
    ax[2].text(1, 1, tr_rad.stats.channel,
                  horizontalalignment='right',
                  verticalalignment='bottom',
                  transform=ax[2].transAxes,
                  size='x-large')
    ax[2].set_xlabel('Seconds after event')
    ax[0].set_title('Network: {}, Station: {},\
 Dist (deg): {}, Depth (km): {} \nStart: {} \nEnd: {}'.format(
                  tr_vert.stats.network,
                  tr_vert.stats.station,
                  round(tr_vert.stats.sac['gcarc'],3),
                  round(tr_vert.stats.sac['evdp'],3),
                  tr_vert.stats.starttime,
                  tr_vert.stats.endtime))
    for ii in ax:
        ii.set_xticks(range(int(time[0]),int(time[-1]),50))
        ii.grid()
        ii.ticklabel_format(style='sci',scilimits=(0,0),axis='y')

    plt.show()



















