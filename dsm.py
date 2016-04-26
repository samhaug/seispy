#!/usr/bin/env python

import numpy as np
import scipy
from matplotlib import pyplot as plt
import subprocess as sp
from seispy import convert

def read_file(fname):

    f = open(fname)
    flist = f.read().strip().split('\n')
    flist = flist[16:-11]
    inc = []
    for ii in np.arange(0,len(flist),6):
        inc.append(flist[ii:ii+6])
    return inc

def grab_layer(f_layer):

    layer_dict = {}

    l1 = f_layer[0].strip().split()
    layer_dict['rmin'] = float(l1[0])
    layer_dict['rmax'] = float(l1[1])

    layer_dict['rhod'] = float(l1[2])
    layer_dict['rhoc'] = float(l1[3])
    layer_dict['rhob'] = float(l1[4])
    layer_dict['rhoa'] = float(l1[5])

    lph = f_layer[1].strip().split()
    layer_dict['phd'] = float(lph[0])
    layer_dict['phc'] = float(lph[1])
    layer_dict['phb'] = float(lph[2])
    layer_dict['pha'] = float(lph[3])

    lpv = f_layer[2].strip().split()
    layer_dict['pvd'] = float(lpv[0])
    layer_dict['pvc'] = float(lpv[1])
    layer_dict['pvb'] = float(lpv[2])
    layer_dict['pva'] = float(lpv[3])

    lsh = f_layer[3].strip().split()
    layer_dict['shd'] = float(lsh[0])
    layer_dict['shc'] = float(lsh[1])
    layer_dict['shb'] = float(lsh[2])
    layer_dict['sha'] = float(lsh[3])

    lsv = f_layer[4].strip().split()
    layer_dict['svd'] = float(lsv[0])
    layer_dict['svc'] = float(lsv[1])
    layer_dict['svb'] = float(lsv[2])
    layer_dict['sva'] = float(lsv[3])

    n = f_layer[5].strip().split()
    layer_dict['nd'] = float(n[0])
    layer_dict['nc'] = float(n[1])
    layer_dict['nb'] = float(n[2])
    layer_dict['na'] = float(n[3])
    layer_dict['qa'] = float(n[4])
    layer_dict['qk'] = float(n[5])

    return layer_dict

def make_layer(inc):
    val_list = []
    for ii in range(0,len(inc)):
        l = grab_layer(inc[ii])
        lrad = l['rmax']-l['rmin']
        r = np.linspace(l['rmin'],l['rmax'])
        a = 6371.
        pv = l['pva']*pow(r/a,3)+l['pvb']*pow(r/a,2)+l['pvc']*(r/a)+l['pvd']
        sv = l['sva']*pow(r/a,3)+l['svb']*pow(r/a,2)+l['svc']*(r/a)+l['svd']
        rho = l['rhoa']*pow(r/a,3)+l['rhob']*pow(r/a,2)+l['rhoc']*(r/a)+l['rhod']
        val_list.append([r,pv,sv,rho])
    return val_list

def make_discont(val_list,thick,amp):
    left_rad = val_list[0][0]
    cent_rad = val_list[1][0]
    right_rad = val_list[2][0]

    extend = (cent_rad[-1]-cent_rad[0]-thick)/2.
    left_rad_end = left_rad[-1]+extend
    right_rad_beg = right_rad[0]-extend


    left_r = np.linspace(left_rad[0],left_rad_end)
    center_r = np.linspace(left_rad_end,right_rad_beg)
    right_r = np.linspace(right_rad_beg,right_rad[-1])

    p_amp_0 =     val_list[0][1][0]
    p_new_amp_1 = val_list[0][1][-1]*(1+amp)
    p_new_amp_2 = val_list[2][1][0]*(1-amp)
    p_amp_2_end = val_list[2][1][-1]

    s_amp_0 =     val_list[0][2][0]
    s_new_amp_1 = val_list[0][2][-1]*(1+amp)
    s_new_amp_2 = val_list[2][2][0]*(1-amp)
    s_amp_2_end = val_list[2][2][-1]

    rho_amp_0 =     val_list[0][3][0]
    rho_new_amp_1 = val_list[0][3][-1]*(1+amp)
    rho_new_amp_2 = val_list[2][3][0]*(1-amp)
    rho_amp_2_end = val_list[2][3][-1]

    p_left_amp = np.linspace(p_amp_0,p_new_amp_1)
    p_cent_amp = np.linspace(p_new_amp_1,p_new_amp_2)
    p_right_amp = np.linspace(p_new_amp_2,p_amp_2_end)

    s_left_amp = np.linspace(s_amp_0,s_new_amp_1)
    s_cent_amp = np.linspace(s_new_amp_1,s_new_amp_2)
    s_right_amp = np.linspace(s_new_amp_2,s_amp_2_end)

    rho_left_amp = np.linspace(rho_amp_0,rho_new_amp_1)
    rho_cent_amp = np.linspace(rho_new_amp_1,rho_new_amp_2)
    rho_right_amp = np.linspace(rho_new_amp_2,rho_amp_2_end)

    plt.plot(left_r,p_left_amp)
    plt.plot(center_r,p_cent_amp)
    plt.plot(right_r,p_right_amp)

    plt.plot(left_r,s_left_amp)
    plt.plot(center_r,s_cent_amp)
    plt.plot(right_r,s_right_amp)

    plt.plot(left_r,rho_left_amp)
    plt.plot(center_r,rho_cent_amp)
    plt.plot(right_r,rho_right_amp)

    return ([left_r/6371.,center_r/6371.,right_r/6371.],
           [p_left_amp,p_cent_amp,p_right_amp],
           [s_left_amp,s_cent_amp,s_right_amp],
           [rho_left_amp,rho_cent_amp,rho_right_amp])

def write_list_inc(r,p,s,rho,f):
    p = np.polyfit(r,p,1)
    s = np.polyfit(r,s,1)
    rho = np.polyfit(r,rho,1)

    f.write('{}  {}  {}  {}  {}  {}'.format(
             str(np.round(r[0]*6371,3)),
             str(np.round(r[-1]*6371,3)),
             str(np.round(rho[1],3)),
             str(np.round(rho[0],3)),
             '0.000',
             '0.000 \n'))
    f.write('                {}  {}  {}  {}'.format(
             str(np.round(p[1],3)),
             str(np.round(p[0],3)),
             '0.000',
             '0.000 \n'))
    f.write('                {}  {}  {}  {}'.format(
             str(np.round(p[1],3)),
             str(np.round(p[0],3)),
             '0.000',
             '0.000 \n'))
    f.write('                {}  {}  {}  {}'.format(
             str(np.round(s[1],3)),
             str(np.round(s[0],3)),
             '0.000',
             '0.000 \n'))
    f.write('                {}  {}  {}  {}'.format(
             str(np.round(s[1],3)),
             str(np.round(s[0],3)),
             '0.000',
             '0.000 \n'))
    f.write('                {}  {}  {}  {}  {}  {}'.format(
             '1.000',
             '0.000',
             '0.000',
             '0.000',
             '312.0',
             '57283.0 \n'))

def plot_list(val_list):
    for ii in val_list:
        plt.plot(ii[0],ii[1])
        plt.plot(ii[0],ii[2])
        plt.plot(ii[0],ii[3])
    plt.show()

def visualize(filename):
    inc = read_file(filename)
    val_list = make_layer(inc)
    plot_list(val_list)

def write_output(thick,amplitude):
    sp.call('rm test',shell=True)
    inc = read_file('input.prem.psv')
    #select list increment to vary depth of discontinuity
    val_list = make_layer(inc[18:21])
    r,p,s,rho = make_discont(val_list,thick,amplitude)
    f = open('test','w+')
    for ii in range(0,len(r)):
        write_list_inc(r[ii],p[ii],s[ii],rho[ii],f)

    f.close()
    filename = str(thick)+'_'+str(amplitude)+'.psv'
    #f.close()
    #f = open(str(filename))
    sp.call('cat head_input test tail_input > '+str(filename),shell=True)

def sac_metadata(st):

   for tr in st:
       tr.stats.sac['evla'] = -8.13
       tr.stats.sac['evlo'] = -71.27
       tr.stats.sac['stla'] = 40.290
       tr.stats.sac['stlo'] = -116.500
       tr.stats.sac['o'] = 0
   st = convert.set_gcarc(st)
   return st








