#!/usr/bin/env python

import numpy as np
import obspy
from mpl_toolkits.basemap import Basemap
from obspy.imaging.beachball import Beach
import matplotlib
from matplotlib import pyplot as plt
import os
import myplot.basemap

def mapplot(**kwargs):
   """
   produce map on which to draw objects
   """
   proj = kwargs.get('projection','eck4')
   lat = kwargs.get('lat_0',45)
   lon = kwargs.get('lon_0',-100)
   res = kwargs.get('resolution','l')

   map = Basemap(height=16700000,width=12000000,
            resolution='l',area_thresh=1000.,projection='omerc',
            lon_0=-100,lat_0=15,lon_2=-120,lat_2=65,lon_1=-50,lat_1=-55)

   #map = Basemap(projection=proj,lat_0=lat,lon_0=lon,resolution=res)
   map.drawcoastlines(linewidth=0.25)
   map.drawcountries(linewidth=0.25)
   #map.fillcontinents(color='olive',lake_color='aqua')
   #map.drawmapboundary(fill_color='darkblue')
   return map

def stat_coord(st):
   """
   Return list of station coordinates for a stream object
   """
   coord_list = []
   lat_list = []
   lon_list = []
   for tr in st:
       lat_list.append(tr.stats.sac['stla'])
       lon_list.append(tr.stats.sac['stlo'])
   coord_list = [lat_list,lon_list]
   return coord_list

def add_source(st, map, **kwargs):
    coord = (st[0].stats.sac['evla'],st[0].stats.sac['evlo'])
    x,y = map(coord[1],coord[0])
    #ax = plt.gca()
    #focmecs = [0.136, -0.591, 0.455, -0.396, 0.046, -0.615]
    #bball = Beach(focmecs, xy=(x, y),width=500000, linewidth=1, alpha=0.85)
    #bball.set_zorder(10)
    #ax.add_collection(bball)
    map.scatter(x,y,c='red',s=90,marker='*',lw=0.2)
    return coord

def add_station(coord_list, map, **kwargs):
   """
   From coordinate list, add stations to map
   """
   mark = kwargs.get('marker','v')
   color = kwargs.get('color','yellow')
   size = kwargs.get('size',10)
   lw = kwargs.get('lw',0.2)
   #map = mapplot(**kwargs)
   x,y = map(coord_list[1],coord_list[0])
   for ii in coord_list:
       map.scatter(x,y,c=color,s=size,marker=mark,lw=lw)

def source_reciever_plot(st, **kwargs):
   """
   Plot source and reciever on map
   """
   save = kwargs.get('save',False)

   topo = kwargs.get('topo',False)
   m = mapplot()
   m.drawparallels(np.arange(-80.,81.,20.),labels=[True])
   m.drawmeridians(np.arange(-180.,181.,20.))
   source_coord = add_source(st,m)
   coord_list = stat_coord(st)
   add_station(coord_list,m)
   for ii in range(0,len(coord_list[0]),2):
       m.drawgreatcircle(source_coord[1],source_coord[0],
       coord_list[1][ii],coord_list[0][ii],c='k',lw=0.3,alpha=0.3)
   title = os.getcwd().split('/')
   ax = plt.gca()
   ax.set_title('{} \n Depth (km): {} '.format(
               title[5],round(st[0].stats.sac['evdp'],3)))
   if topo != False:
       myplot.basemap.drawtopography(m,alpha=0.5,cmap=matplotlib.cm.gray)

   if save != False:
       plt.savefig(save+'/map.pdf',format='pdf')
   if save == False:
       plt.show()







