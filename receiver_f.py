#!/usr/bin/env python

import obspy
import numpy as np

def rotate_ne_rt(st_n, st_e, **kwargs):
    st_r = st_n.copy()
    st_t = st_n.copy()
    for idx, tr in enumerate(st_n):
        r,t = obspy.signal.rotate.rotate_NE_RT(tr.data,
              st_e[idx].data,tr.stats.sac['gcarc'])
        st_r[idx].data = r
        st_t[idx].data = t
