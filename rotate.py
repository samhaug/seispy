#!/usr/bin/env python

from numpy import sin,cos,pi

def rotate_ne_rt(stn,ste):
    '''
    Produce radial and transverse streams from north and east streams

    '''
    if len(stn) != len(ste):
        raise TypeError("North and East streams must have equal no. traces")

    def from_obspy(n, e, ba):
        if len(n) != len(e):
            raise TypeError("North and East component have different length.")
        if ba < 0 or ba > 360:
            raise ValueError("Back Azimuth should be between 0 and 360 degrees.")
        r = e*sin((ba+180)*2*pi/360)+n*cos((ba+180)*2*pi/360)
        t = e*cos((ba+180)*2*pi/360)-n*sin((ba+180)*2*pi/360)
        return r, t

    str = stn.copy()
    stt = stn.copy()
    for idx,tr in enumerate(stn):
        r,t = from_obspy(stn[idx].data,ste[idx].data,stn[idx].stats.baz)
        str[idx].data = r
        str[idx].stats.channel = 'R'
        stt[idx].data = t
        stt[idx].stats.channel = 'T'
    return str,stt
