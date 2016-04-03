#!/usr/bin/env python

import numpy as np
from scipy.stats import laplace
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d


def make_discont(discont_list):

    discont_rho = np.zeros(1145)
    discont_vp = np.zeros(1145)
    discont_vs = np.zeros(1145)
    for ii in discont_list:
        domain = np.linspace(3480,5771,num=1145)
        one = np.exp(-1*(ii['h']-domain[domain<ii['h']])/ii['b'])
        two = np.exp(-1*(domain[domain>=ii['h']]-ii['h'])/ii['b'])

        if ii['sign'] == 'positive':
            one *= -1
        if ii['sign'] == 'negative':
            two *= -1

        thickness = ii['thick']/2.
        y1 = one[-1]
        y2 = two[0]
        m = (y2-y1)/thickness
        trans_x = np.linspace(-1*(thickness/2),thickness/2,num=int(thickness))
        trans_y = m*trans_x
        discont = np.hstack((one,trans_y,two))[0:int(-1*thickness)]
        vs_d = (discont/discont.max())*ii['vs']
        vp_d = (discont/discont.max())*ii['vp']
        rho_d = (discont/discont.max())*ii['rho']
        discont_rho += rho_d
        discont_vp  += vp_d
        discont_vs  += vs_d
    discont_rho += 1
    discont_vp += 1
    discont_vs += 1

    return discont_vp, discont_vs, discont_rho, domain


def interp_mantle(rho_discont,vs_discont,vp_discont,domain):
    #name = '/home/samhaug/PREM_'+str(radius)+'_b_'+str(b)+'_v_'+str(v_discont)+'_rho_'+str(rho_discont)+'_'+str(sign)+'.bm'
    name = '/home/samhaug/PREM_discont.bm'
    prem_um = np.loadtxt('prem_um.dat')
    prem_core = np.loadtxt('prem_core.dat')
    prem = np.loadtxt('prem_lm.dat')
    rad = prem[:,0]/1000
    rho = prem[:,1]
    vpv = prem[:,2]
    vsv = prem[:,3]
    qka = prem[:,4]
    qmu = prem[:,5]
    vph = prem[:,6]
    vsh = prem[:,7]
    eta = prem[:,8]

    f_rho = interp1d(rad,rho)
    f_vpv = interp1d(rad,vpv)
    f_vsv = interp1d(rad,vsv)
    f_vph = interp1d(rad,vph)
    f_vsh = interp1d(rad,vsh)
    f_qka = interp1d(rad,qka)
    f_qmu = interp1d(rad,qmu)
    f_eta = interp1d(rad,eta)

    rho_interp = f_rho(domain)
    vpv_interp = f_vpv(domain)
    vsv_interp = f_vsv(domain)
    vph_interp = f_vph(domain)
    vsh_interp = f_vsh(domain)
    qka_interp = f_qka(domain)
    qmu_interp = f_qmu(domain)
    eta_interp = f_eta(domain)

    #max_index = np.where(discont == discont.max())[0][0]
    rho_d = rho_interp*rho_discont
    vpv_d = vpv_interp*vp_discont
    vsv_d = vsv_interp*vs_discont
    vph_d = vph_interp*vp_discont
    vsh_d = vsh_interp*vs_discont

    new_rho = rho_d
    new_vpv = vpv_d
    new_vsv = vsv_d
    new_vph = vph_d
    new_vsh = vsh_d

    out = np.vstack((np.round(domain*1000,3),
                     np.round(new_rho,3),
                     np.round(new_vpv,3),
                     np.round(new_vsv,3),
                     np.round(qka_interp,3),
                     np.round(qmu_interp,3),
                     np.round(new_vph,3),
                     np.round(new_vsh,3),
                     np.round(eta_interp,3)))

    out = np.flipud(np.transpose(out))
    f = open(name,'w')
    f.write('ANELASTIC T \nANISOTROPIC F \nUNITS m \nCOLUMNS radius rho vpv vsv qka qmu vph vsh eta\n')
    f.close()
    full_txt = open(name,'a')
    full = np.vstack((prem_um,out,prem_core))
    np.savetxt(full_txt,full,fmt='%1.3f')
    full_txt.close()

discont_list = [{'h':5221,'b':25,'sign':'positive','thick':10,'vs':0.03,'vp':0.03,'rho':0.03},
#                {'h':4871,'b':25,'sign':'negative','thick':10,'vs':0.03,'vp':0.03,'rho':0.03},
                {'h':4571,'b':25,'sign':'positive','thick':50,'vs':0.03,'vp':0.03,'rho':0.03}]

discont_vp, discont_vs, discont_rho, domain = make_discont(discont_list)
interp_mantle(discont_rho,discont_vs,discont_vp,domain)


