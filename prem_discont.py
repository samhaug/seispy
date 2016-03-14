#!/usr/bin/env python

import numpy as np
from scipy.stats import laplace
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d


def make_discont(radius,b,sign):

    domain = np.linspace(3480,5771,num=1500)
    one = np.exp(-1*(radius-domain[domain<radius])/b)
    two = np.exp(-1*(domain[domain>=radius]-radius)/b)

    if sign == 'positive':
        one *= -1
    if sign == 'negative':
        two *= -1
    discont = np.hstack((one,two))
    discont = discont/discont.max()
    return discont,domain,b,radius,sign


def interp_mantle(discont,domain,v_discont,rho_discont,b,radius,sign):
    name = '/home/samhaug/PREM_'+str(radius)+'_b_'+str(b)+'_v_'+str(v_discont)+'_rho_'+str(rho_discont)+'_'+str(sign)+'.bm'

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

    max_index = np.where(discont == discont.max())[0][0]
    rho_d = discont*rho_interp[max_index]*rho_discont
    vpv_d = discont*vpv_interp[max_index]*v_discont
    vsv_d = discont*vsv_interp[max_index]*v_discont
    vph_d = discont*vph_interp[max_index]*v_discont
    vsh_d = discont*vsh_interp[max_index]*v_discont

    new_rho = rho_interp + rho_d
    new_vpv = vpv_interp + vpv_d
    new_vsv = vsv_interp + vsv_d
    new_vph = vph_interp + vph_d
    new_vsh = vsh_interp + vsh_d

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

discont, domain,b,radius,sign = make_discont(3681,25,'positive')
interp_mantle(discont,domain,0.03,0.03,b,radius,sign)


