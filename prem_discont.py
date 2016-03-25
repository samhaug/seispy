#!/usr/bin/env python

import numpy as np
from scipy.stats import laplace
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d


def make_discont(discont_list):

    discont_rho = np.zeros(1500)
    discont_v = np.zeros(1500)
    for ii in discont_list:
        domain = np.linspace(3480,5771,num=1500)
        one = np.exp(-1*(ii[0]-domain[domain<ii[0]])/ii[1])
        two = np.exp(-1*(domain[domain>=ii[0]]-ii[0])/ii[1])

        if ii[2] == 'positive':
            one *= -1
        if ii[2] == 'negative':
            two *= -1

        discont = np.hstack((one,two))
        v_d = (discont/discont.max())*ii[3]
        rho_d = (discont/discont.max())*ii[4]
        discont_rho += rho_d
        discont_v  += v_d
    discont_rho += 1
    discont_v += 1

    return discont_v, discont_rho,domain


def interp_mantle(rho_discont,v_discont,domain):
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
    vpv_d = vpv_interp*v_discont
    vsv_d = vsv_interp*v_discont
    vph_d = vph_interp*v_discont
    vsh_d = vsh_interp*v_discont

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

discont_list = [[5521,25,'negative',0.03,0.03],
                [4871,25,'negative',0.03,0.03],
                [4571,25,'positive',0.03,0.03]]

discont_v, discont_rho, domain = make_discont(discont_list)
interp_mantle(discont_rho,discont_v,domain)


