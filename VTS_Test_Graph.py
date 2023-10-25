# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from scipy.interpolate import UnivariateSpline as interpolate

mpl.rcParams['lines.markersize'] = 2
mpl.rcParams["lines.linewidth"] = 0.5

Data = pd.read_csv('data/CBMMD-R_2023-09-27_09-37.txt')

G = 278e9
def surface_resistance(Q):
    Rres = G/Q
    return Rres

def sort(E,Q):
    Q = Q[~np.isnan(Q)]
    E = E[~np.isnan(E)]

    E_indices = np.argsort(E)
    E = E[E_indices]
    Q = Q[E_indices]
    return E, Q

def BCS_resistance(E_4K,E_2K,Q_4K,Q_2K):

    E_4K, Q_4K = sort(E_4K,Q_4K)
    E_2K, Q_2K = sort(E_2K,Q_2K)

    R_4K = surface_resistance(Q_4K)
    R_2K = surface_resistance(Q_2K)

    E_min = max((E_2K.min(),E_4K.min()))
    E_max = min((E_2K.max(),E_4K.max()))

    E = np.linspace(E_min,E_max,num=100)
    
    R_res = interpolate(E_2K,R_2K,s=6e-1*E_2K.shape[0])(E)
    R_BCS = interpolate(E_4K,R_4K,s=6e-1*E_4K.shape[0])(E)-R_res
    return E, R_BCS, R_res


mm = 0.0393701
fig, ax = plt.subplots(3,1,figsize=(86*mm,180*mm),dpi=1200)

ax[0].scatter(Data['Eacc 4.4K Before CBP'],Data['Qo 4.4K Before CBP'],c='m',marker='o',label='As-coated')
ax[0].scatter(Data['Eacc 4.4K After CBP'],Data['Qo 4.4K After CBP'],c='c',marker='^',label='Polished')
ax[0].scatter(Data['Eacc 4.4K After Recoating'],Data['Qo 4.4K After Recoating'],c='y',marker='+',s=20,label='Re-coated')

ax[0].set_yscale('log')
ax[0].set_xlim(0,17)
ax[0].set_ylim(3e8,2e11)
ax[0].grid()

ax[1].scatter(Data['Eacc 2.0K Before CBP'],Data['Qo 2.0K Before CBP'],c='m',marker='o',label='As-coated')
ax[1].scatter(Data['Eacc 2.0K After CBP'],Data['Qo 2.0K After CBP'],c='c',marker='^',label='Polished')
ax[1].scatter(Data['Eacc 2.0K After Recoating'],Data['Qo 2.0K After Recoating'],c='y',marker='+',s=20,label='Re-coated')

ax[1].set_yscale('log')
ax[1].set_xlim(0,17)
ax[1].set_ylim(3e8,2e11)
ax[1].grid()

Data = pd.read_csv('data/CBMMD-R_2023-09-27_09-37.txt')


E, R_BCS, R_res = BCS_resistance(Data['Eacc 4.4K Before CBP'],Data['Eacc 2.0K Before CBP'],Data['Qo 4.4K Before CBP'],Data['Qo 2.0K Before CBP'])
ax[2].plot(E,R_BCS,label='As-coated R$_{BCS}$',c='m',ls='--')
ax[2].plot(E,R_res,label='As-coated R$_{res}$',c='indigo',ls='--')

E, R_BCS, R_res = BCS_resistance(Data['Eacc 4.4K After CBP'],Data['Eacc 2.0K After CBP'],Data['Qo 4.4K After CBP'],Data['Qo 2.0K After CBP'])
ax[2].plot(E,R_BCS,label='Polished R$_{BCS}$',c='c',ls='-')
ax[2].plot(E,R_res,label='Polished R$_{res}$',c='blue',ls='-')

E, R_BCS, R_res = BCS_resistance(Data['Eacc 4.4K After Recoating'],Data['Eacc 2.0K After Recoating'],Data['Qo 4.4K After Recoating'],Data['Qo 2.0K After Recoating'])
ax[2].plot(E,R_BCS,label='Re-coated R$_{BCS}$',c='y',ls='-.')
ax[2].plot(E,R_res,label='Re-coated R$_{res}$',c='darkgoldenrod',ls='-.')

ax[2].set_xlim(0,17)
ax[2].set_ylim(0,40)
ax[2].grid()

ax[0].legend(loc='upper right')
ax[2].legend(loc='upper right', ncol=1, fontsize='xx-small')

ax[2].set_xlabel('Accelerating Field Strength [MV/m]')
ax[0].set_ylabel('Quality Factor')
ax[1].set_ylabel('Quality Factor')
ax[2].set_ylabel('R$_{res}$, R$_{BCS}$(4.4 K) [nOhm]')

ax[0].set_title('Temperature: 4.4 K')
ax[1].set_title('Temperature: 2.0 K')
ax[2].set_title('Surface Resistance')

fig.tight_layout()

plt.savefig('doc/figs/VTS_Test_Graph.png')