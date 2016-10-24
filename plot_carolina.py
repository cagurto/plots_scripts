import os
import sys
from astropy import units as u
from astropy.coordinates import SkyCoord
import numpy as np
import matplotlib.pyplot as plt

#
# if my path is available
mypath = '/Users/ltesti/git_repository'
pathset = True
for p in sys.path:
    if p == mypath:
        pathset = False

if pathset:
    sys.path.append(mypath)

#
import uvdata as uvlt
#import ltutils.uvdata as uvlt



#
#
#plt.rc('text', usetex=True)
def plot_fun(name,uvbinned,xmin=0.,xmax=630.,ymin=-2.,ymax=10.,autoscale=True):
    if autoscale:
        sfac=0.2
        imin = 1000.*uvbinned['imb'].min()
        imax = 1000.*uvbinned['imb'].max()
        ymin = imin-(imax-imin)*sfac
        ymax = imax+(imax-imin)*sfac
        ymin = min(ymin, -ymax)
        rmin = 1000.*uvbinned['reb'].min()
        rmax = 1000.*uvbinned['reb'].max()
        ymin = min(ymin,rmin-(rmax-rmin)*sfac)
        ymax = max(ymax,rmax+(rmax-rmin)*2.*sfac)
    # plt.rc('font', family='serif')
    fig = plt.figure(figsize=(7, 7))
    fig.subplots_adjust(.18,.15,.9,.9,0,0)
    ax0 = plt.subplot2grid((4,1),(0, 0), rowspan=3)
    ax0.errorbar(bin,1000.*uvbinned['reb'],yerr=1000.*uvbinned['ereb'],fmt='o',color='r')
    ax0.plot([0,1000],[0,0],linestyle='dotted',color='k')
    ax0.set_ylim(ymin,ymax)
    ax0.set_xlim(xmin,xmax)
    ax0.set_ylabel(r'Real (mJy)',fontsize=16)
    ax0.set_xticklabels([])
    plt.title(name,fontsize=18)
    ax1 = plt.subplot2grid((4,1),(3, 0), rowspan=1)
    ax1.errorbar(bin,1000.*uvbinned['imb'],yerr=1000.*uvbinned['eimb'],fmt='o',color='r')
    ax1.plot([0,1000],[0,0],linestyle='dotted',color='k')
    ax1.set_ylim(ymin,-ymin)
    ax1.set_xlim(xmin,xmax)
    ax1.set_ylabel(r'Imag. (mJy)',fontsize=16)
    ax1.set_xlabel(r'baseline (k$\lambda$)',fontsize=16)
    fig.savefig(name+'.pdf')


bin = np.linspace(50.,850.,15)
myuv = uvlt.UVDataMS('dummy',('el29_spw_one_phs.ms',tb))
uvbi = myuv.calc_binned(bin)
plot_fun('Elias_2-29_ALMA_Band_6',uvbi, autoscale=True)

