#  new_model_1+3m.py
#  
#
#  Created by Carolina Agurto on 23/05/17.
#
import pylab as plb
import numpy as np
from pylab import *
import scipy.constants as sc
import scipy.stats as ss
import matplotlib.pyplot as plt
import scipy
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp

#SMA DATA 1.3 mm subcompact+extended
data1 = plb.loadtxt('Per50-concat_1mm-f.txt')
uvd1  = sqrt(data1[:,0]**2 + data1[:,1]**2)
wt1 = data1[:,4]
ws1 = sum(data1[:,4])
amp1  = data1[:,2]*data1[:,4].sum() / ws1
sigmadata1 = 1./sqrt(wt1)
nu_1 = 100.

#NOEMA DATA 3 mm subcompact+extended
data3 = plb.loadtxt('Per50-concat_3mm-f.txt')
uvd3  = sqrt(data3[:,0]**2 + data3[:,1]**2)
wt3 = data3[:,4]
ws3 = sum(data3[:,4])
amp3  = data3[:,2]*data3[:,4].sum() / ws3
sigmadata3 = 1./sqrt(wt3)
nu_2 = 230.


#Define the function with uvdistances, flux of the gaussian and flux disk (constant)


#gaussian + constant 1mm

def f1(uv_in, Fd1, Fg1, uv_d):
    fa1 = 0.291*exp(-(uv_in)**2./(2.*uv_d**2.)) + Fd1
    return fa1
#gaussian + constant 3mm

#def f2(uv_in, Fd2, Fg2, uv_d):
#    fa2 = 0.02*exp(-(uv_in)**2./(2.*uv_d**2.)) + Fd2
#    return fa2


def f2(uv_in, Fd2, Fg2, alpha_g, alpha_d, uv_d):
    fa2 = Fg1f1*exp(-(uv_in)**2./(2.*uv_d**2.))*(0.43)**alpha_g+Fd1f1*(0.43)**alpha_d
    return fa2

#guest values for every function 1mm

guessf1 =[0.2, 0.291, 160.]

#guest values for every function 3mm

#guessf2 =[0.02 , 0.06, 300]
guessf2 =[0.02 , 0.06, 3., 3., 300]


popt1,pcov1 = curve_fit(f1, uvd1, amp1, guessf1, sigmadata1, bounds=(0., [0.2, 0.291, 160.]))
popt2,pcov2 = curve_fit(f2, uvd3, amp3, guessf2, sigmadata3, bounds=(0., [0.02, 0.06, 3., 3., 300.]))

#1mm
Fd1f1 = round(popt1[0],3)
Fg1f1 = round(popt1[1],3)
uvd1f1 = round(popt1[2],3)


Fd1f1st = str(Fd1f1)
Fg1f1st = str(Fg1f1)
uvd1f1st = str(uvd1f1)


#3mm
Fd2f2 = round(popt2[0],5)
Fg2f2 = round(popt2[1],5)
al_gf2 = round(popt2[2],5)
al_df2 = round(popt2[3],5)
uvd2f2 = round(popt2[4],5)

Fd2f2st = str(Fd2f2)
Fg2f2st = str(Fg2f2)
al_gf2st = str(al_gf2)
al_df2st = str(al_df2)
uvd2f2st = str(uvd2f2)

#plot models

uv1 = arange(0.,160,10)
uv3 = arange(0.,300,10)

#------------------------------------------------------------------------

#plot data 1mm
plt.plot(uvd1,amp1,'o',label='data')

#plot function 1mm
plt.plot(uv1,f1(uv1,*popt1),'-g',label=r'fit: av=20 Fg='+Fg1f1st+' uvd='+uvd1f1st)
plt.legend()
plt.ylim((0.,1.2))
plt.xlabel('u-v distance (kilo\lambda)')
plt.ylabel('Real amplitude (Jy)')
plt.show()

#------------------------------------------------------------------------

#plot data 3mm
plt.plot(uvd3,amp3,'o',label='data')

#plot function 3mm
plt.plot(uv3,f2(uv3,*popt2),'-g',label=r'fit: Fg='+Fg2f2st+' uvd='+uvd2f2st+' $\alpha$='+al_df2st)
#+' $\beta$='+al_df2st)#+' $\alpha$='+al_gf2st)
plt.legend()
plt.ylim((0.,0.1))
plt.xlabel('u-v distance (kilo\lambda)')
plt.ylabel('Real amplitude (Jy)')
plt.show()

#------------------------------------------------------------------------

#save data

savetxt('gaussianfix_model_test1mm.txt',np.c_[data1[:,0],data1[:,1],f1(uvd1,*popt1),data1[:,3],data1[:,4]])
savetxt('gaussianfix_model_test3mm.txt',np.c_[data3[:,0],data3[:,1],f2(uvd3,*popt2),data3[:,3],data3[:,4]])


