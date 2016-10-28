import astropy.units as u
from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
data1=np.loadtxt('thickice.txt')
data2=np.loadtxt('thinice.txt')
datam=ascii.read('opacity_marco.txt')
#plot wavelenght v/s kappa_abs in cm^2/gram
plt.plot(data1[:,0]*u.micron.to(u.centimeter),data1[:,1],color='c',label ='$Thick \ ice \ mantles \ q=3.5$')
plt.plot(data2[:,0]*u.micron.to(u.centimeter),data2[:,1],color='g',label ='$Thin \ ice \ mantles \ q=3.5$')
plt.plot(datam['lambda'],datam['kappa_q_3_0_porous']*100.,color='r',label='$Porous \ q=3.0$')
plt.plot(datam['lambda'],datam['kappa_q_3_5_porous']*100.,color='r',linestyle='--',label='$Porous \ q=3.5$')
plt.plot(datam['lambda'],datam['kappa_q_3_0_compact']*100.,color='b',label='$Compact \ q=3.0$')
plt.plot(datam['lambda'],datam['kappa_q_3_5_compact']*100.,color='b',linestyle='--',label='$Compact \ q=3.5$')
#plt.plot(datam[:,0],datam[:,1],color='r',label ='$kappaq=3.0 porous$')
plt.legend(loc='lower left')
plt.xscale('log')
plt.yscale('log')
plt.xlabel("$\lambda\; [\mathrm{cm}]$")
plt.ylabel("$\kappa\; [\mathrm{cm}^2/\mathrm{g}]$")
plt.savefig('python_plot.pdf')
plt.show()
