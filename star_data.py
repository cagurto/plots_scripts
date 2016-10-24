__author__ = 'ltesti - 8 Jul 2015'

import math as m
import numpy as np

star_names = {'fttau' : 'fttau','ft tau' : 'fttau', 'dltau': 'dltau', 'dl tau': 'dltau',
              'hd163296' : 'hd163296', 'hd 163296' : 'hd163296',
              'G34' : 'lupus_iii_87', 'LupusIII_87' : 'lupus_iii_87', 'SSTc2dJ160918.1-390453' : 'lupus_iii_87',
              'dummy' : 'dummy', 'Dummy' : 'dummy', 'DUMMY' : 'dummy'}

star_parameters = dict(
       fttau=dict(Name='FT Tau',lum=0.38, temp=5000., mass=0.85, dist=140., inc=23., pa=29.),
       hd163296=dict(Name='HD 163296',lum=37.7, temp=9250., mass=2.47, dist=122., inc=44., pa=133.),
       dltau=dict(Name='DL Tau',lum=1.12, temp=4060., mass=0.7, dist=140., inc=43., pa=144.),
       dummy=dict(Name='Dummy', lum=1.0, temp=5500., mass=1.0, dist=125., inc=30., pa=0.)
                       #dummy=dict(Name='Dummy', lum=1.0, temp=5500., mass=1.0, dist=160., inc=0., pa=0.)  # to change with correct parameters
       )


class StarData:

    def __init__(self, sname):
        if sname.lower() in star_names:
            outkey = star_names[sname.lower()]
            self.name = star_parameters[outkey]['Name']
            self.pa = star_parameters[outkey]['pa']/180.*m.pi
            self.inc = star_parameters[outkey]['inc']/180.*m.pi
            self.mass = star_parameters[outkey]['mass']
            self.llstar = np.log10(star_parameters[outkey]['lum'])
            self.teff = star_parameters[outkey]['temp']
            self.unknown_star = False
        else:
            print("Star {0} does not recognized!".format(sname))
            self.unkown_star = True

    def __str__(self):
        rep = "Star "+self.name+": M/Msun=%6.2f  Log(L/Lsun)=%9.3e  Teff=%8.1f  pa=%5.1fdeg  inc=%4.1fdeg" % (self.mass,self.llstar,self.teff,self.pa*180./m.pi,self.inc*180./m.pi)

        return rep
