# -*- coding: utf-8 -*-
"""
Created on Mon May 14 16:57:31 2018
Namedtuple containing all the same fields and values as are in unitConstants.m

@author: David
"""

#from collections import namedtuple
import numpy as np

#ucdict = {"meter": 1., "second": 1., "kg": 1., "C": 1., "Ohm": 1., "Farad": 1.,
#          "Kelvin": 1., "Joule": 1., "Watt": 1., "radian": 1., "km": 1e3,
#          "cm": 1e-2, "mm": 1e-3, "inch": 2.54e-2, "um": 1e-6, "nm": 1e-9,
#          "pm": 1e-12, "mW": 1e-3, "uW": 1e-6, "nW": 1e-9, "pW": 1e-12,
#          "Ampere": 1, "mA": 1e-3, "uA": 1e-6, "Hz": 1, "kHz": 1e3, 
#          "deg": np.pi/180, "mrad": 1e-3, "urad": 1e-3, "arcsec": np.pi/6.48e5,
          
          
#uctuple = namedtuple('uctuple', sorted(ucdict))
#SIC = uctuple(**ucdict)
#SIC = namedtuple('SIC', SIC._fields + ('km'))

#class unitsConstants(object):
meter = 1.
second = 1.
kg = 1.
C = 1.
Ohm = 1.
Farad = 1.
Kelvin = 1.
Joule = 1.
Watt = 1.
radian = 1.
km = 1e3 * meter
cm = 1e-2 * meter
mm = 1e-3 * meter
inch = 2.54 * cm
um = 1e-3 * mm
nm = 1e-3 * um
pm = 1e-3 * nm
mW = 1e-3 * Watt
uW = 1e-3 * mW
nW = 1e-3 * uW
kW = 1e3 * Watt
MW = 1e3 * kW
Ampere = 1 * C / second
mA = 1e-3 * Ampere
uA = 1e-3 * mA
Hz = 1 / second
kHz = 1e3 * Hz
deg = (np.pi / 180) * radian
mrad = 1e-3 * radian
urad = 1e-3 * mrad
arcsec = deg / 3600
mas = 1e-3 * arcsec
uas = 1e-3 * mas
minute = 60 * second
hour = 3600 * second
day = 24 * hour
usec = 1e-6 * second
ppb = 1e-9
ppt = 1e-12

h_planck = 6.626068e-34 #* meter^2 * kg / second
c_light = 299792458 #* meter / second
jupiterRadius = 69911000 #* meter
earthRadius = 6371000 #* meter
sunAbsMag = 4.83 # What units?
lightyear = 9.4607e15 #* meter
AU = 149597870700 #* meter
pc = AU/arcsec