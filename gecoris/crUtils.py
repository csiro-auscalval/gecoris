#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""crUtils 

Module containing corner reflector utility functions

Copyright (C) 2020 by R.Czikhardt

Email: czikhardt.richard@gmail.com
Last edit: 28.7.2020

This file is part of GECORIS - Geodetic Corner Reflector (In)SAR Toolbox.

    GECORIS is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    GECORIS is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with GECORIS. If not, see <https://www.gnu.org/licenses/>.
"""

import numpy as np
from shapely.geometry import Polygon

# Default trihedrals:
theta0 = 54.7356
omega0 = 45
shape = 'triangular'
flipped = False
aperture = np.array([[1,0,0], [0,0,1], [0,1,0]]) # triangular
mirror = False

# types of CR:
def TRIANGULAR():
    return (theta0,omega0,aperture,flipped,mirror)

def IGRS():
    flipped = True
    return (theta0,omega0,aperture,flipped,mirror)

def SQUARE():
    aperture= np.array([[1,0,0], [1,0,1], [0,0,1], [0,1,1], [0,1,0], [1,1,0]])
    return (theta0,omega0,aperture,flipped,mirror)

def SQUAREFLIPPED():
    aperture= np.array([[1,0,0], [1,0,1], [0,0,1], [0,1,1], [0,1,0], [1,1,0]])
    flipped = True
    return (theta0,omega0,aperture,flipped,mirror)

def TRIANGULARFLIPPED():
    flipped = True
    return (theta0,omega0,aperture,flipped,mirror)

def SQUAREBASED():
    theta0 = 61.42
    aperture= np.array([[1,0,0], [0,0,1], [0,1,0], [1,1,0]])
    return (theta0,omega0,aperture,flipped,mirror)

# TODO: trimmed, notched, etc.

def crTypeSwitch(argument):
    switcher = {
        'IGRS': IGRS,
        'TRIANGULAR': TRIANGULAR,
        'SQUARE': SQUARE,
        'SQUAREBASED': SQUAREBASED,
        'SQUAREFLIPPED': SQUAREFLIPPED
    }
    # Get the function from switcher dictionary
    func = switcher.get(argument, lambda: "Invalid CR type")
    # Execute the function
    return func()

def crRCS0(crType,a0,zenDip,aziDip,zas,wavelength):
    (theta0,omega0,aperture,flipped,mirror) = crTypeSwitch(crType.upper())
    
    # Compute angle in corner reflector reference frame
    azimuthAngle = np.mod(zas[1]*180/np.pi+360,360)
    if azimuthAngle < 180:
        domega = azimuthAngle-90-aziDip
    else:
        domega = azimuthAngle-270+aziDip
    if flipped:
        theta = 90-(zas[0]*180/np.pi - zenDip)
    else:
        theta = zas[0]*180/np.pi - zenDip
    dtheta = theta-theta0
    omega = (omega0+domega)/180*np.pi
    
    # Equivalent area    
    aeq = crAeq(theta/180*np.pi,omega,aperture,mirror)
    
    # RCS in [dB]:
    RCS = 4*np.pi*aeq**2 * a0**4 / wavelength**2
    RCSdB = 10*np.log10(RCS)
    return RCSdB,dtheta,domega

def crAeq(theta,omega,aperture,mirror):
    
    # compute incidence vector n3 and axis n1 and n2 of projection plane
    n1 =  np.array([-np.sin(omega), np.cos(omega), 0])
    n2 =  np.array([-np.sin(omega)*np.cos(theta), -np.cos(omega)*np.cos(theta),
                    np.sin(theta)])
    n3 =  np.array([-np.cos(omega)*np.sin(theta), -np.sin(omega)*np.sin(theta),
                    -np.cos(theta)])
    
    # Project aperture and inverted aperture onto projection plane
    prj_aperture= aperture @ np.c_[n1, n2]
    inv_aperture= -prj_aperture

    polyPrj_Aperture = Polygon(prj_aperture)
    polyInv_aperture = Polygon(inv_aperture)
    # area:
    #x = inv_aperture.T[0]
    #y = inv_aperture.T[1]
    #Ara = 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))
    Ara = polyInv_aperture.area
    
    # intersect aperture with inverted aperture:
    intersectAperture = polyPrj_Aperture.intersection(polyInv_aperture)
    
    # mirror reflecting faces in xy, xz and yz plane, and project, 
    # necessary for complex shapes, but not necessary for simple shapes
    # TODO: finish!
    
    # Compute effective area
    Aeq = intersectAperture.area
    return Aeq
