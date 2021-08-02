#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""geoUtils

Module containing geodetic utility functions

Copyright (C) 2021 by R.Czikhardt

Email: czikhardt.richard@gmail.com
Last edit: 1.7.2021

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
import datetime
import time

speedOfLight = 299792458.0


def plh2xyz(plh):
    """Return a numpy.array([x,y,z]) of ECEF coordinates
    
    Transformation of ellipsoidal to ECEF coordinates using GRS80 ellipsoid.
    
    input: numpy.array([phi,lambda,h]) in [radians]
    """
    
    a = 6378137
    f = 0.003352810681182
    e2 = 2*f - f**2
    
    if plh.ndim > 1:
        plh = plh.T
    
    N = np.divide(a,np.sqrt(1 - e2*np.power(np.sin(plh[0]),2)))    
    X = (N+plh[2])*np.cos(plh[0])*np.cos(plh[1])
    Y = (N+plh[2])*np.cos(plh[0])*np.sin(plh[1])
    Z = (N-N*e2+plh[2])*np.sin(plh[0])
    
    xyz = np.array([X,Y,Z])
    if xyz.ndim > 1:
        xyz = xyz.T
    return xyz


def xyz2plh(xyz):
    """Return a numpy.array([phi,lambda,h]) of GRS80 ellipsoidal coordinates
    in [radians]
    
    Transformation of ECEF coordinates to GRS80 ellipsoidal coordinates using
    iteration method.
    
    input: numpy.array([x,y,z]) 
    """
    
    a = 6378137
    f = 0.003352810681182
    e2 = 2*f - f**2
    
    if xyz.ndim > 1:
        xyz = xyz.T
    
    r = np.sqrt(np.power(xyz[0],2)+np.power(xyz[1],2))
    
    #iteration
    i = 1
    maxiter = 5
    Np = xyz[2]
    while i < maxiter:
        phi = np.arctan((xyz[2] + e2*Np)/r)
        N = a/np.sqrt(1 - e2*np.power(np.sin(phi),2))
        Np = N*np.sin(phi)
        i += 1
    
    plh = np.array([phi,np.arctan2(xyz[1],xyz[0]),r/np.cos(phi)-N])
    if plh.ndim > 1:
        plh = plh.T
    return plh


def xyz2zas(xyz,satVec):
    """Return a np.array([z,a,s]) of zenith,azimuth [rad] and distance [m]
    
    Between satellite state vector and ECEF position.
    
    input: np.array([x,y,z]) ECEF coordinates
    """ 
    
    dx = satVec-xyz    
    plh = xyz2plh(xyz)
    
    if plh.ndim > 1:    
        # unit normal vector:
        n = np.array([np.cos(plh[:,0])*np.cos(plh[:,1]),
                      np.cos(plh[:,0])*np.sin(plh[:,1]),
                      np.sin(plh[:,0])]).T
        ip = n[:,0]*dx[:,0] + n[:,1]*dx[:,1] + n[:,2]*dx[:,2]
            
        s = np.sqrt(np.sum(np.power(dx,2),axis=1)) # distance
        z = np.arccos(ip/s) # zenith angle
        a = np.arctan2(-n[:,1]*dx[:,0] + n[:,0]*dx[:,1], ip*-n[:,2] + dx[:,2])
        zas = np.array([z,a,s]).T
    else:
        dx = dx[0]
        n = np.array([np.cos(plh[0])*np.cos(plh[1]),
                      np.cos(plh[0])*np.sin(plh[1]),
                      np.sin(plh[0])]).T
        ip = n[0]*dx[0] + n[1]*dx[1] + n[2]*dx[2]
            
        s = np.sqrt(np.sum(np.power(dx,2))) # distance
        z = np.arccos(ip/s) # zenith angle
        a = np.arctan2(-n[1]*dx[0] + n[0]*dx[1], ip*-n[2] + dx[2])
        zas = np.array([z,a,s])
    
    return zas


def nev2xyz(nev,plh):
    """Return a numpy.array([dx,dy,dz]) of ECEF coordinate differences
    
    Transformation of local topocentric differences to ECEF coordinate
    differences for given point with ellipsoidal plh coordiantes.
    
    input:  numpy.array([n,e,v]) in [m] 
            numpy.array([phi,lambda,h]) in [radians]
    """
    
    n = np.array([np.cos(plh[0])*np.cos(plh[1]),
                 np.cos(plh[0])*np.sin(plh[1]),
                 np.sin(plh[0])])    
    cphi = np.sqrt(1-n[2]**2)
    if nev.ndim > 1:
        dxyz = np.array([(-n[0]*n[2]*nev[:,0] - n[1]*nev[:,1])/cphi + n[0]*nev[:,2],
                         (-n[1]*n[2]*nev[:,0] + n[0]*nev[:,1])/cphi + n[1]*nev[:,2],
                         cphi*nev[:,0] + n[2]*nev[:,2]])
    else:
        dxyz = np.array([(-n[0]*n[2]*nev[0] - n[1]*nev[1])/cphi + n[0]*nev[2],
                        (-n[1]*n[2]*nev[0] + n[0]*nev[1])/cphi + n[1]*nev[2],
                        cphi*nev[0] + n[2]*nev[2]])
    return dxyz


def xyz2nev(xyz,plh):
    """Return a numpy.array([n,e,v]) of topocentric coordiantes
    
    Transformation of ECEF coordinate differences to local topocentric 
    differences to for given point with ellipsoidal plh coordiantes.
    
    input:  numpy.array([n,e,v]) in [m] 
            numpy.array([phi,lambda,h]) in [radians]
    """
    # compute unit normal vector:
    n = np.array([np.cos(plh[0])*np.cos(plh[1]),
                  np.cos(plh[0])*np.sin(plh[1]),
                  np.sin(plh[0])])
    cphi = np.sqrt(1-n[2]**2)
    if xyz.ndim > 1:
        ip = n[0]*xyz[:,0] + n[1]*xyz[:,1] + n[2]*xyz[:,2]
        nev = np.array([(ip*-n[2] + xyz[:,2])/cphi,
                        (-n[1]*xyz[:,0] + n[0]*xyz[:,1])/cphi,
                        ip])
    else:
        ip = n[0]*xyz[0] + n[1]*xyz[1] + n[2]*xyz[2]
        nev = np.array([(ip*-n[2] + xyz[2])/cphi,
                (-n[1]*xyz[0] + n[0]*xyz[1])/cphi,
                ip])
    return nev


def nev2los(nev,heading,incidence):
    """Return a numpy.array of LOS differences
    
    Transformation of local topocentric coordinate differences (nev) to 
    satellite line-of-sight (LOS) differences for given acquisition heading 
    (w.r.t. north) and incidence angle.
    
    input:  numpy.array([n,e,v]) in [m] 
            heading, incidence in [radians]
    """
    if nev.ndim < 2:
        n = nev[0]
        e = nev[1]
        v = nev[2]
    else:
        n = nev[:,0]
        e = nev[:,1]
        v = nev[:,2]    
    los = (v*np.cos(incidence) - np.sin(incidence)*
           (n*np.cos(heading - 3*np.pi/2) + e*np.sin(heading - 3*np.pi/2)))
    return los


def etrf2itrf(xyzETRF,t):
    """Return a np.array([x,y,z]) in ITRF2014 epoch t
    
    Transformation from ETRF2000 into ITRF2014 at epoch t
    
    input: np.array([x,y,z]) in ETRF2000 epoch t0 = 2010.0
    """    
    
    # ITRF2014 - ETRF2000 transformation parameters
    # at epoch 2010.0 (ETRF Technical Note v.1, Altamimi, 2018)
    t0 = 2010.0
    T0 = np.array([54.7/1000,52.2/1000,-74.1/1000]) # translation
    D0 = 2.12e-9; # scale
    Rx0 = 1.701/1e3/3600*np.pi/180; # rotations
    Ry0 = 10.290/1e3/3600*np.pi/180;
    Rz0 = -16.632/1e3/3600*np.pi/180;
    # rates:
    dT = np.array([0.1,0.1,-1.9])/1000;
    dD = 0.11e-9; 
    dRx = 0.081/1e3/3600*np.pi/180;
    dRy = 0.490/1e3/3600*np.pi/180;
    dRz = -0.792/1e3/3600*np.pi/180;
    
    # interpolate transf. parameters to requested epoch:
    T = T0 + dT*(t-t0);
    D = D0 + dD*(t-t0);
    Rx = Rx0 + dRx*(t-t0);
    Ry = Ry0 + dRy*(t-t0);
    Rz = Rz0 + dRz*(t-t0);
    
    xyzITRF = np.empty(xyzETRF.shape)
    # inverse Helmert transformation:
    if xyzETRF.ndim > 1:
        for xyzE,xyzI in zip(xyzETRF,xyzITRF):
            xyzI[0] = xyzE[0] - T[0] - D*xyzE[0] + Rz*xyzE[1] - Ry*xyzE[2]
            xyzI[1] = xyzE[1] - T[1] - D*xyzE[1] - Rz*xyzE[0] + Rx*xyzE[2]
            xyzI[2] = xyzE[2] - T[2] - D*xyzE[2] + Ry*xyzE[0] - Rx*xyzE[1]
    else:
        xyzITRF[0] = xyzETRF[0] - T[0] - D*xyzETRF[0] + Rz*xyzETRF[1] - Ry*xyzETRF[2]
        xyzITRF[1] = xyzETRF[1] - T[1] - D*xyzETRF[1] - Rz*xyzETRF[0] + Rx*xyzETRF[2]
        xyzITRF[2] = xyzETRF[2] - T[2] - D*xyzETRF[2] + Ry*xyzETRF[0] - Rx*xyzETRF[1]
    
    return xyzITRF

#TODO: def itrf2etrf()

def itrf2itrf(xyz, t, t0 = 2010.0):
    """Return a np.array([x,y,z]) in ITRF2014 epoch t
    
    Transformation from ITRF2014 t0=2010.0 into ITRF2014 at epoch t
    
    input: np.array([x,y,z]) in ITRF2014 epoch t0 = 2010.0
    """   
    # Eurasian tectonic plate rates:
    dRx = 0.081/1e3/3600*np.pi/180;
    dRy = 0.490/1e3/3600*np.pi/180;
    dRz = -0.792/1e3/3600*np.pi/180;
    # rotation matrix:
    Rx = dRx*(t-t0);
    Ry = dRy*(t-t0);
    Rz = dRz*(t-t0);
    
    xyzOut = np.empty(xyz.shape)
    # inverse Helmert transformation:
    if xyz.ndim > 1:
        for xyzE,xyzI in zip(xyz,xyzOut):
            xyzI[0] = xyzE[0] + Rz*xyzE[1] - Ry*xyzE[2]
            xyzI[1] = xyzE[1] - Rz*xyzE[0] + Rx*xyzE[2]
            xyzI[2] = xyzE[2] + Ry*xyzE[0] - Rx*xyzE[1]
    else:
        xyzOut[0] = xyz[0] + Rz*xyz[1] - Ry*xyz[2]
        xyzOut[1] = xyz[1] - Rz*xyz[0] + Rx*xyz[2]
        xyzOut[2] = xyz[2] + Ry*xyz[0] - Rx*xyz[1]
    return xyzOut


def decimalYear(date):
    """Return a decimal year
    input: date in 'YYYYMMDD' string format
    """
    
    def sinceEpoch(date): # returns seconds since epoch
        return time.mktime(date.timetuple())
    s = sinceEpoch

    year = date.year
    startOfThisYear = datetime.datetime(year=year, month=1, day=1)
    startOfNextYear = datetime.datetime(year=year+1, month=1, day=1)

    yearElapsed = s(date) - s(startOfThisYear)
    yearDuration = s(startOfNextYear) - s(startOfThisYear)
    fraction = yearElapsed/yearDuration

    return date.year + fraction

def roundTime(dt=None, roundTo=60):
   """Round a datetime object to any time lapse in seconds
   dt : datetime.datetime object, default now.
   roundTo : Closest number of seconds to round to, default 1 minute.
   Author: Thierry Husson 2012 - Use it as you want but don't blame me.
   """
   if dt == None : dt = datetime.datetime.now()
   seconds = (dt.replace(tzinfo=None) - dt.min).seconds
   rounding = (seconds+roundTo/2) // roundTo * roundTo
   return dt + datetime.timedelta(0,rounding-seconds,-dt.microsecond)

def sod2UTC(seconds):
    """Return UTC string
    input: seconds of day (sod)
    """
    return time.strftime('T%H:%M:%S', time.gmtime(seconds))
    

def orbitFit(orbit, verbose = 0, plotFlag = 0):
    """Return a orbitFit dict
    
    Satellite state vector interpolation using Chebyshev polynomials of 
    7th order (according to DLR recommendations). Function returns Chebyshev 
    polynomial coefficients. Use these to evaluate orbit state at given time
    via function 'orbitVal'.
    
    input: snappy 'orbit' object (as read by 'readMetadata' function)
    """
    
    # parse masterorb: 
    t = np.array([orb.time_mjd for orb in orbit])
    x = np.array([orb.x_pos for orb in orbit])
    y = np.array([orb.y_pos for orb in orbit])
    z = np.array([orb.z_pos for orb in orbit])
    x_vel = np.array([orb.x_vel for orb in orbit])
    y_vel = np.array([orb.y_vel for orb in orbit])
    z_vel = np.array([orb.z_vel for orb in orbit])
    t = np.mod(t,1)*24*3600 # convert t to [seconds] of day
    
    # interpolate orbits using Chebyshev polynomials of 7th order:
    t0 = (min(t) + max(t))/2
    px = t - t0 # time argument px (centered around mid interval)
    cX = np.polynomial.chebyshev.chebfit(px,x,7) # position
    cY = np.polynomial.chebyshev.chebfit(px,y,7)
    cZ = np.polynomial.chebyshev.chebfit(px,z,7)
    cVX = np.polynomial.chebyshev.chebfit(px,x_vel,7) # velocity
    cVY = np.polynomial.chebyshev.chebfit(px,y_vel,7)
    cVZ = np.polynomial.chebyshev.chebfit(px,z_vel,7)
    cAX = np.polynomial.chebyshev.chebder(cVX) # acceleration
    cAY = np.polynomial.chebyshev.chebder(cVY)
    cAZ = np.polynomial.chebyshev.chebder(cVZ)
    
    if verbose:
        # position fit residuals:
        xRes = np.polynomial.chebyshev.chebval(px,cX) - x
        yRes = np.polynomial.chebyshev.chebval(px,cY) - y
        zRes = np.polynomial.chebyshev.chebval(px,cZ) - z
        xStd = np.std(xRes)
        yStd = np.std(yRes)
        zStd = np.std(zRes)
        print(f'Orbit fit position residuals: X {xStd:.4f} m, Y {yStd:.4f} m, Z {zStd:.4f} m. ')
        # velocity residuals:
        vxRes = np.polynomial.chebyshev.chebval(px,np.polynomial.chebyshev.chebder(cX)) - x_vel
        vyRes = np.polynomial.chebyshev.chebval(px,np.polynomial.chebyshev.chebder(cY)) - y_vel
        vzRes = np.polynomial.chebyshev.chebval(px,np.polynomial.chebyshev.chebder(cZ)) - z_vel
        vxStd = np.std(vxRes)
        vyStd = np.std(vyRes)
        vzStd = np.std(vzRes)
        print(f'Orbit fit velocity residuals: vX {vxStd:.4f} m/s, vY {vyStd:.4f} m/s, vZ {vzStd:.4f} m/s. ')
    ### fit using regular polynomials like this:
    #cX2 = np.polynomial.polynomial.polyfit(px,y,4)
    #xRes2 = np.polynomial.polynomial.polyval(px,cX2) - y
    
    if plotFlag:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots() 
        ax.plot(px, xRes, label='xRes')
        ax.plot(px, yRes, label='yRes')
        ax.plot(px, zRes, label='zRes')
        ax.set_ylabel('residuals [m]')
        ax.set_xlabel('t [s]')
        ax.set_title("Orbit fit position residuals")
        ax.legend()
    
    orbitFit = dict()
    orbitFit["t0"] = t0
    orbitFit["cX"] = cX
    orbitFit["cY"] = cY
    orbitFit["cZ"] = cZ
    orbitFit["cVX"] = cVX
    orbitFit["cVY"] = cVY
    orbitFit["cVZ"] = cVZ
    orbitFit["cAX"] = cAX
    orbitFit["cAY"] = cAY
    orbitFit["cAZ"] = cAZ
    
    return orbitFit


def orbitVal(orbitFit, tAzi):
    """Return a orbit state vector and velocity
    
    Eevaluate orbit state vector at given azimuth time using fitted Chebyshev 
    polynomials by function 'orbitFit'.
    
    input: 'orbitFit' dict (as generated by 'orbitFit' function)
            tAz - azimuth time [seconds of day]
    """    
    
    # evalueate position at requested azimuth epoch:
    x = np.polynomial.chebyshev.chebval(tAzi-orbitFit['t0'],orbitFit['cX'])
    y = np.polynomial.chebyshev.chebval(tAzi-orbitFit['t0'],orbitFit['cY'])
    z = np.polynomial.chebyshev.chebval(tAzi-orbitFit['t0'],orbitFit['cZ'])

    # velocity at epoch t
    vX = np.polynomial.chebyshev.chebval(tAzi-orbitFit['t0'],orbitFit['cVX'])
    vY = np.polynomial.chebyshev.chebval(tAzi-orbitFit['t0'],orbitFit['cVY'])
    vZ = np.polynomial.chebyshev.chebval(tAzi-orbitFit['t0'],orbitFit['cVZ'])

    
    satVec = np.array([x,y,z])
    
    satVel = np.array([vX,vY,vZ])
    
    return satVec,satVel


def tropoDelay(xyz,satVec, method='approx'):
    """Return a slant tropospheric delay in [s]
    
    Between satellite state vector and ECEF position.
    
    input: np.array([x,y,z]) ECEF coordinates and satVec
    """     
    
    plh = xyz2plh(xyz)
    zas = xyz2zas(xyz,satVec)
    if plh.ndim > 1:
        h = plh[:,2]
        z = zas[:,0]
    else:
        h = plh[2]
        z = zas[0]
    
    if method == 'approx':
        zenithDelay = h**2/(8.55e7) - h/3411 + 2.41
        slantDelay = zenithDelay/np.cos(z)/speedOfLight;
    else:
        raise('Requested tropo. delay estiamtion method not implemented.')
    return slantDelay


def getSET(xyzP,date,t):
    """
    % function to compute effect of SET in tide-free model on ECEF coordinates
    % of station in ITRF
    
    % according to IERS conventions 2003
    % matlab original kindly provided by Dr. Branislav Habel
    % Department of Theoretical Geodesy, SUT Bratislava, Slovakia
    
    % ------
    % input:
    % - xyzP   - TRF position of reflector
    % - date   - date of acquisition ('YYYYMMDD')
    % - t      - seconds of day for satellite zero-Doppler position
    
    % output:
    % - xyzSET - effect of SET on ECEF coordinates in [m]
    % ------
    
    % part of Geodetic radar reflector toolbox & g-mapit
    % 15.6.2020
    % R. Czikhardt
    %% ------------------------------------------------------------------------
    """
    import numpy as np
    
    X = xyzP[0]
    Y = xyzP[1]
    Z = xyzP[2]
    
    # time to UTC:
    YY = float(date[:4])
    MM = float(date[4:6])
    DD = float(date[6:])
    hh = np.floor(t/3600)
    mm = np.floor((t-hh*3600)/60)
    ss = t-hh*3600-mm*60
    
    if MM <= 2:
        YY -= 1
        MM += 12
    
    A = np.fix(YY/100)
    B = 2 - A + np.fix(A/4)
    juld = (np.fix(365.25*(YY + 4716)) + np.fix(30.6001*(MM + 1)) + DD + B
            - 1524.5 + hh/24 + mm/(24*60) + ss/(24*60*60))
    
    MJD = juld - 2400000.5
    #time = datenum(YY,MM,DD,hh,mm,ss)
    
    #constants:
    deg = np.pi / 180; minu = np.pi / (180*60); sec = np.pi / (180*60*60)
    
    Re  = 6378136.6;                #rovnikovy polomer Zeme,
    GMe = 3.986004418e14;       #geocentricka gravitacna konstanta Zeme,
    GMm = 4.903e12;             #geocentricka gravitacna konstanta Mesiaca,
    GMs = 1.32712442076e20;     #geocentricka gravitacna konstanta Slnka,
    AU  = 1.496e11;             #astronomicka jednotka;
    
    #%GEOCENTRICKE SURADNICE STANOVISKA:   
    fi     = np.arctan2(Z,np.sqrt(X*X + Y*Y))
    lamb = np.arctan2(Y,X)
    
    r  = np.sqrt(X*X + Y*Y + Z*Z)   #velkost geocentrickeho vektora stanoviska,
    r0 = np.array([X, Y, Z])/r             #jednotkovy geocentricky vektor stanoviska,
    
    #% rotation matrix
    Rot = np.array([[-np.sin(fi)*np.cos(lamb), -np.sin(fi)*np.sin(lamb), np.cos(fi)],
                    [-np.sin(lamb), np.cos(lamb), 0],
                    [np.cos(fi)*np.cos(lamb), np.cos(fi)*np.sin(lamb), np.sin(fi)]]) #rotacna matica
     
    #PREVOD MJD:
    DJS = 36525;            #dlzka julianskeho storocia v dnoch,
    DJT = 365250;           #dlzka julianskeho tisicrocia v dnoch,
    JD0 = 2451545;          #juliansky datum v epoche J2000.0,
    JD  = MJD + 2400000.5;  #juliansky datum v okamihu MJD,
    JS  = (JD - JD0)/DJS;   #julianske storocie,
    JM  = (JD - JD0)/DJT;   #julianske tisicrocie
      
    #GREENWICHSKY HVIEZDNY CAS GST(IERS2003, str.47):
    #Argumenty lunisolarnej nutacie(IERS2003, str.48):
    l     = (134.96340251*deg + 1717915923.217800*sec*JS + 
        31.879200*sec*JS**2 + 0.05163500*sec*JS**3 - 
        0.0002447000*sec*JS**4);
    ll    = (357.52910918*deg +  129596581.048100*sec*JS - 
        0.553200*sec*JS**2 + 0.00013600*sec*JS**3 - 0.0000114900*sec*JS**4);
    F     =  (93.27209062*deg + 1739527262.847800*sec*JS - 
        12.751200*sec*JS**2 - 0.00103700*sec*JS**3 + 
        0.0000041700*sec*JS**4);
    D     = (297.85019547*deg + 1602961601.209000*sec*JS - 
        6.370600*sec*JS**2 + 0.00659300*sec*JS**3 - 0.0000316900*sec*JS**4);
    Omega = (125.04455501*deg - 6962890.543100*sec*JS + 
        7.472200*sec*JS**2 + 0.00770200*sec*JS**3 - 0.0000593900*sec*JS**4);
    
    #Argumenty planetarnej nutacie(IERS2003, str.49):
    l_Me = 4.402608842 + 2608.7903141574*JS;
    l_Ve = 3.176146697 + 1021.3285546211*JS;
    l_E = 1.753470314 + 628.3075849991*JS;
    l_Ma = 6.203480913 + 334.0612426700*JS;
    l_Ju = 0.599546497 + 52.9690962641*JS;
    l_Sa = 0.874016757 + 21.3299104960*JS;
    l_Ur = 5.481293872 + 7.4781598567*JS;
    l_Ne = 5.311886287 + 3.8133035638*JS;
    pa = 0.024381750*JS + 0.00000538691*JS**2;
     
    #Zlozka theta(UT1):
    T_u   = JD - 2451545.0;
    theta = 2*np.pi*(0.7790572732640 + 1.00273781191135448*T_u);
    
    #Polynomicka zlozka:
    PZ = (0.014506 + 4612.15739966*JS + 1.39667721*JS**2 -
        0.00009344*JS**3 + 0.00001882*JS**4)*sec;
    
    #Zlozka delta_psi*np.cos(epsilon_A):
    dPcE = 0; #zanedbana vo vypocte
    
    #%Nepolynomicka zlozka:
    #Argumenty lunisolarnej a planetarnej nutacie
    #--------------------------------------------------------------------------------------------------------------
    #     i    C_{s,j})_i     C_{c,j})_i    l    l'   F    D   Om L_Me L_Ve  L_E L_Ma  L_J L_Sa  L_U L_Ne  p_A
    #--------------------------------------------------------------------------------------------------------------
    A = np.array([
         [1,2640.96,-0.39,0,0,0,0,1,0,0,0,0,0,0,0,0,0],
         [2,63.52,-0.02,0,0, 0, 0, 2,0, 0,  0,0,0,0,0,0, 0],
         [3,11.75, 0.01,0,0, 2,-2, 3,0, 0,  0,0,0,0,0,0, 0],
         [4,11.21, 0.01,0,0, 2,-2, 1,0, 0,  0,0,0,0,0,0, 0],
         [5,-4.55, 0.00,0,0, 2,-2, 2,0, 0,  0,0,0,0,0,0, 0],
         [6, 2.02, 0.00,0,0, 2, 0, 3,0, 0,  0,0,0,0,0,0, 0],
         [7, 1.98, 0.00,0,0, 2, 0, 1,0, 0,  0,0,0,0,0,0, 0],
         [8,-1.72, 0.00,0,0, 0, 0, 3,0, 0,  0,0,0,0,0,0, 0],
         [9,-1.41,-0.01,0,1, 0, 0, 1,0, 0,  0,0,0,0,0,0, 0],
        [10,-1.26,-0.01,0,1, 0, 0,-1,0, 0,  0,0,0,0,0,0, 0],
        [11,-0.63, 0.00,1,0, 0, 0,-1,0, 0,  0,0,0,0,0,0, 0],
        [12,-0.63, 0.00,1,0, 0, 0, 1,0, 0,  0,0,0,0,0,0, 0],
        [13, 0.46, 0.00,0,1, 2,-2, 3,0, 0,  0,0,0,0,0,0, 0],
        [14, 0.45, 0.00,0,1, 2,-2, 1,0, 0,  0,0,0,0,0,0, 0],
        [15, 0.36, 0.00,0,0, 4,-4, 4,0, 0,  0,0,0,0,0,0, 0],
        [16,-0.24,-0.12,0,0, 1,-1, 1,0,-8, 12,0,0,0,0,0, 0],
        [17, 0.32, 0.00,0,0, 2, 0, 0,0, 0,  0,0,0,0,0,0, 0],
        [18, 0.28, 0.00,0,0, 2, 0, 2,0, 0,  0,0,0,0,0,0, 0],
        [19, 0.27, 0.00,1,0, 2, 0, 3,0, 0,  0,0,0,0,0,0, 0],
        [20, 0.26, 0.00,1,0, 2, 0, 1,0, 0,  0,0,0,0,0,0, 0],
        [21,-0.21, 0.00,0,0, 2,-2, 0,0, 0,  0,0,0,0,0,0, 0],
        [22, 0.19, 0.00,0,1,-2, 2,-3,0, 0,  0,0,0,0,0,0, 0],
        [23, 0.18, 0.00,0,1,-2, 2,-1,0, 0,  0,0,0,0,0,0, 0],
        [24,-0.10, 0.05,0,0, 0, 0, 0,0, 8,-13,0,0,0,0,0,-1],
        [25, 0.15, 0.00,0,0, 0, 2, 0,0, 0,  0,0,0,0,0,0, 0],
        [26,-0.14, 0.00,2,0,-2, 0,-1,0, 0,  0,0,0,0,0,0, 0],
        [27, 0.14, 0.00,1,0, 0,-2, 1,0, 0,  0,0,0,0,0,0, 0],
        [28,-0.14, 0.00,0,1, 2,-2, 2,0, 0,  0,0,0,0,0,0, 0],
        [29, 0.14, 0.00,1,0, 0,-2,-1,0, 0,  0,0,0,0,0,0, 0],
        [30, 0.13, 0.00,0,0, 4,-2, 4,0, 0,  0,0,0,0,0,0, 0],
        [31,-0.11, 0.00,0,0, 2,-2, 4,0, 0,  0,0,0,0,0,0, 0],
        [32, 0.11, 0.00,1,0,-2, 0,-3,0, 0,  0,0,0,0,0,0, 0],
        [33, 0.11, 0.00,1,0,-2, 0,-1,0, 0,  0,0,0,0,0,0, 0],
        [34,-0.87, 0.00,0,0, 0, 0, 1,0, 0,  0,0,0,0,0,0, 0]]);
    
    ARG_LN = l*A[:,3] + ll*A[:,4] + F*A[:,5] + D*A[:,6] + Omega*A[:,7]
    ARG_PN = l_Me*A[:,8] + l_Ve*A[:,9] + l_E*A[:,10] + l_Ma*A[:,11] + l_Ju*A[:,12] + l_Sa*A[:,13] + l_Ur*A[:,14] + l_Ne*A[:,15] + pa*A[:,16]
    NZ = ((A[:,1]*10**(-6)*sec)*np.sin(ARG_LN+ARG_PN) 
        + (A[:,2]*10**(-6)*sec)*np.cos(ARG_LN+ARG_PN))
    NZ = np.sum(NZ)
    
    GST = (theta + PZ + NZ + dPcE); #vysledny Grenwichsky hviezdny cas v okamihu MJD;
    GST = np.remainder(GST,2*np.pi)
    GST = GST + 2*np.pi if GST < 0 else GST
    
    ##POLOHA SLNKA:  
    #Sklon ekliptiky:     
    ecl = 23*deg + 26*minu + 21.448*sec - 46.8150*sec*JS - 0.00059*sec*JS**2 + 0.001813*sec*JS**3;
    
    #Periodicke cleny pre ekliptikalnu dlzku Zeme:
    L0 = np.array([   175347046 *np.cos(0        +     0         *JM),
              3341656 *np.cos(4.6692568+  6283.0758500 *JM),  
                34894 *np.cos(4.62610  + 12566.15170   *JM),
                 3497 *np.cos(2.7441   +  5753.3849    *JM),   
                 3418 *np.cos(2.8289   +     3.5231    *JM), 
                 3136 *np.cos(3.6277   + 77713.7715    *JM),
                 2676 *np.cos(4.4181   +  7860.4194    *JM),
                 2343 *np.cos(6.1352   +  3930.2097    *JM),
                 1324 *np.cos(0.7425   + 11506.7698    *JM),
                 1273 *np.cos(2.0371   +   529.6910    *JM),
                 1199 *np.cos(1.1096   +  1577.3435    *JM),
                  990 *np.cos(5.233    +  5884.927     *JM),
                  902 *np.cos(2.045    +    26.298     *JM),
                  857 *np.cos(3.508    +   398.149     *JM),
                  780 *np.cos(1.179    +  5223.694     *JM),
                  753 *np.cos(2.533    +  5507.553     *JM),
                  505 *np.cos(4.583    + 18849.228     *JM),
                  492 *np.cos(4.205    +   775.523     *JM),
                  357 *np.cos(2.920    +     0.067     *JM),
                  317 *np.cos(5.849    + 11790.629     *JM),
                  284 *np.cos(1.899    +   796.298     *JM),
                  271 *np.cos(0.315    + 10977.079     *JM),
                  243 *np.cos(0.345    +  5486.778     *JM),
                  206 *np.cos(4.806    +  2544.314     *JM),
                  205 *np.cos(1.869    +  5573.143     *JM),
                  202 *np.cos(2.458    +  6069.777     *JM),
                  156 *np.cos(0.833    +   213.299     *JM),
                  132 *np.cos(3.411    +  2942.463     *JM),
                  126 *np.cos(1.083    +    20.775     *JM),
                  115 *np.cos(0.645    +     0.980     *JM),
                  103 *np.cos(0.636    +  4694.003     *JM),
                  102 *np.cos(0.976    + 15720.839     *JM),
                  102 *np.cos(4.267    +     7.114     *JM),
                   99 *np.cos(6.21     +  2146.17      *JM),
                   98 *np.cos(0.68     +   155.42      *JM),
                   86 *np.cos(5.98     +161000.69      *JM),
                   85 *np.cos(1.30     +  6275.96      *JM),
                   85 *np.cos(3.67     + 71430.70      *JM),
                   80 *np.cos(1.81     + 17260.15      *JM),
                   79 *np.cos(3.04     + 12036.46      *JM),
                   75 *np.cos(1.76     +  5088.63      *JM),
                   74 *np.cos(3.50     +  3154.69      *JM),
                   74 *np.cos(4.68     +   801.82      *JM),
                   70 *np.cos(0.83     +  9437.76      *JM),
                   62 *np.cos(3.98     +  8827.39      *JM),
                   61 *np.cos(1.82     +  7084.90      *JM),
                   57 *np.cos(2.78     +  6286.60      *JM),
                   56 *np.cos(4.39     + 14143.50      *JM),
                   56 *np.cos(3.47     +  6279.55      *JM),
                   52 *np.cos(0.19     + 12139.55      *JM),
                   52 *np.cos(1.33     +  1748.02      *JM),
                   51 *np.cos(0.28     +  5856.48      *JM),
                   49 *np.cos(0.49     +  1194.45      *JM),
                   41 *np.cos(5.37     +  8429.24      *JM),
                   41 *np.cos(2.40     + 19651.05      *JM),
                   39 *np.cos(6.17     + 10447.39      *JM),
                   37 *np.cos(6.04     + 10213.29      *JM),
                   37 *np.cos(2.57     +  1059.38      *JM),
                   36 *np.cos(1.71     +  2352.87      *JM),
                   36 *np.cos(1.78     +  6812.77      *JM),
                   33 *np.cos(0.59     + 17789.85      *JM),
                   30 *np.cos(0.44     + 83996.85      *JM),
                   30 *np.cos(2.74     +  1349.87      *JM),
                   25 *np.cos(3.16     +  4690.48      *JM) ]);
    L1 = np.array([628331966747 *np.cos(0        +     0         *JM),
               206059 *np.cos(2.678235 +  6283.075850  *JM),
                 4303 *np.cos(2.6351   + 12566.1517    *JM),
                  425 *np.cos(1.590    +     3.523     *JM),
                  119 *np.cos(5.796    +    26.298     *JM),
                  109 *np.cos(2.966    +  1577.344     *JM),
                   93 *np.cos(2.59     + 18849.23      *JM),
                   72 *np.cos(1.14     +   529.69      *JM),
                   68 *np.cos(1.87     +   398.15      *JM),
                   67 *np.cos(4.41     +  5507.55      *JM),
                   59 *np.cos(2.89     +  5223.69      *JM),
                   56 *np.cos(2.17     +   155.42      *JM),
                   45 *np.cos(0.40     +   796.30      *JM),
                   36 *np.cos(0.47     +   775.52      *JM),
                   29 *np.cos(2.65     +     7.11      *JM),
                   21 *np.cos(5.34     +     0.98      *JM),
                   19 *np.cos(1.85     +   486.78      *JM),
                   19 *np.cos(4.97     +   213.30      *JM),
                   17 *np.cos(2.99     +  6275.96      *JM),
                   16 *np.cos(0.03     +  2544.31      *JM),
                   16 *np.cos(1.43     +  2146.17      *JM),
                   15 *np.cos(1.21     + 10977.08      *JM),
                   12 *np.cos(2.83     +  1748.02      *JM),
                   12 *np.cos(3.26     +  5088.63      *JM),
                   12 *np.cos(5.27     +  1194.45      *JM),
                   12 *np.cos(2.08     +  4694.00      *JM),
                   11 *np.cos(0.77     +   553.57      *JM),
                   10 *np.cos(1.30     +  6286.60      *JM),
                   10 *np.cos(4.24     +  1349.87      *JM),
                    9 *np.cos(2.70     +   242.73      *JM),
                    9 *np.cos(5.64     +   951.72      *JM),
                    8 *np.cos(5.30     +  2352.87      *JM),
                    6 *np.cos(2.65     +  9437.76      *JM),
                    6 *np.cos(4.67     +  4690.48      *JM) ]);
    L2 = np.array([       52919 *np.cos(0        +     0         *JM),
                 8720 *np.cos(1.0721   +  6283.0758    *JM),   
                  309 *np.cos(0.867    + 12566.152     *JM),
                   27 *np.cos(0.05     +      3.52     *JM),
                   16 *np.cos(5.19     +     26.30     *JM),
                   16 *np.cos(3.68     +    155.42     *JM),
                   10 *np.cos(0.76     +  18849.23     *JM),
                    9 *np.cos(2.06     +  77713.77     *JM),
                    7 *np.cos(0.83     +    775.52     *JM),
                    5 *np.cos(4.66     +   1577.34     *JM),
                    4 *np.cos(1.03     +      7.11     *JM),
                    4 *np.cos(3.44     +   5573.14     *JM),
                    3 *np.cos(5.14     +    769.30     *JM),
                    3 *np.cos(6.05     +   5507.55     *JM),
                    3 *np.cos(1.19     +    242.73     *JM),
                    3 *np.cos(6.12     +    529.69     *JM),
                    3 *np.cos(0.31     +    398.15     *JM),
                    3 *np.cos(2.28     +    553.57     *JM),
                    2 *np.cos(4.38     +   5223.69     *JM),
                    2 *np.cos(3.75     +      0.98     *JM) ]);
    L3 = np.array([         289 *np.cos(5.844    +   6283.076    *JM),
                   35 *np.cos(0        +      0        *JM),
                   17 *np.cos(5.49     +  12566.15     *JM),
                    3 *np.cos(5.20     +    155.42     *JM),
                    1 *np.cos(4.72     +      3.52     *JM),
                    1 *np.cos(5.30     +  18849.23     *JM),
                    1 *np.cos(5.97     +    242.73     *JM) ]);
    L4 = np.array([         114 *np.cos(3.142    +      0        *JM),
                    8 *np.cos(4.13     +   6283.08     *JM),
                    1 *np.cos(3.84     +  12566.15     *JM) ]);
    L5 =             1 *np.cos(3.14     +      0        *JM)  ;
    
    #Periodicke cleny pre ekliptikalnu sirku Zeme:
    B0 = np.array([ 280 *np.cos(3.199 + 84334.662 *JM),
          102 *np.cos(5.422 +  5507.553 *JM),
           80 *np.cos(3.88  +  5223.69  *JM),
           44 *np.cos(3.70  +  2352.87  *JM),
           32 *np.cos(4.00  +  1577.34  *JM) ]);   
    B1 = np.array([   9 *np.cos(3.90  +  5507.55  *JM),
            6 *np.cos(1.73  +  5223.69  *JM) ]);
    
    #Periodicke cleny pre sprievodic Zeme:
    R0 = np.array([ 100013989*np.cos(0         +     0         *JM),
            1670700*np.cos(3.0984635 +  6283.0758500 *JM),
              13956*np.cos(3.05525   + 12566.15170   *JM),
               3084*np.cos(5.1985    + 77713.7715    *JM),
               1628*np.cos(1.1739    +  5753.3849    *JM),
               1576*np.cos(2.8469    +  7860.4194    *JM),
                925*np.cos(5.453     + 11506.770     *JM),
                542*np.cos(4.564     +  3930.210     *JM),
                472*np.cos(3.661     +  5884.927     *JM),
                346*np.cos(0.964     +  5507.553     *JM),
                329*np.cos(5.900     +  5223.694     *JM),
                307*np.cos(0.299     +  5573.143     *JM),
                243*np.cos(4.273     + 11790.629     *JM),
                212*np.cos(5.847     +  1577.344     *JM),
                186*np.cos(5.022     + 10977.079     *JM),
                175*np.cos(3.012     + 18849.228     *JM),
                110*np.cos(5.055     +  5486.778     *JM),
                 98*np.cos(0.89      +  6069.78      *JM), 
                 86*np.cos(5.69      + 15720.84      *JM),
                 86*np.cos(1.27      +161000.69      *JM),
                 65*np.cos(0.27      + 17260.15      *JM),
                 63*np.cos(0.92      +   529.69      *JM),
                 57*np.cos(2.01      + 83996.85      *JM),
                 56*np.cos(5.24      + 71430.70      *JM),
                 49*np.cos(3.25      +  2544.31      *JM),
                 47*np.cos(2.58      +   775.52      *JM),
                 45*np.cos(5.54      +  9437.76      *JM),
                 43*np.cos(6.01      +  6275.96      *JM),
                 39*np.cos(5.36      +  4694.00      *JM),
                 38*np.cos(2.39      +  8827.39      *JM),
                 37*np.cos(0.83      + 19651.05      *JM),
                 37*np.cos(4.90      + 12139.55      *JM),
                 36*np.cos(1.67      + 12036.46      *JM),
                 35*np.cos(1.84      +  2942.46      *JM),
                 33*np.cos(0.24      +  7084.90      *JM),
                 32*np.cos(0.18      +  5088.63      *JM),
                 32*np.cos(1.78      +   398.15      *JM),
                 28*np.cos(1.21      +  6286.60      *JM),
                 28*np.cos(1.90      +  6279.55      *JM),
                 26*np.cos(4.59      + 10447.39      *JM) ]);
    R1 = np.array([    103019*np.cos(1.107490  +  6283.075850  *JM),
               1721*np.cos(1.0644    + 12566.1517    *JM),
                702*np.cos(3.142     +     0         *JM),
                 32*np.cos(1.02      + 18849.23      *JM),
                 31*np.cos(2.84      +  5507.55      *JM),
                 25*np.cos(1.32      +  5223.69      *JM),
                 18*np.cos(1.42      +  1577.34      *JM),
                 10*np.cos(5.91      + 10977.08      *JM),
                  9*np.cos(1.42      +  6275.96      *JM),
                  9*np.cos(0.27      +  5486.78      *JM) ]); 
    R2 = np.array([       359*np.cos(5.7846    +  6283.0758    *JM),
                124*np.cos(5.579     + 12566.152     *JM),
                 12*np.cos(3.14      +     0         *JM),
                  9*np.cos(3.63      + 77713.77      *JM),
                  6*np.cos(1.87      +  5573.14      *JM),
                  3*np.cos(5.47      + 18849.23      *JM) ]);
    R3 = np.array([       145*np.cos(4.273     +  6283.076     *JM),
                  7*np.cos(3.92      + 12566.15      *JM) ]);     
    R4 =           4*np.cos(2.56      +  6283.078     *JM)  ;      
    
    #Heliocentricke ekliptikalne suradnice Zeme:
    LE = (np.sum(L0) + np.sum(L1)*JM + np.sum(L2)*JM**2 + np.sum(L3)*JM**3 + np.sum(L4)*JM**4 + np.sum(L5)*JM**5)*10**-8;
    BE = (np.sum(B0) + np.sum(B1)*JM)*10**-8;
    RE = AU*(np.sum(R0) + np.sum(R1)*JM + np.sum(R2)*JM**2 + np.sum(R3)*JM**3 + np.sum(R4)*JM**4)*10**-8; 
    
    #Geocentricke ekliptikalne suradnice Slnka:
    xg_eclS = -RE * np.cos(BE) * np.cos(LE);
    yg_eclS = -RE * np.cos(BE) * np.sin(LE);
    zg_eclS = -RE * np.sin(BE);
    
    #Geocentricke rovnikove suradnice Slnka:
    xg_S = xg_eclS;
    yg_S = yg_eclS * np.cos(ecl) - zg_eclS * np.sin(ecl);
    zg_S = yg_eclS * np.sin(ecl) + zg_eclS * np.cos(ecl);
    
    RA  = np.arctan2( yg_S, xg_S );        #rektascenzia a deklinacia
    Dec = np.arctan2( zg_S, np.sqrt(xg_S*xg_S + yg_S*yg_S) );
    
    FIS     = Dec ;         #geocentricka sirka,
    lambS = RA - GST;     #geocentricka dlzka vzhladom na Greenwich,
    
    XS = RE * np.cos(FIS)*np.cos(lambS); 
    YS = RE * np.cos(FIS)*np.sin(lambS);
    ZS = RE * np.sin(FIS);
    
    RS  = np.sqrt(XS*XS + YS*YS + ZS*ZS);   #velkost geocentrickeho vektora Slnka,
    R0S = np.array([XS/RS, YS/RS, ZS/RS]);         #jednotkovy geocentricky vektor Slnka,
    
    #%POLOHA MESIACA:
    #Stredna dlzka mesiaca:
    LL = 218.3164591*deg + 481267.88134236*deg*JS - 0.0013268*deg*JS**2 + JS**3/  (538841*deg)  - JS**4/ (65194000*deg);
    #Elongacia Mesiaca: 
    D  = 297.8502042*deg + 445267.11151680*deg*JS - 0.0016300*deg*JS**2 + JS**3/  (545868*deg)  - JS**4/(113065000*deg);
    #Stredna anomalia Slnka:
    M  = 357.5291092*deg +  35999.05029090*deg*JS - 0.0001536*deg*JS**2 + JS**3/(24490000*deg);
    #Stredna anomalia Mesiaca:
    MM = 134.9634114*deg + 477198.86763130*deg*JS + 0.0089970*deg*JS**2 + JS**3/   (69699*deg)  - JS**4/ (14712000*deg);
    #Argument sirky Mesiaca:
    F  =  93.2720993*deg + 483202.01752730*deg*JS - 0.0034029*deg*JS**2 - JS**3/  (3526000*deg) + JS**4/(863310000*deg);
    
    #Doplnkove argumenty:
    A1 = 119.75*deg + 131.849*deg*JS;
    A2 = 53.09*deg + 479264.290*deg*JS;
    A3 = 313.45*deg + 481266.484*deg*JS;
    
    #Excentricita drahy:
    E  = 1 - 0.002516*JS - 0.0000074*JS**2;
    
    #Periodicke cleny pre ekliptikalnu dlzku, sirku a sprievodic:
    #         D   M   MM  F    Sum(L)   Sum(R)    D   M   MM  F    Sum(B)
    PT = np.array([[1,0,0,1,0,6288774,-20905355,0,0,0,1,5128122],
    	[2,2,0,-1,0,1274027,-3699111,0,0,1,1,280602],
    	[3,2,0,0,0,658314,-2955968,0,0,1,-1,277693],
    	[4,0,0,2,0,213618,-569925,2,0,0,-1,173237],
    	[5,0,1,0,0,-185116,48888,2,0,-1,1,55413],
    	[6,0,0,0,2,-114332,-3149,2,0,-1,-1,46271],
    	[7,2,0,-2,0,58793,246158,2,0,0,1,32573],
    	[8,2,-1,-1,0,57066,-152138,0,0,2,1,17198],
    	[9,2,0,1,0,53322,-170733,2,0,1,-1,9266],
    	[10,2,-1,0,0,45758,-204586,0,0,2,-1,8822],
    	[11,0,1,-1,0,-40923,-129620,2,-1,0,-1,8216],
    	[12,1,0,0,0,-34720,108743,2,0,-2,-1,4324],
    	[13,0,1,1,0,-30383,104755,2,0,1,1,4200],
    	[14,2,0,0,-2,15327,10321,2,1,0,-1,-3359],
    	[15,0,0,1,2,-12528,0,2,-1,-1,1,2463],
    	[16,0,0,1,-2,10980,79661,2,-1,0,1,2211],
    	[17,4,0,-1,0,10675,-34782,2,-1,-1,-1,2065],
    	[18,0,0,3,0,10034,-23210,0,1,-1,-1,-1870],
    	[19,4,0,-2,0,8548,-21636,4,0,-1,-1,1828],
    	[20,2,1,-1,0,-7888,24208,0,1,0,1,-1794],
    	[21,2,1,0,0,-6766,30824,0,0,0,3,-1749],
    	[22,1,0,-1,0,-5163,-8379,0,1,-1,1,-1565],
    	[23,1,1,0,0,4987,-16675,1,0,0,1,-1491],
    	[24,2,-1,1,0,4036,-12831,0,1,1,1,-1475],
    	[25,2,0,2,0,3994,-10445,0,1,1,-1,-1410],
    	[26,4,0,0,0,3861,-11650,0,1,0,-1,-1344],
    	[27,2,0,-3,0,3665,14403,1,0,0,-1,-1335],
    	[28,0,1,-2,0,-2689,-7003,0,0,3,1,1107],
    	[29,2,0,-1,2,-2602,0,4,0,0,-1,1021],
    	[30,2,-1,-2,0,2390,10056,4,0,-1,1,833],
    	[31,1,0,1,0,-2348,6322,0,0,1,-3,777],
    	[32,2,-2,0,0,2236,-9884,4,0,-2,1,671],
    	[33,0,1,2,0,-2120,5751,2,0,0,-3,607],
    	[34,0,2,0,0,-2069,0,2,0,2,-1,596],
    	[35,2,-2,-1,0,2048,-4950,2,-1,1,-1,491],
    	[36,2,0,1,-2,-1773,4130,2,0,-2,1,-451],
    	[37,2,0,0,2,-1595,0,0,0,3,-1,439],
    	[38,4,-1,-1,0,1215,-3958,2,0,2,1,422],
    	[39,0,0,2,2,-1110,0,2,0,-3,-1,421],
    	[40,3,0,-1,0,-892,3258,2,1,-1,1,-366],
    	[41,2,1,1,0,-810,2616,2,1,0,1,-351],
    	[42,4,-1,-2,0,759,-1897,4,0,0,1,331],
    	[43,0,2,-1,0,-713,-2117,2,-1,1,1,315],
    	[44,2,2,-1,0,-700,2354,2,-2,0,-1,302],
    	[45,2,1,-2,0,691,0,0,0,1,3,-283],
    	[46,2,-1,0,-2,596,0,2,1,1,-1,-229],
    	[47,4,0,1,0,549,-1423,1,1,0,-1,223],
    	[48,0,0,4,0,537,-1117,1,1,0,1,223],
    	[49,4,-1,0,0,520,-1571,0,1,-2,-1,-220],
    	[50,1,0,-2,0,-487,-1739,2,1,-1,-1,-220],
    	[51,2,1,0,-2,-399,0,1,0,1,1,-185],
    	[52,0,0,2,-2,-381,-4421,2,-1,-2,-1,181],
    	[53,1,1,1,0,351,0,0,1,2,1,-177],
    	[54,3,0,-2,0,-340,0,4,0,-2,-1,176],
    	[55,4,0,-3,0,330,0,4,-1,-1,-1,166],
    	[56,2,-1,2,0,327,0,1,0,1,-1,-164],
    	[57,0,2,1,0,-323,1165,4,0,1,-1,132],
    	[58,1,1,-1,0,299,0,1,0,-1,-1,-119],
    	[59,2,0,3,0,294,0,4,-1,0,-1,115],
    	[60,2,0,-1,-2,0,8752,2,-2,0,1,107]]);
    
    #Suma L a R:
    temp = PT[np.logical_or(PT[:,2] == 1,PT[:,2] == -1)].T
    suma_1L = np.sum(temp[5]*np.sin(temp[1]*D + temp[2]*M*E + temp[3]*MM + temp[4]*F))
    suma_1R = np.sum(temp[6]*np.cos(temp[1]*D + temp[2]*M*E + temp[3]*MM + temp[4]*F))
    temp = PT[np.logical_or(PT[:,2] == 2,PT[:,2] == -2)].T
    suma_2L = np.sum(temp[5]*np.sin(temp[1]*D + temp[2]*M*E*E + temp[3]*MM + temp[4]*F))
    suma_2R = np.sum(temp[6]*np.cos(temp[1]*D + temp[2]*M*E*E + temp[3]*MM + temp[4]*F))
    temp = PT[PT[:,2] == 0].T
    suma_0L = np.sum(temp[5]*np.sin(temp[1]*D + temp[2]*M + temp[3]*MM + temp[4]*F))
    suma_0R = np.sum(temp[6]*np.cos(temp[1]*D + temp[2]*M + temp[3]*MM + temp[4]*F))
    suma_L = suma_0L + suma_1L + suma_2L
    suma_R = suma_0R + suma_1R + suma_2R
    
    #Suma B:
    temp = PT[np.logical_or(PT[:,8] == 1,PT[:,8] == -1)].T
    suma_1B = np.sum(temp[11]*np.sin(temp[7]*D + temp[8]*M*E + temp[9]*MM + temp[10]*F))
    temp = PT[np.logical_or(PT[:,8] == 2,PT[:,8] == -2)].T
    suma_2B = np.sum(temp[11]*np.sin(temp[7]*D + temp[8]*M*E*E + temp[9]*MM + temp[10]*F))
    temp = PT[PT[:,8] == 0].T
    suma_0B = np.sum(temp[11]*np.sin(temp[7]*D + temp[8]*M + temp[9]*MM + temp[10]*F))
    suma_B = suma_0B + suma_1B + suma_2B
    
    #% Geocentricke ekliptikalne suradnice Mesiaca:
    SUMA_L = suma_L + 3958*np.sin(A1) + 1962*np.sin(LL-F) + 318*np.sin(A2);
    SUMA_B = suma_B - 2235*np.sin(LL) + 382*np.sin(A3) + 175*np.sin(A1 - F) + 175*np.sin(A1+F) + 127*np.sin(LL-MM) - 115*np.sin(LL+MM);
    
    LM = LL + (SUMA_L/1000000)*deg;
    BM = (SUMA_B/1000000)*deg;
    RM = (385000.56 + suma_R/1000)*1000;
    
    xg_eclM = RM * np.cos(BM) * np.cos(LM);
    yg_eclM = RM * np.cos(BM) * np.sin(LM);
    zg_eclM = RM * np.sin(BM);
    
    #Geocentricke rovnikove suradnice Mesiaca:
    xg_M = xg_eclM;
    yg_M = yg_eclM * np.cos(ecl) - zg_eclM * np.sin(ecl);
    zg_M = yg_eclM * np.sin(ecl) + zg_eclM * np.cos(ecl);
    
    RA  = np.arctan2( yg_M, xg_M );        #rektascenzia a deklinacia
    Dec = np.arctan2( zg_M, np.sqrt(xg_M*xg_M + yg_M*yg_M) );
    
    FIM     = Dec ;         #geocentricka sirka,
    lambM = RA - GST;     #geocentrick� dlzka vzhladom na Greenwich,
    lambM = lambM + 2*np.pi if lambM < 2*np.pi else lambM
    
    XM = RM * np.cos(FIM)*np.cos(lambM); 
    YM = RM * np.cos(FIM)*np.sin(lambM);
    ZM = RM * np.sin(FIM);
    
    RM  = np.sqrt(XM*XM + YM*YM + ZM*ZM);   #velkost geocentrickeho vektora Mesiaca,
    R0M = np.array([XM/RM, YM/RM, ZM/RM]);         #jednotkovy geocentricky vektor Mesiaca,
       
    #%KROK 1: KOREKCIE POCITANE V  CASOVEJ OBLASTI (IERS2003, str.79)   
    #A)IN-PHASE:
    h_0 = 0.6078; 
    h_2 = -0.0006; 
    l_0 = 0.0847; 
    l_2 = 0.0002;
    h_3 = 0.292 ; 
    l_3 = 0.015;    
    h2  = h_0 + h_2*((3 * np.sin(fi)*np.sin(fi) - 1)/2);
    l2  = l_0 + l_2*((3 * np.sin(fi)*np.sin(fi) - 1)/2);
    
    #Mesiac pre n = 2:
    A2         = h2 * r0 * (3/2*np.dot(R0M,r0)**2 -1/2);
    B2         = l2 *  3 * np.dot(R0M,r0) * (R0M - np.dot(R0M,r0) * r0);
    delta_r_M2 = (GMm * Re**4)/(GMe * RM**3) * (A2 + B2);
    #Slnko pre  n = 2:
    A2         = h2 * r0 * (3/2*np.dot(R0S,r0)**2 -1/2);
    B2         = 3* l2 * np.dot(R0S,r0) * (R0S - np.dot(R0S,r0) * r0);
    delta_r_S2 = (GMs * Re**4)/(GMe * RS**3) * (A2 + B2);
    
    #Mesiac pre n = 3:
    A3         = h_3 * r0 * (5/2 * np.dot(R0M,r0)**3 -3/2*np.dot(R0M,r0));
    B3         = l_3 * (15/2 * np.dot(R0M,r0)**2 - 3/2)*(R0M - np.dot(R0M,r0) * r0);
    delta_r_M3 = (GMm * Re**5)/(GMe * RM**4) * (A3 + B3);
    #Slnko pre  n = 3 => prispevok zanedbatelny
    
    delta_r_M = delta_r_M2 + delta_r_M3;
    delta_r_S = delta_r_S2;
    
    ner_inphase_M = Rot @ delta_r_M;
    ner_inphase_S = Rot @ delta_r_S;
    
    n_inphase_M = ner_inphase_M[0];
    e_inphase_M = ner_inphase_M[1];
    r_inphase_M = ner_inphase_M[2];
    
    n_inphase_S = ner_inphase_S[0];   #zlozka v smere sever;
    e_inphase_S = ner_inphase_S[1];   #zlozka v smere vychod;
    r_inphase_S = ner_inphase_S[2];   #zlozka v radialnom smere;        
    
    #B)PRISPEVOK K PRIECNEJ ZLOZKE:     
    #MESIAC:
    P21M = (3 * XM * ZM)/(RM**2 * np.cos(lambM));
    P22M = (3 * (XM**2 - YM**2))/(RM**2 * np.cos (2*lambM));
    #denne slapy:
    n_cont_d = -0.0012 * np.sin(fi) * (GMm * Re**4)/(GMe * RM**3)* P21M * np.sin(fi) * np.cos(lamb - lambM);
    e_cont_d =  0.0012 * np.sin(fi) * (GMm * Re**4)/(GMe * RM**3)* P21M * np.cos(2*fi) * np.sin(lamb - lambM);
    #poldenne slapy: 
    n_cont_p = -1/2*0.0024 * np.sin(fi) * np.cos(fi) * (GMm * Re**4)/(GMe * RM**3)* P22M * np.cos(2*(lamb - lambM));
    e_cont_p = -1/2*0.0024 * np.sin(fi) * np.cos(fi) * (GMm * Re**4)/(GMe * RM**3)* P22M * np.sin(fi) * np.sin(2*(lamb - lambM));
    
    n_con_M = n_cont_d + n_cont_p;
    e_con_M = e_cont_d + e_cont_p;
    
    #SLNKO:
    P21S = (3 * XS * ZS)/(RS**2 * np.cos(lambS));
    P22S = (3 * (XS**2 - YS**2))/(RS**2 * np.cos (2*lambS));
    #denne slapy:
    n_cont_d = -0.0012 * np.sin(fi) * (GMs * Re**4)/(GMe * RS**3)* P21S * np.sin(fi) * np.cos(lamb - lambS);
    e_cont_d =  0.0012 * np.sin(fi) * (GMs * Re**4)/(GMe * RS**3)* P21S * np.cos(2*fi) * np.sin(lamb - lambS);
    #poldenne slapy: 
    n_cont_p = -1/2*0.0024 * np.sin(fi) * np.cos(fi) * (GMs * Re**4)/(GMe * RS**3)* P22S * np.cos(2*(lamb - lambS));
    e_cont_p = -1/2*0.0024 * np.sin(fi) * np.cos(fi) * (GMs * Re**4)/(GMe * RS**3)* P22S * np.sin(fi) * np.sin(2*(lamb - lambS));
    
    n_con_S = n_cont_d + n_cont_p;
    e_con_S = e_cont_d + e_cont_p;
    
    #Celkovy prispevok casti B) Mesiac + Slnko: 
    n_con = n_con_M + n_con_S;
    e_con = e_con_M + e_con_S;
    
    #C)OUT-OF-PHASE:  
    #MESIAC:
    #denne slapy:
    r_oop_d = -3/4*(-0.0025)*(GMm * Re**4)/(GMe * RM**3)*np.sin(2*FIM)*np.sin(2*fi)*np.sin(lamb - lambM);
    n_oop_d = -3/2*(-0.0007)*(GMm * Re**4)/(GMe * RM**3)*np.sin(2*FIM)*np.cos(2*fi)*np.sin(lamb - lambM);
    e_oop_d = -3/2*(-0.0007)*(GMm * Re**4)/(GMe * RM**3)*np.sin(2*FIM)*np.sin(fi)  *np.cos(lamb - lambM);
    #poldenne slapy: 
    r_oop_p = -3/4*(-0.0022)*(GMm * Re**4)/(GMe * RM**3)*np.cos(FIM)**2*np.cos(fi)**2*np.sin(2*(lamb - lambM));
    n_oop_p =  3/4*(-0.0007)*(GMm * Re**4)/(GMe * RM**3)*np.cos(FIM)**2*np.sin(2*fi)*np.sin(2*(lamb - lambM));
    e_oop_p = -3/4*(-0.0007)*(GMm * Re**4)/(GMe * RM**3)*np.cos(FIM)**2*2*np.cos(fi)*np.cos(2*(lamb - lambM));
    
    r_oop_M = (r_oop_d + r_oop_p);
    n_oop_M = (n_oop_d + n_oop_p);
    e_oop_M = (e_oop_d + e_oop_p);
    
    #SLNKO:
    #denne slapy:
    r_oop_d = -3/4*(-0.0025)*(GMs * Re**4)/(GMe * RS**3)*np.sin(2*FIS)*np.sin(2*fi)*np.sin(lamb - lambS);
    n_oop_d = -3/2*(-0.0007)*(GMs * Re**4)/(GMe * RS**3)*np.sin(2*FIS)*np.cos(2*fi)*np.sin(lamb - lambS);
    e_oop_d = -3/2*(-0.0007)*(GMs * Re**4)/(GMe * RS**3)*np.sin(2*FIS)*np.sin(fi)  *np.cos(lamb - lambS);
    #poldenne slapy: 
    r_oop_p = -3/4*(-0.0022)*(GMs * Re**4)/(GMe * RS**3)*np.cos(FIS)**2*np.cos(fi)**2*np.sin(2*(lamb - lambS));
    n_oop_p =  3/4*(-0.0007)*(GMs * Re**4)/(GMe * RS**3)*np.cos(FIS)**2*np.sin(2*fi)*np.sin(2*(lamb - lambS));
    e_oop_p = -3/4*(-0.0007)*(GMs * Re**4)/(GMe * RS**3)*np.cos(FIS)**2*2*np.cos(fi)*np.cos(2*(lamb - lambS));
    
    r_oop_S = (r_oop_d + r_oop_p);
    n_oop_S = (n_oop_d + n_oop_p);
    e_oop_S = (e_oop_d + e_oop_p);    
    
    #Celkovy prispevok casti C) Mesiac + Slnko: 
    r_oop = r_oop_M + r_oop_S;
    n_oop = n_oop_M + n_oop_S;
    e_oop = e_oop_M + e_oop_S;   
    
    #%KROK 2: KOREKCIE POCITANE VO FREKVENCNEJ OBLASTI (IERS2003, str.79)
    #Delaunayove premenne:
    l     = 134.96340251*deg + 1717915923.217800*sec*JS + 31.879200*sec*JS**2 + 0.05163500*sec*JS**3 - 0.0002447000*sec*JS**4;
    ll    = 357.52910918*deg +  129596581.048100*sec*JS -  0.553200*sec*JS**2 + 0.00013600*sec*JS**3 - 0.0000114900*sec*JS**4;
    F     =  93.27209062*deg + 1739527262.847800*sec*JS - 12.751200*sec*JS**2 - 0.00103700*sec*JS**3 + 0.0000041700*sec*JS**4;
    D     = 297.85019547*deg + 1602961601.209000*sec*JS -  6.370600*sec*JS**2 + 0.00659300*sec*JS**3 - 0.0000316900*sec*JS**4;
    Omega = 125.04455501*deg -    6962890.543100*sec*JS +  7.472200*sec*JS**2 + 0.00770200*sec*JS**3 - 0.0000593900*sec*JS**4;    
    
    #A)DENNE SLAPY - teseralne:
    #Tab.7.5a str.82 m = 1
    A = np.array([
    [13.39866 ,135.655 ,1 ,-2 ,0 ,1 ,0 ,0 ,1 ,0 ,2 ,0 ,2 ,-0.08 ,0.00 ,-0.01 ,0.01],
    [13.94083 ,145.545 ,1 ,-1 ,0 ,0 ,-1 ,0 ,0 ,0 ,2 ,0 ,1 ,-0.10 ,0.00 ,0.00 ,0.00],
    [13.94303 ,145.555 ,1 ,-1 ,0 ,0 ,0 ,0 ,0 ,0 ,2 ,0 ,2 ,-0.51 ,0.00 ,-0.02 ,0.03],
    [14.49669 ,155.655 ,1 ,0 ,0 ,1 ,0 ,0 ,1 ,0 ,0 ,0 ,0 ,0.06 ,0.00 ,0.00 ,0.00],
    [14.91787 ,162.556 ,1 ,1 ,-3 ,0 ,0 ,1 ,0 ,1 ,2 ,-2 ,2 ,-0.06 ,0.00 ,0.00 ,0.00],
    [14.95893 ,163.555 ,1 ,1 ,-2 ,0 ,0 ,0 ,0 ,0 ,2 ,-2 ,2 ,-1.23 ,-0.07 ,0.06 ,0.01],
    [15.03886 ,165.545 ,1 ,1 ,0 ,0 ,-1 ,0 ,0 ,0 ,0 ,0 ,-1 ,-0.22 ,0.01 ,0.01 ,0.00],
    [15.04107 ,165.555 ,1 ,1 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,12.00 ,-0.78 ,-0.67 ,-0.03],
    [15.04328 ,165.565 ,1 ,1 ,0 ,0 ,1 ,0 ,0 ,0 ,0 ,0 ,1 ,1.73 ,-0.12 ,-0.10 ,0.00],
    [15.08214 ,166.554 ,1 ,1 ,1 ,0 ,0 ,-1 ,0 ,-1 ,0 ,0 ,0 ,-0.50 ,-0.01 ,0.03 ,0.00],
    [15.12321 ,167.555 ,1 ,1 ,2 ,0 ,0 ,0 ,0 ,0 ,-2 ,2 ,-2 ,-0.11 ,0.01 ,0.01 ,0.00]]);
    
    NF = l*A[:,8] + ll*A[:,9] + F*A[:,10] + D*A[:,11] + Omega*A[:,12];
    Theta_f = GST + np.pi - NF;
    #Radialna zlozka:
    r_teseralne = np.sum((A[:,13]*np.sin(Theta_f +lamb) + A[:,14]*np.cos(Theta_f + lamb))*np.sin(2*fi))/1000;
    #Transverzalne zlozky:   
    n_teseralne = np.sum((A[:,15]*np.sin(Theta_f +lamb) + A[:,16]*np.cos(Theta_f + lamb))*np.cos(2*fi))/1000; 
    e_teseralne = np.sum((A[:,15]*np.cos(Theta_f +lamb) - A[:,16]*np.sin(Theta_f + lamb))*np.sin(fi))/1000;
    
    #B)ZONALNE SLAPY:
    #Tab.7.5b str.82
    B = np.array([
    [0.00221,55.565,0,0,0,0,1,0,0,0,0,0,1,0.47,0.16,0.23,0.07],
    [0.08214,57.555,0,0,2,0,0,0,0,0,-2,2,-2,-0.20,-0.11,-0.12,-0.05],
    [0.54438,65.455,0,1,0,-1,0,0,-1,0,0,0,0,-0.11,-0.09,-0.08,-0.04],
    [1.09804,75.555,0,2,0,0,0,0,0,0,-2,0,-2,-0.13,-0.15,-0.11,-0.07],
    [1.10024,75.565,0,2,0,0,1,0,0,0,-2,0,-1,-0.05,-0.06,-0.05,-0.03]]);
    
    NF = l*B[:,8] + ll*B[:,9] + F*B[:,10] + D*B[:,11] + Omega*B[:,12];
    Theta_f = - NF;
    #Radialna zlozka:
    r_zonalne = np.sum((3/2 * np.sin(fi)**2 - 1/2)*(B[:,13]*np.cos(Theta_f) + B[:,14]*np.sin(Theta_f)))/1000;
    #Transverzalne zlozky:
    n_zonalne = np.sum((B[:,15]*np.cos(Theta_f) + B[:,16]*np.sin(Theta_f))*np.sin(2*fi))/1000;
    
    #%% VYSLEDNY POHYB STANOVISKA pre "TIDE FREE MODEL" - lokalny suradnicovy system:
    n_tf = (n_inphase_M + n_inphase_S + n_con + n_oop + n_teseralne + n_zonalne)*1000;
    e_tf = (e_inphase_M + e_inphase_S + e_con + e_oop + e_teseralne)*1000;
    r_tf = (r_inphase_M + r_inphase_S + r_oop + r_teseralne + r_zonalne)*1000;        
    
    #VYSLEDNY POHYB STANOVISKA pre "TIDE FREE MODEL" - v geocentrickom terestrickom ss (vzhladom na Greenwich):
    dXYZt_tf = Rot.T @ np.array([n_tf, e_tf, r_tf]);  # v mm
    ner_tf = np.array([n_tf,e_tf,r_tf]);
    
    #Permanentne slapy = deformacia zemskeho povrchu
    #/IERS 2003, Technical Note, str.83/
    
    #Permanentna deformacia v smere sprievodica:
    r_ps = (-0.1206 + 0.0001 * ((3 * np.sin(fi)**2 -1)/2))*((3 * np.sin(fi) * np.sin(fi) -1)/2)*1000;
    
    #Permanentna deformacia v smere sever-juh:
    n_ps = (-0.0252 - 0.0001 * (3 * np.sin(fi)**2 -1)/2 * np.sin(2 * fi))*1000;
    
    #VYSLEDNY POHYB STANOVISKA pre "MEAN TIDE MODEL" - lokalny suradnicovy system:
    n_mt = n_tf + n_ps;
    e_mt = e_tf;
    r_mt = r_tf + r_ps;
    
    #VYSLEDNY POHYB STANOVISKA pre "MEAN TIDE MODEL" - v geocentrickom terestrickom ss (vzhladom na Greenwich):
    dXYZt_mt = Rot.T @ np.array([n_mt, e_mt, r_mt]);  # v mm
    ner_mt = np.array([n_mt,e_mt,r_mt]); 
    
    permanent = np.array([n_ps,0,r_ps]);
    
    #fi = -90:0.1:90;
    #fi = fi*pi/180;
    #P2 = (3*np.sin(fi)**2 - 1)/2;
    #n_ps = (-0.0252 - 0.0001 .* P2 .* np.sin(2*fi))*1000;
    #e_ps = zeros(numel(fi),1);
    #r_ps = ((-0.1206 + 0.0001 * P2) .* P2)*1000; 
    
    xyzSET = dXYZt_tf/1e3; # in [m]
    
    return xyzSET
