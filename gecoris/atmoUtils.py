#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""atmoUtils

Module containing atmospheric modelling utility functions

Copyright (C) 2021 by R.Czikhardt

Email: czikhardt.richard@gmail.com
Last edit: 1.10.2020

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

import os
import re
import requests
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d, interpn, RegularGridInterpolator
from scipy.integrate import cumtrapz
from gecoris import geoUtils, radarUtils


def getERA5(outDir,aoi,year,month,day,time):
    """API form to download ERA5 ECMWF data
    
    for slant delay computation
    
    input: year, month, day, time, outDir, area (N,W,S,E)
    output: ERA5 netCDF saved in outDir
    """
    import cdsapi
        
    outDir += (year + month + day)
    #outDir = 'D:/DELFT/AtmoDelay_experiment/temp.netcdf'
    #aoi = '53.50/6.35/53.00/7.15'
    #year = '2019'
    #month = '03'
    #day = '04'
    #time = '04:00'
    
    ##
    c = cdsapi.Client()
    
    c.retrieve(
        'reanalysis-era5-pressure-levels',
        {
            'product_type':'reanalysis',
            'format':'netcdf',
            'area':aoi,
            'variable':[
                'geopotential','relative_humidity','temperature'
            ],
            'pressure_level':[
                '1','2','3',
                '5','7','10',
                '20','30','50',
                '70','100','125',
                '150','175','200',
                '225','250','300',
                '350','400','450',
                '500','550','600',
                '650','700','750',
                '775','800','825',
                '850','875','900',
                '925','950','975',
                '1000'
            ],
            'year':year,
            'month':month,
            'day':day,
            'time':time
        },
        outDir)

def prepareAtmo(stack,outDir):
    """Prepare ERA5 ECMWF data for stack
    
    input: stack instance, outDir
    output: ERA5 netCDF saved in outDir
    """
    
    # constants:
    hMax = 30e3 # [m]
    hIncr = 20 # [m] height integration step
    R = 6371e3
    
    acqDates = stack.acqDates
    
    # get scene center and sat. position to get required AOI crop:
    metadata = stack.metadata[acqDates[-1]]
    centerPLH = np.array([metadata['centerLat']/180*np.pi,
                         metadata['centerLon']/180*np.pi,
                         metadata['centerH']])
    centerXYZ = geoUtils.plh2xyz(centerPLH)
    _,_,satVec = radarUtils.xyz2t(centerXYZ,metadata) 
    LOS = satVec - centerXYZ
    LOS = LOS/np.linalg.norm(LOS)
    incAngle = geoUtils.xyz2zas(centerXYZ,satVec)[0]
    slantMax = hMax/np.cos(incAngle)
    slantDist = np.arange(0,slantMax,hIncr/np.cos(incAngle))
    xyz_int = np.array([centerXYZ]).T + LOS.T*np.tile(slantDist,(3,1))
    plh_int = geoUtils.xyz2plh(xyz_int.T)
    # scene size:
    latSize = metadata['nAzimuth']*metadata['nBursts']*metadata['azimuthSpacing']/R
    lonSize = metadata['nRange']*metadata['rangeSpacing']/R
    offset = np.ceil(max(lonSize*180/np.pi,latSize*180/np.pi))
    # aoi:
    minLon = np.min(plh_int[:,1])*180/np.pi -offset # [deg]
    maxLon = np.max(plh_int[:,1])*180/np.pi +offset
    minLat = np.min(plh_int[:,0])*180/np.pi -offset
    maxLat = np.max(plh_int[:,0])*180/np.pi +offset
    aoi = str(maxLat) +'/'+ str(minLon) +'/'+ str(minLat) +'/'+ str(maxLon)
    # acq. time same for all dates, round to hour:
    time = geoUtils.roundTime(metadata['acqDate'],3600).strftime('%H:%M')
    
    # loop through acqDates
    for date in acqDates:
        year = date[:4]
        month = date[4:6]
        day = date[6:]        
        netcdf_fname = outDir + date    
        if not os.path.exists(netcdf_fname):
            getERA5(outDir,aoi,year,month,day,time)


def getTropoDelay(xyz_P,xyz_S,date,atmoDir,method='Jolviet'):
    
    from netCDF4 import Dataset
    
    # Constants:
    # Parameters from Saastaamoinen and Jolivet et al 2014, GRL
    Rd = 287.05            # [j/kg/K] specific gas constants dry air
    Rv = 461.495           # [j/kg/K] water vapour
    k1 = 0.776             # [K/Pa]
    k2 = 0.716             # [K/Pa]
    k3 = 3.75e3            # [K^2/Pa]
    # constants:
    hStn = 218 # [m]
    hMax = 30e3 # [m]
    hIncr = 20 # height integration step
    
    #LOS vector
    LOS = xyz_S-xyz_P
    LOS = LOS/np.linalg.norm(LOS)
    
    # compute integration points (P_int) position up to h_max
    incAngle = geoUtils.xyz2zas(xyz_P,xyz_S)[0]
    slantMax = hMax/np.cos(incAngle)
    slantDist = np.arange(0,slantMax,hIncr/np.cos(incAngle))
    xyz_int = np.array([xyz_P]).T + LOS.T*np.tile(slantDist,(3,1))
    plh_int = geoUtils.xyz2plh(xyz_int.T)
    
    #% load netcdf
    netcdf_fname = atmoDir + date
    
    ds = Dataset(netcdf_fname)
    # parse variables:
    lon = ds['longitude'][:].filled()
    lat = ds['latitude'][:].filled()
    Plevel = ds['level'][:].filled()
    Temp = ds['t'][:].filled().squeeze().T # re-order
    Geopot = ds['z'][:].filled().squeeze().T
    Humid = ds['r'][:].filled().squeeze().T
    # prepare meshgrids:
    latGrid, lonGrid, Pressure = np.meshgrid(lat, lon, Plevel)
    lonGrid = lonGrid/180*np.pi
    latGrid = latGrid/180*np.pi
    
    Pressure = np.tile(Plevel,(lon.shape[0],lat.shape[0],1))
    Pressure *= 100 # from [hPa] to [Pa]
    
    #% convert humidity to WaterVapour
    # a.) water vapour saturation pressure:
    # C1 = -7.85951783
    # C2 = 1.84408259
    # C3 = -11.7866497
    # C4 =  22.6807411
    # C5 = -15.9618719
    # C6 = 1.80122502
    # Pc = 220640 # [hPa]
    # Tc = 647.096 # [K]
    # theta = 1 - Temp/Tc
    # tmp = Tc/Temp*(C1*theta + C2*theta**1.5 + C3*theta**3 + 
    #                C4*theta**3.5 + C5*theta**4 + C6*theta**7.5)
    # P_WS = Pc*np.exp(tmp) # [hPa]
    # # b.) partial water vapour pressure from humidity:
    # e = Humid/100*P_WS
    
    #% other variant:
    svpw = 6.1121*np.exp((17.502*(Temp-273.16))/(240.97+Temp-273.16))
    svpi = 6.1121*np.exp((22.587*(Temp-273.16))/(273.86+Temp-273.16))
    tempbound1 = 273.16
    tempbound2 = 250.16
    wgt = (Temp - tempbound2)/(tempbound1 - tempbound2)
    svp = svpi + (svpw - svpi)*wgt**2
    svp[Temp > tempbound1] = svpw[Temp > tempbound1]
    svp[Temp < tempbound2] = svpi[Temp < tempbound2]
    e = Humid/100*svp;
    e *= 100 # from [hPa] to [Pa]
    
    #%
    #fig = plt.figure()
    #ax = fig.add_subplot(111, projection='3d')
    #scat = ax.scatter(lonGrid, latGrid, Pressure, c=Temp.flatten(), alpha=0.5)
    #fig.colorbar(scat, shrink=0.5, aspect=5)
    
    #% prepare height profiles:
    # Convert Geopotential to Geopotential Height and then to Geometric Height
    g0 = 9.80665
    H = Geopot/g0
    # map of g with latitude
    g = 9.80616*(1 - 0.002637*np.cos(2*latGrid) + 0.0000059*(np.cos(2*latGrid))**2)
    # map of Re with latitude
    Rmax = 6378137;
    Rmin = 6356752;
    Re = np.sqrt(1/(((np.cos(latGrid)**2)/Rmax**2) + ((np.sin(latGrid)**2)/Rmin**2)))
    # Calculate Geometric Height, Z
    Z = (H*Re)/(g/g0*Re - H)
    
    #% spline-interpolate in height profiles:       
    hInt = np.flipud(np.arange(hStn,hMax,hIncr))
    
    def spline1d(z,var,z_new):
        f = interp1d(z, var, kind='cubic',bounds_error=False, 
                     fill_value=np.nan)
        return f(z_new)
    
    #e_int = np.apply_along_axis(spline1d,2,Z,e)
    e_int = np.empty(Z.shape[:2]+ hInt.shape)
    P_int = np.empty(Z.shape[:2]+ hInt.shape)
    T_int = np.empty(Z.shape[:2]+ hInt.shape)
    for i,j in np.ndindex(Z.shape[:2]):
        e_int[i,j,:] = spline1d(Z[i,j,:], e[i,j,:], hInt)
        P_int[i,j,:] = spline1d(Z[i,j,:], Pressure[i,j,:], hInt)
        T_int[i,j,:] = spline1d(Z[i,j,:], Temp[i,j,:], hInt)
        
    #% trilinear interpolation   
    lonVec = lon/180*np.pi
    latVec = np.flip(lat/180*np.pi)
    hVec = np.flip(hInt,0)    
    lph_int = np.array([plh_int[:,1],plh_int[:,0],plh_int[:,2]]).T
    
    #
    e_P = interpn((lonVec, latVec, hVec),np.flip(e_int,(0,2)),lph_int,
                  bounds_error=False,fill_value=np.nan)
    T_P = interpn((lonVec, latVec, hVec),np.flip(T_int,(0,2)),lph_int,
                  bounds_error=False,fill_value=np.nan)
    P_P = interpn((lonVec, latVec, hVec),np.flip(P_int,(0,2)),lph_int,
                  bounds_error=False,fill_value=np.nan)
    
    # remove nan values outside (extrapol.)
    slantDist = slantDist[~np.isnan(e_P)]
    e_P = e_P[~np.isnan(e_P)]
    T_P = T_P[~np.isnan(T_P)]
    P_P = P_P[~np.isnan(P_P)]
    
    # integrate along LOS direction
    if method == 'default':
        tmp = k1*(P_P-e_P)/T_P + k2*e_P/T_P + k3*e_P/(T_P**2)
        LOS = np.flipud(cumtrapz(np.flipud(tmp),np.flipud(slantDist)))[0]*1e-6
    elif method == 'Jolviet':
        tmp = ((k2-(Rd*k1/Rv))*e_P/T_P + k3*e_P/(T_P**2))
        LOSw = np.flipud(cumtrapz(np.flipud(tmp),np.flipud(slantDist)))[0]*1e-6
        # g between h_P and h_max (using normal gradient):
        g_m = np.mean(g.flatten()) - 0.3086e-5*hMax/2
        # dry delay (using cosine projection...)
        LOSd = 1/np.cos(incAngle)*(1e-6)*((k1*Rd/g_m)*(P_P - P_P[0]))[-1]
        LOS = LOSw + LOSd
    
    return -LOS # TODO

# IONOSPHERE UTILS:
def parse_map(tecmap, exponent = -1):
    tecmap = re.split('.*END OF TEC MAP', tecmap)[0]
    return np.stack([np.fromstring(l, sep=' ') for l in re.split('.*LAT/LON1/LON2/DLON/H\\n',tecmap)[1:]])*10**exponent
    
def get_tecmaps(filename):
    with open(filename) as f:
        ionex = f.read()
        return [parse_map(t) for t in ionex.split('START OF TEC MAP')[1:]]

def get_tec(tecmap, lat, lon):
    i = round((87.5 - lat)*(tecmap.shape[0]-1)/(2*87.5))
    j = round((180 + lon)*(tecmap.shape[1]-1)/360)
    return tecmap[i,j]

def ionex_filename(year, day, centre, zipped = True):
    return '{}g{:03d}0.{:02d}i{}'.format(centre, day, year % 100, '.Z' if zipped else '')

def ionex_ftp_path(year, day, centre):
    return 'https://cddis.nasa.gov/archive/gnss/products/ionex/{:04d}/{:03d}/{}'.format(year, day, ionex_filename(year, day, centre))

def ionex_local_path(year, day, centre = 'cod', directory = '/tmp', zipped = False):
    return directory + '/' + ionex_filename(year, day, centre, zipped)
    
def download_ionex(year, day, centre = 'cod', output_dir = '/tmp'):
    #wget.download(ionex_ftp_path(year, day, centre), output_dir)
    r = requests.get(ionex_ftp_path(year, day, centre))
    rename = False
    if r.status_code == 404: # use CODE rapid solution if no final
        centre = 'cor'
        rename = True
        r = requests.get(ionex_ftp_path(year, day, centre))
    with open(ionex_local_path(year, day, zipped=True, 
                               directory=output_dir, centre=centre), 'wb') as fd:
        for chunk in r.iter_content(chunk_size=128):
            fd.write(chunk)
    subprocess.call(['gzip', '-d', 
                     ionex_local_path(year, day, zipped=True,
                                      directory=output_dir, centre=centre)])
    if rename:
        subprocess.call(['mv',ionex_local_path(year, day, directory=output_dir, centre=centre),
                              ionex_local_path(year, day, directory=output_dir, centre='cod')])

def getIPP(satVec,xyz):
    # model radius:
    RM = (6371 + 450)*1e3 # km
    LOS = satVec - xyz
    LOS = LOS/np.linalg.norm(LOS)    
    # distance from point to IPP (intersection of LOS and model sphere):
    d = -LOS@xyz + np.sqrt((LOS@xyz)**2 - (np.linalg.norm(xyz))**2 + RM**2)
    IPP = xyz + LOS*d
    return geoUtils.xyz2plh(IPP)

def ionexInterpolator(acqDate, atmoDir = '/tmp/'):
    # parse datetime acqDate:
    year = acqDate.year
    doy = acqDate.timetuple().tm_yday
    # download ionex:
    if not os.path.exists(ionex_local_path(year, doy, directory=atmoDir)):
        download_ionex(year, doy, output_dir=atmoDir)
    tecmaps = get_tecmaps(ionex_local_path(year, doy, directory=atmoDir))
    #% get interpolation grid:
    TECgrid = np.array(tecmaps[:-1])
    lat = np.arange(-87.5,88,2.5)
    lon = np.arange(-180,181,5)
    t = np.arange(24)
    interpolator = RegularGridInterpolator((t,lat,lon), TECgrid)
    return interpolator

def getIonoDelay(xyz, satVec, acqDate, atmoDir = '/tmp/'):
    RM = (6371 + 450)*1e3 # km
    # prepare interpolator:
    interp = ionexInterpolator(acqDate, atmoDir = atmoDir)
    # prepare IPP:
    plhIPP = getIPP(satVec,xyz)
    # acquisition time:
    hh = acqDate.hour+acqDate.minute/60+acqDate.second/3600
    # interpolate TEC:
    TEC_IPP = interp(np.array([hh,
                               plhIPP[0]*180/np.pi,
                               plhIPP[1]*180/np.pi]))
    # mapping function:
    z,_,_ = geoUtils.xyz2zas(xyz,np.array([satVec]))
    MF = 1/np.sqrt(1-(np.linalg.norm(xyz)/RM*np.sin(z))**2)
    # scaling factor for Sentinel-1:
    f = 5.405e9 # Hz
    scale = 0.9
    iono_delay = 40.3e16/(f**2)*TEC_IPP*MF*scale
    return iono_delay[0]


# preparation for VMF3
# for now use supplied:
def parse_supplied_tropodelay(acqDate, file):
    q = np.genfromtxt(file, delimiter=',')
    idx = np.where(q[:,0] == int(acqDate))[0]
    return q[idx,1]
