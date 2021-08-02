# -*- coding: utf-8 -*-
"""batchDownload

Routine to automatically download Sentinel-1 SLC data for specified AOI 
and orbit using aria2c or sentinelsat utilities

Copyright (C) 2021 by R.Czikhardt

Email: czikhardt.richard@gmail.com
Last edit: 25.8.2020

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
###############################################################################

print('Geodetic Corner Reflector (In)SAR Toolbox (GECORIS) v.1.0')
print('Copyright (c) 2021 Richard Czikhardt, czikhardt.richard@gmail.com')
print('Dept. of Theoretical Geodesy, Slovak University of Technology')
print('-----------------------------------------------------------------')
print('License: GNU GPL v3+')
print('-----------------------------------------------------------------')

import sys
import os
import subprocess
import time


def batchDownload(parms):    
    # parse parameters:
    downloadMethod = parms['parameters']['downloadMethod']
    dataPath = parms['parameters']['dataPath']
    user = parms['parameters']['scihubUser']
    password = parms['parameters']['scihubPassword']
    updateSwitch = parms['parameters']['updateSwitch']
    min_lon = str(parms['AOI']['minLongitude'])
    max_lon = str(parms['AOI']['maxLongitude'])
    min_lat = str(parms['AOI']['minLatitude'])
    max_lat = str(parms['AOI']['maxLatitude'])
    aoiFile = parms['AOI']['geojson']
    startDate = str(parms['AOI']['startDate'])
    endDate = str(parms['AOI']['endDate'])
    swaths = parms['AOI']['swaths']
    
    # set endDate to current date by default:
    if not endDate:
        endDate = time.strftime('%Y%m%d')
    
    # construct AOI polygon:
    aoiPolygon = ("polygon(("+min_lon+" "+min_lat+", "+min_lon+" "
                  +max_lat+", "+max_lon+" "+max_lat+", "+max_lon+" "
                  +min_lat+", "+min_lon+" "+min_lat+" ))")
    
    
    # loop through tracks
    for swath in swaths:
        if not os.path.exists(dataPath + swath):
            os.makedirs(dataPath + swath)
        os.chdir(dataPath + swath)
        # get orbit direction
        if swath[0:3] == 'ASC':
            orbit = 'Ascending'
            orbitNumber = swath[3:]
        elif swath[0:3] == 'DSC':
            orbit = 'Descending'
            orbitNumber = swath[3:]
        else:
            raise Exception("Some problem with orbit!")
        
        # Find latest downloaded SLC - set startDate accordingly (+1 day)
        if updateSwitch:
            startDate = str(int(max([slc[17:25] 
                                     for slc in os.listdir(dataPath + swath)
                                         if slc.endswith('.zip')])))
        
        if downloadMethod == 'sentinelsat':
            #%% download variant with sentinelsat
            parms = ("-u "+user+" -p "+password+" -g "+aoiFile+" -s "
                     +startDate+" -e "+endDate
                     +" -d --producttype SLC -q \"orbitdirection="
                     +orbit+",relativeorbitnumber="+orbitNumber+
                     "\" --url \"https://scihub.copernicus.eu/dhus\"")
            # Execute sentinelsat download and wait for it to finish (.call method)
            subprocess.call("sentinelsat "+parms,shell=True)
        elif downloadMethod == 'aria2c':
            #%% download variant with aria2c
            parms = ("--http-auth-challenge=true --http-user=\'"
                     +user+"\' --http-passwd=\'"+password
                     +"\' --continue=true --auto-file-renaming=false"
                     +" \"https://api.daac.asf.alaska.edu/services/search/param?platform=Sentinel-1&intersectsWith="
                     +aoiPolygon+")&start="+startDate[0:4]+"-"+startDate[4:6]
                     +"-"+startDate[6:8]+"T00:00:00" + "&end="+endDate[0:4]+"-"
                     +endDate[4:6]+"-"+endDate[6:8]+
                     "T00:00:00&processingLevel=SLC&relativeOrbit="
                     +orbitNumber+"&output=metalink\"")
            subprocess.call("aria2c "+parms,shell=True)

if __name__ == "__main__":
    # load parameters file:
    parmFile = sys.argv[1]
    with open(parmFile,'r') as inp:
        parameters = eval(inp.read())
    # call function:
    batchDownload(parameters)
