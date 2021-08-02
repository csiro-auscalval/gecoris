#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""batchPrepare

Routine to prepare S1 SLC bursts and extract AOI for SNAP processing, 
optimised for cross-swath

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

import os
import sys
import subprocess
import time
import zipfile
import re
import xml.etree.ElementTree as ET
from shapely import wkt,geometry


def batchPrepare(parameters):
    
    # PARAMETERS:
    # read input parameters:
    snapPath = parameters['SNAP_parameters']['snapPath']
    repoPath = parameters['SNAP_parameters']['repoPath']
    ram = parameters['SNAP_parameters']['ram']
    cores = parameters['SNAP_parameters']['cores']
    workDir = parameters['Stack_parameters']['workDir']
    dataDirs = parameters['Stack_parameters']['dataDirs']
    swaths = parameters['Stack_parameters']['swaths']
    min_lon = parameters['Stack_parameters']['min_lon']
    max_lon = parameters['Stack_parameters']['max_lon']
    min_lat = parameters['Stack_parameters']['min_lat']
    max_lat = parameters['Stack_parameters']['max_lat']
    startDate = parameters['Stack_parameters']['startDate']
    endDate = parameters['Stack_parameters']['endDate'] 
    ram = str(ram)
    cores = str(cores*10)
    
    aoi = ('POLYGON (('+min_lon+' '+min_lat+', '+max_lon+' '+min_lat+', '
           +max_lon+' '+max_lat+', '+min_lon+' '+max_lat+', '+min_lon+' '
           +min_lat+'))')
    
    # set endDate to current date by default:
    if not endDate:
        endDate = time.strftime('%Y%m%d')
    
    delimiter = os.sep
    graph_assembly = delimiter + 'graphs/PrepareGraph_assembly.xml'
    graph_single = delimiter + 'graphs/PrepareGraph_single.xml'
    
    # prepare output stack directory and find out if only update requested
    out_dir = workDir + 'slaves/'
    
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
        update = False
        print('Specified output dir does not exist, creating one.')
    else:
        formerSlaves = [slv[0:8] for slv in os.listdir(out_dir) 
                        if slv.endswith('.dim')]
        update = True
        print('Stack already exists, updating...') 
    
    # load dates and filenames of input SLCs
    img_files = []
    img_date = []
    for in_dir in dataDirs:
        dir_contents = os.listdir(in_dir)
        dir_contents.sort()
        for fname in dir_contents:
            if fname.endswith('.zip'): # test if it is S1 .zip FILE
                temp = fname[17:25] # parse date
    			# test if within requested date period
                if startDate <= temp <= endDate:
                    # test if more than one scene for particular date
                    if img_date[-1:] == [temp]: # -1: collon here in order not to receive error on first iteration    
                        img_files[-1] = ', '.join([img_files[-1], in_dir 
                                                   + delimiter + fname])
                    else:
                        img_files.append(in_dir + delimiter + fname) # list of filenames
                        img_date.append(temp)  
                    #print(temp)
    
    #%% if update, select only new SLCs:
    if update:
        img_date_new = [img for img in img_date if img not in formerSlaves]
        if len(img_date_new)==0:
            print('No new images to process, aborting.')
            sys.exit()
        img_files_new = []
        for img in img_date_new:         
            img_files_new.append(img_files[img_date.index(img)])            
        img_date = img_date_new
        img_files = img_files_new
        #if len(img_date)==1:
        #    img_files = [img_files]
    
    #img_files.sort()
    print('New images to process:')
    print(img_files)
    if not isinstance(img_files,list):
        img_files = [img_files]
    
    #%% Iterate through dates and run prepareGraph in SNAP gpt
    for img in img_date:
        # check if more acquisition dates:
        if len(img_date) > 1:
            inp = img_files[img_date.index(img)].split(',') # index input
        else: # only single SLC update
            # check if single/more SLC per date:
            if isinstance(img_files[0],list):
                inp = img_files[0]
            else:
                inp = img_files
            if ',' in inp[0]:
                inp = inp[0].split(',')
        proc = False
        for sfile in inp:
            sfile = sfile.strip() #remove trailing spaces
            # test if aoi within frame:
            try: # use try if zip.file corrupted
                if inSLCframe(sfile,aoi):
                    proc = True
                    # do the processing
                    for swath in swaths:
                        ofile = '\"' + out_dir + img + '_' + swath + '.dim"'
                        cmd = ('gpt -c '+ram+'m -q '+cores+' -x -e ' + '\"' + repoPath + graph_single + '\" ' + '-PsourceFile1=' + 
                               '\"' + sfile + '\"' + ' -Pswath=\"' + swath + '\"' + ' -Paoi=\"' + aoi + '\"' + ' -PoutputFile=' + ofile)
                        # switch between Linux/Windows:
                        if os.name == 'nt':
                            os.system(cmd)
                        else:
                            subprocess.call([snapPath+cmd],shell=True)
            except: # zip file corrupted
                print('File '+ sfile +' corrupted, skipping.')
                proc = True
        if not proc:
            print('Assembling SLCs...')
            # make assembly:
            sfiles = ",".join(['\"' + e + '\"' for e in inp]) # adds "" and then joins with comma
            for swath in swaths:
                ofile = '\"' + out_dir + img + '_' + swath + '.dim\"'
                cmd = ('gpt -c '+ram+'m -q '+cores+' -x -e ' + '\"' + repoPath + graph_assembly
                   + '\" ' + '-PsourceFiles=' + sfiles + ' -Pswath=\"' + swath + '\"' + ' -Paoi=\"' +
                aoi + '\"' + ' -PoutputFile=' + ofile)
                #print(cmd)
                # switch between Linux/Windows:
                if os.name == 'nt':
                    os.system(cmd)
                else:
                    subprocess.call([snapPath+cmd],shell=True)
                
    print('All SLC images cropped and prepared.')

def inSLCframe(slcFile,wktAOI):
    # open manifest.safe from SLC .zip file:
    try:
        zf = zipfile.ZipFile(slcFile, 'r')
    except:
        raise Exception("Zip file corrupt.")
    manifest = [n for n in zf.namelist() if re.search("manifest.safe",n)]
    f = zf.open(manifest[0])
    
    # parse gml frame from manifest.safe:
    root = ET.fromstring(f.read())
    frame = root[1].find("./metadataObject[@ID='measurementFrameSet']")
    gmlFrame = frame[0][0][0][0][0][0].text 
    
    # covert gml frame to shapely polygon:
    gmlFrame = [[float(i) for i in j.split(',')] for j in gmlFrame.split(' ')]
    frame = geometry.Polygon([i[::-1] for i in gmlFrame]) # change point order for shapely
    
    aoi = wkt.loads(wktAOI) # convert input wkt to shapely polygon
    
    return frame.contains(aoi)
    
if __name__ == "__main__":   
    # load parameters file:
    parmFile = sys.argv[1]
    with open(parmFile,'r') as inp:
        parms = eval(inp.read())
    # call function:
    batchPrepare(parms)