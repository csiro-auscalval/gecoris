#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""ioUtils

Module containing input/output utility functions

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

import glob
import json
import numpy as np
import os
import csv
import copy
# gecoris packages:
from gecoris import classes, geoUtils


def toJSON(cr,outDir):
    """
    Serialize Reflector instance to JSON in specified output directory.
    """
    if isinstance(cr, classes.Stack):
        outFile = outDir + os.sep + 'stack_' + cr.id + '.json'
    else:
        outFile = outDir + os.sep + cr.id + '.json'
    with open(outFile, "w") as write_file:
        json.dump(cr, write_file, cls = classes.crEncoder,indent=2)
  
def fromJSON(inFile):
    """
    Load Reflector instance from specified JSON (full path)
    """     
    with open(inFile, "r") as in_file:
        cr = json.load(in_file, cls = classes.crDecoder)
    return cr

def fromCSV(inFile):
    """Initialise reflector instances from CSV log file
    
    CSV station log file assumes same geo coordinates for both ascending and 
    descending reflectors of the station, given in ETRS-89, ETRF2000, 2010.0
    
    input: CSV station log
    output: list of Reflector instances
    """
    CRs = []
    with open(inFile, "r") as f:
        q = csv.reader(f, delimiter=',')
        for p in q: 
            # skip header
            if p[0] == 'ID':
                continue
            # instantize Coordinates:
            asc = classes.Coordinates(float(p[5])/180*np.pi,
                                      float(p[6])/180*np.pi,
                                      float(p[7]))
            dsc = copy.copy(asc)
            # create instance of reflector:
            if p[1] =='CREF':
                CRs.append(
                    classes.CornerReflector(p[0],asc,dsc,
                                            p[2][:8],
                                            p[3][:8],
                                            p[4][:8],
                                            p[10],
                                            float(p[11]),
                                            float(p[12]),
                                            float(p[13])))
            elif p[1] == 'CAT':
                CRs.append(
                    classes.Transponder(p[0],asc,dsc,
                                        p[2][:8],
                                        p[3][:8],
                                        p[4][:8],
                                        None))
            else:
                print('Not supported type of reflector: '+p[1]+ ', skipping.')
    return CRs

def stacks_fromCSV(inFile):
    """Initialise stack instances from CSV log file
    
    input: CSV stacks log
    output: list of Stack instances
    """
    stacks = []
    with open(inFile, "r") as f:
        q = csv.reader(f, delimiter=',')
        for p in q: 
            # skip header
            if p[1].upper() == 'ID':
                continue
            else:
                # create instance of Stack:
                stacks.append(
                    classes.Stack(p[1],p[2],p[0],p[3],p[4]))
    return stacks

def csv_xyz2plh(inFile,outFile):
    """
    Convert ECEF to ellipsoidal coordinates in CSV file
    """
    with open(inFile, "r") as f:
        with open(outFile,"w") as out_f:
            q = csv.reader(f, delimiter=',')
            out = csv.writer(out_f, delimiter = ',')
            for p in q:
                xyz = np.array([float(p[1]),float(p[2]),float(p[3])])
                plh = geoUtils.xyz2plh(xyz)
                outStr = ([p[0],
                          f"{plh[0]*180/np.pi:.10f}",
                          f"{plh[1]*180/np.pi:.10f}",
                          f"{plh[2]:.10f}"])
                out.writerow(outStr)

def csv_plh2xyz(inFile,outFile):
    """
    Convert ECEF to ellipsoidal coordinates in CSV file
    """
    with open(inFile, "r") as f:
        with open(outFile,"w") as out_f:
            q = csv.reader(f, delimiter=',')
            out = csv.writer(out_f, delimiter = ',')
            for p in q:
                plh = np.array([float(p[6])/180*np.pi,
                                float(p[7])/180*np.pi,
                                float(p[8])])
                xyz = geoUtils.plh2xyz(plh)
                outStr = ([p[0],
                          f"{xyz[0]:.4f}",
                          f"{xyz[1]:.4f}",
                          f"{xyz[2]:.4f}"])
                out.writerow(outStr)


def statsToJSON(cr,outDir):
    """
    Export reflector instance statistics to JSON in specified output directory
    """    
    outFile = outDir + os.sep + cr.id + '_stats.json'
    stats = {}
    for stack in cr.stacks:
        stats[stack['id']] = cr.getStackStats(stack['id'])
    
    with open(outFile, "w") as write_file:
        json.dump(stats, write_file,indent=2)
        

def loadJSON(inJSON):
    """
    Load JSON as dictionary
    """    
    with open(inJSON, "r") as read_file:
          jsonDict = json.load(read_file)
    return jsonDict

def load_station_log(inJSON):
    """Initialise reflector instances from JSON station-log-file
    
    input: full path to JSON station-log-file
    output: list of Reflector instances
    """    
    with open(inJSON, "r") as read_file:
          jsonDict = json.load(read_file)
          stations = []
          for station in jsonDict['stations']:
              if station['ascending']['orientation']:
                  coords = station['ascending']['coordinates']
                  ascending = classes.Coordinates(
                      coords['latitude']/180*np.pi, 
                      coords['longitude']/180*np.pi, 
                      coords['elevation'],
                      coords['CRS'],
                      coords['FRAME'],
                      coords['EPOCH'],
                      coords['EPSG'])
              else:
                  ascending = None
              if station['descending']['orientation']:
                  coords = station['descending']['coordinates']
                  descending = classes.Coordinates(
                      coords['latitude']/180*np.pi, 
                      coords['longitude']/180*np.pi, 
                      coords['elevation'],
                      coords['CRS'],
                      coords['FRAME'],
                      coords['EPOCH'],
                      coords['EPSG'])
              else:
                  descending = None
              if station['type'] == 'CREF':
                  stations.append(
                      classes.CornerReflector(
                          station['id'],
                          ascending,
                          descending,
                          station['installDate'][:8],
                          station['startDate'][:8],
                          station['endDate'][:8],
                          station['geometry']['shape'],
                          station['geometry']['leglength'],
                          station['geometry']['zenithDip'],
                          station['geometry']['azimuthDip']))
              elif station['type'] == 'CAT':
                  stations.append(
                      classes.Transponder(
                          station['id'],
                          ascending,
                          descending,
                          station['installDate'][:8],
                          station['startDate'][:8],
                          station['endDate'][:8],
                          station['RCS0']))
          
    return stations

def load_stacks(inJSON):
    """Initialise stack instances from JSON station-log-file
    
    input: full path to JSON station-log-file
    output: list of Stack instances
    """      
    with open(inJSON, "r") as read_file:
          jsonDict = json.load(read_file)
          stacks = []
          for stack in jsonDict['stacks']:
              stacks.append(
                  classes.Stack(
                          stack['id'],
                          stack['dataPath'],
                          stack['sensor'],
                          stack['subswath'],
                          stack['type']))    
    return stacks


def listHDF(data):
    print('Datasets:')
    for k,v in data.items():
        print(v)
    print('Attributes:')
    for k,v in data.attrs.items():
        print(k,v)
        

def load_stations(station_log, process_dir):
    """Load already analysed stations
    
    input: full path to JSON station-log-file and process dir
    output: list of Station instances and flag = True if analysis performed
    """ 
    if station_log.lower().endswith('.json'):
        stations = load_station_log(station_log)
    elif station_log.endswith('.csv'):
        stations = fromCSV(station_log)
    else:
        print('Unknown station log file format. Use .json or .csv.')
        raise
    # load analysed .json logs:
    refl_flag = True
    try:
        logs = glob.glob(process_dir + os.sep +"*.json")
        for i in range(len(stations)):
            inJSON = [q for q in logs 
                          if stations[i].id+'.json' in q.split(os.sep)[-1]]
            stations[i] = fromJSON(inJSON[0])
    except:
        print('No refl. analysis performed for reflectors in specified dir.')
        refl_flag = False
    return stations, refl_flag


def prepare_insar_parms(station_log, stacks_log, outDir, aoi = 10):

    if station_log.lower().endswith('.json'):
        stations = load_station_log(station_log)
    elif station_log.endswith('.csv'):
        stations = fromCSV(station_log)
    else:
        print('Unknown station log file format. Use .json or .csv.')
        raise
    if stacks_log.lower().endswith('.json'):
        stacks = load_stacks(stacks_log)
    elif stacks_log.lower().endswith('.csv'):
        stacks = stacks_fromCSV(stacks_log)
    else:
        print('Unknown stacks log file format. Use .json or .csv.')
    
    parms = dict()
    padd = aoi/2*1e3/30/3600
    # determine aoi (+- 0.1 deg ~ +- 10 km):
    plh = np.array([s.get_plh()[0] for s in stations])
    parms['min_lat'] = np.min(plh[:,0])*180/np.pi - padd
    parms['max_lat'] = np.max(plh[:,0])*180/np.pi + padd
    parms['min_lon'] = np.min(plh[:,1])*180/np.pi - padd
    parms['max_lon'] = np.max(plh[:,1])*180/np.pi + padd
    parms['aver_h'] = np.mean(plh[:,2])
    # determine start:
    parms['startDate'] = min([s.startDate for s in stations])
    
    parms_all = [dict() for i in range(len(stacks))]
    # determine stacks:
    for stack, par in zip(stacks, parms_all):
        par['stackDir'] = stack.stackDir
        par['stackId'] = stack.id
        par['subswath'] = stack.subswath
        par['outDir'] = outDir + os.sep + '/insar/'
        par.update(parms)
    return parms_all


#def progress_bar(i, n):
    