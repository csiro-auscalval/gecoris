#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Main function to perform radar reflector analysis using GRTx

R. Czikhardt 
29.7.2020
"""

print('Geodetic Corner Reflector (In)SAR Toolbox (GECORIS) v.1.0')
print('Copyright (c) 2021 Richard Czikhardt, czikhardt.richard@gmail.com')
print('Dept. of Theoretical Geodesy, Slovak University of Technology')
print('-----------------------------------------------------------------')
print('License: GNU GPL v3+')
print('-----------------------------------------------------------------')

import sys
import os
import glob
from gecoris import ioUtils, plotUtils, atmoUtils


def main(parms):    
    # unpack:
    stationLog = parms['stationLog']
    stacksLog = parms['stackLog']
    outDir = parms['outDir']
    posFlag = parms['precisePosFlag']
    ovsFactor = parms['ovsFactor']
    atmoFlag = parms['atmoFlag']
    plotFlag = parms['plotFlag']
    #
    print('Initializing reflectors...')
    if stationLog.lower().endswith('.json'):
        stations = ioUtils.load_station_log(stationLog)
    elif stationLog.endswith('.csv'):
        stations = ioUtils.fromCSV(stationLog)
    else:
        print('Unknown station log file format. Use .json or .csv.')
        raise
    # prepare atmo. data dir:
    if atmoFlag:
        atmoDir = outDir + os.sep + 'atmo/'
        if not os.path.exists(atmoDir):
            os.makedirs(atmoDir)
    #
    print('Loading SAR data stacks...')
    if stacksLog.lower().endswith('.json'):
        stacks = ioUtils.load_stacks(stacksLog)
    elif stacksLog.lower().endswith('.csv'):
        stacks = ioUtils.stacks_fromCSV(stacksLog)
    else:
        print('Unknown stacks log file format. Use .json or .csv.')
    # load data:
    for stack in stacks: 
        stack.readData()
        ioUtils.toJSON(stack, parms['outDir']) # save to JSON
        if atmoFlag:
            print('Preparing atmo. models for stack '+stack.id)
            atmoUtils.prepareAtmo(stack, atmoDir)        
    # check if output directory exists and load already processed:
    if not os.path.exists(outDir):
        os.makedirs(outDir)
        logs = []
    else:
        logs = glob.glob(outDir + os.sep +"*.json")
    #
    print(str(len(stations))+ ' reflectors on ' +str(len(stacks))
          + ' stacks to process.')
    # iterate over list indices to modify objects:
    for i in range(len(stations)):
        # check if analysis already performed:
        inJSON = [q for q in logs 
                  if stations[i].id+'.json' in q.split(os.sep)[-1]]
        if inJSON:
            print('Station '+stations[i].id+ ' already analyzed, updating.')
            stations[i] = ioUtils.fromJSON(inJSON[0])
            for stack in stacks:
                if atmoFlag:
                    stations[i].updateStack(stack, ovsFactor=ovsFactor, 
                                            posFlag=posFlag, atmoDir=atmoDir)
                else:
                    stations[i].updateStack(stack, ovsFactor=ovsFactor,
                                            posFlag=posFlag)
        else:
            for stack in stacks:
                if atmoFlag:
                    stations[i].addStack(stack, ovsFactor=ovsFactor, 
                                         posFlag=posFlag, plotFlag=plotFlag,
                                         outDir=outDir+stations[i].id,
                                         atmoDir=atmoDir)
                else:
                    stations[i].addStack(stack, ovsFactor=ovsFactor, 
                                         posFlag=posFlag, plotFlag=plotFlag,
                                         outDir=outDir+stations[i].id)
        stations[i].print_all_timeseries(outDir)
        print('Removing outliers.')
        stations[i].detectOutliers()
        print('Performing RCS and SCR analysis.')
        stations[i].RCSanalysis()
        if posFlag and plotFlag > 0:
            print('Plotting ALE.')
            for stack in stacks:
                if atmoFlag:
                    stations[i].plotALE(stack.id,outDir,atmoDir)
                else:
                    stations[i].plotALE(stack.id,outDir)
            stations[i].plotALE_TS(outDir)
        if plotFlag > 0:
            print('Plotting RCS time series.')
            stations[i].plotRCS(outDir)
        #
        print('Exporting to JSON dumps.')
        stations[i].toJSON(outDir)
        stations[i].statsToJSON(outDir)
    #
    if len(stations) > 2 and plotFlag > 0:
        print('Generating network plots.')
        plotUtils.plotNetworkALE(stations, outDir+os.sep+'network_ALE.png')
        plotUtils.plotNetworkRCS(stations, outDir+os.sep+'network_RCS.png')
        plotUtils.plotNetworkSCR_hz(stations, outDir+os.sep+'network_SCR.png')
    #
    print('Done. Thank you for using GECORIS. Do not forget to reference.')
    print('In case of any inquries, please contact czikhardt.richard@gmail.com')


def parse_parms(parms_file):
    try:
        with open(parms_file,'r') as inp:
            try:
                parms = eval(inp.read())
            except:
                print("Something wrong with parameters file.")
                raise
        return parms
    except:
        print("Specified parameters file not found.")
        raise


if __name__ == "__main__":
    # load parameters:
    if len(sys.argv) > 1:
        parms = parse_parms(sys.argv[1])
        main(parms)
    else:
        print('Not enough input arguments!')
