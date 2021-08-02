#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 16 18:14:25 2021

@author: rc
"""

import sys
import os
from gecoris import insarUtils, ioUtils, plotUtils, insar

def CR_insar(input_parms, decomp = True):

    #%% inputs:
    station_log = input_parms['stationLog']
    stacks_log = input_parms['stackLog']
    process_dir = input_parms['outDir']
    # initialise:
    stations, refl_flag = ioUtils.load_stations(station_log, process_dir)
    
    #%% do insar analysis for all stacks:
    insar_parms = ioUtils.prepare_insar_parms(station_log, stacks_log, process_dir)
    for parm in insar_parms:
        # chcek if insar proc. done already:
        hdf_file = parm['outDir'] + 'insar_' + parm['stackId'] + '.hdf5'
        if not os.path.isfile(hdf_file):
            insar.insar(parm)
        else:
            print('InSAR analysis for stack'+parm['stackId']+'already performed.')
    
    #%% plot individual reflectors:
    for i in range(len(stations)):
        for parm in insar_parms:
            hdf_file = parm['outDir'] + 'insar_' + parm['stackId'] + '.hdf5'
            stations[i].add_insar_ts(hdf_file, parm['stackId'])
        # plot TS:
        stations[i].plot_insar_ts(out_dir = parm['outDir'])
        # TODO: add decomposition
        if decomp:
            stations[i].decomp_insar_ts()
            stations[i].plot_decomp_ts(out_dir = parm['outDir'])
    
#%%

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
    if len(sys.argv) == 2:
        parms = parse_parms(sys.argv[1])
        CR_insar(parms)
    elif len(sys.argv) > 2:
        parms = parse_parms(sys.argv[1])
        CR_insar(parms, decomp=True)
    else:
        print('Not enough input arguments!')