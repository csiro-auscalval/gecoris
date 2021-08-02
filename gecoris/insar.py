#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Main function to perform InSAR analysis using GECORIS

R. Czikhardt 
4.3.2021
"""

import sys
from gecoris import insarUtils, plotUtils, ioUtils

def insar(input_parms):

    #%% DEFAULT PARAMETERS:
    
    parms = {
        # parameters:
        'D_A_thr' : 0.25,       # Norm. Amplitude Dispersion (NAD) threshold on PS sel.
        'model' : 'seasonal',   # func. model to solve ambiguities, 'linear' / 'seasonal'
        'reference' : 'auto',   # 'auto' / 'free' / 'contraint'
        'APS_flag' : 1,         # 0/1 = estimate and remove Atmospheric Phase Screen (APS)
        'plot_flag': 1,         # plots, 0 = none, 1 = sparse, 2 = detailed
    }
    
    # update default parms by input:
    parms.update(input_parms)
        
    #%% prepare insar data in HDF:
    data = insarUtils.prepare_insar_HDF(parms)
        
    #%% or open existing:
    #HDFfile = '/data/gecoris_test/insar/insar_DSC51.hdf5'
    #data = insarUtils.openHDF(HDFfile)
    
    #%% plot PSC
    if parms['plot_flag']:
        plotUtils.plot_psc(data['psc'], parms['outDir'] + 'psc_D_A.png')
    
    #%% create network:
    data['network/arcs'] = insarUtils.createNetwork(data['psc'], data['stack'],
                                                    n_type = 'delaunay',
                                                    plotFlag = parms['plot_flag'])
    
    #%% temporal ambig:
    insarUtils.temporalAmbiguity(data['stack'],
                                 data['psc'],
                                 data['network'], model='seasonal')
    if parms['plot_flag']:
        plotUtils.plot_network_quality(data['network'], outDir = parms['outDir'])
    
    #%% outlier detection:
    insarUtils.remove_outliers(data)
    
    #%% ref. point:
    insarUtils.setRefPoint(data['psc2'], data['network2'], 
                           method='auto') # TODO: change to parms
    # check network conditioning:
    insarUtils.network_cond(data) 
    
    #%% spatial ambiguity:    
    insarUtils.spatialAmbiguity(data['network2'])
        
    #%% network solution
    insarUtils.solve_network(data['network2'])
    
    #%% plot it:
    if parms['plot_flag']:
        plotUtils.plot_network_parms(data['network2'], data['psc2'], 
                                     parms['outDir']+'/parms_1st.png')
    
    #%% APS estimation:
    # TODO: aps_flag
    insarUtils.APS(data['psc2'], data['network2'], data['stack'], data['psc'], 
                   plotFlag = parms['plot_flag'], 
                   apsDir = parms['outDir'] + '/aps/')
    
    #%% 2nd iteration after APS:
    insarUtils.temporalAmbiguity(data['stack'], 
                                 data['psc'], 
                                 data['network2'], model='seasonal')
    if parms['plot_flag']:
        plotUtils.plot_network_quality(data['network2'])
    #%% 2nd outlier removal:
    del data['network']
    data['network'] = data['network2']    
    insarUtils.remove_outliers(data)
    #%% ref. point (former):
    data['network2'].attrs['refIdx'] = insarUtils.ref_coords2idx(
        data['psc2'], 
        data['network2'].attrs['refAz'], 
        data['network2'].attrs['refR'])
    # check network conditioning:
    insarUtils.network_cond(data) 
        
    #%% spatial ambiguity:    
    insarUtils.spatialAmbiguity(data['network2'])
    
    #%% network solution
    #insarUtils.solve_network(data['network2'])
    
    #%% PRECISE NETWORK SOL:
    stack = ioUtils.fromJSON(parms['outDir']+'stack_'+parms['stackId']+'.json')
    insarUtils.getPreciseH2PH(data['psc2'], stack)
    insarUtils.solve_network_precise(data['network2'], data['psc2'], 
                                     model='seasonal')
    
    #%% plot it:
    if parms['plot_flag']:
        plotUtils.plot_network_parms(data['network2'], data['psc2'], 
                                     parms['outDir']+'/parms_2nd_fixh.png')
    #%% remove RPN and re-solve
    insarUtils.remove_RPN(data['network2'])
    insarUtils.solve_network(data['network2'])
    if parms['plot_flag']:
        plotUtils.plot_network_parms(data['network2'], data['psc2'], 
                                     parms['outDir']+'/parms_2nd_fixh.png')
    
    #%% plot it:
    if parms['plot_flag']:
        plotUtils.plot_network_parms(data['network2'], data['psc2'], 
                                     parms['outDir']+'/parms_2nd.png')
    
    #%% export HDF to CSV:
    outCSV = parms['outDir']+'/insar_'+parms['stackId']+'.csv'
    insarUtils.HDF2csv(data, outCSV)
    # -optionally- convert to shapefile:
    #csv2shp(outCSV, outCSV.replace('.csv','.shp')) 
    
    #%% close dataset:
    data.close()


def parse_parms(parms_file):
    if isinstance(parms_file, str):
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
    elif isinstance(parms_file, dict):
        return parms_file
    else:
        print('Wrong input format.')
        raise


if __name__ == "__main__":
    # load parameters:
    if len(sys.argv) > 1:
        parms = parse_parms(sys.argv[1])
        insar(parms)
    else:
        print('Not enough input arguments!')