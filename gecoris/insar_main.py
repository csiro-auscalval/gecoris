#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""insar_main.py

Main template script to perform InSAR time series analysis.

Copyright (C) 2021 by R.Czikhardt

Email: czikhardt.richard@gmail.com
Last edit: 6.7.2021
"""

from gecoris import insarUtils, plotUtils

#%% input parameters:

parms = {
    'stackId' : 's1_DSC124',
    'stackDir' : '/data/GUDS/CR_Prievidza/DSC124/',
    'subswath' : 'IW1',
    'startDate' : '',
    #% define aoi:
    'min_lon' : 18.60,
    'max_lon' : 18.74,
    'min_lat' : 48.70,
    'max_lat' : 48.79,
    'aver_h' : 440,
    # parameters:
    'D_A_thr' : 0.4,       # Norm. Amplitude Dispersion (NAD) threshold on PS sel.
    'oversample_factor' : 16, # oversampling to detect sub-pix position of PS, keep = 1 to skip
    'network' : 'redundant', # network formation strategy: 'redundant'/ 'delaunay'
    'model' : 'seasonal',   # func. model to solve ambiguities, 'linear' / 'seasonal'
    'bounds' : {    # soft estimation bounds 
        'dH' : 20, # [m]
        'vel': 5,  # [mm/year]
        'seasonal': 2 # [mm]
        },
    'reference' : 'auto',   # 'auto' / 'free' / 'contraint'
    'ref_longitude' : None, # ref. point coordinates
    'ref_latitude' : None, # only use in case of 'constraint' solution
    'ref_height' : None, # ellipsoidal
    'APS_flag' : 1,         # 0/1 = estimate and remove Atmospheric Phase Screen (APS)
    'atmoWindow' : 200 ,     # atmo-filter window length in [days]
    'densify_flag': 1,       # densify 1st order network by 2nd order PSc
    'densify_stdThr' : 6,   # [mm], threshold for densification
    'plot_flag': 1,         # plots, 0 = none, 1 = sparse, 2 = detailed
    # outputs:
    'outDir' : '/data/insarsk/insarsk_workshop/HN/'
}
    
#%% prepare insar data in HDF:
data = insarUtils.prepare_insar_HDF(parms)

# or open existing:
#data = insarUtils.openHDF('<path_to_HDF>')

#%% order to 1st/2nd order PSc:
insarUtils.order_psc(data)

#%% plot PSC
if parms['plot_flag']:
    plotUtils.plot_psc(data['psc'], 
                       parms['outDir'] + '/psc1_NAD.png')
    if 'psc_B' in data:
        plotUtils.plot_psc(data['psc_B'], 
                           parms['outDir'] + '/psc2_NAD.png')

#%% create network:
data['network/arcs'] = insarUtils.createNetwork(data['psc'], 
                                                data['stack'],
                                                n_type = parms['network'],
                                                plotFlag = parms['plot_flag'],
                                                outDir = parms['outDir'])

#%% temporal ambiguities:
insarUtils.temporalAmbiguity(data['stack'],
                             data['psc'],
                             data['network'], 
                             model = parms['model'],
                             bounds = parms['bounds'])

if parms['plot_flag']:
    plotUtils.plot_network_quality(data['network'], 
                                   outDir = parms['outDir'])

#%% outlier detection:
insarUtils.remove_outliers(data)

#%% ref. point (datum):
insarUtils.setRefPoint(data['psc2'], 
                       data['network2'], 
                       method = parms['reference'],
                       refLon = parms['ref_longitude'],
                       refLat = parms['ref_latitude'])
# check reference datum conditioning
insarUtils.network_cond(data)

#%% spatial ambiguity:    
insarUtils.spatialAmbiguity(data['network2'])

#%% network solution
insarUtils.solve_network(data['network2'])

if parms['plot_flag']:
    plotUtils.plot_network_parms(data['network2'], 
                                 data['psc2'], 
                                 parms['outDir']+'/parms_1st.png')

#%% APS estimation:
if parms['APS_flag']:
    insarUtils.APS(data['psc2'], 
                   data['network2'], 
                   data['stack'], 
                   data['psc'],
                   atmoWindow = parms['atmoWindow'],
                   plotFlag = parms['plot_flag'],
                   apsDir = parms['outDir'] + '/aps/')
    #% 2nd iteration after APS:
    insarUtils.temporalAmbiguity(data['stack'], 
                                 data['psc'], 
                                 data['network2'], 
                                 model = parms['model'],
                                 bounds = parms['bounds'])
    if parms['plot_flag']:
        plotUtils.plot_network_quality(data['network2'],
                                       outDir = parms['outDir'])
    #% 2nd outlier removal:
    del data['network']
    data['network'] = data['network2']    
    insarUtils.remove_outliers(data)
    #% ref. point (former):
    if parms['reference'] != 'free':
        data['network2'].attrs['refIdx'] = insarUtils.ref_coords2idx(
            data['psc2'], 
            data['network2'].attrs['refAz'], 
            data['network2'].attrs['refR'])
    #% spatial ambiguity:    
    insarUtils.spatialAmbiguity(data['network2'])

#%% precise network solution
insarUtils.getPreciseH2PH(data['psc2'], 
                          parms['outDir'],
                          parms['stackId'])
insarUtils.solve_network_precise(data['network2'], 
                                 data['psc2'], 
                                 model = parms['model'])

#%% remove RPN and re-solve
insarUtils.remove_RPN(data['network2'],
                           plotFlag = parms['plot_flag'],
                           outDir = parms['outDir'])
insarUtils.solve_network_precise(data['network2'], 
                                 data['psc2'], 
                                 model = parms['model'])
if parms['plot_flag']:
    plotUtils.plot_network_parms(data['network2'], 
                                 data['psc2'], 
                                 parms['outDir']+'/parms_2nd.png')

#%% densify network:
if parms['densify_flag']:
    insarUtils.getPreciseH2PH(data['psc_B'], 
                              parms['outDir'],
                              parms['stackId'])
    insarUtils.densify_network(data, 
                               k = 3, 
                               mode_thr = 1, 
                               std_thr = parms['densify_stdThr'])
    
#%% make positioning corrections:
insarUtils.fix_geocoding(data, parms)
    
#%% export to CSV
outCSV = parms['outDir']+'/insar_'+parms['stackId']+'.csv'
insarUtils.HDF2csv(data, outCSV)