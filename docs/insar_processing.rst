.. _insar_processing:

****************
InSAR processing
****************

Main InSAR time series analysis workflow based on the TU Delft geodetic estimation theory.
This module assumes you have already coregistered SLC stacks of your AOI using SNAP, see :ref:`snap_processing`.
All consequent processing steps are contained in the script file :download:`insar_main.py <../gecoris/insar_main.py>`

.. warning:: 

   Currently, this is a limited access feature, see :ref:`intro`.

We'll use two gecoris libraries for processing and plotting:

.. code:: python

    from gecoris import insarUtils, plotUtils

Processing parameters
=====================

Define the processing parameters in ``parms`` dictionary:

.. code:: python

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


.. note::

   Detailed description of the parameters can be found in: :ref:`parms_table`.


Extract data from SNAP processing
=================================

Creating new dataset
--------------------

Extract SNAP-coregistered stack located at ``parms['stackDir']``:

.. code:: python

    data = insarUtils.prepare_insar_HDF(parms)

This will create two binary `HDF5 files <https://www.hdfgroup.org/solutions/hdf5/>`_  and JSON file in your ``parms['outDir']``: 

* ``stack_<stackId>.hdf5`` - containing SLC data stack,

* ``insar_<stackId>.hdf5`` - containing Persistent Scatter candidates (PSc) selected
  using ``parms['D_A']`` threshold,

* ``stack_<stackId>.json`` - containing stack metadata.

If the datastack has been already extracted in previous session, only PSc selection is run.
You can also open an existing HDF datastack:

.. code:: python

    data = insarUtils.openHDF('insar_<stackId>.hdf5')

.. note::

   You can list contents of the currently open HDF datastack with ``data.keys()`` method. 
   
   Full HDF datastack specificiation is in: :ref:`HDF_datastack`.

.. warning::

   Always close the already open datastack with ``data.close()`` method before opening new one.

Updating existing dataset with new SLCs
---------------------------------------

If you have existing HDF datastack in ``parms['stackDir']`` and only want to update it with new SLC images:

.. code:: python

    data = insarUtils.prepare_insar_HDF(parms, update=True)


First-order network
===================

Order PSc
---------

Split PSc into 1st and 2nd order, based on NAD threshold of 0.25. Optional second argument can be used to modify this default value.

.. code:: python

    insarUtils.order_psc(data)

Plot the PSc in geographic coordinates:

.. code:: python

    if parms['plot_flag']:
    	plotUtils.plot_psc(data['psc'], 
                           parms['outDir'] + '/psc1_NAD.png')
    if 'psc_B' in data:
        plotUtils.plot_psc(data['psc_B'], 
                           parms['outDir'] + '/psc2_NAD.png')

.. _create_network:

Create network
--------------

.. code:: python

   data['network/arcs'] = insarUtils.createNetwork(data['psc'], 
                                                   data['stack'],
                                                   n_type = parms['network'],
                                                   plotFlag = parms['plot_flag'],
                                                   outDir = parms['outDir'])

Solve temporal ambiguities
--------------------------

Ambiguities are estimated over all network arcs:

.. code:: python

   insarUtils.temporalAmbiguity(data['stack'],
                                data['psc'],
                                data['network'], 
                                model = parms['model'],
                                bounds = parms['bounds'])

.. note::
   
   The ``bounds`` parameter defines the dictionary of soft estimation bound for unknown parameters
   (residual height, displacement model coefficient), and is site and deformation phenomena-specific.

Overall quality statistics of the estimation can be plotted using:

.. code:: python

   if parms['plot_flag']:
       plotUtils.plot_network_quality(data['network'], 
                                      outDir = parms['outDir'])
                                      
Unreliable (outlying) arcs are then removed automatically, or using specific threshold on the 
standard deviation of residuals (RMSE) as an optional second argument:

.. code:: python

   insarUtils.remove_outliers(data)

.. note::
   
   Outleir removal step also automatically removes points that might become isolated by unreliable arcs removal
   in order to attain the network consistency.

Integrate network spatially
---------------------------

First, define the reference datum:

.. code:: python

   insarUtils.setRefPoint(data['psc2'], 
                          data['network2'], 
                          method = parms['reference'],
                          refLon = parms['ref_longitude'],
                          refLat = parms['ref_latitude'])

The ``method`` parameter can be set to:

- ``'auto'`` - reference point choosen automatically in barycentre and under good coherence conditions

- ``'constraint'`` - reference point selected nearest to the given coordinates (e.g. of a corner reflector)

- ``'free'`` - datum-free network solution (equivalent to taking the overall average as reference)


Secondly, check the network conditioning under the defined reference by:

.. code:: python

   insarUtils.network_cond(data)

.. note::
   
   This step automatically removes isolated networks that might cause ill-conditioning of the estimation.


Finally, spatially integrate ambiguities:

.. code:: python

   insarUtils.spatialAmbiguity(data['network2'])


Solve network
-------------

.. code:: python

   insarUtils.solve_network(data['network2'])
   
Plot estimated parameters per points:

.. code:: python

   if parms['plot_flag']:
       plotUtils.plot_network_parms(data['network2'], 
                                    data['psc2'], 
                                    parms['outDir']+'/parms_1st.png')

Atmospheric phase screen (APS) estimation
=========================================

Estimate and remove APS on the solved first-order network and re-run all the previous steps:

.. code:: python

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

.. warning::
   
   APS estimation is a necesarry step for processing areas larger than approximatelly 10 x 10 km. 
   Note however that APS estimation can be biased if using very small AOI.
   The ``atmoWindow`` parameter is long-wavelength signal filter window length (in days).
   This parameter should reflect the temporal length and sampling of 
   the dataset. Too small values alias other (displacement) signals, 
   whereas too large values result in inefficient filtering and 
   consequently biased APS estimation.


Precise network solution
========================

Precise network solution includes APS correction and refined height-to-phase conversion factors:

.. code:: python

   insarUtils.getPreciseH2PH(data['psc2'], 
                             parms['outDir'],
                             parms['stackId'])
   insarUtils.solve_network_precise(data['network2'], 
                                    data['psc2'], 
                                    model = parms['model'])

Remove reference phase noise (RPN) and solve again:

.. code:: python

   insarUtils.remove_RPN(data['network2'],
                              plotFlag = parms['plot_flag'],
                              outDir = parms['outDir'])
   insarUtils.solve_network_precise(data['network2'], 
                                    data['psc2'], 
                                    model = parms['model'])
                                    
Plot refined network parameters:

.. code:: python

   if parms['plot_flag']:
       plotUtils.plot_network_parms(data['network2'], 
                                    data['psc2'], 
                                    parms['outDir']+'/parms_2nd.png')
                                    
                                    
Network densification
=====================

Densify first-order PS network by second-order PS candidates:

.. code:: python

   if parms['densify_flag']:
       insarUtils.getPreciseH2PH(data['psc_B'], 
                                 parms['outDir'],
                                 parms['stackId'])
       insarUtils.densify_network(data, 
                                  k = 3, 
                                  mode_thr = 1, 
                                  std_thr = parms['densify_stdThr'])

.. note::

   Parameter ``k`` defines number of nearest first-order PS to connect the second-order PSc to,
   ``mode_thr`` is a maximum of ambiguity misclosures in the solution, and 
   ``std_thr`` is a threshold on standard deviation of residuals (RMSE) of the second-order PSc.

Geocoding and export
====================     

Perform geocoding refinement (using estimated residual heights and sub-pixel positions):

.. code:: python

   insarUtils.fix_geocoding(data, parms)
   
Export results to standard CSV:

.. code:: python

   outCSV = parms['outDir']+'/insar_'+parms['stackId']+'.csv'
   insarUtils.HDF2csv(data, outCSV)


To visualize the results, see :ref:`visualization`.

.. note:: 

   It's easy to build your own data exporter. See for example custom modification ``insarUtils.HDF2csv_remotio()``, exporting to https://remotio.space CSV standard instead.

.. warning:: 

   If you desire to repeat the processing with different parameters, it is not necesarry to perform data extraction again.
   Simply call ``insarUtils.resetNetwork(data)`` and continue from :ref:`create_network`.
                     
