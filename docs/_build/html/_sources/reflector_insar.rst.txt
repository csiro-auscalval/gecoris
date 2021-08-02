.. _reflector_insar:

********
CR InSAR
********

This module performs InSAR time series analysis
in a network of corner reflectors plus surrounding coherent point scatterers (PS).
It is following the steps in :ref:`insar_processing`.

This module assumes you have already prepared co-registered SLC stacks of your AOI using SNAP, see :ref:`snap_processing`.

.. note::

   It is recommended to perform InSAR analysis :ref:`insar_processing` of your AOI manually before using the automatic routines in this module.


Three-tracks example
--------------------

Initialize already processed reflectors:

.. code:: python

   from gecoris import ioUtils

   station_log = '/data/GUDS/CR_Prievidza/CR_Prievidza.csv'
   stacks_log = '/data/GUDS/CR_Prievidza/input_stacks_Prievidza.csv'
   process_dir = '/data/GUDS/CR_Prievidza/gecoris'

   stations, refl_flag = ioUtils.load_stations(station_log, process_dir)

Define full paths to the results of InSAR analysis:

.. code:: python

   hdf_files = [
       '/data/GUDS/CR_VYCHOD/insar_Presov/insar_ASC175.hdf5',
       '/data/GUDS/CR_VYCHOD/insar_Presov/insar_DSC51.hdf5',
       '/data/GUDS/CR_VYCHOD/insar_Presov/insar_DSC124.hdf5']

   stack_ids = [
       's1_asc175',
       's1_dsc51'.
       's1_dsc124]

Extract InSAR time series for reflectors:

.. code:: python

   for i in range(len(stations)):
       for stackId, hdf in zip(stack_ids, hdf_files):
           stations[i].add_insar_ts(hdf, stackId)

Finally, perform LOS decomposition into vertical and east-west components and plot the time series:

.. code:: python

   decomp = True
   for i in range(len(stations)):
       # plot TS:
       stations[i].plot_insar_ts(out_dir = '/data/GUDS/CR_Prievidza/insar_Prievidza')
       if decomp:
           stations[i].decomp_insar_ts()
           stations[i].plot_decomp_ts(
               out_dir = '/data/GUDS/CR_Prievidza/insar_Prievidza')


.. figure:: _static/HRD-KU-1_s1_asc175_insar.png
    :align: center
    :figwidth: 600px
    
    Line-of-sight (LOS) displacement time series of double-reflector from three overlapping tracks.

.. figure:: _static/HRD-KU-1_insar_decomp.png
    :align: center
    :figwidth: 600px
    
    Vertical and east-west displacement time series of double-reflector transformed from LOS measurements from 
    three overlapping tracks.


Automatic routines
------------------

One should first use :ref:`reflector_monitoring` module to perform reflectivity time series analysis of the network of CR.

Then, all required parameters are automatically parsed from respective ``gecoris.parms`` and analysis is performed for all available stacks using:

.. code:: shell

   python gecoris/CR_insar.py gecoris.parms


.. warning::

   It is highly recommended to modify default ``insar.parms`` parameters in a second iteration for your AOI's specifics.


All data and results of InSAR processing are sequentially stored in binary HDF and human-readable CSV formats.


