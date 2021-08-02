.. _snap_processing:

***************
SNAP processing
***************

This module contains Python wrappers for SNAP to extract Sentinel-1 bursts containing your Area of Interest (AOI) and co-register them into a stack of images within a common reference geometry.
Sentinel-1 data of your AOI should be already downloaded and accessible locally or via the network filesystem. For automatic download, see :ref:`sentinel_download`.

Burst extraction
----------------

For each stack covering your AOI, prepare input parameters file :download:`snap.parms <../templates/snap.parms>`:

.. code:: python
   
   { 
   "SNAP_parameters" : {
       "snapPath" : "/home/rc/sw/snap/bin/",          # path to SNAP binaries
       "repoPath" : "/home/rc/sw/gecoris/gecoris/",   # path to gecoris
       "ram" : 100000,                                # [MB]
       "cores" : 64                                   #
   },
   "Stack_parameters" : {
       "workDir" : '/home/rc/CR_Partizanske/DSC124/', # processing directory
       "dataDirs" : ('/home/rc/SNT1/DSC124/',),       # tuple of paths to SLC data
       "swaths" : ('IW1',),                           # tuple of all usable swaths
       "min_lon" : '18.34054',                        # Specify  AOI corners:
       "max_lon" : '18.34055',
       "min_lat" : '48.63017',
       "max_lat" : '48.63018',
       "startDate" : '20190730',
       "endDate" : '',                                # leave blank to use all SLC
       "masterDate" : ''                              # leave blank to select optimal
   }}
   
This wrapper calls SNAP to extract corresponding bursts and download precise orbits:

.. code:: shell

   python gecoris/batchPrepare.py snap.parms
   
Extracted SLC bursts are stored in ``slaves`` directory under your ``workDir`` specified in input parameters.

.. note::

   If your AOI spans multiple subswaths, you must specify all of them in ``swaths`` parameter as a comma-separated tuple.


Co-registration
---------------

For :ref:`insar_processing` and :ref:`scr_prediction` co-registered image stacks are required.

First, you must choose the *master* (geometric reference) image 
either automatically, minimizing temporal and geometric baselines:

.. code:: shell
   
   python gecoris/getMaster.py snap.parms
   
Or it can be selected manually and placed in the ``master`` directory under your ``workDir``.

Co-registration is then performed by:

.. code:: shell
   
   python gecoris/batchCoreg.py snap.parms


.. note::

   Current implementation uses single-master configuration, 
   geometry-based co-registration and ESD correction. Sub-swath wise.
