.. _scr_prediction:

***************
SCR predicition
***************

This module performs Signal-to-Clutter ratio (SCR) prediction over candidate site for corner reflector installation. 

First, prepare ``coreg`` stacks using :ref:`snap_processing`.

SCR prediction is done in a crop around specified geodetic coordinates. It is driven by input in :download:`planning.json <../templates/planning.json>`:

.. code:: json

   {
     "id": "PEM2",
     "longitude": 18.340544,
     "latitude": 48.630174,
     "elevation": 250.952,
     "cropSize": 200,
     "RCS": 30,
     "stackDir": "/data/CR_Partizanske/DSC124/",
     "oversamplingFactor": 16,
     "outTiff": "/data/CR_Partizanske/SCR/SCR_DSC124.tiff"
   }

.. note::

   ``cropSize`` is radius in metres. ``RCS`` is expected reflector's Radar Cross Section in [dBm2]. 
   Expected RCS of chosen reflector type can be simulated using tools in the ``gecoris.crUtils`` module.


SCR simulation is performed by:

.. code:: shell

   python gecoris/scrSimul.py planning.json
   
Result is a GeoTIFF with predicted SCR over your AOI.

Example output
--------------

.. figure:: _static/predictedSCR_example.png
    :align: center
    :figwidth: 600px

    Maps of predicted SCR, given 1 m inner-leg-length square trihedral reflector (30 dBm2 RCS), 
    computed on the 1-year of Sentinel-1 time series over specific landslide area in Slovakia.


