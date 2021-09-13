.. _visualization:

******************
QGIS visualization
******************

Spatial visualization of InSAR results and time series plots are possible via `QGIS <https://qgis.org/en/site/>`_ plugin `insar2qgis <https://bitbucket.org/memorid/insar2qgis/src/master/>`_.
Please follow plugin installation instructions in its README file.

Currently, the plugin ingests only ``.shp`` files.
GECORIS contains utility function ``csv2shp.py``, which can be used to convert standard ``.csv`` output from :ref:`insar_processing` to a shapefile:

.. code:: shell

   python csv2shp.py <results>.csv <results>.shp
