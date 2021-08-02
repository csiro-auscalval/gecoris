.. _sentinel_download:

************************
Sentinel-1 data download
************************

Download Sentinel-1 SLC data in batches.
Sample input parameters file :download:`download.parms <../templates/download.parms>`:

.. code:: python

   { 
   "parameters" : {
       "downloadMethod" : "aria2c",	# "aria2c" / "sentinelsat"
       "dataPath" : "/data/SNT1/",
       "scihubUser" : "gmapit",
       "scihubPassword" : "gmapit.2020",
       "updateSwitch" : False  		# True/False
   },
   "AOI" : {				# either as corners of polygon or .geojson file
       "minLongitude" : 18.34054,
       "maxLongitude" : 18.34055,
       "minLatitude" : 48.63017,
       "maxLatitude" : 48.63018,
       "geojson" : "/data/Slovakia.json",
       "startDate" : '20200501',
       "endDate" : '',			# current date by default if left blank
       "swaths" : ('DSC124','ASC73')     # tuple of all swaths
   }}


Call download routine in shell session:

.. code:: shell

   python gecoris/batchDownload.py download.parms

.. note::

   Either `aria2c <https://aria2.github.io/>`_ or `sentinelsat <https://sentinelsat.readthedocs.io/en/stable/>`_ must be installed.
