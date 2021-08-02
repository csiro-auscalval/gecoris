### Data preparation & coregistration module (SNAPpy):

This module contains Python wrappers for SNAP to download Sentinel-1 data, extract Area of Interest (AOI) and co-register stack of images.



#### 1. Download Sentinel-1 SLC data

Sample input parameters file [*download.parms*](./templates/download.parms):

```Python
{ 
"parameters" : {
    "downloadMethod" : "aria2c",	# "aria2c" / "sentinelsat"
    "dataPath" : "/data/SNT1/",
    "scihubUser" : "gmapit",
    "scihubPassword" : "gmapit.2020",
    "updateSwitch" : False  		# True/False
},
"AOI" : {							# either as corners of polygon or .geojson file
    "minLongitude" : 18.34054,
    "maxLongitude" : 18.34055,
    "minLatitude" : 48.63017,
    "maxLatitude" : 48.63018,
    "geojson" : "/data/Slovakia.json",
    "startDate" : '20200501',
    "endDate" : '',					# current date by default if left blank
    "swaths" : ('DSC124','ASC73')   # tuple of all swaths
}}
```

**You have to have [aria2c](https://aria2.github.io/) or [sentinelsat](https://sentinelsat.readthedocs.io/en/stable/) installed.*

Call routine:

```bash
python gecoris/batchDownload.py download.parms
```



#### 2. Prepare AOI bursts

For each AOI stack, prepare input parameters file (template in [*snap.parms*](./templates/snap.parms)):

```python
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
```

Prepare extracted AOI SLC bursts with precise orbits in `slaves` directory:

```bash
python gecoris/batchPrepare.py snap.parms
```

Prepared SLC bursts in `slaves` directory are ready for time series analysis by [Reflector (network) monitoring module.](./reflector_monitoring_module.md)



#### 3. -optional- Co-registration

For [SCR prediction](./planning_module.md) or further [InSAR analysis](./InSAR_module.md), co-registered image stacks (i.e. transformed to common reference image) are required.

Choose optimal master (geometric reference) and prepare `master` directory:

```bash
python gecoris/getMaster.py snap.parms
```

Co-registration:

```bash
python gecoris/batchCoreg.py snap.parms
```

*- using single-master configuration, geometry-based co-registration and ESD correction. Sub-swath wise.*

