************
Introduction
************

Geodetic Corner Reflector (In)SAR Toolbox (GECORIS) is a Python framework for integrated analysis of satellite SAR time series of natural and artificial radar reflectors, such as corner reflectors or radar transponders. 
GECORIS can:

- perform geodetic positioning of artificial reflectors;
- assess the clutter level of a particular site before reflector installation;
- estimate the Radar Cross Section (RCS) time series to track reflector’s performance and detect outliers, e.g. due to debris accumulation or damage;
- estimate the Signal-to-Clutter Ratio (SCR) to predict the positioning precision and the InSAR phase variance;
- estimate the InSAR displacement time series of the reflector network.

GECORIS is partly built on the `Sentinels Application Platform (SNAP) <https://github.com/senbox-org/>`_. Currently supported are Sentinel-1 SLC measurements. 

GECORIS is described in this `paper <https://doi.org/10.3390/rs13050926>`_.


Copyright
---------

*(c) 2021 by Richard Czikhardt*

*Department of Theoretical Geodesy and Geoinformatics, Slovak University of Technology*

Support
-------

.. warning:: 

   THIS IS A RESEARCH CODE PROVIDED TO YOU "AS IS" WITH NO WARRANTIES OF CORRECTNESS. USE AT YOUR OWN RISK.

.. note:: 
   
   If you would like a support and debugging, please contact me via czikhardt.richard@gmail.com.


How to cite
-----------

Czikhardt R, van der Marel H, Papco J. GECORIS: An Open-Source Toolbox for Analyzing Time Series of Corner Reflectors in InSAR Geodesy. Remote Sensing. 2021; 13(5):926. https://doi.org/10.3390/rs13050926



License
-------

GNU GPL v3+
