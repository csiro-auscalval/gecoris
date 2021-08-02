# GECORIS
Geodetic Corner Reflector (In)SAR Toolbox - a Python framework for analyzing satellite SAR time series of artificial radar reflectors, such as corner reflectors or radar transponders. 
GECORIS can:

* perform geodetic positioning of artificial reflectors;
* assess the clutter level of a particular site before reflector installation;
* estimate the Radar Cross Section (RCS) time series to track reflector’s performance and detect outliers, e.g. due to debris accumulation or damage;
* estimate the Signal-to-Clutter Ratio (SCR) to predict the positioning precision and the InSAR phase variance;
* estimate the InSAR displacement time series of the reflector network.

GECORIS is partly built on the [Sentinels Application Platform (SNAP)](https://github.com/senbox-org/). Currently supported are Sentinel-1 SLC measurements. 

GECORIS is described in this [paper](https://doi.org/10.3390/rs13050926).

>*Copyright (c) 2021 by Richard Czikhardt*

> *Dept. of Theoretical Geodesy and Geoinformatics, Slovak University of Technology*



## Notice

THIS IS A RESEARCH CODE PROVIDED TO YOU "AS IS" WITH NO WARRANTIES OF CORRECTNESS. USE AT YOUR OWN RISK.

If you would like a support/debugging, please contact me via czikhardt.richard@gmail.com

Contributions are much appreciated!

*Richard Czikhardt* 



## [Installation](./installation.md)



## Usage

- ### [Data preparation & coregistration module (SNAPpy)](./data_preparation_module.md)


- ### [Planning module](./planning_module.md)

- ### [Reflector (network) monitoring module](./reflector_monitoring_module.md)

- ### [InSAR module](./InSAR_module.md)



## Cite this work

Czikhardt R, van der Marel H, Papco J. GECORIS: An Open-Source Toolbox for Analyzing Time Series of Corner Reflectors in InSAR Geodesy. Remote Sensing. 2021; 13(5):926. https://doi.org/10.3390/rs13050926



## License

GNU GPL v3+
