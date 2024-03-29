# Installation

GECORIS is tested on GNU/Linux and partly on MS Windows 10. Recommended minimum setup:

* **Python 3.6**
* **SNAP 8.0**

We recommend using [conda](https://docs.conda.io/en/latest/index.html) environment to install gecoris along with SNAP and all Python dependencies. 

### 1. Install [miniconda](https://docs.conda.io/en/latest/miniconda.html)

### 2. Create Python environment 

For SNAP, Python =<3.6 is required:

```bash
conda create -n gecoris python=3.6
conda activate gecoris
```

### 3. Install Python dependencies

```bash
conda install numpy scipy matplotlib pandas shapely requests scikit-image h5py=2.10
conda install -c conda-forge gdal netcdf4
```

### 4. Install [SNAP](http://step.esa.int/main/download/snap-download/) including Sentinel-1 toolbox

*Note for GNU/Linux users:*

Install additional dependency:

```bash
sudo apt-get install libgfortran3 libgfortran5 python3-distutils
```

`libgfortran3` was discontinued on Ubuntu 20 based distros. Must be downloaded and installed manually from [here](https://packages.ubuntu.com/bionic/amd64/libgfortran3/download) with additional [gcc dependency](https://packages.ubuntu.com/bionic/amd64/gcc-6-base/download).

### 5. Configure SNAP with Python 

```bash
cd <snap-install-dir>/bin
./snappy-conf <python-exe>
```
where `<snap-install-dir>` is SNAP installation directory and `<python-exe>` is full path to conda Python binaries (`<miniconda-install-dir>/envs/gecoris/python` on Linux).

Install `snappy` package:

```bash
cd <snappy-dir>
python setup.py install
```

where `<snappy-dir>` is package directory (`$HOME/.snap/snap-python/snappy` on Linux).

More detailed instructions can be found [here](https://senbox.atlassian.net/wiki/spaces/SNAP/pages/50855941/Configure+Python+to+use+the+SNAP-Python+snappy+interface).

### 6. Download development version of gecoris and setup

```bash
git clone https://memorid@bitbucket.org/memorid/gecoris.git
```

Edit PATH environmental variables:

- on GNU/Linux, append to `.bashrc` file:

```bash
export gecoris_HOME=~/sw/gecoris
export PYTHONPATH=${PYTHONPATH}:${gecoris_HOME}
export PATH=${PATH}:${gecoris_HOME}/gecoris
```

- on WINDOWS in command prompt:

```
setx gecoris_HOME C:\sw\gecoris /m
setx PYTHONPATH "%PYTHONPATH%;%gecoris_HOME%" /m
setx PATH "%PATH%;%gecoris_HOME%\gecoris" /m
```


### Additional requirements

- For Sentinel-1 SLC download routines: [aria2c](https://aria2.github.io/) or [sentinelsat](https://sentinelsat.readthedocs.io/en/stable/)



For precise atmospheric delay corrections in positioning module

* CDS API is needed to auto-download ERA5 ECMWF data:

```
conda install -c conda-forge cdsapi
```

​		Then, follow instructions to setup:

​		https://cds.climate.copernicus.eu/api-how-to

- We use IONEX data from CDDIS. Due to new NASA data access policy, one needs to prepare `.netrc` file in his home dir (`~`):
```
machine urs.earthdata.nasa.gov login <your_earthdata_login> password <your_earthdata_pass>
```

Follow instructions here: 

https://wiki.earthdata.nasa.gov/display/EL/How+To+Access+Data+With+cURL+And+Wget