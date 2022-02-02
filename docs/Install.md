
## Download the tool
   * git clone https://github.com/WGLab/NanoMod

## Install dependent packages

### Option 1: 
  Users can simple run the commands below to install the dependent packages
  ```
  cd NanoMod
  conda env create -f env.py27nanomod.yml
  ```
  
### Option 2:
  Optionally, you can also install the dependent packages below by yourself.

  Prerequisites:

	* Python 2.7
	* python packages:
	        + h5py
		+ scipy
		+ numpy
	* BWA MEM
	* HDF5
	* rpy2
	* R package:
	        + ggplot2
	        + gridExtra
	        + scales


You can run the commands below line by line to manually set up the conda environment and install the required packages:
```
conda create -n py27nanomod python=2.7.14
conda activate py27nanomod 
conda install -c conda-forge readline=6.2 
conda install -c bioconda samtools==1.9 bwa==0.7.17
conda install -c r rpy2==2.7.0 r-ggplot2==2.1.0 r-gridExtra==2.0.0 r-scales==0.3.0
conda install -c anaconda h5py==2.7.1 hdf5==1.10.1 scipy==1.2.1 numpy==1.15.4
conda install -c biobuilds r-cowplot==0.9.1
```

## Test:
   To run `NanoMod`, users need to `source activate py27nanomod` and then go to the NanoMod folder.

	* The script to be run is bin/NanoMod.py: `python bin/NanoMod.py`
	
## Usage:
 For how to use them, please refer to [Usage](https://github.com/WGLab/NanoMod/blob/master/docs/Usage.md)

