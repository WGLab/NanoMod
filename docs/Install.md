
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

## Test:
   To run `NanoMod`, users need to `source activate py27nanomod` and then go to the NanoMod folder.

	* The script to be run is bin/NanoMod.py: `python bin/NanoMod.py`
	
## Usage:
 For how to use them, please refer to [Usage](https://github.com/WGLab/NanoMod/blob/master/docs/Usage.md)

