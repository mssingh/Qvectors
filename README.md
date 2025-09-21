# Qvectors

### About

Simple package to take ERA5 data from gadi and plot Q-vectors and PV diagnostics for a given region. 
The script takes a range of times and makes an mp4 as well as providing each frame as a .png.

Works for the Australian region, but probably some bugs remain when trying to shift the box.

### Instructions

To run the package, clone the repository into a directory on `gadi`. Then run the python script
`plot_map.py`.

You will need the following dependencies:

1) The CLEX/Weather of the 21st CEntury conda environment

	a) Ensure you are a member of `hh5`

	b) Type the following commands
```
		$ module use /g/data/hh5/public/modules
		$ module load conda/analysis3-22.10 
```

2) Access to ERA5 datasets

	a) Ensure you are a member of `rt52`

	b) Ensure you are a member of `uc16`


Ensure that you request `storage=gdata/hh5+gdata/rt52+gdata/uc16` when asking for resources on `gadi`


Once you have run the code, it should produce some `.png` maps in the `Figures/` directory, and some 
`.mp4` animations in the `Animations/` directory. 

Take a look at `plot_map.py` to see how you can change the time and region for plotting. 

Unfortuantely, real time data is not available, so you have to plot historic data. I am working 
on getting access to ECMWF forecasts so I can make the plots with forecast data.



