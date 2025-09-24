# Qvectors

## About

Simple package to take ERA5 data from gadi and plot Q-vectors and PV diagnostics for a given region. 
The script takes a range of times and makes an mp4 as well as providing each frame as a .png.

Works for the Australian region, but probably some bugs remain when trying to shift the box.

## Requirements

To run the package, clone the repository into a directory on `gadi`. Then run the python script
`plot_map.py`.

You will need the following dependencies:

1) The CLEX/Weather of the 21st CEntury conda environment

	Ensure you are a member of `xp65`

2) Access to ERA5 datasets

	a) Ensure you are a member of `rt52`

	b) Ensure you are a member of `uc16`


## Instructions

You can use the scripts either by logging into gadi via ssh and running an interactive job on the
supercomputer, or using the Australian Research Environment (ARE).

	a) Log in to gadi using your username and password:

```
		ssh -Y uername@gadi.nci.org.au
```

	b) Type the following commands:
```
		module use /g/data/xp65/public/modules
		module load conda/analysis3
```

	c) Clone the repository into a directory on gadi, and change into that directory :
```
		git clone git@github.com:mssingh/Qvectors.git Qvectors
		cd Qvectors
```

	d) Begin an interactive job on gadi:
```
		qsub -I -X -P k10  -q normal -l storage=storage=gdata/xp65+gdata/rt52+gdata/uc16 -l walltime=8:00:00,mem=12000MB -l ncpus=1
```

		You may have to replace `k10` with a different project.

	e) You may now run the python script `plot_map.py` using
```
		python3 plot_map.py
```

You can edit the script using the editor of your choice. You may also use an interactive
python environment such as ipython.

You could also run the code in a JupyterLab instance on the ARE, or through VisualStudio.
Feel free to do this if you like, and you are comfortable with gadi.

Once you have run the code, it should produce some `.png` maps in the `Figures/` directory, and some 
`.mp4` animations in the `Animations/` directory. 

Take a look at `plot_map.py` to see how you can change the time and region for plotting. 

Unfortuantely, real time data is not available, so you have to plot historic data. I am working 
on getting access to ECMWF forecasts so I can make the plots with forecast data.



