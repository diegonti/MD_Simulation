# Python Scripts for Statistics and Visualization (in parallel)

Here are all the scripts used by the program to compute the statistics of the output results and to visualize them. Its the same version as the one in `serial/` but now it used the builtin `multiprocessing` module to paralelize the work among cores. Its set up in a way that each core will take care of une obverbable data and make its correspondent plots and stats.

## `stats.py` : 

Script to compute the statistics (average, $\langle x \rangle$, standard deviation (STD, $\sigma$) and autocorrelation time, $\tau$) of the raw results of the simulation. Uses the Block Average method to compute the stats in blocks, while the STD vs. block size ($m$) is fitted to a exponential function $\sigma(m) = a-be^{-m/\tau}$, where the $\sigma$ and $\tau$ are gathered.

Returns an output file, `name_stats.py`, with the statistical parameters ($\langle x \rangle$, $\sigma$ and $\tau$)), diffusion coeffitient ($D$) and fitting parameters ($a, b, c$) and also plots of the statistichal error ($\sigma$) vs. the block size ($m$).

The script accepts several flags for running:
```
$ python3 stats.py [-h] [-ip input_path] [-op output_path] [-s start] [-f final]
```
```
[-ip input_path] : Is mandatory and requires the path of the simulation log file. 

[-op output_path] : Folder where the plots are saved. Optional, defaults to `plots/`.

[-s start] : Start frame (int). Frame from which the output data is considered. Optional, defaults to the first frame.

[-f final] : Final frame (int). Frame up to which the output data is considered. Optional, defaults to the last frame.
```

## `visualization.py`

Script for plotting the results of the output data from the MD runs. Plots Energies, Temperature, Pressure, MSD, Total momentum and RDF.

The script can be executed with the same flags as before:  
```
$ python3 visualization.py [-h] [-t] [-ip input_path] [-op output_path] [-s start] [-f final]
```
Altough now this one is used if a GIF file of the trajectory wants to be created:
```
[-t] To create a trajectory.gif file of the simulation trajectory (for large trajectories it may take some time, VMD may be better).
```

<br>

## Help and Requeriments

Both files will give more information about the flags and usage with the flag `-h`:
```
$ python3 stats.py -h
$ python3 visualization.py -h
```
The requeriments needed for using this scripts are: `numpy`, `scypy`, `matplotlib`, `ase`, which are present in the `requirements.txt`, and can be installed with:
```
$ pip install -r scripts/requeriments.txt
or $ pip install module    (for each module)
```

For mor information/help, contact [Diego Ontiveros](https://github.com/diegonti).