"""
Script to compute the statistics (mean, std and autocorrelation time) of the raw results of the simulation.
Uses the Block Average method to compute the average <x> and the std \sigma.
The STD vs. block size is fitted to a exponential function f(x) = a-b*exp(-x/tau).
Returns an output file with the stats and plots of the statistichal error (STD) vs. the block size.

Use with: $ python3 stats.py -ip input_path -op output_path -s start -f final
although if any arguments is not present chooses from the default 
path, start,finish = ./plots/, 0,-1 

Diego Ontiveros
"""
import os
import warnings
import time as cpu_time
from argparse import ArgumentParser, Namespace
import multiprocessing as mp

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
# warnings.filterwarnings("ignore", category=RuntimeWarning)
# warnings.filterwarnings("ignore", category=sp.optimize.OptimizeWarning)

to = cpu_time.time()

###################### FUNCTION DEFINITION ######################

def blockAverage(data, maxBlockSize=None):
	"""Computes the block average of a timeseries "x", and 
	provides error bounds for the estimated mean <x>. 

	Parameters
	--------------------
	`data`  : Time series array of an observable X
	`maxBlocksize` : Maximum number of observations/block

	Returns
	--------------------
	`m_points` : Points used (observations/block) for computing averages
	`blockVar` : Variances for each blocksize array 
	`blockMean` : Variances for each blocksize array 	
	"""
 
	Nobs = len(data)           # total number of observations in data
	minBlockSize = 1           # min 1 observation/block
 
	if maxBlockSize is None: maxBlockSize = int(Nobs/4)   # max: 4 blocks (otherwise can't calc variance)
  
	# m_points = 2**n until being less of the inputed maxblocksize
	power = np.arange(int(np.log(maxBlockSize)/np.log(2)))
	m_points = 2**power

	NumBlocks = len(m_points)   				# total number of block sizes

	blockMean = np.zeros(NumBlocks)             # mean for each blocksize
	blockVar  = np.zeros(NumBlocks)             # variance associated with each blockSize


	# Loop for all considered blocksizes (m)
	for k,m in enumerate(m_points):

		Nblock    = int(Nobs/m)               # Number of blocks
		obsProp   = np.zeros(Nblock)          # Container for parcelling block 

		# Loop to chop datastream into blocks and take average
		for i in range(1,Nblock+1):
			
			i1 = (i-1) * m
			i2 =  i1 + m
			obsProp[i-1] = np.mean(data[i1:i2])

		blockMean[k] = np.mean(obsProp)
		blockVar[k]  = np.var(obsProp)/(Nblock - 1)
	
	return m_points, blockVar, blockMean


def plotBlockAverage(m_points,blockVar,blockMean,fit_params=None,save=True,save_name="BlockAverage.jpg",label="x"):
	"""
	Returns a plot of the square root of the block variance and the block means as a function of block size m.

	Parameters
	---------------
	`m_points` : array with the number of samples/block used.
	`blockVar` : array with the block variances gathered.
	`blockMean` : array with the block means gathered.
	`fit_params`: (optional) optimized parameters of the fitting function a-b*exp(-m/t).
	`save`: (optional) Boolean to save or not the image. Defauls to True.
	`save_name`: (optional) File name of the saved image.
	`label` : (optional) Formatted label of the obserbable. In form: Observable(Units).
	"""
	# Getting label and units for the axis labels
	observable_label = label.split("(")[0]
	observable_units = label.split(")")[0].split("(")[-1]

	fig,(ax1,ax2) = plt.subplots(1,2, figsize = (10,5))

	# Plotting the STD points for each block (m)
	ax1.scatter(m_points, np.sqrt(blockVar),marker="x",c="k",lw=1)

	# Plots the fitting to the STD points
	if fit_params is not None:
		x = np.linspace(m_points[0],m_points[-1])
		ax1.plot(x,fit_function(x,*fit_params),"k:",label="fit")

	# Labels and Legend
	ax1.set_xlabel('m')
	ax1.set_ylabel(f'$\sigma_m$ ({observable_label}) ({observable_units})')
	ax1.legend()

	# Plotting the mean points for each block with the STD as errorbars
	ax2.errorbar(m_points, blockMean, yerr=np.sqrt(blockVar),fmt="k.-",capsize=5)
	ax2.set_ylabel(f'$\\langle${observable_label}$\\rangle$ ({observable_units})')
	ax2.set_xlabel('m')

	fig.tight_layout()

	if save: fig.savefig(opath+save_name,dpi=600)


def fit(fit_function,m_points,blockSTD):
	"""Fitting the computed block statistichal errors to the fit_function given"""
	params,cov = sp.optimize.curve_fit(fit_function,m_points,blockSTD)
	return params,cov

def fit_function(x,a,b,tau):
	"""Fit function for the statistichal error points."""
	return a-b*np.exp(-x/tau)


def write(line,file):
	"""Prints and writes line in the terminal and output file."""
	print(line,flush=True)
	with open(file,"a",encoding="utf-8") as outFile: outFile.write(line+"\n")

def title(text,before=15,after=15,separator="-",head=2,tail=1):
    """Prints text in a title-style."""
    separator = str(separator)
    line = "\n"*head+separator*before+text+separator*after+"\n"*tail
    return line


def individualStats(input_data): 
	"""The individual work that an individual worker will do. Function that will be paralelized."""

	xi,data_labels_i,data_labels_f_i = input_data
	n_tot = len(xi)
	
	# Calculate Block Averages (Binning)
	m_points_i,blockVar_i,blockMean_i = blockAverage(xi[start:finish],maxBlockSize=int(n_tot/100))

	# Calulate fitting for the Block_STD computed
	params_i, cov_i = fit(fit_function,m_points_i,np.sqrt(blockVar_i))      # Fitting --> Getting optimized params
	x_m = np.linspace(m_points_i[0],m_points_i[-1])                         # linspace for the fitted function (smoother)
	fitted_sigma = fit_function(x_m,*params_i)                              # values for the fitted function
	# params.append(params_i)
	
	# Computing statistichal parameters
	mean = blockMean_i.mean()                                               # Calculation of the total mean (average of BockMeans)
	std = fitted_sigma[-1]                                                  # Correlated STD will be the one at large BlockSizes (m) (last values of fitted function, plateau)
	tau = params_i[-1]														# Autocorrelation time (from fitting)
	
	# Writing and Plotting results
	observable_label = data_labels_i.split("(")[0].strip()
	write(data_format.format(data_labels_i,mean,std,tau),file=outFile)
	plotBlockAverage(m_points_i,blockVar_i,blockMean_i,fit_params=params_i,label=data_labels_f_i, save_name=f"BlockAverage_{observable_label}.png")

	return params_i

###################### MAIN PROGRAM ######################

if __name__=="__main__":
	print()

	# User input management
	parser = ArgumentParser(description="Script to compute the statistics (mean, std and autocorrelation time) of the raw results of the simulation.\
	Uses the Block Average method to compute the average <x> and the std \sigma. \
	The STD vs. block size is fitted to a exponential function f(x) = a-b*exp(-x/tau).\
	Returns an output file with the stats and plots of the statistichal error (STD) vs. the block size.")
	parser.add_argument("-ip","--ipath",help="Input path (str). File name where the simulation data is.", type=str)
	parser.add_argument("-op","--opath",help="Output path (str). Folder name where the plots/output will be created. Defaults to './plots/'.",default="./plots/", type=str)
	parser.add_argument("-s","--start",help="Start frame (int). Frame from which the output data is considered. Defaults to the first frame.",default=None, type=int)
	parser.add_argument("-f","--final",help="Final frame (int). Frame up to which the output data is considered. Defaults to the last frame.",default=None, type=int)
	parser.add_argument("-t","--None",help="Just for mantaining the same arguments between scripts.",default=None, type=int)

	args: Namespace = parser.parse_args()
	ipath,opath,start,finish = args.ipath,args.opath,args.start,args.final
	if not opath.endswith("/"): opath += "/"

	sim_name = ipath.split("_log")[0]

	try: os.mkdir(opath)
	except FileExistsError: pass

	# Loading data from test files
	dataT = np.loadtxt(ipath,skiprows=1)     			# Thermodynamic data
	data = dataT.T                          			# Each parameter in a column
	t,E,Epot,Ekin,Tinst,P,MSD,p = data      			# Getting each parameter

	# Defining data labels 
	data_labels = ["E (kJ/mol)","Epot (kJ/mol)","Ekin (kJ/mol)","T (K)","P (Pa)","MSD (A^2)","Pt (kg*m/s)"]
	data_labels_f = ["E(kJ/mol)","$E_{pot}$(kJ/mol)","$E_{kin}$(kJ/mol)","T(K)","P(Pa)","MSD($\\AA^2$)","$P_T$(kg*m/s)"]

	# Creating output file
	outFile = sim_name+"_stats.log"
	if os.path.exists(outFile): os.remove(outFile)

	# Formatted table columns and title
	titles_format =" {:^23} | " * (4 + 0)
	data_format =" {:^23} | " + " {:^23.15e} | " * (3 + 0)
	write("",file=outFile)
	write("Statistichal parameters for each computed observable:\n",file=outFile)
	write(titles_format.format("Observable","Average, <x>", "STD, \u03C3", "Autocorrelation time, \u03C4"),file=outFile)
	write(title("",before=23*4,head=0,tail=0),outFile)

	# Paralel Statistics
	input_data = zip(data[1:],data_labels,data_labels_f)
	pool = mp.Pool(len(data_labels))
	params = pool.map(individualStats,input_data)
	pool.close()
	pool.join()


	# Difussion Coeffitient (linear fit of MSD to y= ax + b)
	(a,b),D_residual,extra1,extra2,extra3 = np.polyfit(t[start:finish],MSD[start:finish],deg=1,full=True)
	D = a/6  * 1e12 /1e20 # (Diffusion coeffitient)
	write(f"\nEstimated Difussion Coefitient (A^2/ps): D = {D}",file=outFile)

	# Writing fitting parameters
	write("\n",file=outFile)
	write("Fitting parameters for block STD fitting to the function f(x) = a-b*exp(-x/tau):\n",file=outFile)
	write(titles_format.format("Observable","a", "b", "\u03C4"),file=outFile)
	write(title("",before=23*4,head=0,tail=0),outFile)
	for i,param in enumerate(params):
		write(data_format.format(data_labels[i], *param),file=outFile)


	# Ending
	tf = cpu_time.time()
	write(f"\n\nProcess finished in {tf-to:.2f}s.\n",file=outFile)
	write("Have a nice day :) \n",file=outFile)
