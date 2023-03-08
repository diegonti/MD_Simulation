"""
Script for plotting the results of the output data from the MD runs.
Plots Energies, Temperature, Pressure, MSD and RDF.

Use: $ python3 visualization.py path start finish
although if any arguments is not present chooses from the default 
path, start,finish = ./plots/, 0,-1 

Diego Ontiveros
"""

import os
import sys
import time as cpu_time

import numpy as np
import matplotlib.pyplot as plt
to = cpu_time.time()


def parseUserArguments(arguments):
    """Parses the user input for the script."""

    path = "./plots/"
    start,finish = None,None

    for i,arg in enumerate(arguments):
        try: arguments[i] = int(arg)
        except ValueError: pass

    # When no inputs are given, stay with default
    if len(arguments) == 0: pass

    # When 1 argument is given see if its the path or start
    elif len(arguments) == 1:
        try: start = int(arguments[0])
        except ValueError: path = arguments[0]

    # When 2 arguments are given see if they are path+s or s+f
    elif len(arguments) == 2:
        if isinstance(arguments[0],str): path,start = arguments
        elif isinstance(arguments[1],str): start,path = arguments
        else: start,finish = arguments

    # When 3 arguments are given see where the path is
    elif len(arguments) == 3:
        if isinstance(arguments[0],str): path,start,finish = arguments
        elif isinstance(arguments[0],int): start, finish, path = arguments
    
    return path,start,finish


def makePlot(t,lines:list[np.ndarray],colours:list[str],labels:list[str],file_name:str,start=None,finish=None,save=True,**kwargs):
    """General function to make a plot of output line (or lines) and save the image.

    Parameters
    ----------
    `t` : x array of the plot. np.array.
    `lines` : List with the arrays of the lines to plot. When only one line, accepts a np.array.
    `colours` : List with the colours of each line. When only one line, accepts a string.
    `labels` : List with the labels of each line. When only one line, accepts a string.
    `file_name` : Name of the output image.
    `start` : Optional. Index from which iteration to start the plot. By default None.
    `finish` : Optional. Index at which iteration the plot ends. By default None.
    `save` : Optional. Bool to save or not the figure. By default True.
    `**kwargs`: Optional. General plot wkargs (xylim,figsize,xylabels).

    Returns
    -------
    `fig` : Drawn Figure.
    `ax` : Drawn Axis.
    """

    if not isinstance(lines,list): lines, colours, labels = [lines],[colours],[labels]

    fig,ax = plt.subplots(figsize=kwargs.get("figzise",(7,7)))

    for line,colour,label in zip(lines,colours,labels):
        plt.plot(t[start:finish],line[start:finish],c=colour,lw=0.5,label=label)

    ax.set_xlim(kwargs.get("xlim",(t[start:finish][0],t[start:finish][-1])))
    ax.set_ylim(kwargs.get("ylim",(None,None)))
    ax.set_xlabel(kwargs.get("xlabel","Time (ps)"))
    ax.set_ylabel(kwargs.get("ylabel",None))
    ax.legend()
    fig.tight_layout()
    if save: fig.savefig(path+file_name,dpi=600)

    return fig,ax


############################### MAIN PROGRAM ###############################

if __name__ == "__main__":

    # User input management
    arguments = sys.argv[1:]
    path,start,finish = parseUserArguments(arguments)
    try: os.mkdir(path)
    except FileExistsError: pass
    print("\nMaking Plots...")
    
    # Loading data from files
    dataT = np.loadtxt("output.dat",skiprows=0)     			# Thermodynamic data
    data = dataT.T                          					# Each parameter in a column
    t,E,Epot,Ekin,Tinst,P,MSD,p = data      					# Getting each parameter
    
    # dataRDF = np.loadtxt("RDF_out.dat").T
    # r,RDF = dataRDF

    # Energies Plot
    makePlot(t,[Ekin,Epot,E],["r","b","k"],["$E_{kin}$","$E_{pot}$","$E$"],
            file_name="energies.png",
            xlabel="Time (ps)", ylabel="Energy (kJ/mol)")

    # Temperature Plot
    makePlot(t,Tinst,"r","$T_{inst}$",
            file_name="temperature.png",
            xlabel="Time (ps)", ylabel="Temperature (K)")

    # Pressure Plot
    makePlot(t,P,"purple","P",
            file_name="pressure.png",
            xlabel="Time (ps)",ylabel="Pressure (Pa)")

    # MSD Plot
    start,finish = None,None
    figMSD,axMSD = makePlot(t,MSD,"green","MSD",
            file_name="MSD.png", save=False,
            xlabel="Time (ps)", ylabel="MSD ($\AA^2$)")
    a,b = np.polyfit(t[start:finish],MSD[start:finish],deg=1)     # linear fit to ax + b
    D = a/6  * 1e12 /1e20 # (Diffusion coeffitient)
    #! Show D in plot?

    axMSD.plot(t[start:finish],a*t[start:finish]+b,"k:",alpha=0.75)
    figMSD.savefig(path+"MSD.png",dpi=600)

    # RDF Plot
    # makePlot(r,RDF,"red","RDF",
    #          file_name="RDF.png",
    #          xlabel="RDF",ylabel="r ($\AA$)")
    # 


    # Trajectory #!(make this user input?)
    make_trajectory = False
    if make_trajectory:
        from ase import io
        if os.path.exists("trajectory.xyz"):
            try:
                frames = io.read("trajectory.xyz", index=":")
                io.write(path+"trajectory.gif", frames, interval=50)
                print(f"Trajectory.gif was created.")

            except:
                print("ase could not read the trajectory.xyz file. Make sure the format is correct (as .xyz).")
        else: print("Structure trajectory not found. Make sure the file was created.")

tf = cpu_time.time()
print(f"\nProcess finished in {tf-to:.2f}s.\n")