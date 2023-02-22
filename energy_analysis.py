from readinputs import *

import matplotlib.pyplot as plt
import numpy as np
import os

#get log filename 
input_prefix = "repeats"+str(lattice_repeats)+"_a"+str(lattice_spacing)+"_R"+",".join([str(ele) for ele in radii.split(',')])+"_fracpos"+str(fraction_positive)+"_debye"+str(debye_length)+"_brushL"+",".join([str(ele) for ele in brush_lengths.split(',')])+"_charges-"+",".join([str(ele) for ele in surface_potentials.split(',')])+"_seed"+str(seed)
logfile= input_prefix + ".run.200000002.log"
print("file being analyzed:", logfile)
#get cluster size
clusterfile = input_prefix + '.clusteranalysis.txt'

if os.path.exists(logfile):

    f = open(logfile, 'r')
    data = f.readlines()

    timesteps = []
    energies = []
    temps = []

    for line in data[1:]:
        d = line.split('\t')
        timesteps.append(int(d[0]))
        energies.append(float(d[1]))
        temps.append(float(d[2].split('\n')[0]))

    #Plotting the desired quantities as a function of simulation timestep
    
    ##Potential energy vs timestep
    fig, ax = plt.subplots(figsize=(10,8))
    ax.plot(timesteps, energies)
    ax.set_xlabel('Timestep')
    ax.set_ylabel ('Potential Energy')
    if os.path.exists(clusterfile):
        cluster_size= float(open(clusterfile,'r').readlines()[0])
        ax.set_title(f"Potential Energy vs Timestep, fitness score: {cluster_size}")
    else:
        ax.set_title('Potential Energy vs Timestep')

    savename = os.path.splitext(logfile)[0]+'.energy.graph.jpg'
    plt.savefig(savename)
    plt.show()
    plt.close()
    
    ##Temperature vs timestep
    fig, ax = plt.subplots(figsize=(10,8))
    
    ax.plot(timesteps, temps)
    ax.set_xlabel('Timestep')
    ax.set_ylabel('Temperature')
    ax.set_title('Temperature vs Timestep')
    
    savename = os.path.splitext(logfile)[0]+'.temp.graph.jpg'
    plt.savefig(savename)
    plt.show()
    plt.close()
    