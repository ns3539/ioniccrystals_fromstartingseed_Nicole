from readinputs import * 

import numpy as np
import argparse
import os
import sys

#GA fitness function
def fitness_func(cluster):
    fitness = cluster 
    return fitness


outputprefix="repeats"+str(lattice_repeats)+"_a"+str(lattice_spacing)+"_R"+",".join([str(ele) for ele in radii.split(',')])+"_fracpos"+str(fraction_positive)+"_debye"+str(debye_length)+"_brushL"+",".join([str(ele) for ele in brush_lengths.split(',')])+"_charges-"+",".join([str(ele) for ele in surface_potentials.split(',')])+"_seed"+str(seed)
analysisfile = outputprefix+'.clusteranalysis.txt'

print("cluster file being analyzed:", analysisfile)

if os.path.exists(analysisfile):
    cluster_size= float(open(analysisfile,'r').readlines()[0])
    fitness_score = fitness_func(cluster_size)
    
    radii = [float(r) for r in radii.split(',')]
    
    brush_lengths = [float(x) for x in brush_lengths.split(',')]
    
    surface_potentials = [float(y) for y in surface_potentials.split(',')]

    save_data=[lattice_repeats, lattice_spacing]
    
    for r in radii:
        save_data.append(r)
        
    save_data.append(fraction_positive)
    save_data.append(debye_length)
    for x in brush_lengths:
        save_data.append(x)
    for y in surface_potentials:
        save_data.append(y)

    print("parameters:", save_data)
    print("fitness score:", fitness_score)
    
    save_data.append(fitness_score)
    
    save_data = np.array(save_data, dtype=object)
    
    
    #convert to numpy array to save data
    #save_data = np.array([params, fitness_score], dtype=object)

    #basedir = os.path(analysisfile).split('sc')
    #basedir = os.path.basename()
    basedir = os.getcwd().split('GA_generation')
    
    print(basedir)
    generation=str(basedir[1][0])
    print("generation:", generation)
    
    
    fitnessfile=os.path.join(basedir[0],f"generation{generation}.fitness.csv")
    
    os.makedirs(os.path.dirname(fitnessfile),exist_ok=True)
    #output_file_fitness = 'generation'+ generation + '.fitness.txt'
    #print(output_file_fitness)
    
    f= open(fitnessfile, 'a')
    f.write(str(save_data) + '\n')
    #np.savetxt(f, save_data) 
    f.close()
        
    #with open(file_data, 'a') as f:
        #f.write(str([params, fitness_score]) + '\n')
else:

    print("path does not exist")
    print("either cluster file has not yet been generated, or there was an error in generating it")
    
    failurefile=os.path.join(basedir[0],f"generation{generation}.failure.csv")
    os.makedirs(os.path.dirname(failurefile),exist_ok=True)
    #output_file_fitness = 'generation'+ generation + '.fitness.txt'
    #print(output_file_fitness)
    
    f= open(failurefile, 'a')
    f.write(str(outprefix) + '\n')
    #np.savetxt(f, save_data) 
    f.close()