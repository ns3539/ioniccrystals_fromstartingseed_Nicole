#from readinputs import *

import numpy as np
import argparse
import os
import sys


######################## supporting functions for GA algorithm ########################
# adapted from https://towardsdatascience.com/genetic-algorithm-implementation-in-python-5ab67bb124a6
# modified to impose parameter limits

def select_mating_pool(pop, fitness, num_parents):
    # Selecting the best individuals in the current generation as parents for producing the offspring of the next generation.
    parents = np.empty((num_parents, pop.shape[1]))
    for parent_num in range(num_parents):
        max_fitness_idx = np.where(fitness == np.max(fitness))
        max_fitness_idx = max_fitness_idx[0][0]
        parents[parent_num, :] = pop[max_fitness_idx, :]
        fitness[max_fitness_idx] = -99999999999
    return parents

def crossover(parents, offspring_size):
    offspring = np.empty(offspring_size)
    # The point at which crossover takes place between two parents. Usually it is at the center.
    crossover_point = np.uint8(offspring_size[1]/2)

    for k in range(offspring_size[0]):
        # Index of the first parent to mate.
        parent1_idx = k%parents.shape[0]
        # Index of the second parent to mate.
        parent2_idx = (k+1)%parents.shape[0]
        # The new offspring will have its first half of its genes taken from the first parent.
        offspring[k, 0:crossover_point] = parents[parent1_idx, 0:crossover_point]
        # The new offspring will have its second half of its genes taken from the second parent.
        offspring[k, crossover_point:] = parents[parent2_idx, crossover_point:]
    return offspring

def mutation(offspring_crossover, num_weights):
    # Mutation randomly selects one of the genes to change and applies a random change within appropriately specified limits
    for offspring in range(offspring_crossover.shape[0]):
        # The random value to be added to the gene.
        index = np.random.randint(0, int(num_weights)) #randomly select a gene index
        random_value = np.around(np.random.uniform(-1.0, 1.0, 1), 2) #random change to add to/subtract from gene value, rounded to 2 decimal places
        
        mutation = offspring_crossover[offspring, index] + random_value #proposed change
        
        test_score = test_mutation(mutation, index) #test the change
        if test_score == 0:
            #acccept the change
            offspring_crossover[offspring, index] = mutation#offspring_crossover[offspring, index] + random_value
        elif test_score == 1:
            #reject the change
            offspring_crossover[offspring, index] = offspring_crossover[offspring, index]
        #offspring_crossover[idx, 0] = offspring_crossover[idx, 0] + random_value
    return offspring_crossover
    
def test_mutation(proposed_mutation, parameter_index):
    #Function to test if proposed mutation is valid within appropriate limits. 
    #Like a fitness function, this is specific to each optimization problem and the nature of the parameters
    
    #lattice_spacing
    if parameter_index == 0: 
        if proposed_mutation>=7 and proposed_mutation<=10:
            return 0
        else:
            return 1
        
    #negative radius
    if parameter_index == 1:
        if proposed_mutation>=99 and proposed_mutation<=111:
            return 0
        else:
            return 1
    
    #positive radius
    if parameter_index == 2:
        if proposed_mutation>=80 and proposed_mutation<=90:
            return 0
        else:
            return 1
    
    #debye length
    if parameter_index == 3:
        if proposed_mutation>=4 and proposed_mutation<=6:
            return 0
        else:
            return 1
        

####################################################################################

inputfile = 'inputs.yaml'
f=open(inputfile)
data=f.readlines()
gen = int(data[32].split(':')[1][1])
print(gen)

#get params and fitness functions

filename=f"generation{gen}.fitness.csv"
f = open(filename)
data = f.readlines()
print("accessing fitness scores for generation", gen)

#to turn string back into list
new_data = []
for i in data:
    new_data.append(i.split('[')[1].split(']')[0])
    
list_data=[]
for i in new_data:
    list_data.append([float(ele) for ele in i.split(' ')])
    
print(list_data)

#define variables
gen_size = len(list_data) #number of input combinations
print("population size:", gen_size)
num_parents_mating=int(gen_size/5)  #number of solutions to be selected as parents in the next generation (the division by 5 is arbitrary)
print("number of solutions being selected as parents:", num_parents_mating)
num_weights=4  #number of parameters being varied

######indices for reference######

# lattice_repeats = 0
# lattice_spacing = 1
# radiusN = 2
# radiusP = 3
# frac_positive = 4
# debye_length = 5
# brush_lengths = 6, 7
# surface_potentials = 8, 9


input_params=[]
fitness = []

for i in list_data:
    params = [i[1], i[2], i[3], i[5]]
    fitness_score = i[10]
    input_params.append(params)
    fitness.append(fitness_score)
    
input_params = np.array(input_params)
 
# Selecting the best parents in the population for mating.
parents = select_mating_pool(input_params, fitness, 
                                    num_parents_mating)
                                    
print("parents:", parents)

# Generating next generation using crossover.

offspring_crossover = crossover(parents, offspring_size=(len(input_params)-parents.shape[0], num_weights))

# Adding some variations to the offsrping using mutation.
offspring_mutation = mutation(offspring_crossover, num_weights)# Creating the new population based on the parents and offspring.

input_params[0:parents.shape[0], :] = parents
input_params[parents.shape[0]:, :] = offspring_mutation

print("optimized parameters:", input_params) #return optimized inputs to feed into next generation
print("size of next generation:", len(input_params))

#Save output file containing analysis data in the same directory as simulation trajectory (for reference/records)

output_filename=f"generation{gen}.GAanalysis.csv"

#d=np.array([input_params],dtype=object)
#print(d)
#d.dump(output_file_GAanalysis)

#os.makedirs(os.path.dirname(output_filename),exist_ok=True)
    
f= open(output_filename, 'w')
f.write(str(input_params) + '\n')
#np.savetxt(f, save_data) 
f.close()

#increase generation number
gen +=1

#Load parameters into .sh script to run MD simulations for next generation

f = open('vary_parameters_run.sh')
lines = f.readlines()
savelines = lines[0:44]


#get values for unchanged parameters

lattice_repeats = int(list_data[0][0])
frac_positive = float(list_data[0][4])
brush_lengths = str(list_data[0][6]) +',' + str(list_data[0][7])
surface_potentials = str(list_data[0][8]) +',' + str(list_data[0][9])
seed = 1


#updated parameters 

for i in input_params:
    
    ls = float(i[0])
    rn = str(i[1])
    rp = str(i[2])
    dl = str(i[3])
    r = rn+','+rp

    #savelines.append(f'submit_simulation ${lattice_repeats} ${ls} $"{r}" ${frac_positive} ${dl} $"{brush_lengths}" $"{surface_potentials}" ${seed}\n')
    new_job = f'submit_simulation {lattice_repeats} {ls} "{r}" {frac_positive} {dl} "{brush_lengths}" "{surface_potentials}" {seed} {gen}\n'
    if new_job not in savelines:
        savelines.append(new_job)
    
savefile = f'vary_parameters_run_future.sh'
with open(savefile, 'w') as f:
    for line in savelines:
        f.write(line)
        
#prepare other .sh files

#cluster analysis
f = open('run_cluster_analysis.sh')
lines = f.readlines()
savelines = lines[0:46]

#get values for unchanged parameters

lattice_repeats = int(list_data[0][0])
frac_positive = float(list_data[0][4])
brush_lengths = str(list_data[0][6]) +',' + str(list_data[0][7])
surface_potentials = str(list_data[0][8]) +',' + str(list_data[0][9])
seed = 1


for i in input_params:
    
    ls = float(i[0])
    rn = str(i[1])
    rp = str(i[2])
    dl = str(i[3])
    r = rn+','+rp
    
    new_job = f'analyze_runs {lattice_repeats} {ls} "{r}" {frac_positive} {dl} "{brush_lengths}" "{surface_potentials}" {seed} {gen}\n'
    if new_job not in savelines:
        savelines.append(new_job)
    
    #savelines.append(f'analyze_runs {lattice_repeats} {ls} "{r}" {frac_positive} {dl} "{brush_lengths}" "{surface_potentials}" {seed} {gen}\n')

savefile = f'run_cluster_analysis_future.sh'
with open(savefile, 'w') as f:
    for line in savelines:
        f.write(line)
        
#energy analysis
f = open('run_energy_analysis.sh')
lines = f.readlines()
savelines = lines[0:45]

#get values for unchanged parameters

lattice_repeats = int(list_data[0][0])
frac_positive = float(list_data[0][4])
brush_lengths = str(list_data[0][6]) +',' + str(list_data[0][7])
surface_potentials = str(list_data[0][8]) +',' + str(list_data[0][9])
seed = 1

for i in input_params:
    
    ls = float(i[0])
    rn = str(i[1])
    rp = str(i[2])
    dl = str(i[3])
    r = rn+','+rp
    
    new_job = f'make_graphs {lattice_repeats} {ls} "{r}" {frac_positive} {dl} "{brush_lengths}" "{surface_potentials}" {seed} {gen}\n'
    if new_job not in savelines:
        savelines.append(new_job)
    #savelines.append(f'make_graphs {lattice_repeats} {ls} "{r}" {frac_positive} {dl} "{brush_lengths}" "{surface_potentials}" {seed} {gen}\n')

savefile = f'run_energy_analysis_future.sh'
with open(savefile, 'w') as f:
    for line in savelines:
        f.write(line)
        
# fitness function file
f = open('run_GA_fitness.sh')
lines = f.readlines()
savelines = lines[0:47]

#get values for unchanged parameters

lattice_repeats = int(list_data[0][0])
frac_positive = float(list_data[0][4])
brush_lengths = str(list_data[0][6]) +',' + str(list_data[0][7])
surface_potentials = str(list_data[0][8]) +',' + str(list_data[0][9])
seed = 1


for i in input_params:
    
    ls = float(i[0])
    rn = str(i[1])
    rp = str(i[2])
    dl = str(i[3])
    r = rn+','+rp

    #new_job = f'submit_simulation {lattice_repeats} {ls} "{r}" {frac_positive} {dl} "{brush_lengths}" "{surface_potentials}" {seed} {gen}\n'
    #if new_job not in savelines:
        #savelines.append(new_job)
    savelines.append(f'get_fitness {lattice_repeats} {ls} "{r}" {frac_positive} {dl} "{brush_lengths}" "{surface_potentials}" {seed} {gen}\n')

savefile = f'run_GA_fitness_future.sh'
with open(savefile, 'w') as f:
    for line in savelines:
        f.write(line)
        
# update generation in input file

inputfile = 'inputs.yaml'
f=open(inputfile, 'r')
data=f.readlines()
f.close()

g=open(inputfile, 'w')
data[32] = f'gen: {gen}'
for line in data:
    g.write(line)

       
#GA algorithm - only need to update generation number

f = open('GA_algorithm.sh')
lines = f.readlines()

lines[22] = f'for gen in {gen}; do\n'

savefile = f'GA_algorithm.sh' #overwrite the original file since only one line is changing 
with open(savefile, 'w') as f:
    for num, line in enumerate(lines):
        f.write(line)
