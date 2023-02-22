import numpy as np
import sys
import re
import os
import os.path
import argparse
import yaml
import copy


with open ('inputs.yaml') as f:
        data = yaml.load(f,Loader=yaml.FullLoader)


parser = argparse.ArgumentParser()
parser.add_argument('--lattice_repeats', action= "store",nargs ="*", default=10, dest = "lattice_repeats", type = int)
parser.add_argument('--lattice_spacing', action= "store",nargs ="*", default=5.0, dest = "lattice_spacing", type = float)
parser.add_argument('--lattice_type',action= "store",nargs ="*", default="sc", dest = "lattice_type", type = str)
parser.add_argument('--radii',action="store", nargs="*", default='65,150',dest="radii",type=str)
parser.add_argument('--fraction_positive', action= "store", nargs="*", default=0.5, dest="fraction_positive", type=float)
parser.add_argument('--debye_length',action="store", nargs="*", default=5.0, dest="debye_length",type=float)
parser.add_argument('--brush_lengths',action="store", nargs="*", default='10.0,10.0', dest="brush_lengths",type=str)
parser.add_argument('--surface_potentials',action="store", nargs="*", default='50.0,50.0', dest="surface_potentials",type=str)
parser.add_argument('--seed',action="store", nargs="*", default=10, dest="seed",type=int)
parser.add_argument('--gen',action="store", nargs="*", default=0, dest="gen",type=int)
print("")
print("Adjusted Values")
print("--------------------")
results = parser.parse_args()
print("lattice_repeats:",results.lattice_repeats)
print("lattice_spacing:",results.lattice_spacing)
print("lattice_type:",results.lattice_type)
print("radii:",results.radii)
print("fraction_positive:",results.fraction_positive)
print("debye_length:",results.debye_length)
print("brush_lengths:",results.brush_lengths)
print("surface_potentials:",results.surface_potentials)
print("seed:",results.seed)
print("generation:", results.gen)

result_list = []

for gen in results.gen:
    for lattice_repeats in results.lattice_repeats:
        for lattice_spacing in results.lattice_spacing:
            for radii in results.radii: 
                for fraction_positive in results.fraction_positive: 
                    for debye_length in results.debye_length:
                        for brush_lengths in results.brush_lengths:
                            for surface_potentials in results.surface_potentials:
                                for seed in results.seed:
    
                                    lattice_type=results.lattice_type[0]
                                    new_data = copy.copy(data)                        
                                    new_data['lattice_repeats']=lattice_repeats
                                    new_data['lattice_spacing']=lattice_spacing
                                    new_data['lattice_type']=lattice_type
                                    new_data['radii']=radii
                                    new_data['fraction_positive']=fraction_positive
                                    new_data['debye_length']=debye_length
                                    new_data['brush_lengths']=brush_lengths
                                    new_data['surface_potentials']=surface_potentials
                                    new_data['seed']=seed
                                    new_data['gen']=gen
    
                                    path=os.path.join(os.getcwd(),"GA_generation"+str(gen),str(lattice_type),"repeats"+str(lattice_repeats),"latticespacing"+str(lattice_spacing),"radii"+str(radii),"fraction_positive"+str(fraction_positive),"debye_length"+str(debye_length),"brush_lengths"+str(brush_lengths),"surface_potentials"+str(surface_potentials),"seed"+str(seed))
    
                                    file_data=os.path.join(path,"inputs.yaml")
                                    os.makedirs(os.path.dirname(file_data),exist_ok=True)
            
                                    with open (file_data, 'w') as outfile_data:
                                        yaml.dump(new_data, outfile_data, default_flow_style= False)
                                    print("Updated Yaml Input File generated.")
    
                                    print("*******************************************************************************************")
    
                                    print("Full path to run directory:")
                                    print(os.path.dirname(file_data))

