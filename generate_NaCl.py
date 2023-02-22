#!/usr/bin/python
import argparse
import numpy as np
import os

parser=argparse.ArgumentParser()
parser.add_argument("--radiusP",default=65.0,help="Radius type 1 (default: %(default)s)",type=float)
parser.add_argument("--radiusN",default=150.0,help="Radius type 2 (default: %(default)s)",type=float)
parser.add_argument("-B","--brush_length",help="Brush length (default: %(default)s)",default=10.0,type=float)
parser.add_argument("--fudge",default=0,help="Multiple of brush length to add to the spacing (default:0)",type=float)
parser.add_argument("-r","--repeats",help="Number of repeats in each direction (default: %(default)s)",default=5,type=int)
parser.add_argument("--outdir",default="NaCl")
args = parser.parse_args()
locals().update(vars(args))

approx_min = 2*radiusN+2*radiusP+4*brush_length
#approx_min = radiusN+radiusP+2*brush_length
#approx_min=467.85

large_size = 1.0
max_radius = max(radiusN,radiusP)
min_radius=min(radiusN,radiusP)
small_size = min(radiusN,radiusP)/max_radius

bl_reduced= brush_length/(max_radius)
lattice_const = (approx_min) / max_radius 


particle_types=['N','P']
"""
if not lattice_const*np.sqrt(2)/2. > 2*large_size + 2*bl_reduced:
    print("Warning: lattice_const*sqrt(2)/2 is too small compared to 2R_max + 2*brush",lattice_const*np.sqrt(2)/2.,2*(large_size+bl_reduced))
else:
    print("lattice_const*sqrt(2)/2 is good compared to 2R_max+2*brush",lattice_const*np.sqrt(2)/2.,2*(large_size+bl_reduced))
"""
from ase.lattice.cubic import SimpleCubicFactory
class NaClFactory(SimpleCubicFactory):
    "A factory for creating NaCl"

    bravais_basis = [[0, 0, 0], [0, 0, 0.5], [0, 0.5, 0], [0, 0.5, 0.5],
                     [0.5, 0, 0], [0.5, 0, 0.5], [0.5, 0.5, 0],
                     [0.5, 0.5, 0.5]]
    element_basis = (0, 1, 1, 0, 1, 0, 0, 1)
    #element_basis = (1,0,0,1,0,1,1,0)

fact = NaCl = Rocksalt = NaClFactory()


atoms = fact(directions=[[1,0,0], [0,1,0], [0,0,1]],
            size=(repeats,repeats,repeats), symbol=['N','P'],  
            latticeconstant=lattice_const)

if outdir:
    atoms.write(os.path.join(outdir,"NaCl_rN%.0f_rP%.0f_r%i_B%.0f_lc%.2f_fudge%.1f.xyz"%(radiusN,radiusP,repeats,brush_length,lattice_const*max_radius,fudge)))

