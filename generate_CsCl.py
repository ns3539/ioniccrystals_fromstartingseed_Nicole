#!/usr/bin/python
import argparse
import numpy as np
import os

parser=argparse.ArgumentParser()
parser.add_argument("--radiusN",default=105.0,help="Radius type 1 (default: %(default)s)",type=float)
parser.add_argument("--radiusP",default=85.0,help="Radius type 2 (default: %(default)s)",type=float)
parser.add_argument("-B","--brush_length",help="Brush length (default: %(default)s)",default=10,type=float)
parser.add_argument("-r","--repeats",help="Number of repeats in each direction (default: %(default)s)",default=5,type=int)
parser.add_argument("--outdir",default="CsCl")
args = parser.parse_args()
locals().update(vars(args))

#approx_min = radiusN+radiusP+2*brush_length
#large_size = 1.0
#small_size = min(radiusN,radiusP)/max(radiusN,radiusP)
#touching_dist = approx_min/(max(radiusN,radiusP))
#lattice_const = 2*touching_dist/np.sqrt(3)

approx_min=241.09
max_radius=max(radiusN,radiusP)
lattice_const=approx_min/max_radius

from ase.lattice.cubic import SimpleCubicFactory
class CsClFactory(SimpleCubicFactory):
    "A factory for creating CsCl"
    bravais_basis = [[0, 0, 0], [0.5, 0.5, 0.5], ]
    element_basis = (0,1)

fact = CsCl = Rocksalt = CsClFactory()


atoms = fact(directions=[[1,0,0], [0,1,0], [0,0,1]],
            size=(repeats,repeats,repeats), symbol=['N','P'], 
            latticeconstant=lattice_const)

atoms.write(os.path.join(outdir,"CsCl_rN%.0f_rP%.0f_r%i_lc%.2f_B%.0f.xyz"%(radiusN,radiusP,repeats,lattice_const*max_radius,brush_length)))
