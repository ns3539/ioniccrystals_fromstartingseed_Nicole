import time
import hoomd
import hoomd.md
import numpy as np
import os
import sys
import re
import os.path
import argparse
import yaml
from copy import copy
import pickle
from gsd import hoomd as gsd
import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
import random
import glob
#import seaborn as sns
#from numpy.lib.stride_tricks import sliding_window_view
#sns.set_context("talk", font_scale=1.0)


def round_to_int(no):
    result=[]
    nint=int(no)
    if(np.abs(nint-no)<0.5):
        res=int(np.floor(no))
    else:
        res=int(np.ceil(no))
    return res

def sample_spherical(npoints, ndim=3):
#simple points on a sphere from here https://stackoverflow.com/questions/33976911/generate-a-random-sample-of-points-distributed-on-the-surface-of-a-unit-sphere
    vec = np.random.randn(ndim, npoints)
    vec /= np.linalg.norm(vec, axis=0)
    return vec.reshape((-1,3))

def sphere_fibonacci_grid_points (ng,r):
# Output, real xgB(3,ng): the grid points.
  phi = ( 1.0 + np.sqrt ( 5.0 ) ) / 2.0
  theta = np.zeros ( ng )
  sphi = np.zeros ( ng )
  cphi = np.zeros ( ng )
  for i in range ( 0, ng ):
    i2 = 2 * i - ( ng - 1 )
    theta[i] = 2.0 * np.pi * float ( i2 ) / phi
    sphi[i] = float ( i2 ) / float ( ng )
    cphi[i] = np.sqrt ( float ( ng + i2 ) * float ( ng - i2 ) ) / float ( ng )
  xg = np.zeros ( ( ng, 3 ) )
  for i in range ( 0, ng ) :
    xg[i,0] = cphi[i] * np.sin ( theta[i] ) * r
    xg[i,1] = cphi[i] * np.cos ( theta[i] ) * r
    xg[i,2] = sphi[i] * r
  return xg

def read_xyz_file(seed_file, particle_type_list):
    fh = open(seed_file,'r')
    line = fh.readline()
    config = []
    type_dict = {}
    type_list = []
    n_atoms = int(line)
    #skip a blank line
    line = fh.readline()
    line = fh.readline()
    while line:
        type, x, y, z = line.split()
        #if not type in type_dict:
        #    type_dict[type] = len(type_dict)
        #type_list.append(type_dict[type])
        type_list.append(particle_type_list.index(type))
        config.append(np.array((x,y,z)).astype(float))
        line = fh.readline()
    config = np.array(config)
    assert len(config) == n_atoms, "xyz reader- Number of atoms doesn't match the number of lines read in"
    return config, type_list

def screened_potential_shifted_force(r,rmin,rmax,radius_sum,steric_prefactor,electrostatic_prefactor,H,d):

    #d=debye length,H=2*brush_length 

    #a and b are set so that the force and energy goes to zero at rmax
    #make assumption that rmax > H

    #separate distances for electrostatic and hard sphere

    l = r - radius_sum
    lcut = rmax - radius_sum

    a = electrostatic_prefactor/d*np.exp(-lcut/d)
    b = -electrostatic_prefactor*np.exp(-lcut/d) - lcut*a

    p_term1 = electrostatic_prefactor*np.exp(-l/d) + a*l + b
    p_term2 = steric_prefactor*(28*((H/l)**.25-1) + (20./11)*(1-(l/H)**2.75)+ 12*(l/H-1))

    if type(p_term2) == np.ndarray: p_term2[np.isnan(p_term2)]=0

    potential = p_term1 + p_term2*(l<H)

    # F = -dU/dr
    f_term1 = electrostatic_prefactor/d*np.exp(-l/d) - a
    f_term2 = -steric_prefactor/H*(12-7*(H/l)**(1.25)-5*(l/H)**(1.75))
    if type(f_term2) == np.ndarray: f_term2[np.isnan(f_term2)]=0
    force = f_term1 + (f_term2)*(l<H)

    return potential, force

def flat_list(alist): # [[1,2],[3,4]] --> [1,2,3,4]
    flatlist=list()
    for sublist in alist:
        for i in sublist:
            flatlist.append(i)
    return(flatlist)

def remove_duplicates(alist):
    res = []
    for i in alist:
        if i not in res:
            res.append(i)

    return res

