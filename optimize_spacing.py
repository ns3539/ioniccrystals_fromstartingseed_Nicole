#!/usr/bin/python
import argparse
import numpy as np
import os
import sys

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
    
    return potential,force

parser=argparse.ArgumentParser()
parser.add_argument("--radiusN",default=65.0,help="Radius type 1 (default: %(default)s)",type=float)
parser.add_argument("--radiusP",default=150.0,help="Radius type 2 (default: %(default)s)",type=float)
parser.add_argument("-B","--brush_length",help="Brush length (default: %(default)s)",default=10.0,type=float)
parser.add_argument("--fudge",default=0,help="Multiple of brush length to add to the spacing (default:0)",type=float)
parser.add_argument("--debye_length",default=5.25,help="Debye length in nm",type=float)
parser.add_argument("-r","--repeats",help="Number of repeats in each direction (default: %(default)s)",default=5,type=int)
parser.add_argument("--outdir",default="NaCl")
args = parser.parse_args()
locals().update(vars(args))


#approx_min = 2*radiusN+2*radiusP+4*brush_length
approx_min = radiusN+radiusP+2*brush_length
large_size = 1.0
max_radius = max(radiusN,radiusP)
small_size = min(radiusN,radiusP)/max_radius

bl_reduced= brush_length/max_radius
#lattice_const = (approx_min + fudge*brush_length) / max_radius 

#lattice_const=((2*large_size + 2*bl_reduced) + (fudge/max_radius))

lattice_const=approx_min/max_radius

if not lattice_const*np.sqrt(2)/2. > 2*large_size + 2*bl_reduced:
    print("Warning: lattice_const*sqrt(2)/2 is too small compared to 2R_max + 2*brush",lattice_const*np.sqrt(2)/2.,2*(large_size+bl_reduced))
else:
    print("lattice_const*sqrt(2)/2 is good compared to 2R_max+2*brush",lattice_const*np.sqrt(2)/2.,2*(large_size+bl_reduced))

particle_types=['N','P']
brush_lengths=[brush_length,brush_length]
radii_particles=[radiusN,radiusP]
surface_potentials=[50.0,50.0]
dielectric_constant=80.0
mV_to_kBT = 25.7
joule_to_kBT = 4.11e-21
permitivity = 8.85e-12 #Farad/M
brush_density=0.09
#debye_length=5.5

electrostatic_prefactors={}
steric_prefactors={}

electrostatic_constant=(2*np.pi*dielectric_constant*permitivity/joule_to_kBT)
surface_potentials=np.array(surface_potentials)
surface_potentials[0] = surface_potentials[0]*-1

print(surface_potentials)
surface_potentials_in_V=surface_potentials/1000

for i in range(len(particle_types)):
    for j in range(i,len(particle_types)):
        particlei=particle_types[i]
        particlej=particle_types[j]
        pair=(particlei,particlej)
        key=pair
        effectiveradius=np.round(2./(1./radii_particles[i] + 1./radii_particles[j]),3)
        effectiveradius_in_m=1e-9*effectiveradius
        electrostatic_prefactors[key]=electrostatic_constant*effectiveradius_in_m*surface_potentials_in_V[i]*surface_potentials_in_V[j]

        repulsionradius=(radii_particles[i]+ radii_particles[j])/2.
        sum_brush_length=brush_lengths[i]+brush_lengths[j]
        steric_prefactors[key]=16*np.pi*repulsionradius*(sum_brush_length**2)*(brush_density**(3./2))/35


rate=0.00001
precision = 0.000001 #This tells us when to stop the algorithm
previous_step_size = 1 #
max_iters = 10000000 # maximum number of iterations
minima_list=[]

for i in range(len(particle_types)):
    typei = particle_types[i]

    for j in range(i,len(particle_types)):
        typej = particle_types[j]

        key=(typei,typej)

        if(typei==typej):
            radius_sum = 2*radii_particles[particle_types.index(typei)]
        else:
            radius_sum = radii_particles[particle_types.index(typei)]+radii_particles[particle_types.index(typej)]

        sum_brush_lengths=brush_lengths[particle_types.index(typei)] + brush_lengths[particle_types.index(typej)]

        #print("*******************************************************************************************************************************************************************************")
        #print("potential for pair "+str(key)+" with radius_sum "+str(radius_sum))
        #potential_range = radius_sum + 20*debye_length

        pot_min=1.00005*radius_sum
        pot_max = radius_sum+20*debye_length
        dpot = (pot_max-pot_min)/5000.
        test_range = np.arange(pot_min,pot_max+dpot,dpot)
        if(typei!=typej):
            key=(typei,typej)
            pot,force = screened_potential_shifted_force(test_range,rmin=pot_min, rmax=pot_max, radius_sum=radius_sum, steric_prefactor=steric_prefactors[key],electrostatic_prefactor=electrostatic_prefactors[key],H=sum_brush_lengths,d=debye_length)
            minima_list.append((i,j,key,np.min(pot),test_range[np.argmin(pot)]))

            cur_x = pot_min # The algorithm starts at x=pot_min
            iters = 0 #iteration counter
            while previous_step_size > precision and iters < max_iters:
                print("***********************************************************************************************")
                prev_x = cur_x #Store current x value in prev_x
                prev_pot,prev_force=screened_potential_shifted_force(prev_x,rmin=pot_min, rmax=pot_max, radius_sum=radius_sum, steric_prefactor=steric_prefactors[key],electrostatic_prefactor=electrostatic_prefactors[key],H=sum_brush_lengths,d=debye_length)
                cur_x = cur_x - rate * (-prev_force) #Grad descent
                previous_step_size = abs(cur_x - prev_x) #Change in x
                iters = iters+1 #iteration count
                print("Iteration:",iters,"\n Value of spacing is:",cur_x) #Print iterations
                cur_pot,cur_force=screened_potential_shifted_force(cur_x,rmin=pot_min, rmax=pot_max, radius_sum=radius_sum, steric_prefactor=steric_prefactors[key],electrostatic_prefactor=electrostatic_prefactors[key],H=sum_brush_lengths,d=debye_length)
                print("Energy at this value of spacing:",cur_pot)
    
            print("The local minimum occurs at", cur_x)

print("i,j,key,potential_min, potential_min_location")
for min_info in minima_list:
    print(min_info)


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
    atoms.write(os.path.join(outdir,"NaCl_rN%.0f_rP%.0f_r%i_B%.0f_fudge%.1f.xyz"%(radiusN,radiusP,repeats,brush_length,fudge)))

