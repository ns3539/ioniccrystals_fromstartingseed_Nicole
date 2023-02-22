from readinputs import *

mV_to_kBT = 25.7
joule_to_kBT = 4.11e-21

surface_potentials=[float(i) for i in surface_potentials.split(',')]
surface_potentials[0] = surface_potentials[0]*-1
print("Charges (in mV) are:")
print(surface_potentials)

num_particle_types=len(surface_potentials)

radii_particles=[float(i) for i in radii.split(',')]
print(radii_particles)

#fractions=[float(i) for i in fractions.split(',')]

brush_lengths=[float(i) for i in brush_lengths.split(',')]

#assert len(fractions)==2,"fraction list should have a size of 2, one value for the number of rigid bodies in the lattice and the other for the free positive particles"
assert len(radii_particles)==2, "2 radius values for the 2 types of particles"
assert len(brush_lengths)==2, "2 brush length values for the 2 types of particles"
assert len(surface_potentials)==2, "2 surface potential values for the 2 types of particles"
assert surface_potentials[0]<0, "First particle in the list should be the negative particle"
assert surface_potentials[1]>0, "Next particle should be the positive particle"

num_negative_types=0
num_positive_types=0

#particle types are put in a list below

for i in range(len(surface_potentials)):
    value=surface_potentials[i]
    if(value<0):
        num_negative_types+=1
    else:
        num_positive_types+=1

base_types=['P','N']
particle_types=[]

for j in range(num_negative_types):
    if(num_negative_types==1):
        particle_types.append(base_types[1])
    else:
        particle_types.append(base_types[1]+str(j+1))
for k in range(num_positive_types):
    if(num_positive_types==1):
        particle_types.append(base_types[0])
    else:
        particle_types.append(base_types[0]+str(k+1))

if(num_positive_types==0 or num_negative_types==0):
    num_ionic_radii=0
else:
    num_ionic_radii=num_positive_types*num_negative_types


print("Particle types:",particle_types)

ionic_radii={}


radius_dict={}

if(num_ionic_radii!=0):
    for i in range(len(particle_types)):
        for j in range(i,len(particle_types)):
            particlei=particle_types[i]
            particlej=particle_types[j]
            pair=(particlei,particlej)
            key=pair
                #use harmonic average from Derjaguin approximation: https://en.wikipedia.org/wiki/Derjaguin_approximation
            er= np.round(2./(1./radii_particles[i] + 1./radii_particles[j]),3)
            radius_dict[key]=er
            if(particlei[0]!=particlej[0]):
                ionic_radii[key]=er

max_radius = np.max(radii_particles)
radius_list=radii_particles + list(ionic_radii.values())
print("Radius of N, P and ionic radii(s): ",radius_list)

permitivity = 8.85e-12 #Farad/M
electrostatic_constant=(2*np.pi*dielectric_constant*permitivity/joule_to_kBT)

surface_potentials=np.array(surface_potentials)
surface_potentials_in_V=surface_potentials/1000

electrostatic_prefactors={}
steric_prefactors={}

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

brush2=sum_brush_length
radiusN=radii_particles[0]
radiusP=radii_particles[1]

print("radiusN:",radiusN)
print("radiusP:",radiusP)

print("Electrostatic prefactors are: ",electrostatic_prefactors)
print("Steric prefactors are: ",steric_prefactors)

outputprefix="repeats"+str(lattice_repeats)+"_a"+str(lattice_spacing)+"_R"+",".join([str(ele) for ele in radii_particles])+"_fracpos"+str(fraction_positive)+"_debye"+str(debye_length)+"_brushL"+",".join([str(ele) for ele in brush_lengths])+"_charges"+",".join([str(ele) for ele in surface_potentials])+"_seed"+str(seed)

