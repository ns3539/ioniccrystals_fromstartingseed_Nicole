from snapshot_creation import *

##Creating neighborlist

#nl = hoomd.md.nlist.cell(r_buff=0.8*ionic_radius)
nl = hoomd.md.nlist.cell()
nl.reset_exclusions(exclusions = ['1-2', 'body','constraint'])

##Implementing gravity

if gravity>0:
    typeN = hoomd.group.type(type='N')
    typeP = hoomd.group.type(type='P')
    fgrav = -gravity
    gravity_force_P = hoomd.md.force.constant(fvec=[0,0,fgrav],group=typeP)
    #gravity_force_N = hoomd.md.force.constant(fvec=[0,0,fgrav*massN],group=N)

##Implementing closed walls on the 6 faces of the simulation box to prevent particles going out of the box

if closedwalls=='True':
    print("Turning on closed walls along x/y/z with max_z:",mod_maxZ)

    #create 6 repulsive walls
    upper_wall_x = hoomd.md.wall.plane(origin=(mod_maxX,0,0),normal=(-1,0,0),inside=True)
    lower_wall_x = hoomd.md.wall.plane(origin=(-mod_maxX,0,0),normal=(1,0,0),inside=True)
    upper_wall_y = hoomd.md.wall.plane(origin=(0,mod_maxY,0),normal=(0,-1,0),inside=True)
    lower_wall_y = hoomd.md.wall.plane(origin=(0,-mod_maxY,0),normal=(0,1,0),inside=True)
    upper_wall_z = hoomd.md.wall.plane(origin=(0,0,mod_maxZ),normal=(0,0,-1),inside=True)
    lower_wall_z = hoomd.md.wall.plane(origin=(0,0,-mod_maxZ),normal=(0,0,1),inside=True)

    wall_group = hoomd.md.wall.group(upper_wall_x,lower_wall_x,upper_wall_y,lower_wall_y,upper_wall_z,lower_wall_z)
    wall_force = hoomd.md.wall.slj(wall_group,r_cut=max_radius*(2.0**(1.0/6.0)))
    wall_force.force_coeff.set('N',epsilon=1,sigma=radiusN,alpha=0,r_cut=radiusN*(2.0**(1.0/6.0)))
    wall_force.force_coeff.set('P',epsilon=1,sigma=radiusP,alpha=0,r_cut=radiusP*(2.0**(1.0/6.0)))

##Implementing a charged attractive wall at the bottom of the simulation box

# if chargewall!=0 and not closedwalls:
if chargewall!=0:
    eps=np.abs(chargewall)
    z_chgwall = -mod_maxZ + 0.5*max_radius
    #upper_wall = hoomd.md.wall.plane(origin=(0,0,0),normal=(0,0,-1),inside=True)
    lower_wall = hoomd.md.wall.plane(origin=(0,0,z_chgwall),normal=(0,0,1),inside=True)
    wall_group = hoomd.md.wall.group(lower_wall)

    if chargewall < 0:
        attractive_wall_force = hoomd.md.wall.lj(wall_group,r_cut=radiusP*5)
        attractive_wall_force.force_coeff.set('N',epsilon=0,sigma=radiusN*2)
        attractive_wall_force.force_coeff.set('P',epsilon=eps,sigma=radiusP*2)

    if chargewall > 0:
        attractive_wall_force = hoomd.md.wall.lj(wall_group,r_cut=radiusN*5)
        attractive_wall_force.force_coeff.set('N',epsilon=eps,sigma=radiusN*2)
        attractive_wall_force.force_coeff.set('P',epsilon=0,sigma=radiusP*2)


table = hoomd.md.pair.table(width=5000,nlist=nl)
minima_list=[]

for i in range(len(particle_types)):
    typei = particle_types[i]

    for j in range(i,len(particle_types)):
        typej = particle_types[j]

        key=(typei,typej)
        my_radius = radius_dict[key]

        if(typei==typej):
            radius_sum = 2*radius_dict[key]
        else:
            radius_sum = radii_particles[particle_types.index(typei)]+radii_particles[particle_types.index(typej)]

        sum_brush_lengths=brush_lengths[particle_types.index(typei)] + brush_lengths[particle_types.index(typej)]

        print("*******************************************************************************************************************************************************************************")
        print("Setting potential for pair "+str(key)+" with radius_sum "+str(radius_sum))

        potential_range = radius_sum + 20*debye_length
        table.pair_coeff.set(typei,typej, func=screened_potential_shifted_force,rmin=1.00005*radius_sum, rmax=potential_range,coeff=dict(radius_sum=radius_sum, steric_prefactor=steric_prefactors[key],electrostatic_prefactor=electrostatic_prefactors[key],H=sum_brush_lengths,d=debye_length),)

        pot_min=1.00005*radius_sum
        pot_max = radius_sum+20*debye_length
        dpot = (pot_max-pot_min)/5000.
        test_range = np.arange(pot_min,pot_max+dpot,dpot)
        pot,force = screened_potential_shifted_force(test_range,rmin=pot_min, rmax=pot_max, radius_sum=radius_sum, steric_prefactor=steric_prefactors[key],electrostatic_prefactor=electrostatic_prefactors[key],H=sum_brush_lengths,d=debye_length)
        np.savetxt(outputprefix+".pairpotential_%i,%i.txt"%(i,j), np.concatenate( (test_range.reshape(-1,1),pot.reshape(-1,1),force.reshape(-1,1)),axis=-1),header="Distance                   Potential                   Force" )
        minima_list.append((i,j,key,np.min(pot),test_range[np.argmin(pot)]))

print("i,j,key,potential_min, potential_min_location")
for min_info in minima_list:
    print(min_info)


