from define_quantities import *

################################## Set up the simulation to be restartable if progressfile is found ##################################


progressfile = outputprefix+'.progress.txt'
if os.path.exists(progressfile):
    continuesim=True
    prevsteps = int(open(progressfile,'r').readlines()[-1])
    inputfile = outputprefix+'.run.%i.gsd'%prevsteps
    totalsteps = prevsteps + nsteps
else:
    continuesim = False
    prevsteps = 0
    totalsteps = prevsteps + nsteps

outputgsdfile = outputprefix+'.run.%i.gsd'%totalsteps
outputlogfile = outputprefix+'.run.%i.log'%totalsteps

######################################################### Snapshot creation ##########################################################


if(continuesim==True):
    system=hoomd.init.read_gsd(inputfile,frame=-1)
    snapshot = system.take_snapshot()
    particle_type_list = system.particles.types
    print("Particle types:",particle_type_list)
    #half of the box lengths in x, y, z directions
    max_z = snapshot.box.Lz/2.
    max_x = snapshot.box.Lx/2.
    max_y = snapshot.box.Ly/2.
    print("Box size:")
    print(system.box)

    #NEVER EXPAND THE BOX ON RESTART. BOX TO BE EXPANDED ONLY FOR FIRST RUN 

    mod_maxX = system.box.Lx/2.
    mod_maxY = system.box.Ly/2.
    mod_maxZ = system.box.Lz/2.

else:
    print("Lattice_type: ",lattice_type)
    if lattice_type == "fcc":
        lattice = hoomd.lattice.fcc
    elif lattice_type == "bcc":
        lattice = hoomd.lattice.bcc
    elif lattice_type == "sc":
        lattice = hoomd.lattice.sc
    else:
        print("Lattice type %s is not supported (only fcc, bcc, sc)"%lattice_type)
        sys.exit(1)


    if(radiusN>radiusP):
        latticesite='N'
        satellite='P'
        radius_latticesite=radiusN
        radius_satellite=radiusP
    else:
        latticesite='P'
        satellite='N'
        radius_latticesite=radiusP
        radius_satellite=radiusN
    

    print("Lattice site is:",latticesite)

    orbit_distance = (radiusP+radiusN+brush2)*orbit_factor
    print("Putting satellite particles at a distance of %f from center particles"%orbit_distance)
    num_central = lattice_repeats**3
    num_satellite_per_center = int((1-fraction_positive)/fraction_positive)
    print("num_satellite_per_center:",num_satellite_per_center)


    center_unit_cell = lattice(a=lattice_spacing*2*radius_latticesite)
    center_in_uc = center_unit_cell.N
    satellite_in_uc = center_in_uc*num_satellite_per_center
    print("center_in_uc:",center_in_uc)
    print("satellite_in_uc:",satellite_in_uc)

    total_in_uc = satellite_in_uc + center_in_uc

    print("Total number of particles in each unit cell:",total_in_uc)
    print("Total lattice spacing:",lattice_spacing*2*radius_latticesite)

    #scale mass if not set
    if massN is None:
        massN = massP/(radiusP/radiusN)**3
        #massN=massP
        print("Setting mass N: %f (mass P: %f)"%(massN,massP))

    mass_dict = {}
    for i in range(len(particle_types)):
        if(particle_types[i][0]=='P'):
            mass_dict[particle_types[i]]=massP
        elif(particle_types[i][0]=='N'):
            mass_dict[particle_types[i]]=massN


    particle_positions = center_unit_cell.position

    all_positions = []
    all_masses = []
    all_types = []
    all_diameters = []

    for i in range(center_in_uc):
        xyz = particle_positions[i]
        all_positions.append(xyz)
        all_masses.append(mass_dict[latticesite])
        all_types.append(latticesite)
        all_diameters.append(2*radius_latticesite)

        #sattelite_positions = sample_spherical(num_satellite_per_center)*orbit_distance
        sattelite_positions = sphere_fibonacci_grid_points (num_satellite_per_center,orbit_distance)
        for xyz_s in sattelite_positions:
            all_positions.append(xyz_s+xyz)
            all_masses.append(massN)
            all_types.append(satellite)
            all_diameters.append(2*radius_satellite)


    uc = hoomd.lattice.unitcell(N = total_in_uc,
                    a1 = center_unit_cell.a1,
                    a2 = center_unit_cell.a2,
                    a3 = center_unit_cell.a3,
                    dimensions = center_unit_cell.dimensions,
                    position = all_positions,
                    type_name = all_types,
                    mass = all_masses,
                    diameter = np.array(all_diameters),)

    snapshot = uc.get_snapshot()
    snapshot.replicate(lattice_repeats,lattice_repeats,lattice_repeats)
    snapshot.particles.types=[latticesite,satellite]

    print(snapshot.particles.types)

    N=len(snapshot.particles.typeid)
    
    print("Total number of lattice sites:",N)

    particle_type_list=snapshot.particles.types


    print(snapshot.particles.position)
        
    #add seed
    if seed_file is not None:
        if os.path.exists(seed_file):
            seed_xyz, seed_particle_types = read_xyz_file(seed_file,particle_types)
            # divided by max diameter, have to scale back up
            if scale_by == "max_diameter":
                seed_positions = seed_xyz * 2*max_radius
                print("Warning: scaling xyz by max diameter")
            elif scale_by == "max_radius":
                seed_positions = seed_xyz * max_radius
                print("Warning: scaling xyz by max radius")
            else:
                print("Warning: not scaling xyz by anything")
                seed_positions = seed_xyz
        else:
            print("Path does not exist")
            seed_positions = None
    else:
        print("Seed file is None")
        seed_positions = None

    #half of the box lengths in x, y, z directions
    max_z = snapshot.box.Lz/2.
    max_x = snapshot.box.Lx/2.
    max_y = snapshot.box.Ly/2.
    
    print("Seed positions:")
    print(seed_positions)

    system=hoomd.init.read_snapshot(snapshot)

    print("Box size:")
    print(system.box)


    print("max_x,max_y,max_z")
    print(max_x,max_y,max_z)


    if(closedwalls=='True'):
        hoomd.update.box_resize(Lx = 2.*max_x+8*max_radius, Ly = 2.*max_y+8*max_radius, Lz=2.*max_z+8*max_radius, period=None,scale_particles=False)
        print("Updated box size in presence of closed walls:")
        print(system.box)
    
    mod_maxX = system.box.Lx/2.
    mod_maxY = system.box.Ly/2.
    mod_maxZ = system.box.Lz/2.

    
    """
    ptypes=[]
    typeids=[]

    for i in range(len(system.particles)):
        ptypes.append(system.particles[i].type)
        typeids.append(system.particles[i].typeid)

    particle_type_list = system.particles.types
    print(particle_type_list)
    """
    #print(seed_particle_types)

    print("Number of particles in system before adding any seed crystal: ",len(system.particles))

    print("Removing lattice particles that overlap with particles of the seed crystal.........")
    if seed_positions is not None:
        first_idx = 0
        for pidx, pxyz in enumerate(seed_positions):
            particle_type_id = seed_particle_types[pidx]
            particle_type = particle_types[particle_type_id]
            t = system.particles.add(particle_type)
            if first_idx ==0: first_idx = int(t)
            system.particles[t].position = pxyz
            system.particles[t].diameter = radii_particles[particle_type_id]*2
            #print(particle_type,particle_type_id,system.particles[t].diameter)
            system.particles[t].mass = mass_dict[particle_type]
            
        print("Seed particles assigned their attributes.......")
        print("Number of particles in system after adding the seed crystal: ",len(system.particles))

        # search for overlaps
        tags_to_remove = []
         
        box_size=np.array([system.box.Lx,system.box.Ly,system.box.Lz])
        
        print("Start time:")
        start = time.process_time()
        print(start)
        
        for pidx, pxyz in enumerate(seed_positions):
            particle_type_id = seed_particle_types[pidx]
            radius=radius_list[particle_type_id]

            systemtypeids=np.array([p.typeid for p in system.particles])
            systempositions=np.array([p.position for p in system.particles])
            systemradii=np.array([0.5*(p.diameter) for p in system.particles])
            systemtags=np.array([p.tag for p in system.particles])
            
            cutoff=np.ones(systemtypeids.shape[0])
            cutoff*=((systemtypeids==particle_type_id)*(2*systemradii*1.3)**2)+((systemtypeids!=particle_type_id)*(((radius+systemradii)*1.3)**2))
            #print(cutoff,"     ",systemradii,"     ",radius)             
            dr=pxyz-systempositions
            #dr=system.box.min_image(dr)
            dr = dr - box_size*np.floor(dr/box_size+0.5)
            dr2_mag=np.array([np.dot(k,k) for k in dr])

            condition1=(dr2_mag<cutoff)
            condition2=(systemtags<first_idx)
            netcondition=(condition1*condition2)
            #print("Seed particle ",pidx," No of system particles: ",netcondition.shape[0])

            locations=np.array(np.where(netcondition==True)[0])
            tags=[systemtags[i] for i in locations if len(locations)!=0]
            if(len(tags)!=0):
                tags_to_remove.append(tags)

        """
        
        for pidx, pxyz in enumerate(seed_positions):
            particle_type_id = seed_particle_types[pidx]
            for p in system.particles:
                pid = p.typeid
                if particle_type_id == pid:
                    cutoff = (2*radius_list[pid]*1.3)**2
                else:
                    cutoff = ((radius_list[particle_type_id]+radius_list[pid])*1.3)**2
                dr = pxyz - p.position
                dr = dr - box_size*np.floor(dr/box_size+0.5)
                #dr = system.box.min_image(dr)
                dr2_mag = np.dot(dr,dr)
                if dr2_mag < cutoff and p.tag<first_idx: tags_to_remove.append(p.tag)
        """

        
        print("Time taken for this seed overlap removal part to execute: ",time.process_time() - start)
        
        tags_to_remove_1=flat_list(tags_to_remove)
        tags_to_remove=remove_duplicates(tags_to_remove_1)
        
        print("**************************************************************************************************")
        print("Removing %i particles that overlap with seed"%len(tags_to_remove))
        print(tags_to_remove)
        print("**************************************************************************************************")
        for t in set(tags_to_remove):
            system.particles.remove(t)
    
        print("No of particles in system after removing overlaps:",len(system.particles))
        
