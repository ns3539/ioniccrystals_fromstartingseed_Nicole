from add_interactions import *

if temperature.find(':')>0:
    temp_list = np.array( temperature.split(':') ,dtype=float)
    assert len(temp_list)%2==0 and len(temp_list)>=0,"must give an even number of temperature arguments"
    temp_list = [tuple(x) for x in temp_list.reshape((-1,2))]
    temperature = hoomd.variant.linear_interp(points=temp_list)
else:
    temperature = float(temperature)

all = hoomd.group.all()


hoomd.analyze.log(filename=outputlogfile,
                  quantities=['potential_energy', 'temperature'],
                  period=dump_frequency//10,
                  overwrite=True)

if mode.lower() == "minimize":
    fire = hoomd.md.integrate.mode_minimize_fire(dt=dt,Etol=1e-7,min_steps=1000)
    nve=hoomd.md.integrate.nve(group=hoomd.group.all())
    hoomd.dump.gsd(outputprefix+'_minimize.seedcrystal.gsd', period=500, group=all, overwrite=True,dynamic=['attribute','momentum','topology'])
    while not(fire.has_converged()):
        hoomd.run(1000000)
    sys.exit()

elif mode.lower() == "langevin":
    print("mode:",mode)
    hoomd.md.integrate.mode_standard(dt=dt)
    ld=hoomd.md.integrate.langevin(group=all, kT=temperature, seed=seed)
    print("Using gamma value as:", gamma)
    for particle_type in particle_type_list:
        ld.set_gamma(particle_type,gamma)


elif mode.upper() == "NVT":
    hoomd.md.integrate.mode_standard(dt=dt)
    integrator = hoomd.md.integrate.nvt(group=all, kT=temperature, tau=1/gamma)
    integrator.randomize_velocities(seed=seed)

elif mode.upper() == "NPT":
    hoomd.md.integrate.mode_standard(dt=dt)
    integrator1 = hoomd.md.integrate.nvt(group=all, kT=temperature, tau=1/gamma)
    eq = hoomd.dump.gsd(outputprefix+'.eq.gsd', period=dump_frequency, group=all, overwrite=True)
    integrator1.randomize_velocities(seed=seed)
    hoomd.run(100000)
    eq.disable()
    integrator1.disable()
    integrator = hoomd.md.integrate.npt(group=all, kT=1.0, tau=100*dt, tauP=100*dt, P=1e-8)
    #integrator.randomize_velocities(seed=seed)

else:
    print("Mode '%s' not supported"%mode)
    sys.exit(1)


hoomd.dump.gsd(outputgsdfile, period=dump_frequency, group=all, overwrite=True,dynamic=['attribute','momentum','topology'])
hoomd.run(nsteps)

#very end of simulation, assuming successful end
fh = open(progressfile,'a')
fh.write("%i\n"%totalsteps)
fh.close()



