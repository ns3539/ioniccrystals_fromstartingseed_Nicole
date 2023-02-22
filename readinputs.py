from functions import *
"""
parser=argparse.ArgumentParser()
parser.add_argument("--gpu",default=False,action="store_true")
parser.add_argument("--closedwalls",default=False,help="Turn on walls perp to x/y/z (default: False)",action="store_true")
parser.add_argument("--gravity",default=0.5,help="Turn on downward gravity constant force in z by setting to a positive number (default: 0)",type=float)
parser.add_argument("--chargewall",default=0,help="Turn on wall perp to z attracting -(>0) or +(<0) particles (default: False)",type=float)
parser.add_argument("-R","--radii",help="Radii of the different types of particles as a string",type=str)
parser.add_argument("-B","--brush_lengths",help="Brush lengths as a string",type=str)
parser.add_argument("-s","--brush_density",help="Brush density (default: %(default)s)",default=0.09,type=float)
parser.add_argument("-d","--debye_length",help="Debye length (default: %(default)s)",default=4,type=float)
parser.add_argument("--surface_potentials",help="Surface potentials in mV as a string",type=str)
parser.add_argument("--dielectric_constant",help="Dielectric constant (default: %(default)s)",default=80,type=float)
parser.add_argument("--gamma",help="Drag coefficient (default: %(default)s)",default=0.01,type=float)
parser.add_argument("--massP",help="Mass of a positive particle",default=1.0,type=float)
parser.add_argument("--massN",help="Mass of a negative particle",default=None,type=float)
parser.add_argument("--seed",help="Random seed (default: %(default)s)",default=1,type=int)
parser.add_argument("--fraction_positive",help="Fraction positive charge (default: %(default)s)",default=0.1,type=float)
parser.add_argument("-a","--lattice_spacing",help="Lattice spacing in terms of positive particle diameter (default: %(default)s)",default=1.5,type=float)
parser.add_argument("--lattice_repeats",default=5,type=int,help="times to repliacte the system in each direction (default: %(default)s)") #no of times we want to replicate in one direction
parser.add_argument("--lattice_type",default="bcc",help="Lattice type (bcc, sc, fcc) (default: %(default)s)")
parser.add_argument("--orbit_factor",default=1.3,type=float,help="Factor beyond sum of radii and brush at which to start N particles (default: %(default)s)")
parser.add_argument("--dt",help="Simulation time step (default: %(default)s)",type=float,default=0.1)
parser.add_argument("--seed_file",help="xyz file specifying some seed coordinates",default=None)
parser.add_argument("-i","--inputfile",default=None,help="Input file (gsd file, optional)",type=str)
parser.add_argument("--scale_by",default="max_diameter",help="Scale xyz by this (default: %(default)s), options: max_diameter, max_radius")
parser.add_argument("-o","--outputprefix",help="Output prefix (required)",type=str,required=True)
parser.add_argument("-n","--nsteps",help="Number of steps to run",type=int,required=True)
parser.add_argument("-T","--temperature",help="Temperature at which to run, or list of tuples specifying annealing schedule",default="1.0")
parser.add_argument("-m","--mode",help="Integrator/dynamics scheme. Allowed: Minimize, Langevin, NVT (default: %(default)s)",default="Langevin")
args = parser.parse_args()

locals().update(vars(args))
"""


####################################################### Loading the yaml input file #################################################

with open ('inputs.yaml') as f:
        data = yaml.load(f,Loader=yaml.FullLoader)

####################################################### Reading general info ########################################################

generalinputdict=data # create a dictionary for inputs


gpu=str(generalinputdict['gpu'])
closedwalls=str(generalinputdict['closedwalls'])
gravity=float(generalinputdict['gravity'])
chargewall=float(generalinputdict['chargewall'])
radii=str(generalinputdict['radii'])
brush_lengths=str(generalinputdict['brush_lengths'])
brush_density=float(generalinputdict['brush_density'])
debye_length=float(generalinputdict['debye_length'])
surface_potentials=str(generalinputdict['surface_potentials'])
dielectric_constant=float(generalinputdict['dielectric_constant'])
seed_file=generalinputdict['seed_file']
gen=int(generalinputdict['gen'])
print("generation:", gen)


gamma=float(generalinputdict['gamma'])
massP=float(generalinputdict['massP'])
if(generalinputdict['massN']!=None):
    massN=float(generalinputdict['massN'])
else:
    massN=generalinputdict['massN']
seed=int(generalinputdict['seed'])
dt=float(generalinputdict['dt'])
dumptime=int(generalinputdict['dumptime'])


fraction_positive=float(generalinputdict['fraction_positive'])
lattice_spacing=float(generalinputdict['lattice_spacing'])
lattice_repeats=int(generalinputdict['lattice_repeats'])
lattice_type=str(generalinputdict['lattice_type'])
orbit_factor=float(generalinputdict['orbit_factor'])


scale_by=str(generalinputdict['scale_by'])
outputprefix=str(generalinputdict['outputprefix'])
inputfile=str(generalinputdict['inputfile'])

nsteps=int(generalinputdict['nsteps'])
temperature=str(generalinputdict['temperature'])
mode=str(generalinputdict['mode'])

if gpu=='True':
    hoomd.context.initialize("--mode=gpu")
    print("Running on the GPU")
else:
    hoomd.context.initialize("--mode=cpu")

np.random.seed(seed)
dump_frequency=int(dumptime/dt)


