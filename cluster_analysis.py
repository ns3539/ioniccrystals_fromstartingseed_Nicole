from readinputs import *

from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist

import numpy as np
import gsd.hoomd
import argparse
import os
import sys

#functions for analysis to obtain cluster sizes 
def distance(r,box_size):
    r=r-box_size*np.floor(r/box_size+0.5)
    return(np.sqrt(r[0]**2+r[1]**2+r[2]**2))

def split_into_clusters(link_mat,thresh,n):
   c_ts=n
   clusters={}
   for row in link_mat:
      if row[2] < thresh:
          n_1=int(row[0])
          n_2=int(row[1])

          if n_1 >= n:
             link_1=clusters[n_1]
             del(clusters[n_1])
          else:
             link_1= [n_1]

          if n_2 >= n:
             link_2=clusters[n_2]
             del(clusters[n_2])
          else:
             link_2= [n_2]

          link_1.extend(link_2)
          clusters[c_ts] = link_1
          c_ts+=1
      else:
          return clusters


#get gsd filename for trajectory
inputprefix="repeats"+str(lattice_repeats)+"_a"+str(lattice_spacing)+"_R"+",".join([str(ele) for ele in radii.split(',')])+"_fracpos"+str(fraction_positive)+"_debye"+str(debye_length)+"_brushL"+",".join([str(ele) for ele in brush_lengths.split(',')])+"_charges-"+",".join([str(ele) for ele in surface_potentials.split(',')])+"_seed"+str(seed)
progressfile = inputprefix+'.progress.txt'
print(progressfile)
if os.path.exists(progressfile):
    initial_steps= int(open(progressfile,'r').readlines()[0])
    final_steps = int(open(progressfile,'r').readlines()[-1])
    
    initial_filename = inputprefix+'.run.%i.gsd'%initial_steps
    final_filename = inputprefix+'.run.%i.gsd'%final_steps
    
    print('file for initial cluster data:', initial_filename)
    print('file for final cluster data:', final_filename)
    
    #analysis for initial seed
    traj = gsd.hoomd.open(initial_filename, "rb")
    numsnap = len(traj)
    gsdfile=os.path.basename(initial_filename)
    
    #firstframe=0
    #lastframe=numsnap-1
    
    box=traj[0].configuration.box[:3]
    
    cutoff=230 #choose the cutoff according to the size of the particles forming the cluster
    #clusters_save=[]
    
    # get largest cluster in the first frame of initial file 
    snap=traj[0]
    timestep=snap.configuration.step
    #timestep_list.append(timestep)
    points = snap.particles.position
    pairwise_distances=pdist(points,metric='euclidean')
    z = linkage(pairwise_distances, method='single',metric='euclidean')
    N=snap.particles.position.shape[0]
    clustering = split_into_clusters(z,cutoff,N)
    cluster_num = 0
    cluster_sizes=[]
    for cluster in clustering:
        cluster_num+=1
        cluster_sizes.append(len(clustering[cluster]))
    #noofclusters_allframes.append(len(cluster_sizes))
    #index=cluster_sizes.index(maxclustersize)
    if clustering==None:
        print("No clusters in first frame! Most likely your distance cutoff is too high")
    #largestclustersize_allframes.append(maxclustersize)
    else:
        maxclustersize_first=np.max(cluster_sizes)
        print('size of largest cluster in first frame:', maxclustersize_first)
        
        #get largest cluster in the last frame at end of simulation
        
        traj = gsd.hoomd.open(final_filename, "rb")
        numsnap = len(traj)
        gsdfile=os.path.basename(final_filename)
    
        #firstframe=0
        lastframe=numsnap-1
        
        box=traj[0].configuration.box[:3]
        snap=traj[int(lastframe)]
        timestep=snap.configuration.step
        #timestep_list.append(timestep)
        points = snap.particles.position
        pairwise_distances=pdist(points,metric='euclidean')
        z = linkage(pairwise_distances, method='single',metric='euclidean')
        N=snap.particles.position.shape[0]
        clustering = split_into_clusters(z,cutoff,N)
        cluster_num = 0
        cluster_sizes=[]
        for cluster in clustering:
            cluster_num+=1
            cluster_sizes.append(len(clustering[cluster]))
        maxclustersize_last=np.max(cluster_sizes)
        #noofclusters_allframes.append(len(cluster_sizes))
        #index=cluster_sizes.index(maxclustersize)
        if clustering==None:
            print("No clusters in last frame! Most likely your distance cutoff is too high")
        #largestclustersize_allframes.append(maxclustersize)
        else:
            print('size of largest cluster in last frame:', maxclustersize_last)
            
            #clusters_save.append(maxclustersize_last_max/max_clustersize_first)
            relative_size = maxclustersize_last/maxclustersize_first
            print("ratio of cluster sizes:", relative_size)
        
        
    #Save output file containing analysis data in the same directory as simulation trajectory
    #output_file_clusteranalysis=os.path.splitext(filename)[0]+'.clusteranalysis.txt'
    output_file_clusteranalysis = inputprefix + '.clusteranalysis.txt'
    #d=np.array([clusters_save],dtype=object)
    #print(d)
    #d.dumps(output_file_clusteranalysis)
    
    with open(output_file_clusteranalysis, 'w') as f:
        #f.write(str(clusters_save))
        f.write(str(relative_size))


else:
    print("file path does not exist")
