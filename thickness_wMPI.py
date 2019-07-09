"""
Created on Sun Jul 07 13:38:20 2019

@author: amaan
"""
"""
IMPORT MODULES
"""


from mpi4py import MPI
import numpy as np
import time
import datetime
import json
import cPickle as pickle

import get
from stl import mesh
import matplotlib.pyplot as plt
import csv
import sys

"""
CALCULATE OBSTRUCTIONS
"""
def calc_OBS_para(shells,rukawat,src,elem,elem_r,elem_noruk,elem_noruk_t,thetas,n,i):
    for ii in range(i*n,(i+1)*n):
                    # replace this with elem
                    # get an element p
        p = get.p(shells.x[ii,:],shells.y[ii,:],shells.z[ii,:])
                        
                        # get coordinates of center of the element
        c_p = get.center(p['p1'],p['p2'],p['p3'])
                        
                        # Get Normal to element
        N_p = get.normal(p['p1'],p['p2'],p['p3'])
                        
                        # get ray from source
                        
                        #ray = get.sourcevector(src['NO'],c_p)
        ray = get.sourcevector(src['NO'],c_p)
                        
                        # get angle between ray and element
        m = 0
        n = 0
        ray_angle = get.angle(N_p,ray['d'])
                        # ray angle less than 90 implies source faces the face
        if ray_angle < 90:  
                            
                            # Check for obstruction with evaporator assembly and centerstack-stack
                            # loop through all elements in rukawat.
                    #        print('now checking against all rukawat', i)
            for j in range(0,elem_r):
                                # get an element from rukawat
                rwat = get.p(rukawat.x[j,:],rukawat.y[j,:],rukawat.z[j,:])
                                # get center of rukawat element
                c_rwat = get.center(rwat['p1'],rwat['p2'],rwat['p3'])
                                # get normal to obstruction/rukawat
                N_rwat = get.normal(rwat['p1'],rwat['p3'],rwat['p2'])
                                # Find Q, the point along ray that lies on the plane of rwat
                Q = get.rwatQ(N_rwat,ray,c_rwat)
                                # Check if Q belongs to rwat and if yes is t_Q < t
                OBS = get.obstruction(rwat,Q,ray['t'],N_rwat,ray['d'])          
                if OBS == True:
                                #print("Obstruction found in rukawat",j,ii)
                    break
                else:
                    m = m+1
                            
                                
                            # check for obstruction from shell elements
            if OBS == False:
                    #            print('now checking against all shell elements', i)   
                for k in [u for u in range(0,elem) if u!=ii]:
                    rwat = get.p(shells.x[k,:],shells.y[k,:],shells.z[k,:])
                    c_rwat = get.center(rwat['p1'],rwat['p2'],rwat['p3'])
                    N_rwat = get.normal(rwat['p1'],rwat['p3'],rwat['p2'])
                    Q = get.rwatQ(N_rwat,ray,c_rwat)
                    OBS = get.obstruction(rwat,Q,ray['t'],N_rwat,ray['d'])            
                    if OBS == True:
                                    #print("Obstruction found in SHELLS break",ii,k)
                        break
                    else:
                        n = n+1
                                
            if OBS == False:
                #print "no obstruction found in element",ii,m,n 
                elem_noruk = np.append(elem_noruk,ii)
                thetas = np.append(thetas,ray_angle)
                elem_noruk_t = np.append(elem_noruk_t,ray['t'])
    return {'i':elem_noruk, 'thetas':thetas, 't':elem_noruk_t}       
    

"""
MAIN SCRIPT
"""
comm = MPI.COMM_WORLD
nprocs = comm.Get_size()
rank = comm.Get_rank()

with open('input') as ip:
    ip_lines = ip.readlines()

shells = mesh.Mesh.from_file(ip_lines[0].rstrip("\n"))
rukawat = mesh.Mesh.from_file(ip_lines[1].rstrip("\n"))
r_qcm = np.float(ip_lines[2].rstrip("\n"))
t_qcm = np.float(ip_lines[3].rstrip("\n"))

# Only master writes to Diagnostics
if rank == 0:
	diagnostics = open("DIAGNOSTICS","a+")	
	diagnostics.write("\r\n")
	diagnostics.write("Run Commenced at "+str((datetime.datetime.now()))+"\r\n")
src = get.source()
#elem = 50
elem = np.size(shells.areas)        # total number of shell elements
elem_r = np.size(rukawat.areas)     # number of elements in rukawat
start = time.time()	

if rank == 0:
	diagnostics.write("Number of cores %s\r\n" %nprocs)
	
n = elem/nprocs
iterable = list(range(nprocs))
OUTPUT = dict(i=np.array([]),thetas=np.array([]),t = np.array([]))
	
if rank == 0:
	diagnostics.write("Number of elements -  %s\r\n" %elem)
	diagnostics.write("Number of Obstructions - %s\r\n" %elem_r)

# CALL IT!
OP = calc_OBS_para(shells,rukawat,src,np.size(shells.areas),elem_r,OUTPUT['i'],OUTPUT['t'],OUTPUT['thetas'],int(n),rank)

# Barrier till all processes finish
comm.barrier()
print "Obstruction Calc barrier crossed by process %s" %rank

"""
OUTPUT HANDLING 
"""

# Send data to master
if rank != 0:
	req = comm.isend(OP, dest=0,tag=rank)
	req.wait()
print "OP Send barrier crossed by process %s" %rank
# At master piece data together
if rank == 0:
	OUTPUT = dict(i=np.array([]),thetas=np.array([]),t = np.array([]))
	OUTPUT['i'] = np.append(OUTPUT['i'],OP['i'])
	OUTPUT['thetas'] = np.append(OUTPUT['thetas'],OP['thetas'])
	OUTPUT['t'] = np.append(OUTPUT['t'],OP['t'])
	print "Output dictionary initialized by process %s" %rank
	for i in range(1,nprocs):
		req = comm.irecv(source=i,tag=i)
		OP_wrkr = req.wait()
		OUTPUT['i'] = np.append(OUTPUT['i'],OP_wrkr['i'])
		OUTPUT['thetas'] = np.append(OUTPUT['thetas'],OP_wrkr['thetas'])
		OUTPUT['t'] = np.append(OUTPUT['t'],OP_wrkr['t'])
	with open('OUTPUT.DAT','w+') as OP:
    		op_writer = csv.writer(OP, delimiter='\t')
    		for i in range(0,np.size(OUTPUT['i'])):
        		op_writer.writerow([OUTPUT['i'][i],OUTPUT['thetas'][i],OUTPUT['t'][i]])
	print("Execution time --- %s seconds ---\r\n" % (time.time() - start)) 
	diagnostics.write("Execution time --- %s seconds ---\r\n" % (time.time() - start))
	diagnostics.write("Number of elements with no obstruction - %s \r\n" %np.size(OUTPUT['i']))
	diagnostics.write("Ending run at "+str((datetime.datetime.now()))+"\r\n")
	diagnostics.close()
	t_i = t_qcm*np.cos(OUTPUT['thetas']*np.pi/180)*(r_qcm**2)/(OUTPUT['t']**2)
	print "thickness acquired"

# PLOT
	i = 0
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	#for e in range(0,np.size(shells.areas)):
	for e in OUTPUT['i']:
   	 # get the element with no obstruction
    		e = int(e)
    		x_p = list(shells.x[e,:])
   		y_p = list(shells.y[e,:])
    		z_p = list(shells.z[e,:])
#    get.patch(ax,x_p,y_p,z_p,1,0,1)
    		get.patch(ax,x_p,y_p,z_p,t_i[i],np.min(t_i),np.max(t_i))
    		i = i+1
	ax.scatter(src['NO'][0],src['NO'][1],src['NO'][2],s=50,c='r')
	get.view(ax,3)
	ax.set_zlim(-400,400)
	ax.set_ylim(-500,500)
	ax.set_xlim(-500,500)
	plt.show()



