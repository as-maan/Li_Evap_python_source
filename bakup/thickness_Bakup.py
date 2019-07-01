# -*- coding: utf-8 -*-
"""
Created on Sun Jun 16 13:38:20 2019

@author: amaan
"""
import numpy as np
import get
from stl import mesh
import matplotlib.pyplot as plt
import multiprocessing as mp
import time
from functools import partial


    
def calc_OBS(shells,rukawat,src,elem,elem_r,elem_noruk,elem_noruk_t,thetas,i):
    #for ii in range(i*n,(i+1)*n):
        # replace this with elem
        # get an element p
    p = get.p(shells.x[i,:],shells.y[i,:],shells.z[i,:])
        
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
            N_rwat = get.normal(rwat['p1'],rwat['p2'],rwat['p3'])
                # Find Q, the point along ray that lies on the plane of rwat
            Q = get.rwatQ(N_rwat,ray,c_rwat)
                # Check if Q belongs to rwat and if yes is t_Q < t
            OBS = get.obstruction(rwat,Q,ray['t'],N_rwat)          
            if OBS == True:
    #                print('Obstruction found in rukawat, break',i)
                break
            else:
                m = m+1
            
                
            # check for obstruction from shell elements
        if OBS == False:
    #            print('now checking against all shell elements', i)
                
            for k in [u for u in range(0,elem) if u!=i]:
                rwat = get.p(shells.x[k,:],shells.y[k,:],shells.z[k,:])
                c_rwat = get.center(rwat['p1'],rwat['p2'],rwat['p3'])
                N_rwat = get.normal(rwat['p1'],rwat['p2'],rwat['p3'])
                Q = get.rwatQ(N_rwat,ray,c_rwat)
                OBS = get.obstruction(rwat,Q,ray['t'],N_rwat)            
                if OBS == True:
    #                    print('Obstruction found in shells, break',i)
                    break
                else:
                    n = n+1
                
        if OBS == False:
            print "no obstruction found in element",i,m,n 
            elem_noruk = np.append(elem_noruk,i)
            thetas = np.append(thetas,ray_angle)
            elem_noruk_t = np.append(elem_noruk_t,ray['t'])

    return {'i':elem_noruk, 'thetas':thetas, 't':elem_noruk_t}        
        
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
                N_rwat = get.normal(rwat['p1'],rwat['p2'],rwat['p3'])
                                # Find Q, the point along ray that lies on the plane of rwat
                Q = get.rwatQ(N_rwat,ray,c_rwat)
                                # Check if Q belongs to rwat and if yes is t_Q < t
                OBS = get.obstruction(rwat,Q,ray['t'],N_rwat)          
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
                    N_rwat = get.normal(rwat['p1'],rwat['p2'],rwat['p3'])
                    Q = get.rwatQ(N_rwat,ray,c_rwat)
                    OBS = get.obstruction(rwat,Q,ray['t'],N_rwat)            
                    if OBS == True:
                                    #print("Obstruction found in SHELLS break",ii,k)
                        break
                    else:
                        n = n+1
                                
            if OBS == False:
                print "no obstruction found in element",ii,m,n 
                elem_noruk = np.append(elem_noruk,ii)
                thetas = np.append(thetas,ray_angle)
                elem_noruk_t = np.append(elem_noruk_t,ray['t'])
    return {'i':elem_noruk, 'thetas':thetas, 't':elem_noruk_t}       
    
    


# import shells, obstruction and source
shells = mesh.Mesh.from_file('/media/amaan/Storage/Work/PPPL/Analysis/Li_Evap/shells/CAD_Source/dummy_shells_2.stl')
rukawat = mesh.Mesh.from_file('/media/amaan/Storage/Work/PPPL/Analysis/Li_Evap/shells/CAD_Source/dummy_Obs.stl')
#print("Geometry acquired")
src = get.source()

elem = np.size(shells.areas)        # total number of shell elements
elem_r = np.size(rukawat.areas)     # number of elements in rukawat

r_qcm = 1.0
t_qcm = 1.0  


start = time.time()

# Parallel Execution
num_cores = 4
n = elem/num_cores
iterable = list(range(num_cores))
pool = mp.Pool(num_cores)
OUTPUT = dict(i=np.array([]),thetas=np.array([]),t = np.array([]))
#print("Number of elements", elem)
#print("Number of Obstructions", elem_r)
func = partial(calc_OBS_para,shells,rukawat,src,np.size(shells.areas),elem_r,OUTPUT['i'],OUTPUT['t'],OUTPUT['thetas'],n)
OP = pool.map(func,iterable)
pool.close()
pool.join()
# Stitch outputs together

for k in iterable:
    OUTPUT['i'] = np.append(OUTPUT['i'],OP[k]['i'])
    OUTPUT['thetas'] = np.append(OUTPUT['thetas'],OP[k]['thetas'])
    OUTPUT['t'] = np.append(OUTPUT['t'],OP[k]['t'])

        
"""
# Series Execution
elem_noruk = np.array([])
elem_noruk_t = np.array([])
thetas = np.array([])
for i in range(0,elem):
    OUTPUT = calc_OBS(shells,rukawat,src,elem,elem_r,elem_noruk,elem_noruk_t,thetas,i)
    elem_noruk = OUTPUT['i']
    thetas = OUTPUT['thetas']
    elem_noruk_t = OUTPUT['t']
 
"""
print("Execution time --- %s seconds ---" % (time.time() - start))         



    
t_i = t_qcm*np.cos(OUTPUT['thetas']*np.pi/180)*(r_qcm**2)/(OUTPUT['t']**2)

# make plot
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

get.view(ax,3)

ax.set_zlim(-2,2)
ax.set_ylim(-2,4)
ax.set_xlim(0,10)
plt.show()

'''
ax.set_zlim(-400,400)
ax.set_ylim(-500,500)
ax.set_xlim(-500,500)
plt.show()
'''    

