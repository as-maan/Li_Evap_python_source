# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 17:46:43 2019

@author: amaan
"""
import numpy as np
import math
import matplotlib as mpl
import statistics
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from pycallgraph import PyCallGraph
from pycallgraph.output import GraphvizOutput

def source():
#    NO = np.array([0.0,0.0,0.0])
    NO = 25.4*np.array([-3.3601,-16.8923,0.3908])
    FG = 25.4*np.array([3.6954,18.5778,0.3908])    
    return {'NO':NO, 'FG':FG}

def p(x,y,z):
    [x1,x2,x3] = x
    [y1,y2,y3] = y
    [z1,z2,z3] = z
    p1 = [x1,y1,z1]
    p2 = [x2,y2,z2]
    p3 = [x3,y3,z3]
    return p1,p2,p3
    
def center(p1,p2,p3):
    x = [p1[0],p2[0],p3[0]]
    y = [p1[1],p2[1],p3[1]]
    z = [p1[2],p2[2],p3[2]]
    c = [sum(x)/float(len(x)),sum(y)/float(len(y)),sum(z)/float(len(z))]
    #print c
    return c

def normal(N):
    n = -N/np.sqrt(N.dot(N))
    #n = np.cross(np.array(p3)-np.array(p1),np.array(p2)-np.array(p1))/np.linalg.norm(np.cross(np.array(p3)-np.array(p1),np.array(p2)-np.array(p1)))
    return n

def sourcevector(p_s,c_p):
    t = np.sqrt((c_p-p_s).dot(c_p-p_s))

#    t = np.linalg.norm(c_p-p_s)
    d = (c_p-p_s)/t
    return p_s,d,t

def angle(N_p,d):
    norm_N = np.sqrt(N_p.dot(N_p))
#    norm_N = np.linalg.norm(N_p)
    norm_d = np.sqrt(d.dot(d))
#    norm_d = np.linalg.norm(d)
    theta = (180/math.pi)*(math.acos(np.dot(N_p,d)/(norm_N*norm_d)))

    return theta

def rwatQ(N_rwat,p_s,d,t,c_rwat):
    d_scalar = np.dot(N_rwat,c_rwat)
    t_Q = (d_scalar-np.dot(N_rwat,p_s))/(np.dot(N_rwat,d))
    p_Q = p_s + t_Q*d
    return p_Q,t_Q       

def obstruction(rwat1,rwat2,rwat3,p_Q,t_Q,t,N_rwat,d):

	cond0 = np.dot(-N_rwat,d)

	a = [x1 - x2 for (x1, x2) in zip(rwat2, rwat1)]
	b = [x1 - x2 for (x1, x2) in zip(p_Q.tolist(),rwat1)]
	n = [a[1]*b[2] - a[2]*b[1],a[2]*b[0] - a[0]*b[2],a[0]*b[1] - a[1]*b[0]]

	cond1 = np.dot(np.array(n),N_rwat)

	a = [x1 - x2 for (x1, x2) in zip(rwat3,rwat2)]
	b = [x1 - x2 for (x1, x2) in zip(p_Q.tolist(),rwat2)]
	n = [a[1]*b[2] - a[2]*b[1],a[2]*b[0] - a[0]*b[2],a[0]*b[1] - a[1]*b[0]]
  
	cond2 = np.dot(np.array(n),N_rwat)

	a = [x1 - x2 for (x1, x2) in zip(rwat1,rwat3)]
	b = [x1 - x2 for (x1, x2) in zip(p_Q.tolist(),rwat3)]
	n = [a[1]*b[2] - a[2]*b[1],a[2]*b[0] - a[0]*b[2],a[0]*b[1] - a[1]*b[0]]

	cond3 = np.dot(np.array(n),N_rwat)

	cond4 = t - t_Q 
	cond  = [cond4,cond0,cond1,cond2,cond3]
	if all(i>=0 for i in cond) == True:
		OBS = True
	else:
		OBS = False
	return OBS
        
def patch(ax, x, y, z, v, vmin, vmax, cmap_name='viridis'):
    cmap = mpl.cm.get_cmap(cmap_name)               # Get colormap by name
    c = cmap(mpl.colors.Normalize(vmin, vmax)(v))   # Normalize value and get color
    pc = Poly3DCollection([list(zip(x,y,z))])       # Create PolyCollection from coords
    pc.set_facecolor(c)                             # Set facecolor to mapped value
    pc.set_edgecolor('none')                           # Set edgecolor to black
    ax.add_collection3d(pc)                         # Add PolyCollection to axes
    return pc

def view(ax, code):
    if code == 2: #view(2) sets the default two-dimensional view, az = 0, el = 90.
        ax.view_init(0, 180)     # (args are reversed from MATLAB)

    if code == 3: #view(3) sets the default three-dimensional view, az = –37.5, el = 30.
        ax.view_init(-30, -37.5)

    if code == 1: #view(3) sets the default three-dimensional view, az = –37.5, el = 30.
        ax.view_init(0, 0)


