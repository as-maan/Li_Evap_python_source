# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 17:46:43 2019

@author: amaan
"""
import numpy as np
import math
import matplotlib as mpl
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
def source():
    NO = np.array([0.0,0.0,0.0])
#    NO = 25.4*np.array([-3.3601,-16.8923,0.3908])
    FG = 25.4*np.array([3.6954,18.5778,0.3908])    
    return {'NO':NO, 'FG':FG}

def p(x,y,z):
    [x1,x2,x3] = x
    [y1,y2,y3] = y
    [z1,z2,z3] = z
    p1 = np.array([x1,y1,z1])
    p2 = np.array([x2,y2,z2])
    p3 = np.array([x3,y3,z3])
    return {'p1':p1, 'p2':p2, 'p3':p3}
    
def center(p1,p2,p3):
    x = [p1[0],p2[0],p3[0]]
    y = [p1[1],p2[1],p3[1]]
    z = [p1[2],p2[2],p3[2]]
    c = [np.average(x),np.average(y),np.average(z)]
    return c

def normal(p1,p2,p3):
    n = np.cross(p2-p1,p3-p1)/np.linalg.norm(np.cross(p2-p1,p3-p1))
    return n

def sourcevector(p_s,c_p):
    t = np.linalg.norm(c_p-p_s)
    d = (c_p-p_s)/t
    return {'p_s':p_s, 'd':d, 't':t}

def angle(N_p,d):
    theta = (180/math.pi)*(math.acos(np.dot(N_p,d)/(np.linalg.norm(N_p)*np.linalg.norm(d))))
    return theta

def rwatQ(N_rwat,R,c_rwat):
    d = np.dot(N_rwat,c_rwat)
    t_Q = (d-np.dot(N_rwat,R['p_s']))/(np.dot(N_rwat,R['d']))
    p_Q = R['p_s'] + t_Q*R['d']
    return {'p_Q':p_Q, 't_Q':t_Q}       

def obstruction(rwat,Q,t,N_rwat):
    cond1 = np.dot(np.cross((rwat['p2'] - rwat['p1']),(Q['p_Q']-rwat['p1'])),N_rwat)
    cond2 = np.dot(np.cross((rwat['p3'] - rwat['p2']),(Q['p_Q']-rwat['p2'])),N_rwat)
    cond3 = np.dot(np.cross((rwat['p1'] - rwat['p3']),(Q['p_Q']-rwat['p3'])),N_rwat)
    cond4 = t - Q['t_Q'] 
    cond  = [cond1,cond2,cond3,cond4]
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
        ax.view_init(90, 0)     # (args are reversed from MATLAB)

    if code == 3: #view(3) sets the default three-dimensional view, az = –37.5, el = 30.
        ax.view_init(30, -37.5)

