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
import datetime
from functools import partial
import csv
import sys

    
# import shells, obstruction and source


shells = mesh.Mesh.from_file('/media/amaan/Storage/Work/PPPL/Analysis/Li_Evap/shells/fortran_implementation/CAD_Source/shell.stl')
rukawat = mesh.Mesh.from_file('/media/amaan/Storage/Work/PPPL/Analysis/Li_Evap/shells/fortran_implementation/CAD_Source/Li_Evap_Export_Obs.STL')

elem = np.size(shells.areas)


with open("/media/amaan/Storage/Work/PPPL/Analysis/Li_Evap/shells/fortran_implementation/CAD_Source/Geom/shells.x.DAT","w") as X:
    for i in range(0,elem):
        X.write(str(shells.x[i][0]) + "\t" + str(shells.x[i][1]) + "\t" + str(shells.x[i][2]) + "\n")
        
with open("/media/amaan/Storage/Work/PPPL/Analysis/Li_Evap/shells/fortran_implementation/CAD_Source/Geom/shells.y.DAT","w") as X:
    for i in range(0,elem):
        X.write(str(shells.y[i][0]) + "\t" + str(shells.y[i][1]) + "\t" + str(shells.y[i][2]) + "\n")

with open("/media/amaan/Storage/Work/PPPL/Analysis/Li_Evap/shells/fortran_implementation/CAD_Source/Geom/shells.z.DAT","w") as X:
    for i in range(0,elem):
        X.write(str(shells.z[i][0]) + "\t" + str(shells.z[i][1]) + "\t" + str(shells.z[i][2]) + "\n")

with open("/media/amaan/Storage/Work/PPPL/Analysis/Li_Evap/shells/fortran_implementation/CAD_Source/Geom/shells.norm.DAT","w") as X:
    for i in range(0,elem):
        X.write(str(shells.normals[i][0]) + "\t" + str(shells.normals[i][1]) + "\t" + str(shells.normals[i][2]) + "\n")
elem = np.size(rukawat.areas)
with open("/media/amaan/Storage/Work/PPPL/Analysis/Li_Evap/shells/fortran_implementation/CAD_Source/Geom/ruk.x.DAT","w") as X:
    for i in range(0,elem):
        X.write(str(rukawat.x[i][0]) + "\t" + str(rukawat.x[i][1]) + "\t" + str(rukawat.x[i][2]) + "\n")
        
with open("/media/amaan/Storage/Work/PPPL/Analysis/Li_Evap/shells/fortran_implementation/CAD_Source/Geom/ruk.y.DAT","w") as X:
    for i in range(0,elem):
        X.write(str(rukawat.y[i][0]) + "\t" + str(rukawat.y[i][1]) + "\t" + str(rukawat.y[i][2]) + "\n")

with open("/media/amaan/Storage/Work/PPPL/Analysis/Li_Evap/shells/fortran_implementation/CAD_Source/Geom/ruk.z.DAT","w") as X:
    for i in range(0,elem):
        X.write(str(rukawat.z[i][0]) + "\t" + str(rukawat.z[i][1]) + "\t" + str(rukawat.z[i][2]) + "\n")

with open("/media/amaan/Storage/Work/PPPL/Analysis/Li_Evap/shells/fortran_implementation/CAD_Source/Geom/ruk.norm.DAT","w") as X:
    for i in range(0,elem):
        X.write(str(rukawat.normals[i][0]) + "\t" + str(rukawat.normals[i][1]) + "\t" + str(rukawat.normals[i][2]) + "\n")