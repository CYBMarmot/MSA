import numpy as np
import math as m
import pandas as pd
import matplotlib.pyplot as plt

#Function Verified
def local3DFrame(A,E,G,J,Iz,Iy,L):
    ek = np.array([[A*E/L,0,0,0,0,0,-A*E/L,0,0,0,0,0],[0,12*E*Iz/(L**3),0,0,0,6*E*Iz/(L**2),0,-12*E*Iz/(L**3),0,0,0,6*E*Iz/(L**2)],[0,0,12*E*Iy/(L**3),0,-6*E*Iy/(L**2),0,0,0,-12*E*Iy/(L**3),0,-6*E*Iy/(L**2),0],[0,0,0,J*G/L,0,0,0,0,0,-J*G/L,0,0],[0,0,-6*E*Iy/(L**2),0,4*E*Iy/L,0,0,0,6*E*Iy/(L**2),0,2*E*Iy/L,0],[0,6*E*Iz/(L**2),0,0,0,4*E*Iz/L,0,-6*E*Iz/(L**2),0,0,0,2*E*Iz/L],[-A*E/L,0,0,0,0,0,A*E/L,0,0,0,0,0],[0,-12*E*Iz/(L**3),0,0,0,-6*E*Iz/(L**2),0,12*E*Iz/(L**3),0,0,0,-6*E*Iz/(L**2)],[0,0,-12*E*Iy/(L**3), 0,6*E*Iy/(L**2),0,0,0,12*E*Iy/(L**3),0,6*E*Iy/(L**2),0],[0,0,0,-J*G/L,0,0,0,0,0,J*G/L,0,0],[0,0, -6*E*Iy/(L**2),0,2*E*Iy/L,0,0,0,6*E*Iy/(L**2), 0, 4*E*Iy/L,0],[0,6*E*Iz/(L**2), 0, 0, 0, 2*E*Iz/L, 0, -6*E*Iz/(L**2), 0, 0, 0, 4*E*Iz/L]])
    
    return ek

def global3DFrame(A,E,G,J,Iz,Iy,x1,y1,z1,x2,y2,z2,roll):
    xdif = x2-x1
    ydif = y2-y1
    zdif = z2-z1
    len3 = np.sqrt(xdif**2+ydif**2+zdif**2)
    rollang = np.radians(roll)
    xp = np.array([xdif/len3,ydif/len3,zdif/len3])
    yunit = np.array([0,0,1])
    ybar = np.cross(yunit,xp)

    RotMatr

    
    








