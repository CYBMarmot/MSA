import numpy as np
import math as m
import pandas as pd
import matplotlib.pyplot as plt

#Function Verified
def klocal(x1,y1,x2,y2,A,E):
    lenx = x2-x1
    leny = y2-y1
    lenT = np.sqrt(lenx**2+leny**2)
    C = lenx/lenT
    S = leny/lenT
    stiffScalar = E*A/lenT
    ForceCoefMatr = np.array([[C**2,S*C,-C**2,-S*C],[S*C,S**2,-S*C,-S**2],[-C**2,-S*C,C**2,S*C],[-S*C,-S**2,S*C,S**2]])
    kCo = ForceCoefMatr
    Lk = stiffScalar * kCo
    return Lk, kCo

#Function Verified
def kglobal(coords, areas, moduli, con):
    mem_num = np.shape(areas)[1]
    con_num = np.shape(con)[0]
    total_con = np.max(con)
    dof = 2*total_con
    Gk = np.zeros((dof,dof))

    for i in range(con_num):
        n1 = con[i,0]
        n2 = con[i,1]
        x1 = coords[n1-1,0]
        y1 = coords[n1-1,1]
        x2 = coords[n2-1,0]
        y2 = coords[n2-1,1]
        A = areas[0,i]
        E = moduli[0,i]
        ka = klocal(x1,y1,x2,y2,A,E)[0]   
        indexvec = np.array([[2*n1-2,2*n1-1,2*n2-2,2*n2-1]])
        indexLen = np.shape(indexvec)[1]
        for j in range(indexLen):
            for k in range(indexLen):
                row = indexvec[0,j]
                col = indexvec[0,k]
                Gk[row,col] = Gk[row,col] + ka[j,k]
    return Gk

#Function verified
def orgElements(P,kG):
    vecsize = np.shape(kG)[0]
    P_num = np.shape(P)[0]
    f_orig = np.zeros((vecsize,1))
    for i in range(P_num):
        node = P[i,0]
        if P[i,1] == 1:
            f_orig[2*node-2] = P[i,2]
        if P[i,1] == 2:
            f_orig[2*node-1] = P[i,2]
    return f_orig

#Function Verified
def dumShift(kG,Bc,P):
    f_orig = orgElements(P,kG)
    f1 = np.copy(f_orig)
    itera = np.shape(Bc)[0]
    for i in range(itera):
        node = Bc[i,0]
        if Bc[i,1] == 1:
            index = 2*node-2
        if Bc[i,1] == 2:
            index = 2*node-1
        kCol = np.transpose(np.array([kG[:,index]]))
        cdisp = Bc[i,2]
        f1 = f1 - kCol*cdisp
    return f1, f_orig

#Function Verified
def kglobal2f2(Bc,P,coords, areas, moduli, con):
    kG = kglobal(coords, areas, moduli, con)
    f1 = dumShift(kG,Bc,P)[0]
    f_orig = dumShift(kG,Bc,P)[1]
    f2 = np.copy(f1)
    itera = np.shape(Bc)[0]
    itera1 = np.shape(kG)[0]
    kG2 = np.copy(kG)
    for i in range(itera):
        node = Bc[i,0]
        if Bc[i,1] == 1:
            index = 2*node-2
        if Bc[i,1] == 2:
            index = 2*node-1  
        for j in range(itera1):
            kG2[j,index] = 0
            kG2[index,j] = 0
        kG2[index,index] = kG[index,index]
        rdisp = Bc[i,2]
        f2[index] = rdisp*kG[index,index]

    return kG2, f2, kG, f_orig

#Function Verified
def solveTruss(Bc,P,coords, areas, moduli, con):
    mats = kglobal2f2(Bc,P,coords, areas, moduli, con)
    kG = mats[2]
    kG2 = mats[0]
    f2 = mats[1]
    
    disp = np.linalg.inv(kG2) @ f2
    Force = kG  @ disp

    veclnum = np.shape(areas)[1]
    pAxial = np.zeros((veclnum,1))
    fAxial = np.zeros((veclnum,1))
    memID = np.zeros((veclnum,1))
    lengths = np.zeros((veclnum,1))
    con_num = np.shape(con)[0]
    
    for i in range(con_num):
        n1 = con[i,0]
        n2 = con[i,1]
        x1 = coords[n1-1,0]
        y1 = coords[n1-1,1]
        x2 = coords[n2-1,0]
        y2 = coords[n2-1,1]
        xdif = x2-x1
        ydif = y2 - y1
        lens = np.sqrt(xdif**2+ydif**2)
        A = areas[0,i]
        E = moduli[0,i]
        ka = klocal(x1,y1,x2,y2,A,E)[0]   
        indexvec = np.array([[2*n1-2,2*n1-1,2*n2-2,2*n2-1]])
        indexLen = np.shape(indexvec)[1]
        TabE = np.zeros((indexLen,1))
        for j in range(indexLen):
            TabE[j,0] = disp[indexvec[0,j]]
        TabF = ka @ TabE
        
        N = (TabE[0]**2 + TabE[1]**2)**0.5
        c1 = np.absolute(np.sin(np.arctan(ydif/xdif)))
        c2 = np.absolute(np.cos(np.arctan(ydif/xdif)))
        if (c1 < 0.0000001):
            Nc = TabF[2]/c2
        else:
            Nc = TabF[3]/c1
            
        memID[i,0] = i
        lengths[i,0] = lens
        pAxial[i,0] = Nc
        fAxial[i,0] = np.absolute(Nc/A)

    Rtable1 = np.concatenate((lengths, pAxial, fAxial), axis=1)
    RtablePres1 = pd.DataFrame(Rtable1, columns=['Length', 'Force','Stress'])

    supports = np.array([np.unique(Bc[:,1])])
  

    itera3 = np.shape(supports)[1]
    xreact = np.zeros((itera3,1))
    yreact = np.zeros((itera3,1))
    nodeID = np.zeros((itera3,1))

    for k in range(itera3):
        nodenum = supports[0,k]
        index = np.array([[2*nodenum-2,2*nodenum-1]])
        xreact[k,0] = Force[index[0,0],0]
        yreact[k,0] = Force[index[0,1],0]
        nodeID[k,0] = nodenum

    Rtable2 = np.concatenate((nodeID, xreact, yreact), axis=1)
    RtablePres2 = pd.DataFrame(Rtable2, columns=['Node', 'Fx','Fy']) 

    itera4 = np.shape(coords)[0]
    x_loc = np.transpose(np.array([coords[:,0]]))
    y_loc = np.transpose(np.array([coords[:,1]]))
    dispX = np.zeros((itera4,1))
    dispY = np.zeros((itera4,1))
    id3 = np.zeros((itera4,1))

    for l in range(itera4):
        nodeval = l+1
        id3[l,0] = nodeval
        dispX[l,0] = disp[2*nodeval-2,0]
        dispY[l,0] = disp[2*nodeval-1,0]

    Rtable3 = np.concatenate((id3,x_loc,y_loc,dispX, dispY), axis=1)
    RtablePres3 = pd.DataFrame(Rtable3, columns=['Node','GX', 'GY',' GX Displacement', 'GY Displacement']) 
        
    return Rtable1, RtablePres1, Rtable2, RtablePres2, Rtable3, RtablePres3, dispX, dispY,lengths,fAxial

# Function verified, this is now working well.
def plotResults(con,coords,scale1,dispX,dispY):

    M = np.shape(con)[0]
    plt.figure(figsize=(6, 6))

    dispX = np.copy(dispX)*scale1
    dispY = np.copy(dispY)*scale1
    
    for i in range(M):
        fnode = con[i,0]
        snode = con[i,1]
        x1 = coords[fnode-1,0]
        x1d = coords[fnode-1,0] + dispX[fnode-1,0] 
        x2 = coords[snode-1,0]
        x2d =  coords[snode-1,0] + dispX[snode-1,0]
        y1 = coords[fnode-1,1]
        y1d = coords[fnode-1,1] + dispY[fnode-1,0]
        y2 = coords[snode-1,1]
        y2d = coords[snode-1,1] + dispY[snode-1,0]
        xs = np.array([[x1],[x2]])
        ys = np.array([[y1],[y2]])
        xs1 = np.array([[x1d],[x2d]])
        ys1 = np.array([[y1d],[y2d]])
        plt.plot(xs, ys, color = 'purple', linestyle = 'solid', markersize = 8)
        plt.plot(xs1, ys1 , color ='red', linestyle = 'dashed',markersize = 8)

    plt.title("Truss Structure")
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.grid(True)
    plt.axis('equal')
    plt.show()
    
    
        

    

    
    








