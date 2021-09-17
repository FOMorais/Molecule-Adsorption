import os
from os import listdir
from os.path import isfile, join
import sys
from glob import glob
import pandas as pd
import time
import numpy as np
from scipy.spatial import distance
from scipy.linalg import svdvals
import random


def new_dihedrall(p0,p1,p2,p3):
    """Praxeolitic formula
    1 sqrt, 1 cross product"""
    p0 = np.array(p0)
    p1 = np.array(p1)
    p2 = np.array(p2)
    p3 = np.array(p3)

    b0 = -1.0*(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= np.linalg.norm(b1)

    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x))

def wiki_dihedral(p0,p1,p2,p3):
    """formula from Wikipedia article on "Dihedral angle"; formula was removed
    from the most recent version of article (no idea why, the article is a
    mess at the moment) but the formula can be found in at this permalink to
    an old version of the article:
    https://en.wikipedia.org/w/index.php?title=Dihedral_angle&oldid=689165217#Angle_between_three_vectors
    uses 1 sqrt, 3 cross products
    """
    p0 = np.array(p0)
    p1 = np.array(p1)
    p2 = np.array(p2)
    p3 = np.array(p3)

    b0 = -1.0*(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    b0xb1 = np.cross(b0, b1)
    b1xb2 = np.cross(b2, b1)

    b0xb1_x_b1xb2 = np.cross(b0xb1, b1xb2)

    y = np.dot(b0xb1_x_b1xb2, b1)*(1.0/np.linalg.norm(b1))
    x = np.dot(b0xb1, b1xb2)

    return np.degrees(np.arctan2(y, x))



file1 = open('table_feat.txt', 'w')

baseFolder = glob("*.xyz")
baseFolder.sort()
for i in baseFolder:

    filename = str(i)         #just place the name of the xyz file here
    atoms = []
    coordinates = []
    xyz = open(filename)
    n_atoms = int(xyz.readline())
    empty = xyz.readline()
    for line in xyz:
        atom,x,y,z = line.split()
        atoms.append(atom)
        coordinates.append([float(x), float(y), float(z)])
    xyz.close()


    if not (str(atoms[0]) == 'H'):
        i=atoms.index('H')
        coordinates[i],coordinates[0],atoms[i],atoms[0]=coordinates[0],coordinates[i],atoms[0],atoms[i]

    if not (str(atoms[1]) == 'C'):
        ii=atoms.index('C')
        coordinates[ii],coordinates[1],atoms[ii],atoms[1]=coordinates[1],coordinates[ii],atoms[1],atoms[ii]

    if not (str(atoms[2]) == 'H'):
        iii=[i for i, n in enumerate(atoms) if n == 'H'][1]
        coordinates[iii],coordinates[2],atoms[iii],atoms[2]=coordinates[2],coordinates[iii],atoms[2],atoms[iii]

    if not (str(atoms[3]) == 'H'):
        iv=[i for i, n in enumerate(atoms) if n == 'H'][2]
        coordinates[iv],coordinates[3],atoms[iv],atoms[3] = coordinates[3],coordinates[iv],atoms[3],atoms[iv]

    if not (str(atoms[4]) == 'H'):
        try:
           v=[i for i, n in enumerate(atoms) if n == 'H'][3]
           coordinates[v],coordinates[4],atoms[v],atoms[4] = coordinates[4],coordinates[v],atoms[4],atoms[v]
        except IndexError:
            pass

#print(atoms)
#print(coordinates)


#print("Number of atoms:  %d" % n_atoms)


    if atoms[0]=='H':           
        a=np.array(coordinates[0])
    if atoms[1]=='C':               
        b=np.array(coordinates[1])  
    if atoms[2]=='H':               
        c=np.array(coordinates[2])  
    if atoms[3]=='H':               
        d=np.array(coordinates[3])  
    if atoms[4]=='H':               
        e=np.array(coordinates[4])  
        #print(i)
       
        bc = c - b
        be = e - b

        cosine_angle = np.dot(bc, be) / (np.linalg.norm(bc) * np.linalg.norm(be))
        angle = np.arccos(cosine_angle)

        a5=np.degrees(angle)

        bd = d - b
        bc = e - b

        cosine_angle = np.dot(bd, be) / (np.linalg.norm(bd) * np.linalg.norm(be))
        angle = np.arccos(cosine_angle)

        a6=np.degrees(angle)

        ba = a - b
        bc = c - b

        cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
        angle = np.arccos(cosine_angle)

        a1=np.degrees(angle)

        ba = a - b
        bd = d - b

        cosine_angle = np.dot(ba, bd) / (np.linalg.norm(ba) * np.linalg.norm(bd))
        angle = np.arccos(cosine_angle)

        a2=np.degrees(angle)

        ba = a - b
        be = e - b

        cosine_angle = np.dot(ba, be) / (np.linalg.norm(ba) * np.linalg.norm(be))
        angle = np.arccos(cosine_angle)

        a3=np.degrees(angle)

        bd = d - b
        bc = c - b

        cosine_angle = np.dot(bd, bc) / (np.linalg.norm(bd) * np.linalg.norm(bc))
        angle = np.arccos(cosine_angle)

        a4=np.degrees(angle)
        HCHav=(a1+a2+a3+a4+a5+a6)/6.
        #print("Mean HCH angle for CH4: "+str(HCHav))
    #print((a1+a2+a3+a4+a5+a6)/6.)
    
        n=n_atoms-5
    #print("We have "+str(n)+" "+atoms[5]+" atoms.")
    
        indexC = atoms.index('C')
        posC=coordinates[indexC]
        I=[]
        for x in range(n):
            a=x+5  #I is the vector with the TM8 positions
            I.append(coordinates[a])
    #print(I)
        R=[]
        V1=[]
        V2=[]
        V3=[]
        V4=[]
        #RR=[]
        """
        for x in range(len(I)):
            R.append(distance.euclidean(posC, I[x])) #calculating the euclidan distances C-TMx 
        
        """
        posH1=coordinates[0]
        posH2=coordinates[2]
        posH3=coordinates[3]
        posH4=coordinates[4]
        
        
        
        
        for x in range(len(I)):
            R.append(distance.euclidean(posC  , I[x]))
            V1.append(distance.euclidean(posH1, I[x])) #calculating the euclidan distances H1-TMx
            V2.append(distance.euclidean(posH2, I[x])) #calculating the euclidan distances H2-TMx
            V3.append(distance.euclidean(posH3, I[x])) #calculating the euclidan distances H3-TMx
            V4.append(distance.euclidean(posH4, I[x])) #calculating the euclidan distances H4-TMx
			
        H=[]
        
        H.append(V1)
        H.append(R)
        H.append(V2)
        H.append(V3)
        H.append(V4)
        H=np.array(H)   #MATRIZ COM TODAS AS DISTANCIAS H,C - TM
        #print(H)
        indmin = np.unravel_index(np.argmin(H, axis=None), H.shape)  #PEGANDO O VALOR MINIMO DESSA MATRIZ
        #print(indmin)                                               #POSICAO DO VALOR MINIMO NA MATRIZ
        #print(indmin[0])                                            #COM ISSO DESCOBRIMOS DE QUAL H ESSA DISTANCIA VEIO
        #print(indmin[1])                                            #COM ISSO DESCOBRIMOS QUAL TM FOI RESPONSAVEL POR ESSA DMIN
        TM_close=I[indmin[1]] 
        #print("tm close")                             #COORDENADAS DO TM MAIS PROXIMO DA MOLECULA
        #print(TM_close) 
        
        dH=[] 
        dHs=[posH1,posH2,posH3,posH4]
                                          
        dH.append(distance.euclidean(posH1, TM_close))
        dH.append(distance.euclidean(posH2, TM_close))
        dH.append(distance.euclidean(posH3, TM_close))
        dH.append(distance.euclidean(posH4, TM_close))
        #print(dH)
        #print(dHs)
        dH=np.array(dH)
        dHmin = np.argmin(dH)
        #print("h mais prox do tm close")
        #print(dHmin)
        H_close=dHs[dHmin]
        #print(H_close)
        
        dHmax = np.argmax(dH)
        #print("h mais longe do tm close")
        #print(dHmax)
        H_far=dHs[dHmax]
        #print(H_far)
        
        #print("distancia minima")
        DHcTM=distance.euclidean(H_close, TM_close)
        DCtm=distance.euclidean(posC  , TM_close)
        #print("distancia max")
        DHfTM=distance.euclidean(H_far  , TM_close)      
        
        Hdied=wiki_dihedral(H_far,posC,H_close,TM_close)
        
        ba = np.array(H_close) - np.array(posC)
        bd = np.array(TM_close) - np.array(posC)

        cosine_angle = np.dot(ba, bd) / (np.linalg.norm(ba) * np.linalg.norm(bd))
        angle = np.arccos(cosine_angle)

        HCTM=np.degrees(angle)
        
        R=[]
        for x in range(len(I)):
            R.append(distance.euclidean(TM_close, I[x]))
        m = min(y for y in R if y > 0)
        coff=m*1.25 #calculating the cutoff radius from m multiplied by a factor f
        H=[]
        
        for k in range(len(R)):
            if (R[k]<=coff and  (R[k]>0)):
                H.append(k)
        E=[]
        for j in H:
            j=int(j)
            E.append(R[j])

        CN=len(H)
        dav=sum(E)/float(len(E))
        
        
        #####do C
        R=[]
        for x in range(len(coordinates)):
            R.append(distance.euclidean(posC, coordinates[x]))
        m = min(y for y in R if y > 0)
        coff=m*1.25 #calculating the cutoff radius from m multiplied by a factor f
        H=[]
        
        for k in range(len(R)):
            if (R[k]<=coff and  (R[k]>0)):
                H.append(k)
        E=[]
        for j in H:
            j=int(j)
            E.append(R[j])

        CNC=len(H)
        #davC=sum(E)/float(len(E))
        
        #print(str(filename)+"\t"+str(round(HCHav,5))+"\t"+str(round(DHcTM,5))+"\t"+str(round(DHfTM,5))+"\t"+str(round(DCtm,5))+"\t"+str(round(abs(Hdied),5))+"\t"+str(CN)+"\t"+str(round(abs(dav),5))+"\t"+str(round(abs(HCTM),5))+"\n")
        file1.write(str(filename)+"\t"+str(round(HCHav,4))+"\t"+str(round(DHcTM,4))+"\t"+str(round(DHfTM,4))+"\t"+str(round(DCtm,4))+"\t"+str(round(abs(Hdied),4))+"\t"+str(CN)+"\t"+str(round(abs(dav),4))+"\t"+str(round(abs(HCTM),4))+"\t"+str(round(abs(CNC),4))+"\n")

    else:

        #print(filename)
        
        ba = a - b
        bc = c - b

        cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
        angle = np.arccos(cosine_angle)

        a1=np.degrees(angle)

        ba = a - b
        bd = d - b

        cosine_angle = np.dot(ba, bd) / (np.linalg.norm(ba) * np.linalg.norm(bd))
        angle = np.arccos(cosine_angle)

        a2=np.degrees(angle)

        bd = d - b
        bc = c - b

        cosine_angle = np.dot(bd, bc) / (np.linalg.norm(bd) * np.linalg.norm(bc))
        angle = np.arccos(cosine_angle)

        a4=np.degrees(angle)
        
        HCHav=(a1+a2+a4)/3.
        #print("Mean HCH angle for CH3: "+str(HCHav))
       
        n=n_atoms-4
    
        indexC = atoms.index('C')
        posC=coordinates[indexC]
        I=[]
        for x in range(n):
            a=x+4  #I is the vector with the TM8 positions
            I.append(coordinates[a])
        R=[]
        V1=[]
        V2=[]
        V3=[]

        posH1=coordinates[0]
        posH2=coordinates[2]
        posH3=coordinates[3]
        
        
        
        
        for x in range(len(I)):
            R.append(distance.euclidean(posC  , I[x]))
            V1.append(distance.euclidean(posH1, I[x])) #calculating the euclidan distances H1-TMx
            V2.append(distance.euclidean(posH2, I[x])) #calculating the euclidan distances H2-TMx
            V3.append(distance.euclidean(posH3, I[x])) #calculating the euclidan distances H3-TMx
			
        H=[]
        
        H.append(V1)
        H.append(R)
        H.append(V2)
        H.append(V3)
        H=np.array(H)   #MATRIZ COM TODAS AS DISTANCIAS H,C - TM
        #print(H)
        indmin = np.unravel_index(np.argmin(H, axis=None), H.shape)  #PEGANDO O VALOR MINIMO DESSA MATRIZ
        #print(indmin)                                               #POSICAO DO VALOR MINIMO NA MATRIZ
        #print(indmin[0])                                            #COM ISSO DESCOBRIMOS DE QUAL H ESSA DISTANCIA VEIO
        #print(indmin[1])                                            #COM ISSO DESCOBRIMOS QUAL TM FOI RESPONSAVEL POR ESSA DMIN
        TM_close=I[indmin[1]] 
                                                                    #COORDENADAS DO TM MAIS PROXIMO DA MOLECULA
        
        
        dH=[] 
        dHs=[posH1,posH2,posH3]
                                          
        dH.append(distance.euclidean(posH1, TM_close))
        dH.append(distance.euclidean(posH2, TM_close))
        dH.append(distance.euclidean(posH3, TM_close))

        dH=np.array(dH)
        dHmin = np.argmin(dH)
        H_close=dHs[dHmin]
        
        dHmax = np.argmax(dH)
        H_far=dHs[dHmax]
        DHcTM=distance.euclidean(H_close, TM_close)
        DCtm=distance.euclidean(posC  , TM_close)
        DHfTM=distance.euclidean(H_far  , TM_close)      
        
        Hdied=wiki_dihedral(H_far,posC,H_close,TM_close)
        
        
        ba = np.array(H_close) - np.array(posC)
        bd = np.array(TM_close) - np.array(posC)

        cosine_angle = np.dot(ba, bd) / (np.linalg.norm(ba) * np.linalg.norm(bd))
        angle = np.arccos(cosine_angle)

        HCTM=np.degrees(angle)
        
        R=[]
        for x in range(len(I)):
            R.append(distance.euclidean(TM_close, I[x]))
        m = min(y for y in R if y > 0)
        coff=m*1.25 #calculating the cutoff radius from m multiplied by a factor f
        H=[]
        
        for k in range(len(R)):
            if (R[k]<=coff and  (R[k]>0)):
                H.append(k)
        E=[]
        for j in H:
            j=int(j)
            E.append(R[j])

        CN=len(H)
        dav=sum(E)/float(len(E))
        #####do C
        R=[]
        for x in range(len(coordinates)):
            R.append(distance.euclidean(posC, coordinates[x]))
        m = min(y for y in R if y > 0)
        coff=m*1.25 #calculating the cutoff radius from m multiplied by a factor f
        H=[]
        
        for k in range(len(R)):
            if (R[k]<=coff and  (R[k]>0)):
                H.append(k)
        E=[]
        for j in H:
            j=int(j)
            E.append(R[j])

        CNC=len(H)
        #davC=sum(E)/float(len(E))
        #print(CNC, davC)


        #print(str(filename)+"\t"+str(round(HCHav,5))+"\t"+str(round(DHcTM,5))+"\t"+str(round(DHfTM,5))+"\t"+str(round(DCtm,5))+"\t"+str(round(abs(Hdied),5))+"\t"+str(CN)+"\t"+str(round(abs(dav),5))+"\t"+str(round(abs(HCTM),5))+"\n")
        file1.write(str(filename)+"\t"+str(round(HCHav,4))+"\t"+str(round(DHcTM,4))+"\t"+str(round(DHfTM,4))+"\t"+str(round(DCtm,4))+"\t"+str(round(abs(Hdied),4))+"\t"+str(CN)+"\t"+str(round(abs(dav),4))+"\t"+str(round(abs(HCTM),4))+"\t"+str(round(abs(CNC),4))+"\n")
        

       
file1.close()
