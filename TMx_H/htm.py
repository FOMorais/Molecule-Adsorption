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
    if atoms[0]=='H':           
        a=np.array(coordinates[0])    
        n=n_atoms-1
        I=[]
        for x in range(n):
            a=x+1  #I is the vector with the TM8 positions
            I.append(coordinates[a])   
        V1=[]
        posH1=coordinates[0]                          
        for x in range(len(I)):
            V1.append(distance.euclidean(posH1, I[x])) #calculating the euclidan distances H1-TMx			
        V1.sort()
        V1=V1[:4]                  
        
        H=[]
        for x in range(len(I)):
            H.append(distance.euclidean(posH1  , I[x]))

        H=np.array(H)   #MATRIZ COM TODAS AS DISTANCIAS H - TM
        indmin = np.unravel_index(np.argmin(H, axis=None), H.shape)  #PEGANDO O VALOR MINIMO DESSA MATRIZ
        #print(indmin)                                               #POSICAO DO VALOR MINIMO NA MATRIZ
        #print(indmin[0])                                            #COM ISSO DESCOBRIMOS DE QUAL H ESSA DISTANCIA VEIO
        #print(indmin[1])                                            #COM ISSO DESCOBRIMOS QUAL TM FOI RESPONSAVEL POR ESSA DMIN
        TM_close=I[indmin[0]] 
        print(TM_close)                                                            #COORDENADAS DO TM MAIS PROXIMO DA MOLECULA
        
        R=[]
        for x in range(len(I)):
            R.append(distance.euclidean(TM_close, I[x]))
        print(R)
        m = min(y for y in R if y > 0)
        coff=m*1.25 #calculating the cutoff radius from m multiplied by a factor f
        H=[]
        print(coff)
        for k in range(len(R)):
            if (R[k]<=coff and  (R[k]>0)):
                H.append(k)
        E=[]
        for j in H:
            j=int(j)
            E.append(R[j])

        CN=len(H)
        dav=sum(E)/float(len(E))        
        

        
        R=[]
        for x in range(len(I)):
            R.append(distance.euclidean(posH1, I[x]))
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

        CNH=len(H)
        davH=sum(E)/float(len(E))        

        
        
        
        
        
        
                               
        print(str(filename)+"\t"+str(round(V1[0],5))+"\t"+str(round(V1[1],5))+"\t"+str(round(CN,5))+"\t"+str(round(dav,5))+"\t"+str(round(CNH,5))+"\t"+str(round(davH,5))+"\n")
        file1.write(str(filename)+"\t"+str(round(V1[0],5))+"\t"+str(round(V1[1],5))+"\t"+str(round(CN,5))+"\t"+str(round(dav,5))+"\t"+str(round(CNH,5))+"\t"+str(round(davH,5))+"\n")

       
file1.close()
