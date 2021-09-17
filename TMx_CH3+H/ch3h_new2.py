#Felipe Orlando Morais

import numpy as np
from scipy.spatial import distance
from glob import glob
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


    filename = str(i)
    #print(filename)         #just place the name of the xyz file here
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
            print("ATENCAO")
            pass



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
    for x in range(len(I)):
        R.append(distance.euclidean(posC, I[x])) #calculating the euclidan distances C-TMx 
#print(R)
    #print("The smaller distance C-TM is: "+str(min(R)))

#now for the four Hidrogens


    posH1=coordinates[0]
    V=[]
    I1=[]
    I2=[]
    I3=[]
    I4=[]
    for x in range(len(I)):
        I1.append(distance.euclidean(posH1, I[x]))
        V.append(distance.euclidean(posH1, I[x])) #calculating the euclidan distances H1-TMx 

    posH2=coordinates[2]
    for x in range(len(I)):
        I2.append(distance.euclidean(posH2, I[x]))
        V.append(distance.euclidean(posH2, I[x])) #calculating the euclidan distances H2-TMx 

    posH3=coordinates[3]
    for x in range(len(I)):
        I3.append(distance.euclidean(posH3, I[x]))
        V.append(distance.euclidean(posH3, I[x])) #calculating the euclidan distances H3-TMx 

    posH4=coordinates[4]
    for x in range(len(I)):
        I4.append(distance.euclidean(posH4, I[x]))
        V.append(distance.euclidean(posH4, I[x])) #calculating the euclidan distances H4-TMx 

#print("The smaller distance H-TM is: "+str(min(V)))
    if (min(V)) in I1:
        Hads=np.array(posH1)
        a=np.array(posH2)
        b=np.array(posC)
        c=np.array(posH3)
        d=np.array(posH4)
    if (min(V)) in I2:
        Hads=np.array(posH2)
        a=np.array(posH1)
        b=np.array(posC)
        c=np.array(posH3)
        d=np.array(posH4)
    if (min(V)) in I3:
        Hads=np.array(posH3)
        a=np.array(posH1)
        b=np.array(posC)
        c=np.array(posH2)
        d=np.array(posH4)
    if (min(V)) in I4:
        Hads=np.array(posH4)
        a=np.array(posH1)
        b=np.array(posC)
        c=np.array(posH2)
        d=np.array(posH3)


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

    a3=np.degrees(angle)
    HCHav=(a1+a2+a3)/3.
    # print("Mean HCH angle for CH3: "+str(HCHav))
    CHads=distance.euclidean(posC, Hads)

    #print("Distance C-Hads: "+str(CHads))

    
    #for x in range(len(I)):
    #    H.append(distance.euclidean(Hads, I[x])) #calculating the euclidan distances Hads-TMx
    #print("The smaller distance Hads-TM is: "+str(min(H)))
    
    H=[]
    #print(a)  #H1
    #print(c)  #H2
    #print(d)  #H3
    R=[]
    V1=[]
    V2=[]
    V3=[]
    
    for x in range(len(I)):
            R.append(distance.euclidean(b, I[x]))
            V1.append(distance.euclidean(a, I[x])) #calculating the euclidan distances H1-TMx
            V2.append(distance.euclidean(c, I[x])) #calculating the euclidan distances H2-TMx
            V3.append(distance.euclidean(d, I[x])) #calculating the euclidan distances H3-TMx
			
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
    dHs=[a,c,d]
                                          
    dH.append(distance.euclidean(a, TM_close))
    dH.append(distance.euclidean(c, TM_close))
    dH.append(distance.euclidean(d, TM_close))
    dH=np.array(dH)
    
    dHmin = np.argmin(dH)
    H_close=dHs[dHmin]
        
    dHmax = np.argmax(dH)
    H_far=dHs[dHmax]
    
    DHcTM=distance.euclidean(H_close, TM_close)
    DCtm= distance.euclidean(b , TM_close)
    DHfTM=distance.euclidean(H_far  , TM_close)      
        
    Hdied=wiki_dihedral(H_far,posC,H_close,TM_close)
    
    #####H-C-TM
    ba = np.array(H_close) - np.array(posC)
    bd = np.array(TM_close) - np.array(posC)

    cosine_angle = np.dot(ba, bd) / (np.linalg.norm(ba) * np.linalg.norm(bd))
    angle = np.arccos(cosine_angle)

    HCTM=np.degrees(angle)
    
    #COORDENACAO TM PROX  
    print(TM_close)  
    R=[]
    for x in range(len(I)):
            R.append(distance.euclidean(TM_close, I[x]))
    #print(R)  
         
    m = min(y for y in R if y > 0)
    #print(m) 
    coff=m*1.25 #calculating the cutoff radius from m multiplied by a factor f
    H=[]
    #print(coff)    
    for k in range(len(R)):
            if (R[k]<=coff and  (R[k]>0)):
                    H.append(k)
    E=[]
    for j in H:
            j=int(j)
            E.append(R[j])

    CNTM=len(H)
    davTM=sum(E)/float(len(E))
    
    #COORDENACAO Hads   
    R=[]
    for x in range(len(I)):
            R.append(distance.euclidean(Hads, I[x]))
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

    CNHads=len(H)
    davHads=sum(E)/float(len(E))
    
    
    
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
    
    print(str(filename)+"\t"+str(round(HCHav,5))+"\t"+str(round(DHcTM,5))+"\t"+str(round(DHfTM,5))+"\t"+str(round(DCtm,5))+"\t"+str(round(abs(CHads),5))+"\t"+str(round(abs(Hdied),5))+"\t"+str(round(abs(HCTM),5))+"\t"+str(round(abs(CNTM),5))+"\t"+str(round(abs(davTM),5))+"\t"+str(round(abs(CNHads),5))+"\t"+str(round(abs(davHads),5))+"\t"+str(round(abs(CNC),5))+"\n")
    file1.write(str(filename)+"\t"+str(round(HCHav,5))+"\t"+str(round(DHcTM,5))+"\t"+str(round(DHfTM,5))+"\t"+str(round(DCtm,5))+"\t"+str(round(abs(CHads),5))+"\t"+str(round(abs(Hdied),5))+"\t"+str(round(abs(HCTM),5))+"\t"+str(round(abs(CNTM),5))+"\t"+str(round(abs(davTM),5))+"\t"+str(round(abs(CNHads),5))+"\t"+str(round(abs(davHads),5))+"\t"+str(round(abs(CNC),5))+"\n")
file1.close()
