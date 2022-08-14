import mdtraj as md
import numpy as np
import math
import pickle
from Bio.PDB import PDBParser
from Bio import PDB
import pandas as pd

parser = PDBParser(PERMISSIVE=1)

#with open("list_structures.temp") as file:
    #lines=file.readlines()
    #pdbLIG = parser.load(file=lines[0].split()[1])
structure = parser.get_structure("lig", "lig.pdb")
vecs=[]
crd=57.29577951
h,w = 1000001, 5;
puck=np.zeros((h, 200))
traj=md.load('lig.pdb')
nres=traj.n_residues
sugt=[2,8,4,5,1,7,3,9,0,6]
ang=[]
for model in structure:
        for residue in model.get_residues():
                residue_list = str(residue.get_full_id()[3][1])
                Co1 = residue["C1'"].get_vector()
                Co2 = residue["C2'"].get_vector()
                Co3 = residue["C3'"].get_vector()
                Co4 = residue["C4'"].get_vector()
                Co5 = residue["O4'"].get_vector()
                coord=[Co1,Co2,Co3,Co4,Co5]
            
                angle = PDB.calc_dihedral(Co1,Co2,Co3,Co4)
                angle2 = PDB.calc_dihedral(Co2,Co3,Co4,Co5)
                angle3 = PDB.calc_dihedral(Co3,Co4,Co5,Co1)
                angle4 = PDB.calc_dihedral(Co4,Co5,Co1,Co2)
                angle5 = PDB.calc_dihedral(Co5,Co1,Co2,Co3)
                angles=[angle,angle2,angle3,angle4,angle5]
                ang.append(angles)
                np.savetxt('ang.dat',ang,delimiter='  ', fmt='%f')
df= pd.read_csv('ang.dat', skiprows=0, delimiter="\s+",header=None)

num=[0,1,2,3,4]
cosdata=[]
sindata=[]
for i in range(0,5):
    cos=[math.cos(144*i/crd)]
    sin=[math.sin(144*i/crd)]
    
    cosdata += cos
    sindata += sin


df_trans=df.T
for x in df_trans:
    get_a=lambda x: x*cosdata*crd
    get_b=lambda x:x*sindata*crd

df2=df_trans.apply(get_a)
df3=df_trans.apply(get_b)
Total = df2.sum()
Total2=df3.sum()
a=Total*0.4
b=Total2*(-0.4)
np.savetxt('a.dat',a,delimiter='  ', fmt='%f')
np.savetxt('b.dat',b,delimiter='  ', fmt='%f')
mula = a.mul(a, axis=0)
mulb=b.mul(b, axis=0)
amp=(mula+mulb)**(1/2)

    




all = pd.concat([a, b, amp], axis=1)
amp_pos=all.apply(lambda x: x[2] if x[2]>30 else 0, axis=1)
all = pd.concat([a, b, amp, amp_pos], axis=1)

cp=all.apply(lambda x: x[0]/x[3] if x[3]>0 else 0, axis=1)
all = pd.concat([a, b, amp, amp_pos,cp], axis=1)
sp=all.apply(lambda x: x[1]/x[3] if x[3]>0 else 0, axis=1)
all = pd.concat([a, b, amp, amp_pos,cp, sp], axis=1)
pha=all.apply(lambda x: ((math.cos(x[4])**(-1)))*crd if x[3]>0 else 0, axis=1)
all = pd.concat([a, b, amp, amp_pos, cp, sp, pha], axis=1)
pha2=all.apply(lambda x: -x[6] if x[5]<0 else 0, axis=1)
all = pd.concat([a, b, amp, amp_pos, cp, sp, pha,pha2], axis=1)

phase=all.apply(lambda x: int(x[7]+360) if x[7]<0 else x[6], axis=1)
all = pd.concat([a, b, amp, amp_pos, cp, sp, pha,pha2,phase], axis=1)
sugar=all.apply(lambda x: int(x[8]/36), axis=1)
all = pd.concat([a, b, amp, amp_pos, cp, sp, pha,pha2,phase,sugar], axis=1)
print(all)
print(sugar.values)
sugar1=[]
sugar1.append(sugar.values)
np.savetxt('pucker.dat',sugar1,delimiter='  ', fmt='% 4d')

#print(all[0])
#pha=np.arctan(df[3]/df[4])
#print(pha)



#phase=amp.where(amp<0, lambda a: a/amp)
#print(phase)

#amp=m666666ath.sqrt(a.values*a.values)
#print(amp)