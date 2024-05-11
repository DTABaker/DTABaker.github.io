#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import sys
sys.path.append('/home/dtbaker/.local/lib/python3.9/site-packages')
from tabulate import tabulate
import Bio
from Bio.PDB import *
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# In[2]:


yourPDB=input('Please input your PDB id. Ex: 4jpp \n')
pdbfile=PDBList()
pdbf= pdbfile.retrieve_pdb_file(yourPDB,pdir="/home/dtbaker/jupyter_notebook/pdbfiles/")
print(pdbf)
print(yourPDB)


# In[3]:


test1=pdbf.split(".")[1]
print(test1)
pdb_id = pdbf.split(".")[0]
print(pdb_id)
if test1 == "cif":
    parser = MMCIFParser()
    structure = parser.get_structure(pdb_id, pdbf)
    print(type(structure))
    for model in structure:
        print(model)
        for chain in model:
            print(chain)
            for residue in chain:
                print(residue)
                #for atom in residue:
                    #print(atom)
elif test1 == "pdb":
    parser = PDBParser(PERMISSIVE = True, QUIET = False)
    structure = parser.get_structure(pdb_id, pdbf)
    print(type(structure))
    for model in structure:
        print(model)
        for chain in model:
            print(chain)
            for residue in chain:
                print(residue)
                #for atom in residue:
                    #print(atom)
elif test1 == "mmtf":
    parser = MMTFParser()
    structure = MMTFParser.get_structure("PDB/pdbf")

elif test1 == "ent":
    pqr_parser = PDBParser(PERMISSIVE=1, is_pqr=True)
    structure = parser.get_structure(pdb_id, pdbf, is_pqr=True)


# In[4]:


print(structure.header["name"])


# In[5]:


print(structure.header["resolution"])
#print(structure.header["release_date"] or structure.header["deposition_date"])


# In[6]:


keywords = structure.header.keys()
print(keywords)


# In[7]:


print(pdbf)
print(pdb_id)

#This line removes the directory location subtext to pdbfile-for my directory
id_w_oDir=pdbf.split("/")[5]
print(id_w_oDir)


# In[8]:


#This confirms the file type of cif and prints some common searched cells
if test1 == "cif":
    mmcif_dict = MMCIF2Dict(pdbf)
    for key, value in mmcif_dict.items():
        if "?" not in value:
            if "_entity." in key:
                print(key, value)
            if "_cell." in key:
                print(key, value)
    #solv_cont= mmcif_dict["_expt1_crystal.density_percent_sol"]


# In[9]:


#Two ways of selecting particular things of a chain
res_list=Selection.unfold_entities(chain, "R")
atom_list=Selection.unfold_entities(chain, "A")

print(type(structure))
print(chain)
#print(chain2)
print(residue)

#Used to count number of specific parts 
print(len(res_list))
print(len(atom_list))
print("\n")

#shows polypeptides present
ppb = CaPPBuilder()
for pp in ppb.build_peptides(structure):
    print(pp.get_sequence())


# In[10]:


#prints the number of chains in a PDB model (Dr.Wang code)
num_chain=0
for model in structure:
    chain_info=dict()
    for chain in model:
        num_chain=num_chain +1
          
            
print(num_chain)           


# In[11]:


#prints the number of chains in a PDB model 
num_chain=0
chainIdsList=[]
for model in structure:
    chain_info=dict()
    for chain in model:
        chainIdsList=chain.id
        num_chain=num_chain +1
        print(chainIdsList)
        
#Below does the same as the for loop more efficiently          
#print(len(list(structure.get_chains())))
print("The number of chains:",num_chain)


# In[12]:


#Dr.Wang's code on extracting pdb info
for model in structure:
    chain_info=dict()
    for chain in model:
        print("chain id is", chain.id, chain.get_id)
        chain_info[chain.id]=dict()
        chain_info[chain.id]["chain ID"]=chain.id
        chain_info[chain.id]["total residues"]=len(list(chain.get_residues()))
        AA_res=0
        AA_wat=0
        AA_non=0
        AA_non_names = []
        
        for residue in chain:
            if residue.id[0]==' ':
                AA_res += 1
            elif residue.id[0]=='W':
                AA_wat += 1
            else:
                AA_non += 1
                print(residue.resname)
                AA_non_names.append(residue.resname)
            chain_info[chain.id]["Regular _AA number"]=AA_res
            chain_info[chain.id]["Water Number"]=AA_wat
            chain_info[chain.id]["Non-standard AA number"]=AA_non
            
            if AA_non !=0:
                chain_info[chain.id]["names of Non-standard AA"]=AA_non_names
chain_info


# In[13]:


print(list(chain))
print("\n")
print("len list of chain", len(list(chain)))
print("len list of model", len(list(model)))
print("len of list of residue", len(list(residue)))


# In[14]:


import pandas as pd
chainDF = pd.DataFrame(chain_info)
chainDF


# In[15]:


#Counts the number of chains
print(len(list(structure.get_chains())))


# In[16]:


#Calculates the CA distances within a single chain, pairwise Residues
def calculate_ca_distances_within_chain(chain):
    ca_atoms = [residue['CA'] for residue in chain if 'CA' in residue]
    distances = []
    for i in range(len(ca_atoms)):
        for j in range(i + 1, len(ca_atoms)):
            distance = ca_atoms[i] - ca_atoms[j]
            distances.append({'Chain:':chain.id, 'Residue_i':i+1, 'Residue_j':j+1, 'Distance':distance})
            
    return pd.DataFrame(distances)


#Calculates the CA distance between two chain, pairwise Residues
def calculate_ca_distances_between_chains(chain1, chain2):
    ca_atoms1 = [residue['CA'] for residue in chain1 if 'CA' in residue]
    ca_atoms2 = [residue['CA'] for residue in chain2 if 'CA' in residue]
    distances1 = []
    for i, atom1 in enumerate(ca_atoms1):
        for j, atom2 in enumerate(ca_atoms2):
            distance1 = atom1 - atom2
            distances1.append({'Chain 1:': chain1.id, 'Residue_i':i+1, 'Chain 2:': chain2.id, 'Residue_j':j+1, 'Distance':distance1})
    return pd.DataFrame(distances1)


# In[17]:


#Chain identifiers for easier calling.
chain1=structure[0]["A"]
chain2=structure[0]["B"]
chain3=structure[0]["C"]
chain4=structure[0]["D"]
chain5=structure[0]["E"]


# In[18]:


#Calling of functions and saved. Below tabbed out code outputs the save to a CSV file.

distanceAA=calculate_ca_distances_within_chain(chain1)
#distanceAA.to_csv('WithinChainOutput.csv', index=False)
distanceBB=calculate_ca_distances_within_chain(chain2)
distanceCC=calculate_ca_distances_within_chain(chain3)
distanceDD=calculate_ca_distances_within_chain(chain4)
distanceEE=calculate_ca_distances_within_chain(chain5)



distanceAB=calculate_ca_distances_between_chains(chain1, chain2)
distanceBC=calculate_ca_distances_between_chains(chain2, chain3)
distanceCD=calculate_ca_distances_between_chains(chain3, chain4)
distanceDE=calculate_ca_distances_between_chains(chain4, chain5)

#distanceAB.to_csv('BetweenChainOutput.csv', index=False)


# In[19]:


distanceAB


# In[20]:


def calculate_ca_distances(structure):
    ca_distances = []
    chains = list(structure.get_chains())
    
    # Loop through each pair of chains without repeating pairs
    for i, chain_i in enumerate(chains):
        for j, chain_j in enumerate(chains):
            if i < j:
                ca_atoms_i = [atom for atom in chain_i.get_atoms() if atom.get_name() == 'CA']
                ca_atoms_j = [atom for atom in chain_j.get_atoms() if atom.get_name() == 'CA']
                
                # Calculate distances for all CA atom pairs between two chains
                for atom_i in ca_atoms_i:
                    for atom_j in ca_atoms_j:
                        distance = atom_i - atom_j
                        ca_distances.append({
                            'Chain_i': chain_i.id, 'Residue_i': atom_i.get_parent().id[1],
                            'Chain_j': chain_j.id, 'Residue_j': atom_j.get_parent().id[1],
                            'Distance': distance
                        })

    # Convert the list of distances to a DataFrame
    return pd.DataFrame(ca_distances)


distance_df = calculate_ca_distances(structure)
#distance_df.to_csv('DfOutput.csv', index=False)
print(distance_df)


# In[21]:


plt.figure(figsize=(12, 8))

plt.hist(distanceAA['Distance'], bins=25, color='blue', label='Distance AA', histtype='step', linewidth=2, alpha=0.7)
plt.hist(distanceBB['Distance'], bins=25, color='red', label='Distance BB', histtype='step', linewidth=2, alpha=0.7)
plt.hist(distanceCC['Distance'], bins=25, color='green', label='Distance CC', histtype='step', linewidth=2, alpha=0.7)
plt.hist(distanceDD['Distance'], bins=25, color='cyan', label='Distance DD', histtype='step', linewidth=2, alpha=0.7)
plt.hist(distanceEE['Distance'], bins=25, color='purple', label='Distance EE', histtype='step', linewidth=2, alpha=0.7)


plt.legend(title='Chains:')


plt.xlabel('Cα Distance (Å)')
plt.ylabel('Frequency')
plt.title('Cα Distances for each chain in 4jpp')

plt.show()


# In[22]:


plt.figure(figsize=(12, 8))

plt.hist(distanceAB['Distance'], bins=25, color='blue', label='Distance AB', histtype='step', linewidth=1.5, alpha=0.5)
plt.hist(distanceBC['Distance'], bins=25, color='red', label='Distance BC', histtype='step', linewidth=1.5, alpha=0.5)
plt.hist(distanceCD['Distance'], bins=25, color='green', label='Distance CD', histtype='step', linewidth=1.5, alpha=0.5)
plt.hist(distanceDE['Distance'], bins=25, color='cyan', label='Distance DE', histtype='step', linewidth=1.5, alpha=0.5)



plt.legend(title='Chains:')


plt.xlabel('Cα Distance (Å)')
plt.ylabel('Frequency')
plt.title('Cα Distances for each chain pair in 4jpp')

plt.show()


# In[ ]:




