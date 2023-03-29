#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 23 19:29:58 2019
Measure the overlap energy of all residues normalized by number of contacts
-All non-bonded atoms, non-backbone inter residue overlaps
@author: agrigas115
"""
import sys
import numpy as np
import itertools

backbone_atoms_list = ['N', 'CA', 'C', 'O', 'OXT', 'H1', 'H2', 'H3', 'H',
                       'HA']

bonded_list = ['NCA', 'NC', 'NO', 'CAC', 'CAO', 'CO',
               'CAH1', 'CAH2', 'CAH3', 'CH1', 'CH2', 'CH3',
               'OH1', 'OH2', 'OH3', 'CAH', 'CH', 'OH',
               'NHA', 'CHA', 'OHA', 'H1H2', 'H1H3', 'H2H3',
               'HAH', 'HAH1', 'HAH2', 'HAH3', #backbone-backbone interactions
               
               'CACB', 'CBCG1', 'CBCG2', 
               'NH1', 'NH2', 'NH3', 'CAHA', 'CBHB', 'CG1HG11', 
               'CG1HG12', 'CG1HG13', 'CG2HG21', 'CG2HG22', 
               'CG2HG23', #VAL
               
               'CBOG', 'NH', 'CBHB2', 'CBHB3', 'OGHG', #SER
               
               'CGCD1', 'CGCD2', 'CD1CE1', 'CD2CE2', 'CE1CZ',
               'CE2CZ', 'CZOH', 'OHHH', 'CD1HD1', 'CD2HD2',
               'CE1HE1', 'CE2HE2', 'CZHH', #TYR
               
               'CGOD1', 'CGOD2', #ASP
               
               'CAHA2', 'CAHA3', #GLY
               
               'CGND1', 'CGCD2', 'ND1CE1', 'CD2NE2', 'CD2HD2',
               'CE1HE1', 'CE1NE2', #HIS
               
               'CZHZ', #PHE
               
               'CGCD1', 'CGCD2', 'CGHG', 'CD1HD11', 'CD1HD12', 'CD1HD13',
               'CD2HD21', 'CD2HD22', 'CD2HD23', #LEU
               
               'CBOG1', 'CHCG2', 'OG1HG1', 'CG2HG21', 'CG2HG22', 
               'CG2HG23', #THR
               
               'CBCG', 'CGCD', 'CDCE', 'CENZ', 'CBHB2', 'CBHB3',
               'CGHG2', 'CGHG3', 'CDHD2', 'CDHD3', 'CEHE2', 'CEHE3',
               'NZHZ1', 'NZHZ2', 'NZHZ3',  #LYS
               
               'CG1CD1', 'CD1HD11', 'CD1HD12', 'CD1HD13', #ILE
               
               'CGND2', 'ND2HD21', 'ND2HD22', #ASN
               
               'CDNE', 'NECZ', 'CZNH1', 'CZNH2', 'NEHE', 
               'NH1HH11', 'NH1HH12', 'NH2HH21', 'NH2HH22', #ARG
               
               'CDOE1', 'CDOE2', #GLU
               
               'CDNE2', 'NE2HE21', 'NE2HE22', #GLN
               
               'NCD', 'CDHD2', 'CDHD3', #PRO
               
               'CBHB1', 'CBHB2', 'CBHB3', #ALA
               
               'CBSG', 'SGHG', #CYS
               
               'CGSD', 'SDCE', 'CDHE1', 'CDHE2', 'CDHE3', 'CEHE1', #MET
               
               'CGCD1', 'CD1NE1', 'NE1CE2', 'CE2CD2', 'CE2CZ2',
               'CD2CG', 'CD2CE3', 'CE3CZ3', 'CZ3CH2', 'CH2CZ2', 
               'CD1HD1', 'NE1HE1', 'CZ2HZ2', 'CH2HH2', 'CZ3HZ3',
               'CE3HE3', #TRP
               
               'COXT' #Terminus, i dunno Jennifer's code should remove this atom type
                       #but i had a case where it has this atom type
               ]
#List of all bonded atoms

def U(r, rad_i, rad_j):
    sigma_ij = rad_i+rad_j
    if sigma_ij-r > 0:
        return (1.0/72.0)*(1.0-(sigma_ij/r)**6.0)**2.0
    else:
        return 0.0
#repulsive lennard-jones potential

decoy = str(sys.argv[1])
decoy = '../pdb/'+decoy

# read in command line argument for decoy structure

overlaps_list = []
measured_list = []
name = decoy[7:-4]

energy_list = []
print('Working on ' + decoy[7:])
with open(decoy, 'r') as f:
    info = f.read()
sourcelines = info.splitlines()
decoy_data = []
for line in sourcelines:
    decoy_data.append([int(line.split()[2]), float(line.split()[3]), float(line.split()[4]),
                       float(line.split()[5]), float(line.split()[6]), line.split()[0]])
#reading in the coordinate data for the decoy structure and saving a list
#of the [resid, radius, x, y, z]

previous_resid = decoy_data[0][0]

total_energy = 0
overlap_id_list = []
total_overlap_dict = {}
last_resid = decoy_data[-1][0]
last_resid_list = []
for entry in decoy_data:
    if entry[0] == last_resid:
        last_resid_list.append(entry)
for pair in itertools.permutations(range(0, len(decoy_data)), 2):
    atom_1 = decoy_data[pair[0]]
    atom_2 = decoy_data[pair[1]]
    
    current_resid = atom_1[0]
    
    resid_1 = atom_1[0]
    resid_2 = atom_2[0]

    atomtype_1 = atom_1[5]
    atomtype_2 = atom_2[5]
    if not atomtype_1 in backbone_atoms_list: #first atom is sidechain
        if (resid_1 != resid_2 #Are you on different residues?
        and atomtype_1+atomtype_2 != 'NC' 
        and atomtype_2+atomtype_1 != 'NC'): #Are you the only atoms bonded between residues?
            x_1 = float(atom_1[2])
            y_1 = float(atom_1[3])
            z_1 = float(atom_1[4])
            x_2 = float(atom_2[2])
            y_2 = float(atom_2[3])
            z_2 = float(atom_2[4])
            
            r=np.sqrt((x_2-x_1)**2+(y_2-y_1)**2+(z_2-z_1)**2)
            
            rad_1 = float(atom_1[1])
            rad_2 = float(atom_2[1])
            
            overlap_energy = U(r, rad_1, rad_2)
            
            if current_resid == previous_resid:
                if overlap_energy > 0:
                    overlaps_list.append([atomtype_1+'_'+str(resid_1)+atomtype_2+'_'+str(resid_2), overlap_energy])
                    total_energy += overlap_energy
                    overlap_id_list.append(atom_2[0])
                else:
                    total_energy += 0

            else:
                if len((overlap_id_list)) == 0:
                    energy_list.append([previous_resid, 0])
                    total_energy = overlap_energy
                    overlaps_list.append([str(previous_resid), '0'])
                    previous_resid = current_resid
                else:
                    energy_list.append([previous_resid,total_energy/len(set(overlap_id_list))])
                    total_energy = overlap_energy
                    previous_resid = current_resid
                    overlap_id_list = []
            
# =============================================================================
#                 elif (resid_1 == resid_2 #Are you on the same residue?
#                 and atomtype_2 in backbone_atoms_list #second atom is backbone
#                 and not (atomtype_1 in ['CB', 'HB', 'HB2', 'HB3'])
#                 and not (atomtype_1+atomtype_2 in bonded_list) 
#                 and not (atomtype_2+atomtype_1 in bonded_list) #Are you bonded?
#                 and atomtype_1+atomtype_2 != atomtype_1+atomtype_1 
#                 and atomtype_1+atomtype_2 != atomtype_2+atomtype_2 #Are you self overlap?
#                 and not (atomtype_1+atomtype_2+str(current_resid) in intra_measured_list)): #Have you been measured yet?
#                     intra_measured_list.append(atomtype_2+atomtype_1+str(current_resid))
#                     #print(current_resid)
#                     #print(atomtype_1+atomtype_2)
#                     
#                     x_1 = float(atom_1[2])
#                     y_1 = float(atom_1[3])
#                     z_1 = float(atom_1[4])
#                     x_2 = float(atom_2[2])
#                     y_2 = float(atom_2[3])
#                     z_2 = float(atom_2[4])
#                     
#                     r=np.sqrt((x_2-x_1)**2+(y_2-y_1)**2+(z_2-z_1)**2)
#                     
#                     rad_1 = float(atom_1[1])
#                     rad_2 = float(atom_2[1])
#                     
#                     overlap_energy = U(r, rad_1, rad_2)
#                     #print(overlap_energy)
#                     if overlap_energy > 0:
#                         intra_overlaps_list.append([atomtype_1+'_'+str(resid_1)+atomtype_2+'_'+str(resid_2), overlap_energy])
#                     if current_resid == previous_resid:
#                         if overlap_energy > 0:
#                             total_energy += overlap_energy
#                             overlap_id_list.append(atom_2[0])
#                         else:
#                             total_energy += 0
#         
#                     else:
#                         if len(set(overlap_id_list)) == 0:
#                             energy_list.append([previous_resid, 0])
#                             total_energy = overlap_energy
#                             previous_resid = current_resid
#                         else:
#                             energy_list.append([previous_resid,total_energy/len(set(overlap_id_list))])
#                             total_energy = overlap_energy
#                             previous_resid = current_resid
#                             overlap_id_list = [atom_2[0]]
# =============================================================================
### Same thing but just for the last residue because how I iter over
#### misses the last residue
overlap_id_list = []
total_energy = 0.0
for atom_1 in last_resid_list:
    resid_1 = atom_1[0]
    atomtype_1 = atom_1[5]
    
    x_1 = float(atom_1[2])
    y_1 = float(atom_1[3])
    z_1 = float(atom_1[4])
    for atom_2 in decoy_data:
        resid_2 = atom_2[0]
        atomtype_2 = atom_2[5]
        if not atomtype_1 in backbone_atoms_list: #first atom is sidechain
            if (resid_1 != resid_2 #Are you on different residues?
            and atomtype_1+atomtype_2 != 'NC' 
            and atomtype_2+atomtype_1 != 'NC'): #Are you the only atoms bonded between residues?
                x_2 = float(atom_2[2])
                y_2 = float(atom_2[3])
                z_2 = float(atom_2[4])
                
                r=np.sqrt((x_2-x_1)**2+(y_2-y_1)**2+(z_2-z_1)**2)
                
                rad_1 = float(atom_1[1])
                rad_2 = float(atom_2[1])
                
                overlap_energy = U(r, rad_1, rad_2)
                total_energy += overlap_energy
                if overlap_energy > 0.0 :
                    overlaps_list.append([atomtype_1+'_'+str(resid_1)+atomtype_2+'_'+str(resid_2), overlap_energy])
                    overlap_id_list.append(resid_2)
            
# =============================================================================
#                 elif (resid_1 == resid_2 
#                   and atomtype_2 in backbone_atoms_list #second atom is backbone
#                   and not (atomtype_1 in ['CB', 'HB', 'HB2', 'HB3'])
#                   and not (atomtype_1+atomtype_2 in bonded_list) 
#                   and not (atomtype_2+atomtype_1 in bonded_list) 
#                   and atomtype_1+atomtype_2 != atomtype_1+atomtype_1 
#                   and atomtype_1+atomtype_2 != atomtype_2+atomtype_2):
#                     x_1 = float(atom_1[2])
#                     y_1 = float(atom_1[3])
#                     z_1 = float(atom_1[4])
#                     x_2 = float(atom_2[2])
#                     y_2 = float(atom_2[3])
#                     z_2 = float(atom_2[4])
#                     
#                     r=np.sqrt((x_2-x_1)**2+(y_2-y_1)**2+(z_2-z_1)**2)
#                     
#                     rad_1 = float(atom_1[1])
#                     rad_2 = float(atom_2[1])
#                     
#                     overlap_energy = U(r, rad_1, rad_2)
#                     if overlap_energy > 0:
#                         intra_overlaps_list.append([atomtype_1+atomtype_2, overlap_energy])
#                     if current_resid == previous_resid:
#                         if overlap_energy > 0:
#                             total_energy += overlap_energy
#                             overlap_id_list.append(atom_2[0])
#                         else:
#                             total_energy += 0
# =============================================================================

if len((overlap_id_list)) == 0:
    energy_list.append([last_resid, 0])
    overlaps_list.append([str(last_resid), '0'])
else:
    energy_list.append([last_resid,total_energy/len(set(overlap_id_list))])
        
with open('../overlaplist_data/'+decoy[7:-4]+'_overlaplist.txt', 'w') as f:
    for i in range(0, len(energy_list)):
        f.write(str(energy_list[i][0])+'\t'+str(energy_list[i][1]))
        f.write('\n')
#Write overlap data

with open('../contactlist_data/'+decoy[7:-4]+'_contactlist.txt', 'w') as f:
    for i in range(0, len(overlaps_list)):
        line_format = '{:<13}  {:<13}'.format(overlaps_list[i][0], overlaps_list[i][1])
        f.write(line_format)
        f.write('\n')
#Write a file of the atoms with intra-residue overlaps and greater
#energy than 0
        
