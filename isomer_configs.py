#!/usr/bin/python
"""
Created on Thu Oct 13 13:15:04 2016

@author: nenian charles
"""

from pymatgen.io.vasp import Poscar
from pymatgen.io.cif import CifParser
import os, itertools, os.path

"""User please provide some basic info about the units you are searching for """
cations_indx = [] #If you know the indices of the focal cation, list them
cati = 'Na' # If you don't have the indices provide the name of the species
ani = 'O' # provide the name of the anion specie you are searching for
cutoff = 2.5 #Angstroms

f = open("isomer_conformations.txt", 'w')
f.close()


#Search for derivative structure folders
for dictionary in next(os.walk('.'))[1]:
    f = open("isomer_conformations.txt", 'a')    
    os.chdir(dictionary)
    f.write("{}\n".format(dictionary))
    f.write("structure #\tIsomer conformation\n")
    f.write("-------------------------------------------------------------------------------\n")
    structures_list = [file for file in next(os.walk('.'))[2] if file.endswith(".vasp")]
    #print len(structures_list)
    if len(structures_list) == 0:
        structures_list = [file for file in next(os.walk('.'))[2] if file.endswith(".cif")]
        if len(structures_list) == 0:
            print('no structures found')
            pass
        else:
            file_format = 'cif'
    else:
        file_format = 'poscar'
        
    isomer_count = {'cis':0, 'trans':0, 'fac':0, 'mer':0, 'mixed':0}
    for entry in structures_list:
        struct_name = entry.replace("{}-".format(dictionary), '')
        struct_name = struct_name.replace('.vasp','')
        #print struct_name
        if file_format == 'poscar':
            poscar = Poscar.from_file(entry)
            struct = poscar.structure
        elif file_format == 'cif':
            parser = CifParser(entry)
            struct = parser.get_structures()[0]
        #cations_indx = [0,3,5,6]
        if len(cations_indx) == 0:            
            cations_indx = struct.indices_from_symbol(cati)
            
        anions_indx = struct.indices_from_symbol(ani)
    
        cation_neighs = dict()
        ordering_rank = []
        for cation in cations_indx:
            myneigh = []
            for anion in anions_indx:
                distance = struct.get_distance(cation,anion)
                if distance <= cutoff:
                    myneigh.append(anion)
            ordering_rank.append(len(myneigh))
    
            cation_neighs[cation] = myneigh
    
        ordering_rank = list(set(ordering_rank))
        if len(ordering_rank) > 1:
            print ("somethings wrong")
            break
        else:
            ordering_rank = ordering_rank[0]
        isomer_conformations = []    
        if ordering_rank == 1:
            print ("anion order doesn't support cis, trans, fac or mer")
            pass
        if ordering_rank == 2:
            site_angles = []
            for cation in cations_indx:
                angle = struct.get_angle(cation_neighs[cation][0], cation, cation_neighs[cation][1])
                #site_angles.append(angle)
                if angle == 90:
                    isomer_conformations.append('cis')
                elif angle == 180 or angle == 0:
                    isomer_conformations.append('trans')
            if len(set(isomer_conformations)) == 1:
                isomer = list(set(isomer_conformations))[0]
                f.write("{}\t\t{}\n".format(struct_name, isomer))
                isomer_count[isomer] += 1
            elif len(set(isomer_conformations)) > 1:
                f.write("{}\t\tmixed [cis and trans]\n".format(struct_name))
                isomer_count['mixed'] += 1
            #print site_angles    
        if ordering_rank == 3:
            for cation in cations_indx:
                site_angles = []
                for anion_pair in itertools.combinations(cation_neighs[cation],2):
                    angle = struct.get_angle(anion_pair[0], cation, anion_pair[1])
                    site_angles.append(angle)

                character = [ang == 90 for ang in site_angles]
                if all(character):
                    isomer_conformations.append('fac')
                else:
                    isomer_conformations.append('mer')
                
            if len(set(isomer_conformations)) == 1:
                isomer = list(set(isomer_conformations))[0]
                f.write("{}\t\t{}\n".format(struct_name, isomer))
                isomer_count[isomer] += 1
                #f.write("{}\t\t{}\n".format(struct_name, list(set(isomer_conformations))[0]))
            elif len(set(isomer_conformations)) > 1:
                f.write("{}\t\tmixed [fac and mer]\n".format(struct_name))
                isomer_count['mixed'] += 1
    
    f.write('\n')
    f.write('{}\n'.format(isomer_count))
    f.write('\n')
    f.close()        
    os.chdir("..")
