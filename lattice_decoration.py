# -*- coding: utf-8 -*-
"""
Created on Sat Sep 10 12:23:35 2016

@author: nenian
"""

# import numpy and scipy etc
import itertools
import random
import os 
import multiprocessing
import numpy as np
import scipy as sp
from pymatgen import Structure
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer 
from pymatgen.io.vasp import Poscar
from apply_distortions import match_atoms

def count_crystal_system(crystal_system_stats, space_grp, tot):
    if space_grp <= 2:
        crystal_system_stats['Triclinic'] += tot
    elif space_grp <= 15:
        crystal_system_stats['Monoclinic'] += tot
    elif space_grp <= 74:
        crystal_system_stats['Orthorhombic'] += tot
    elif space_grp <= 142:
        crystal_system_stats['Tetragonal'] += tot
    elif space_grp <= 167:
        crystal_system_stats['Trigonal'] += tot
    elif space_grp <= 194:
        crystal_system_stats['Hexagonal'] += tot
    elif space_grp <= 230:
        crystal_system_stats['Cubic'] += tot
    
    return crystal_system_stats

def bravais_collect(collect_structures, total_structs):
    f = open("ORDERING_OUTPUT.txt", 'a')
    crystal_system_stats = {'Triclinic':0,
                                'Monoclinic':0,
                                'Orthorhombic':0,
                                'Tetragonal':0,
                                 'Trigonal':0,
                                'Hexagonal':0,
                                'Cubic':0}
                                
    f.write("Ordering stats by space group.\n")    

    for spgrp in collect_structures:  
        f.write("space_group # {}\t {}\n".format(spgrp,len(collect_structures[spgrp])) )
        crystal_system_stats = count_crystal_system(crystal_system_stats, spgrp, len(collect_structures[spgrp]))
            
    f.write("\n")        
    f.write("The total # of unique structures is {}\n\n".format(total_structs)) 
    f.write("Ordering stats by crystal system.\n")
    for Brav_lat in crystal_system_stats:
        f.write("{}\t {}\n".format(Brav_lat, crystal_system_stats[Brav_lat]))
    
    f.write("-------------------------------------------------------------------------------------\n\n")
    f.close() 
        
#def random_product(*args, **kwds):
#    "Random selection from itertools.product(*args, **kwds)"
#    pools = map(tuple, args) * kwds.get('repeat', 1)
#    return tuple(random.choice(pool) for pool in pools)
def random_product(*args, repeat=1):
    "Random selection from itertools.product(*args, **kwds)"
    pools = [tuple(pool) for pool in args] * repeat
    return tuple(random.choice(pool) for pool in pools)

def random_combination(iterable, r):
    "Random selection from itertools.combinations(iterable, r)"
    pool = tuple(iterable)
    n = len(pool)
    indices = sorted(random.sample(range(n), r))
    return tuple(pool[i] for i in indices)


def compare_structures(strt_list_to_compare):
    matcher = StructureMatcher()
    strt_holder = []
    for grp in matcher.group_structures(strt_list_to_compare):
        strt_holder.append(grp[0])

    return strt_holder

def declutter_list(tup):
    spgrp, structures_list = tup
    decluttered = compare_structures(structures_list)

    return (int(spgrp), decluttered)

def declutter_dict(structure_dict, nprocs):
    
    pool = multiprocessing.Pool(processes=nprocs)
    
    reduced_struct_dict = dict(pool.map(declutter_list, [(space_group,structure_dict[space_group]) for space_group in structure_dict]))
    total_structs = 0
    for spgrp in reduced_struct_dict:
        total_structs += len(reduced_struct_dict[spgrp])
        
    return reduced_struct_dict, total_structs
    

def write_to_file(invars, collect_structures):
    f = open("{}/ORDERING_OUTPUT.txt".format(invars.run_dir), 'a')
    f.write("Writing new structures to disk\n")
    f.write("-------------------------------------------------------------------------------------\n\n")  
    f.close()
    for entry in collect_structures:
        f = open("{}/ORDERING_OUTPUT.txt".format(invars.run_dir), 'a')
        folder_name = "space_group-No-{}".format(entry)
        
        if os.path.isdir('{}/{}'.format(invars.run_dir, folder_name)):
            f.write("{} directory is already present\n".format(folder_name))
            f.write("Going in to look for files\n")
            current_dir = '{}/{}'.format(invars.run_dir, folder_name)
            
        else:
            f.write("Creating {} directory\n".format(folder_name))
            os.mkdir(folder_name)
            current_dir = "{}/{}".format(invars.run_dir, folder_name)
            
        f.close()
        for i, new_strt in enumerate(collect_structures[entry]):
            
            str_name = "{}-str-{}.vasp".format(folder_name, (i+1))
            
            if os.path.isfile("{}/{}".format(current_dir, str_name)):
                f = open("{}/ORDERING_OUTPUT.txt".format(invars.run_dir), 'a')
                f.write("{} file is already present\n".format(str_name))
                f.close()
                pass
            else:
                
                w = Poscar(new_strt)
                w.write_file("{}/{}".format(current_dir, str_name))
                
            if invars.calc_distort == '.FASLE.':
                pass      
            
            elif invars.calc_distort == '.TRUE.' \
            and os.path.isdir('{}/distortions-str-{}'.format(current_dir, (i+1))) is False:
                
                f = open("{}/ORDERING_OUTPUT.txt".format(invars.run_dir), 'a')
                f.write("Evaluating distortions for {}\n".format(str_name))
                f.write("-------------------------------------------------------------------------------------\n")
                f.write('creating distortions-{} directory\n\n'.format(i+1))
                os.mkdir('{}/distortions-str-{}'.format(current_dir,(i+1)))
                 
                dist_dir = '{}/distortions-str-{}'.format(current_dir,(i+1))
                f.close()                                     
                
                #print distortion_files
                match_atoms(new_strt, invars, dist_dir)
                
            elif invars.calc_distort == '.TRUE.' \
            and os.path.isdir('{}/distortions-str-{}'.format(current_dir, (i+1))) is True:
                f = open("{}/ORDERING_OUTPUT.txt".format(invars.run_dir), 'a')
                f.write('distortions-str-{} directory already present\nGoing in...\n\n'.format(i+1))
                
                dist_dir = '{}/distortions-str-{}'.format(current_dir, (i+1))               
                f.close()
                match_atoms(new_strt, invars, dist_dir)


def collect(invars, parent, collect_structures, new_species, count_structures):        
            
    tiling = []
    atom_list = np.array([i.symbol for i in parent.species])
    '''If vacancy ordering is requested...
    Use boolean mask to create vacancies in the cell '''     
    if invars.ion_sub == 'Vac':
        
        mask = np.ones(parent.num_sites, dtype=bool)
        new_species = list(new_species)
        mask[new_species] = False
        tiling = atom_list[mask]
        
        reduced_frac_coords = parent.frac_coords[mask]
        tmp_struct = Structure(parent.lattice.matrix, tiling, reduced_frac_coords)
        
        '''Else, order new atoms by label'''
    else:    
        for i, specie in enumerate(atom_list):
            if i in new_species:
                tiling.append(invars.ion_sub)
            else:
                tiling.append(str(specie))
            
        tmp_struct = Structure(parent.lattice.matrix, tiling, parent.frac_coords)
            
    space_group = int(SpacegroupAnalyzer(tmp_struct, symprec=invars.symprec).get_space_group_number())
            
    # The user may not want to collect structures with P1 symmetry 
    if invars.P1_evaluate == '.FALSE.' and space_group == 1:
        pass
    # Bias the search for unique structures to a few select space groups
    
    elif invars.search_bias != None and space_group not in invars.search_bias:
        pass
    else:
        count_structures += 1
        if space_group in collect_structures:
            
            # Limit bin size to speed up code
            if len(collect_structures[space_group]) <= invars.maxbinsize:
                collect_structures[space_group].append(tmp_struct)
            else:
                count_structures -= 1
        else: 
            collect_structures[space_group] = [tmp_struct]
            
    return collect_structures, count_structures

def anion_combo_cal(invars, parent, i):

    NN = parent.get_neighbors(parent.sites[i], invars.cutoff_radius,\
    include_index=True)
            
    # Calculate how many ions will be replaced
    combo_order = int(invars.sub_ratio*len(NN))            
            
    # We only need the index of the neigh ion            
    neighs_indxs = [neigh[2] for neigh in NN]                
                
    combos = [subset for subset in itertools.combinations(neighs_indxs,combo_order)]
                               
    return combos
        
''' Routine to decorate cation sublattice '''
def cation_decorate(invars, parent):
    os.chdir(invars.run_dir)
    f = open("ORDERING_OUTPUT.txt", 'a')
    f.write("Statistics for cation ordering:\n\n")
    f.write("-------------------------------------------------------------------------------------\n\n")    
    f.close()
    num_cations_to_order = parent.composition.as_dict()[invars.center_atom]
    combo_order = int(invars.sub_ratio*num_cations_to_order) 
    
    collect_structures = dict()
    count_structures = 0
    
    for new_cations in itertools.combinations(parent.indices_from_symbol(invars.center_atom), combo_order):
        
        collect_structures, count_structures = collect(invars, parent, collect_structures, new_cations, count_structures)        
    
    return collect_structures
            
''' Routine to decorate lattice by combination '''
def anions_decorate(invars, parent):
    os.chdir(invars.run_dir)
    f = open("ORDERING_OUTPUT.txt", 'a')
    f.write("Statistics for anion ordering:\n\n")
    f.write("-------------------------------------------------------------------------------------\n\n")
    
    #Find the nearest anion neighbors of the focal cation
    allcombos = []
    if invars.cation_list == None: 
        
        for i, specie in enumerate(parent.species):
            if str(specie.symbol) == invars.center_atom:
                combos = anion_combo_cal(invars, parent, i)
                allcombos.append(combos)
                
    else:
        for i in invars.cation_list:
            
            combos = anion_combo_cal(invars, parent, i)    
            allcombos.append(combos)
            
    total_combinations = len(combos)**parent.composition.get(invars.center_atom)
    f.write("Total combinations to evaluate: {}\n".format(total_combinations))


    # initialize collection   
    collect_structures = dict()
    count_structures = 0    
    
    """Exhaustive search: Slow and wasteful for large systems"""
    if invars.randsearch == '.FALSE.':
        
        for ligand_combos in itertools.product(*allcombos):
            # convert tuple of tuples to list  

            new_ligands = [item  for sublist in ligand_combos for item in sublist]
            
            collect_structures, count_structures = collect(invars, parent, collect_structures, new_ligands, count_structures)            
       
    else:
        """Random search: Much faster. Drastically improves scalability of code"""
        count = 0
        
        # If user does not specify the maximum number of iterations use 
        # total combo length        
        if invars.numiter == None:
            invars.numiter = total_combinations        
            
        while count_structures <= invars.maxrand and count <= invars.numiter:
            count += 1 
            ligand_combos = random_product(*allcombos)
            
            # convert tuple of tuples to list  
            new_ligands = [item  for sublist in ligand_combos for item in sublist]
            
            collect_structures, count_structures = collect(invars, parent, collect_structures, new_ligands, count_structures)
            #if count_structures % 100 == 0:
                #print count_structures


    f.close()
    return collect_structures

''' Routine to decorate lattice by combination for sites manually 
provided by user'''   
def unrestricted_decorate(invars, parent):
    os.chdir(invars.run_dir)
    f = open("ORDERING_OUTPUT.txt", 'a')
    f.write("Statistics for anion ordering:\n\n")
    f.write("-------------------------------------------------------------------------------------\n\n")
    f.close()    
    
    combo_order = int(len(invars.ions_list_to_order)*invars.sub_ratio)
    #print(combo_order)
    collect_structures = dict()
    count_structures = 0
    
    if invars.randsearch == '.FALSE.':    
        for new_ions in itertools.combinations(invars.ions_list_to_order, combo_order):
        
            collect_structures, count_structures = collect(invars, parent, collect_structures, new_ions, count_structures)        
    
        return collect_structures
        
    else:
        count = 0
        
        # If user does not specify the maximum number of iterations use 
        # total combo length        
        if invars.numiter == None:
            invars.numiter = sp.special.comb(len(invars.ions_list_to_order), combo_order)        
            print(invars.numiter) 
        while count_structures <= invars.maxrand and count <= invars.numiter:
            count += 1 
            new_ions = random_combination(invars.ions_list_to_order, combo_order)
            #print(new_ions)
            # convert tuple of tuples to list  
            #new_ligands = [item  for sublist in ligand_combos for item in sublist]
            
            collect_structures, count_structures = collect(invars, parent, collect_structures, new_ions, count_structures)        
            
    return collect_structures
