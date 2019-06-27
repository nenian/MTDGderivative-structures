import os
import numpy as np
from pymatgen import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp import Poscar
import spglib
#from ordering_lattice import compare_structures

''' Rotate cell by 90 degrees to ensure that all possible
    cases of distortion + cell decoration '''

def X_cell_rotation():    
    angle = np.deg2rad(90)
    Rot_X = [[1,0,0],
             [0,np.cos(angle),-np.sin(angle)],
             [0,np.sin(angle), np.cos(angle)]]
    
    return Rot_X
    
def Y_cell_rotation():
    angle = np.deg2rad(90) 
    Rot_Y = [[np.cos(angle), 0, np.sin(angle)],
              [0, 1, 0],
              [-np.sin(angle), 0, np.cos(angle)]]
    
    return Rot_Y

def Z_cell_rotation():
    angle = np.deg2rad(90)    
    Rot_Z = [[np.cos(angle), -np.sin(angle), 0],
              [np.sin(angle), np.cos(angle), 0],
              [0, 0, 1]]    
    
    return Rot_Z


def match_atoms(str1, invars, dist_dir):
    f = open("{}/ORDERING_OUTPUT.txt".format(invars.run_dir), 'a')
    distortion_files = os.listdir(invars.distortion_directory)
    f.write("distortion\t\tSpace group #\n")
    f.write("-------------------------------------------------------------------------------------\n")
    for filename in distortion_files:
        poscar = Poscar.from_file(("{}/{}").format(invars.distortion_directory,filename))
        str2 = poscar.structure
    
        new_coords = []
        for i, posi1 in enumerate(str1.cart_coords):
            distance = 1e100
            coord_holder = []
            bounds = np.array(str1.lattice.abc)
            for j, posi2 in enumerate(str2.cart_coords):
                min_dists = np.min(np.dstack(((np.array(posi1) - np.array(posi2)) % bounds, (np.array(posi2) - np.array(posi1)) % bounds)), axis = 2)
                dists = np.sqrt(np.sum(min_dists ** 2, axis = 1))
               
                if dists < distance:
                   distance = dists
                   coord_holder = str2.frac_coords[j]

            new_coords = np.append(new_coords, coord_holder)
        new_coords = new_coords.reshape(str1.num_sites,3)    
        # Arrange atoms in a list that matches the list of the a0a0a0 structure
        distorted_structure = Structure(str1.lattice, str1.species, new_coords)
        space_group = int(SpacegroupAnalyzer(distorted_structure, invars.symprec).get_space_group_number())
        str_name = "{}-SG-{}.vasp".format(filename, space_group)
        if os.path.isfile(str_name):
            f.write("{}\t\t\t{}\n".format(filename,space_group))          
            
        else:
            f.write("{}\t\t\t{}\n".format(filename,space_group))
            w = Poscar(distorted_structure)
            w.write_file("{}/{}".format(dist_dir,str_name))
            
        #f.write("\n")
    f.write("-------------------------------------------------------------------------------------\n")
    f.write("-------------------------------------------------------------------------------------\n")
    f.close()
    #lattice = Lattice.from_parameters(a=str1.basis_vectors[0][0], b=str1.basis_vectors[1][1], c=str1.basis_vectors[2][2], alpha=90,beta=90,gamma=90)
    
    #rot_X_new_coords = np.dot(new_coords, X_cell_rotation())
    #print "New Coords original"
    #print new_coords
    #print "rotated in X new Coords"
    #print rot_X_new_coords
    #all_permut_of_str.append(Structure(lattice, str1.species, rot_X_new_coords))
    
    #all_permut_of_str.append(structure_x_rot)
    
    #rot_Y_new_coords = np.dot(new_coords, Y_cell_rotation())
    #all_permut_of_str.append(Structure(lattice, str1.species, rot_Y_new_coords))
    
    #rot_Z_new_coords = np.dot(new_coords, Z_cell_rotation())
    
    #all_permut_of_str.append(Structure(lattice, str1.species, rot_Z_new_coords))
    
    #print len(all_permut_of_str) 
    
    #all_permut_of_str = compare_structures(all_permut_of_str)
    
    '''
    for i, structure in enumerate(all_permut_of_str):
        space_group = SpacegroupAnalyzer(structure).get_spacegroup_number()
        #print("{} rotation has str {} space group {}".format(filename, i, space_group))
        str_name = "{}-str-{}-SG{}.vasp".format(filename, i, space_group)
    
        w = Poscar(structure)
        w.write_file(str_name)
    return new_coords
    ''' 
