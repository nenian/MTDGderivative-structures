# -*- coding: utf-8 -*-
"""
Created on Fri Sep  9 10:24:10 2016

@author: nenian
"""

import configparser
from io import StringIO
import os, datetime

class Conf:

                
    def __init__(self, filename):

        def errormessage(error):
            
            if error == 1:
                message = "\t\t\t     ___ _ __ _ __ ___  _ __ \n\
                            / _ \ '__| '__/ _ \| '__|\n\
                            |  __/ |  | | | (_) | |   \n\
                             \___|_|  |_|  \___/|_|\n\n\
                Calculation mode invalid or absent.\n\
                Current modes are:\n\
                0 = apply distortions only\n\
                1 = cation substitution\n\
                2 = anion substitution \n\n\
                Note, if you are trying to do post-process plotting\n\
                plese set plotting = .TRUE. in the input file."
                
                return message
                
            elif error == 2:
                message = "\t\t\t     ___ _ __ _ __ ___  _ __ \n\
                            / _ \ '__| '__/ _ \| '__|\n\
                            |  __/ |  | | | (_) | |   \n\
                             \___|_|  |_|  \___/|_|\n\n\
                \nPlease specify check your input.\n\
                For calculation mode > 0 you must enter\n\
                focal cation, list of anions and cations,\n\
                The ion to be subbed in, rato of ions after subs \n\n"
                
                return message
                
            elif error == 3:
                message= "\t\t\t     ___ _ __ _ __ ___  _ __ \n\
                            / _ \ '__| '__/ _ \| '__|\n\
                            |  __/ |  | | | (_) | |   \n\
                             \___|_|  |_|  \___/|_|\n\n\
                             \nPlease provide the full path to the directory where the distortion casts are stored\n\n"
                
                return message
            
            elif error == 4:
                message = "\t\t\t     ___ _ __ _ __ ___  _ __ \n\
                            / _ \ '__| '__/ _ \| '__|\n\
                            |  __/ |  | | | (_) | |   \n\
                             \___|_|  |_|  \___/|_|\n\n\
                             Log file not found!! Please check your work directory.\n\n"
                return message
        
        def sort_number_lists(str_list):
            num_list = []
            for x in map(str.strip,str_list.split(',')):
                
                if x.isdigit():
                    num_list.append(int(x))
                    continue
                if x[0] == '<':
                    num_list.extend(range(1,int(x[1:])+1))
                    continue
                if x[0] == '>':
                    num_list.extend(range(int(x[1:]), 231))
                    continue
                if '-' in x:
                    xr = list(map(str.strip,x.split('-')))
                    num_list.extend(range(int(xr[0]),int(xr[1])+1))
                    continue 
                else:
                    errormessage(2)
            return num_list
            
        parser = configparser.ConfigParser()

        with open(filename) as stream:
            stream = StringIO("[top]\n" + stream.read())  #Trick ConfigParser 
                                                          #with fake section header
            parser.readfp(stream)
       
        #Set defaults 
        # Indicate if P1 structures will be considered
        # for systems with many combinations exculding P1
        # speeds up the code considerably
        # Default = FALSE
        # Where should the code run/store files?
        try:
            self.run_dir = parser.get('top', 'SAVEIN')
        except configparser.NoOptionError:
            self.run_dir = os.getcwd()
            
        '''Switch over to work dir and open log file for calculation '''
        os.chdir(self.run_dir)

        '''Check to see if log file/ output is present'''
        if os.path.isfile("ORDERING_OUTPUT.txt"):
            log_file = ".TRUE."
        else:
            log_file = ".FALSE."
        #Restart mode is false by default
        try:
            self.restart = parser.get('top', 'restart')
            if self.restart == ".TRUE." and log_file == ".TRUE.":
            
                f = open('ORDERING_OUTPUT.txt', 'a')
                f.write("\n")
                f.write("Ordering code executed on {:%b %d, %Y at %H:%M:%S}\n".format(datetime.datetime.now()))
                f.write("-------------------------------------------------------------------------------------\n\n")
                f.write("RESTART: {}\n".format(self.restart))
            elif self.restart == ".TRUE." and log_file == ".FALSE.":
                f = open('ORDERING_OUTPUT.txt', 'w')
                f.write("Ordering code executed on {:%b %d, %Y at %H:%M:%S}\n".format(datetime.datetime.now()))
                f.write("-------------------------------------------------------------------------------------\n\n")
                f.write("Input variables for this run:\n\n")
                f.write("-------------------------------------------------------------------------------------\n\n")
                error = errormessage(4)
                f.write("{}\n".format(error))
                self.calc_mode = "ERROR"
            
        except configparser.NoOptionError:
            self.restart = ".FALSE."
            f = open('ORDERING_OUTPUT.txt', 'w')
            f.write("Ordering code executed on {:%b %d, %Y at %H:%M:%S}\n".format(datetime.datetime.now()))
            f.write("-------------------------------------------------------------------------------------\n\n")
            f.write("Input variables for this run:\n\n")
            f.write("-------------------------------------------------------------------------------------\n\n")
            f.write("RESTART: {}\n".format(self.restart))
        


        #Plotting mode is false by default           
        try:
            self.plotting = parser.get('top', 'plotting')
        except configparser.NoOptionError:
            self.plotting = ".FALSE."
        
        try:
            self.symprec = float(parser.get('top', 'symprec'))
        except configparser.NoOptionError:
            self.symprec = 0.001
        
        try:
            self.calc_mode = int(parser.get('top','calc_mode'))
            if self.calc_mode == 1:
                self.calc_mode_descrption = "Cation ordering"
            elif self.calc_mode == 2:
                self.calc_mode_descrption = "Anion ordering"
            elif self.calc_mode == 3:
                self.calc_mode_descrption = "Unrestricted"
            
        except configparser.NoOptionError:
            if self.restart == '.TRUE.' and  self.plotting == '.TRUE.':
                print("in the right place")
                self.calc_mode = "post-processing"
                self.calc_mode_descrption = "Plotting"
            else:
                print("in the wrong place")
                error = errormessage(1)
                f.write("{}".format(error))
                self.calc_mode = "ERROR"
            #sys.exit()            
        f.write("Calculation mode: {} -- {}\n".format(self.calc_mode, self.calc_mode_descrption))    
        
        try:
            self.P1_evaluate = parser.get('top', 'P1_EVAL')
        except configparser.NoOptionError:
            self.P1_evaluate = '.FALSE.'
        
        
        try:
            self.randsearch = parser.get('top', 'rand_search')
        except configparser.NoOptionError:
            self.randsearch = '.FALSE.'
            
        try:
            self.maxrand = int(parser.get('top', 'max_rand'))
        except configparser.NoOptionError:
            self.maxrand = 2000
            
        try:
            self.nprocs = int(parser.get('top', 'nprocs'))
        except configparser.NoOptionError:
            self.nprocs = 1
            
        try:
            self.maxbinsize = int(parser.get('top','maxbin'))
        except configparser.NoOptionError:
            self.maxbinsize = 1e6
        
        try:
            self.numiter = int(parser.get('top','numiter'))
        except configparser.NoOptionError:
            self.numiter = None
            
        # Apply distortions after ordering?
        # Default = FALSE
        try:
            self.calc_distort = parser.get('top', 'distortions')
        except configparser.NoOptionError:
            self.calc_distort = '.FALSE.'        
        
        # If the user selects to calculate distortions
        # they must also provide the location of the 
        # folder with the distortion files
        if self.calc_distort == '.TRUE.':
            try:
                self.distortion_directory = parser.get('top', 'dist_der')
            except configparser.NoOptionError:
                error = errormessage(3)
                f.write("{}".format(error))
                self.calc_mode = "ERROR"
        
        try:
            self.calc_mode = int(parser.get('top','calc_mode'))
            
        except configparser.NoOptionError:
            errormessage(1)
        
        try:
            self.cation_list = parser.get('top', 'focalcats').split()
            self.cation_list = [int(x) for x in self.cation_list]
        except configparser.NoOptionError:
            self.cation_list = None
        
        # Search Bias
        try:
            
            bias = parser.get('top', 'search_bias')
            
            self.search_bias = sort_number_lists(bias)
                    
            #Ensure no duplicates
            self.search_bias = list(set(self.search_bias))
        except configparser.NoOptionError:
            self.search_bias = None    
                
        # Ions list for unrestricted search
        try:
            str_of_ions = parser.get('top', 'orderions')
            
            
            self.ions_list_to_order = []
            for x in map(str.strip,str_of_ions.split(',')):
                
                if x.isdigit():
                    self.ions_list_to_order.append(int(x))
                    continue
                if x[0] == '<':
                    self.ions_list_to_order.extend(range(0,int(x[1:])))
                    continue
                if '-' in x:
                    xr = list(map(str.strip,x.split('-')))
                    self.ions_list_to_order.extend(range(int(xr[0]),int(xr[1])+1))
                    continue 
                else:
                    errormessage(2)
                       
            
            self.ions_list_to_order = list(set(self.ions_list_to_order))
            
        except configparser.NoOptionError:
            self.ions_list_to_order = None
            
        #print(self.ions_list_to_order)
        
        if self.calc_mode == 1 or self.calc_mode == 2 or self.calc_mode == 3:
            try:
                self.center_atom=parser.get('top', 'cenat')
                self.cations = parser.get('top','cations').split()
                self.anions =parser.get('top', 'anions').split()
                self.ion_sub=parser.get('top' ,'ionsub')
                
                #sub_ratio can be entered as a float of a fraction x/y
                try: 
                    self.sub_ratio=float(parser.get('top', 'subratio'))
                except ValueError:
                    frac = parser.get('top', 'subratio')
                    x = frac.split('/')
                    self.sub_ratio = float(x[0])/float(x[1])
                    
            except configparser.NoOptionError:
                
                error = errormessage(2)
                f.write("{}".format(error))
                self.calc_mode = "ERROR"
                
        elif self.calc_mode < 0 or self.calc_mode is not "post-processing":
            error = errormessage(1)
            f.write("{}".format(error))
            self.calc_mode = "ERROR"
        else:
            self.center_atom=None
                
            self.cations = None
            self.anions = None
            self.ion_sub = None
            self.sub_ratio = None
        
        if self.calc_mode == 2:
            try:
                self.cutoff_radius=float(parser.get('top', 'cutoff'))
            except configparser.NoOptionError:
                f.write("\nFor anion substitution you must provide a cutoff\n\
                radius in order to compute the nearest neighbor atoms.\n\
                Using default cutoff_radius = 3.2 Angstroms\n\n")
                self.cutoff_radius= 3.2
                
        else:
            self.cutoff_radius= None
                #errormessage(3)
            
        '''Write Input data to log file'''
        if self.calc_mode == "ERROR" or self.restart == ".TRUE.":
            pass
        else:
            f.write("Cations: {}\n".format(self.cations))
            f.write("Anions: {}\n".format(self.anions))
            f.write("Focal cation: {}\n".format(self.center_atom))  
            if self.cation_list is not None:
                f.write("User specified focal cation {}\n".format(self.cation_list))
            f.write("Substituted ion: {}\n".format(self.ion_sub))
            f.write("Substitution ratio: {}\n".format(self.sub_ratio))
            f.write("Radius for nearest neighbor calculation: {}\n".format(self.cutoff_radius))
            f.write("Evaluate structures with P1 symmetry: {}\n".format(self.P1_evaluate))
            f.write("Apply distortions: {}\n".format(self.calc_distort))
            if self.calc_distort == ".TRUE.":
                f.write("Distortion 'casts' located in {}\n".format(self.distortion_directory))
            f.write("-------------------------------------------------------------------------------------\n\n")
        f.close()
                
