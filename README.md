# RMINDD
RMINDD - Radiation-Matter Interaction and Damage calculation using Nuclear Data
RMINDD is a package of computer programs written in Python programming language. It can read the recent basic evaluated nuclear data libraries such as ENDF/B-VII.1, ENDF/B-VIII.0, TENDL-2023, etc. to estimate the quantities that determine the primary radiation damage in structural materials such as Fe, Ni, Cr, Si, W, etc. used in nuclear applications. The basic features of this program so far developed are distributed here for academic and scientific use. Future developments are possible and will be available here once developed and tested. This README file also be updated as when required.

Please cite the following article when you use this program:
U. Saha, "Primary radiation damage due to neutron interactions using inexplicit evaluated nuclear data: a case study in isotopes of tungsten using ENDF/B-VIII.0 and TENDL-2019", Pramana – J. Phys. (2024) 98:5, https://doi.org/10.1007/s12043-023-02682-2

The input cards / keywords / commands must be given through an input file. The keywords are listed in the distributed file "Input_RMINDD.txt". Any of the modules can be repeatedly run for same or different ENDF/B file for same or different isotope as per requirement as long as it is numbered properly in creasing order in the sequence in which it is run. For example, if run module EngdepU first with Module_number = 1 and again run it after two modules (say, RecedU and TransmU), then before starting with the module specific keywords for EngdepU number it as Module_number = 4, and so on. Given below is a sample with complete list of keywords (available at present).

## Input and data for module EngdepU of RMINDD

Module_number = 1
Module_name = "EngdepU"

Element_isotope = 
Raw_ENDF6_file = 
Preprocessed_ENDF6_file = 
MAT_num = 
Num_reaction = 
Reaction_num = 
Atom_displ_model = 
Threshold_Ed = 
b_arcdpa = 
c_arcdpa =
Output_filenames =  
## multigrouping yes=1 or no=0
Multigroup = 
Energy_group_type_index = 
## input neutron spectrum yes=1 or no=0
## if yes, then there should be Energy group limits in 'Energy-GroupLimits.txt' file and
## the neutron spectrum in 'NeutronSpectrum.txt' file
Input_n_spectrum =  
Num_MT_to_multigroup = 
MTs_to_multigroup = 

## Input and data for module RecedU of RMINDD

Module_number = 2
Module_name = "RecedU"

Element_isotope = 
Raw_ENDF6_file = 
Preprocessed_ENDF6_file = 
MAT_num = 
Num_reaction = 
Reaction_num = 
Num_group_limits = 
Num_fine_en_points = 
Energy_group_type_index = 


## Input and data for module TransmU of RMINDD

Module_number = 3
Module_name = "TransmU"

Element_isotope = 
Raw_ENDF6_file = 
Preprocessed_ENDF6_file = 
MAT_num = 
Transmgas_group_file = 
Transmnucl_group_file = 
Transmgas_point_file = 
Transmnucl_point_file = 


## Input and data for module EngdepU of RMINDD

Module_number = 1
Module_name = "EngdepU"

Element_isotope = 
Raw_ENDF6_file = 
Preprocessed_ENDF6_file = 
MAT_num = 
Num_reaction = 
Reaction_num = 
Atom_displ_model = 
Threshold_Ed = 
b_arcdpa = 
c_arcdpa =
Output_filenames =  
## multigrouping yes=1 or no=0
Multigroup = 
Energy_group_type_index = 
## input neutron spectrum yes=1 or no=0
## if yes, then there should be Energy group limits in 'Energy-GroupLimits.txt' file and
## the neutron spectrum in 'NeutronSpectrum.txt' file
Input_n_spectrum =  
Num_MT_to_multigroup = 
MTs_to_multigroup = 

## Input and data for module CombinU of RMINDD
## CombinU requires output files from EngdepU for its calculations

Module_number = 5
Module_name = "CombinU"

Files_directory = 
Num_elements_target = 
Elements_target = 
## give in the same order as given elements
Element_stoichiometries = 	
Num_isotopes_total = 
## directories with these names will be searched for required files
Isotopes_evaluated = 
## to be given in the same order as the isotopes are given
Percent_abundances_all = 
## better to arrange data files according to the order of name of elements in the target 
Recoil_damage_energy_file = 
Recoil_damage_energy_file = 
Element_Ed_target = 
Element_Ed_target = 
Ed_below_no_displacement = 
Ed_below_no_displacement = 
Atom_dpaXS_type =  		## "MD-based", "NRT-based", "Both"
## there should be as many damage efficiency file as the number of target elements
Damage_efficiency_file = 
Damage_efficiency_file =
