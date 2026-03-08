''' 
>> Code: RMINDD - (Radiation-Matter Interaction and Damage calculation using Nuclear Data)
>> Perform: Calcuation of metrics of neutron radiation damage in a material using ENDF-6 files
>> Module: ReadU.py -- reads and checks provided input data
>> Author: Uttiyoarnab Saha
>> Version and Date: 1.0 and 01/07/2022
>> Last modified: 01/07/2022, Kolkata
>> Update: 01/07/2022
>> Major changes: 

=========================================================================================
'''

import sys

'''
Reads the MAT number and name of isotope present on the evaluation.
Prints MAT number and isotope name provided and those present on the
evaluation to putput file.
Checks if MAT numbers are the same, if not, raises an exception.
Called when EngdepU, RecedU and TransmU modules are called for computations.
'''

def compareMATNumbers (raw_ENDF6_file, MAT_num):
	ifile_rawENDF = open (raw_ENDF6_file, 'r')
	ifile_rawENDF.readline()
	data_rawENDF = ifile_rawENDF.readline().split()
	#mat = data_rawENDF[6]
	if (len(data_rawENDF) == 7):
		mat = data_rawENDF[5][-4:]
	if (len(data_rawENDF) == 8):
		mat = data_rawENDF[6]	# mat = data_rawENDF[5][-4:]

	for i in range(3):
		ifile_rawENDF.readline()
	data_rawENDF = ifile_rawENDF.readline().split()
	iso = data_rawENDF[0] + data_rawENDF[1]
	ifile_rawENDF.close()

	print( f'Evaluation on {raw_ENDF6_file} is: ', iso, file = ofile_outRMINDD)
	print('', file = ofile_outRMINDD)
	print( 'Isotope provided on input file is:' element_isotope_name, file = ofile_outRMINDD)

	## Check if given mat number matches with ENDF mat number
	print(f'MAT on {raw_ENDF6_file} is {mat}, MAT given in input {MAT_num}', file = ofile_outRMINDD)

	if (mat != MAT_num):
		print( 'Error', file = ofile_outRMINDD)
		print( 'Material does not exist on tape01', file = ofile_outRMINDD)
		print('', file = ofile_outRMINDD)
		raise Exception ('Input MAT mismatch: Material does not exist on tape01')

'''
Read the input file provided by user. Check the inputs for consistencies.
'''

def readCheckInputFile(inpRMINDD, ofile_outRMINDD):
	ifile_inpRMINDD = open(inpRMINDD, 'r')
	ifile_inpRMINDD.seek(0, 2)
	eof = ifile_inpRMINDD.tell()
	ifile_inpRMINDD.seek(0, 0)

	while (ifile_inpRMINDD.tell() != eof):
		line = ifile_inpRMINDD.readline()
		# blank line
		if (line.strip(' ') == '\n'):
			continue
		data1 = line.split()
		words = []

		for i in range(len(data1)):
			data1[i] = data1[i].strip(' ')
			data1[i] = data1[i].strip(',')
			data1[i] = data1[i].strip('"')
			data1[i] = data1[i].strip("'")
			words.append(data1[i])

		## comment line
		if (words[0] == '#'):
			continue

		if (words[0] == 'Module_name'):
			global module_name = words[2]

		## CombinU input cards

		if (words[0] == 'Num_elements_target'):
			num_elements_target = int(words[2])

		if (words[0] == 'Elements_target'):
			element_Ed = []
			element_Ed_bnd = []
			global element_stoich = [0]*num_elements_target
			global element_recdamen = []
			global elements_target = [0]*num_elements_target
			for i in range(num_elements_target):
				elements_target[i] = words[2+i]
			if (len(elements_target) != num_elements_target):
				raise Exception ('Total number of elements in target is not equal to number of elements given!')

		if (words[0] == 'Element_stoichiometries'):
			for i in range(num_elements_target):
				element_stoich[i] = float(words[2+i])

		if (words[0] == 'Files_directory'):
			global files_dir = words[3]

		if (words[0] == 'Element_Ed_target'):
			element_Ed.append(words[2])
			globals()[f'Ed_{words[2]}'] = float(words[3])

		if (words[0] == 'Ed_below_no_displacement'):
			element_Ed_bnd.append(words[2])
			globals()[f'Ed_bnd_{words[2]}'] = float(words[3])

		if (words[0] == 'Num_isotopes_total'):
			num_isotopes_total = int(words[2])
			global isotopes_evaluated = [0]*num_isotopes_total

		if (words[0] == 'Isotopes_evaluated'):
			for i in range(num_isotopes_total):
				isotopes_evaluated[i] = words[2+i]
			if (len(isotopes_evaluated) != num_isotopes_total):
				raise Exception ('Total number of isotopes evaluated is not equal to total number of isotopes!')

		if (words[0] == 'Percent_abundances_all'):
			global percent_abundances_all = [0]*num_isotopes_total
			for i in range(num_isotopes_total):
				percent_abundances_all[i] = float(words[2+i])/100
			if (len(percent_abundances_all) != num_isotopes_total):
				raise Exception ('Number of abundances given is not equal to the total number of isotopes!')

		if (words[0] == 'Recoil_damage_energy_file'):
			element_recdamen.append(words[2])
			globals()[f'ifile_Rec_dam_en_{words[2]}'] = words[3]

		if (words[0] == 'Atom_dpaXS_type'):
			atom_dpaXS_type = words[2]
			if (atom_dpaXS_type == 'Both' or atom_dpaXS_type == 'MD-based'):
				global element_dameff = []

		if (words[0] == 'Damage_efficiency_file'):
			element_dameff.append(words[2])
			globals()[f'ifile_Dam_eff_{words[2]}'] = words[3]

		## EngdepU (or RecedU) input cards

		if (words[0] == 'Element_isotope'):
			element_isotope_name = words[2]

		if (words[0] == 'Raw_ENDF6_file'):
			global raw_ENDF6_file = words[2]

		if (words[0] == 'Preprocessed_ENDF6_file'):
			global preprocessed_ENDF6_file = words[2]

		if (words[0] == 'MAT_num'):
			MAT_num = int(words[2])

		if (words[0] == 'Num_reaction'):
			global num_reac = int(words[2])
			global num_reac_array = [0]*num_reac

		if (words[0] == 'Reaction_num'):
			for i in range(2, num_reac+2):
				num_reac_array[i-2] = int(words[i])
			if (num_reac > 1):
				for value in num_reac_array:
					if (value == 7):
						print("Please do only partials or only total. In total all partials will also be done.")
						sys.exit()

		if (words[0] == 'Atom_displ_model'):
			global atom_displ_model = words[2]

		if (words[0] == 'Threshold_Ed'):
			global threshold_Ed = float(words[2])

		if (words[0] == 'b_arcdpa'):
			global b_arcdpa = float(words[2])

		if (words[0] == 'c_arcdpa'):
			global c_arcdpa = float(words[2])

		if (words[0] == 'Multigroup'):
			global multigroup = words[2]

		if (words[0] == 'Energy_group_type_index'):
			global en_group_type = int(words[2])

		if (words[0] == 'Input_n_spectrum'):
			global input_n_spec = int(words[2])

		if (words[0] == 'Num_MT_to_multigroup'):
			global num_MT_multigroup = int(words[2])
			if (num_MT_multigroup > 0):
				global num_MT_group_array = [0]*num_MT_multigroup

		if (words[0] == 'MTs_to_multigroup' and num_MT_multigroup > 0):
			if (len(words) != num_MT_multigroup + 2):
				raise Exception('Number of MTs given for multigrouping is not equal to what is mentioned.')
				sys.exit()
			## store partial MTs to multigroup
			for i in range(2, num_MT_multigroup+2):
				num_MT_group_array[i-2] = int(words[i])

		if (words[0] == 'Num_group_limits'):
			global num_group_limits = int(words[2])

		if (words[0] == 'Num_fine_en_points'):
			global num_fine_en_points = int(words[2])

		if (words[0] == 'Num_partial_reac_to_sum'):
			global num_partial_reac_tosum = int(words[2])
			if (num_partial_reac_tosum > 0):
				global partial_reac_tosum = [0]*num_partial_reac_tosum

		if (words[0] == 'Partial_reac_to_sum' and num_partial_reac_tosum > 0):
			if (len(words) != num_partial_reac_tosum + 2):
				print('Number of partial reactions given to sum is not equal to what is mentioned.')
				sys.exit()
			for i in range(2, num_partial_reac_tosum + 2):
				partial_reac_tosum[i-2] = int(words[i])

	ifile_inpRMINDD.close()

	## checks for different modules

	if (module_name == "EngdepU" or module_name == "RecedU"):
		compareMATNumbers (raw_ENDF6_file, MAT_num)

		# Check if num_reac is within 1 and 7
		if (1 > num_reac or num_reac > 7):
			print( 'Error', file = ofile_outRMINDD)
			print( 'wrong reaction index; please follow the list', file = ofile_outRMINDD)
			raise Exception('Indices of reactions to compute must be between 1 and 7')
		
	if (module_name == "CombinU"):
		if (len(element_recdamen) != num_elements_target):
			raise Exception ('Recoil and damage energy data for all elements in target are not given!')
		if (len(element_dameff) != num_elements_target):
			raise Exception ('Damage efficiency data for all elements in target are not given!')
		if (len(element_Ed) != num_elements_target):
			raise Exception ('Threshold lattice displacement energy data for all elements in target are not given!')
		if (len(element_Ed_bnd) != num_elements_target):
			raise Exception ('Ed below which no displacement occurs in target for all elements in target are not given!')
	
		for element in element_recdamen:
			if (element not in elements_target):
				raise Exception ('Element mismatch between target and damage energy data provided!')
		for element in element_dameff:
			if (element not in elements_target):
				raise Exception ('Element mismatch between target and damage efficiency data provided!')
		for element in element_Ed:
			if (element not in elements_target):
				raise Exception ('Element mismatch between target and threshold displacement energy data provided!')
		for element in element_Ed_bnd:
			if (element not in elements_target):
				raise Exception ('Element mismatch between target and Ed_bnd energy data provided!')
