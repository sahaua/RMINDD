''' 
Code: RMINDD - (Radiation Matter Interaction and Damage calculation using Nuclear Data)
Perform: Calcuation of metrics of neutron radiation damage in a material using ENDF-6 files
Module: ReadU.py -- reads and checks provided input data
Author: Uttiyoarnab Saha
Version and Date: 1.0 and 01/07/2022
Last modified: 01/07/2022, Kolkata
Update: 01/07/2022
Major changes: 

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

def compareMATNumbers (ofile_outRMINDD, raw_ENDF6_file, element_isotope_name, MAT_num):
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
	print( 'Isotope provided on input file is:', element_isotope_name, file = ofile_outRMINDD)

	## Check if given mat number matches with ENDF mat number
	print(f'MAT on {raw_ENDF6_file} is {mat}, MAT given in input {MAT_num}', file = ofile_outRMINDD)

	if (mat != MAT_num):
		print( 'Error', file = ofile_outRMINDD)
		print( f'Material does not exist on {raw_ENDF6_file}', file = ofile_outRMINDD)
		print('', file = ofile_outRMINDD)
		raise Exception (f'Input MAT mismatch: Material does not exist on {raw_ENDF6_file}')

'''
Read the input file provided by user. Check the inputs for consistencies.
'''

def readCheckInputFile(inpRMINDD, ofile_outRMINDD):

	dict_input_file_variables = {}

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

		if (words[0] == 'Module_number'):
			module_number = int(words[2])
			dict_input_file_variables[module_number] = {}

		if (words[0] == 'Module_name'):
			dict_input_file_variables[module_number]['module_name'] = words[2]

		## CombinU input cards

		if (words[0] == 'Num_elements_target'):
			dict_input_file_variables[module_number]['num_elements_target'] = int(words[2])

		if (words[0] == 'Elements_target'):
			dict_input_file_variables[module_number]['element_Ed'] = []
			dict_input_file_variables[module_number]['element_Ed_bnd'] = []
			dict_input_file_variables[module_number]['element_stoich'] = [0]*dict_input_file_variables[module_number]['num_elements_target']
			dict_input_file_variables[module_number]['element_recdamen'] = []
			dict_input_file_variables[module_number]['elements_target'] = [0]*dict_input_file_variables[module_number]['num_elements_target']
			for i in range(dict_input_file_variables[module_number]['num_elements_target']):
				dict_input_file_variables[module_number]['elements_target'][i] = words[2+i]
			if (len(dict_input_file_variables[module_number]['elements_target']) != dict_input_file_variables[module_number]['num_elements_target']):
				raise Exception ('Total number of elements in target is not equal to number of elements given!')

		if (words[0] == 'Element_stoichiometries'):
			for i in range(dict_input_file_variables[module_number]['num_elements_target']):
				dict_input_file_variables[module_number]['element_stoich'][i] = float(words[2+i])

		if (words[0] == 'Files_directory'):
			dict_input_file_variables[module_number]['files_dir'] = words[2]

		if (words[0] == 'Element_Ed_target'):
			dict_input_file_variables[module_number]['element_Ed'].append(words[2])
			dict_input_file_variables[module_number][f'Ed_{words[2]}'] = float(words[3])

		if (words[0] == 'Ed_below_no_displacement'):
			dict_input_file_variables[module_number]['element_Ed_bnd'].append(words[2])
			dict_input_file_variables[module_number][f'Ed_bnd_{words[2]}'] = float(words[3])

		if (words[0] == 'Num_isotopes_total'):
			dict_input_file_variables[module_number]['num_isotopes_total'] = int(words[2])
			dict_input_file_variables[module_number]['isotopes_evaluated'] = [0]*dict_input_file_variables[module_number]['num_isotopes_total']

		if (words[0] == 'Isotopes_evaluated'):
			for i in range(dict_input_file_variables[module_number]['num_isotopes_total']):
				dict_input_file_variables[module_number]['isotopes_evaluated'][i] = words[2+i]
			if (len(dict_input_file_variables[module_number]['isotopes_evaluated']) != dict_input_file_variables[module_number]['num_isotopes_total']):
				raise Exception ('Total number of isotopes evaluated is not equal to total number of isotopes!')

		if (words[0] == 'Percent_abundances_all'):
			dict_input_file_variables[module_number]['percent_abundances_all'] = [0]*dict_input_file_variables[module_number]['num_isotopes_total']
			for i in range(dict_input_file_variables[module_number]['num_isotopes_total']):
				dict_input_file_variables[module_number]['percent_abundances_all'][i] = float(words[2+i])/100
			if (len(dict_input_file_variables[module_number]['percent_abundances_all']) != dict_input_file_variables[module_number]['num_isotopes_total']):
				raise Exception ('Number of abundances given is not equal to the total number of isotopes!')

		if (words[0] == 'Recoil_damage_energy_file'):
			dict_input_file_variables[module_number]['element_recdamen'].append(words[2])
			dict_input_file_variables[module_number][f'ifile_Rec_dam_en_{words[2]}'] = words[3]

		if (words[0] == 'Atom_dpaXS_type'):
			dict_input_file_variables[module_number]['atom_dpaXS_type'] = words[2]
			if (dict_input_file_variables[module_number]['atom_dpaXS_type'] == 'Both' or dict_input_file_variables[module_number]['atom_dpaXS_type'] == 'MD-based'):
				dict_input_file_variables[module_number]['element_dameff'] = []

		if (words[0] == 'Damage_efficiency_file'):
			dict_input_file_variables[module_number]['element_dameff'].append(words[2])
			dict_input_file_variables[module_number][f'ifile_Dam_eff_{words[2]}'] = words[3]

		## EngdepU or RecedU or TransmU input cards

		if (words[0] == 'Element_isotope'):
			dict_input_file_variables[module_number]['element_isotope_name'] = words[2]

		if (words[0] == 'Raw_ENDF6_file'):
			dict_input_file_variables[module_number]['raw_ENDF6_file'] = words[2]

		if (words[0] == 'Preprocessed_ENDF6_file'):
			dict_input_file_variables[module_number]['preprocessed_ENDF6_file'] = words[2]

		if (words[0] == 'MAT_num'):
			dict_input_file_variables[module_number]['MAT_num'] = words[2]

		if (words[0] == 'Num_reaction'):
			dict_input_file_variables[module_number]['num_reac'] = int(words[2])
			dict_input_file_variables[module_number]['num_reac_array'] = [0]*dict_input_file_variables[module_number]['num_reac']

		if (words[0] == 'Reaction_num'):
			for i in range(2, dict_input_file_variables[module_number]['num_reac']+2):
				dict_input_file_variables[module_number]['num_reac_array'][i-2] = int(words[i])
			if (dict_input_file_variables[module_number]['num_reac'] > 1):
				for value in dict_input_file_variables[module_number]['num_reac_array']:
					if (value == 7):
						print("Please do only partials or only total. In total all partials will also be done.")
						sys.exit()

		if (words[0] == 'Atom_displ_model'):
			dict_input_file_variables[module_number]['atom_displ_model'] = words[2]

		if (words[0] == 'Threshold_Ed'):
			dict_input_file_variables[module_number]['threshold_Ed'] = float(words[2])

		if (words[0] == 'b_arcdpa'):
			dict_input_file_variables[module_number]['b_arcdpa'] = float(words[2])

		if (words[0] == 'c_arcdpa'):
			dict_input_file_variables[module_number]['c_arcdpa'] = float(words[2])

		if (words[0] == 'Output_filenames'):
			dict_input_file_variables[module_number]['output_filenames'] = words[2]

		if (words[0] == 'Multigroup'):
			dict_input_file_variables[module_number]['multigroup'] = int(words[2])

		if (words[0] == 'Energy_group_type_index'):
			dict_input_file_variables[module_number]['en_group_type'] = int(words[2])

		if (words[0] == 'Input_n_spectrum'):
			dict_input_file_variables[module_number]['input_n_spec'] = int(words[2])

		if (words[0] == 'Num_MT_to_multigroup'):
			dict_input_file_variables[module_number]['num_MT_multigroup'] = int(words[2])
			if (dict_input_file_variables[module_number]['num_MT_multigroup'] > 0):
				dict_input_file_variables[module_number]['num_MT_group_array'] = [0]*dict_input_file_variables[module_number]['num_MT_multigroup']

		if (words[0] == 'MTs_to_multigroup' and dict_input_file_variables[module_number]['num_MT_multigroup'] > 0):
			if (len(words) != dict_input_file_variables[module_number]['num_MT_multigroup'] + 2):
				raise Exception('Number of MTs given for multigrouping is not equal to what is mentioned.')
				sys.exit()
			## store partial MTs to multigroup
			for i in range(2, dict_input_file_variables[module_number]['num_MT_multigroup']+2):
				dict_input_file_variables[module_number]['num_MT_group_array'][i-2] = int(words[i])

		if (words[0] == 'Num_group_limits'):
			dict_input_file_variables[module_number]['num_group_limits'] = int(words[2])

		if (words[0] == 'Num_fine_en_points'):
			dict_input_file_variables[module_number]['num_fine_en_points'] = int(words[2])

		if (words[0] == 'Num_partial_reac_to_sum'):
			dict_input_file_variables[module_number]['num_partial_reac_tosum'] = int(words[2])
			if (dict_input_file_variables[module_number]['num_partial_reac_tosum'] > 0):
				dict_input_file_variables[module_number]['partial_reac_tosum'] = [0]*dict_input_file_variables[module_number]['num_partial_reac_tosum']

		if (words[0] == 'Partial_reac_to_sum' and dict_input_file_variables[module_number]['num_partial_reac_tosum'] > 0):
			if (len(words) != dict_input_file_variables[module_number]['num_partial_reac_tosum'] + 2):
				print('Number of partial reactions given to sum is not equal to what is mentioned.')
				sys.exit()
			for i in range(2, dict_input_file_variables[module_number]['num_partial_reac_tosum'] + 2):
				dict_input_file_variables[module_number]['partial_reac_tosum'][i-2] = int(words[i])

		## TransmU input cards
		if (words[0] == 'Transmgas_group_file'):
			dict_input_file_variables[module_number]['transmgas_group_file'] = words[2]
	
		if (words[0] == 'Transmnucl_group_file'):
			dict_input_file_variables[module_number]['transmnucl_group_file'] = words[2]
	
		if (words[0] == 'Transmgas_point_file'):
			dict_input_file_variables[module_number]['transmgas_point_file'] = words[2]

		if (words[0] == 'Transmnucl_MF5_point_file'):
			dict_input_file_variables[module_number]['transmnucl_MF5_point_file'] = words[2]

		if (words[0] == 'Transmnucl_net_group_file'):
			dict_input_file_variables[module_number]['transmnucl_net_group_file'] = words[2]

	ifile_inpRMINDD.close()

	## checks for different modules

	if (dict_input_file_variables[module_number]['module_name'] == "EngdepU" or dict_input_file_variables[module_number]['module_name'] == "RecedU" or 
	dict_input_file_variables[module_number]['module_name'] == "TransmU"):
		compareMATNumbers (ofile_outRMINDD, dict_input_file_variables[module_number]['raw_ENDF6_file'], \
		dict_input_file_variables[module_number]['element_isotope_name'], dict_input_file_variables[module_number]['MAT_num'])

		# Check if num_reac is within 1 and 7
		if ((dict_input_file_variables[module_number]['module_name'] != "TransmU") and 
		(1 > dict_input_file_variables[module_number]['num_reac'] or dict_input_file_variables[module_number]['num_reac'] > 7)):
			print( 'Error', file = ofile_outRMINDD)
			print( 'wrong reaction index; please follow the list', file = ofile_outRMINDD)
			raise Exception('Indices of reactions to compute must be between 1 and 7')

	if (dict_input_file_variables[module_number]['module_name'] == "CombinU"):
		if (len(dict_input_file_variables[module_number]['element_recdamen']) != dict_input_file_variables[module_number]['num_elements_target']):
			raise Exception ('Recoil and damage energy data for all elements in target are not given!')
		if (len(dict_input_file_variables[module_number]['element_dameff']) != dict_input_file_variables[module_number]['num_elements_target']):
			raise Exception ('Damage efficiency data for all elements in target are not given!')
		if (len(dict_input_file_variables[module_number]['element_Ed']) != dict_input_file_variables[module_number]['num_elements_target']):
			raise Exception ('Threshold lattice displacement energy data for all elements in target are not given!')
		if (len(dict_input_file_variables[module_number]['element_Ed_bnd']) != dict_input_file_variables[module_number]['num_elements_target']):
			raise Exception ('Ed below which no displacement occurs in target for all elements in target are not given!')
		
		for element in dict_input_file_variables[module_number]['element_recdamen']:
			if (element not in dict_input_file_variables[module_number]['elements_target']):
				raise Exception ('Element mismatch between target and damage energy data provided!')
		for element in dict_input_file_variables[module_number]['element_dameff']:
			if (element not in dict_input_file_variables[module_number]['elements_target']):
				raise Exception ('Element mismatch between target and damage efficiency data provided!')
		for element in dict_input_file_variables[module_number]['element_Ed']:
			if (element not in dict_input_file_variables[module_number]['elements_target']):
				raise Exception ('Element mismatch between target and threshold displacement energy data provided!')
		for element in dict_input_file_variables[module_number]['element_Ed_bnd']:
			if (element not in dict_input_file_variables[module_number]['elements_target']):
				raise Exception ('Element mismatch between target and Ed_bnd energy data provided!')


	return(dict_input_file_variables)
