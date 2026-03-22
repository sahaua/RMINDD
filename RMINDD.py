''' 
>> Code: RMINDD - (Radiation-Matter Interaction and Damage calculation using Nuclear Data)
>> Perform: Calcuation of metrics of neutron radiation damage in a material using ENDF-6 files
>> Author: Uttiyoarnab Saha
>> Version and Date: 1.0 and 01/07/2022
>> Last modified: 01/07/2022, Kolkata
>> Update: 01/07/2022
>> Major changes: 

=========================================================================================
'''

import numpy, datetime, sys
from time import process_time

## import calculation modules
import UtilsU, ReadU, RecedU, EngdepU, CombinU, TransmU

def printIndexesforReactions (ofile_outRMINDD):
	print("1 = n,g" , file = ofile_outRMINDD)
	print("2 = n,n" , file = ofile_outRMINDD)
	print("3 = n,n'" , file = ofile_outRMINDD)
	print("4 = n,xn" , file = ofile_outRMINDD)
	print("5 = n,particle", file = ofile_outRMINDD)
	print("6 = n,anything", file = ofile_outRMINDD)
	print("7 = total", file = ofile_outRMINDD)
	print('', file = ofile_outRMINDD)

	print('------------------------------------------------', file = ofile_outRMINDD)
	print('', file = ofile_outRMINDD)

start_time = process_time()

day_execution = datetime.date.today()
time_execution = datetime.datetime.now().strftime('%H:%M:%S')

outRMINDD = "Output_RMINDD.txt"
ofile_outRMINDD = open (outRMINDD, 'a')

print('~~~~ RMINDD ~~~~', file = ofile_outRMINDD)
print(day_execution, ' ', time_execution, '\n', file = ofile_outRMINDD)

## command line input for input file name
inpRMINDD = sys.argv[1]

## all data read in from the input file
dict_input_file_variables = ReadU.readCheckInputFile(inpRMINDD, ofile_outRMINDD)

num_of_module_calls = len(dict_input_file_variables)

for i in range(num_of_module_calls):
	## unique energy array and number of energy points
	if (dict_input_file_variables[mod_num]['module_name'] in ['EngdepU', 'RecedU', 'TransmU']):
		(NPt, Etu) = UtilsU.uqce(ofile_outRMINDD, ifile_preprocessedENDF6)

	'''
	The purpose of EngdepU is to compute dpa and heating cross sections due to
	interactions by neutrons in the isotope of a material. It also produces basic
	cross sections and recoil energies. Both point and multigrouped dpa and heating
	can be obtained.
	'''
	if (dict_input_file_variables[mod_num]['module_name'] == "EngdepU"):
		print('~~ RMINDD-EngdepU ~~')
		print('~~ RMINDD-EngdepU ~~', file = ofile_outRMINDD)
		print(':Messages for you:', file = ofile_outRMINDD)
		print('--------------------', file = ofile_outRMINDD)
		print('', file = ofile_outRMINDD)

		printIndexesforReactions(ofile_outRMINDD)

		print('The computed dpa and heating cross sections can', file = ofile_outRMINDD)
		print('be found in files:', file = ofile_outRMINDD)
		print('ndpa--.txt and nheat--.txt and', file = ofile_outRMINDD)
		print('ndpagrouped--.txt and nheatgrouped--.txt' , file = ofile_outRMINDD)
		print('for each reaction', file = ofile_outRMINDD)
		print('', file = ofile_outRMINDD)
		print('Total from CPO reactions is in: ....3001.txt', file = ofile_outRMINDD)
		print('Total from (n, xn) reactions is in: ....1601.txt', file = ofile_outRMINDD)
		print('Total from (n, anything) inexplicit reaction data is in: ....5001.txt', file = ofile_outRMINDD)

		ifile_rawENDF6 = open(dict_input_file_variables[mod_num]['raw_ENDF6_file'], 'r')
		ifile_preprocessedENDF6 = open(dict_input_file_variables[mod_num]['preprocessed_ENDF6_file'], 'r')

		EngdepU.controlAllReactionsHeatingDPA (ofile_outRMINDD, ifile_rawENDF6, ifile_preprocessedENDF6, \
		dict_input_file_variables[mod_num]['input_n_spec'], dict_input_file_variables[mod_num]['num_reac_array'], dict_input_file_variables[mod_num]['num_reac'], \
		NPt, Etu, dict_input_file_variables[mod_num]['atom_displ_model'], dict_input_file_variables[mod_num]['threshold_Ed'], \
		dict_input_file_variables[mod_num]['b_arcdpa'], dict_input_file_variables[mod_num]['c_arcdpa'], \
		dict_input_file_variables[mod_num]['multigroup'], dict_input_file_variables[mod_num]['en_group_type'])

		ifile_preprocessedENDF6.close()
		ifile_rawENDF6.close()

		if (dict_input_file_variables[mod_num]['num_MT_multigroup'] > 0):
			print( 'Multigroup partial reactions .....')
			for i in range (dict_input_file_variables[mod_num]['num_MT_multigroup']):
				print( 'MT = ', dict_input_file_variables[mod_num]['num_MT_group_array'][i])


	'''
	The purpose of RecedU is to produce the PKA spectra induced by reactions of energetic
	neutrons with the isotopes of material. The PKA spectra are produced in a energy 
	group-to-group matrix format. Reaction-wise data, partial sums as well as total PKA
	spectra can be produced.
	'''
	if (dict_input_file_variables[mod_num]['module_name'] == "RecedU"):
		print('~~ RMINDD-RecedU ~~')
		print('~~ RMINDD-RecedU ~~', file = ofile_outRMINDD)
		print(':Messages for you:', file = ofile_outRMINDD)
		print('--------------------', file = ofile_outRMINDD)
		print('', file = ofile_outRMINDD)

		printIndexesforReactions(ofile_outRMINDD)

		print('The computed PKA spectra can be found in files:', file = ofile_outRMINDD)
		print('PKA-MATRICES.txt -- each reaction', file = ofile_outRMINDD)
		print('n-allPKAspectra.txt -- sum total', file = ofile_outRMINDD)

		ifile_rawENDF6 = open(dict_input_file_variables[mod_num]['raw_ENDF6_file'], 'r')
		ifile_preprocessedENDF6 = open(dict_input_file_variables[mod_num]['preprocessed_ENDF6_file'], 'r')

		RecedU.FINE_ENERGY_CALL_REAC (ofile_outRMINDD, ifile_rawENDF6, ifile_preprocessedENDF6, dict_input_file_variables[mod_num]['input_n_spec'], \
		dict_input_file_variables[mod_num]['element_isotope_name'], dict_input_file_variables[mod_num]['en_group_type'], \
		dict_input_file_variables[mod_num]['num_group_limits'], dict_input_file_variables[mod_num]['num_fine_en_points'], \
		dict_input_file_variables[mod_num]['num_reac_array'])

		ifile_preprocessedENDF6.close()
		ifile_rawENDF6.close()

		if (dict_input_file_variables[mod_num]['num_partial_reac_tosum'] > 0):
			ofile1001 = open('n-sum-partialsPKAspectra.txt', 'a')
			print(dict_input_file_variables[mod_num]['element_isotope_name'], file = ofile1001)
			dsdt = numpy.zeros((dict_input_file_variables[mod_num]['num_group_limits']-1, dict_input_file_variables[mod_num]['num_group_limits']-1))
			dsdt = RecedU.ALLSUM (dict_input_file_variables[mod_num]['num_group_limits'], dict_input_file_variables[mod_num]['partial_reac_tosum'])
			print('The sum of recoil nuclei energy spectra for given partial reactions', file = ofile1001)
			for it in range (dict_input_file_variables[mod_num]['num_group_limits']-1):
				print (['{:.6E}'.format(dsdt[it][jt]) for jt in range (dict_input_file_variables[mod_num]['num_group_limits']-1)], file = ofile1001)
			ofile1001.close()

	'''
	The purpose of TransmU is to find the neutron induced gas and transmutation nuclide production cross sections in the
	given isotope of the target element.
	'''
	if (dict_input_file_variables[mod_num]['module_name'] == "TransmU"):
		print('~~ RMINDD-TransmU ~~')
		print('~~ RMINDD-TransmU ~~', file = ofile_outRMINDD)
		print(':Messages for you:', file = ofile_outRMINDD)
		print('--------------------', file = ofile_outRMINDD)
		print('', file = ofile_outRMINDD)

		TransmU.ActivationGasProduction (ofile_outRMINDD, ifile_rawENDF6, ifile_preprocessedENDF6, \
		NPt, Etu, dict_input_file_variables[mod_num])

	''' 
	The purpose of CombinU is to find neutron induced dpa cross sections in the
	multi-element target material. It should be run only after having the required quantities
	for each isotope in the target material calculated using the EngdepU module.
	'''
	if (dict_input_file_variables[mod_num]['module_name'] == "CombinU"):
		print( '~~ RMINDD-CombinU ~~')
		print( '~~ RMINDD-CombinU ~~', file = ofile_outRMINDD)
		print( ':Messages for you:', file = ofile_outRMINDD)
		print('--------------------', file = ofile_outRMINDD)
		print('', file = ofile_outRMINDD)
	
		CombinU.combineXSMultiElementTarget (ofile_outRMINDD, dict_input_file_variables[mod_num])

	## end of looping over as many modules given in input.

stop_time = process_time()
total_time = stop_time - start_time
print('', file = ofile_outRMINDD)
print( 'Total time taken:', file = ofile_outRMINDD)
print(total_time, file = ofile_outRMINDD)

ofile_outRMINDD.close()
