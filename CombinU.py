## Module CombinU.py of the RMINDD code
## Uttiyoarnab Saha

## ------------------------------------------------------------------------


import numpy
import re
import glob
import os, sys
import matplotlib.pyplot as plt

## Separate the symbol and mass number from names of directories 

def separateSymbolMassNumber(word):
	if not isinstance(word, str):
		raise ValueError("Input must be a string.")

	# Extract symbols and mass numbers separately
	symbol = ''.join(re.findall(r'[A-Za-z]', word))
	mass_num = ''.join(re.findall(r'\d', word))
	return (symbol, mass_num)

#=======Multigroup=======*

## Multigroup, according to requirement, the point dpa and heating
## cross sections into the chosen neutron energy group structure.

def groupmulti(igtype, NP, E, sdpa):

	if (igtype==0):
		ifile = open('Energy-GroupLimits.txt', 'r')
		Ngl = int(ifile.readline().split()[-1])
		Eg = numpy.zeros(Ngl); gsdpa = numpy.zeros(Ngl)
		for i in reversed(range(Ngl)):
			Eg[i] = float(ifile.readline().split()[0])
		ifile.close()

	if (igtype==7):
		Ngl = 101
		nre = Ngl-1
		Eg = numpy.zeros(Ngl); gsdpa = numpy.zeros(Ngl)
		Eg = engrp7()

	if (igtype==8):
		Ngl = 101
		nre = Ngl-1
		Eg = numpy.zeros(Ngl); gsdpa = numpy.zeros(Ngl)
		Eg = engrp8()

	fi = numpy.ones(Ngl)

	ifg = 0
	for i in range (Ngl-1):
		if (Eg[i]<=E[0] and E[1]<=Eg[i+1]):
			ifg = i
			break

	for i in range (ifg, Ngl-1):
		Eg1 = Eg[i]
		Eg2 = Eg[i+1]
		Nsect = 10000
		h = (Eg2-Eg1)/Nsect
		bcs = 0; dhval = 0; dhcs = 0; denominator = 0; deno1 = 0; deno2 = 0
		bcs1 = 0; bcs2 = 0; dhval1 = 0; dhval2 = 0; dhcs1 = 0; dhcs2 = 0
		for j in range (Nsect):
			t = Eg1 + j*h

			flux1 = numpy.interp(t, Eg, fi)	#srchintrp3 (Eg,fi,Ngl,t)
			dhcs1 = numpy.interp(t+h, E, sdpa) * flux1
			deno1 = flux1

			flux2 = numpy.interp(t+h, Eg, fi)	#srchintrp3 (Eg,fi,Ngl,t+h)
			dhcs2 = numpy.interp(t+h, E, sdpa) * flux2
			deno2 = flux2

			dhcs = dhcs + (h/2)*(dhcs1 + dhcs2)
			denominator = denominator + (h/2)*(deno1 + deno2)

		if (dhcs != 0 and denominator != 0):
			gsdpa[i] = dhcs/denominator

	return(Ngl,Eg,gsdpa)

def engrp7():
	Ngl = 101
	Eg = [1.00E-04, 1.00E-03, 1.00E-02, 2.30E-02, 5.00E-02, 7.60E-02,
	1.10E-01, 1.70E-01, 2.50E-01, 3.80E-01, 5.50E-01, 8.40E-01, 1.28E+00,
	1.90E+00, 2.80E+00, 4.20E+00, 6.30E+00, 9.20E+00, 1.30E+01, 2.10E+01,
	3.00E+01, 4.50E+01, 6.90E+01, 1.00E+02, 1.30E+02, 1.70E+02, 2.20E+02,
	2.80E+02, 3.60E+02, 4.50E+02, 5.70E+02, 7.60E+02, 9.60E+02, 1.28E+03,
	1.60E+03, 2.00E+03, 2.70E+03, 3.40E+03, 4.50E+03, 5.50E+03, 7.20E+03,
	9.20E+03, 1.20E+04, 1.50E+04, 1.90E+04, 2.50E+04, 3.20E+04, 4.00E+04,
	5.20E+04, 6.60E+04, 8.80E+04, 1.10E+05, 1.30E+05, 1.60E+05, 1.90E+05,
	2.20E+05, 2.50E+05, 2.90E+05, 3.20E+05, 3.60E+05, 4.00E+05, 4.50E+05,
	5.00E+05, 5.50E+05, 6.00E+05, 6.60E+05, 7.20E+05, 7.80E+05, 8.40E+05,
	9.20E+05, 1.00E+06, 1.20E+06, 1.40E+06, 1.60E+06, 1.80E+06, 2.00E+06,
	2.30E+06, 2.60E+06, 2.90E+06, 3.30E+06, 3.70E+06, 4.10E+06, 4.50E+06,
	5.00E+06, 5.50E+06, 6.00E+06, 6.70E+06, 7.40E+06, 8.20E+06, 9.00E+06,
	1.00E+07, 1.10E+07, 1.20E+07, 1.30E+07, 1.40E+07, 1.50E+07, 1.60E+07,
	1.70E+07, 1.80E+07, 1.90E+07, 2.00E+07]
	
	return(Eg)

def engrp8():
	Ngl = 101
	Eg = [1.00E-05,1.00E-03,1.00E-02,2.30E-02,5.00E-02,7.60E-02,
	1.10E-01,1.70E-01,2.50E-01,3.80E-01,5.50E-01,8.40E-01,
	1.28E+00,1.90E+00,2.80E+00,4.20E+00,6.30E+00,9.20E+00,
	1.30E+01,2.10E+01,3.00E+01,4.50E+01,6.90E+01,1.00E+02,
	1.30E+02,1.70E+02,2.20E+02,2.80E+02,3.60E+02,4.50E+02,
	5.70E+02,7.60E+02,9.60E+02,1.28E+03,1.60E+03,2.00E+03,
	2.70E+03,3.40E+03,4.50E+03,5.50E+03,7.20E+03,9.20E+03,
	1.20E+04,1.50E+04,1.90E+04,2.50E+04,3.20E+04,4.00E+04,
	5.20E+04,6.60E+04,8.80E+04,1.10E+05,1.30E+05,1.60E+05,
	1.90E+05,2.20E+05,2.50E+05,2.90E+05,3.20E+05,3.60E+05,
	4.00E+05,4.50E+05,5.00E+05,5.50E+05,6.00E+05,6.60E+05,
	7.20E+05,7.80E+05,8.40E+05,9.20E+05,1.00E+06,1.20E+06,
	1.40E+06,1.60E+06,1.80E+06,2.00E+06,2.30E+06,2.60E+06,
	2.90E+06,3.30E+06,3.70E+06,4.10E+06,4.50E+06,5.00E+06,
	5.50E+06,6.00E+06,6.70E+06,7.40E+06,8.20E+06,9.00E+06,
	1.00E+07,1.10E+07,1.20E+07,1.30E+07,1.40E+07,1.50E+07,
	1.60E+07,1.70E+07,1.80E+07,1.90E+07,2.00E+07]

	return(Eg)


def readInputFile(ifilename):
	ifile = open(ifilename, 'r')
	ifile.seek(0, 2)			# go to the end of file
	eof = ifile.tell()			# get the end-of-file position
	ifile.seek(0, 0)			# go to the start of file
	
	while (ifile.tell() != eof):
		line = ifile.readline()
		# blank line
		if (line.strip(' ') == '\n'):
			continue
		data1 = line.split()
		words = []

		for i in range(len(data1)):
			data1[i] = data1[i].strip(' ')
			data1[i] = data1[i].strip(',')
			data1[i] = data1[i].strip('"')
			words.append(data1[i])

		# comment line
		if (words[0] == '#'):
			continue

		if (words[0] == 'Num_elements_target'):
			num_elements_target = int(words[2])

		if (words[0] == 'Elements_target'):
			element_Ed = []
			element_Ed_bnd = []
			element_stoich = [0]*num_elements_target
			element_recdamen = []
			elements_target = [0]*num_elements_target
			for i in range(num_elements_target):
				elements_target[i] = words[2+i]
			if (len(elements_target) != num_elements_target):
				raise Exception ('Total number of elements in target is not equal to number of elements given!')

		if (words[0] == 'Element_stoichiometries'):
			for i in range(num_elements_target):
				element_stoich[i] = float(words[2+i])

		if (words[0] == 'Element_Ed_target'):
			element_Ed.append(words[2])
			globals()[f'Ed_{words[2]}'] = float(words[3])

		if (words[0] == 'Ed_below_no_displacement'):
			element_Ed_bnd.append(words[2])
			globals()[f'Ed_bnd_{words[2]}'] = float(words[3])

		if (words[0] == 'Num_isotopes_total'):
			num_isotopes_total = int(words[2])
			isotopes_evaluated = [0]*num_isotopes_total

		if (words[0] == 'Isotopes_evaluated'):
			for i in range(num_isotopes_total):
				isotopes_evaluated[i] = words[2+i]
			if (len(isotopes_evaluated) != num_isotopes_total):
				raise Exception ('Total number of isotopes evaluated is not equal to total number of isotopes!')

		if (words[0] == 'Percent_abundances_all'):
			percent_abundances_all = [0]*num_isotopes_total
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
				element_dameff = []

		if (words[0] == 'Damage_efficiency_file'):
			element_dameff.append(words[2])
			globals()[f'ifile_Dam_eff_{words[2]}'] = words[3]

	ifile.close()
	
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

	return (elements_target, isotopes_evaluated, element_recdamen, element_dameff, percent_abundances_all, element_stoich)

lib = sys.argv[1]
ifilename = sys.argv[2]
(elements_target, isotopes_evaluated, element_recdamen, element_dameff, \
percent_abundances_all, element_stoich) = readInputFile(ifilename)

## have to test for some universal value for threshold lattice displacement energy in the target
## testing for target SiC -- according to the results of tests, may have to change in inputs that are taken! 
checkvalue_Ed = 2.0/(1/35 + 1/20)

## Read and store damage energy and damage efficiency data

for element in element_dameff:
	ifile = open (globals()[f'ifile_Dam_eff_{element}'], 'r')
	ifile.readline()
	ifile.readline()
	ifile.readline()

	globals()[f'num_data_{element}_ref'] = int(ifile.readline().split()[0])
	globals()[f'T_dam_{element}_ref'] = numpy.zeros(globals()[f'num_data_{element}_ref'])
	globals()[f'Dam_eff_{element}_ref'] = numpy.zeros(globals()[f'num_data_{element}_ref'])

	for i in range(globals()[f'num_data_{element}_ref']):
		data = ifile.readline().split()
		globals()[f'T_dam_{element}_ref'][i] = float(data[0]) * 1000.0
		globals()[f'Dam_eff_{element}_ref'][i] = float(data[1])
	ifile.close()


## Find damage energy from SRIM-2013 for the target
## Read and store data for recoil energy, damage energy and number of vacancies
## obtained from SRIM-2013

for element in element_recdamen:
	ifile = open(globals()[f'ifile_Rec_dam_en_{element}'], 'r')
	ifile.readline()
	line = ifile.readline()
	globals()[f'num_data_{element}_srim'] = int(line.split()[1])
	globals()[f'E_R_{element}_srim'] = []; globals()[f'T_dam_{element}_srim'] = []; globals()[f'num_vac_{element}_srim'] = []

	for i in range(globals()[f'num_data_{element}_srim']):
		line = ifile.readline()
		data = line.split()
		globals()[f'E_R_{element}_srim'].append(float(data[0]))
		globals()[f'T_dam_{element}_srim'].append(float(data[-2]))
		globals()[f'num_vac_{element}_srim'].append(float(data[-1]))
	ifile.close()


## Using reaction-wise isotopic cross sections and recoil energies
## from evaluated nuclear data to find damage energies and number 
## of displacements using the damage efficiency data, whcih are
## computed using the EngdepU module.

## The data are expected to be in directories named according to the isotopes of elements,
## i.e., as per the input isotopes_evaluated in the input file such as Si28, C12, Ni58, etc.

## First extracting the required data into containers

## Total n-interaction cross sections in isotopes

os.chdir(lib)
for isotope in isotopes_evaluated:
	os.chdir(isotope)
	ifile = open('E8.1.nheat1.txt')
	globals()[f'num{isotope}'] = int(ifile.readline().split()[0])
	globals()[f'En{isotope}_MT1'] = numpy.zeros(globals()[f'num{isotope}'])
	globals()[f'XS{isotope}_MT1'] = numpy.zeros(globals()[f'num{isotope}'])
	
	for i in range(globals()[f'num{isotope}']):
		data = ifile.readline().split()
		globals()[f'En{isotope}_MT1'][i] = float(data[0])
		globals()[f'XS{isotope}_MT1'][i] = float(data[1])
	ifile.close()
	os.chdir('../')


## Loop through all files in each isotope directory

for isotope in isotopes_evaluated:
	os.chdir(isotope)
	try:
		(element, mass_num) = separateSymbolMassNumber(isotope)
	except ValueError as e:
		print(f"Error: {e}")

	## Match files with a specific pattern
	globals()[f'file_pattern{isotope}'] = "nheat*.txt"
	globals()[f'files{isotope}'] = glob.glob(globals()[f'file_pattern{isotope}'])
	globals()[f'MTreac{isotope}'] = []

	for file in globals()[f'files{isotope}']:
		MTreac = file.split('.')[0][5:]
		with open(file, 'r') as ifile:
			num = int(ifile.readline().split()[0])

			globals()[f'En{isotope}_MT{MTreac}'] = numpy.zeros(num)
			globals()[f"XS1{isotope}_MT{MTreac}"] = numpy.zeros(num)
			globals()[f"Er1{isotope}_MT{MTreac}"] = numpy.zeros(num)
			for i in range(num):
				data = ifile.readline().split()
				globals()[f'En{isotope}_MT{MTreac}'][i] = float(data[0])
				globals()[f'XS1{isotope}_MT{MTreac}'][i] = float(data[1])
				globals()[f'Er1{isotope}_MT{MTreac}'][i] = float(data[2])
		ifile.close()
		globals()[f'MTreac{isotope}'].append(MTreac)

		globals()[f'Tdam{isotope}_MT{MTreac}'] = numpy.zeros(globals()[f'num{isotope}'])
		globals()[f'DamEff{isotope}_MT{MTreac}'] = numpy.zeros(globals()[f'num{isotope}'])
		globals()[f'dpa{isotope}_MT{MTreac}'] = numpy.zeros(globals()[f'num{isotope}'])
		globals()[f'NRTdpa{isotope}_MT{MTreac}'] = numpy.zeros(globals()[f'num{isotope}'])
		globals()[f'Er{isotope}_MT{MTreac}'] = numpy.interp(globals()[f'En{isotope}_MT1'], globals()[f'En{isotope}_MT{MTreac}'], globals()[f'Er1{isotope}_MT{MTreac}'])
		globals()[f'XS{isotope}_MT{MTreac}'] = numpy.interp(globals()[f'En{isotope}_MT1'], globals()[f'En{isotope}_MT{MTreac}'], globals()[f'XS1{isotope}_MT{MTreac}'])
		globals()[f'Tdam{isotope}_MT{MTreac}'] = numpy.interp (globals()[f'Er{isotope}_MT{MTreac}'], globals()[f'E_R_{element}_srim'], globals()[f'T_dam_{element}_srim'])
		globals()[f'numvacSi28_MT{MTreac}'] = numpy.interp (globals()[f'Tdam{isotope}_MT{MTreac}'], globals()[f'T_dam_{element}_srim'], globals()[f'num_vac_{element}_srim'])
		globals()[f'DamEff{isotope}_MT{MTreac}'] = numpy.interp (globals()[f'Tdam{isotope}_MT{MTreac}'], globals()[f'T_dam_{element}_ref'], globals()[f'Dam_eff_{element}_ref'])
		for i in range(globals()[f'num{isotope}']):
			if (globals()[f'Tdam{isotope}_MT{MTreac}'][i] < globals()[f'Ed_{element}']):
				globals()[f'dpa{isotope}_MT{MTreac}'][i] = 0.0
				globals()[f'NRTdpa{isotope}_MT{MTreac}'][i] = 0.0
			if (globals()[f'Tdam{isotope}_MT{MTreac}'][i] >= globals()[f'Ed_{element}'] and globals()[f'Tdam{isotope}_MT{MTreac}'][i] < 2*globals()[f'Ed_{element}']/0.8):
				globals()[f'dpa{isotope}_MT{MTreac}'][i] = 1.0 * globals()[f"XS{isotope}_MT{MTreac}"][i]
				globals()[f'NRTdpa{isotope}_MT{MTreac}'][i] = 1.0 * globals()[f"XS{isotope}_MT{MTreac}"][i]
			if (globals()[f'Tdam{isotope}_MT{MTreac}'][i] >= 2*globals()[f'Ed_{element}']/0.8):
				globals()[f'dpa{isotope}_MT{MTreac}'][i] = 0.8/(2*globals()[f'Ed_{element}']) * globals()[f'DamEff{isotope}_MT{MTreac}'][i] * globals()[f'Tdam{isotope}_MT{MTreac}'][i] \
				*  globals()[f"XS{isotope}_MT{MTreac}"][i]
				globals()[f'NRTdpa{isotope}_MT{MTreac}'][i] = 0.8/(2*globals()[f'Ed_{element}']) * globals()[f'Tdam{isotope}_MT{MTreac}'][i] * globals()[f"XS{isotope}_MT{MTreac}"][i]
	os.chdir('../')

## Calculate dpa cross section in the target

## Summing over reaction channels in each isotope
## total dpa cross section in the target

for isotope in isotopes_evaluated:
	try:
		(element, mass_num) = separateSymbolMassNumber(isotope)
	except ValueError as e:
		print(f"Error: {e}")
	globals()[f'dpa{isotope}_MT1not102'] = numpy.zeros(globals()[f'num{isotope}'])
	globals()[f'NRTdpa{isotope}_MT1not102'] = numpy.zeros(globals()[f'num{isotope}'])
	for MT in globals()[f'MTreac{isotope}']:
		if (MT != 102):
			for i in range(globals()[f'num{isotope}']):
				if (globals()[f'Tdam{isotope}_MT{MT}'][i] >= globals()[f'Ed_bnd_{element}']):
					globals()[f'dpa{isotope}_MT1not102'][i] = globals()[f'dpa{isotope}_MT1not102'][i] + globals()[f'dpa{isotope}_MT{MT}'][i]
					globals()[f'NRTdpa{isotope}_MT1not102'][i] = globals()[f'NRTdpa{isotope}_MT1not102'][i] + globals()[f'NRTdpa{isotope}_MT{MT}'][i]


## make unique En grid out of En's of all isotopes

all_elements = ''.join(str(item) for item in elements_target)
globals()[f'En{all_elements}'] = numpy.concatenate([globals()[f'En{isotope}_MT1'] for isotope in isotopes_evaluated])
globals()[f'En{all_elements}_unique'] = numpy.unique(globals()[f'En{all_elements}'])
num_data = len(globals()[f'En{all_elements}_unique'])

## Find all quantities to the unique energy grid

for isotope in isotopes_evaluated:
	globals()[f'XS{isotope}_MT102unique'] = numpy.interp(globals()[f'En{all_elements}_unique'], globals()[f'En{isotope}_MT102'], globals()[f'XS{isotope}_MT102'])

	globals()[f'XS{isotope}_MT1unique'] = numpy.interp(globals()[f'En{all_elements}_unique'], globals()[f'En{isotope}_MT1'], globals()[f'XS{isotope}_MT1'])

	globals()[f'XS{isotope}_MT1not102unique'] = numpy.zeros(num_data)
	for i in range(num_data):
		globals()[f'XS{isotope}_MT1not102unique'][i] = globals()[f'XS{isotope}_MT1unique'][i] - globals()[f'XS{isotope}_MT102unique'][i]

	globals()[f'dpa{isotope}_MT102unique'] = numpy.interp(globals()[f'En{all_elements}_unique'], globals()[f'En{isotope}_MT102'], globals()[f'dpa{isotope}_MT102'])

	globals()[f'NRTdpa{isotope}_MT102unique'] = numpy.interp(globals()[f'En{all_elements}_unique'], globals()[f'En{isotope}_MT102'], globals()[f'NRTdpa{isotope}_MT102'])

	globals()[f'dpa{isotope}_MT1not102unique'] = numpy.interp(globals()[f'En{all_elements}_unique'], globals()[f'En{isotope}_MT1'], globals()[f'dpa{isotope}_MT1not102'])

	globals()[f'NRTdpa{isotope}_MT1not102unique'] = numpy.interp(globals()[f'En{all_elements}_unique'], globals()[f'En{isotope}_MT1'], globals()[f'NRTdpa{isotope}_MT1not102'])


	## Abundances of Si isotopes for finding the cross sections
	## and recoil energies

	globals()[f'ab{isotope}'] = percent_abundances_all[isotopes_evaluated.index(isotope)]

for element in elements_target:
	globals()[f'XS{element}_MT102unique'] = numpy.zeros(num_data)
	globals()[f'dpa{element}_MT102unique'] = numpy.zeros(num_data)
	globals()[f'NRTdpa{element}_MT102unique'] = numpy.zeros(num_data)

	globals()[f'XS{element}_MT1not102unique'] = numpy.zeros(num_data)
	globals()[f'dpa{element}_MT1not102unique'] = numpy.zeros(num_data)
	globals()[f'NRTdpa{element}_MT1not102unique'] = numpy.zeros(num_data)


for i in range (num_data):
	for element in elements_target:
		for isotope in isotopes_evaluated:
			try:
				(element_iso, mass_num) = separateSymbolMassNumber(isotope)
			except ValueError as e:
				print(f"Error: {e}")

			if (element == element_iso):
				globals()[f'XS{element}_MT102unique'][i] = globals()[f'XS{element}_MT102unique'][i] + globals()[f'ab{isotope}']*globals()[f'XS{isotope}_MT102unique'][i]

				globals()[f'dpa{element}_MT102unique'][i] = globals()[f'dpa{element}_MT102unique'][i] + globals()[f'ab{isotope}']*globals()[f'dpa{isotope}_MT102unique'][i]

				globals()[f'NRTdpa{element}_MT102unique'][i] = globals()[f'NRTdpa{element}_MT102unique'][i] + globals()[f'ab{isotope}']*globals()[f'NRTdpa{isotope}_MT102unique'][i]

				globals()[f'XS{element}_MT1not102unique'][i] = globals()[f'XS{element}_MT1not102unique'][i] + globals()[f'ab{isotope}']*globals()[f'XS{isotope}_MT1not102unique'][i]
				
				globals()[f'dpa{element}_MT1not102unique'][i] = globals()[f'dpa{element}_MT1not102unique'][i] + globals()[f'ab{isotope}']*globals()[f'dpa{isotope}_MT1not102unique'][i]

				globals()[f'NRTdpa{element}_MT1not102unique'][i] = globals()[f'NRTdpa{element}_MT1not102unique'][i] + globals()[f'ab{isotope}']*globals()[f'NRTdpa{isotope}_MT1not102unique'][i]



## Combine the dpa cross sections in elements

globals()[f'dpa{all_elements}'] = numpy.zeros(num_data)
globals()[f'NRTdpa{all_elements}'] = numpy.zeros(num_data)
#globals()[f'simple_dpa{all_elements}'] = numpy.zeros(num_data)
#globals()[f'simple_NRTdpa{all_elements}'] = numpy.zeros(num_data)


# print data to output files
ofile = open('CombinU-XS.out', 'w')
print(num_data, file = ofile)

for element in elements_target:
	globals()[f'ofile1{element}'] = open('CombinU-Probabilities-'+element+'.out', 'w')
	print(num_data, '(n,g) / (n, other)', file = globals()[f'ofile1{element}'])

for i in range(num_data):
	denominator_1not102 = 0
	for element in elements_target:
		denominator_1not102 = denominator_1not102 + globals()[f'XS{element}_MT1not102unique'][i]
	for element in elements_target:
		globals()[f'f{element}1not102'] = element_stoich[elements_target.index(element)]*globals()[f'XS{element}_MT1not102unique'][i]/denominator_1not102
	denominator102 = 0
	for element in elements_target:
		denominator102 = denominator102 + globals()[f'XS{element}_MT102unique'][i]
	for element in elements_target:
		globals()[f'f{element}102'] = element_stoich[elements_target.index(element)]*globals()[f'XS{element}_MT102unique'][i]/denominator102

	for element in elements_target:
		globals()[f'dpa{all_elements}'][i] = globals()[f'dpa{all_elements}'][i] + globals()[f'f{element}102']*globals()[f'dpa{element}_MT102unique'][i] + globals()[f'f{element}1not102']*globals()[f'dpa{element}_MT1not102unique'][i]
		globals()[f'NRTdpa{all_elements}'][i] = globals()[f'NRTdpa{all_elements}'][i] + globals()[f'f{element}102']*globals()[f'NRTdpa{element}_MT102unique'][i] + globals()[f'f{element}1not102']*globals()[f'NRTdpa{element}_MT1not102unique'][i]
		
		#globals()[f'simple_dpa{all_elements}'][i] = globals()[f'simple_dpa{all_elements}'][i] + element_stoich[elements_target.index(element)] * globals()[f'dpa{element}_unique'][i]
		#globals()[f'simple_NRTdpa{all_elements}'][i] = globals()[f'simple_NRTdpa{all_elements}'][i] + element_stoich[elements_target.index(element)] * globals()[f'NRTdpa{element}_unique'][i]

	print(globals()[f'En{all_elements}_unique'][i], globals()[f'NRTdpa{all_elements}'][i], globals()[f'dpa{all_elements}'][i], file = ofile)

	for element in elements_target:
		print(globals()[f'En{all_elements}_unique'][i], globals()[f'f{element}102']/elements_stoich[elements_target.index(element)], globals()[f'f{element}1not102']/elements_stoich[elements_target.index(element)], file = globals()[f'ofile1{element}'])

for element in elements_target:
	globals()[f'ofile1{element}'].close()

## multigroup dpa cross sections

(Ngl, Eg, globals()[f'grouped_dpa{all_elements}']) = groupmulti(7, num_data, globals()[f'En{all_elements}_unique'], globals()[f'dpa{all_elements}'])
(Ngl, Eg, globals()[f'grouped_NRTdpa{all_elements}']) = groupmulti(7, num_data, globals()[f'En{all_elements}_unique'], globals()[f'NRTdpa{all_elements}'])

print(Ngl, file = ofile)
for i in range(Ngl):
	print(Eg[i], globals()[f'grouped_NRTdpa{all_elements}'][i], globals()[f'grouped_dpa{all_elements}'][i], file = ofile)

ofile.close()
