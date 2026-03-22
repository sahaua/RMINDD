''' 
Code: RMINDD - (Radiation-Matter Interaction and Damage calculation using Nuclear Data)
Perform: Calculation of metrics of neutron radiation damage in a material using ENDF-6 files
Module: UtilsU.py -- compilation some useful functions used in other modules
Author: Uttiyoarnab Saha
Version and Date: 1.0 and 01/07/2022
Last modified: 01/07/2022, Kolkata
Update: 01/07/2022
Major changes: 

=========================================================================================
'''

import numpy

'''
#======= Make the required unique common energy =======*

It creates an unique energy array out of the MF3 MT1 energy points
given in the pre-processed ENDF-6 file which is used as the common 
energy array for partial and total dpa and heating cross section 
computation. This energy array is called unique because any repetition
(may occur at the dense resonances regions) in the pre-processed energy 
points are found out and removed, so that each energy point is present
only once.
'''

def uqce (ofile_outRMINDD, ifile_preprocessedENDF6):
	## Et=Energy array in MT=1, Etu=unique of Et
	## extraction of total energy points

	ifile_preprocessedENDF6.seek(0, 0)

	while True:
		line = ifile_preprocessedENDF6.readline()
		if (line == ''):
			break  
		data = eachLineInfo(line)
		MAT = int(data[6]); MF = int(data[7]); MT  = int(data[8])
		if (MAT != -1):
			if (MF == 3):
				if (MT == 1):
					line = ifile_preprocessedENDF6.readline() 
					data = eachLineInfo(line)
					QM = float(data[0]); QI =  float(data[1]); NR = int(data[4]); NPt = int(data[5])
					Et = [0]*NPt
					LR = int(ifile_preprocessedENDF6.readline().split()[1])
					i = 0
					while (i < NPt):
						line = ifile_preprocessedENDF6.readline()
						data = eachLineInfo(line)
						for j in range(0,5,2):
							if (data[j] != ''):
								Et[i] = float(data[j])
								i += 1
							else:
								i += 1
								break
		else:
			break

	print('', file = ofile_outRMINDD)
	print(NPt,' Total cross sections energy points', file = ofile_outRMINDD)

	## make unique common energy
	Etu = numpy.array(Et)
	Etu = numpy.unique(Etu)
	NPt = len(Etu)

	print('', file = ofile_outRMINDD)
	print(NPt,' Unique total cross sections energy points', file = ofile_outRMINDD)

	return(NPt, Etu)

'''
The following four functions helps to read the formatted data in different sections of
the ENDF-6 file. For specific data, the corresponding function is called during reading data.

#========The data in each line are explicitly extracted=======*
'''

def eachLineInfo(line):
	data = [0]*9
	for i in range(6):
		s = ''						# 6 data each of 11 places
		for char in range(i*11,(i+1)*11):
			s = s + line[char]
		data[i] = s.lstrip(' ')

	s = ''
	for char in range(66,70):			# MAT data of 4 places
		s = s + line[char]
	data[6] = s.lstrip(' ')
	# this is done because some ENDF-6 files only give a blank space here in the first line (creates problem)
	if ((data[6]) == ''):
		data[6] = 1

	s = ''
	for char in range(70,72):			# MF data of 2 places
		s = s + line[char]
	data[7] = s.lstrip(' ')

	s = ''
	for char in range(72,75):			# MT data of 3 places
		s = s + line[char]
	data[8] = s.lstrip(' ')

	for i in range(6):
		putE = 1
		for x in data[i]:
			if (x == 'E' or x == 'e'):
				putE = 0
				break
		if (putE == 1):
			data[i] = "E-".join(data[i].split('-'))
			data[i] = "E+".join(data[i].split('+'))
			data[i] = data[i].lstrip('E')

	return(data)

#========The data in lines of different types are explicitly extracted=======*

def lineType1Info(line):
	data = eachLineInfo(line)
	dataV1 = 0; dataV2 = 0; dataV3 = 0; dataV4 = 0; dataV5 = 0; dataV6 = 0
	dataV7 = 0; dataV8 = 0; dataV9 = 0
	iflspace = 0
	for element in data:
		if (element == ''):
			iflspace = 1
			break
	if (iflspace == 0):
		dataV1 = float(data[0]); dataV2 = float(data[1]); dataV3 = int(data[2])
		dataV4 = int(data[3]); dataV5 = int(data[4]); dataV6 = int(data[5])
		dataV7 = int(data[6]); dataV8 = int(data[7]); dataV9 = int(data[8])
	return(dataV1,dataV2,dataV3,dataV4,dataV5,dataV6,dataV7,dataV8,dataV9)

def lineType2Info(line):
	data = eachLineInfo(line)
	dataV1 = float(data[0]); dataV2 = float(data[1]); dataV3 = int(data[2])
	dataV4 = int(data[3]); dataV5 = int(data[4]); dataV6 = int(data[5])
	dataV7 = int(data[6]); dataV8 = int(data[7]); dataV9 = int(data[8])
	return(dataV1,dataV2,dataV3,dataV4,dataV5,dataV6,dataV7,dataV8,dataV9)

def lineType3Info(filehandle,numdata,numvariables):
	# numvariables (1 / 2) denotes no. of variables the given data has to be read into
	if (numvariables == 2):
		xdata = [0]*numdata
		ydata = [0]*numdata
	if (numvariables == 1):
		xdata = [0]*numdata

	i = 0
	if (numvariables == 2):
		while (i < numdata):
			line = filehandle.readline()
			data = eachLineInfo(line)
			# run_limit variable is introduced because in some files 
			#(e.g. ENDF/B-VIII.0, Si28
			#  1.000000+0 1.000000+0          0          2          1          21425 6 51
			#	2          2          0          0          0          01425 6 51)
			# extra data (more than 2) are given in the interpolation ranges specification line!
			run_limit = 6
			if (2*numdata <= run_limit):
				run_limit = 2*numdata
			for j in range(0,run_limit,2):
				if (data[j] != ''):
					xdata[i] = float(data[j])
					ydata[i] = float(data[j+1])
					i += 1
				else:
					i += 1
					break
	if (numvariables == 1):
		while (i < numdata):
			line = filehandle.readline()
			data = eachLineInfo(line)
			for j in range(6):
				if (data[j] != ''):
					xdata[i] = float(data[j])
					i += 1
				else:
					i += 1
					break

	if (numvariables == 2):
		return(xdata, ydata)
	if (numvariables == 1):
		return(xdata)

#========The Data from Output File of PKA Matrix in Specific Format==========*

def lineType4Info(filehandle,nre,numdata):
	xdata = [0]*numdata
	i = 0
	while (i < numdata):
		line = filehandle.readline()
		line = line.lstrip('[')
		line = line.rstrip(']')
		data = line.split()
		for num in range(len(data)):
			if (num == len(data) - 1):
				data[num] = float(data[num].lstrip("'").rstrip("']"))
			elif (num == 0):
				data[num] = float(data[num].lstrip("['").rstrip("',"))
			else:
				data[num] = float(data[num].lstrip("'").rstrip("',"))
		for j in range(len(data)):
			#print(i,j)
			xdata[i] = float(data[j])
			i += 1
	return(xdata)

'''
From a value of Z get the name of the corresponding element. 
'''
def elementFromZValue(Z):
	periodic_table = {
	1: "H", 2: "He", 3: "Li", 4: "Be", 5: "B", 6: "C", 7: "N", 8: "O", 9: "F", 10: "Ne", 11: "Na", 12: "Mg", 13: "Al",
	14: "Si", 15: "P", 16: "S", 17: "Cl", 18: "Ar", 19: "K", 20: "Ca", 21: "Sc", 22: "Ti", 23: "V", 24: "Cr", 25: "Mn",
	26: "Fe", 27: "Co", 28: "Ni", 29: "Cu", 30: "Zn", 31: "Ga", 32: "Ge", 33: "As", 34: "Se", 35: "Br", 36: "Kr", 37: "Rb",
	38: "Sr", 39: "Y", 40: "Zr", 41: "Nb", 42: "Mo", 43: "Tc", 44: "Ru", 45: "Rh", 46: "Pd", 47: "Ag", 48: "Cd", 49: "In",
	50: "Sn", 51: "Sb", 52: "Te", 53: "I", 54: "Xe", 55: "Cs", 56: "Ba", 57: "La", 58: "Ce", 59: "Pr", 60: "Nd", 61: "Pm",
	62: "Sm", 63: "Eu", 64: "Gd", 65: "Tb", 66: "Dy", 67: "Ho", 68: "Er", 69: "Tm", 70: "Yb", 71: "Lu", 72: "Hf", 73: "Ta",
	74: "W", 75: "Re", 76: "Os", 77: "Ir", 78: "Pt", 79: "Au", 80: "Hg", 81: "Tl", 82: "Pb", 83: "Bi", 84: "Po", 85: "At",
	86: "Rn", 87: "Fr", 88: "Ra", 89: "Ac", 90: "Th", 91: "Pa", 92: "U", 93: "Np", 94: "Pu", 95: "Am", 96: "Cm", 97: "Bk",
	98: "Cf", 99: "Es", 100: "Fm", 101: "Md", 102: "No", 103: "Lr", 104: "Rf", 105: "Db", 106: "Sg", 107: "Bh", 108: "Hs",
	109: "Mt", 110: "Ds", 111: "Rg", 112: "Cn", 113: "Uut", 114: "Fl", 115: "Uup", 116: "Lv", 117: "Uus", 118: "Uuo"
	}
	return(periodic_table[Z])
