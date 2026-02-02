#	>> Code: RMINDD (Radiation Matter Interaction from Nuclear Data and Damage)
#	>> Perform: Calcuation of Effects of Radiation in Matter from ENDF-6 files
#	>> Author: Dr. Uttiyoarnab Saha
#	>> Version and Date: 1.0 and 01/07/2022
#	>> Last modified: 01/07/2022, Kolkata
#	>> Update: 01/07/2022
#	>> Major changes: 
#
#=========================================================================================

# This program ReadU.py is a part of RMINDD.py. It reads input file and the ENDF-6 nuclear
# data files.

import sys

def readInputFile(input_file_name):	
	ifile = open(input_file_name, 'r')
	ifile.seek(0, 2)			# go to the end of file
	eof = ifile.tell()			# get the end-of-file position
	ifile.seek(0, 0)			# go to the start of file

	while (ifile.tell() != eof):
		line = ifile.readline()
		# blank line
		if (line.strip(' ') == '\n'):
			continue
		data1 = line.split()
		data = []

		for i in range(len(data1)):
			data.append(data1[i].strip(' '))

		# comment line
		if (data[0] == '#'):
			continue

		if (data[0] == 'Num_modules'):
			num_modules = int(data[2])

		if (data[0] == 'Module_name'):
			module_name = data[2]

		if (data[0] == 'Element_isotope'):
			el_iso_name = data[2]

		if (data[0] == 'Raw_ENDF6_file'):
			raw_ENDF6_file = data[2]

		if (data[0] == 'Processed_ENDF6_file'):
			processed_ENDF6_file = data[2]

		if (data[0] == 'MAT_num'):
			MAT_num = int(data[2])

		if (data[0] == 'Num_reac'):
			num_reac = int(data[2])
			num_reac_array = [0]*num_reac

		if (data[0] == 'Reaction_num'):
			for i in range(2, num_reac+2):
				num_reac_array[i-2] = int(data[i])
			if (num_reac > 1):
				for value in num_reac_array:
					if (value == 7):
						print("Please do only partials or only total. In total all partials will also be done.")
						sys.exit()

		if (data[0] == 'Atom_displ_model'):
			atom_displ_model = data[2]

		if (data[0] == 'Threshold_lattice_disp_en'):
			threshold_latt_disp_en = float(data[2])

		if (data[0] == 'b_arcdpa'):
			b_arcdpa = float(data[2])

		if (data[0] == 'c_arcdpa'):
			c_arcdpa = float(data[2])

		if (data[0] == 'Multigroup'):
			multigroup = data[2]

		if (data[0] == 'Energy_group_type_index'):
			en_group_type = int(data[2])

		if (data[0] == 'Input_n_spectrum'):
			input_n_spec = int(data[2])

		if (data[0] == 'Input_en_group'):
			input_en_group = int(data[2])

		if (data[0] == 'Num_MT_to_multigroup'):
			num_MT_multigroup = int(data[2])
			if (num_MT_multigroup > 0):
				num_MT_group_array = [0]*num_MT_multigroup

		if (data[0] == 'MTs_to_multigroup' and num_MT_multigroup > 0):
			if (len(data) != num_MT_multigroup + 2):
				print('Number of MTs given for multigrouping is not equal to what is mentioned.')
				sys.exit()
			for i in range(2, num_MT_multigroup+2):
				num_MT_group_array[i-2] = int(data[i])

		if (data[0] == 'Num_group_limits'):
			num_group_limits = int(data[2])

		if (data[0] == 'Num_fine_en_points'):
			num_fine_en_points = int(data[2])

		if (data[0] == 'Num_partial_reac_to_sum'):
			num_partial_reac_tosum = int(data[2])
			if (num_partial_reac_tosum > 0):
				partial_reac_tosum = [0]*num_partial_reac_tosum

		if (data[0] == 'Partial_reac_to_sum' and num_partial_reac_tosum > 0):
			if (len(data) != num_partial_reac_tosum + 2):
				print('Number of partial reactions given to sum is not equal to what is mentioned.')
				sys.exit()
			for i in range(2, num_partial_reac_tosum + 2):
				partial_reac_tosum[i-2] = int(data[i])

	ifile.close()

readInputFile('inputfile.txt')

# The data in one line are explicitly extracted
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
	
	s = ''
	for char in range(70,72):			# MF data of 2 places
		s = s + line[char]
	data[7] = s.lstrip(' ')
	
	s = ''
	for char in range(72,75):			# MT data of 3 places
		s = s + line[char]
	data[8] = s.lstrip(' ')
	
	for i in range(6):
		data[i] = "E-".join(data[i].split('-'))
		data[i] = "E+".join(data[i].split('+'))
		data[i] = data[i].lstrip('E')
	
	return(data)
# ======================================================================

# The data in lines of different types of ENDF-6 file formats 
# are explicitly extracted.

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

def lineType3Info (filehandle,numdata,numvariables):
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
			for j in range(0,6,2):
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

# ======================================================================

global E_MT1
# 1 = elastic scattering
global MT_values1 = [2]
# 2 = inelastic scattering 
global MT_values2 = [4] + numpy.linspace(51, 91, num = 41).astype(int)
# 3 = (n, CPO) reactions
global MT_values3 = [11, 22, 23, 24, 25, 28, 29, 30, 32, 33, 34, 35, 36, 41, \
42, 44, 45, 103, 104, 105, 106, 107, 108, 109, 111, 112, 113, 114, 115, 116, 117]
MT_values3 = MT_values3 + numpy.linspace(600, 849, num = 250).astype(int)
# 4 = (n, xn) reactions
global MT_values4 = [16, 17, 37]
# 5 = (n, g) reaction
global MT_values5 = [102]
# 6 = (n, anything) reactions
global MT_values6 = [5]
# 7 = all reactions (total)
global MT_values7 = MT_values1 + MT_values2 + MT_values3 + MT_values4 + MT_values5 + MT_values6

# Function readFile3 reads the required cross sections and stores them 
# in separate arrays, using the processed ENDF-6 file
# It will also return the general unique unionised energy array

def readFile3 (processed_ENDF6_file, num_reac_array):
	ifile = open (processed_ENDF6_file, 'r')
	while True:
		line = ifile.readline()
		data = eachLineInfo(line)
		MAT = int(data[6]); MF = int(data[7]); MT  = int(data[8])
		if (MAT != -1):
			if (MF == 3):
				if (MT == 1):
					line = ifile.readline() 
					data = eachLineInfo(line)
					QM = float(data[0]); QI =  float(data[1])
					NR = int(data[4]); globals()[f'NP_MT{MT}'] = int(data[5])
					globals()[f'E_MT{MT}'] = [0]*globals()[f'NP_MT{MT}']; globals()[f'sig_MT{MT}'] = [0]*globals()[f'NP_MT{MT}']
					LR = int(ifile.readline().split()[1])
					(globals()[f'E_MT{MT}'], globals()[f'sig_MT{MT}']) = lineType3Info(ifile,globals()[f'NP_MT{MT}'],2)
				
				# store energy and cross sections of all reactions needed
				for l in sorted(MT_values7):
					if (MT == l):
						line = ifile.readline()
						data = eachLineInfo(line)
						QM = float(data[0]); QI =  float(data[1])
						NR = int(data[4]); globals()[f'NP_MT{l}'] = int(data[5])
						globals()[f'E_MT{l}'] = [0]*globals()[f'NP_MT{l}']; globals()[f'sig_MT{l}'] = [0]*globals()[f'NP_MT{l}']
						LR = int(ifile.readline().split()[1])
						(globals()[f'E_MT{l}'], globals()[f'sig_MT{l}']) = lineType3Info(ifile,globals()[f'NP_MT{l}'],2)
		else:
			break
	ifile.close()

	# make unique common energy

	global Etu = numpy.array(E_MT1)
	Etu = numpy.unique(Etu)
	global NPt = len(Etu)

# File 4 contains secondary energy angle data, required to be read using
# the raw ENDF-6 file for elastic scattering, inelastic scattering and
# discrete CPO reactions.

def readFile4(raw_ENDF6_file, num_reac_array):
	ifile = open (processed_ENDF6_file, 'r')
	while True:
		line = ifile.readline()
		data = eachLineInfo(line)
		MAT = int(data[6]); MF = int(data[7]); MT  = int(data[8])
		if (MF == 3 and MT == 0):
			ifile.readline()
			line = ifile.readline()
			data = eachLineInfo(line)
			iflspace = 0
			for element in data:
				if (element == ''):
					iflspace = 1
					break
			if (iflspace == 0):
				ZA = float(data[0]); AWR = float(data[1]); l1 = int(data[2]);\
				LTT = int(data[3]); NK = int(data[4]); l2 = int(data[5]);\
				MAT = int(data[6]); MF = int(data[7]); MT = int(data[8])
		if (MAT != -1):
			if (MF == 4):
				# for elastic scattering reaction ....
				if 1 or 7 in num_reac_array:
					if (MT == 2):
						ifile.readline()
						if (LTT == 0):
							alc = [0]*65
						if (LTT == 3 or LTT == 1):
							line = ifile.readline()
							data = eachLineInfo(line)
							NE1 = int(data[5])
							# Legendre polynomial coefficients
							ifile.readline()
							al = numpy.zeros((NE1,65)); EL = [0]*NE1; alc = [0]*65
							for i in range(NE1):
								line = ifile.readline()
								data = eachLineInfo(line)
								T = float(data[0]); EL[i] = float(data[1]); NL = int(data[4])
								NL = NL+1
								al[i][0] = 1
								temporary = [0]*NL
								temporary = lineType3Info(ifile,NL,1)
								for j, value in enumerate(temporary, 1):
									al[i][j] = value
						if (LTT == 3 or LTT == 2):
							line = ifile.readline()
							data = eachLineInfo(line)
							NE2 = int(data[5])
							# Tabulated Probability
							ifile.readline()
							Enf = [0]*NE2; NPr = [0]*NE2; cdata = numpy.zeros((NE2,201))
							fdata = numpy.zeros((NE2,201)); ftotal = numpy.zeros((NE2,64))
							for i in range(NE2):
								line = ifile.readline()
								data = eachLineInfo(line)
								T = float(data[0]); Enf[i] = float(data[1])
								line = ifile.readline()
								data = eachLineInfo(line)
								NPr[i] = int(data[0])
								temporary1 = [0]*201
								temporary2 = [0]*201
								(temporary1, temporary2) = lineType3Info(ifile,NPr[i],2)
								for j, (value1, value2) in enumerate(zip(temporary1, temporary2), 0):
									cdata[i][j] = value1
									fdata[i][j] = value2
				# for inelastic scattering reaction ....
				if 2 or 7 in num_reac_array:
					for l in MT_values2:
						if (MT == l):
							if (MT < 91):
								if4d = 1
							if (MT == 91):
								if4c = 1
							LTT = LTTv
							ifile.readline()
							line = ifile.readline()
							(ZA,AWR,LI,LCT,L2,L3,MAT,MF,MT) = lineType1Info(line)
							if (LTT == 1  and  LI == 0):
								line = ifile.readline()
								(C1,C2,L1,L2,NR,NE4,MAT,MF,MT) = lineType2Info(line)
								(NBT, INTr) = lineType3Info(ifile,NR,2)
								for i in range(NE4):
									line = ifile.readline()
									data = eachLineInfo(line)
									c1 = float(data[0]); En4[i] = float(data[1]); NL4 = int(data[4])
									NL4 = NL4 + 1
									al4[i][0] = 1
									temporary = [0]*NL4
									temporary = lineType3Info (ifile,NL4,1)
									for j, value in enumerate(temporary, 1):
										al4[i][j] = value
				
				# for (n, CPO) discrete level reactions
				if 3 or 7 in num_reac_array:
					n_cpo_ifl4 = 0
					n_cpo_ifspad4 = 0 #flag for secondary particle angular data
					n_cpo_ifspad4al = 0 # flag for secondary particle angular data in 'al' coefficients
					n_cpo_ifspad4muf = 0 # flag for secondary particle angular data in 'mu,f' form

					for l in MT_values3:
						if (MT == l):
							n_cpo_ifl4 = 1; ifspad4 = 1
							ifile.readline()
							if (LTT == 3 or LTT == 1):
								n_cpo_ifspad4al = 1
								line = ifile.readline()
								data = eachLineInfo(line)
								n_cpo_NE1 = int(data[5])
								## Legendre polynomial coefficients
								ifile.readline()
								n_cpo_al4 = [[0]*65]*n_cpo_NE1; n_cpo_EL = [0]*n_cpo_NE1
								for i in range (n_cpo_NE1):
									line = ifile.readline()
									data = eachLineInfo(line)
									T = float(data[0]); n_cpo_EL[i] = float(data[1]); n_cpo_NL4 = int(data[4])
									n_cpo_NL4 = n_cpo_NL4 + 1
									n_cpo_al4[i][0] = 1
									temporary = [0]*NL4
									temporary = lineType3Info (ifile,NL4,1)
									for j, value in enumerate (temporary, 1):
										n_cpo_al4[i][j] = value

							if (LTT == 3 or LTT == 2):
								n_cpo_ifspad4muf = 1
								n_cpo_NE2 = int(ifile.readline().split()[5])
								## Probability
								ifile.readline()
								n_cpo_Enf = [0]*n_cpo_NE2; n_cpo_NPr = [0]*n_cpo_NE2; n_cpo_cdata4 = numpy.zeros((n_cpo_NE2,201))
								n_cpo_fdata4 = numpy.zeros((n_cpo_NE2,201)); n_cpo_ftotal= numpy.zeros((n_cpo_NE2,64))
								n_cpo_fmuE = numpy.zeros((NP,64)); fpr = [0]*64
								for i in range (NE2):
									line = ifile.readline()
									data = eachLineInfo(line)
									T = float(data[0]); Enf[i] = float(data[1]) 
									NPr[i] = int(ifile.readline().split()[0])
									temporary1 = [0]*NPr[i]
									temporary2 = [0]*NPr[i]
									(temporary1,temporary2) = lineType3Info(ifile,NPr[i],2)
									for j, (value1, value2) in enumerate(zip(temporary1, temporary2), 0):
										cdata4[i][j] = value1
										fdata4[i][j] = value2
		else:
			break
	ifile102.close()

# File 5 contains secondary energy data, required to be read using
# the raw ENDF-6 file for inelastic scattering, (n, xn) reactions

def readFile5(raw_ENDF6_file, num_reac_array):
	ifile = open (processed_ENDF6_file, 'r')
		while True:
			line = ifile.readline()
			if (line == ''):
				break
			data = eachlineinfo(line)
			MAT = int(data[6]); MF = int(data[7]); MT = int(data[8])	
			if (MT == 0):
				if (MF == 5):
					line = ifile.readline()
					(ZAv,AWRv,L1,L2,NK5,L3,MAT,MF,MT) = line_type1_info(line)
				if (MF < 5):
					ifile.readline()
					line = ifile.readline()
					(ZAv,AWRv,L1,L2,NK5,L3,MAT,MF,MT) = line_type1_info(line)
				if (MAT != -1):
					if (MF == 5):
						# for inelastic scattering reaction ....
						if 2 or 7 in num_reac_array:
							ift5c = 0
							if5c = 0
							if (MT == 91):
									if5c = 1
									NK5 = NKv
									for NSS in range(NK5):
										ifile.readline()
										line = ifile.readline()
										(C1,C2,L1,LF[NSS],NR,NP5[NSS],MAT,MF,MT) = line_type2_info(line)
										(NBT, INTr) = line_type3_info(ifile,NR,2)
										line_type3_info(ifile,NP5[NSS],2)
										if (LF[NSS] == 1):
											ift5c = 1
											ifile.readline()
											line = ifile.readline()
											(C1,C2,L1,L2,NR,NE5[NSS],MAT,MF,MT) = line_type2_info(line)
											(NBT, INTr) = line_type3_info(ifile,NR,2)
											for i in range (NE5[NSS]):
												ifile.readline()
												line = ifile.readline()
												(C1,En5[NSS][i],L1,L2,NR,NF5[NSS][i],MAT,MF,MT) = line_type2_info(line)
												(NBT, INTr) = line_type3_info(ifile,NR,2)
												temporary1 = [0]*NF5[NSS][i]
												temporary2 = [0]*NF5[NSS][i]
												(temporary1,temporary2) = line_type3_info(ifile,NF5[NSS][i],2)
												for j, (value1, value2) in enumerate(zip(temporary1, temporary2), 0):
													Enp5[NSS][i][j] = value1
													f5[NSS][i][j] = value2

										if (LF[NSS] == 9):
											ifile.readline()
											line = ifile.readline()
											(C1,C2,L1,L2,NR,NE5[NSS],MAT,MF,MT) = line_type2_info(line)
											(NBT, INTr) = line_type3_info(ifile,NR,2)
											temporary1 = [0]*NE5[NSS]
											temporary2 = [0]*NE5[NSS]
											(temporary1,temporary2) = line_type3_info(ifile,NE5[NSS],2)
											for j, (value1, value2) in enumerate(zip(temporary1, temporary2), 0):
												En5[NSS][N] = value1
												tht[NSS][N] = value2

						# for (n, xn) reactions ....
						if 4 or 7 in num_reac_array:
							ift5 = 0
							if5 = 0
							for l in MT_values4:
								if (MT == l):
									ZA = ZAv
									AWR = AWRv
									if5 = 1
									for NSS in range (NK5):
										line = ifile.readline()
										(C1,C2,L1,LF[NSS],NR,NP5[NSS],MAT,MF,MT) = line_type2_info(line)
										(NBT, INTr) = line_type3_info(ifile,NR,2)
										temporary1 = [0]*NP5[NSS]
										temporary2 = [0]*NP5[NSS]
										(temporary1,temporary2) = line_type3_info(ifile,NP5[NSS],2)
										for N, (value1, value2) in enumerate(zip(temporary1, temporary2), 0):
											Eint[NSS][N] = value1
											p[NSS][N] = value2
										if (LF[NSS] == 1):
											line = ifile.readline()
											(C1,C2,L1,L2,NR,NE5[NSS],MAT,MF,MT) = line_type2_info(line)
											(NBT, INTr) = line_type3_info(ifile,NR,2)
										for i in range (NE5[NSS]):
											line = ifile.readline()
											(C1,En5[NSS][i],L1,L2,NR,NF5[NSS][i],MAT,MF,MT) = line_type2_info(line)
											(NBT, INTr) = line_type3_info(ifile,NR,2)
											temporary1 = [0]*NF5[NSS][i]
											temporary2 = [0]*NF5[NSS][i]
											(temporary1,temporary2) = line_type3_info(ifile,NF5[NSS][i],2)
											for j, (value1, value2) in enumerate(zip(temporary1, temporary2), 0):
												Enp5[NSS][i][j] = value1
												f5[NSS][i][j] = value2
										if (LF[NSS] == 9):
											ift5 = 1
											U = C1
											line = ifile.readline()
											(C1,C2,L1,L2,NR,NE5[NSS],MAT,MF,MT) = line_type2_info(line)
											(NBT, INTr) = line_type3_info(ifile,NR,2)
											temporary1 = [0]*NE5[NSS]
											temporary2 = [0]*NE5[NSS]
											(temporary1,temporary2) = line_type3_info(ifile,NE5[NSS],2)
											for N, (value1, value2) in enumerate(zip(temporary1, temporary2), 0):
												En5[NSS][N] = value1
												tht[NSS][N] = value2
			else:
				break
		ifile.close()
def readFile6():

def readFile12():

def readFile15():
