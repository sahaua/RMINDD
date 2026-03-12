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
		data = eachlineinfo(line)
		MAT = int(data[6]); MF = int(data[7]); MT  = int(data[8])
		if (MAT != -1):
			if (MF == 3):
				if (MT == 1):
					line = ifile_preprocessedENDF6.readline() 
					data = eachlineinfo(line)
					QM = float(data[0]); QI =  float(data[1]); NR = int(data[4]); NPt = int(data[5])
					Et = [0]*NPt
					LR = int(ifile_preprocessedENDF6.readline().split()[1])
					i = 0
					while (i < NPt):
						line = ifile_preprocessedENDF6.readline()
						data = eachlineinfo(line)
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
