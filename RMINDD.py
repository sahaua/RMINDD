## >> Code: RMINDD - (Radiation-Matter Interaction and Damage calculation using Nuclear Data)
## >> Perform: Calcuation of metrics of neutron radiation damage in an isotope of a material using ENDF-6 files
## >> Author: Uttiyoarnab Saha
## >> Version and Date: 1.0 and 01/07/2022
## >> Last modified: 01/07/2022, Kolkata
## >> Update: 01/07/2022
## >> Major changes: 
##
## =========================================================================================

import numpy, datetime, sys
from time import process_time
## import calculation modules
import RecedU, EngdepU #use TransmU

day_execution = datetime.date.today()
time_execution = datetime.datetime.now().strftime('%H:%M:%S')

outRMINDD = "Output_RMINDD.txt"
ofile_outRMINDD = open (outRMINDD, 'a')

print('~~~~ RMINDD ~~~~', file = ofile_outRMINDD)
print(day_execution, ' ', time_execution, '\n', file = ofile_outRMINDD)

## command line input for input file name
inpRMINDD = sys.argv[1]

ifile_inpRMINDD = open (inpRMINDD, 'r')

ntasks = int(ifile_inpRMINDD.readline().split()[0])

if (ntasks < 1 or ntasks > 3):
	print('**** Give Valid Value for Number of Tasks ****')

if (1 <= ntasks and ntasks <= 3):
	start_time = process_time()
	for itasks in range (ntasks):
		cmetric = ifile_inpRMINDD.readline().split()[0]

	# Perform tasks based on the results required from the relevant module
	# First one is for energy deposition to estimate dpa and heating
		if (cmetric == 'EngdepU'):

			ifile52 = open ("tape01", 'r')

			print( '~~ RMINDD-EngdepU ~~')
			print( '~~ RMINDD-EngdepU ~~', file = ofile_outRMINDD)
			print( ':Messages for you:', file = ofile_outRMINDD)
			print('--------------------', file = ofile_outRMINDD)
			print('', file = ofile_outRMINDD)
		
			ifile52.readline()
			mat = ifile52.readline().split()[5][-4:] 				# read mat number in ENDF = 52

			for i in range(3):
				ifile52.readline()
			data = ifile52.readline().split()
			iso = data[0] + data[1]
			ifile52.close()
	
			print( 'Evaluation on tape01', file = ofile_outRMINDD)
			print( iso, file = ofile_outRMINDD)
			print('', file = ofile_outRMINDD)
			matg = ifile_inpRMINDD.readline().split()[0] 			# given mat number in input = 50
			
			# Check if given mat number matches with ENDF mat number
			
			##
			print(mat, matg, file = ofile_outRMINDD)
			
			if (mat != matg):
				print( 'Error', file = ofile_outRMINDD)
				print( 'Material does not exist on tape01', file = ofile_outRMINDD)
				print('', file = ofile_outRMINDD)
				break
			
			nreac = int(ifile_inpRMINDD.readline().split()[0])   		# given number of reactions in input = 50
			nra = [0]*nreac
			data = ifile_inpRMINDD.readline().split()
			for i in range(nreac):
				nra[i] = int(data[i])
			data = ifile_inpRMINDD.readline().split()
			mdisp = int(data[0]); Ed = float(data[1]); bad = float(data[2]); cad = float(data[3]) 	#  dpa model, displ. energy, arc dpa parameters(2)
			data = ifile_inpRMINDD.readline().split()
			mgyn = int(data[0]); insp = int(data[1]) 					# multigrouping yes=1 or no=0 and given input spectrum yes = 1, no = 0
			if (mgyn == 1):
				igtype = int(ifile_inpRMINDD.readline().split()[0]) 		# type of group
				NpMTtg = int(ifile_inpRMINDD.readline().split()[0]) 		# number of partial MTs to group
									
				# applies only for total, (n, xn), total (n, CPO) in nra
				if (NpMTtg > 0):					# store partial MTs to group
					nMTprtg = [0]*NpMTtg
					data = ifile_inpRMINDD.readline().split()
					for i in range(NpMTtg):
						nMTprtg[i] = int(data[i])

			# Check if nreac is within 1 and 7 and proceed in main calculatiun (uqce), else print error
			if (1 <= nreac and nreac <= 7):
			
				EngdepU.uqce (ofile_outRMINDD,insp,nra,nreac,mdisp,Ed,bad,cad,mgyn,igtype)
			#call extractdata(nra(i))
			#call filerecds()
				if (NpMTtg > 0):
					print( 'Multigroup partial reactions .....')
					for i in range (NpMTtg):
						print( 'MT = ',nMTprtg[i])
						#call groupmulti (insp,nMTprtg(i),igtype,1)
						#call groupmulti (insp,nMTprtg(i),igtype,2)
			else:
				print( 'Error', file = ofile_outRMINDD)
				print( 'wrong reaction index; please follow the list', file = ofile_outRMINDD)
				print( "1 = n,g" , file = ofile_outRMINDD)
				print( "2 = n,n" , file = ofile_outRMINDD)
				print( "3 = n,n'" , file = ofile_outRMINDD)
				print( "4 = n,xn" , file = ofile_outRMINDD)
				print( "5 = n,particle", file = ofile_outRMINDD)
				print( "6 = n,anything", file = ofile_outRMINDD)
				print( "7 = total", file = ofile_outRMINDD)
				print('', file = ofile_outRMINDD)
			
			print('------------------------------------------------', file = ofile_outRMINDD)
			print('', file = ofile_outRMINDD)
			print('The computed dpa and heating cross sections can', file = ofile_outRMINDD)
			print('be found in files:', file = ofile_outRMINDD)
			print('ndpa--.txt and nheat--.txt and', file = ofile_outRMINDD)
			print('ndpagrouped--.txt and nheatgrouped--.txt' , file = ofile_outRMINDD)
			print('for each reaction', file = ofile_outRMINDD)
			print('', file = ofile_outRMINDD)
			print('Total from CPO reactions is in: ....3001.txt', file = ofile_outRMINDD)
			print('Total from (n, xn) reactions is in: ....1601.txt', file = ofile_outRMINDD)
			print('Total from (n, anything) reactions is in: ....5001.txt', file = ofile_outRMINDD)
			
			stop_time = process_time()
			total_time = stop_time - start_time
			print('', file = ofile_outRMINDD)
			print( 'Total time taken:', file = ofile_outRMINDD)
			print(total_time, file = ofile_outRMINDD)

		# Second one is to calculate recoil energy distribution

		if (cmetric == 'RecedU'):

			ifile102 = open("tape01", 'r')

			print( '~~ RMINDD-RecedU ~~')
			print( ':Messages for you:', file = ofile_outRMINDD)
			print('--------------------', file = ofile_outRMINDD)
			print('', file = ofile_outRMINDD)
	  
			ifile102.readline()
			mat = ifile102.readline().split()[5][-4:] 				# read mat number in ENDF = 52
			matg = ifile_inpRMINDD.readline().split()[0] 			# given mat number in input = 50
			
			for i in range(3):
				ifile102.readline()
			data = ifile102.readline().split()
			iso = data[0] + data[1]
			ifile102.close()
			
			print( 'Evaluation on tape01', file = ofile_outRMINDD)
			print( iso, file = ofile_outRMINDD)
			print('', file = ofile_outRMINDD)
			
			# Check if given mat number matches with ENDF mat number
			
			##
			print(mat, matg, file = ofile_outRMINDD)
			
			if (mat != matg):
				print( 'Error', file = ofile_outRMINDD)
				print( 'Material does not exist on tape01', file = ofile_outRMINDD)
				print('', file = ofile_outRMINDD)
				break

			# nrct = Number of reactions to calculate
			# nrcta = Array of reaction indices for which to calculate [1=elastic,
			# 2=inelastic, 3=remaining particle thresholds, 4=(n,xn), 5=(n,g), 
			# 6=(n,anything-transmuted nuclei), 7=sum (total)]
			# nrg = Number of energy group limits
			# nbpoints = Number of energy points to add inside each group
			# igtype = which group? 
			# 1 = VITAMIN-J 175, 2 = 26 group, 3 = 33 group, 4 = 238 group,
			# 5 = 198 group,6 = 709 group, 7 = 640 group, 8 = 100 group,
			# 9 = 47 group, 10 = DLC-2 100 group, 11 = 229 group, 12 = 229 group.

		
			eliso = ifile_inpRMINDD.readline().strip('\n')
			nrct = int(ifile_inpRMINDD.readline().split()[0])
			if (nrct > 0):
				nrcta = [0]*nrct
				line = ifile_inpRMINDD.readline()
				for i in range(nrct):
					nrcta[i] = int(line.split()[i])
			if (nrct > 1):
				for value in nrcta:
					if (value == 7):
						print("Please do only partials or only total. In total all partials will also be done.")
						break
			line = ifile_inpRMINDD.readline()
			igtype = int(line.split()[0])
			nrg  = int(line.split()[1])
			nbpoints = int(line.split()[2])
			nre = nrg - 1
			insp = int(ifile_inpRMINDD.readline().split()[0])
			num_partial_reac_tosum = int(ifile_inpRMINDD.readline().split()[0])
			if (num_partial_reac_tosum > 0):
				partial_reac_tosum = [0]*num_partial_reac_tosum
				line = ifile_inpRMINDD.readline()
				for i in range(num_partial_reac_tosum):
					partial_reac_tosum[i] = int(line.split()[i])
			
			if (1 <= nrct and nrct <= 6):
				for i in range (nrct):
					RecedU.FINE_ENERGY_CALL_REAC(ofile_outRMINDD,insp,eliso,igtype,nrg,nbpoints,nrcta[i],nrcta)

			if (num_partial_reac_tosum > 0):
				ofile1001 = open('n-sum-partialsPKAspectra.txt', 'a')
				print(eliso, file = ofile1001)
				dsdt = numpy.zeros((nre,nre))
				dsdt = RecedU.ALLSUM (nrg,partial_reac_tosum)
				print('The sum of recoil nuclei energy spectra for given partial reactions', file = ofile1001)
				for it in range (nre):
					print (['{:.6E}'.format(dsdt[it][jt]) for jt in range (nre)], file = ofile1001)
				ofile1001.close()
			
			print('------------------------------------------------', file = ofile_outRMINDD)
			print('', file = ofile_outRMINDD)
			print('The computed PKA spectra can be found in files:', file = ofile_outRMINDD)
			print('PKA-MATRICES.txt -- each reaction', file = ofile_outRMINDD)
			print('n-allPKAspectra.txt -- sum total', file = ofile_outRMINDD)
		
			stop_time = process_time()
			total_time = stop_time - start_time
			print('', file = ofile_outRMINDD)
			print( 'Total time taken:', file = ofile_outRMINDD)
			print(total_time, file = ofile_outRMINDD)

ofile_outRMINDD.close()

ifile_inpRMINDD.close()
