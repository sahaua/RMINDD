#	>> Code: RMINDD (Radiation Matter Interaction from Nuclear Data and Damage)
#	>> Perform: Calcuation of Effects of Radiation in Matter from ENDF-6 files
#	>> Author: Dr. Uttiyoarnab Saha
#	>> Version and Date: 1.0 and 01/07/2022
#	>> Last modified: 01/07/2022, Kolkata
#	>> Update: 01/07/2022
#	>> Major changes: 
#
#=========================================================================================

import RecedU
import EngdepU
import numpy
#use TransmU

from time import process_time
	
ifile_inpRMINDD = open ("Input_RMINDD.txt", 'r')
	
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

			ofile51 = open ("Output_RadEMC-EngdepU.txt", 'a')
			ifile52 = open ("tape01", 'r')

			print( '~~ RadEMC-EngdepU ~~')
			print( ':Messages for you:', file = ofile51)
			print('--------------------', file = ofile51)
			print('', file = ofile51)
		
			ifile52.readline()
			mat = ifile52.readline().split()[5][-4:] 				# read mat number in ENDF = 52

			for i in range(3):
				ifile52.readline()
			data = ifile52.readline().split()
			iso = data[0] + data[1]
			ifile52.close()
	
			print( 'Evaluation on tape01', file = ofile51)
			print( iso, file = ofile51)
			print('', file = ofile51)
			matg = ifile_inpRMINDD.readline().split()[0] 			# given mat number in input = 50
			
			# Check if given mat number matches with ENDF mat number
			
			##
			print(mat, matg, file = ofile51)
			
			if (mat != matg):
				print( 'Error', file = ofile51)
				print( 'Material does not exist on tape01', file = ofile51)
				print('', file = ofile51)
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
			
				EngdepU.uqce (insp,nra,nreac,mdisp,Ed,bad,cad,mgyn,igtype)
			#call extractdata(nra(i))
			#call filerecds()
				if (NpMTtg > 0):
					print( 'Multigroup partial reactions .....')
					for i in range (NpMTtg):
						print( 'MT = ',nMTprtg[i])
						#call groupmulti (insp,nMTprtg(i),igtype,1)
						#call groupmulti (insp,nMTprtg(i),igtype,2)
			else:
				print( 'Error', file = ofile51)
				print( 'wrong reaction index; please follow the list', file = ofile51)
				print( "1 = n,g" , file = ofile51)
				print( "2 = n,n" , file = ofile51)
				print( "3 = n,n'" , file = ofile51)
				print( "4 = n,xn" , file = ofile51)
				print( "5 = n,particle", file = ofile51)
				print( "6 = n,anything", file = ofile51)
				print( "7 = total", file = ofile51)
				print('', file = ofile51)
			
			print('------------------------------------------------', file = ofile51)
			print('', file = ofile51)
			print('The computed dpa and heating cross sections can', file = ofile51)
			print('be found in files:', file = ofile51)
			print('ndpa--.txt and nheat--.txt and', file = ofile51)
			print('ndpagrouped--.txt and nheatgrouped--.txt' , file = ofile51)
			print('for each reaction', file = ofile51)
			print('', file = ofile51)
			print('Total from CPO reactions is in: ....3001.txt', file = ofile51)
			print('Total from (n, xn) reactions is in: ....1601.txt', file = ofile51)
			print('Total from (n, anything) reactions is in: ....5001.txt', file = ofile51)
			
			stop_time = process_time()
			total_time = stop_time - start_time
			print('', file = ofile51)
			print( 'Total time taken:', file = ofile51)
			print(total_time, file = ofile51)
			
			ofile51.close()

		# Second one is to calculate recoil energy distribution

		if (cmetric == 'RecedU'):
			
			ofile51 = open ("Output_RadEMC-RecedU.txt", 'a')
			ifile102 = open("tape01", 'r')

			print( '~~ RadEMC-RecedU ~~')
			print( ':Messages for you:', file = ofile51)
			print('--------------------', file = ofile51)
			print('', file = ofile51)
	  
			ifile102.readline()
			mat = ifile102.readline().split()[5][-4:] 				# read mat number in ENDF = 52
			matg = ifile_inpRMINDD.readline().split()[0] 			# given mat number in input = 50
			
			for i in range(3):
				ifile102.readline()
			data = ifile102.readline().split()
			iso = data[0] + data[1]
			ifile102.close()
			
			print( 'Evaluation on tape01', file = ofile51)
			print( iso, file = ofile51)
			print('', file = ofile51)
			
			# Check if given mat number matches with ENDF mat number
			
			##
			print(mat, matg, file = ofile51)
			
			if (mat != matg):
				print( 'Error', file = ofile51)
				print( 'Material does not exist on tape01', file = ofile51)
				print('', file = ofile51)
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
					RecedU.FINE_ENERGY_CALL_REAC(insp,eliso,igtype,nrg,nbpoints,nrcta[i],nrcta)

			if (num_partial_reac_tosum > 0):
				ofile1001 = open('n-sum-partialsPKAspectra.txt', 'a')
				print(eliso, file = ofile1001)
				dsdt = numpy.zeros((nre,nre))
				dsdt = RecedU.ALLSUM (nrg,partial_reac_tosum)
				print('The sum of recoil nuclei energy spectra for given partial reactions', file = ofile1001)
				for it in range (nre):
					print (['{:.6E}'.format(dsdt[it][jt]) for jt in range (nre)], file = ofile1001)
				ofile1001.close()
			
			print('------------------------------------------------', file = ofile51)
			print('', file = ofile51)
			print('The computed PKA spectra can be found in files:', file = ofile51)
			print('PKA-MATRICES.txt -- each reaction', file = ofile51)
			print('n-allPKAspectra.txt -- sum total', file = ofile51)
		
			stop_time = process_time()
			total_time = stop_time - start_time
			print('', file = ofile51)
			print( 'Total time taken:', file = ofile51)
			print(total_time, file = ofile51)

			ofile51.close()

	ifile_inpRMINDD.close()
