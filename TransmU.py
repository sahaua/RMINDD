''' 
Code: RMINDD - (Radiation-Matter Interaction and Damage calculation using Nuclear Data)
Perform: Calcuation of metrics of neutron radiation damage in a material using ENDF-6 files
Module: TransmU.py -- finds activation and gasproduction cross sections
Author: Uttiyoarnab Saha
Version and Date: 1.0 and 01/07/2022
Last modified: 01/07/2022, Kolkata
Update: 01/07/2022
Major changes: 

=========================================================================================
'''

import numpy
import math
import os

def ActivationGasProduction (ofile_outRMINDD,ifile_rawENDF6,ifile_preprocessedENDF6,eliso,en_group_type,input_n_spec,
transmgas_group_file, transmnucl_group_file, transmgas_point_file, transmnucl_MF5_point_file, 
transmnucl_net_group_file, NPt, Etu):

## Calculation of gas production and total activation cross section
## due to charged particle production reactions and (n,xn) and (n,g)
## reactions of neutron.

	tnEyld = numpy.zeros((200,200)); tnYld = numpy.zeros((200,200))
	tnNBTp = numpy.zeros((200,20)); tnINTrp = numpy.zeros((200,20))
	tnNyld = [0]*200; tnZp = [0]*200

	Eyld = numpy.zeros((5,200)); Yld = numpy.zeros((5,200))
	Nyld = [0]*5; iflMTtppr = [0]*5
	NBTp = numpy.zeros((5,20)); INTrp = numpy.zeros((5,20))
	NBTpp = [0]*20; NBTpd = [0]*20; NBTpt = [0]*20; NBTp3He = [0]*20; NBTpal = [0]*20
	INTrpp = [0]*20; INTrpd = [0]*20; INTrpt = [0]*20; INTrp3He = [0]*20; INTrpal = [0]*20

	## Assuming maximum 200 transmuted isotopes from MF6 MT5 and 
	## remaining are explicit transmutation reactions

	Zvaltrack = [0]*230; Avaltrack = [0]*230; cntrack = [0]*230

	ofile_transmgas_group = open(transmgas_group_file, 'w')
	ofile_transmnucl_group = open(transmnucl_group_file, 'w')
	ofile_transmgas_point = open(transmgas_point_file, 'w')
	ofile_transmnucl_MF5_point = open(transmnucl_MF5_point_file, 'w')
	ofile_transmnucl_net_group = open(transmnucl_net_group_file, 'w')

	print('E/(n,xp)/(n,xd)/(n,xt)/(n,x3He)/(n,xa)/(n,act.)', file = ofile_transmgas_group)
	print('Energy(eV) - Cross section(barns)', file = ofile_transmnucl_group)
	print('E/(n,xp)/(n,xd)/(n,xt)/(n,x3He)/(n,xa)/(n,act.)', file = ofile_transmgas_point)
	print('Energy(eV) - Cross section(barns)', file = ofile_transmnucl_MF5_point)
	print('Energy (eV) / Cross-section (barns)', file = ofile_transmnucl_net_group)

	nrab = 37
	MTnum = [0]*42; iflag = [0]*42

	MTnum = [5,11,16,17,22,23,24,25,28,29,30,32,33,34,35,36,37,41, \
	42,44,45,102,103,104,105,106,107,108,109,110,111,112,113,114, \
	115,116,117,203,204,205,206,207]
	
	## Read and Store the Values of Z and A of the target nucleus
	
	ifile_preprocessedENDF6.seek(0, 0)
	while True:
		line = ifile_preprocessedENDF6.readline()
		if (line == ''):
			break
		data = eachlineinfo(line)
		MAT = int(data[6]); MF = int(data[7]); MT  = int(data[8])
		if ( MT == 0 and MF == 0 ):
			line = ifile_preprocessedENDF6.readline()
			data = eachlineinfo(line)
			iflspace = 0
			for element in data:
				if (element == ''):
					iflspace = 1
					break
			if (iflspace == 0):
				ZAv = float(data[0]); AWRv = float(data[1]); L0 = int(data[2])
				L1 = int(data[3]); NKv = int(data[4]); L2 = int(data[5])
				MAT = int(data[6]); MF = int(data[7]); MT = int(data[8])

				Ztarget = int(ZAv)/1000
				Atarget =  mod(int(ZAv),1000)
		else:
			break

	## REMEMBER to SEND to this program the unique ENERGY ARRAY and total number of points
	## extracted from FILE 1 ........<<<<  
	
	print ( NPt, ' Total cross section energy points', file = ofile_outRMINDD)
	
	for i in range (NP):
		if (Etu[i] >= 20.0E+06):
			ifull = i
			break
	NP = ifull
	print ( NP, ' Energy points up to 20 MeV', file = ofile_outRMINDD)

	print ( 'Reading Data .....' )
	
	## Gas producion cross sections MT = 203, ...., 207	
	## may be given sometimes; nreac = 38 to 42 => these MTs
	
	
	sig5tot = [0]*NPt; sigparttot = [0]*NPt; sppoint = [0]*NPt; sdpoint = [0]*NPt
	strpoint = [0]*NPt; s3Hepoint = [0]*NPt; sapoint = [0]*NPt; sigtpoint = [0]*NPt
	
	for nreac in range(38, 43):
		ifile_preprocessedENDF6.seek(0, 0)
	
		iflag[nreac] = 0
		iflMTtppr[nreac-38] = 0
	
		MTi = MTnum[nreac]
	
		## extraction of cross sections
		NP1 = 0
		while (True):
			line = ifile_preprocessedENDF6.readline()
			if (line == ''):
				break
			data = eachlineinfo(line)
			MAT = int(data[6]); MF = int(data[7]); MT  = int(data[8])
			if (MT == 0 and MF != 0):
				line = ifile_preprocessedENDF6.readline()
				(ZAv,AWRv,L0,L1,NKv,L2,MAT,MF,MT) = line_type1_info(line)
			if (MAT != -1):
				if (MF == 3):
					if (MT == MTi):
						iflag[nreac] = 1
						iflMTtppr[nreac-38] = 1
						print('MT = ', MTi, '.....')
						print('Gas production MT found: ', MTi, file = ofile_outRMINDD)
						ZA = ZAv
						AWR = AWRv
						line = ifile_preprocessedENDF6.readline()
						data = eachlineinfo(line)
						QM = float(data[0]); QI =  float(data[1])
						NR = int(data[4]); NP1 = int(data[5])
						E1 = numpy.zeros(NP1); sig1 = numpy.zeros(NP1)
						ifile_preprocessedENDF6.readline()
						(E1, sig1) = line_type3_info(ifile_preprocessedENDF6,NP1,2)
			else:
				break
	
		if (iflag[nreac] == 1):
			if (MTi in [203, 204, 205, 206, 207]):
			if (MTi == 203 or MTi == 204 or MTi == 205 or MTi == 206 or MTi == 207):
				globals()[f'E{MTi}'] = []; globals()[f'sig{MTi}'] = []
				globals()[f'E{MTi}'] = numpy.append(globals()[f'E{MTi}'], E1)
				globals()[f'sig{MTi}'] = numpy.append(globals()[f'sig{MTi}'], sig1)

	iflp103 = 0
	iflp104 = 0
	iflp105 = 0
	iflp106 = 0
	iflp107 = 0

	## Calculation of gas production and activation from individual
	## CPO and other reactions; nreac = 1 to nrab(=37) => these MTs.

	for nreac in range(nrab):
		ifile_preprocessedENDF6.seek(0, 0)

		iflag[nreac] = 0
		MTi = MTnum[nreac]

		## extraction of cross sections
		NP1 = 0
		while(True):
			line = ifile_preprocessedENDF6.readline()
			if (line == ''):
				break
			data = eachlineinfo(line)
			MAT = int(data[6]); MF = int(data[7]); MT  = int(data[8])
			if (MT == 0 and MF != 0):
				line = ifile_preprocessedENDF6.readline()
				(ZAv,AWRv,L0,L1,NKv,L2,MAT,MF,MT) = line_type1_info(line)
			if (MAT != -1):
				if (MF == 3):
					if (MT == MTi):
						iflag[nreac] = 1
						if (MT == 103) iflp103 = 1
						if (MT == 104) iflp104 = 1
						if (MT == 105) iflp105 = 1
						if (MT == 106) iflp106 = 1
						if (MT == 107) iflp107 = 1
						print( 'MT = ', MTi, '.....')
						ZA = ZAv
						AWR = AWRv
						line = ifile_preprocessedENDF6.readline()
						data = eachlineinfo(line)
						QM = float(data[0]); QI =  float(data[1])
						NR = int(data[4]); NP1 = int(data[5])
						E1 = numpy.zeros(NP1); sig1 = numpy.zeros(NP1)
						ifile_preprocessedENDF6.readline()
						(E1, sig1) = line_type3_info(ifile_preprocessedENDF6,NP1,2)
			else:
				break

		if (iflag[nreac] == 1):
			print('', file = ofile_outRMINDD)
			print(f'Activation and/or Gas production MT = {MTi} given', file = ofile_outRMINDD)
			globals()[f'E{MTi}'] = []; globals()[f'sig{MTi}'] = []
			globals()[f'E{MTi}'] = numpy.append(globals()[f'E{MTi}'], E1)
			globals()[f'sig{MTi}'] = numpy.append(globals()[f'sig{MTi}'], sig1)
			
	## The calculation from MT = 600 to 849 is performed when 
	## MT = 103, ....., 107 are not given. These are performed in
	## the subroutine "adddiscnth".

	print('', file = ofile_outRMINDD)

	iflMTpr = 0

	if (iflp103 == 0):
		E103 = numpy.zeros(NPt); sig103 = numpy.zeros(NPt)
		sig103 = adddiscnth(iflMTpr, 103, Etu, NPt)
		if (iflMTpr == 1):
			print('Discrete + Continuum (n, p) data represented', file = ofile_outRMINDD)
			print("From discrete (n,p) .....")
			print('MT = ', 103, '.....')
			iflag[23] = 1

	iflMTpr = 0

	if (iflp104 == 0):
		E104 = numpy.zeros(NPt); sig104 = numpy.zeros(NPt)
		sig104 = adddiscnth(iflMTpr, 104, Etu, NPt)
		if (iflMTpr == 1):
			print('Discrete + Continuum (n, d) data represented', file = ofile_outRMINDD)
			print("From discrete (n,d) .....")
			print('MT = ', 104, '.....')
			iflag[24] = 1

	iflMTpr = 0

	if (iflp105 == 0):
		E105 = numpy.zeros(NPt); sig105 = numpy.zeros(NPt)
		sig105 = adddiscnth(iflMTpr, 105, Etu, NPt)
		if (iflMTpr == 1):
			print('Discrete + Continuum (n, t) data represented', file = ofile_outRMINDD)
			print("From discrete (n,t) .....")
			print('MT = ', 105, '.....')
			iflag[25] = 1

	iflMTpr = 0

	if (iflp106 == 0):
		E106 = numpy.zeros(NPt); sig106 = numpy.zeros(NPt)
		sig106 = adddiscnth(iflMTpr, 106, Etu, NPt)
		if (iflMTpr == 1):
			print('Discrete + Continuum (n, 3He) data represented', file = ofile_outRMINDD)
			print("From discrete (n,3He) .....")
			print('MT = ', 106, '.....')
			iflag[26] = 1

	iflMTpr = 0

	if (iflp107 == 0):
		E107 = numpy.zeros(NPt); sig107 = numpy.zeros(NPt)
		sig107 = adddiscnth(iflMTpr, 107, Etu, NPt)
		if (iflMTpr == 1):
			print('Discrete + Continuum (n, a) data represented', file = ofile_outRMINDD)
			print("From discrete (n,a) .....")
			print('MT = ', 107, '.....')
			iflag[27] = 1

	## The yields for the production of charged particles are collected
	## from MF = 6, MT = 5 in order to compute the respective gas production 
	## contributions from cross sections in MF = 3, MT = 5.

	(NBTp,INTrp,Nyld,Eyld,Yld) = gtYMf6Mt5 (tnNBTp,tnINTrp,tnNyld,tnEyld,tnYld,tnZp,itnt)

	Eyldp = []; Eyldd = []; Eyldtr = []; Eyld3He = []; EyldHe = []; Yldp = []
	Yldd = []; Yldtr = []; Yld3He = []; YldHe = []; NBTpp = []; NBTpd = []; NBTpt = []
	NBTp3He = []; NBTpal = []; INTrpp = []; INTrpd = []; INTrpt = []; INTrp3He = []; INTrpal = []

	Eyldp = numpy.append(Eyldp, Eyld[0][0:Nyld[0]])
	Eyldd = numpy.append(Eyldd, Eyld[1][0:Nyld[1]])
	Eyldtr = numpy.append(Eyldtr, Eyld[2][0:Nyld[2]])
	Eyld3He = numpy.append(Eyld3He, Eyld[3][0:Nyld[3]])
	EyldHe = numpy.append(EyldHe, Eyld[4][0:Nyld[4]])

	Yldp = numpy.append(Yldp, Yld[0][0:Nyld[0]])
	Yldd = numpy.append(Yldd, Yld[1][0:Nyld[1]])
	Yldtr = numpy.append(Yldtr, Yld[2][0:Nyld[2]])
	Yld3He = numpy.append(Yld3He, Yld[3][0:Nyld[3]])
	YldHe = numpy.append(YldHe, Yld[4][0:Nyld[4]])

	NBTpp = numpy.append(NBTpp, NBTp[0][:])
	NBTpd = numpy.append(NBTpd, NBTp[1][:])
	NBTpt = numpy.append(NBTpt, NBTp[2][:])
	NBTp3He = numpy.append(NBTp3He, NBTp[3][:])
	NBTpal = numpy.append(NBTpal, NBTp[4][:])
	
	INTrpp = numpy.append(INTrpp, INTrp[0][:])
	INTrpd = numpy.append(INTrpd, INTrp[1][:])
	INTrpt = numpy.append(INTrpt, INTrp[2][:])
	INTrp3He = numpy.append(INTrp3He, INTrp[3][:])
	INTrpal = numpy.append(INTrpal, INTrp[4][:])
	
	print('', file = ofile_outRMINDD)
	print('Multigroup activation cross sections asked in:', file = ofile_outRMINDD)
	print('', file = ofile_outRMINDD)
	
	if (en_group_type == 0):
		print( '0 -- User-defined energy groups', file = ofile_outRMINDD)
		ifile = open('Energy-GroupLimits.txt', 'r')
		Ngl = int(ifile.readline().split()[0])
		for i in reversed(range(Ngl)):
			Eg[i] = float(ifile.readline().split()[0])
		ifile.close()
	if (en_group_type == 1):
		(Eg,Ngl) = engrp1()
		print( '1 -- VITAMIN-J 175 groups', file = ofile_outRMINDD)
	if (en_group_type == 2):
		(Eg,Ngl) = engrp2()
		print( '2 -- 238 groups', file = ofile_outRMINDD)
	if (en_group_type == 3):
		(Eg,Ngl) = engrp3()
		print( '3 -- 198 groups', file = ofile_outRMINDD)
	if (en_group_type == 4):
		(Eg,Ngl) = engrp4()
		print( '4 -- 33 groups', file = ofile_outRMINDD)
	if (en_group_type == 5):
		(Eg,Ngl) = engrp5()
		print( '5 -- 26 groups', file = ofile_outRMINDD)
	if (en_group_type == 6):
		(Eg,Ngl) = engrp6()
		print( '6 -- DLC-2 100 groups', file = ofile_outRMINDD)
	if (en_group_type == 7):
		(Eg,Ngl) = engrp7()
		print( '7 -- 198 groups to 229 groups (200 MeV)', file = ofile_outRMINDD)
	if (en_group_type == 8):
		(Eg,Ngl) = engrp8()
		print( '8 -- 198 groups to 229 groups (150 MeV)', file = ofile_outRMINDD)
	if (en_group_type == 9):
		(Eg,Ngl) = engrp9()
		print( '9 -- 616 groups (from DEMO-HCPB-FW)', file = ofile_outRMINDD)

	print('', file = ofile_outRMINDD)
	print('Activation cross sections can be found in file:', file = ofile_outRMINDD)
	print('GasProDatatest', file = ofile_outRMINDD)
	print('', file = ofile_outRMINDD)

	## Multigrouping the required cross sections in the given energy
	## group structure and then adding into separate collections for 
	## particular gas production species and activation.

	print('Grouping: ', Ngl-1, 'groups .....')

	gsp = numpy.zeros(Ngl); gsd = numpy.zeros(Ngl); gstr = numpy.zeros(Ngl); gs3He = numpy.zeros(Ngl) 
	gsa = numpy.zeros(Ngl)
	gsig5 = []
	gsigt = numpy.zeros(Ngl)
	tnYldg = numpy.zeros(Ngl)
	sigtrack = numpy.zeros((Ngl,230))

	Yldpg = terpol(NBTpp,INTrpp,Nyld[0],Eyldp,Yldp,Ngl,Eg)
	Ylddg = terpol(NBTpd,INTrpd,Nyld[1],Eyldd,Yldd,Ngl,Eg)
	Yldtrg = terpol(NBTpt,INTrpt,Nyld[2],Eyldtr,Yldtr,Ngl,Eg)
	Yld3Heg = terpol(NBTp3He,INTrp3He,Nyld[3],Eyld3He,Yld3He,Ngl,Eg)
	YldHeg = terpol(NBTpal,INTrpal,Nyld[4],EyldHe,YldHe,Ngl,Eg)

	for nreac in range(42):
		MTi = MTnum(nreac)
		if (iflag[nreac] == 1):
			NPtg = len(globals()[f'E{MTi}'])
			if (MTi == 5):
				sig5tot = terpolAPR(2,NPtg,globals()[f'E{MTi}'],globals()[f'sig{MTi}'],NPt,Etu)
			else:
				sigparttot = terpolAPR(2,NPtg,globals()[f'E{MTi}'],globals()[f'sig{MTi}'],NPt,Etu)

			## Cross section gets multigrouped based on n-spectrum here .....
			gsig = groupmulti (input_n_spec,globals()[f'E{MTi}'],globals()[f'sig{MTi}'],NPtg,Eg,Ngl) 

			## cross sections from lumped MT = 5 data
			if (MTi == 5):
				gsig5 = numpy.append(gsig5, gsig)
				for i in range(Ngl):
					gsp[i] = gsp[i] + Yldpg[i]*gsig[i]
					gsd[i] = gsd[i] + Ylddg[i]*gsig[i]
					gstr[i] = gstr[i] + Yldtrg[i]*gsig[i]
					gs3He[i] = gs3He[i] + Yld3Heg[i]*gsig[i]
					gsa[i] = gsa(i) + YldHeg[i]*gsig[i]
					gsigt[i] = gsigt[i] + gsig[i]

				## segregations of point cross sections
				Yldpg = terpol(NBTpp,INTrpp,Nyld[0],Eyldp,Yldp,NPt,Etu) 
				Ylddg = terpol(NBTpd,INTrpd,Nyld[1],Eyldd,Yldd,NPt,Etu)
				Yldtrg = terpol(NBTpt,INTrpt,Nyld[2],Eyldtr,Yldtr,NPt,Etu)
				Yld3Heg = terpol(NBTp3He,INTrp3He,Nyld[3],Eyld3He,Yld3He,NPt,Etu)
				YldHeg = terpol(NBTpal,INTrpal,Nyld[4],EyldHe,YldHe,NPt,Etu)

				for i in range(NPt):
					sppoint[i] = sppoint[i] + Yldpg[i]*sig5tot[i]
					sdpoint[i] = sdpoint[i] + Ylddg[i]*sig5tot[i]
					strpoint[i] = strpoint[i] + Yldtrg[i]*sig5tot[i]
					s3Hepoint[i] = s3Hepoint[i] + Yld3Heg[i]*sig5tot[i]
					sapoint[i] = sapoint[i] + YldHeg[i]*sig5tot[i]
					sigtpoint[i] = sigtpoint[i] + sig5tot[i]
			## enf of if MTi = 5
			
			if (iflMTtppr[0] == 0):
				if (MTi==28 or MTi==41 or MTi==42 or MTi==44 or MTi==45 or \
				MTi==103 or MTi==111 or MTi==112 or MTi==115 or MTi==116):
					for i in range(Ngl):
						gsp[i] = gsp[i] + gsig[i]
						gsigt[i] = gsigt[i] + gsig[i]
					for i in range(NPt):
						sppoint[i] = sppoint[i] + sigparttot[i]
						sigtpoint[i] = sigtpoint[i] + sigparttot[i]
					end do

			if (MTi == 203):
				for i in range(Ngl):
					gsp[i] = gsig[i]
					gsigt[i] = gsigt[i] + gsig[i]
				for i in range(NPt):
					sppoint[i] = sigparttot[i]
					sigtpoint[i] = sigtpoint[i] + sigparttot[i]

			if (iflMTtppr[1] == 0):
				if (MTi==11 or MTi==32 or MTi==35 or MTi==104 or \
				MTi==114 or MTi==115 or MTi==117):
					for i in range(Ngl):
						gsd[i] = gsd(i) + gsig[i]
						gsigt[i] = gsigt[i] + gsig[i]
					for i in range(NPt):
						sdpoint[i] = sdpoint[i] + sigparttot[i]
						sigtpoint[i] = sigtpoint[i] + sigparttot[i]

			if (MTi == 204):
				for i in range(Ngl):
					gsd[i] = gsig[i]
					gsigt[i] = gsigt[i] + gsig[i]
				for i in range(NPt):
					sdpoint[i] = sigparttot[i]
					sigtpoint[i] = sigtpoint[i] + sigparttot[i]

			if (iflMTtppr[2] == 0):
				if (MTi==33 or MTi==36 or MTi==105 or MTi==113 or MTi==116):
					for i in range(Ngl):
						gstr[i] = gstr[i] + gsig[i]
						gsigt[i] = gsigt[i] + gsig[i]
					for i in range(NPt):
						strpoint[i] = strpoint[i] + sigparttot[i]
						sigtpoint[i] = sigtpoint[i] + sigparttot[i]

			if (MTi == 205):
				for i in range(Ngl):
					gstr[i] = gsig[i]
					gsigt[i] = gsigt[i] + gsig[i]
				for i in range(NPt):
					strpoint[i] = sigparttot[i]
					sigtpoint[i] = sigtpoint[i] + sigparttot[i]

			if (iflMTtppr[3] == 0):
				if (MTi==34 or MTi==106):
					for i in range(Ngl):
						gs3He[i] = gs3He[i] + gsig[i]
						gsigt[i] = gsigt[i] + gsig[i]
					for i in range(NPt):
						s3Hepoint[i] = s3Hepoint[i] + sigparttot[i]
						sigtpoint[i] = sigtpoint[i] + sigparttot[i]

			if (MTi == 206):
				for i in range(Ngl):
					gs3He[i] = gsig[i]
					gsigt[i] = gsigt[i] + gsig[i]
				for i in range(NPt):
					s3Hepoint[i] = sigparttot[i]
					sigtpoint[i] = sigtpoint[i] + sigparttot[i]
				
			if (iflMTtppr[4] == 0):
				if (MTi==22 or MTi==23 or MTi==24 or MTi==25 or MTi==29 or \
				MTi==30 or MTi==35 or MTi==36 or MTi==45 or MTi==107 or \
				MTi==108 or MTi==109 or MTi==112 or MTi==113 or MTi==114 or	MTi==117):
					for i in range(Ngl):
						gsa[i] = gsa[i] + gsig[i]
						gsigt[i] = gsigt[i] + gsig[i]
					for i in range(NPt):
						sapoint[i] = sapoint[i] + sigparttot[i]
						sigtpoint[i] = sigtpoint[i] + sigparttot[i]

			if (MTi == 207):
				for i in range(Ngl):
					gsa[i] = gsig[i]
					gsigt[i] = gsigt[i] + gsig[i]
				for i in range(NPt):
					sapoint[i] = sigparttot[i]
					sigtpoint[i] = sigtpoint[i] + sigparttot[i]
					end do

			if (MTi == 16 or MTi == 17 or MTi == 37 or MTi == 102):
				for i in range(Ngl):
					gsigt[i] = gsigt[i] + gsig[i]
				for i in range(NPt):
					sigtpoint[i] = sigtpoint[i] + sigparttot[i]

			print('MT = ', MTi)

			if (MTi != 5):
				if (MTi==11 or MTi==33 or MTi==42):
					Zval = Ztarget-1 
					Aval = Atarget-3
					Zvaltrack[0] = Zval
					Avaltrack[0] = Aval
					cntrack[0] = cntrack[0]+1
					for k in range(Ngl):
						sigtrack[k][0] = sigtrack[k][0] + gsig[k]
				if (MTi == 16):
					Zval = Ztarget
					Aval = Atarget-1
					Zvaltrack[1] = Zval
					Avaltrack[1] = Aval
					cntrack[1] = cntrack[1]+1
					for k in range(Ngl):
						sigtrack[k][1] = sigtrack[k][1] + gsig[k]
				if (MTi == 17):
					Zval = Ztarget
					Aval = Atarget-2
					Zvaltrack[2] = Zval
					Avaltrack[2] = Aval
					cntrack[2] = cntrack[2]+1
					for k in range(Ngl):
						sigtrack[k][2] = sigtrack[k][2] + gsig[k]
				if (MTi == 22):
					Zval = Ztarget-2
					Aval = Atarget-4
					Zvaltrack[3] = Zval
					Avaltrack[3] = Aval
					cntrack[3] = cntrack[3]+1
					for k in range(Ngl):
						sigtrack[k][3] = sigtrack[k][3] + gsig[k]
				if (MTi == 23):
					Zval = Ztarget-6
					Aval = Atarget-12
					Zvaltrack[4] = Zval
					Avaltrack[4] = Aval
					cntrack[4] = cntrack[4]+1
					for k in range(Ngl):
						sigtrack[k][4] = sigtrack[k][4] + gsig[k]
				if (MTi == 24):
					Zval = Ztarget-2
					Aval = Atarget-5
					Zvaltrack[5] = Zval
					Avaltrack[5] = Aval
					cntrack[5] = cntrack[5]+1
					for k in range(Ngl):
						sigtrack[k][5] = sigtrack[k][5] + gsig[k]
				if (MTi == 25):
					Zval = Ztarget-2
					Aval = Atarget-6
					Zvaltrack[6] = Zval
					Avaltrack[6] = Aval
					cntrack[6] = cntrack[6]+1
					for k in range(Ngl):
						sigtrack[k][6] = sigtrack[k][6] + gsig[k]
				if (MTi == 28 or MTi == 104):
					Zval = Ztarget-1
					Aval = Atarget-1
					Zvaltrack[7] = Zval
					Avaltrack[7] = Aval
					cntrack[7] = cntrack[7]+1
					for k in range(Ngl):
						sigtrack[k][7] = sigtrack[k][7] + gsig[k]
				if (MTi == 29):
					Zval = Ztarget-4
					Aval = Atarget-8
					Zvaltrack[8] = Zval
					Avaltrack[8] = Aval
					cntrack[8] = cntrack[8]+1
					for k in range(Ngl):
						sigtrack[k][8] = sigtrack[k][8] + gsig[k]
				if (MTi == 30): 
					Zval = Ztarget-4
					Aval = Atarget-9
					Zvaltrack[9] = Zval
					Avaltrack[9] = Aval
					cntrack[9] = cntrack[9]+1
					for k in range(Ngl):
						sigtrack[k][9] = sigtrack[k][9] + gsig[k]
				if (MTi == 32 or MTi == 41 or MTi == 105):
					Zval = Ztarget-1
					Aval = Atarget-2
					Zvaltrack[10] = Zval
					Avaltrack[10] = Aval
					cntrack[10] = cntrack[10]+1
					for k in range(Ngl):
						sigtrack[k][10] = sigtrack[k][10] + gsig[k]
				if (MTi == 34 or MTi == 107 or MTi == 116):
					Zval = Ztarget-2
					Aval = Atarget-3
					Zvaltrack[11] = Zval
					Avaltrack[11] = Aval
					cntrack[11] = cntrack[11]+1
					for k in range(Ngl):
						sigtrack[k][11] = sigtrack[k][11] + gsig[k]
				if (MTi == 35 or MTi == 113):
					Zval = Ztarget-5
					Aval = Atarget-10
					Zvaltrack[12] = Zval
					Avaltrack[12] = Aval
					cntrack[12] = cntrack[12]+1
					for k in range(Ngl):
						sigtrack[k][12] = sigtrack[k][12] + gsig[k]
				if (MTi == 36):
					Zval = Ztarget-5
					Aval = Atarget-11
					Zvaltrack[13] = Zval
					Avaltrack[13] = Aval
					cntrack[13] = cntrack[13]+1
					for k in range(Ngl):
						sigtrack[k][13] = sigtrack[k][13] + gsig[k]
				if (MTi == 37):
					Zval = Ztarget
					Aval = Atarget-3
					Zvaltrack[14] = Zval
					Avaltrack[14] = Aval
					cntrack[14] = cntrack[14]+1
					for k in range(Ngl):
						sigtrack[k][14] = sigtrack[k][14] + gsig[k]
				if (MTi == 44 or MTi == 106 or MTi == 115):
					Zval = Ztarget-2
					Aval = Atarget-2
					Zvaltrack[15] = Zval
					Avaltrack[15] = Aval
					cntrack[15] = cntrack[15]+1
					for k in range(Ngl):
						sigtrack[k][15] = sigtrack[k][15] + gsig[k]
				if (MTi == 45 or MTi == 117):
					Zval = Ztarget-3
					Aval = Atarget-5
					Zvaltrack[16] = Zval
					Avaltrack[16] = Aval
					cntrack[16] = cntrack[16]+1
					for k in range(Ngl):
						sigtrack[k][16] = sigtrack[k][16] + gsig[k]
				if (MTi == 102):
					Zval = Ztarget
					Aval = Atarget+1
					Zvaltrack[17] = Zval
					Avaltrack[17] = Aval
					cntrack[17] = cntrack[17]+1
					for k in range(Ngl):
						sigtrack[k][17] = sigtrack[k][17] + gsig[k]
				if (MTi == 103):
					Zval = Ztarget-1
					Aval = Atarget
					Zvaltrack[18] = Zval
					Avaltrack[18] = Aval
					cntrack[18] = cntrack[18]+1
					for k in range(Ngl):
						sigtrack[k][18] = sigtrack[k][18] + gsig[k]
				if (MTi == 108):
					Zval = Ztarget-4
					Aval = Atarget-7
					Zvaltrack[19] = Zval
					Avaltrack[19] = Aval
					cntrack[19] = cntrack[19]+1
					for k in range(Ngl):
						sigtrack[k][19] = sigtrack[k][19] + gsig[k]
				if (MTi == 109):
					Zval = Ztarget-6
					Aval = Atarget-11
					Zvaltrack[20] = Zval
					Avaltrack[20] = Aval
					cntrack[20] = cntrack[20]+1
					for k in range(Ngl):
						sigtrack[k][20] = sigtrack[k][20] + gsig[k]
				if (MTi == 111):
					Zval = Ztarget-2
					Aval = Atarget-1
					Zvaltrack[21] = Zval
					Avaltrack[21] = Aval
					cntrack[21] = cntrack[21]+1
					for k in range(Ngl):
						sigtrack[k][21] = sigtrack[k][21] + gsig[k]
				if (MTi == 112):
					Zval = Ztarget-3
					Aval = Atarget-4
					Zvaltrack[22] = Zval
					Avaltrack[22] = Aval
					cntrack[22] = cntrack[22]+1
					for k in range(Ngl):
						sigtrack[k][22] = sigtrack[k][22] + gsig[k]
				if (MTi == 114):
					Zval = Ztarget-5
					Aval = Atarget-9
					Zvaltrack[23] = Zval
					Avaltrack[23] = Aval
					cntrack[23] = cntrack[23]+1
					for k in range(Ngl):
						sigtrack[k][23] = sigtrack[k][23] + gsig[k]
				print('MT = ', MTi, 'Z = ', Zval, 'A = ', Aval, file = ofile_transmnucl_group)
				for k in range(Ngl):
					print(Eg[k], gsig[k], file = ofile_transmnucl_group)
				print('-----------------', file = ofile_transmnucl_group)

	print(eliso, Ngl, file = ofile_transmgas_group)
	for k in range(Ngl):
		print(Eg[k],gsp[k],gsd[k],gstr[k],gs3He[k],gsa[k],gsigt[k], file = ofile_transmgas_group)
	print('-------------------', file = ofile_transmgas_group)

	## Print the Point Cross Sections (Energies corresponding
	## to the MF=3 MT=1) of Production of Light Charged Particles
	
	print(eliso, NPt, file = ofile_transmgas_point)
	for k in range(NPt):
		print(Etu[k],sppoint[k],sdpoint[k],strpoint[k],s3Hepoint[k],sapoint[k],sigtpoint[k], file = ofile_transmgas_point)
	print('-------------------', file = ofile_transmgas_point)

	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	## Find and Print the Multigrouped cross sections of Transmuted Nuclides

	for itn in range(itnt):
		itnZ = tnZp[itn]/1000
		itnA =  math.fmod(tnZp[itn],1000)

		tnYldg = terpol(tnNBTp[itn][0:19],tnINTrp[itn][0:19],tnNyld[itn],tnEyld[itn][0:199],tnYld[itn][0:199],Ngl,Eg)

		iflsametrnsN = 0
		for i in range(230):
			if (Zvaltrack[i]==itnZ and Avaltrack[i]==itnA):
				cntrack[i] = cntrack[i] + 1
				iflsametrnsN = 1
				for k in range(Ngl):
					sigtrack[k][i] = sigtrack[k][i] + tnYldg[k]*gsig5[k]
				break

		if (iflsametrnsN == 0):
			for i in range(230):
				if (Zvaltrack[i]==0 and Avaltrack[i]==0):
					Zvaltrack[i] = itnZ
					Avaltrack[i] = itnA
					cntrack[i] = cntrack[i] + 1
					for k in range(Ngl):
						sigtrack[k][i] = sigtrack[k][i] + tnYldg[k]*gsig5[k]
					break

		print(itn,eliso,itnZ,itnA, file = ofile_transmnucl_group)
		for k in range(Ngl):
			print(Eg[k], tnYldg[k]*gsig5[k], file = ofile_transmnucl_group)
		print('-----------------', file = ofile_transmnucl_group)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	print('List of Transmuted Nuclides Produced and Number of Times Each is Found', file = ofile_transmnucl_group)
	print('Sl. No.Z/A/Counts', file = ofile_transmnucl_group)

	## Find Net Number of Transmuted Nuclides
	for i in reversed(range(230)):
		if (Zvaltrack[i] != 0 and Avaltrack[i] != 0 and cntrack[i] != 0):
			itrnstot = i
			break

	for i in range(itrnstot+1):
		print(i, Zvaltrack[i], Avaltrack[i], cntrack[i], file = ofile_transmnucl_group)
	print('------------------'

	for i in range(itrnstot):
		print( i, eliso, Zvaltrack[i], Avaltrack[i], file = ofile_transmnucl_net_group)
		for k in range(Ngl):
			print(Eg[k], sigtrack[k][i], file = ofile_transmnucl_net_group)
		print('------------------', file = ofile_transmnucl_net_group)

	## Find and Print the Point Cross Sections (Energies corresponding
	## to the MF=3 MT=1) of Transmuted Nuclides

	for itn in range(itnt):
		itnZ = tnZp[itn]/1000
		itnA =  math.fmod(tnZp[itn],1000)
		tnYldg = terpol(tnNBTp[itn][0:19],tnINTrp[itn][0:19],tnNyld[itn], \
				tnEyld[itn][0:199],tnYld[itn][0:199],NPt,Etu)
		print(itn, eliso, itnZ, itnA, file = ofile_transmnucl_MF5_point)
		for k in range(NPt):
			if (tnYldg[k] != 0):
				kstart = k-1
				break
		for  k in range(kstart, NPt):
			print(Etu[k], tnYldg[k]*sig5tot[k], file = ofile_transmnucl_MF5_point)
		print('------------------', file = ofile_transmnucl_MF5_point)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	ofile_transmgas_group.close()
	ofile_transmnucl_group.close()
	ofile_transmgas_point.close()
	ofile_transmnucl_MF5_point.close()
	ofile_transmnucl_net_group.close()
		
##=======================================================================


## To energy multigroup cross sections depending on neutron spectrum
	
def groupmulti (input_n_spec,E1p,sig,NPtg,Eg,Ngl):
	gsig = numpy.zeros(Ngl)

	if (input_n_spec == 1):
		fi = numpy.zeros(Ngl)
		ifile407 = open('NeutronSpectrum.txt', 'r')
		nre = int(ifile407.readline().split()[0])
		for i in reversed(range(nre)):
			fi[i] = float(ifile407.readline().split()[0])
		ifile407.close()

	ifg = 0
	for i in range(Ngl-1):
		if (Eg[i]<=E1p[0] and E1p[0]<=Eg[i+1]):
			ifg = i
			break

	for i in range (ifg, Ngl-1):
		Eg1 = Eg[i]
		Eg2 = Eg[i+1]
		Nsect = 10000
		h = (Eg2-Eg1)/Nsect
		bcs = 0; denominator = 0; deno1 = 0; deno2 = 0
		bcs1 = 0; bcs2 = 0
		for j in range (Nsect):
			t = Eg1 + j*h
			if (t < 0.1):
				L = 1
			if (0.1 <= t and t < 820.3e+03):
				L = 2
			if (t >= 820.3e+03):
				L = 3
			if (input_n_spec == 0):
				flux1 = spectrum (t,L)
			if (input_n_spec == 1):
				flux1 = numpy.interp(t, Eg, fi)

			bcs1 = numpy.interp(t+h, Etu, sig) * flux1
			deno1 = flux1

			if (t+h < 0.1):
				L = 1
			if (0.1 <= t+h and t+h < 820.3e+03):
				L = 2
			if (t+h >= 820.3e+03):
				L = 3
			if (input_n_spec == 0):
				flux2 = spectrum (t+h,L)
			if (input_n_spec == 1):
				flux2 = numpy.interp(t+h, Eg, fi)

			bcs2 = numpy.interp(t+h, Etu, sig) * flux2
			deno2 = flux2

			bcs = bcs + (h/2)*(bcs1 + bcs2)
			denominator = denominator + (h/2)*(deno1 + deno2)

		if (bcs != 0 and denominator != 0):
			gsig[i] = bcs/denominator

	return(gsig)

##=======================================================================

## To interpolate into required energy points using the specific
## scheme in the range.
	
def terpol (NRin,Intrin,n1,E1,Y1,n2,E2):
	for i in range(n2):
		for j in range(n1):
			if (E2[i] == E1[j]):
				Y2[i] = Y1[j]
				break
			if (E1[j] < E2[i] and E2[i] < E1[j+1]):
				diff1 = E2[i] - E1[j]
				diff2 = E2[i] - E2[i-1]

				for k in range(20):
					if (j <= NRin[k]):
						intflg = Intrin[k]
						break

				if (diff1 <= diff2):
					Y2[i] = terpolin(intflg,E2[i],E1[j],E1[j+1],Y1[j],Y1[j+1])

				if (diff2 < diff1):
					Y2[i] = terpolin(intflg,E2[i],E2[i-1],E1[j+1],Y2[i-1],Y1[j+1])
				break
	return(Y2)
## ===============================================================
## Schemes of interpolation
	
def terpolin (intflg,x,x1,x2,y1,y2):
	interpscheme = intflg
	if (interpscheme == 1 or interpscheme == 11 or interpscheme == 21):
		y = y1
	if (interpscheme == 2 or interpscheme == 12 or interpscheme == 22):
		y = y1+(y2-y1)*(x-x1)/(x2-x1)
	if (interpscheme == 3 or interpscheme == 13 or interpscheme == 23):
		y = y1+(y2-y1)*(log(x/x1))/log(x2/x1)
	if (interpscheme == 4 or interpscheme == 14 or interpscheme == 24):
		y = y1*exp((x-x1)*log(y2/y1)/(x2-x1))
	if (interpscheme == 5 or interpscheme == 15 or interpscheme == 25):
		y = y1*exp(log(x/x1)*log(y2/y1)/log(x2/x1))

	return(y)

## ===============================================================

## Interpolation of Arrays via Particular Rule 

def terpolAPR (iprule,n1,x1,y1,n2,x2):
	for i in range(n2):
		for j in range(n1):
			if (x2[i] == x1[j]):
				y2[i] = y1[j]
				break
			if (x1[j] < x2[i] and x2[i] < x1[j+1]):
				diff1 = x2[i] - x1[j]		
				diff2 = x2[i] - x2[i-1]
				if (diff1 <= diff2):
					y2[i] = terpolin(iprule,x2[i],x1[j],x1[j+1],y1[j],y1[j+1])
				if (diff2 < diff1):
					y2[i] = terpolin(iprule,x2[i],x2[i-1],x1[j+1],y2[i-1],y1[j+1])
				break
	return(y2)

#=======Get the yields of particle production=======*
	
# Collect the yields of five species of light charged particles
# from File 6 MT = 5.
	
def gtYMf6Mt5():
	Eyld = [[0]*200]*5; Yld = [[0]*200]*5
	iZp = [0]*5; iAp = [0]*5; Nyld = [0]*5
	NBTp = [[0]*20]*5; INTrp = [[0]*20]*5

	iZp[0] = 1001
	iZp[1] = 1002
	iZp[2] = 1003
	iZp[3] = 2003
	iZp[4] = 2004

	iAp[0] = 1
	iAp[1] = 2
	iAp[2] = 3
	iAp[3] = 3
	iAp[4] = 4

	ifile_rawENDF6.seek(0, 0)
	ifile = ifile_rawENDF6

	while (True):
		line = ifile.readline()
		data = eachlineinfo(line)
		MAT = int(data[6]); MF = int(data[7]); MT  = int(data[8])
		if (MAT == -1):
			ifl6mt5 = 0
			break
		if (MT == 0):
			if (MF == 6):
				line = ifile.readline()
				data = eachlineinfo(line)
				ZAv = float(data[0]); AWRv = float(data[1]); l1 = int(data[2])
				LCT = int(data[3]); NKv = int(data[4]); l2 = int(data[5])
				MAT = int(data[6]); MF = int(data[7]); MT = int(data[8])

			if (MF < 6):
				line = ifile.readline()
				data = eachlineinfo(line)
				ZAv = float(data[0]); AWRv = float(data[1]); l1 = int(data[2])
				LCT = int(data[3]); NKv = int(data[4]); l2 = int(data[5])
				MAT = int(data[6]); MF = int(data[7]); MT = int(data[8])

		if (MAT != -1):
			if (MF == 6):
				if (MT == 5):
					ifl6mt5 = 1	
					NK = NKv
					for NSS in range (NK):
						line = ifile.readline()
						data = eachlineinfo(line)
						ZAP = float(data[0]); AWP = float(data[1]); LIP = int(data[2])
						LAW = int(data[3]); NR6 = int(data[4]); NP6 = int(data[5])
						MAT = int(data[6]); MF = int(data[7]); MT = int(data[8])

						iflgp = 0
						if(int(ZAP) == iZp[0] or int(ZAP) == iZp[1] or int(ZAP) == iZp[2] \
						or int(ZAP) == iZp[3] or int(ZAP) == iZp[4]):
							iflgp = 1
						
						if (int(ZAP) == iZp[0] and int(ceiling(AWP)) == iAp[0]):
							ip = 0
						if (int(ZAP) == iZp[1] and int(ceiling(AWP)) == iAp[1]):
							ip = 1
						if (int(ZAP) == iZp[2] and int(ceiling(AWP)) == iAp[2]):
							ip = 2
						if (int(ZAP) == iZp[3] and int(ceiling(AWP)) == iAp[4]):
							ip = 3
						if (int(ZAP) == iZp[4] and int(ceiling(AWP)) == iAp[5]):
							ip = 4
		
						if (iflgp == 0):
							N = 0
							while (N < NR6):
								line = ifile.readline()
								data = eachlineinfo(line)
								for i in range(0,6,2):
									NBT = int(data[i])
									INTr = int(data[i+1])
									N += 1
							N = 0
							while (N < NP6):
								line = ifile.readline()
								data = eachlineinfo(line)
								for i in range(0,6,2):
									E = float(data[i])
									Y = float(data[i+1])
									N += 1

						if (iflgp == 1):
							Nyld[ip] = NP6
							N = 0
							while (N < NR6):
								line = ifile.readline()
								data = eachlineinfo(line)
								for i in range(0,6,2):
									NBTp[ip][N] = int(data[i])
									INTrp[ip][N] = int(data[i+1])
									N += 1
							N = 0
							while (N < NP6):
								line = ifile.readline()
								data = eachlineinfo(line)
								for i in range(0,6,2):
									Eyld[ip][N] = float(data[i])
									Yld[ip][N] = float(data[i+1])
									N += 1

						if (LAW != 0):	# added this condition for calculating from JEFF-3.3
							line = ifile.readline()
							data = eachlineinfo(line)
							c1 = float(data[0]); c2 = float(data[1]); l3 = int(data[2])
							l4 = int(data[3]); NR6 = int(data[4]); NE6 = int(data[5])
							MAT = int(data[6]); MF = int(data[7]); MT = int(data[8])

							N = 0
							while (N < NR6):
								line = ifile.readline()
								data = eachlineinfo(line)
								for i in range(0,6,2):
									NBT = int(data[i])
									INTr = int(data[i+1])
									N += 1
							for i in range(NE6):
								line = ifile.readline()
								data = eachlineinfo(line)
								c1 = float(data[0]); En = float(data[1]); ND = int(data[2])
								NA = int(data[3]); NW = int(data[4]); NEP = int(data[5])
								MAT = int(data[6]); MF = int(data[7]); MT = int(data[8])
								NS = int(data[9])

								J1 = 0
								while (J1 < NW):
									line = ifile.readline()
									data = eachlineinfo(line)
									for k in range(6):
										Bvall[J1] = float(data[k])
										J1 += 1

			else:
				break

	return (NBTp,INTrp,Nyld,Eyld,Yld)
	
#=======To find for availble MT=======*
	
	# The presence of a particular reaction cross section (MT in MF3) is
	# searched from the information about evaluation given in File 1 of raw
	# ENDF-6 file.
	
def FindMT(MTfind, ifile_rawENDF6):
	MFs = numpy.zeros(1000); MTs = numpy.zeros(1000)
	(nfiles,MFs,MTs) = file1(ifile_rawENDF6)

	iflMTpr = 0
	for i in range(nfiles):
		if (MFs[i]==3 and MTs[i]==MTfind):
			iflMTpr = 1
			break
	return(iflMTpr)

#=======Get MFs and MTs from File 1 of raw data=======*
 	
## The directory of Files (MFs) and the corresponding Sections (MTs)
## given in the evaluation are read from File 1 and returned.
	
def file1(ifile_rawENDF6): 
	## maximum of NXC = 350 (ENDF-102), 
	## but deviates for Mn55 ENDF/B-VII.1, so changed to 1000 
	ifile_rawENDF6.seek(0, 0)
	ifile = ifile_rawENDF6
	ifile.readline()
	line = ifile.readline()
	data = eachlineinfo(line)
	ZA = float(data[0]); AWR = float(data[1])
	LRP = int(data[2]); LFI = int(data[3]); NLIB = int(data[4])
	NMOD = int(data[5]); MAT = int(data[6]); MF = int(data[7])
	MT = int(data[8])

	line = ifile.readline()
	data = eachlineinfo(line)
	ELIS = float(data[0]); STA = float(data[1]);
	LIS = int(data[2]); LISO = int(data[3]); num = int(data[4]);
	NFOR = int(data[5]); MAT = int(data[6]); MF = int(data[7])
	MT = int(data[8])

	line = ifile.readline()
	data = eachlineinfo(line)
	AWI = float(data[0]); EMAX = float(data[1])
	LREL = int(data[2]); num = int(data[3]); NSUB = int(data[4]);
	NVER = int(data[5]); MAT = int(data[6]); MF = int(data[7])
	MT = int(data[8])

	line = ifile.readline()
	data = eachlineinfo(line)
	TEMP = float(data[0]); c2 = float(data[1])
	LDRV = int(data[2]); num = int(data[3]); NWD = int(data[4])
	NXC = int(data[5]); MAT = int(data[6]); MF = int(data[7])
	MT = int(data[8])
	# NWD = number of descriptive records (lines)
	# NXC = number of (MF,MT,NC,MOD) lines
	for i in range(NWD):
		ifile.readline()

	MFs = numpy.zeros(1000); MTs = numpy.zeros(1000)

	for i in range(NXC):
		line = ifile.readline()
		data = eachlineinfo(line)
		blnk = data[0]; blnk = data[1]
		MFs[i] = int(data[2]); MTs[i] = int(data[3])
		NCn = int(data[4]); MODn = int(data[5]); MAT = int(data[6])
		MF = int(data[7]); MT = int(data[8])

	nfiles = NXC

	return(nfiles,MFs,MTs)
## =====================================================================================

## Calculation of total (n,p), (n,d), (n,tr), (n,3He) and (n,a)
## cross sections from MT = 600 to 849.
		
def adddiscnth(ifile_preprocessedENDF6,iflMTpr,MTnth,Etu,NPt):
	signth = numpy.zeros(NPt)
	iflMTpr = 0
	if (MTnth == 103):
		MTi = 600
		MTimax = 649
	if (MTnth == 104):
		MTi = 650
		MTimax = 699
	if (MTnth == 105):
		MTi = 700
		MTimax = 749
	if (MTnth == 106):
		MTi = 750
		MTimax = 799
	if (MTnth == 107):
		MTi = 800
		MTimax = 849

	ifile_preprocessedENDF6.seek(0, 0)

	while (True):
		NP1 = 0
		iflp = 0
		if (MTi <= MTimax):
			line = ifile_preprocessedENDF6.readline()
			if (line == ''):
				break
			data = eachlineinfo(line)
			MAT = int(data[6]); MF = int(data[7]); MT  = int(data[8])
			if (MT == 0 and MF != 0):
				line = ifile_preprocessedENDF6.readline()
				(ZAv,AWRv,L0,L1,NKv,L2,MAT,MF,MT) = line_type1_info(line)
			if (MAT != -1):
				if (MF == 3):
					if (MT == MTi or MT == MTimax):
						if (MT == MTimax):
							MTi = MT
						iflp = 1
						iflMTpr = 1
						write(*,*) 'MT = ', MTi, '.....'
						ZA = ZAv
						AWR = AWRv
						read(500,2) QM, QI, LR, NR, NP1
						allocate(E1(NP1),sig1(NP1))
						read(500,*) 
						read(500,3) (E1(j), sig1(j), j = 1, NP1)
						MTi = MTi + 1
					end if
				end if
			else:
				break
	
			if (iflp == 1):
				for i in range(NPt):
					if (Etu[i] >= E1[0]):
						istart = i
						break
				for i in range(istart-1, NPt):
					for j in range(NP1):
						if (Etu[i] == E1[j]):
							signth[i] = signth[i] + sig1[j]
							break
						if (Etu[i] > E1[j] and Etu[i] < E1[j+1]):
							zf = crstd(Etu[i],E1[j],E1[j+1],sig1[j],sig1[j+1])
							signth[i] = signth[i] + zf
							break
			else:
				break

	return(signth)

##==================================================================

## The weighting spectrum for multigrouping of cross sections
		
def spectrum (En,L):
 	match L:
	case 1:
		fi = En * math.exp(-En/0.0253)
		return(fi)
	case 2:
		fi = 1/En
		return(fi)
	case 3:
		fi = math.sqrt(En) * math.exp(-En/(1.415e+6))
		return(fi)
	case 4:
		fi = 1
		return(fi)
	case _:
		fi = 0
		return(fi)

def crstd(x,x1,x2,y1,y2):
	y = y1 + ((y2-y1)*(x-x1)/(x2-x1))
	return(y)

##=======================================================================

## Energy Group Structures

def engrp1():
	Ngl = 176
	Eg = [1.000E-05,1.000E-01,4.140E-01,5.320E-01,6.830E-01,
	8.760E-01,1.130E+00,1.450E+00,1.860E+00,2.380E+00,	
	3.060E+00,3.930E+00,5.040E+00,6.480E+00,8.320E+00,	
	1.070E+01,1.370E+01,1.760E+01,2.260E+01,2.900E+01,	
	3.730E+01,4.790E+01,6.140E+01,7.890E+01,1.010E+02,	
	1.300E+02,1.670E+02,2.140E+02,2.750E+02,3.540E+02,	
	4.540E+02,5.830E+02,7.490E+02,9.610E+02,1.230E+03,	
	1.580E+03,2.030E+03,2.250E+03,2.490E+03,2.610E+03,	
	2.750E+03,3.040E+03,3.350E+03,3.710E+03,4.310E+03,	
	5.530E+03,7.100E+03,9.120E+03,1.060E+04,1.170E+04,	
	1.500E+04,1.930E+04,2.190E+04,2.360E+04,2.420E+04,	
	2.480E+04,2.610E+04,2.700E+04,2.850E+04,3.180E+04,	
	3.430E+04,4.090E+04,4.630E+04,5.250E+04,5.660E+04,	
	6.740E+04,7.200E+04,7.950E+04,8.250E+04,8.650E+04,	
	9.800E+04,1.110E+05,1.170E+05,1.230E+05,1.290E+05,	
	1.360E+05,1.430E+05,1.500E+05,1.580E+05,1.660E+05,	
	1.740E+05,1.830E+05,1.930E+05,2.020E+05,2.130E+05,	
	2.240E+05,2.350E+05,2.470E+05,2.730E+05,2.870E+05,	
	2.950E+05,2.970E+05,2.990E+05,3.020E+05,3.340E+05,	
	3.690E+05,3.880E+05,4.080E+05,4.500E+05,4.980E+05,	
	5.230E+05,5.500E+05,5.780E+05,6.080E+05,6.390E+05,	
	6.720E+05,7.070E+05,7.430E+05,7.810E+05,8.210E+05,	
	8.630E+05,9.070E+05,9.620E+05,1.000E+06,1.110E+06,	
	1.160E+06,1.220E+06,1.290E+06,1.350E+06,1.420E+06,	
	1.500E+06,1.570E+06,1.650E+06,1.740E+06,1.830E+06,	
	1.920E+06,2.020E+06,2.120E+06,2.230E+06,2.310E+06,	
	2.350E+06,2.370E+06,2.390E+06,2.470E+06,2.590E+06,	
	2.730E+06,2.870E+06,3.010E+06,3.170E+06,3.330E+06,	
	3.680E+06,4.070E+06,4.490E+06,4.720E+06,4.970E+06,	
	5.220E+06,5.490E+06,5.770E+06,6.070E+06,6.380E+06,	
	6.590E+06,6.700E+06,7.050E+06,7.410E+06,7.790E+06,	
	8.190E+06,8.610E+06,9.050E+06,9.510E+06,1.000E+07,	
	1.050E+07,1.110E+07,1.160E+07,1.220E+07,1.250E+07,	
	1.280E+07,1.350E+07,1.380E+07,1.420E+07,1.460E+07,	
	1.490E+07,1.570E+07,1.650E+07,1.690E+07,1.730E+07,
	1.960E+07]	

	return(Eg, Ngl)

def engrp2():
	Ngl = 239
	Eg = [0.1000E-04,0.1000E-03,0.5000E-03,0.7500E-03,0.1000E-02,
	0.1200E-02,0.1500E-02,0.2000E-02,0.2500E-02,0.3000E-02,
	0.4000E-02,0.5000E-02,0.7500E-02,0.1000E-01,0.2530E-01,
	0.3000E-01,0.4000E-01,0.5000E-01,0.6000E-01,0.7000E-01,
	0.8000E-01,0.9000E-01,0.1000E+00,0.1250E+00,0.1500E+00,
	0.1750E+00,0.2000E+00,0.2250E+00,0.2500E+00,0.2750E+00,
	0.3000E+00,0.3250E+00,0.3500E+00,0.3750E+00,0.4000E+00,
	0.4500E+00,0.5000E+00,0.5500E+00,0.6000E+00,0.6250E+00,
	0.6500E+00,0.7000E+00,0.7500E+00,0.8000E+00,0.8500E+00,
	0.9000E+00,0.9250E+00,0.9500E+00,0.9750E+00,0.1000E+01,
	0.1010E+01,0.1020E+01,0.1030E+01,0.1040E+01,0.1050E+01,
	0.1060E+01,0.1070E+01,0.1080E+01,0.1090E+01,0.1100E+01,
	0.1110E+01,0.1120E+01,0.1130E+01,0.1140E+01,0.1150E+01,
	0.1175E+01,0.1200E+01,0.1225E+01,0.1250E+01,0.1300E+01,
	0.1350E+01,0.1400E+01,0.1450E+01,0.1500E+01,0.1590E+01,
	0.1680E+01,0.1770E+01,0.1860E+01,0.1940E+01,0.2000E+01,
	0.2120E+01,0.2210E+01,0.2300E+01,0.2380E+01,0.2470E+01,
	0.2570E+01,0.2670E+01,0.2770E+01,0.2870E+01,0.2970E+01,
	0.3000E+01,0.3050E+01,0.3150E+01,0.3500E+01,0.3730E+01,
	0.4000E+01,0.4750E+01,0.5000E+01,0.5400E+01,0.6000E+01,
	0.6250E+01,0.6500E+01,0.6750E+01,0.7000E+01,0.7150E+01,
	0.8100E+01,0.9100E+01,0.1000E+02,0.1150E+02,0.1190E+02,
	0.1290E+02,0.1375E+02,0.1440E+02,0.1510E+02,0.1600E+02,
	0.1700E+02,0.1850E+02,0.1900E+02,0.2000E+02,0.2100E+02,
	0.2250E+02,0.2500E+02,0.2750E+02,0.3000E+02,0.3125E+02,
	0.3175E+02,0.3325E+02,0.3375E+02,0.3460E+02,0.3550E+02,
	0.3700E+02,0.3800E+02,0.3910E+02,0.3960E+02,0.4100E+02,
	0.4240E+02,0.4400E+02,0.4520E+02,0.4700E+02,0.4830E+02,
	0.4920E+02,0.5060E+02,0.5200E+02,0.5340E+02,0.5900E+02,
	0.6100E+02,0.6500E+02,0.6750E+02,0.7200E+02,0.7600E+02,
	0.8000E+02,0.8200E+02,0.9000E+02,0.1000E+03,0.1080E+03,
	0.1150E+03,0.1190E+03,0.1220E+03,0.1860E+03,0.1925E+03,
	0.2075E+03,0.2100E+03,0.2400E+03,0.2850E+03,0.3050E+03,
	0.5500E+03,0.6700E+03,0.6830E+03,0.9500E+03,0.1150E+04,
	0.1500E+04,0.1550E+04,0.1800E+04,0.2200E+04,0.2290E+04,
	0.2580E+04,0.3000E+04,0.3740E+04,0.3900E+04,0.6000E+04,
	0.8030E+04,0.9500E+04,0.1300E+05,0.1700E+05,0.2500E+05,
	0.3000E+05,0.4500E+05,0.5000E+05,0.5200E+05,0.6000E+05,
	0.7300E+05,0.7500E+05,0.8200E+05,0.8500E+05,0.1000E+06,
	0.1283E+06,0.1500E+06,0.2000E+06,0.2700E+06,0.3300E+06,
	0.4000E+06,0.4200E+06,0.4400E+06,0.4700E+06,0.4995E+06,
	0.5500E+06,0.5730E+06,0.6000E+06,0.6700E+06,0.6790E+06,
	0.7500E+06,0.8200E+06,0.8611E+06,0.8750E+06,0.9000E+06,
	0.9200E+06,0.1010E+07,0.1100E+07,0.1200E+07,0.1250E+07,
	0.1317E+07,0.1356E+07,0.1400E+07,0.1500E+07,0.1850E+07,
	0.2354E+07,0.2479E+07,0.3000E+07,0.4304E+07,0.4800E+07,
	0.6434E+07,0.8187E+07,0.1000E+08,0.1284E+08,0.1384E+08,
	0.1455E+08,0.1568E+08,0.1733E+08,0.2000E+08] 	
	
	return(Eg, Ngl)

def engrp3():
	Ngl = 199
	Eg = [1.0000E-05,5.0000E-04,2.0000E-03,5.0000E-03,1.0000E-02,
	1.4500E-02,2.1000E-02,3.0000E-02,4.0000E-02,5.0000E-02,
	7.0000E-02,1.0000E-01,1.2500E-01,1.5000E-01,1.8400E-01,
	2.2500E-01,2.7500E-01,3.2500E-01,3.6680E-01,4.1399E-01,
	5.0000E-01,5.3158E-01,6.2506E-01,6.8256E-01,8.0000E-01,
	8.7643E-01,1.0000E+00,1.0400E+00,1.0800E+00,1.1253E+00,
	1.3000E+00,1.4450E+00,1.8554E+00,2.3824E+00,3.0590E+00,
	3.9279E+00,5.0435E+00,6.4760E+00,8.3153E+00,1.0677E+01,
	1.3710E+01,1.7604E+01,2.2603E+01,2.9023E+01,3.7266E+01,
	4.7851E+01,6.1442E+01,7.8893E+01,1.0130E+02,1.3007E+02,
	1.6702E+02,2.1445E+02,2.7536E+02,3.5357E+02,4.5400E+02,
	5.8295E+02,7.4852E+02,9.6112E+02,1.2341E+03,1.5846E+03,
	2.0347E+03,2.2487E+03,2.4852E+03,2.6126E+03,2.7465E+03,
	3.0354E+03,3.3546E+03,3.7074E+03,4.3074E+03,5.5308E+03,
	7.1017E+03,9.1188E+03,1.0595E+04,1.1709E+04,1.5034E+04,
	1.9305E+04,2.1875E+04,2.3579E+04,2.4176E+04,2.4788E+04,
	2.6058E+04,2.7000E+04,2.8501E+04,3.1828E+04,3.4307E+04,
	4.0868E+04,4.6309E+04,5.2475E+04,5.6562E+04,6.7379E+04,
	7.1998E+04,7.9499E+04,8.2503E+04,8.6517E+04,9.8037E+04,
	1.1109E+05,1.1679E+05,1.2277E+05,1.2907E+05,1.3569E+05,
	1.4264E+05,1.4996E+05,1.5764E+05,1.6573E+05,1.7422E+05,
	1.8316E+05,1.9255E+05,2.0242E+05,2.1280E+05,2.2371E+05,
	2.3518E+05,2.4724E+05,2.7324E+05,2.8725E+05,2.9452E+05,
	2.9721E+05,2.9849E+05,3.0197E+05,3.3373E+05,3.6883E+05,
	3.8774E+05,4.0762E+05,4.5049E+05,4.9787E+05,5.2340E+05,
	5.5023E+05,5.7844E+05,6.0810E+05,6.3928E+05,6.7206E+05,
	7.0651E+05,7.4274E+05,7.8082E+05,8.2085E+05,8.6294E+05,
	9.0718E+05,9.6164E+05,1.0026E+06,1.1080E+06,1.1648E+06,
	1.2246E+06,1.2874E+06,1.3534E+06,1.4227E+06,1.4957E+06,
	1.5724E+06,1.6530E+06,1.7377E+06,1.8268E+06,1.9205E+06,
	2.0190E+06,2.1225E+06,2.2313E+06,2.3069E+06,2.3457E+06,
	2.3653E+06,2.3852E+06,2.4660E+06,2.5924E+06,2.7253E+06,
	2.8651E+06,3.0119E+06,3.1664E+06,3.3287E+06,3.6788E+06,
	4.0657E+06,4.4933E+06,4.7237E+06,4.9659E+06,5.2205E+06,
	5.4881E+06,5.7695E+06,6.0653E+06,6.3763E+06,6.5924E+06,
	6.7032E+06,7.0469E+06,7.4082E+06,7.7880E+06,8.1873E+06,
	8.6071E+06,9.0484E+06,9.5123E+06,1.0000E+07,1.0513E+07,
	1.1052E+07,1.1618E+07,1.2214E+07,1.2523E+07,1.2840E+07,
	1.3499E+07,1.3840E+07,1.4191E+07,1.4550E+07,1.4918E+07,
	1.5683E+07,1.6487E+07,1.6905E+07,1.9640E+07]

	return(Eg, Ngl)

def engrp4():
	Ngl = 34
	Eg = [1.000e-05,1.000e-01,5.400e-01,4.000e+00,8.315e+00,
	1.371e+01,2.260e+01,4.017e+01,6.790e+01,9.166e+01,1.486e+02,
	3.043e+02,4.540e+02,7.485e+02,1.230e+03,2.030e+03,3.355e+03,
	5.531e+03,9.119e+03,1.503e+04,2.479e+04,4.087e+04,6.738e+04,
	1.111e+05,1.832e+05,3.020e+05,4.979e+05,8.209e+05,1.353e+06,
	2.231e+06,3.679e+06,6.065e+06,1.000e+07,1.964e+07]

	return(Eg, Ngl)

def engrp5():
	Ngl = 27
	Eg = [1.0E-05,0.0253E+00,0.4642E+00,1.0E+00,2.1544E+00,
	4.6416E+00,1.0E+01,2.1544E+01,4.6416E+01,1.0E+02,2.1544E+02,
	4.6416E+02,1.0E+03,2.1544E+03,4.6416E+03,1.0E+04,2.1544E+04,
	4.6416E+04,1.0E+05,2.0E+05,0.4E+06,0.8E+06,1.4E+06,2.5E+06, 
	4.0E+06,6.5E+06,1.05E+07]

	return(Eg, Ngl)

def engrp6():
	Ngl = 101
	Eg = [1.0000E-05,4.1399E-01,5.3158E-01,6.8256E-01,8.7642E-01,
	1.1254E+00,1.4450E+00,1.8554E+00,2.3824E+00,3.0590E+00,
	3.9279E+00,5.0435E+00,6.4760E+00,8.3153E+00,1.0677E+01,
	1.3710E+01,1.7603E+01,2.2603E+01,2.9023E+01,3.7267E+01,
	4.7851E+01,6.1442E+01,7.8893E+01,1.0130E+02,1.3007E+02,
	1.6702E+02,2.1445E+02,2.7536E+02,3.5358E+02,4.5400E+02,
	5.8295E+02,7.4852E+02,9.6112E+02,1.2341E+03,1.5846E+03,
	2.0347E+03,2.6126E+03,3.3546E+03,4.3074E+03,5.5308E+03,
	7.1017E+03,9.1188E+03,1.1709E+04,1.5034E+04,1.9305E+04,
	2.4788E+04,3.1828E+04,4.0868E+04,5.2475E+04,6.7379E+04,
	8.6517E+04,1.1109E+05,1.2277E+05,1.3569E+05,1.4996E+05,
	1.6573E+05,1.8316E+05,2.0242E+05,2.2371E+05,2.4724E+05,
	2.7324E+05,3.0197E+05,3.3373E+05,3.6883E+05,4.0762E+05,
	4.5049E+05,4.9787E+05,5.5023E+05,6.0810E+05,6.7206E+05,
	7.4274E+05,8.2085E+05,9.0718E+05,1.0026E+06,1.1080E+06,
	1.2246E+06,1.3534E+06,1.4957E+06,1.6530E+06,1.8268E+06,
	2.0190E+06,2.2313E+06,2.4660E+06,2.7253E+06,3.0119E+06,
	3.3287E+06,3.6788E+06,4.0657E+06,4.4933E+06,4.9659E+06,
	5.4881E+06,6.0653E+06,6.7032E+06,7.4082E+06,8.1873E+06,
	9.0484E+06,1.0000E+07,1.1052E+07,1.2214E+07,1.3500E+07,
	1.5000E+07]

	return(Eg, Ngl)

def engrp7():		## EXTENDED FROM 198 GROUP STRUCTURE
	Ngl = 229
	Eg = [1.0000E-05,5.0000E-04,2.0000E-03,5.0000E-03,1.0000E-02,
	1.4500E-02,2.1000E-02,3.0000E-02,4.0000E-02,5.0000E-02,
	7.0000E-02,1.0000E-01,1.2500E-01,1.5000E-01,1.8400E-01,
	2.2500E-01,2.7500E-01,3.2500E-01,3.6680E-01,4.1399E-01,
	5.0000E-01,5.3158E-01,6.2506E-01,6.8256E-01,8.0000E-01,
	8.7643E-01,1.0000E+00,1.0400E+00,1.0800E+00,1.1253E+00,
	1.3000E+00,1.4450E+00,1.8554E+00,2.3824E+00,3.0590E+00,
	3.9279E+00,5.0435E+00,6.4760E+00,8.3153E+00,1.0677E+01,
	1.3710E+01,1.7604E+01,2.2603E+01,2.9023E+01,3.7266E+01,
	4.7851E+01,6.1442E+01,7.8893E+01,1.0130E+02,1.3007E+02,
	1.6702E+02,2.1445E+02,2.7536E+02,3.5357E+02,4.5400E+02,
	5.8295E+02,7.4852E+02,9.6112E+02,1.2341E+03,1.5846E+03,
	2.0347E+03,2.2487E+03,2.4852E+03,2.6126E+03,2.7465E+03,
	3.0354E+03,3.3546E+03,3.7074E+03,4.3074E+03,5.5308E+03,
	7.1017E+03,9.1188E+03,1.0595E+04,1.1709E+04,1.5034E+04,
	1.9305E+04,2.1875E+04,2.3579E+04,2.4176E+04,2.4788E+04,
	2.6058E+04,2.7000E+04,2.8501E+04,3.1828E+04,3.4307E+04,
	4.0868E+04,4.6309E+04,5.2475E+04,5.6562E+04,6.7379E+04,
	7.1998E+04,7.9499E+04,8.2503E+04,8.6517E+04,9.8037E+04,
	1.1109E+05,1.1679E+05,1.2277E+05,1.2907E+05,1.3569E+05,
	1.4264E+05,1.4996E+05,1.5764E+05,1.6573E+05,1.7422E+05,
	1.8316E+05,1.9255E+05,2.0242E+05,2.1280E+05,2.2371E+05,
	2.3518E+05,2.4724E+05,2.7324E+05,2.8725E+05,2.9452E+05,
	2.9721E+05,2.9849E+05,3.0197E+05,3.3373E+05,3.6883E+05,
	3.8774E+05,4.0762E+05,4.5049E+05,4.9787E+05,5.2340E+05,
	5.5023E+05,5.7844E+05,6.0810E+05,6.3928E+05,6.7206E+05,
	7.0651E+05,7.4274E+05,7.8082E+05,8.2085E+05,8.6294E+05,
	9.0718E+05,9.6164E+05,1.0026E+06,1.1080E+06,1.1648E+06,
	1.2246E+06,1.2874E+06,1.3534E+06,1.4227E+06,1.4957E+06,
	1.5724E+06,1.6530E+06,1.7377E+06,1.8268E+06,1.9205E+06,
	2.0190E+06,2.1225E+06,2.2313E+06,2.3069E+06,2.3457E+06,
	2.3653E+06,2.3852E+06,2.4660E+06,2.5924E+06,2.7253E+06,
	2.8651E+06,3.0119E+06,3.1664E+06,3.3287E+06,3.6788E+06,
	4.0657E+06,4.4933E+06,4.7237E+06,4.9659E+06,5.2205E+06,
	5.4881E+06,5.7695E+06,6.0653E+06,6.3763E+06,6.5924E+06,
	6.7032E+06,7.0469E+06,7.4082E+06,7.7880E+06,8.1873E+06,
	8.6071E+06,9.0484E+06,9.5123E+06,1.0000E+07,1.0513E+07,
	1.1052E+07,1.1618E+07,1.2214E+07,1.2523E+07,1.2840E+07,
	1.3499E+07,1.3840E+07,1.4191E+07,1.4550E+07,1.4918E+07,
	1.5683E+07,1.6487E+07,1.6905E+07,1.9640E+07,
	2.10E+07,2.40E+07,2.60E+07,2.80E+07,3.00E+07,3.20E+07,
	3.50E+07,3.80E+07,4.00E+07,4.50E+07,5.00E+07,5.50E+07,
	6.00E+07,7.00E+07,8.00E+07,9.00E+07,1.00E+08,1.10E+08,
	1.20E+08,1.30E+08,1.40E+08,1.50E+08,1.60E+08,1.70E+08,
	1.80E+08,1.90E+08,1.95E+08,1.98E+08,1.99E+08,2.00E+08]

	return(Eg, Ngl)

def engrp8():	# EXTENDED FROM 198 GROUP STRUCTURE
	Ngl = 229
	Eg = [1.0000E-05,5.0000E-04,2.0000E-03,5.0000E-03,1.0000E-02,
	1.4500E-02,2.1000E-02,3.0000E-02,4.0000E-02,5.0000E-02,
	7.0000E-02,1.0000E-01,1.2500E-01,1.5000E-01,1.8400E-01,
	2.2500E-01,2.7500E-01,3.2500E-01,3.6680E-01,4.1399E-01,
	5.0000E-01,5.3158E-01,6.2506E-01,6.8256E-01,8.0000E-01,
	8.7643E-01,1.0000E+00,1.0400E+00,1.0800E+00,1.1253E+00,
	1.3000E+00,1.4450E+00,1.8554E+00,2.3824E+00,3.0590E+00,
	3.9279E+00,5.0435E+00,6.4760E+00,8.3153E+00,1.0677E+01,
	1.3710E+01,1.7604E+01,2.2603E+01,2.9023E+01,3.7266E+01,
	4.7851E+01,6.1442E+01,7.8893E+01,1.0130E+02,1.3007E+02,
	1.6702E+02,2.1445E+02,2.7536E+02,3.5357E+02,4.5400E+02,
	5.8295E+02,7.4852E+02,9.6112E+02,1.2341E+03,1.5846E+03,
	2.0347E+03,2.2487E+03,2.4852E+03,2.6126E+03,2.7465E+03,
	3.0354E+03,3.3546E+03,3.7074E+03,4.3074E+03,5.5308E+03,
	7.1017E+03,9.1188E+03,1.0595E+04,1.1709E+04,1.5034E+04,
	1.9305E+04,2.1875E+04,2.3579E+04,2.4176E+04,2.4788E+04,
	2.6058E+04,2.7000E+04,2.8501E+04,3.1828E+04,3.4307E+04,
	4.0868E+04,4.6309E+04,5.2475E+04,5.6562E+04,6.7379E+04,
	7.1998E+04,7.9499E+04,8.2503E+04,8.6517E+04,9.8037E+04,
	1.1109E+05,1.1679E+05,1.2277E+05,1.2907E+05,1.3569E+05,
	1.4264E+05,1.4996E+05,1.5764E+05,1.6573E+05,1.7422E+05,
	1.8316E+05,1.9255E+05,2.0242E+05,2.1280E+05,2.2371E+05,
	2.3518E+05,2.4724E+05,2.7324E+05,2.8725E+05,2.9452E+05,
	2.9721E+05,2.9849E+05,3.0197E+05,3.3373E+05,3.6883E+05,
	3.8774E+05,4.0762E+05,4.5049E+05,4.9787E+05,5.2340E+05,
	5.5023E+05,5.7844E+05,6.0810E+05,6.3928E+05,6.7206E+05,
	7.0651E+05,7.4274E+05,7.8082E+05,8.2085E+05,8.6294E+05,
	9.0718E+05,9.6164E+05,1.0026E+06,1.1080E+06,1.1648E+06,
	1.2246E+06,1.2874E+06,1.3534E+06,1.4227E+06,1.4957E+06,
	1.5724E+06,1.6530E+06,1.7377E+06,1.8268E+06,1.9205E+06,
	2.0190E+06,2.1225E+06,2.2313E+06,2.3069E+06,2.3457E+06,
	2.3653E+06,2.3852E+06,2.4660E+06,2.5924E+06,2.7253E+06,
	2.8651E+06,3.0119E+06,3.1664E+06,3.3287E+06,3.6788E+06,
	4.0657E+06,4.4933E+06,4.7237E+06,4.9659E+06,5.2205E+06,
	5.4881E+06,5.7695E+06,6.0653E+06,6.3763E+06,6.5924E+06,
	6.7032E+06,7.0469E+06,7.4082E+06,7.7880E+06,8.1873E+06,
	8.6071E+06,9.0484E+06,9.5123E+06,1.0000E+07,1.0513E+07,
	1.1052E+07,1.1618E+07,1.2214E+07,1.2523E+07,1.2840E+07,
	1.3499E+07,1.3840E+07,1.4191E+07,1.4550E+07,1.4918E+07,
	1.5683E+07,1.6487E+07,1.6905E+07,1.9640E+07,
	2.10E+07,2.40E+07,2.60E+07,2.80E+07,3.00E+07,3.20E+07,
	3.50E+07,3.80E+07,4.00E+07,4.50E+07,5.00E+07,5.50E+07,
	6.00E+07,6.50E+07,7.00E+07,7.50E+07,8.00E+07,8.50E+07,
	9.00E+07,9.50E+07,1.00E+08,1.05E+08,1.10E+08,1.20E+08,
	1.25E+08,1.30E+08,1.35E+08,1.40E+08,1.45E+08,1.50E+08]

	return(Eg, Ngl)

def engrp9():	# 616 energy groups DEMO-HCPB-FW
	Ngl = 617
	Eg = [1.00E-05,1.05E-05,1.10E-05,1.15E-05,1.20E-05,1.26E-05,
	1.32E-05,1.38E-05,1.45E-05,1.51E-05,1.58E-05,1.66E-05,1.74E-05,
	1.82E-05,1.91E-05,2.00E-05,2.09E-05,2.19E-05,2.29E-05,2.40E-05,
	2.51E-05,2.63E-05,2.75E-05,2.88E-05,3.02E-05,3.16E-05,3.31E-05,
	3.47E-05,3.63E-05,3.80E-05,3.98E-05,4.17E-05,4.37E-05,4.57E-05,
	4.79E-05,5.01E-05,5.25E-05,5.50E-05,5.75E-05,6.03E-05,6.31E-05,
	6.61E-05,6.92E-05,7.24E-05,7.59E-05,7.94E-05,8.32E-05,8.71E-05,
	9.12E-05,9.55E-05,1.00E-04,1.05E-04,1.10E-04,1.15E-04,1.20E-04,
	1.26E-04,1.32E-04,1.38E-04,1.45E-04,1.51E-04,1.58E-04,1.66E-04,
	1.74E-04,1.82E-04,1.91E-04,2.00E-04,2.09E-04,2.19E-04,2.29E-04,
	2.40E-04,2.51E-04,2.63E-04,2.75E-04,2.88E-04,3.02E-04,3.16E-04,
	3.31E-04,3.47E-04,3.63E-04,3.80E-04,3.98E-04,4.17E-04,4.37E-04,
	4.57E-04,4.79E-04,5.01E-04,5.25E-04,5.50E-04,5.75E-04,6.03E-04,
	6.31E-04,6.61E-04,6.92E-04,7.24E-04,7.59E-04,7.94E-04,8.32E-04,
	8.71E-04,9.12E-04,9.55E-04,1.00E-03,1.05E-03,1.10E-03,1.15E-03,
	1.20E-03,1.26E-03,1.32E-03,1.38E-03,1.45E-03,1.51E-03,1.58E-03,
	1.66E-03,1.74E-03,1.82E-03,1.91E-03,2.00E-03,2.09E-03,2.19E-03,
	2.29E-03,2.40E-03,2.51E-03,2.63E-03,2.75E-03,2.88E-03,3.02E-03,
	3.16E-03,3.31E-03,3.47E-03,3.63E-03,3.80E-03,3.98E-03,4.17E-03,
	4.37E-03,4.57E-03,4.79E-03,5.01E-03,5.25E-03,5.50E-03,5.75E-03,
	6.03E-03,6.31E-03,6.61E-03,6.92E-03,7.24E-03,7.59E-03,7.94E-03,
	8.32E-03,8.71E-03,9.12E-03,9.55E-03,1.00E-02,1.05E-02,1.10E-02,
	1.15E-02,1.20E-02,1.26E-02,1.32E-02,1.38E-02,1.45E-02,1.51E-02,
	1.58E-02,1.66E-02,1.74E-02,1.82E-02,1.91E-02,2.00E-02,2.09E-02,
	2.19E-02,2.29E-02,2.40E-02,2.51E-02,2.63E-02,2.75E-02,2.88E-02,
	3.02E-02,3.16E-02,3.31E-02,3.47E-02,3.63E-02,3.80E-02,3.98E-02,
	4.17E-02,4.37E-02,4.57E-02,4.79E-02,5.01E-02,5.25E-02,5.50E-02,
	5.75E-02,6.03E-02,6.31E-02,6.61E-02,6.92E-02,7.24E-02,7.59E-02,
	7.94E-02,8.32E-02,8.71E-02,9.12E-02,9.55E-02,1.00E-01,1.05E-01,
	1.10E-01,1.15E-01,1.20E-01,1.26E-01,1.32E-01,1.38E-01,1.45E-01,
	1.51E-01,1.58E-01,1.66E-01,1.74E-01,1.82E-01,1.91E-01,2.00E-01,
	2.09E-01,2.19E-01,2.29E-01,2.40E-01,2.51E-01,2.63E-01,2.75E-01,
	2.88E-01,3.02E-01,3.16E-01,3.31E-01,3.47E-01,3.63E-01,3.80E-01,
	3.98E-01,4.17E-01,4.37E-01,4.57E-01,4.79E-01,5.01E-01,5.25E-01,
	5.50E-01,5.75E-01,6.03E-01,6.31E-01,6.61E-01,6.92E-01,7.24E-01,
	7.59E-01,7.94E-01,8.32E-01,8.71E-01,9.12E-01,9.55E-01,1.00E+00,
	1.05E+00,1.10E+00,1.15E+00,1.20E+00,1.26E+00,1.32E+00,1.38E+00,
	1.45E+00,1.51E+00,1.58E+00,1.66E+00,1.74E+00,1.82E+00,1.91E+00,
	2.00E+00,2.09E+00,2.19E+00,2.29E+00,2.40E+00,2.51E+00,2.63E+00,
	2.75E+00,2.88E+00,3.02E+00,3.16E+00,3.31E+00,3.47E+00,3.63E+00,
	3.80E+00,3.98E+00,4.17E+00,4.37E+00,4.57E+00,4.79E+00,5.01E+00,
	5.25E+00,5.50E+00,5.75E+00,6.03E+00,6.31E+00,6.61E+00,6.92E+00,
	7.24E+00,7.59E+00,7.94E+00,8.32E+00,8.71E+00,9.12E+00,9.55E+00,
	1.00E+01,1.05E+01,1.10E+01,1.15E+01,1.20E+01,1.26E+01,1.32E+01,
	1.38E+01,1.45E+01,1.51E+01,1.58E+01,1.66E+01,1.74E+01,1.82E+01,
	1.91E+01,2.00E+01,2.09E+01,2.19E+01,2.29E+01,2.40E+01,2.51E+01,
	2.63E+01,2.75E+01,2.88E+01,3.02E+01,3.16E+01,3.31E+01,3.47E+01,
	3.63E+01,3.80E+01,3.98E+01,4.17E+01,4.37E+01,4.57E+01,4.79E+01,
	5.01E+01,5.25E+01,5.50E+01,5.75E+01,6.03E+01,6.31E+01,6.61E+01,
	6.92E+01,7.24E+01,7.59E+01,7.94E+01,8.32E+01,8.71E+01,9.12E+01,
	9.55E+01,1.00E+02,1.05E+02,1.10E+02,1.15E+02,1.20E+02,1.26E+02,
	1.32E+02,1.38E+02,1.45E+02,1.51E+02,1.58E+02,1.66E+02,1.74E+02,
	1.82E+02,1.91E+02,2.00E+02,2.09E+02,2.19E+02,2.29E+02,2.40E+02,
	2.51E+02,2.63E+02,2.75E+02,2.88E+02,3.02E+02,3.16E+02,3.31E+02,
	3.47E+02,3.63E+02,3.80E+02,3.98E+02,4.17E+02,4.37E+02,4.57E+02,
	4.79E+02,5.01E+02,5.25E+02,5.50E+02,5.75E+02,6.03E+02,6.31E+02,
	6.61E+02,6.92E+02,7.24E+02,7.59E+02,7.94E+02,8.32E+02,8.71E+02,
	9.12E+02,9.55E+02,1.00E+03,1.05E+03,1.10E+03,1.15E+03,1.20E+03,
	1.26E+03,1.32E+03,1.38E+03,1.45E+03,1.51E+03,1.58E+03,1.66E+03,
	1.74E+03,1.82E+03,1.91E+03,2.00E+03,2.09E+03,2.19E+03,2.29E+03,
	2.40E+03,2.51E+03,2.63E+03,2.75E+03,2.88E+03,3.02E+03,3.16E+03,
	3.31E+03,3.47E+03,3.63E+03,3.80E+03,3.98E+03,4.17E+03,4.37E+03,
	4.57E+03,4.79E+03,5.01E+03,5.25E+03,5.50E+03,5.75E+03,6.03E+03,
	6.31E+03,6.61E+03,6.92E+03,7.24E+03,7.59E+03,7.94E+03,8.32E+03,
	8.71E+03,9.12E+03,9.55E+03,1.00E+04,1.05E+04,1.10E+04,1.15E+04,
	1.20E+04,1.26E+04,1.32E+04,1.38E+04,1.45E+04,1.51E+04,1.58E+04,
	1.66E+04,1.74E+04,1.82E+04,1.91E+04,2.00E+04,2.09E+04,2.19E+04,
	2.29E+04,2.40E+04,2.51E+04,2.63E+04,2.75E+04,2.88E+04,3.02E+04,
	3.16E+04,3.31E+04,3.47E+04,3.63E+04,3.80E+04,3.98E+04,4.17E+04,
	4.37E+04,4.57E+04,4.79E+04,5.01E+04,5.25E+04,5.50E+04,5.75E+04,
	6.03E+04,6.31E+04,6.61E+04,6.92E+04,7.24E+04,7.59E+04,7.94E+04,
	8.32E+04,8.71E+04,9.12E+04,9.55E+04,1.00E+05,1.05E+05,1.10E+05,
	1.15E+05,1.20E+05,1.26E+05,1.32E+05,1.38E+05,1.45E+05,1.51E+05,
	1.58E+05,1.66E+05,1.74E+05,1.82E+05,1.91E+05,2.00E+05,2.09E+05,
	2.19E+05,2.29E+05,2.40E+05,2.51E+05,2.63E+05,2.75E+05,2.88E+05,
	3.02E+05,3.16E+05,3.31E+05,3.47E+05,3.63E+05,3.80E+05,3.98E+05,
	4.17E+05,4.37E+05,4.57E+05,4.79E+05,5.01E+05,5.25E+05,5.50E+05,
	5.75E+05,6.03E+05,6.31E+05,6.61E+05,6.92E+05,7.24E+05,7.59E+05,
	7.94E+05,8.32E+05,8.71E+05,9.12E+05,9.55E+05,1.00E+06,1.05E+06,
	1.10E+06,1.15E+06,1.20E+06,1.26E+06,1.32E+06,1.38E+06,1.45E+06,
	1.51E+06,1.58E+06,1.66E+06,1.74E+06,1.82E+06,1.91E+06,2.00E+06,
	2.09E+06,2.19E+06,2.29E+06,2.40E+06,2.51E+06,2.63E+06,2.75E+06,
	2.88E+06,3.02E+06,3.16E+06,3.31E+06,3.47E+06,3.63E+06,3.80E+06,
	3.98E+06,4.17E+06,4.37E+06,4.57E+06,4.79E+06,5.01E+06,5.25E+06,
	5.50E+06,5.75E+06,6.03E+06,6.31E+06,6.61E+06,6.92E+06,7.24E+06,
	7.59E+06,7.94E+06,8.32E+06,8.71E+06,9.12E+06,9.55E+06,1.00E+07,
	1.05E+07,1.10E+07,1.15E+07,1.20E+07,1.26E+07,1.32E+07,1.38E+07,
	1.45E+07,1.51E+07,1.58E+07,1.66E+07,1.74E+07,1.82E+07,1.91E+07,
	2.00E+07,2.00E+07]

	return(Eg, Ngl)