#	>> Code: RMINDD.py, module file: RecedU.py
#	>> Perform: Recoil energy distributions of nuclei from basic ENDF-6 files
#	>> Author: Dr. Uttiyoarnab Saha
#	>> Version and Date: 1.0 and 25/03/2021
#	>> Last modified: 25/03/2021, Kolkata
#	>> Update: 25/03/2021
#	>> Major changes:
#
#======================================================================================================

import numpy
import math
import os

#=======Linear interpolation=======*		
        
def crstd(x,x1,x2,y1,y2):
	y = y1 + ((y2-y1)*(x-x1)/(x2-x1))
	return(y)

#=========================================

# Schemes of interpolation

def TERPOLIN (intflg,x,x1,x2,y1,y2):
	intsch = intflg
	if (intsch == 1 or intsch == 11 or intsch == 21):
		y = y1
	if (intsch == 2 or intsch == 12 or intsch == 22):
		y = y1+(y2-y1)*(x-x1)/(x2-x1)
	if (intsch == 3 or intsch == 13 or intsch == 23):
		y = y1+(y2-y1)*(log(x/x1))/log(x2/x1)
	if (intsch == 4 or intsch == 14 or intsch == 24):
		y = y1*math.exp((x-x1)*log(y2/y1)/(x2-x1))
	if (intsch == 5 or intsch == 15 or intsch == 25):
		y = y1*math.exp(log(x/x1)*log(y2/y1)/log(x2/x1))

	return(y)

# ===============================================================

# To interpolate into required energy points using the specific
# scheme in the range.

def TERPOL (NRin,Intrin,n1,E1,Y1,n2,E2):
	Y2 = [0]*n2
	for i in range (n2):
		for j in range (n1):
			if (E2[i] == E1[j]):
				Y2[i] = Y1[j]
				break
			if (E1[j] < E2[i] and E2[i] < E1[j+1]):
				diff1 = E2[i] - E1[j]
				diff2 = E2[i] - E2[i-1]
				for k in range (len(NRin)):		#20
					if (j <= NRin[k]):
						intflg = Intrin[k]
						break
				if (diff1 <= diff2):
					Y2[i] = TERPOLIN(intflg,E2[i],E1[j],E1[j+1],Y1[j],Y1[j+1])
				if (diff2 < diff1):
					Y2[i] = TERPOLIN(intflg,E2[i],E2[i-1],E1[j+1],Y2[i-1],Y1[j+1])
				break
	return(Y2)

# ===============================================================
	
# Interpolation of Arrays via Particular Rule 

def TERPOLAPR (iprule,n1,x1,y1,n2,x2):
	y2 = [0]*n2
	for i in range (n2):
		for j in range (n1):
			if (x2[i] == x1[j]):
				y2[i] = y1[j]
				break
			if (x1[j] < x2[i] and x2[i] < x1[j+1]):
				diff1 = x2[i] - x1[j]		
				diff2 = x2[i] - x2[i-1]
				if (diff1 <= diff2):
					y2[i] = TERPOLIN(iprule,x2[i],x1[j],x1[j+1],y1[j],y1[j+1])
				if (diff2 < diff1):
					y2[i] = TERPOLIN(iprule,x2[i],x2[i-1],x1[j+1],y2[i-1],y1[j+1])
				break
	return(y2)

#===============================================================

#========The data in each line are explicitly extracted=======*

def eachlineinfo(line):	
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

#========The data in lines of different types of ENDF-6 file are explicitly extracted=======*

def line_type1_info(line):
	data = eachlineinfo(line)
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

def line_type2_info(line):
	data = eachlineinfo(line)
	dataV1 = float(data[0]); dataV2 = float(data[1]); dataV3 = int(data[2])
	dataV4 = int(data[3]); dataV5 = int(data[4]); dataV6 = int(data[5])
	dataV7 = int(data[6]); dataV8 = int(data[7]); dataV9 = int(data[8])
	return(dataV1,dataV2,dataV3,dataV4,dataV5,dataV6,dataV7,dataV8,dataV9)

def line_type3_info (filehandle,numdata,numvariables):
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
			data = eachlineinfo(line)
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
			data = eachlineinfo(line)
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

def line_type4_info(filehandle,nre,numdata):
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

# =======Interpolate cross section to unique common energy=======*
	
	# When the basic / dpa / heating cross sections from partial reactions
	# have to be computed corresponding to the unique energy points from
	# from their individual energy arrays, linear interpolation is performed
	# to find the cross sections corresponding to unique energy points.

def trptuqce (E,s1,Etu):
	s2 = [0]*len(Etu)
	s2 = numpy.interp(Etu, E, s1)
	return(s2)

#=======Gauss-Legendre Quadrature Points and Weights=======*
	
	# GL quadrature 64 points and weights used often for integrature and sometimes
	# as the base angles for which the angular distributions are computed from the
	# given data by interpolation
	
def GQ():
	# ARGUEMENT POINTS FOR THE GAUSS QUADRATURE	
	xabc = numpy.array([-9.99305042E-01,-9.96340117E-01,-9.91013371E-01,
	-9.83336254E-01,-9.73326828E-01,-9.61008800E-01,-9.46411375E-01,-9.29569172E-01,-9.10522137E-01,
	-8.89315446E-01,-8.65999398E-01,-8.40629296E-01,-8.13265315E-01,-7.83972359E-01,-7.52819907E-01,
	-7.19881850E-01,-6.85236313E-01,-6.48965471E-01,-6.11155355E-01,-5.71895646E-01,-5.31279464E-01,
	-4.89403146E-01,-4.46366017E-01,-4.02270158E-01,-3.57220158E-01,-3.11322872E-01,-2.64687162E-01,
	-2.17423644E-01,-1.69644420E-01,-1.21462819E-01,-7.29931218E-02,-2.43502927E-02,2.43502927E-02, 
	7.29931218E-02,1.21462819E-01,1.69644420E-01,2.17423644E-01,2.64687162E-01,3.11322872E-01,  
	3.57220158E-01,4.02270158E-01,4.46366017E-01,4.89403146E-01,5.31279464E-01,5.71895646E-01,  
	6.11155355E-01,6.48965471E-01,6.85236313E-01,7.19881850E-01,7.52819907E-01,7.83972359E-01,  
	8.13265315E-01,8.40629296E-01,8.65999398E-01,8.89315446E-01,9.10522137E-01,9.29569172E-01,  
	9.46411375E-01,9.61008800E-01,9.73326828E-01,9.83336254E-01,9.91013371E-01,9.96340117E-01,9.99305042E-01])

	# WEIGHTS CORRESPONDING TO THE ARGUEMENT POINTS		
	wg = numpy.array([1.78328072E-03,4.14703326E-03,6.50445797E-03,8.84675983E-03,1.11681395E-02,1.34630479E-02,
	1.57260305E-02,1.79517158E-02,2.01348232E-02,2.22701738E-02,2.43527026E-02,2.63774697E-02,2.83396726E-02,
	3.02346571E-02,3.20579284E-02,3.38051618E-02,3.54722133E-02,3.70551285E-02,3.85501532E-02,3.99537411E-02,
	4.12625632E-02,4.24735151E-02,4.35837245E-02,4.45905582E-02,4.54916279E-02,4.62847966E-02,4.69681828E-02,
	4.75401657E-02,4.79993886E-02,4.83447622E-02,4.85754674E-02,4.86909570E-02,4.86909570E-02,4.85754674E-02,
	4.83447622E-02,4.79993886E-02,4.75401657E-02,4.69681828E-02,4.62847966E-02,4.54916279E-02,4.45905582E-02,
	4.35837245E-02,4.24735151E-02,4.12625632E-02,3.99537411E-02,3.85501532E-02,3.70551285E-02,3.54722133E-02,
	3.38051618E-02,3.20579284E-02,3.02346571E-02,2.83396726E-02,2.63774697E-02,2.43527026E-02,2.22701738E-02,
	2.01348232E-02,1.79517158E-02,1.57260305E-02,1.34630479E-02,1.11681395E-02,8.84675983E-03,6.50445797E-03,
	4.14703326E-03,1.78328072E-03])
	
	return(xabc, wg)


#=======Add all pka spectra=======*
	
# Finds the sum of the PKA spectra from partial reactions.
	
def ALLSUM (nrg,nrcta):

	nre = nrg - 1
	numdata = nre*nre
	dsdt = numpy.zeros((nre,nre)); pkaind = numpy.zeros((nre,nre))
	ifile = open('ToAddAll.txt', 'r')

	# maximum number of pka spectra matrices (nsp) to add
	nsp = 6
	if (len(nrcta) == 1):	# total has been sought
		print("n,all .....")
		for i in range(nsp):
			ifile.readline()
			temporary = [0]*numdata
			temporary = line_type4_info(ifile,nre,numdata)
			k1 = 0
			for j in range (nre):
				for k in range (nre):
					pkaind[j][k] = temporary[k1]
					k1 += 1
			for j in range (nre):
				for k in range (nre):
					dsdt[j][k] = dsdt[j][k] + pkaind[j][k]

	if (len(nrcta) > 1): 	# sum of some partials has been sought
		print("Sum of sought partial reactions .....")
		for i in range(nsp):
			flag = 0
			for j in range(len(nrcta)):
				if (i+1 == nrcta[j]):
					flag = 1
					break
			if (flag == 0):
				ifile.readline()
				temporary = [0]*numdata
				temporary = line_type4_info(ifile,nre,numdata)
			if (flag == 1):
				ifile.readline()
				temporary = [0]*numdata
				temporary = line_type4_info(ifile,nre,numdata)
				k1 = 0
				for j in range (nre):
					for k in range (nre):
						pkaind[j][k] = temporary[k1]
						k1 += 1
				for j in range (nre):
					for k in range (nre):
						dsdt[j][k] = dsdt[j][k] + pkaind[j][k]
	return(dsdt)

	ifile.close()
	#os.remove('ToAddAll.txt')

#=======Make the required mesh and call reactions=======*

	# Creates the finer energy array from the chosen energy group
	# structure by introducing more number of equispaced energy points
	# inside each group. This fine energy array of increased length is
	# used for both the incident neutrons and recoil nuclei.
	# After creating this finer array it also calls the reaction specific
	# subroutines according to the input.

def FINE_ENERGY_CALL_REAC (insp,eliso,igtype,nrg,nbpoints,irct,nrcta):

	nre = nrg - 1
	nbge = nrg + (nbpoints*nre)
	fct = 1.0/(nbpoints+1)
	ret = [0]*nbge; erg = [0]*nrg
	if (igtype == 0):
		ifile = open('Energy-GroupLimits.txt', 'r')
		nrg = int(ifile.readline().split()[0])
		for i in range (nre,-1,1):
			erg[i] = float(ifile.readline().split()[0])
		ifile.close()

	if (igtype==1): 
		erg = engrp1()
	if (igtype==2):
		erg = engrp2()
	if (igtype==3):
		erg = engrp3()
	if (igtype==4):
		erg = engrp4()
	if (igtype==5):
		erg = engrp5()
	if (igtype==6):
		erg = engrp6()
	if (igtype==7):
		erg = engrp7()
	if (igtype==8):
		erg = engrp8()
	if (igtype==9):
		erg = engrp9()
	if (igtype==10):
		erg = engrp10()
	if (igtype==11):
		erg = engrp11()
	if (igtype==12):
		erg = engrp12()

	j = 0
	for i in range (nre):
		x = fct*(erg[i+1]-erg[i])
		for k in range (nbpoints+1):
			ret[j] = k*x + erg[i]
			j = j + 1
	ret[nbge-1] = erg[nrg-1]

	if (irct == 1):
		PKAS_ELASTIC (insp,eliso,ret,nbge,nbpoints,nre,igtype)
	if (irct == 2):
		PKAS_INELASTIC (insp,eliso,ret,nbge,nbpoints,nre,igtype)
	if (irct == 3):
		CONTROL_nCPO (insp,eliso,ret,nbge,nbpoints,nre,igtype)
	if (irct == 4):
		CONTROL_nxn (insp,eliso,ret,nbge,nbpoints,nre,igtype)
	if (irct == 5):
		PKAS_ng (insp,eliso,ret,nbge,nbpoints,nre,igtype)
	if (irct == 6):
		PKAS_redtnMF6MT5 (insp,eliso,ret,nbge,nbpoints,nre,igtype)
	if (irct == 7):
		ofile1001 = open('n-allPKAspectra.txt', 'a')
		print(eliso, file = ofile1001)
		PKAS_ELASTIC (insp,eliso,ret,nbge,nbpoints,nre,igtype)
		PKAS_INELASTIC (insp,eliso,ret,nbge,nbpoints,nre,igtype)
		CONTROL_nCPO (insp,eliso,ret,nbge,nbpoints,nre,igtype)
		CONTROL_nxn (insp,eliso,ret,nbge,nbpoints,nre,igtype)
		PKAS_ng (insp,eliso,ret,nbge,nbpoints,nre,igtype)
		PKAS_redtnMF6MT5 (insp,eliso,ret,nbge,nbpoints,nre,igtype)
		# for the total of contributions from all reactions
		dsdt = numpy.zeros((nre,nre))
		dsdt = ALLSUM (nrg,nrcta)
		print('calculated total PKA spectra .....')
		print('writing .....')
		print('The total recoil nuclei energy spectra', file = ofile1001)
		for it in range (nre):
			print (['{:.6E}'.format(dsdt[it][jt]) for jt in range (nre)], file = ofile1001)
		ofile1001.close()

#=======Integrate in the mesh within groups=======*

	# The PKA spectra values obtained in the finer energy points are
	# integrated within each of the original energy groups separately
	# to find the PKA spectra in that group structure.
	
def GROUP_INTEG (ret,nbge,nbpoints,nre,dsgmdT):

	pkamatr = numpy.zeros((nre,nre))
	
	print( 'group .....')
	l = 0														# l -- controls incident group in pka matrix
	for iine in range (0, nbge-(nbpoints+1)+1, (nbpoints+1)):		# iine -- do for all group incident energies
		k = 0													# k -- controls recoil group in pka matrix
		for iree in range (0, nbge-(nbpoints+1)+1, (nbpoints+1)):  # iree -- do for all group recoil energies
			s1 = 0
			s2 = 0
			for j1 in range (iine, (iine + nbpoints+1)):
				dE = ret[j1+1] - ret[j1]
				s = 0
				s3 = 0
				for j in range (iree, (iree + nbpoints+1)):
					dEp = ret[j+1] - ret[j]
					s = s + dsgmdT[j1][j] * dEp
					s3 = s3 + dEp

				if (s3 != 0):
					s = s/s3
				s1 = s1 + s*dE
				s2 = s2 + dE

			if (s2 != 0):
				pkamatr[l][k] = abs(s1/s2)
			k = k + 1

		l = l + 1
	return(pkamatr)
		
#=======Normalize the PKA spectrum to cross section=======*

	# The PKA spectra obtained from GROUP_INTEG for each original energy
	# group are normalised such that when they are integrated over the
	# whole energy range of the recoil nucleus they yield the multigrouped
	# basic neutron cross section for a particular incident neutron energy
	# group.
	
def PKAS_NORM_CROSSSEC (insp,MTtg,igtype,nrg,dsgmdT):
	
	nre = nrg - 1
	erg = [0]*nrg
	gsig = [0]*nre
	
	print('normalize .....')
				
	if (igtype == 0):
		ifile = open('Energy-GroupLimits.txt', 'r')
		nrg = int(ifile.readline().split()[0])
		for i in range (nre,-1,1):
			erg[i] = float(ifile.readline().split()[0])
		ifile.close()

	if (igtype==1): 
		erg = engrp1()
	if (igtype==2):
		erg = engrp2()
	if (igtype==3):
		erg = engrp3()
	if (igtype==4):
		erg = engrp4()
	if (igtype==5):
		erg = engrp5()
	if (igtype==6):
		erg = engrp6()
	if (igtype==7):
		erg = engrp7()
	if (igtype==8):
		erg = engrp8()
	if (igtype==9):
		erg = engrp9()
	if (igtype==10):
		erg = engrp10()
	if (igtype==11):
		erg = engrp11()
	if (igtype==12):
		erg = engrp12()

	gsig = GROUPMULTI (insp,igtype,MTtg,nrg)

	for i in range (nre):
		if (gsig[i] != 0):
			sumnorm = 0
			for j in range (nre):
				sumnorm = sumnorm + dsgmdT[i][j]*(erg[j+1]-erg[j])
			for j in range (nre):
				if (sumnorm != 0):
					dsgmdT[i][j] = dsgmdT[i][j]*gsig[i]/sumnorm

	return(dsgmdT)

#=======Print the PKA matrices=======*

	# The PKA spectra from partial reactions are printed in text
	# file in matrix form under the heading of the particular MT 
	# value.
	
def PRINTPKAS (eliso,MTtg,nre,pkamatr):
	ofile500 = open('PKA-MATRICES.txt', 'a')
	print (eliso, file = ofile500)
	print ( 'writing .....')
		#if (MTtg==2) then
		#open(unit=500,file='n-nPKAspectra.txt',position='append')
		#end if
		#if (MTtg==4) then
		#open(unit=500,file="n-n'PKAspectra.txt",position='append')
		#end if
		#if (MTtg==3) then
		#open(unit=500,file='n-thPKAspectra.txt',position='append')
		#end if
		#if (MTtg==16) then
		#open(unit=500,file='n-2nPKAspectra.txt',position='append')
		#end if
	print('MTPKA = ', MTtg, 'PKA spectra', file = ofile500)
	for i in range (nre):
		print (['{:.6E}'.format(pkamatr[i][j]) for j in range (nre)], file = ofile500)
	
	ofile500.close()

#=======Contributions from thresholds=======*

	# The controlling of the calculation of PKA spectra due to (n,CPO)
	# reactions is performed here. Corresponding to each of this type 
	# of reaction the subroutine PKAS_nCPO is called for the calculation.
			
def CONTROL_nCPO (insp,eliso,ret,nbge,nbpoints,nre,igtype):

	print("n,CPO .....")
	pka3l = numpy.zeros((nbge,nbge)); tmpkaMT103 = numpy.zeros((nbge,nbge))
	tmpkaMT104 = numpy.zeros((nbge,nbge)); tmpkaMT105 = numpy.zeros((nbge,nbge))
	tmpkaMT106 = numpy.zeros((nbge,nbge)); tmpkaMT107 = numpy.zeros((nbge,nbge))
	pka3 = numpy.zeros((nre,nre)); pkaMT103 = numpy.zeros((nre,nre)); pkaMT104 = numpy.zeros((nre,nre))
	pkaMT105 = numpy.zeros((nre,nre)); pkaMT106 = numpy.zeros((nre,nre))
	pkaMT107 = numpy.zeros((nre,nre)); pkamatr = numpy.zeros((nre,nre))

	MTnum = [0]*281   		# 31 + 250
	nrg = nre + 1
	
	MTnum[0]=11; MTnum[1]=22;MTnum[2]=23;MTnum[3]=24
	MTnum[4]=25; MTnum[5]=28;MTnum[6]=29;MTnum[7]=30
	MTnum[8]=32; MTnum[9]=33;MTnum[10]=34;MTnum[11]=35
	MTnum[12]=36; MTnum[13]=41;MTnum[14]=42;MTnum[15]=44
	MTnum[16]=45; MTnum[17]=103;MTnum[18]=104;MTnum[19]=105
	MTnum[20]=106;MTnum[21]=107;MTnum[22]=108;MTnum[23]=109
	MTnum[24]=111;MTnum[25]=112;MTnum[26]=113;MTnum[27]=114
	MTnum[28]=115;MTnum[29]=116;MTnum[30]=117

	for i in range (31, 281):
		MTnum[i] = i + 569
	
	ifdisc103=0; ifdisc104=0; ifdisc105=0; ifdisc106=0; ifdisc107=0
	ifdisdata103=0; ifdisdata104=0; ifdisdata105=0; ifdisdata106=0
	ifdisdata107=0
		
	for j in range (281):
		iflMTpr = 0
		MTfind = MTnum[j]	
		iflMTpr = FindMT(MTfind)
		if (iflMTpr == 1):
			(pka3l,ifdpd,iflpresent) = PKAS_nCPO (MTnum[j],j,ret,nbge)
			if (iflpresent == 1):
				print ('n-threshold: MT = ', MTnum[j])
				MTc = MTnum[j]			
				if (MTc==103):
					tmpkaMT103 = pka3l
				if (MTc==104):
					tmpkaMT104 = pka3l
				if (MTc==105):
					tmpkaMT105 = pka3l
				if (MTc==106):
					tmpkaMT106 = pka3l
				if (MTc==107):
					tmpkaMT107 = pka3l
				if (MTc != 103 and MTc != 104 and MTc != 105 and MTc != 106 and MTc != 107):
					if (MTc < 600):
						pkamatr = GROUP_INTEG(ret,nbge,nbpoints,nre,pka3l)
						pkamatr = PKAS_NORM_CROSSSEC (insp,MTc,igtype,nrg,pkamatr)
						PRINTPKAS (eliso,MTc,nre,pkamatr)
						for isum in range (nre):
							for jsum in range (nre):
								pka3[isum][jsum] = pka3[isum][jsum] + pkamatr[isum][jsum]
			
			
				if (600 <= MTc and MTc <= 649):
					ifdisdata103 = 1
					if (ifdpd == 1):
						ifdisc103 = 1
						pkamatr = GROUP_INTEG (ret,nbge,nbpoints,nre,pka3l)
						pkamatr = PKAS_NORM_CROSSSEC (insp,MTc,igtype,nrg,pkamatr)
						PRINTPKAS (eliso,MTc,nre,pkamatr)
						for isum in range (nre):
							for jsum in range (nre):
								pkaMT103[isum][jsum] = pkaMT103[isum][jsum] + pkamatr[isum][jsum]
					
					if (ifdpd == 0 and ifdisc103 == 1):
						pkamatr = GROUP_INTEG (ret,nbge,nbpoints,nre,pka3l)
						pkamatr = PKAS_NORM_CROSSSEC (insp,MTc,igtype,nrg,pkamatr)
						PRINTPKAS (eliso,MTc,nre,pkamatr)		
						for isum in range (nre):
							for jsum in range (nre):
								pkaMT103[isum][jsum] = pkaMT103[isum][jsum] + pkamatr[isum][jsum]
					
					if (ifdpd == 0 and ifdisc103 == 0):
						pkamatr = GROUP_INTEG(ret,nbge,nbpoints,nre,tmpkaMT103)
						pkamatr = PKAS_NORM_CROSSSEC (insp,MTc,igtype,nrg,pkamatr)
						pkaMT103 = pkamatr 
	
					if (MTc == 649):
						print ('n-threshold: MT = ', 103)
						PRINTPKAS (eliso,103,nre,pkaMT103)
				#-------------------
				
				if (650 <= MTc and MTc <= 699):
					ifdisdata104 = 1
					if (ifdpd == 1):
						ifdisc104 = 1
						pkamatr = GROUP_INTEG (ret,nbge,nbpoints,nre,pka3l)
						pkamatr = PKAS_NORM_CROSSSEC (insp,MTc,igtype,nrg,pkamatr)        
						PRINTPKAS (eliso,MTc,nre,pkamatr)
						for isum in range (nre):
							for jsum in range (nre):
								pkaMT104[isum][jsum] = pkaMT104[isum][jsum] + pkamatr[isum][jsum]
	
					if (ifdpd == 0 and ifdisc104 == 1):
						pkamatr = GROUP_INTEG(ret,nbge,nbpoints,nre,pka3l)
						pkamatr = PKAS_NORM_CROSSSEC (insp,MTc,igtype,nrg,pkamatr)        
						PRINTPKAS (eliso,MTc,nre,pkamatr)		
						for isum in range (nre):
							for jsum in range (nre):
								pkaMT104[isum][jsum] = pkaMT104[isum][jsum] + pkamatr[isum][jsum]
	
					if (ifdpd == 0 and ifdisc104 == 0):
						pkamatr = GROUP_INTEG(ret,nbge,nbpoints,nre,tmpkaMT104)
						pkamatr = PKAS_NORM_CROSSSEC (insp,MTc,igtype,nrg,pkamatr)
						pkaMT104 = pkamatr
	
					if (MTc == 699):
						print ('n-threshold: MT = ', 104)
						PRINTPKAS (eliso,104,nre,pkaMT104)
				#-------------------
	
				if (700 <= MTc and MTc <= 749):
					ifdisdata105 = 1
					if (ifdpd == 1):
						ifdisc105 = 1
						pkamatr = GROUP_INTEG(ret,nbge,nbpoints,nre,pka3l)
						pkamatr = PKAS_NORM_CROSSSEC (insp,MTc,igtype,nrg,pkamatr)
						PRINTPKAS (eliso,MTc,nre,pkamatr)
						for isum in range (nre):
							for jsum in range (nre):
								pkaMT105[isum][jsum] = pkaMT105[isum][jsum] + pkamatr[isum][jsum]
					
					if (ifdpd == 0 and ifdisc105 == 1):
						pkamatr = GROUP_INTEG(ret,nbge,nbpoints,nre,pka3l)
						pkamatr = PKAS_NORM_CROSSSEC (insp,MTc,igtype,nrg,pkamatr)
						PRINTPKAS (eliso,MTc,nre,pkamatr)
						for isum in range (nre):
							for jsum in range (nre):
								pkaMT105[isum][jsum] = pkaMT105[isum][jsum] + pkamatr[isum][jsum]
		
					if (ifdpd == 0 and ifdisc105 == 0):
						pkamatr = GROUP_INTEG(ret,nbge,nbpoints,nre,tmpkaMT105)	
						pkamatr = PKAS_NORM_CROSSSEC (insp,MTc,igtype,nrg,pkamatr)
						pkaMT105 = pkamatr 
					
					if (MTc==749):
						print ('n-threshold: MT = ', 105)
						PRINTPKAS (eliso,105,nre,pkaMT105)
				#-------------------
			
				if (750 <= MTc and MTc <= 799):
					ifdisdata106 = 1
					if (ifdpd == 1):
						ifdisc106 = 1
						pkamatr = GROUP_INTEG(ret,nbge,nbpoints,nre,pka3l)
						pkamatr = PKAS_NORM_CROSSSEC (insp,MTc,igtype,nrg,pkamatr)
						PRINTPKAS (eliso,MTc,nre,pkamatr)
						for isum in range (nre):
							for jsum in range (nre):
								pkaMT106[isum][jsum] = pkaMT106[isum][jsum] + pkamatr[isum][jsum]
					
					if (ifdpd == 0 and ifdisc106 == 1):
						pkamatr = GROUP_INTEG(ret,nbge,nbpoints,nre,pka3l)
						pkamatr = PKAS_NORM_CROSSSEC (insp,MTc,igtype,nrg,pkamatr)
						PRINTPKAS (eliso,MTc,nre,pkamatr)
						for isum in range (nre):
							for jsum in range (nre):
								pkaMT106[isum][jsum] = pkaMT106[isum][jsum] + pkamatr[isum][jsum]
					
					if (ifdpd == 0 and ifdisc106 == 0):
						pkamatr = GROUP_INTEG(ret,nbge,nbpoints,nre,tmpkaMT106)
						pkamatr = PKAS_NORM_CROSSSEC (insp,MTc,igtype,nrg,pkamatr)		
						pkaMT106 = pkamatr
					
					if (MTc == 799):
						print ('n-threshold: MT = ', 106)
						PRINTPKAS (eliso,106,nre,pkaMT106)
				#-------------------
	
				if (800 <= MTc and MTc <= 849):
					ifdisdata107 = 1
					if (ifdpd == 1):
						ifdisc107 = 1
						pkamatr = GROUP_INTEG(ret,nbge,nbpoints,nre,pka3l)
						pkamatr = PKAS_NORM_CROSSSEC (insp,MTc,igtype,nrg,pkamatr)
						PRINTPKAS (eliso,MTc,nre,pkamatr)
						for isum in range (nre):
							for jsum in range (nre):
								pkaMT107[isum][jsum] = pkaMT107[isum][jsum] + pkamatr[isum][jsum]
	
					if (ifdpd == 0 and ifdisc107 == 1):
						pkamatr = GROUP_INTEG(ret,nbge,nbpoints,nre,pka3l)
						pkamatr = PKAS_NORM_CROSSSEC (insp,MTc,igtype,nrg,pkamatr)
						PRINTPKAS (eliso,MTc,nre,pkamatr)
						for isum in range (nre):
							for jsum in range (nre):
								pkaMT107[isum][jsum] = pkaMT107[isum][jsum] + pkamatr[isum][jsum]
	
					if (ifdpd == 0 and ifdisc107 == 0):
						pkamatr = GROUP_INTEG(ret,nbge,nbpoints,nre,tmpkaMT107)
						pkamatr = PKAS_NORM_CROSSSEC (insp,MTc,igtype,nrg,pkamatr)		
						pkaMT107 = pkamatr
					
					if (MTc==849):
						print ('n-threshold: MT = ', 107)
						PRINTPKAS (eliso,107,nre,pkaMT107)
				#-------------------
	
	if (ifdisdata103 == 0):
		iflpr = prflag(nbge,tmpkaMT103)
		if (iflpr == 1):
			pkamatr = GROUP_INTEG(ret,nbge,nbpoints,nre,tmpkaMT103)
			pkamatr = PKAS_NORM_CROSSSEC (insp,103,igtype,nrg,pkamatr)
			pkaMT103 = pkamatr
			PRINTPKAS (eliso,103,nre,pkaMT103)

	if (ifdisdata104 == 0):
		iflpr = prflag(nbge,tmpkaMT104)
		if (iflpr == 1):
			pkamatr = GROUP_INTEG(ret,nbge,nbpoints,nre,tmpkaMT104)
			pkamatr = PKAS_NORM_CROSSSEC (insp,104,igtype,nrg,pkamatr)
			pkaMT104 = pkamatr
			PRINTPKAS (eliso,104,nre,pkaMT104)

	if (ifdisdata105 == 0):
		iflpr = prflag(nbge,tmpkaMT105)
		if (iflpr == 1):
			pkamatr = GROUP_INTEG(ret,nbge,nbpoints,nre,tmpkaMT105)
			pkamatr = PKAS_NORM_CROSSSEC (insp,105,igtype,nrg,pkamatr)
			pkaMT105 = pkamatr
			PRINTPKAS (eliso,105,nre,pkaMT105)
		
	if (ifdisdata106 == 0):
		iflpr = prflag(nbge,tmpkaMT106)
		if (iflpr == 1):
			pkamatr = GROUP_INTEG(ret,nbge,nbpoints,nre,tmpkaMT106)
			pkamatr = PKAS_NORM_CROSSSEC (insp,106,igtype,nrg,pkamatr)
			pkaMT106 = pkamatr
			PRINTPKAS (eliso,106,nre,pkaMT106)

	if (ifdisdata107 == 0):
		iflpr = prflag(nbge,tmpkaMT107)
		if (iflpr == 1):
			pkamatr = GROUP_INTEG(ret,nbge,nbpoints,nre,tmpkaMT107)
			pkamatr = PKAS_NORM_CROSSSEC (insp,107,igtype,nrg,pkamatr)
			pkaMT107 = pkamatr
			PRINTPKAS (eliso,107,nre,pkaMT107)
	
	for isum in range (nre):
		for jsum in range (nre):
			pka3[isum][jsum] = pka3[isum][jsum] + pkaMT103[isum][jsum] + \
			pkaMT104[isum][jsum] + pkaMT105[isum][jsum] + pkaMT106[isum][jsum]+ pkaMT107[isum][jsum]

	print ('n-threshold: MT = ', 3001, '(Total (n, CPO))')
	PRINTPKAS (eliso,3001,nre,pka3)
# ------------------------------------------------------------------		
	ofile1000 = open ('ToAddAll.txt', 'a')
	print('', file = ofile1000)
	for i in range (nre):
		print (['{:.6E}'.format(pka3[i][j]) for j in range (nre)], file = ofile1000)

	ofile1000.close()

#=======Contributions from (n, xn) reactions=======*

	# The controlling of the calculation of PKA spectra due to (n,xn)
	# reactions is performed here. Corresponding to each of this type 
	# of reaction the subroutine PKAS_nxn is called for the calculation.
	# (n,2n), (n,3n) and (n,4n) only.
	
def CONTROL_nxn (insp,eliso,ret,nbge,nbpoints,nre,igtype):

	print("n,xn .....")
	pkaxn1l = numpy.zeros((nbge,nbge))
	pkaxn1 = numpy.zeros((nre,nre)); pkamatr = numpy.zeros((nre,nre))
	MTnum = [0]*3
	nrg = nre + 1
	MTnum[0] = 16; MTnum[1] = 17; MTnum[2] = 37
	
	for i in range (3):
		iflMTpr = 0
		MTfind = MTnum[i]
		iflMTpr = FindMT(MTfind)
		if (iflMTpr == 1):
			(pkaxn1l, iflpresent) = PKAS_nxn (MTnum[i],ret,nbge,nbpoints,nre,igtype)
			if (iflpresent == 1):
				MTc = MTnum[i]
				print ('n, xn:: MT = ', MTc)
				pkamatr = GROUP_INTEG(ret,nbge,nbpoints,nre,pkaxn1l)
				pkamatr = PKAS_NORM_CROSSSEC (insp,MTc,igtype,nrg,pkamatr)
				PRINTPKAS (eliso,MTc,nre,pkamatr)
				for isum in range (nre):
					for jsum in range (nre):
						pkaxn1[isum][jsum] = pkaxn1[isum][jsum] + pkamatr[isum][jsum]

	PRINTPKAS (eliso,1601,nre,pkamatr)
	
	ofile1000 = open ('ToAddAll.txt', 'a')
	print ('', file = ofile1000)
	
	for i in range (nre):
		print (['{:.6E}'.format(pkaxn1[i][j]) for j in range (nre)], file = ofile1000)
	
	ofile1000.close()

#=======To find for availble MT=======*
	
	# Finds if a particular reaction is given in the present evaluation
	# from the directory available in File 1.
	
def FindMT (MTfind):
	MFs = [0]*1000; MTs = [0]*1000
	(nfiles,MFs,MTs) = FILE1()
	iflMTpr = 0
	for i in range (nfiles):
		if (MFs[i] == 3 and MTs[i] == MTfind):
			iflMTpr = 1
			break
	return(iflMTpr)

#=======Present flag=======*

	# Finds if the array of the PKA specta matrix contains
	# any non-zero element.
	
def prflag (n,arr):
	iflpr = 0
	for i in range (n):
		for j in range (n):
			if (arr[i][j] != 0):
				iflpr = 1
				break
	return(iflpr)

# ===========Kernel of Elastic Scattering===========*

	# Calculates the kernel for elastic scattering
	# from Legendre expansion coefficients.

def SPKAEL (alc1,NLa,ret,nbge,A,Z,En,y):
	Pl1 = [0]*NLa
	val2 = [0]*nbge
	Emx = y*En
	for i in range (nbge):
		if (Emx != 0):
			xabc = 1 - 2*ret[i]/Emx
			if (abs(xabc) <= 1.0):	 	#-1.0<=xabc .and. 
				Pl1[0] = 1
				Pl1[1] = xabc
				for j in range (2, NLa):
					Pl1[j] = (((2*(j-1)+1)*xabc*Pl1[j-1])-((j-1)*Pl1[j-2]))/j
				p1 = 0
				for l in range (NLa):
					p1 = p1 + (((2*l)+1)*Pl1[l]*alc1[l])/2
				val2[i] = p1

	return(val2)
				
	# Calculates the kernel for elastic scattering
	# from tabulated mu vs. f(mu,E) data.
		
def SPKAEL1 (xgaussq,fmuE,ret,nbge,A,Z,En,y):
	val2 = [0]*nbge
	Emx = y*En
	for i in range (nbge):
		if (Emx != 0):
			xabc = 1 - 2*ret[i]/Emx
		for k in range (64):	#64-1
			if (xabc == xgaussq[k]):
				val2[i] = fmuE[k][i]
				break
			if(xabc > xgaussq[k] and xabc <= xgaussq[k+1]):
				x = xabc
				x1 = xgaussq[k]
				x2 = xgaussq[k+1]
				y1 = fmuE[k][i]
				y2 = fmuE[k+1][i]
				if(y1 != 0 and (x2-x1) != 0):
					val2[i] = y1*math.exp((x-x1)*math.log(y2/y1)/(x2-x1))
				break
	return(val2)

#=======Uniue energy Array=========*

# Making unique energy array from MF3 MT1
	
def UQCE():
	# Et=Energy array in MT=1, Etu=unique of Et
	# extraction of total energy points
	
	ifile104 = open ('tape02', 'r')
	while True:
		line = ifile104.readline()  
		data = eachlineinfo(line)
		MAT = int(data[6]); MF = int(data[7]); MT  = int(data[8])
		if (MAT != -1):
			if (MF == 3):
				if (MT == 1):
					line = ifile104.readline() 
					data = eachlineinfo(line)
					QM = float(data[0]); QI =  float(data[1]); NR = int(data[4]); NPt = int(data[5])
					Et = [0]*NPt
					LR = int(ifile104.readline().split()[1])
					temporary1 = [0]*NPt
					temporary2 = [0]*NPt
					(temporary1,temporary2) = line_type3_info(ifile104,NPt,2)
					for N, (value1, value2) in enumerate(zip(temporary1, temporary2), 0):
						Et[N] = value1
						s = value2
		else:
			break
	ifile104.close()

	# make unique common energy
	
	Etu = numpy.array(Et)
	Etu = numpy.unique(Etu)
	NPt = len(Etu)
	
	return(Etu, NPt)
#---------------------------------------------------------------------------

# Often unterpolation function

def crstd(x,x1,x2,y1,y2):
	if((x2-x1) != 0):
		y = y1+((y2-y1)*(x-x1)/(x2-x1))
	return(y)

#=======Elastic PKA Spectra=======*

	# Calculation of PKA spectra due to the elastic scattering 
	# interaction of neutron.
		
def PKAS_ELASTIC (insp,eliso,ret,nbge,nbpoints,nre,igtype):

	print ("n,n .....")
	sret = [0]*nbge
	dsgmdT = numpy.zeros((nbge,nbge))
	nquad = 64
	xgaussq = [0]*nquad
	wgaussq = [0]*nquad
	pkamatr = numpy.zeros((nre,nre))
	nrg = nre + 1
	ofile101 = open("Output_RadEMC-RecedU.txt", 'a')

	print(' Elastic scattering (MT = 2) PKA spectra', file = ofile101)
	print('-------------------------------------------', file = ofile101)
# --------------------------------------------------------------------------
	# extraction of total energy grid
	ifile = open ("tape02", 'r')
	while True:
		line = ifile.readline()
		if (line == ''):
			break
		data = eachlineinfo(line)
		MAT = int(data[6]); MF = int(data[7]); MT  = int(data[8])
		if (MAT != -1):
			if (MF == 3):
				if (MT == 1):
					line = ifile.readline()
					data = eachlineinfo(line)
					QM = float(data[0]); QI =  float(data[1])
					NR = int(data[4]); NPt = int(data[5])
					Et = [0]*NPt; Etu = [0]*NPt; siget = [0]*NPt
					alfull = numpy.zeros((nbge,65)); fmuE = numpy.zeros((64,nbge))
					LR = int(ifile.readline().split()[1])
					(Et, sigt) = line_type3_info(ifile,NPt,2)
		else:
			break
	ifile.close()

	print('', file = ofile101)
	print(NPt, ' Energy points in total cross section', file = ofile101)

	# extraction of cross sections
		
	ifile = open ("tape02", 'r')
	while True:
		line = ifile.readline()
		if (line == ''):
			break
		data = eachlineinfo(line)
		MAT = int(data[6]); MF = int(data[7]); MT  = int(data[8])

		if (MAT != -1):
			if (MF == 3):
				if (MT == 2):
					line = ifile.readline()
					data = eachlineinfo(line)
					QM = float(data[0]); QI =  float(data[1])
					NR = int(data[4]); NP = int(data[5])
					E = [0]*NP; sig = [0]*NP
					LR = int(ifile.readline().split()[1])
					(E, sig) = line_type3_info(ifile,NP,2)
		else:
			break
	ifile.close()

	print ('', file = ofile101)
	print (NP, 'Energy points in elastic scattering', file = ofile101)
		
 	# extraction of Legendre polynomial coefficients and Probability
	
	ifile = open ("tape01", 'r')
	while True:
		line = ifile.readline()
		if (line == ''):
			break
		data = eachlineinfo(line)
		MAT = int(data[6]); MF = int(data[7]); MT  = int(data[8])	
		if (MF == 3 and MT == 0):		#.and. NS==99999
			line = ifile.readline()
			data = eachlineinfo(line)
			MAT = int(data[6]); MF = int(data[7]); MT  = int(data[8])
			line = ifile.readline()
			(ZA,AWR,l1,LTT,NK,l2,MAT,MF,MT) = line_type1_info(line)

		if (MAT != -1):
			if (MF == 4):
				if (MT == 2):
					ifile.readline()
					if (LTT == 0):
						alc = [0]*65
					if (LTT == 3 or LTT == 1):
						line = ifile.readline()
						data = eachlineinfo(line)
						NE1 = int(data[5])
						# Legendre polynomial coefficients
						ifile.readline()
						NLarray = [0]*NE1; al = numpy.zeros((NE1,65)); EL = [0]*NE1; alc = [0]*65
						for i in range (NE1):
							line = ifile.readline()
							data = eachlineinfo(line)
							T = float(data[0]); EL[i] = float(data[1]); NL = int(data[4])
							NLarray[i] = NL + 1
							al[i][0] = 1
							temporary = [0]*NLarray[i]
							temporary = line_type3_info(ifile,NLarray[i]-1,1)
							for j, value in enumerate(temporary, 1):
								al[i][j] = value

					if (LTT == 3 or LTT == 2):
						line = ifile.readline()
						data = eachlineinfo(line)
						NE2 = int(data[5])
						# Tabulated Probability
						ifile.readline()
						Enf = [0]*NE2; NPr = [0]*NE2; cdata = numpy.zeros((200,NE2))
						fdata = numpy.zeros((200,NE2)); ftotal = numpy.zeros((64,NE2))
						for i in range (NE2):
							line = ifile.readline()
							data = eachlineinfo(line)
							T = float(data[0]); Enf[i] = float(data[1])
							line = ifile.readline()
							data = eachlineinfo(line)
							NPr[i] = int(data[0])
							temporary1 = [0]*NPr[i]
							temporary2 = [0]*NPr[i]
							(temporary1,temporary2) = line_type3_info(ifile,NPr[i],2)
							for j, (value1, value2) in enumerate(zip(temporary1, temporary2), 0):
								cdata[j][i] = value1
								fdata[j][i] = value2
		else:
			break

	ifile.close()

# --------------------------------------------------------------------------
	A = AWR
	Z = int(ZA/1000)
	y = 4*A/((A+1)**2)
	(xgaussq, wgaussq) = GQ()                                

	# Make unique common energy in MF=3, MT=1

	Etu = numpy.array(Et)
	Etu = numpy.unique(Etu)
	NPt = len(Etu)
	print('', file = ofile101)
	print(NPt,' Unique total energy points', file = ofile101)

	# Finding (n,n) cross section in the unique energy array

	siget = trptuqce (E, sig, Etu)

	# Finding cross sections corresponding to the fine energy group array

	sret = trptuqce (Etu, siget, ret)

	print('', file = ofile101)

	if(LTT == 3):
		print('Legendre coefficients and tabulated probability', file = ofile101)
		print('data representations of angular distribution', file = ofile101)

	if(LTT == 1):
		print('Legendre coefficients representation of angular distribution', file = ofile101)

	if(LTT == 2):
		print('Tabulated probability data representation of angular distribution', file = ofile101)
		
	# call GROUPMULTI(igtype,2,nbge,sret)

	# LOG- LINEAR INTERPOLATION BETWEEN MU AND F(MU,E)
		
	if (LTT == 3 or LTT == 2):
		for k in range (NE2):
			for i in range (64):
				for j in range (int(NPr[k])):
					if (xgaussq[i] == cdata[j][k]):
						ftotal[i][k] = fdata[j][k]
						break
					if(xgaussq[i] > cdata[j][k] and xgaussq[i] < cdata[j+1][k]):
						x = xgaussq[i]
						x1 = cdata[j][k]
						x2 = cdata[j+1][k]
						y1 = fdata[j][k]
						y2 = fdata[j+1][k]
						if (y1 != 0 and (x2-x1) != 0):
							ftotal[i][k] = y1*math.exp((x-x1)*math.log(y2/y1)/(x2-x1))
						break
 		
		# TO GET THE F(MU,E) FOR THE FULL INCIDENT ENERGY RANGE

		for i in range (nbge):
			for j in range (64):
				fmuE[j][i] = 0.5

		for i in range (nbge):
			for j in range (NE2-1):
				if (ret[i] == Enf[j]):
					for k in range (64):
						fmuE[k][i] = ftotal[k][j]
					break
				if (ret[i] > Enf[j] and ret[i] <= Enf[j+1]):
					x = ret[i]
					x1 = Enf[j]
					x2 = Enf[j+1]
					for k in range (64):
						y1 = ftotal[k][j]
						y2 = ftotal[k][j+1]
						if ((x2-x1) != 0):
							fmuE[k][i] = y1 + ((y2-y1)*(x-x1)/(x2-x1))
					break

	# TO GET THE AL(E) FOR THE FULL ENERGY RANGE

	for i in range (nbge):
		alfull[i][0] = 1
		for k in range (1, 65):
			alfull[i][k] = 0

	for i in range (nbge):
		for j in range (NE1-1):
			if (ret[i] == EL[j]):
				for k in range (1, 65):
					alfull[i][k] = al[j][k]
				break
			if (ret[i] > EL[j] and ret[i] <= EL[j+1]):
				diff1 = ret[i] - EL[j]
				diff2 = ret[i] - ret[i-1]
				if (diff1 <= diff2):
					for k in range (1, 65):
						x = ret[i]
						x1 = EL[j]
						x2 = EL[j+1]
						y1 = al[j][k]
						y2 = al[j+1][k]
						if((x2-x1) != 0):
							alfull[i][k] = y1 + ((y2-y1)*(x-x1)/(x2-x1))

				if (diff2 < diff1):
					for k in range (1, 65):
						x = ret[i]
						x1 = ret[i-1]
						x2 = EL[j+1]
						y1 = alfull[i-1][k]
						y2 = al[j+1][k]
						if((x2-x1) != 0):
							alfull[i][k] = y1 + ((y2-y1)*(x-x1)/(x2-x1))

				for k in range (1, 65):
					if (alfull[i][k] == 0):
						alfull[i][k] = alfull[i-1][k]
				break
# --------------------------------------------------------------------------
	# val2 array-- gets the kernel k(E,ER)
	# dsgmdT array-- stores the dsigma/dER
	
	for i in range (nbge):
		if ((LTT == 3 and ret[i] <= EL[NE1-1]) or LTT == 1 or LTT == 0):
			for k in range (65):
				alc[k] = alfull[i][k]
			val2 = [0]*nbge
			val2 = SPKAEL (alc,65,ret,nbge,A,Z,ret[i],y)
			for ia in range (nbge):
				dsgmdT[i][ia] = abs(sret[i]*val2[ia])

		if((LTT == 3 and ret[i] > Enf[0]) or LTT == 2):
			val2 = SPKAEL1 (xgaussq,fmuE,ret,nbge,A,Z,ret[i],y)
			for ia in range (nbge):
				dsgmdT[i][ia] = abs(sret[i]*val2[ia])
	
	pkamatr = GROUP_INTEG (ret,nbge,nbpoints,nre,dsgmdT)	
	pkamatr = PKAS_NORM_CROSSSEC (insp,2,igtype,nrg,pkamatr)
	PRINTPKAS (eliso,2,nre,pkamatr)

	ofile1000 = open ('ToAddAll.txt', 'a')
	print('', file = ofile1000)
	for i in range (nre):
		print (['{:.6E}'.format(pkamatr[i][j]) for j in range (nre)], file = ofile1000)
	
	ofile1000.close()
	ofile101.close()

	# PKAS_ELASTIC function completes here
# ============================================================================

	# Calculates the kernel for inelastic scattering
	# from Legendre expansion coefficients.
	
def SPKAINEL(alc1,NLa,ret,nbge,A,Z,En,Q,y):
	Pl1 = [0]*NLa
	val2 = [0]*nbge
	U = (A+1)*(-Q)/A
	R = numpy.sqrt(1-U/En)
	Emx = y*En
	for i in range (nbge):
		xabc = (1 - (ret[i]/Emx) +  R*R)/(2*R)
		if (-1.0 <= xabc and xabc <= 1.0):
			Pl1[0] = 1
			Pl1[1] = xabc
			for j in range (2,NLa):
				Pl1[j] = (((2*(j-1)+1)*xabc*Pl1[j-1])-((j-1)*Pl1[j-2]))/j
			p1 = 0
			for l in range (NLa):
				p1 = p1 + (((2*l)+1)*Pl1[l]*alc1[l])/2
			val2[i] = p1

	return(val2)

#=======Inelastic PKA Spectra=======*
	
	# Calcualtion of PKA spectra due to inelastic scattering 
	# interaction of neutron.
	
def PKAS_INELASTIC (insp,eliso,ret,nbge,nbpoints,nre,igtype):

	print ("n,n' .....")
	val2 = [0]*nbge; cr2 = [0]*nbge
	dsgmdT = numpy.zeros((nbge,nbge))
	dsgmdTl = numpy.zeros((41,nbge,nbge))
	sret = numpy.zeros((nbge,41))
	pkamatr = numpy.zeros((nre,nre)); pka2 = numpy.zeros((nre,nre))
	EL = numpy.zeros((500,41))
	alc = [0]*65; QI = [0]*41; mta = [0]*41; NP = [0]*41; NE1 = [0]*41
	NLarr = numpy.zeros((500,41))
	xabc = [0]*64; wg = [0]*64; iflpresent = [0]*41

	# for file4
	NE4 = [0]*41
	NL4 = numpy.zeros((41,1000))
	En4 = numpy.zeros((41,1000))
	al4 = numpy.zeros((41,1000,65))

	# for file5
	NK5 = [0]*41
	NP5 = numpy.zeros((41,5)); NE5 = numpy.zeros((41,5)); LF = numpy.zeros((41,5))
	NF5 = numpy.zeros((41,5,100))
	Eint = numpy.zeros((41,5,50)); p = numpy.zeros((41,5,50))
	En5 = numpy.zeros((41,5,100)); tht = numpy.zeros((41,5,100))
	Enp5 = numpy.zeros((5,100,500)); f5 = numpy.zeros((5,100,500))
	
	# for file6
	fiso = [0]*3; NK = [0]*41
	LAW = numpy.zeros((41,3)); NP6 = numpy.zeros((41,3)); LG = numpy.zeros((41,3))
	LEP = numpy.zeros((41,3)); NE6 = numpy.zeros((41,3))
	NEP = numpy.zeros((41,3,100)); NL = numpy.zeros((41,3,100))
	NBT = [0]*50; INTr = [0]*50
	Eint6 = numpy.zeros((41,3,100)); yi = numpy.zeros((41,3,100))
	En = numpy.zeros((41,3,100))
	fr6 = numpy.zeros((1000,nbge))
	Enp = numpy.zeros((3,100,1000)); f = numpy.zeros((3,100,1000))
	al6 = numpy.zeros((41,3,100,65))
	
	nrg = nre + 1
	
	ofile101 = open("Output_RadEMC-RecedU.txt", 'a')
	print('', file = ofile101)
	print(' Inelastic scattering (MT = 4) PKA spectra', file = ofile101)
	print('-------------------------------------------', file = ofile101)

	ifile = open ("tape01", 'r')
	ifile.readline()
	line = ifile201.readline()
	(ZA,AWR,L0,L1,L2,L3,MAT,MF,MT) = line_type1_info(line)
	ifile.close()

	Z = int(ZA/1000)
	A = AWR

	y = A/((A+1)**2)
	n1 = 1/(A+1)
	n2 = A/(A+1)

	print( 'Discrete Levels and Continuum .....')

	m = 51
	for i in range (41):
		mta[i] = m
		m = m + 1

	ifile = open ("tape02", 'r')
	l = 0
	while True:
		line = ifile.readline()
		data = eachlineinfo(line)
		MAT = int(data[6]); MF = int(data[7]); MT  = int(data[8])
		if (MAT != -1):
			if (MF == 3):
				if (MT == 1):
					line = ifile.readline()
					data = eachlineinfo(line)
					QMt = float(data[0]); QIt =  float(data[1])
					LR = int(data[3]); NR = int(data[4]); NPt = int(data[5])
					ifile.readline()
					Et = [0]*NPt; sigt = [0]*NPt; Etu = [0]*NPt; sigit = numpy.zeros((NPt,41))
					E = numpy.zeros((NPt,41)); sig = numpy.zeros((NPt,41))
					sigl = numpy.zeros((NPt,41)); alfull = numpy.zeros((NPt,65))
					(Et, sigt) = line_type3_info(ifile,NPt,2)

				if (MT == mta[l] or MT == 91):
					if (MT == 91):
						l = 40
					iflpresent[l] = 1
					line = ifile.readline()
					data = eachlineinfo(line)
					QM = float(data[0]); QI[l] =  float(data[1])
					LR = int(data[3]); NR = int(data[4]); NP[l] = int(data[5])
					ifile.readline()
					temporary1 = [0]*NP[l]
					temporary2 = [0]*NP[l]
					(temporary1,temporary2) = line_type3_info(ifile,NP[l],2)
					for i, (value1, value2) in enumerate(zip(temporary1, temporary2), 0):
						E[i][l] = value1
						sig[i][l] = value2
					if (l != 40):
						l = l + 1
		else:
			break
	ifile.close()

	if4d = 0
	if4c = 0
	ifile = open ("tape01", 'r')
	l = 0
	while True:
		line = ifile.readline()
		if (line == ''):
			break
		data = eachlineinfo(line)
		MAT = int(data[6]); MF = int(data[7]); MT  = int(data[8])
		if (MT == 0):
			line = ifile.readline()
			(ZA,AWR,L1,LTTv,L2,L3,MAT,MF,MT) = line_type1_info(line)
		if (MAT != -1):
			if (MF == 4):
				if (MT == mta[l] or MT == 91):
					if (MT < 91):
						if4d = 1
					if (MT == 91):
						l = 40
						if4c = 1
					LTT = LTTv
					line = ifile.readline()
					(ZA,AWR,LI,LCT,L2,L3,MAT,MF,MT) = line_type1_info(line)
					if (LTT == 1 and LI == 0):
						line = ifile.readline()
						(C1,C2,L1,L2,NR,NE4[l],MAT,MF,MT) = line_type2_info(line)
						(NBT, INTr) = line_type3_info(ifile,NR,2)
						for i in range (int(NE4[l])):
							line = ifile.readline()
							(c1,En4[l][i],LT,l2,NL4[l][i],l3,MAT,MF,MT) = line_type2_info(line)
							NL4[l][i] = NL4[l][i] + 1
							al4[l][i][0] = 1
							temporary = [0]*int(NL4[l][i])
							temporary = line_type3_info (ifile,int(NL4[l][i])-1,1)
							for j, value in enumerate (temporary, 1):
								al4[l][i][j] = value
					if (l != 40):
						l = l+1
		else:
			break
	ifile.close()

	if5c = 0
	ifile = open ("tape01", 'r')
	while True:
		line = ifile.readline()
		if (line == ''):
			break
		data = eachlineinfo(line)
		MAT = int(data[6]); MF = int(data[7]); MT  = int(data[8])
		if (MT == 0):
			line = ifile.readline()
			(ZA,AWR,L1,L2,NKv,L3,MAT,MF,MT) = line_type1_info(line)
		if (MAT != -1):
			if (MF == 5):
				if (MT == 91):
					if5c = 1
					l = 40
					NK5[l] = NKv
					for NSS in range (int(NK5[l])):
						line = ifile.readline()
						(C1,C2,L1,LF[l][NSS],NR,NP5[l][NSS],MAT,MF,MT) = line_type2_info(line)
						if (LF[l][NSS] == 1):
							line = ifile.readline()
							(C1,C2,L1,L2,NR,NE5[l][NSS],MAT,MF,MT) = line_type2_info(line)
							(NBT, INTr) = line_type3_info(ifile,NR,2)
							for i in range (int(NE5[l][NSS])):
								line = ifile.readline()
								(C1,En5[l][NSS][i],L1,L2,NR,NF5[l][NSS][i],MAT,MF,MT) = line_type2_info(line)
								(NBT, INTr) = line_type3_info(ifile,NR,2)
								temporary1 = [0]*int(NF5[l][NSS][i])
								temporary2 = [0]*int(NF5[l][NSS][i])
								(temporary1,temporary2) = line_type3_info(ifile,int(NF5[l][NSS][i]),2)
								for j, (value1, value2) in enumerate(zip(temporary1, temporary2), 0):
									Enp5[NSS][i][j] = value1
									f5[NSS][i][j] = value2
						if (LF[l][NSS] == 9):
							(NBT, INTr) = line_type3_info(ifile,NR,2)
							temporary1 = [0]*int(NP5[l][NSS])
							temporary2 = [0]*int(NP5[l][NSS])
							(temporary1,temporary2) = line_type3_info(ifile,int(NP5[l][NSS]),2)
							for j, (value1, value2) in enumerate(zip(temporary1, temporary2), 0):
								Eint[l][NSS][j] = value1
								p[l][NSS][j] = value2
							line = ifile.readline()
							(C1,C2,L1,L2,NR,NE5[l][NSS],MAT,MF,MT) = line_type2_info(line)
							(NBT, INTr) = line_type3_info(ifile,NR,2)
							temporary1 = [0]*int(NE5[l][NSS])
							temporary2 = [0]*int(NE5[l][NSS])
							(temporary1,temporary2) = line_type3_info(ifile,int(NE5[l][NSS]),2)
							for j, (value1, value2) in enumerate(zip(temporary1, temporary2), 0):
								En5[l][NSS][j] = value1
								tht[l][NSS][j] = value2
		else:
			break
	ifile.close()

	if6d = 0
	if6c = 0
	if (if4d == 0 or if4c == 0 or if5c == 0):
		ifile = open ("tape01", 'r')
		l = 0
		while True:
			line = ifile.readline()
			if (line == ''):
				break
			data = eachlineinfo(line)
			MAT = int(data[6]); MF = int(data[7]); MT  = int(data[8])
			if (MT == 0):	#.and. NS==99999
				line = ifile.readline()
				(ZAv,AWRv,l1,LCT,NKv,l2,MAT,MF,MT) = line_type1_info(line)
			if (MAT != -1):
				if (MF == 6):
					if (MT == mta[l] or MT == 91):
						if (MT < 91):
							if6d = 1
						if (MT == 91):
							l = 40
							if6c = 1
						NK[l] = NKv
						for NSS in range (int(NK[l])):
							line = ifile.readline()
							(ZAP,AWP,LIP,LAW[l][NSS],NR,NP6[l][NSS],MAT,MF,MT) = line_type1_info(line)
							(NBT, INTr) = line_type3_info(ifile,NR,2)
							temporary1 = [0]*int(NP6[l][NSS])
							temporary2 = [0]*int(NP6[l][NSS])
							(temporary1,temporary2) = line_type3_info(ifile,int(NP6[l][NSS]),2)
							for j, (value1, value2) in enumerate(zip(temporary1, temporary2), 0):
								Eint6[l][NSS][j] = value1
								yi[l][NSS][j] = value2

							if (LAW[l][NSS] == 2):
								line = ifile.readline()
								(c1,c2,l3,l4,NR,NE6[l][NSS],MAT,MF,MT) = line_type2_info(line)
								(NBT, INTr) = line_type3_info(ifile,NR,2)
								for i in range (int(NE6[l][NSS])):
									line = ifile.readline()
									(c1,En[l][NSS][i],LG[l][NSS],l2,NW,NL[l][NSS][i],MAT,MF,MT) = line_type2_info(line)
									NL[l][NSS][i] = NL[l][NSS][i] + 1
									al6[l][NSS][i][0] = 1
									temporary = [0]*int(NL[l][NSS][i])
									temporary = line_type3_info (ifile,int(NL[l][NSS][i])-1,1)
									for j, value in enumerate (temporary, 1):
										al6[l][NSS][i][j] = value 

							if (LAW[l][NSS] == 1):
								line = ifile.readline()
								(c1,c2,LG[l][NSS],LEP[l][NSS],NR,NE6[l][NSS],MAT,MF,MT) = line_type2_info(line)
								(NBT, INTr) = line_type3_info(ifile,NR,2)
								for i in range (int(NE6[l][NSS])):
									line = ifile.readline()
									(c1,En[l][NSS][i],ND,NA,NW,NEP[l][NSS][i],MAT,MF,MT) = line_type2_info(line)
									if (NA != 0):
										if (LG[l][NSS] == 2 or LG[l][NSS] == 1):
											fiso[NSS] = 1
											Ball = [0]*NW
											temporary = [0]*NW
											temporary = line_type3_info (ifile,NW,1)
											for j1, value in enumerate(temporary, 0):
												Ball[j1] = value
											K1 = 0
											K2 = 1
											for j in range (int(NEP[l][NSS][i])):
												Enp[NSS][i][j] = Ball[K1]
												f[NSS][i][j] = Ball[K2]
												K1 = K1+NA+2
												K2 = K2+NA+2

									if (LG[l][NSS] == 1 and NA == 0):
										fiso[NSS] = 1
										temporary1 = [0]*int(NEP[l][NSS][i])
										temporary2 = [0]*int(NEP[l][NSS][i])
										(temporary1,temporary2) = line_type3_info(ifile,int(NEP[l][NSS][i]),2)
										for j, (value1, value2) in enumerate(zip(temporary1, temporary2), 0):
											Enp[NSS][i][j] = value1
											f[NSS][i][j] = value2
						if (l != 40):
							l = l+1
			else:
				break
	ifile.close()

	(xabc, wg) = GQ()

	print('', file = ofile101)
	print(NPt,' Energy points in total cross section', file = ofile101)

	# Make unique common energy in MF=3, MT=1

	Etu = numpy.array(Et)
	Etu = numpy.unique(Etu)
	NPt = len(Etu)
	print('', file = ofile101)
	print(NPt,' Unique total energy points', file = ofile101)
	
	# Finding (n,n') cross section in the unique energy array
	
	for l in range (41):
		if (iflpresent[l] == 1):
			k = 0
			for i in range (NPt):
				for j in range (k, int(NP[l])):
					if (sig[j][l] != 0):
						if (Etu[i] == E[j][l]):
							sigit[i][l] = sig[j][l]
							if (E[j][l] == E[j+1][l]):
								sigit[i][l] = sig[j+1][l]
								k = j + 1
							else:
								k = j
							break
						if (Etu[i] > E[j][l] and Etu[i] < E[j+1][l]):
							sigit[i][l] = crstd(Etu[i],E[j][l],E[j+1][l],sig[j][l],sig[j+1][l])
							k = j
							break
		

	# Finding cross sections corresponding to the fine energy group array

	for l in range (41):
		if (iflpresent[l] == 1):
			for i in range (nbge):
				for j in range (NPt-1):
					if (sigit[j][l] != 0):
						if (ret[i] == Etu[j]):
							sret[i][l] = sigit[j][l]
							break
						if (ret[i] > Etu[j] and ret[i] <= Etu[j+1]):
							if (ret[i] == Etu[j+1]):
								sret[i][l] = sigit[j+1][l]
							else:
								sret[i][l] = crstd(ret[i],Etu[j],Etu[j+1],sigit[j][l],sigit[j+1][l])
							break


	for l in range (40):
		if (iflpresent[l] == 1):
			print('Discrete MT = ', mta[l])

			# TO GET THE AL(E) FOR THE FULL ENERGY RANGE
			if (if6d == 1):
				for i in range (nbge):
					alfull[i][0] = 1
					for k in range(1, 65):
						alfull[i][k] = 0

				for i in range (nbge):
					for j in range (int(NE6[l][0])):
						if (ret[i] ==  En[l][0][j]):
							for k in range(1, 65):
								alfull[i][k] = al6[l][0][j][k]

						if (ret[i]>En[l][0][j] and ret[i]<En[l][0][j+1]):
							diff1 = ret[i] - En[l][0][j]
							diff2 = ret[i] - ret[i-1]
							if (diff1 <= diff2):
								for k in range(1, 65):
									x = ret[i]
									x1 = En[l][0][j]
									x2 = En[l][0][j+1]
									y1 = al6[l][0][j][k]
									y2 = al6[l][0][j+1][k]
									alfull[i][k] = y1 + ((y2-y1)*(x-x1)/(x2-x1))
							if (diff2 < diff1):
								for k in range(1, 65):
									x = ret[i]
									x1 = ret[i-1]
									x2 = En[l][0][j+1]
									y1 = alfull[i-1][k]
									y2 = al6[l][0][j+1][k]
									alfull[i][k] = y1 + ((y2-y1)*(x-x1)/(x2-x1))

							for k in range(1, 65):
								if (alfull[i][k] == 0):
									alfull[i][k] = alfull[i-1][k]
							break

			if (if4d == 1):
				for i in range (nbge):
					alfull[i][0] = 1
					for k in range(1, 65):
						alfull[i][k] = 0

				if (LTT == 1 and LI == 0):
					for i in range (nbge):
						for j in range (int(NE4[l])):	
							if (ret[i] ==  En4[l][j]):
								for k in range (1, 65):
									alfull[i][k] = al4[l][j][k]
								break
							if (ret[i] > En4[l][j] and ret[i] < En4[l][j+1]):
								diff1 = ret[i] - En4[l][j]
								diff2 = ret[i] - ret[i-1]
								if (diff1 <= diff2):
									for k in range(1, 65):
										x = ret[i]
										x1 = En4[l][j]
										x2 = En4[l][j+1]
										y1 = al4[l][j][k]
										y2 = al4[l][j+1][k]
										alfull[i][k] = y1 + ((y2-y1)*(x-x1)/(x2-x1))
								if (diff2 < diff1):
									for k in range(1, 65):
										x = ret[i]
										x1 = ret[i-1]
										x2 = En4[l][j+1]
										y1 = alfull[i-1][k]
										y2 = al4[l][j+1][k]	
										alfull[i][k] = y1 + ((y2-y1)*(x-x1)/(x2-x1))

								for k in range(1, 65):
									if (alfull[i][k] == 0):
										alfull[i][k] = alfull[i-1][k]
								break 

 	# PKA SPECTRA FROM DISCRETE INELASTIC REACTION

			for i in range (nbge):
				Ex = ret[i]
				Q = QI[l]
				for k in range(1, 65):
					alc[k] = alfull[i][k]
				if (sret[i][l] != 0):
					val2 = SPKAINEL(alc,65,ret,nbge,A,Z,Ex,Q,y)
				for ia in range (nbge):
					dsgmdTl[l][i][ia] = abs(val2[ia])

			if(if4d == 1):
				print('Discrete angular distribution', file = ofile101)
				print('represented by Legendre coefficients in File 4', file = ofile101)
			if(if6d == 1):
				print('Discrete angular distribution', file = ofile101) 
				print('represented by Legendre coefficients in File 6', file = ofile101)

	l = 40

	if (iflpresent[l] == 1):
		print('Continuum MT = ', mta[l])

	for i in range (nbge):
		for j in range (nbge):
			dsgmdTl[l][i][j] = 0
		
	if (NK[l] == 3):
		Nrc = 1		# since arrays and lists are indexed from 0 - (n-1)
	if (NK[l] == 2):
		Nrc = 0
	if (NK[l] == 1):
		Nrc = 0

	if (if6c == 1 and Nrc != 0):
		print('', file = ofile101)
		print('Continuum reaction recoil energy distribution given in File 6', file = ofile101)

	# Interpolation for the recoil energies in finer array

		for i in range (int(NE6[l][Nrc])):
			for j in range (nbge):
				for k in range (int(NEP[l][Nrc][i])-1):
					if (Enp[Nrc][i][k] == ret[j]):
						fr6[i][j] = f[Nrc][i][k]
						break
					if (Enp[Nrc][i][k] < ret[j] and ret[j] < Enp[Nrc][i][k+1]):
						if (LEP[l][Nrc] == 1):
							fr6[i][j] = f[Nrc][i][k]
						if (LEP[l][Nrc] == 2):
							diff1 = ret[j] - Enp[Nrc][i][k]
							diff2 = ret[j] - ret[j-1]
							if (diff1 <= diff2):
								x = ret[j]
								x1 = Enp[Nrc][i][k]
								x2 = Enp[Nrc][i][k+1]
								y1 = f[Nrc][i][k]
								y2 = f[Nrc][i][k+1]
								fr6[i][j] = crstd(x,x1,x2,y1,y2)
							if (diff2 < diff1):
								x = ret[j]
								x1 = ret[j-1]
								x2 = Enp[Nrc][i][k+1]
								y1 = fr6[i][j-1]
								y2 = f[Nrc][i][k+1]
								fr6[i][j] = crstd(x,x1,x2,y1,y2)
						break

	# PKA SPECTRA FROM CONTINUUM INELASTIC REACTION
	# Interpolation for the neutron energies in finer array

		for i in range (nbge):
			if (sret[i][l] != 0):
				for j in range (int(NE6[l][Nrc])):
					if (ret[i] == En[l][Nrc][j]):
						for k in range (nbge):
							dsgmdTl[l][i][k] = fr6[j][k]
						break
					if (En[l][Nrc][j] < ret[i] and ret[i] < En[l][Nrc][j+1]):
						diff1 = ret[i] - En[l][Nrc][j]
						diff2 = ret[i] - ret[i-1]
						if (diff1 <= diff2):
							for k in range (nbge):
								x = ret[i]
								x1 = En[l][Nrc][j]	
								x2 = En[l][Nrc][j+1]
								y1 = fr6[j][k]
								y2 = fr6[j+1][k]
								dsgmdTl[l][i][k] = crstd(x,x1,x2,y1,y2)
						if (diff2 < diff1):
							for k in range (nbge):
								x = ret[i]
								x1 = ret[i-1]
								x2 = En[l][Nrc][j+1]
								y1 = dsgmdTl[l][i-1][k]
								y2 = fr6[j+1][k]
								dsgmdTl[l][i][k] = crstd(x,x1,x2,y1,y2)
						break

	# Finding PKA spectra using neutron evaporation model and 
	# isotropic emission of recoil nucleus assumption, for no recoil data
	
	Qs = abs(QI[0])
	Th = Qs/n2

	if (if6c == 0 or Nrc == 0):
		print('', file = ofile101)
		print('Continuum reaction recoil energy distribution', file = ofile101)
		print('not given in File 6, calculated with neutron', file = ofile101)
		print('evaporation model and isotropic emission of', file = ofile101) 
		print('recoil nucleus assumption', file = ofile101)
		for i in range (nbge):
			if (sret[i][l] != 0):
				f0 = sqrt(1 + Th/ret[i])
				f1 = 1 + Th/(2*ret[i])
				T1 = (0.5*4*y*ret[i])*(f1+f0)
				T2 = (0.5*4*y*ret[i])*(f1-f0)
				empmx = n2*(Qs+(n2*ret[i]))   #ret(i) - Th - T2
				empmn = 0				#(ret(i) - Th - T1)/1.8
				thn = sqrt(ret[i]/A)*3.22e+03
				for j in range (nbge):
					if (ret[j] <= T1):		#T2<=ret(j).and.
						sf = abs(Tinteg2 (empmx,empmn,thn,A,n1,n2))
						dsgmdTl[l][i][j] = sf

	MTtg = 51
	for l in range (41):
		for i in range (nbge):
			for j in range (nbge):
				dsgmdTl[l][i][j] = sret[i][l] * dsgmdTl[l][i][j]
				dsgmdT[i][j] = dsgmdTl[l][i][j]
		iflnz = prflag(nbge,dsgmdT)
		if (iflnz == 1):
			pkamatr = GROUP_INTEG (ret,nbge,nbpoints,nre,dsgmdT)
			pkamatr = PKAS_NORM_CROSSSEC (insp,MTtg,igtype,nrg,pkamatr)
			PRINTPKAS (eliso,MTtg,nre,pkamatr)
			for i in range (nre):
				for j in range (nre):
					pka2[i][j] = pka2[i][j] + pkamatr[i][j]
		MTtg = MTtg + 1

	PRINTPKAS (eliso,4,nre,pka2)

	ofile1000 = open ('ToAddAll.txt', 'a')
	print('', file = ofile1000)
	for i in range (nre):
		print (['{:.6E}'.format(pka2[i][j]) for j in range (nre)], file = ofile1000)


	ofile1000.close()
	ofile101.close()

	# PKAS_INELASTIC function completes here
# ==========================================================================

#======= Kernel for discrete (n,CPO) reactions ==========
	
	# Calculates the kernel of reaction for discrete level (n,p)
	# (n,d), etc. reactions from the Legendre expansion coefficients.
		
def SPKATH (alc1,NLa,ret,nbge,A,En,Q,bta):
	Pl1 = [0]*NLa
	val2 = [0]*nbge
	U = (A+1)*abs(Q)/A
	R = 0.0	

	# condition added for not calculatingin case of U < En to avoid negative values in sqrt

	if ((A*(A+1-bta)/bta) != 0 and bta != 0 and (1-U/En) != 0 and En != 0 and U < En):
		R = numpy.sqrt((A*(A+1-bta)/bta)*(1-U/En))
	R1 = bta*R/(A+1-bta)
	if (R1 != 0):		# condition added for not calculatingin case of U < En
		for i in range (nbge):
			xabc = (1+R1*R1-(ret[i]*(A+1)**2/(En*(A+1-bta))))/(2*R1)
			if (abs(xabc) <= 1.0):
				Pl1[0] = 1
				Pl1[1] = xabc
				for j in range (2, NLa):
					Pl1[j] = (((2*(j-1)+1)*xabc*Pl1[j-1])-((j-1)*Pl1[j-2]))/j
				p1 = 0
				for l in range (NLa):
					p1 = p1 + (((2*l)+1)*Pl1[l]*alc1[l])/2
				val2[i] = p1
	return(val2)
#=====================================================

#=======Remaining threshold CPO PKA Spectra=======* 
    
	# Calculation of PKA spectra due to the charged particle emission
	# reactions of neutron.
	
def PKAS_nCPO (MTi,lpr,ret,nbge):
	beta = [0]*281; n1 = [0]*281; ze = [0]*281; cbe = [0]*281
	xabc = [0]*4; wg = [0]*4
	ifl4 = [0]*281; ifl6 = [0]*281
	MTdthl = [0]*250
	alc = [0]*65
	sret = [0]*nbge; val2 = [0]*nbge
	dsgmdT = numpy.zeros((nbge,nbge))

	# for file6
	LAW = [0]*5; NP6 = [0]*5; LG = [0]*5; LEP = [0]*5; NE6 = [0]*5
	NEP = numpy.zeros((5,200)); NL = numpy.zeros((5,200))
	NBT = [0]*50; INTr = [0]*50
	Eint6 = numpy.zeros((5,300)); yi = numpy.zeros((5,300))
	En = numpy.zeros((5,400))
	Enp = numpy.zeros((5,400,1000)); f = numpy.zeros((5,400,1000))
	al6 = numpy.zeros((5,400,65))
	cdata = numpy.zeros((5,400,100)); fdata = numpy.zeros((5,400,100))
	#REAL,DIMENSION(1000,nbge):: fr6
	ZAP = [0]*5; AWP = [0]*5

	ifdpd = 0

	for i in range (250):
		MTdthl[i] = i+600

	# beta array stores the A' values of the emitted light charged particles
	# in a reaction, which is more massive.
	# The array positions correspond to the (n, CPO) reaction MTs

	beta[0]=2;beta[1]=4;beta[2]=4;beta[3]=4;beta[4]=4;beta[5]=1
	beta[6]=4;beta[7]=4;beta[8]=2;beta[9]=3;beta[10]=3;beta[11]=4
	beta[12]=4;beta[13]=1;beta[14]=1;beta[15]=1;beta[16]=4
	beta[17]=1;beta[18]=2;beta[19]=3;beta[20]=3;beta[21]=4
	beta[22]=4;beta[23]=4;beta[24]=1;beta[25]=4;beta[26]=4
	beta[27]=4;beta[28]=2;beta[29]=3;beta[30]=4

	# these are for the MTs represented with discrete + continuum 
	# cross sections
	
	for i in range(31,81):
		beta[i] = 1
	for i in range(81,131):
		beta[i] = 2
	for i in range(131,181):
		beta[i] = 3
	for i in range(181,231):
		beta[i] = 3
	for i in range(231,281):
		beta[i] = 4
	
	# ze array stores the atomic numbers of the emitted light charged
	# particles corresponding to those which are given in array beta
	
	ze[0]=1;ze[1]=2;ze[2]=2;ze[3]=2;ze[4]=2;ze[5]=1;ze[6]=2
	ze[7]=2;ze[8]=1;ze[9]=1;ze[10]=2;ze[11]=2;ze[12]=2;ze[13]=1
	ze[14]=1;ze[15]=1;ze[16]=2;ze[17]=1;ze[18]=1;ze[19]=1;
	ze[20]=2;ze[21]=2;ze[22]=2;ze[23]=2;ze[24]=1;ze[25]=2;
	ze[26]=2;ze[27]=2;ze[28]=1;ze[29]=1;ze[30]=2

	for i in range(31,81):
		ze[i] = 1
	for i in range(81,131):
		ze[i] = 1
	for i in range(131,181):
		ze[i] = 1
	for i in range(181,231):
		ze[i] = 2
	for i in range(231,281):
		ze[i] = 2

	# ARGUEMENT POINTS FOR THE GAUSS QUADRATURE
	xabc = (-0.86114,-0.33998,0.33998,0.86114)
	#WEIGHTS CORRESPONDING TO THE ARGUEMENT POINTS
	wg = (0.34785,0.65215,0.65215,0.34785)
	
	ofile101 = open("Output_RadEMC-RecedU.txt", 'a')
	print('', file = ofile101)
	print(' CPO reaction (MT = ',MTi,') PKA spectra', file = ofile101)
	print('-------------------------------------------', file = ofile101)
	
	iflpresent = 0
	NP = 0
	ifile = open ("tape02", 'r')

	# extraction of cross sections
	while True:
		line = ifile.readline()
		if (line == ''):
			break
		data = eachlineinfo(line)
		MAT = int(data[6]); MF = int(data[7]); MT  = int(data[8])
		if ( MT == 0 and MF != 0):	# .and. NS==99999
			line = ifile.readline()
			data = eachlineinfo(line)
			MAT = int(data[6]); MF = int(data[7]); MT  = int(data[8])
			if ( MT == 0 and MF == 0):
				line = ifile.readline()
			(ZAv,AWRv,L0,L1,NKv,L2,MAT,MF,MT) = line_type1_info(line)
		if (MAT != -1):
			if (MF == 3):
				if (MT == MTi):
					iflpresent = 1
					ifdthl = 0
					for i in range (250):
						if (MT == MTdthl[i]):
							ifdthl = 1
							break
					ZA = ZAv
					AWR = AWRv
					line = ifile.readline()
					data = eachlineinfo(line)
					QM = float(data[0]); QI =  float(data[1])
					LR = int(data[3]); NR = int(data[4]); NP = int(data[5])
					Eall = [0]*NP; sall = [0]*NP; fr6 = numpy.zeros((NP,nbge))
					if (ifdthl == 1):
						alfull = numpy.zeros((nbge,65))
					ifile.readline()
					(Eall, sall) = line_type3_info(ifile,NP,2)
		else:
			break
	ifile.close()

#only if MF3 cross sections are available then do the following

	if (iflpresent == 1):
		Z = int(ZA/1000)
		A = AWR

		ifspad4 = 0 	# flag for secondary particle angular data
		ifspad4al = 0 	# flag for secondary particle angular data in 'al' coefficients
		ifspad4muf = 0 	# flag for secondary particle angular data in 'mu,f' form

		ifile = open ("tape01", 'r')

		while True:
			line = ifile.readline()
			if (line == ''):
				break
			data = eachlineinfo(line)
			MAT = int(data[6]); MF = int(data[7]); MT  = int(data[8])
			if (MF == 3 and MT==0):
				line = ifile.readline()
				data = eachlineinfo(line)
				MAT = int(data[6]); MF = int(data[7]); MT  = int(data[8])
				line = ifile.readline()
				(ZA,AWR,l1,LTT,NK,l2,MAT,MF,MT) = line_type1_info(line)
			if (MAT != -1):
				if (MF == 4):
					if (MT == MTi):
						ifl4[lpr] = 1
						for i in range (250):
							if (MT == MTdthl[i]):
								ifspad4 = 1
								break
						ifile.readline()
						if (LTT == 3 or LTT == 1):
							ifspad4al = 1
							line = ifile.readline()
							data = eachlineinfo(line)
							NE1 = int(data[5])
						# Legendre polynomial coefficients
							ifile.readline()
							al4 = numpy.zeros((NE1,65)); EL = numpy.zeros((NE1))
							for i in range (NE1):
								line = ifile.readline()
								data = eachlineinfo(line)
								T = float(data[0]); EL[i] = float(data[1]); NL4 = int(data[4])
								NL4 = NL4 + 1
								al4[i][0] = 1
								temporary = [0]*NL4
								temporary = line_type3_info (ifile,NL4-1,1)
								for j, value in enumerate (temporary, 1):
									al4[i][j] = value

						if (LTT == 3 or LTT == 2):
							ifspad4muf = 1
							line = ifile.readline()
							data = eachlineinfo(line)
							NE2 = int(data[5])
							# Probability
							ifile.readline()
							Enf = [0]*NE2; NPr = [0]*NE2; cdata4 = numpy.zeros((NE2,201))
							fdata4 = numpy.zeros((NE2,201)); ftotal = numpy.zeros((NE2,64))
							fmuE = numpy.zeros((NP,64)); fpr = [0]*64
							for i in range (NE2):
								line = ifile.readline()
								data = eachlineinfo(line)
								T = float(data[0]); Enf[i] = float(data[1])
								line = ifile.readline()
								data = eachlineinfo(line)
								NPr[i] = int(data[0])
								temporary1 = [0]*NPr[i]
								temporary2 = [0]*NPr[i]
								(temporary1,temporary2) = line_type3_info(ifile,NPr[i],2)
								for j, (value1, value2) in enumerate(zip(temporary1, temporary2), 0):
									cdata4[i][j] = value1
									fdata4[i][j] = value2
			else:
				break
		ifile.close()

		iflr = 0  		#flag to know presence of recoil data in MF = 6
		ifspad6 = 0		#flag for secondary particle angular data
		iffdlcd = 0 	#flag for fdata linear with cdata (LAW=2)
		iflgfdlcd = 0 	#flag for log(fdata) linear with cdata (LAW=2)

		ifile = open ("tape01", 'r')
		while True:
			line = ifile.readline()
			if (line == ''):
				break
			data = eachlineinfo(line)
			MAT = int(data[6]); MF = int(data[7]); MT  = int(data[8])
			if (MT == 0):
				line = ifile.readline()
				(ZAv,AWRv,l1,LCT,NKv,l2,MAT,MF,MT) = line_type1_info(line)
			if (MAT != -1):
				if (MF == 6):
					if (MT == MTi):
						ifl6[lpr] = 1
						NK = NKv
						for NSS in range (NK):
							line = ifile.readline()
							(ZAP[NSS],AWP[NSS],LIP,LAW[NSS],NR6,NP6[NSS],MAT,MF,MT) = line_type1_info(line)
							if (AWP[NSS] > beta[lpr]):
								iflr = 1
								irs = NSS
							(NBT, INTr) = line_type3_info(ifile,NR6,2)
							temporary1 = [0]*NP6[NSS]
							temporary2 = [0]*NP6[NSS]
							(temporary1,temporary2) = line_type3_info(ifile,NP6[NSS],2)
							for N, (value1, value2) in enumerate(zip(temporary1, temporary2), 0):
								Eint6[NSS][N] = value1
								yi[NSS][N] = value2

							if (LAW[NSS] == 2):
								if(int(math.ceil(AWP[NSS])) == beta[lpr]):
									ifspad6 = 1
									isps = NSS
								line = ifile.readline()
								(c1,c2,l3,l4,NR6,NE6[NSS],MAT,MF,MT) = line_type2_info(line)
								(NBT, INTr) = line_type3_info(ifile,NR6,2)
								for i in range (int(NE6[NSS])):
									line = ifile.readline()
									(c1,En[NSS][i],LG[NSS],l2,NW,NL[NSS][i],MAT,MF,MT) = line_type2_info(line)
									if (LG[NSS] == 0):
										NL[NSS][i] = NL[NSS][i] + 1
										al6[NSS][i][0] = 1
										temporary = [0]*int(NL[NSS][i])
										temporary = line_type3_info (ifile,int(NL[NSS][i])-1,1)
										for j, value in enumerate (temporary, 1):
											al6[NSS][i][j] = value
									if (LG[NSS] > 0):
										if (LG[NSS] == 12):
											iffdlcd = 1
										if (LG[NSS] == 14):
											iflgfdlcd = 1
										temporary1 = [0]*int(NL[NSS][i])
										temporary2 = [0]*int(NL[NSS][i])
										(temporary1,temporary2) = line_type3_info(ifile,int(NL[NSS][i]),2)
										for j, (value1, value2) in enumerate(zip(temporary1, temporary2), 0):
											cdata[NSS][i][j] = value1
											fdata[NSS][i][j] = value2
							if (LAW[NSS] == 1):
								line = ifile.readline()
								(c1,c2,LG[NSS],LEP[NSS],NR6,NE6[NSS],MAT,MF,MT) = line_type2_info(line)
								(NBT, INTr) = line_type3_info(ifile,NR6,2)
								for i in range (int(NE6[NSS])):
									line = ifile.readline()
									(c1,En[NSS][i],ND,NA,NW,NEP[NSS][i],MAT,MF,MT) = line_type2_info(line)
									if (NA != 0):
										if (LG[NSS] == 1 or LG[NSS] == 2):
											Ball = [0]*NW
											temporary = [0]*NW
											temporary = line_type3_info (ifile,NW,1)
											for j1, value in enumerate(temporary, 0):
												Ball[j1] = value
											K1 = 0
											K2 = 1
											for j in range (int(NEP[NSS][i])):
												Enp[NSS][i][j] = Ball[K1]
												f[NSS][i][j] = Ball[K2]
												K1 = K1+NA+2
												K2 = K2+NA+2

									if (LG[NSS] == 1 and NA == 0):
										temporary1 = [0]*int(NEP[NSS][i])
										temporary2 = [0]*int(NEP[NSS][i])
										(temporary1,temporary2) = line_type3_info(ifile,int(NEP[NSS][i]),2)
										for j, (value1, value2) in enumerate(zip(temporary1, temporary2), 0):
											Enp[NSS][i][j] = value1
											f[NSS][i][j] = value2
			else:
				break
		ifile.close()

		n1[lpr] = (A+1-beta[lpr])/(A+1)
		cbe[lpr] = (1.029e+6*ze[lpr]*Z)/(beta[lpr]**(1/3)+A**(1/3))

		n2 = A/(A+1)
		n3 = 1/(A+1)

		if (LR == 0):
			Q = abs(QI)
			U = Q/n2

		for i in range (nbge):
			if (ret[i] >= Eall[0]):
				for j in range (NP-1):
					if (sall[j] > 0):
						if (ret[i] == Eall[j]):
							sret[i] = sall[j]
							break
						if (Eall[j] < ret[i] and ret[i] <= Eall[j+1]):
							if (ret[i] == Eall[j+1]):
								sret[i] = sall[j+1]
							else:
								x = ret[i]
								x1 = Eall[j]
								x2 = Eall[j+1]
								y1 = sall[j]
								y2 = sall[j+1]
								if((x2-x1) != 0):
									sret[i] = crstd(x,x1,x2,y1,y2)
							break

# (n,p), (n,d), (n,t), (n,3He), (n,a) discrete level scattering
# leaving the residual nucleus in different excited states can be 
# given in -->> 
# MT = 600..-649 (continuum) for (n,p)
# MT = 650..-699 (continuum) for (n,d)
# MT = 700..-749 (continuum) for (n,t)
# MT = 750..-799 (continuum) for (n,3He)
# MT = 800..-849 (continuum) for (n,a)
  
# The following part calculates the recoil data in these cases 
# only if the angular data for discrete scatterings are given 
# (LAW=2 of MF6). For continuum (i.e., 649, 699, etc.) recoil
# data is computed if continuum distributions are present.
# In other cases, the one particle recoil approximation and idotropic
# emission of recoil nucleus are assumed.

		if (ifspad4 == 1):
			print('', file = ofile101)
			print('Discrete level distribution data given in File 4', file = ofile101)
			for i in range (nbge):
				alfull[i][0] = 1
				for k in range (1, 65):
					alfull[i][k] = 0
			for i in range (nbge):
				if (sret[i] != 0):
					for j in range (NE1):
						if (ret[i] ==  EL[j]):
							for k in range (1, 65):
								alfull[i][k] = al4[j][k]
							break
						if (EL[j] < ret[i] and ret[i] < EL[j+1]):
							diff1 = ret[i] - EL[j]
							diff2 = ret[i] - ret[i-1]
							if (diff1 <= diff2):
								for k in range (1, 65):
									x = ret[i]
									x1 = EL[j]	
									x2 = EL[j+1]
									y1 = al4[j][k]
									y2 = al4[j+1][k]
									if((x2-x1) != 0):
										alfull[i][k] = crstd(x,x1,x2,y1,y2)
							if (diff2 < diff1):
								for k in range (1, 65):
									x = ret[i]
									x1 = ret[i-1]
									x2 = EL[j+1]	
									y1 = alfull[i-1][k]
									y2 = al4[j+1][k]
									if((x2-x1) != 0):
										alfull[i][k] = crstd(x,x1,x2,y1,y2)
							for k in range (1, 65):
								if (alfull[i][k] == 0):
									alfull[i][k] = alfull[i-1][k]
							break

			for i in range (nbge):
				if (sret[i] != 0):
					for k in range (65):
						alc[k] = alfull[i][k]
					val2 = SPKATH(alc,65,ret,nbge,A,ret[i],QI,beta[lpr])
					for ia in range (nbge):
						dsgmdT[i][ia] = abs(val2[ia])

		if (ifspad6 == 1):
			print('', file = ofile101)
			print('Discrete level distribution data given in File 6', file = ofile101)
	
			for i in range (nbge):
				alfull[i][0] = 1
				for k in range (1, 65):
					alfull[i][k] = 0
			for i in range (nbge):
				if (sret[i] != 0):
					for j in range (int(NE6[isps])):
						if (ret[i] ==  En[isps][j]):
							for k in range (65):
								alfull[i][k] = al6[isps][j][k]
							break
						if (En[isps][j] < ret[i] and ret[i] < En[isps][j+1]):
							diff1 = ret[i] - En[isps][j]
							diff2 = ret[i] - ret[i-1]
							if (diff1 <= diff2):
								for k in range (1, 65):
									x = ret[i]
									x1 = En[isps][j]	
									x2 = En[isps][j+1]
									y1 = al6[isps][j][k]
									y2 = al6[isps][j+1][k]
									if((x2-x1) != 0):
										alfull[i][k] = crstd(x,x1,x2,y1,y2)
							if (diff2 < diff1):
								for k in range (1, 65):
									x = ret[i]
									x1 = ret[i-1]
									x2 = En[isps][j+1]
									y1 = alfull[i-1][k]
									y2 = al6[isps][j+1][k]
									if((x2-x1) != 0):
										alfull[i][k] = crstd(x,x1,x2,y1,y2)
							for k in range (1, 65):
								if (alfull[i][k] == 0):
									alfull[i][k] = alfull[i-1][k]
							break
			for i in range (nbge):
				if (sret[i] != 0):
					for k in range (65):
						alc[k] = alfull[i][k]
					val2 = spkath(alc,65,ret,nbge,A,ret[i],QI,beta[lpr])
					for ia in range (nbge):
						dsgmdT[i][ia] = abs(val2[ia])

		if (iflr == 1 and ifspad6 == 0):
			print('', file = ofile101)
			print('Continuum reaction recoil data given in File 6', file = ofile101)
			for i in range (int(NE6[irs])):
				for j in range (nbge):
					for k in range (int(NEP[irs][i]-1)):
						if(Enp[irs][i][k] == ret[j]):
							fr6[i][j] = f[irs][i][k]
							break
						if(Enp[irs][i][k] < ret[j] and ret[j] <= Enp[irs][i][k+1]):
							if (LEP[irs] == 1):
								fr6[i][j] = f[irs][i][k]
							if (LEP[irs] == 2):
								diff1 = ret[j] - Enp[irs][i][k]		
								diff2 = ret[j] - ret[j-1]
								if (diff1 <= diff2):
									x = ret[j]
									x1 = Enp[irs][i][k]
									x2 = Enp[irs][i][k+1]
									y1 = f[irs][i][k]
									y2 = f[irs][i][k+1]
									if((x2-x1) != 0):
										fr6[i][j] = crstd(x,x1,x2,y1,y2)
								if (diff2 < diff1):
									x = ret[j]
									x1 = ret[j-1]
									x2 = Enp[irs][i][k+1]
									y1 = fr6[i][j-1]
									y2 = f[irs][i][k+1]
									if((x2-x1) != 0):
										fr6[i][j] = crstd(x,x1,x2,y1,y2)
							break			
			for i in range (nbge):
				if (sret[i] != 0):
					for j in range (int(NE6[irs])):
						if(ret[i] == En[irs][j]):
							for k in range (nbge):
								dsgmdT[i][k] = fr6[j][k]
							break
						if (En[irs][j] < ret[i] and ret[i] < En[irs][j+1]):
							diff1 = ret[i] - En[irs][j]
							diff2 = ret[i] - ret[i-1]
							if (diff1 <= diff2):
								for k in range (nbge):
									x = ret[i]
									x1 = En[irs][j]
									x2 = En[irs][j+1]
									y1 = fr6[j][k]
									y2 = fr6[j+1][k]
									if((x2-x1) != 0):
										dsgmdT[i][k] = crstd(x,x1,x2,y1,y2)
							if (diff2 < diff1):
								for k in range (nbge):
									x = ret[i]
									x1 = ret[i-1]
									x2 = En[irs][j+1]
									y1 = dsgmdT[i-1][k]
									y2 = fr6[j+1][k]
									if((x2-x1) != 0):
										dsgmdT[i][k] = crstd(x,x1,x2,y1,y2)
							break

		if (iflr == 0 and ifspad6 == 0 and ifspad4 == 0):
			print('', file = ofile101)
			print('No discrete/continuum reaction data for recoil', file = ofile101) 
			print('are given in File 4/6. Calculation performed', file = ofile101)
			print('assuming isotropic emissopn of recoil nucleus', file = ofile101)

			for i in range (nbge):
				if (sret[i] != 0):
					Estar = n1[lpr]*ret[i]
					availE = QI + (n2*ret[i])
					if (availE < cbe[lpr]):
						Ea = availE
					if (cbe[lpr] < availE):
						Ea  = cbe[lpr]
					f1 = Estar
					f2 = 2*sqrt(beta[lpr]*Estar*Ea)
					f3 = beta[lpr]*Ea
					tmax = n3*(f1+f2+f3)
					tmin = n3*(f1-f2+f3)
					for j in range (nbge):
						if (ret[j] >= tmin and ret[j] <= tmax):
							dsgmdT[i][j] = 1.0/(tmax-tmin)

		for i in range (nbge):
			for j in range (nbge):
				if(dsgmdT[i][j] != 0 and sret[i] != 0):
					dsgmdT[i][j] = dsgmdT[i][j]*sret[i]

	# call PKAS_NORM_CROSSSEC(insp,nbge,ret,sret,dsgmdT)
		
		if (600 <= MTi and MTi <= 849 and (ifspad6 == 1 or ifspad4 == 1)):
			ifdpd = 1
		
	# **** the above is done only if MF3 for that MT is present ****

	ofile101.close()

	return(dsgmdT,ifdpd,iflpresent)

 	# PKAS_nCPO function completes here
# =====================================================================================	
		
#=======Remaining threshold PKA Spectra=======* 

	# Calculation of PKA spectra due to the charged particle emission
	# reactions of neutron. The PKA spectra estimated is due to all 
	# reactions given in inexplicit way.
	
def PKAS_redtnMF6MT5 (insp,eliso,ret,nbge,nbpoints,nre,igtype):

	print("n,anything .....")
	erg = [0]*(nre+1)
	gsig = [0]*(nre+1)
	sret = [0]*nbge; sret1 = [0]*nbge
	dsgmdT = numpy.zeros((nbge,nbge))
	pkamatr = numpy.zeros((nre,nre)); pkaMF6MT5 = numpy.zeros((nre,nre))

	# for file6
	Ball = [0]*1000
	fr6 = numpy.zeros((1000,nbge))

	ofile101 = open("Output_RadEMC-RecedU.txt", 'a')
	print(' (n, anything) (heavy recoil nuclei) spectra', file = ofile101)
	print('-----------------------------------------------', file = ofile101)
# --------------------------------------------------------------------------
	iflpresent = 0
	NP = 0
	ifile = open ("tape02", 'r')

	# extraction of cross sections
	while True:
		line = ifile.readline()
		if (line == ''):
			break
		data = eachlineinfo(line)
		MAT = int(data[6]); MF = int(data[7]); MT  = int(data[8])
		if ( MT==0 and MF != 0):	# .and. NS==99999
			line = ifile.readline()
			(ZA,AWR,L0,L1,NKv,L2,MAT,MF,MT) = line_type1_info(line)
		if (MAT != -1):
			if (MF == 3):
				if (MT == 5):
					iflpresent = 1
					line = ifile.readline()
					data = eachlineinfo(line)
					QM = float(data[0]); QI =  float(data[1])
					LR = int(data[3]); NR = int(data[4]); NP = int(data[5])
					Eall = [0]*NP; sall = [0]*NP
					ifile.readline()
					(Eall, sall) = line_type3_info(ifile,NP,2)
		else:
			break
	ifile.close()

# only if MF3 cross sections are available then do the following

	if (iflpresent == 1):
		Z = int(ZA/1000)
		A = AWR

		ifile = open ("tape01", 'r')
		while True:
			line = ifile.readline()
			if (line == ''):
				break
			data = eachlineinfo(line)
			MAT = int(data[6]); MF = int(data[7]); MT  = int(data[8])
			if (MAT == -1):
				break
			if (MT == 0):
				if (MF == 6):
					line = ifile.readline()
					(ZAv,AWRv,l1,LCT,NKv,l2, MAT, MF, MT) = line_type1_info(line)
				if (MF < 6):
					ifile.readline()
					line = ifile.readline()
					(ZAv,AWRv,l1,LCT,NKv,l2,MAT,MF,MT) = line_type1_info(line)
			if (MAT != -1):
				if (MF == 6):
					if (MT == 5):
						NK = NKv
						LAW = [0]*NK; LG = [0]*NK; LEP = [0]*NK; Enp = numpy.zeros((NK,400,1000))
						f = numpy.zeros((NK,400,1000)); Eint6 = numpy.zeros((NK,400))
						Yi = numpy.zeros((NK,400)); NBT = numpy.zeros((NK,20)); INTr = numpy.zeros((NK,20))
						ZAP = [0]*NK; AWP = [0]*NK; Nyld = [0]*NK; NEP = numpy.zeros((NK,400))
						NL = numpy.zeros((NK,400)); NE6 = [0]*NK; En = numpy.zeros((NK,400))
						break
		ifile.close() 

		ifile = open ("tape01", 'r')
		while True:
			line = ifile.readline()
			if (line == ''):
				break
			data = eachlineinfo(line)
			MAT = int(data[6]); MF = int(data[7]); MT  = int(data[8])
			if (MAT == -1):
				break
			if (MT == 0):
				if (MF == 6):
					line = ifile.readline()
					(ZAv,AWRv,l1,LCT,NKv,l2, MAT, MF, MT) = line_type1_info(line)
				if (MF < 6):
					ifile.readline()
					line = ifile.readline()
					(ZAv,AWRv,l1,LCT,NKv,l2, MAT, MF, MT) = line_type1_info(line)
			if (MAT != -1):
				if (MF == 6):
					if (MT == 5):
						NK = NKv
						for NSS in range (NK):
							line = ifile.readline()
							(ZAP[NSS],AWP[NSS],LIP,LAW[NSS],NR6,Nyld[NSS],MAT,MF,MT) = line_type1_info(line)
							temporary1 = [0]*NR6
							temporary2 = [0]*NR6
							(temporary1,temporary2) = line_type3_info(ifile,NR6,2)
							for N, (value1, value2) in enumerate(zip(temporary1, temporary2), 0):
								NBT[NSS][N] = value1
								INTr[NSS][N] = value2
							temporary1 = [0]*Nyld[NSS]
							temporary2 = [0]*Nyld[NSS]
							(temporary1,temporary2) = line_type3_info(ifile,Nyld[NSS],2)
							for N, (value1, value2) in enumerate(zip(temporary1, temporary2), 0):
								Eint6[NSS][N] = value1
								Yi[NSS][N] = value2
							if (LAW[NSS] == 1):
								line = ifile.readline()
								(c1,c2,LG[NSS],LEP[NSS],NR6,NE6[NSS],MAT,MF,MT) = line_type2_info(line)
								(NBTp, INTrp) = line_type3_info(ifile,NR6,2)
								for i in range (int(NE6[NSS])):
									line = ifile.readline()
									(c1,En[NSS][i],ND,NA,NW,NEP[NSS][i],MAT,MF,MT) = line_type2_info(line)
									if (NA != 0):
										if (LG[NSS] == 1 or LG[NSS] == 2):
											temporary = [0]*NW
											temporary = line_type3_info (ifile,NW,1)
											for J1, value in enumerate(temporary, 0):
												Ball[J1] = value
											K1 = 0
											K2 = 1
											for j in range (int(NEP[NSS][i])):
												Enp[NSS][i][j] = Ball[K1]
												f[NSS][i][j] = Ball[K2]
												K1 = K1+NA+2
												K2 = K2+NA+2
									if (LG[NSS] == 1 and NA == 0):
										temporary1 = [0]*int(NEP[NSS][i])
										temporary2 = [0]*int(NEP[NSS][i])
										(temporary1,temporary2) = line_type3_info(ifile,int(NEP[NSS][i]),2)
										for j, (value1, value2) in enumerate(zip(temporary1, temporary2), 0):
											Enp[NSS][i][j] = value1
											f[NSS][i][j] = value2
			else:
				break
		ifile.close()

		for i in range (nbge):
			if (ret[i] >= Eall[0]):
				for j in range (NP-1):
					if (sall[j] > 0):
						if (ret[i] == Eall[j]):
							sret1[i] = sall[j]
							break
						if (Eall[j] < ret[i] and ret[i] <= Eall[j+1]):
							if (ret[i] == Eall[j+1]):
								sret1[i] = sall[j+1]
							else:
								x = ret[i]
								x1 = Eall[j]
								x2 = Eall[j+1]
								y1 = sall[j]
								y2 = sall[j+1]
								sret1[i] = crstd(x,x1,x2,y1,y2)
							break
		
		tnYldg = [0]*nbge; fr6 = numpy.zeros((max(NE6[:]),nbge))

		for NSS in range (NK):
			if (math.ceil(AWP[NSS]) > 4 and math.ceil(ZAP[NSS]//1000) > 2):
				itnZA = int(ZAP[NSS])
				print(itnZA)
				tnYldg = TERPOL(NBT[NSS][:],INTr[NSS][:],Nyld[NSS],Eint6[NSS][:],Yi[NSS][:],nbge,ret)

				for i in range(int(NE6[NSS])):
					for j in range(nbge):
						for k in range(int(NEP[NSS][i])-1):
							if (Enp[NSS][i][k] == ret[j]):
								fr6[i][j] = f[NSS][i][k]
								break
							if (Enp[NSS][i][k] < ret[j] and ret[j] < Enp[NSS][i][k+1]):
								if (LEP[NSS] == 1):
									fr6[i][j] = f[NSS][i][k]
								if (LEP[NSS] == 2):
									diff1 = ret[j] - Enp[NSS][i][k]
									diff2 = ret[j] - ret[j-1]
									if (diff1 <= diff2):
										x = ret[j]
										x1 = Enp[NSS][i][k]
										x2 = Enp[NSS][i][k+1]
										y1 = f[NSS][i][k]
										y2 = f[NSS][i][k+1]
										if((x2-x1) != 0):
											fr6[i][j] = crstd(x,x1,x2,y1,y2)
									if (diff2 < diff1):
										x = ret[j]
										x1 = ret[j-1]
										x2 = Enp[NSS][i][k+1]
										y1 = fr6[i][j-1]
										y2 = f[NSS][i][k+1]
										if((x2-x1) != 0):
											fr6[i][j] = crstd(x,x1,x2,y1,y2)
								break
				for i in range (nbge):
					if (sret1[i] != 0):
						for j in range (int(NE6[NSS])):
							if (ret[i] == En[NSS][j]):
								for k in range (nbge):
									dsgmdT[i][k] = fr6[j][k]
								break
							if (En[NSS][j] < ret[i] and ret[i] < En[NSS][j+1]):
								diff1 = ret[i] - En[NSS][j]
								diff2 = ret[i] - ret[i-1]
								if (diff1 <= diff2):
									for k in range (nbge):
										x = ret[i]
										x1 = En[NSS][j]
										x2 = En[NSS][j+1]
										y1 = fr6[j][k]
										y2 = fr6[j+1][k]
										if((x2-x1) != 0):
											dsgmdT[i][k] = crstd(x,x1,x2,y1,y2)
								if (diff2 < diff1):
									for k in range (nbge):
										x = ret[i]
										x1 = ret[i-1]
										x2 = En[NSS][j+1]
										y1 = dsgmdT[i-1][k]
										y2 = fr6[j+1][k]
										if((x2-x1) != 0):
											dsgmdT[i][k] = crstd(x,x1,x2,y1,y2)
								break
				for i in range (nbge):
					sret[i] = tnYldg[i]*sret1[i]
					for j in range (nbge):
						if(dsgmdT[i][j] != 0 and sret[i] != 0):
							dsgmdT[i][j] =sret[i] * dsgmdT[i][j]

				pkamatr = GROUP_INTEG (ret,nbge,nbpoints,nre,dsgmdT)
				print('normalize .....')
				if (igtype==1): 
					erg = engrp1()
				if (igtype==2):
					erg = engrp2()
				if (igtype==3):
					erg = engrp3()
				if (igtype==4):
					erg = engrp4()
				if (igtype==5):
					erg = engrp5()
				if (igtype==6):
					erg = engrp6()
				if (igtype==7):
					erg = engrp7()
				if (igtype==8):
					erg = engrp8()
				if (igtype==9):
					erg = engrp9()
				if (igtype==10):
					erg = engrp10()
				if (igtype==11):
					erg = engrp11()
				if (igtype==12):
					erg = engrp12()

				nrg = nre + 1
				gsig = TERPOLAPR (2,nbge,ret,sret,nrg,erg)
				for i in range (nre):
					if (gsig[i] != 0):
						sumnorm = 0
						for j in range (nre):
							sumnorm = sumnorm + pkamatr[i][j]*(erg[j+1]-erg[j])
						for j in range (nre):
							if (sumnorm != 0):
								pkamatr[i][j] = pkamatr[i][j]*gsig[i]/sumnorm

				PRINTPKAS (eliso,itnZA,nre,pkamatr)
        
				for i in range (nre):
					for j in range (nre):
						pkaMF6MT5[i][j] = pkaMF6MT5[i][j] + pkamatr[i][j] # sum of MF6 MT5 heavy recoil nuclei
				
		print('Total inexplicit')
		# print sum of all heavy recoil transmuted species from lumped data
		PRINTPKAS (eliso,6501,nre,pkaMF6MT5) 
	
			# for all the heavy recoil nuclei
		# for all the sections in MF=6 MT=5
	# **** the above is done only if MF3 for that MT is present ****

	ofile101.close()

# ----------------------------------------------------------------------------
	ofile1000 = open ('ToAddAll.txt', 'a')
	print('', file = ofile1000)
	for i in range (nre):
		print (['{:.6E}'.format(pkaMF6MT5[i][j]) for j in range (nre)], file = ofile1000)

	ofile1000.close()
# ----------------------------------------------------------------------------		

 	# REDTNMF6MT5 function completes here
# ==========================================================================
		
#=======(n,xn) reactions PKA spectra=======*

	# Calculation of PKA spectra due to (n,2n), (n,3n) and (n,4n)
	# reactions of neutron.
		
def PKAS_nxn (MTi,ret,nbge,nbpoints,nre,igtype):

	cr4 = [0]*nbge
	pka4 = numpy.zeros((nbge,nbge))

	#for file 6
	LAW = [0]*4; NP6 = [0]*4; LG = [0]*4; LEP = [0]*4; NE6 = [0]*4
	NEP = numpy.zeros((4,200)); NL = numpy.zeros((4,200))
	NBT = [0]*50; INTr = [0]*50
	Eint6 = numpy.zeros((4,50)); yi = numpy.zeros((4,50))
	En = numpy.zeros((4,200))
	Enp = numpy.zeros((4,200,5000)); f = numpy.zeros((4,200,5000))
	al6 = numpy.zeros((4,100,65))
	Ball = [0]*1000; fr6 = numpy.zeros((1000,5000))

	nrg = nre + 1

	ofile101 = open("Output_RadEMC-RecedU.txt", 'a')
	print('', file = ofile101)
	print( ' (n, xn) reaction (MT = ',MTi,') PKA spectra', file = ofile101)
	print('-----------------------------------------------', file = ofile101)

	ifile = open ("tape01", 'r')
	iflpresent = 0

	ifile.readline()
	line = ifile.readline()
	(ZA,AWR,l1,LTT,NK,l2, MAT, MF, MT) = line_type1_info(line)
	ifile.close()

	ifile = open ("tape02", 'r')
	while True:
		line = ifile.readline()
		if (line == ''):
			break
		data = eachlineinfo(line)
		MAT = int(data[6]); MF = int(data[7]); MT  = int(data[8])
		if (MAT != -1):
			if (MF == 3):
				if (MT == MTi):
					iflpresent = 1
					line = ifile.readline()
					data = eachlineinfo(line)
					QM = float(data[0]); QI =  float(data[1])
					LR = int(data[3]); NR = int(data[4]); NP = int(data[5])
					E = [0]*NP; sig = [0]*NP
					ifile.readline()
					(E, sig) = line_type3_info(ifile,NP,2)
		else:
			break
	ifile.close()

	if (iflpresent == 1):
		if6 = 0
		ifile = open ("tape01", 'r')

		while True:
			line = ifile.readline()
			if (line == ''):
				break
			data = eachlineinfo(line)
			MAT = int(data[6]); MF = int(data[7]); MT  = int(data[8])
			if (MT == 0):
				if (MF == 6):
					line = ifile.readline()
					(ZAv,AWRv,l1,LCT,NKv,l2, MAT, MF, MT) = line_type1_info(line)
				if (MF < 6):
					ifile.readline()
					line = ifile.readline()
					(ZAv,AWRv,l1,LCT,NKv,l2, MAT, MF, MT) = line_type1_info(line)
			if (MAT != -1):
				if (MF == 6):
					if (MT == MTi):
						ZA = ZAv
						AWR = AWRv
						if6 = 1
						NK = NKv
						for NSS in range (NK):
							line = ifile.readline()
							(ZAP,AWP,LIP,LAW[NSS],NR,NP6[NSS],MAT,MF,MT) = line_type1_info(line)
							temporary1 = [0]*NR
							temporary2 = [0]*NR
							(temporary1,temporary2) = line_type3_info(ifile,NR,2)
							for N, (value1, value2) in enumerate(zip(temporary1, temporary2), 0):
								NBT[N] = value1
								INTr[N] = value2
							temporary1 = [0]*NP6[NSS]
							temporary2 = [0]*NP6[NSS]
							(temporary1,temporary2) = line_type3_info(ifile,NP6[NSS],2)
							for N, (value1, value2) in enumerate(zip(temporary1, temporary2), 0):
								Eint6[NSS][N] = value1
								yi[NSS][N] = value2
							if (LAW[NSS] == 2):
								line = ifile.readline()
								(c1,c2,l3,l4,NR,NE6[NSS],MAT,MF,MT) = line_type2_info(line)
								temporary1 = [0]*NR
								temporary2 = [0]*NR
								(temporary1,temporary2) = line_type3_info(ifile,NR,2)
								for N, (value1, value2) in enumerate(zip(temporary1, temporary2), 0):
									NBT[N] = value1
									INTr[N] = value2
								for i in range (int(NE6[NSS])):
									line = ifile.readline()
									(c1,En[NSS][i],LG[NSS],l2,NW,NL[NSS][i], MAT,MF,MT) = line_type2_info(line)
									NL[NSS][i] = NL[NSS][i] + 1
									al6[NSS][i][0] = 1
									temporary = [0]*int(NL[NSS][i])
									temporary = line_type3_info (ifile,int(NL[NSS][i])-1,1)
									for j, value in enumerate (temporary, 1):
										al6[NSS][i][j] = value
							if (LAW[NSS] == 1):
								line = ifile.readline()
								(c1,c2,LG[NSS],LEP[NSS],NR,NE6[NSS],MAT,MF,MT) = line_type2_info(line)
								temporary1 = [0]*NR
								temporary2 = [0]*NR
								(temporary1,temporary2) = line_type3_info(ifile,NR,2)
								for N, (value1, value2) in enumerate(zip(temporary1, temporary2), 0):
									NBT[N] = value1
									INTr[N] = value2
									# 		read(601,*)
									# 		read(601,*)
								for i in range (int(NE6[NSS])):
									line = ifile.readline()
									(c1,En[NSS][i],ND,NA,NW,NEP[NSS][i],MAT,MF,MT) = line_type2_info(line)
									if (NA != 0):
										if (LG[NSS] == 1 or LG[NSS] == 2):
											Ball = [0]*NW
											temporary = [0]*NW
											temporary = line_type3_info (ifile,NW,1)
											for j1, value in enumerate(temporary, 0):
												Ball[j1] = value
											K1 = 0
											K2 = 1
											for j in range (int(NEP[NSS][i])):
												Enp[NSS][i][j] = Ball[K1]
												f[NSS][i][j] = Ball[K2]
												K1 = K1+NA+2
												K2 = K2+NA+2
									if (LG[NSS] == 1 and NA == 0):
										temporary1 = [0]*int(NEP[NSS][i])
										temporary2 = [0]*int(NEP[NSS][i])
										(temporary1,temporary2) = line_type3_info(ifile,int(NEP[NSS][i]),2)
										for j, (value1, value2) in enumerate(zip(temporary1, temporary2), 0):
											Enp[NSS][i][j] = value1
											f[NSS][i][j] = value2
			else:
				break
		ifile.close()

		for i in range (nbge):
			for j in range (NP-1):
				if (sig[j] > 0):
					if (ret[i] == E[j]):
						cr4[i] = sig[j]
						break
					if (ret[i] > E[j] and ret[i] <= E[j+1]):
						if (ret[i] == E[j+1]):
							cr4[i] = sig[j+1]
						else:
							x = ret[i]
							x1 = E[j]
							x2 = E[j+1]
							y1 = sig[j]
							y2 = sig[j+1]
							if((x2-x1) != 0):
								cr4[i] = crstd(x,x1,x2,y1,y2)
						break

		A = AWR
		Z = int(ZA/1000)
		y = 4*A/((A+1)**2)

		A1 = A-1
		n1 = 1/(A+1)
		n2 = A/(A+1)
		if (LR == 0):
			QM = abs(QI)
			Th = QM / n2

		if (NK == 3):
			Nrc = 1		# since arrays and lists indexed from 0 - (n-1)
		if (NK == 2):
			Nrc = 0
		if (NK == 1):
			Nrc = 0

		if (if6 == 1 and Nrc != 0):
			print('', file = ofile101)
			print('Distribution of recoil nucleus is given in File 6', file = ofile101)

			for i in range (int(NE6[Nrc])):
				for j in range (nbge):
					for k in range (int(NEP[Nrc][i]-1)):
						if(Enp[Nrc][i][k] == ret[j]):
							fr6[i][j] = f[Nrc][i][k]
							break
						if(Enp[Nrc][i][k] < ret[j] and ret[j] <= Enp[Nrc][i][k+1]):
							if (LEP[Nrc] == 1):
								fr6[i][j] = f[Nrc][i][k]
							if (LEP[Nrc] == 2):
								diff1 = ret[j] - Enp[Nrc][i][k]
								diff2 = ret[j] - ret[j-1]
								if (diff1 <= diff2):
									x = ret[j]
									x1 = Enp[Nrc][i][k]
									x2 = Enp[Nrc][i][k+1]
									y1 = f[Nrc][i][k]
									y2 = f[Nrc][i][k+1]
									if((x2-x1) != 0):
										fr6[i][j] = crstd(x,x1,x2,y1,y2)
								if (diff2 < diff1):
									x = ret[j]
									x1 = ret[j-1]
									x2 = Enp[Nrc][i][k+1]	
									y1 = fr6[i][j-1]
									y2 = f[Nrc][i][k+1]	
									if((x2-x1) != 0):
										fr6[i][j] = crstd(x,x1,x2,y1,y2)
							break
			for i in range (nbge):
				if (cr4[i] != 0):
					for j in range (int(NE6[Nrc])):
						if(ret[i] == En[Nrc][j]):
							for k in range (nbge):
								pka4[i][k] = fr6[j][k]
							break
						if(En[Nrc][j] < ret[i] and ret[i] < En[Nrc][j+1]):
							diff1 = ret[i] - En[Nrc][j]
							diff2 = ret[i] - ret[i-1]
							if (diff1 <= diff2):
								for k in range (nbge):
									x = ret[i]
									x1 = En[Nrc][j]	
									x2 = En[Nrc][j+1]
									y1 = fr6[j][k]
									y2 = fr6[j+1][k]
									if((x2-x1) != 0):
										pka4[i][k] = crstd(x,x1,x2,y1,y2)
							if (diff2 < diff1):
								for k in range (nbge):
									x = ret[i]
									x1 = ret[i-1]
									x2 = En[Nrc][j+1]	
									y1 = pka4[i-1][k]
									y2 = fr6[j+1][k]
									if((x2-x1) != 0):
										pka4[i][k] = crstd(x,x1,x2,y1,y2)
							break

	# In case of no recoil nucleus energy distribution data, 
	# evaporation model for emitted neutrons and isotropic recoil
	# nucleus emission is assumed. This is done only in case of 
	# (n,2n) reaction.

		if (MTi == 16 and (if6 == 0 or Nrc == 0)):
			print('', file = ofile101)
			print('Distribution of recoil nucleus is not given in', file = ofile101)
			print('File 6. Calculations are performed for (n, 2n)', file = ofile101)
			print('reaction assuming evaporation model and', file = ofile101)
			print('isotropic emission of recoil nucleus', file = ofile101)
	
			for i in range (nbge):
				if (cr4[i] != 0):
					Enn = ret[i]
					tht = sqrt(ret[i]/A)*3.22e+03
					Em1max = ret[i] - Th
					ul = Em1max
					Tlavg = (n1*n2*Enn) + (n1*ul/n2)
					T1=n1*A*ul/(n2*(A-1))+(A-1)*Tlavg/A+2*sqrt(n1*Tlavg*ul/n2)
					T2=n1*A*ul/(n2*(A-1))+(A-1)*Tlavg/A-2*sqrt(n1*Tlavg*ul/n2)
					emn = 0		#Em1max - T1)/2
					#Em1max = Em1max - T2
					Tx = Tinteg3(A,n1,n2,tht,Em1max,emn)
					for j in range (nbge):
						Er = ret[j]
						if (T2 <= Er and Er <= T1):
							pka4[i][j] = Tx
		
		for i in range (nbge):
			for j in range (nbge):
				pka4[i][j] = cr4[i]*pka4[i][j]
	
	# for the iflpresent if statement
	
	return(pka4, iflpresent)

	ofile101.close()

# PKAS_nxn function completes here
# ==========================================================================


#=======Kernel by assuming evaporation model=======*

	# The evaporation model for the emission of neutrons is assumed
	# and the kernel is computed only in case of (n,2n) reaction,
	# not in case of (n,3n) and (n,4n) reactions.
		
def Tinteg2 (Em2max,emn,tht,A,n1,n2):
	IEEm = tht**2*(1-((1+Em2max/tht)*math.exp(-Em2max/tht)))
	ll = emn
	ul = Em2max
	Nsect = 100
	h = (ul-ll)/Nsect
	s = 0
	for j in range (Nsect):
		r = ll + (j-1)*h
		f1 = r*math.exp(-r/tht) 
		f2 = (r+h)*math.exp(-(r+h)/tht)
		s = s + (h/2) * (f1 + f2)
	return(s/IEEm)

def Tinteg3 (A,n1,n2,tht,Em1max,emn):
	IE = tht**2 * (1 - ((1 + Em1max/tht)*math.exp(-Em1max/tht)))
	ll = emn
	ul = Em1max
	Nsect = 100
	h = (ul-ll)/Nsect
	s = 0
	for j in range (Nsect):
		r = ll + (j-1)*h
		Em2max1 = Em1max - r
		f1 = Tinteg2(Em2max1,0.0,tht,A,n1,n2)
		f11 = r*math.exp(-r/tht) * f1
		Em2max2 = Em1max - (r+h)
		f2 = Tinteg2(Em2max2,0.0,tht,A,n1,n2)
		f22 = (r+h)*math.exp(-(r+h)/tht) * f2
		s = s + (h/2) * (f11 + f22)
	return(s/IE)

#=======(n,g) reaction PKA spectra=======*

	# Calculation of PKA spectra due to (n,g) reaction of neutron.

def PKAS_ng (insp,eliso,ret,nbge,nbpoints,nre,igtype):

	print("n,g .....")
	ntm = [0]*nbge; ntmd = [0]*nbge; ncontd = [0]*nbge
	dsgmdT = numpy.zeros((nbge,nbge))
	sret = [0]*nbge
	pkamatr = numpy.zeros((nre,nre))
	eu = numpy.zeros((nbge,1000)); gu = numpy.zeros((nbge,1000)); eudisc = numpy.zeros((nbge,1000))
	gudisc = numpy.zeros((nbge,1000))
	
	# file 12
	NBT = [0]*20; INTr = [0]*20
	# File 12 -- only the LO = 1 (MULTIPLICITIES) option is coded
	# File 12 -- LO = 2 (TRANSITION PROBABILITY ARRAYS) option
	# is not coded.
	
	#file 6
	NPg6 = [0]*2000; ND6 = [0]*2000; NBT6 = [0]*20; INTr6 = [0]*20
	Eg6 = numpy.zeros((2000,1000)); g6 = numpy.zeros((2000,1000))
	E6 = [0]*2000; Y6 = [0]*2000; En6 = [0]*2000

	#file 15
	NBT15a = [0]*20; INTr15a = [0]*20; NPg15 = [0]*200
	Eg15 = numpy.zeros((200,200)); g15 = numpy.zeros((200,200)); E15 = [0]*100; Y15 = [0]*100; En15 = [0]*200
#------------------------------------------------------------
	nrg = nre + 1
	
	ofile101 = open("Output_RadEMC-RecedU.txt", 'a')
	print('', file = ofile101)
	print(' (n, g) reaction PKA spectra', file = ofile101)
	print('-----------------------------------------------', file = ofile101)

#------------------------------------------------------------
		
	ifl6 = 0
	ifile = open ("tape01", 'r')
	while True:
		line = ifile.readline()
		if (line == ''):
			break
		data = eachlineinfo(line)
		MAT = int(data[6]); MF = int(data[7]); MT  = int(data[8])
		if (MT == 0):
			line = ifile.readline()
			(ZAv,AWRv,L0,LCT,NKv,L2,MAT,MF,MT) = line_type1_info(line)
		if (MAT != -1):
			if (MF == 6):
				if (MT == 102):
					ifl6 = 1
					NK6 = NKv
					line = ifile.readline()
					(C1,C2,LIP,LAW,NR,NP6,MAT,MF,MT) = line_type2_info(line)
					temporary1 = [0]*NR
					temporary2 = [0]*NR
					(temporary1,temporary2) = line_type3_info(ifile,NR,2)
					for N, (value1, value2) in enumerate(zip(temporary1, temporary2), 0):
						NBT[N] = value1
						INTr[N] = value2
					temporary1 = [0]*NP6
					temporary2 = [0]*NP6
					(temporary1,temporary2) = line_type3_info(ifile,NP6,2)
					for N, (value1, value2) in enumerate(zip(temporary1, temporary2), 0):
						E6[N] = value1
						Y6[N] = value2
					Y6tot = [0]*nbge
					# Interpolating to find the total yields corresponding to energy points ret
					for i in range (nbge):
						for j in range (NP6):
							if (ret[i] == E6[j]):
								Y6tot[i] = Y6[j]
								break
							if (E6[j] < ret[i] and ret[i] < E6[j+1]):
								for k in range (NR):
									if (j <= NBT[k]):
										intflg = INTr[k]
										break
								Y6tot[i] = TERPOLIN(intflg,ret[i],E6[j],E6[j+1],Y6[j],Y6[j+1])
								break
					line = ifile.readline()
					(C1,C2,LANG,LEP,NR6,NE6,MAT,MF,MT) = line_type2_info(line)
					temporary1 = [0]*NR6
					temporary2 = [0]*NR6
					(temporary1,temporary2) = line_type3_info(ifile,NR6,2)
					for N, (value1, value2) in enumerate(zip(temporary1, temporary2), 0):
						NBT6[N] = value1
						INTr6[N] = value2
					for i in range (NE6):
						line = ifile.readline()
						(C1,En6[i],ND6[i],NA,NW,NPg6[i],MAT,MF,MT) = line_type2_info(line)
						temporary1 = [0]*NPg6[i]
						temporary2 = [0]*NPg6[i]
						(temporary1,temporary2) = line_type3_info(ifile,NPg6[i],2)
						for N, (value1, value2) in enumerate(zip(temporary1, temporary2), 0):
							Eg6[i][N] = value1
							g6[i][N] = value2
		else:
			break
	ifile.close()

	ifl12 = 0
	ifile = open ("tape01", 'r')
	while True:
		line = ifile.readline()
		if (line == ''):
			break
		data = eachlineinfo(line)
		MAT = int(data[6]); MF = int(data[7]); MT  = int(data[8])

		if (MF <= 12 and MT == 0):
			if (MAT != -1):
				line = ifile.readline()
				(ZAv,AWRv,Lo,L1,NKv,L2,MAT,MF,MT) = line_type1_info(line)
		if (MF == 0):
			if (MAT != -1):
				line = ifile.readline()
				(ZAv,AWRv,Lo,L1,NKv,L2,MAT,MF,MT) = line_type1_info(line)
		if (MAT != -1):
			if (MF == 12):
				if (MT == 102):
					ifl12 = 1
					NK = NKv
					AWR = AWRv
					ZA = ZAv
					NPn12 = [0]*NK; En12 = numpy.zeros((NK,200)); Yk12 = numpy.zeros((NK,200))
					EGk = [0]*NK; ESk = [0]*NK; LP = [0]*NK; NBTk12 = numpy.zeros((NK,20)); INTrk12 = numpy.zeros((NK,20))
					NRk12 = [0]*NK
					line = ifile.readline()
					(C1,C2,L1,L2,NR,NP12,MAT,MF,MT) = line_type2_info(line)
					E12 = [0]*NP12; Y12 = [0]*NP12
					temporary1 = [0]*NR
					temporary2 = [0]*NR
					(temporary1,temporary2) = line_type3_info(ifile,NR,2)
					for N, (value1, value2) in enumerate(zip(temporary1, temporary2), 0):
						NBT[N] = value1
						INTr[N] = value2
					temporary1 = [0]*NP12
					temporary2 = [0]*NP12
					(temporary1,temporary2) = line_type3_info(ifile,NP12,2)
					for N, (value1, value2) in enumerate(zip(temporary1, temporary2), 0):
						E12[N] = value1
						Y12[N] = value2
					Y12tot = [0]*nbge

					# Interpolating to find the total yields corresponding to energy
					# points in ret energy array

					for i in range (nbge):
						for j in range (NP12):
							if (ret[i] == E12[j]):
								Y12tot[i] = Y12[j]
								break
							if (E12[j] < ret[i] and ret[i] < E12[j+1]):
								for k in range (NR):
									if (j <= NBT[k]):
										intflg = INTr[k]
										break
								Y12tot[i] = TERPOLIN(intflg,ret[i],E12[j],E12[j+1],Y12[j],Y12[j+1])
								break

					if (NK == 1):
						for N in range (NP12):
							Yk12[NK][N] = Y12[N]
					if (NK > 1):
						for i in range (NK):
							line = ifile501.readline()
							(EGk[i],ESk[i],LP[i],LF,NRk12[i],NPn12[i],MAT,MF,MT) = line_type2_info(line)
							temporary1 = [0]*NRk12[i]
							temporary2 = [0]*NRk12[i]
							(temporary1,temporary2) = line_type3_info(ifile,NRk12[i],2)
							for N, (value1, value2) in enumerate(zip(temporary1, temporary2), 0):
								NBTk12[i][N] = value1
								INTrk12[i][N] = value2

							temporary1 = [0]*NPn12[i]
							temporary2 = [0]*NPn12[i]
							(temporary1,temporary2) = line_type3_info(ifile,NPn12[i],2)
							for N, (value1, value2) in enumerate(zip(temporary1, temporary2), 0):
								En12[i][N] = value1
								Yk12[i][N] = value2

					# Yield of discrete and continuum g for all neutron energy 
					# (ret energy array) by interpolation

					if (ifl12 == 1):
						Yk12tot = numpy.zeros((NK,nbge))
						if (NK == 1):
							for j in range (nbge):
								Yk12tot[NK][j] = Y12tot[j]
						if (NK > 1):
							for i in range (NK):
								for j in range (nbge):
									for k in range (int(NPn12[i])):
										if (ret[j] == En12[i][k]):
											Yk12tot[i][j] = Yk12[i][k]
											break
										if (En12[i][k] < ret[j] and ret[j] < En12[i][k+1]):
											for l in range (int(NRk12[i])):
												if (k <= NBTk12[i][l]):
													intflg = INTrk12[i][l]
													break
											Yk12tot[i][j] = TERPOLIN(intflg,ret[j],En12[i][k],En12[i][k+1],Yk12[i][k],Yk12[i][k+1])
											break
		else:
			break
	ifile.close()

	ifl15 = 0
	ifile = open ("tape01",'r')
	if (ifl6 == 0):
		while True:
			line = ifile.readline()
			if (line == ''):
				break
			data = eachlineinfo(line)
			MAT = int(data[6]); MF = int(data[7]); MT  = int(data[8])
			if (MT == 0):
				line = ifile.readline()
				(ZAv,AWRv,Lo,L1,NCv,L2,MAT,MF,MT) = line_type1_info(line)
			if (MAT != -1):
				if (MF == 15):
					if (MT == 102):
						ifl15 = 1
						NC = NCv
						line = ifile.readline()
						(C1,C2,L1,LF,NR,NP15,MAT,MF,MT) = line_type2_info(line)
						temporary1 = [0]*NR
						temporary2 = [0]*NR
						(temporary1,temporary2) = line_type3_info(ifile,NR,2)
						for N, (value1, value2) in enumerate(zip(temporary1, temporary2), 0):
							NBT[N] = value1
							INTr[N] = value2
						temporary1 = [0]*NP15
						temporary2 = [0]*NP15
						(temporary1,temporary2) = line_type3_info(ifile,NP15,2)
						for N, (value1, value2) in enumerate(zip(temporary1, temporary2), 0):
							E15[N] = value1
							Y15[N] = value2
						line = ifile.readline()
						(C1,C2,L1,L2,NR15a,NE15,MAT,MF,MT) = line_type2_info(line)
						temporary1 = [0]*NR15a
						temporary2 = [0]*NR15a
						(temporary1,temporary2) = line_type3_info(ifile,NR15a,2)
						for N, (value1, value2) in enumerate(zip(temporary1, temporary2), 0):
							NBT15a[N] = value1
							INTr15a[N] = value2
						for i in range (NE15):
							line = ifile.readline()
							(C1,En15[i],L1,L2,NR,NPg15[i],MAT,MF,MT) = line_type2_info(line)
							temporary1 = [0]*NR
							temporary2 = [0]*NR
							(temporary1,temporary2) = line_type3_info(ifile,NR,2)
							for N, (value1, value2) in enumerate(zip(temporary1, temporary2), 0):
								NBT[N] = value1
								INTr[N] = value2
							temporary1 = [0]*NPg15
							temporary2 = [0]*NPg15
							(temporary1,temporary2) = line_type3_info(ifile,NPg15,2)
							for N, (value1, value2) in enumerate(zip(temporary1, temporary2), 0):
								Eg15[N] = value1
								g15[N] = value2
			else:
				break
	ifile.close()

#-----------------------------------------------------------------------------
	ifile = open ("tape02", 'r')
	ifile.readline()
	line = ifile.readline()
	(ZA,AWR,L0,L1,L2,L3,MAT,MF,MT) = line_type1_info(line)
	while True:
		line = ifile.readline()
		if (line == ''):
			break
		data = eachlineinfo(line)
		MAT = int(data[6]); MF = int(data[7]); MT  = int(data[8])
		if (MAT != -1):
			if (MF == 3):
				if (MT == 102):
					line = ifile.readline()
					data = eachlineinfo(line)
					QM = float(data[0]); QI =  float(data[1])
					LR = int(data[3]); NR = int(data[4]); NP = int(data[5])
					ifile.readline()
					E = [0]*NP; sig = [0]*NP
					temporary1 = [0]*NP
					temporary2 = [0]*NP
					(temporary1,temporary2) = line_type3_info(ifile,NP,2)
					for i, (value1, value2) in enumerate(zip(temporary1, temporary2), 0):
						E[i] = value1
						sig[i] = value2
		else:
			break
	ifile.close()

	print('', file = ofile101)
	print('Number of data given for (n, g) cross section is ',NP, file = ofile101)
#-------------------------------------------------------------
	Z = int(ZA/1000)
	A1 = AWR + 1
	A1f = 1/A1
	if (LR == 0):
		QM = abs(QI)

	emc2 = 939.512E+06 
	tm = emc2*A1
	rtm = 1/tm  

	(Etu, NPt) = UQCE()
	siget = trptuqce(E,sig,Etu)

	#allocate(Etu(NPt),siget(NPt))

	print('', file = ofile101)
	print(NPt,' Unique total energy points', file = ofile101)

	# Finding cross sections corresponding to the fine energy array

	for i in range (nbge):
		for j in range (NPt-1):
			if (ret[i] == Etu[j]):
				sret[i] = siget[j]
				break
			if (ret[i] > Etu[j] and ret[i] <= Etu[j+1]):
				if (ret[i] == Etu[j+1]):
					sret[i] = siget[j+1]
				else:
					x = ret[i]
					x1 = Etu[j]
					x2 = Etu[j+1]
					y1 = siget[j]
					y2 = siget[j+1]
					if((x2-x1) != 0):
						sret[i] = crstd(x,x1,x2,y1,y2)
				break

	if (ifl6 == 1):
		print('', file = ofile101)
		print('Emitted photon data are given in File 6', file = ofile101)
		print('', file = ofile101)
		print('Neutron energy/ discrete gamma/ continuum gamma', file = ofile101)
		for i in range (0,NE6+5,5):
			print('', file = ofile101)
			print(En6[i],' ',ND6[i],' ',NPg6[i]-ND6[i], file = ofile101)
	if (ifl12 == 1):
		print('', file = ofile101)
		print('Emitted photon data are given in File 12', file = ofile101)

	if (ifl15 == 1):
		print('', file = ofile101)
		print('Emitted photon data are given in File 15 ', file = ofile101)
		
	# ENERGY OF THE EMITTED DISCRETE PHOTON FROM FILE 12

	if ((ifl12 == 1 and NK > 1) or (ifl12 == 1 and NK == 1 and ifl15 == 0)):
		EGkp = numpy.zeros((NK,nbge))
		NKd = NK
		if (ifl15 == 1):
			NKd = NK-1
		for i in range (NKd):
			for j in range (nbge):
				if (LP[i] == 2):
					EGkp[i][j] = EGk[i] + (AWR*ret[j]/(AWR+1))
				else:
					EGkp[i][j] = EGk[i]
		
		print('', file = ofile101)
		print('Total number of emitted gamma sections', file = ofile101)
		print(NKd,' discrete and',NK-NKd,' continuum', file = ofile101)
		
	# For each neutron energy, Total yield over all NK contributions
	# must be normalized to Y12tot. Total yields are given only
	# when NK > 1
		
	if (ifl12 == 1 and NK > 1):
		for j in range (nbge):
			sY12 = 0
			for i in range (NK):
				sY12 = sY12 + Yk12tot[i][j]
			for i in range (NK):
				Yk12tot[i][j] = Yk12tot[i][j]*Y12tot[j]/sY12
		 		
 	# AVERAGE OF (Egamma)SQUARE FROM FILE 6 AND FILE 15


	if (ifl15 == 1):
		Eg15t = numpy.zeros((nbge,1000)); g15t = numpy.zeros((nbge,1000))
		for i in range (nbge):
			if (Yk12tot[NK][i] != 0):
				for j in range (NE15):
					if(ret[i] == En15[j]):
						for k in range (int(NPg15[j])):
							Eg15t[i][k] = Eg15[j][k]
							g15t[i][k] = g15[j][k]
						ntm[i] = NPg15[j]
						break
					if (En15[j] < ret[i] and ret[i] < En15[j+1]):
						for k1 in range (NR15a):
							if (j <= NBT15a[k1]):
								iplaw = INTr15a[k1]
								break
						#if (Etu(i)<1.0e+6) iplaw = 1
						diff1 = ret[i] - En15[j]		
						diff2 = ret[i] - ret[i-1]
						if (diff1 <= diff2):
							for k in range (int(NPg15[j])):
								x = ret[i]
								x1 = En15[j]
								x2 = En15[j+1]
								y1 = Eg15[j][k]
								y2 = Eg15[j+1][k]
								y11 = g15[j][k]
								y22 = g15[j+1][k]
								Eg15t[i][k] = TERPOLIN(iplaw,x,x1,x2,y1,y2)
								g15t[i][k] = TERPOLIN(iplaw,x,x1,x2,y11,y22)
						if (diff2 < diff1):
							for k in range (int(NPg15[j])):
								x = ret[i]
								x1 = ret[i-1]
								x2 = En15[j+1]
								y1 = Eg15t[i-1][k]
								y2 = Eg15[j+1][k]
								y11 = g15t[i-1][k]
								y22 = g15[j+1][k]
								Eg15t[i][k] = TERPOLIN(iplaw,x,x1,x2,y1,y2)
								g15t[i][k] = TERPOLIN(iplaw,x,x1,x2,y11,y22)
						ntm[i] = NPg15[j]
						break
		
		for i in range (nbge):
			if (Yk12tot[NK][i] != 0):
				for j in range (int(ntm[i])):
					eu[i][j] = Eg15t[i][j]*Eg15t[i][j]*rtm/2
					gu[i][j] = Yk12tot[NK][i] * g15t[i][j] 

	if (ifl6 == 1):
		Eg6t = numpy.zeros((nbge,1000)); g6t = numpy.zeros((nbge,1000))
		for i in range (nbge):
			for j in range (NE6):
				if(ret[i] == En6[j]):
					for k in range (int(NPg6[j])):
						Eg6t[i][k] = Eg6[j][k]
						g6t[i][k] = g6[j][k]
					ntm[i] = NPg6[j]
					ntmd[i] = ND6[j]
					break
				# histogram interpolations are used to validate 
				# recoil energy (e.g. in W182 from ENDF/B-VII.1)
				if (En6[j] < ret[i] and ret[i] < En6[j+1]):
					for k1 in range (NR6):
						if (j <= NBT6[k1]):
							iplaw = INTr6[k1]
							break
					#if (ret[i] < 1.0e+6):
					#	iplaw = 1
					diff1 = ret[i] - En6[j]	
					diff2 = ret[i] - ret[i-1]
					if (diff1 <= diff2):
						for k in range (int(NPg6[j])):
							x = ret[i]
							x1 = En6[j]
							x2 = En6[j+1]
							y1 = Eg6[j][k]
							y2 = Eg6[j+1][k]
							y11 = g6[j][k]
							y22 = g6[j+1][k]
							Eg6t[i][k] = TERPOLIN(iplaw,x,x1,x2,y1,y2)
							g6t[i][k] = TERPOLIN(iplaw,x,x1,x2,y11,y22)
					if (diff2 < diff1):
						for k in range (int(NPg6[j])):
							x = ret[i]
							x1 = ret[i-1]
							x2 = En6[j+1]
							y1 = Eg6t[i-1][k]
							y2 = Eg6[j+1][k]
							y11 = g6t[i-1][k]
							y22 = g6[j+1][k]
							Eg6t[i][k] = TERPOLIN(iplaw,x,x1,x2,y1,y2)
							g6t[i][k] = TERPOLIN(iplaw,x,x1,x2,y11,y22)
					ntm[i] = NPg6[j]
					ntmd[i] = ND6[j]
					break
		
		for i in range (nbge):
			ncontd[i] = ntm[i]-ntmd[i]
			ndisc = ntmd[i]
			for j in range (ndisc):
				eudisc[i][j] = Eg6t[i][j]*Eg6t[i][j]*rtm/2
				gudisc[i][j] = g6t[i][j]
			for j in range (ndisc+1, ntm[i]):
				eu[i][j] = Eg6t[i][j]*Eg6t[i][j]*rtm/2
				gu[i][j] = g6t[i][j]

	if ((ifl12 == 1 and NK > 1) or (ifl12 == 1 and NK == 1 and ifl15 == 0)):
		NKd = NK
		if (ifl15 == 1):
			NKd = NK-1
		for i in range (nbge):
			for j in range (NKd):
				eudisc[i][j] = EGkp[j][i]*EGkp[j][i]*rtm/2
				gudisc[i][j] = Yk12tot[j][i]

	# Complete set of recoil energies and their distributions from discrete and continuum
	# As obtained from file by inerpolating for nbge incident neutron energies and 
	# nbge recoil atom energies

	if (ifl12 == 1):
		for i in range (nbge):
			for j in range (nbge):
				for k in range (NKd):
					if (ret[j] == eudisc[i][k]):
						dsgmdT[i][j] = gudisc[i][k]
						break
					# given gamma energies are in decreasing order of magnitude
					if (eudisc[i][k] > ret[j] and ret[j] > eudisc[i][k+1]):
						x = ret[j]
						x1 = eudisc[i][k]
						x2 = eudisc[i][k+1]
						y1 = gudisc[i][k]
						y2 = gudisc[i][k+1]
						if((x2-x1) != 0):
							y = crstd(x,x1,x2,y1,y2)
						dsgmdT[i][j] = y
						break
	if (ifl15 == 1):
		for i in range (nbge):
			s15 = 0
			s15 = s15 + sum(gu[i][:])
			s15norm = Yk12tot[NK][i]/s15
			for j in range (nbge):
				for k in range (int(ntm[i])):
					if (ret[j] == eu[i][k]):
						dsgmdT[i][j] = dsgmdT[i][j] + gu[i][k]*s15norm
						break
					if (eu[i][k] < ret[j] and ret[j] < eu[i][k+1]):
						x = ret[j]
						x1 = eu[i][k]
						x2 = eu[i][k+1]
						y1 = gu[i][k]
						y2 = gu[i][k+1]
						y = TERPOLIN(1,x,x1,x2,y1,y2)
						dsgmdT[i][j] = dsgmdT[i][j] + y*s15norm
						break
	if (ifl6 == 1):
		for i in range (nbge):
			s6 = 0
			s6 = s6 + sum(gudisc[i][:]) + sum(gu[i][:])
			for j in range (nbge):
				for k in range (int(ntmd[i])):
					if (ret[j] == eudisc[i][k]):
						dsgmdT[i][j] = gudisc[i][k]
						break
					# given gamma energies are in decreasing order of magnitude
					if (eudisc[i][k] > ret[j] and ret[j] > eudisc[i][k+1]):
						x = ret[j]
						x1 = eudisc[i][k]
						x2 = eudisc[i][k+1]
						y1 = gudisc[i][k]
						y2 = gudisc[i][k+1]
						y = TERPOLIN(LEP,x,x1,x2,y1,y2)
						dsgmdT[i][j] = y
						break
				for k in range (int(ncontd[i])):
					if (ret[j] == eu[i][k]):
						dsgmdT[i][j] = dsgmdT[i][j] + gu[i][k]
						break
					if (eu[i][k] < ret[j] and ret[j] < eu[i][k+1]):
						x = ret[j]
						x1 = eu[i][k]
						x2 = eu[i][k+1]
						y1 = gu[i][k]
						y2 = gu[i][k+1]
						y = TERPOLIN(LEP,x,x1,x2,y1,y2)
						dsgmdT[i][j] = dsgmdT[i][j] + y
						break
				if(s6 == 0):
					dsgmdT[i][j] = dsgmdT[i][j] * Y6tot[i]
				if(s6 != 0):
					dsgmdT[i][j] = dsgmdT[i][j] * Y6tot[i] / s6

	for i in range (nbge):
		for j in range (nbge):
			dsgmdT[i][j] = sret[i]*dsgmdT[i][j]

	pkamatr = GROUP_INTEG(ret,nbge,nbpoints,nre,dsgmdT)
	pkamatr = PKAS_NORM_CROSSSEC (insp,102,igtype,nrg,pkamatr)
	PRINTPKAS (eliso,102,nre,pkamatr)

	ofile101.close()
#-----------------------------------------------------------------------------
	ofile1000 = open ('ToAddAll.txt', 'a')
	print('', file = ofile1000)
	for i in range (nre):
		print (['{:.6E}'.format(pkamatr[i][j]) for j in range (nre)], file = ofile1000)
	ofile1000.close()
	
# PKAS_ng function completes here

#------------------------------------------------------------------------------

#=======File 1 information=======*
	## The directory of Files (MFs) and the corresponding Sections (MTs)
	## given in the evaluation are read from File 1 and returned.
	
def FILE1 ():
	MFs = [0]*1000; MTs = [0]*1000 
		# maximum of NXC = 350 (ENDF-102), 
		# but deviates for Mn55 ENDF/B-VII.1, so changed to 1000 
	ifile = open('tape01','r')
	ifile.readline()
	line = ifile.readline()
	(ZA,AWR,LRP,LFI,NLIB,NMOD,MAT,MF,MT) = line_type1_info(line)
	line = ifile.readline()
	(ELIS,STA,LIS,LISO,num,NFOR,MAT,MF,MT) = line_type1_info(line)
	line = ifile.readline()
	(AWI,EMAX,LREL,num,NSUB,NVER,MAT,MF,MT) = line_type1_info(line)
	line = ifile.readline()
	(TEMP,c2,LDRV,num,NWD,NXC,MAT,MF,MT) = line_type1_info(line)
	for i in range (NWD):
		ifile.readline()
	for i in range (NXC):
		line = ifile.readline()
		data = eachlineinfo(line)
		blnk = data[0]; blnk = data[1]
		MFs[i] = int(data[2]); MTs[i] = int(data[3])
		NCn = int(data[4]); MODn = int(data[5]); MAT = int(data[6])
		MF = int(data[7]); MT = int(data[8])

	nfiles = NXC
	ifile.close()

	return (nfiles,MFs,MTs)

#******************************** WEIGHT SPECTRUM GENERATOR **************************		
 		
	# Spectrum to weight and perform the multigrouping.
 	
def spectrum(En,L):
	if (L == 1):
		fi = En * math.exp(-En/0.0253)
	if (L == 2):
		fi = 1/En
	if (L == 3):
		fi = math.sqrt(En) * math.exp(-En/(1.415e+6))
	if (L == 4):
		fi = 1
	else:
		fi = 0
	return(fi)

#=======Energy multigrouping=======*		
		
def GROUPMULTI (insp,igtype,mttg,nrg):
	gsdpa = [0]*nrg; Eg = [0]*nrg
	Ngl = nrg
	
	if (insp == 1):
		ifile407 = open('NeutronSpectrum.txt', 'r')
		nre = int(ifile407.readline().split()[-1])
		fi = numpy.zeros(nrg)
		for i in reversed(range(nre)):
			fi[i] = float(ifile407.readline().split()[0])
		fi[-1] = fi[-2]
		ifile407.close()
	
	ifile54 = open ("tape02", 'r')
	while True:
		line = ifile54.readline()
		if (line == ''):
			break
		data = eachlineinfo(line)
		MAT = int(data[6]); MF = int(data[7]); MT  = int(data[8])
		if (MAT != -1):
			if (MF == 3):
				if (MT == mttg):
					line = ifile54.readline()
					data = eachlineinfo(line)
					QM = float(data[0]); QI =  float(data[1]); NR = int(data[4]); NP = int(data[5])
					E = [0]*NP; sdpa = [0]*NP
					LR = int(ifile54.readline().split()[3])
					(E,sdpa) = line_type3_info(ifile54,NP,2)
		else:
			break
	ifile54.close()

	if (igtype == 0):
		ifile406 = open('Energy-GroupLimits.txt', 'r')
		Ngl = int(ifile406.readline().split()[-1])
		Eg = [0]*Ngl; gsdpa = [0]*Ngl
		for i in reversed(range(Ngl)):
			Eg[i] = float(ifile406.readline().split()[0])
		ifile406.close()
	
	if (igtype==1):
		Ngl = 176
		nre = Ngl-1
		Eg = numpy.zeros(Ngl); gsdpa = numpy.zeros(Ngl)
		Eg = engrp1()

	if (igtype==2):
		Ngl = 27
		nre = Ngl-1
		Eg = numpy.zeros(Ngl); gsdpa = numpy.zeros(Ngl)
		Eg = engrp2()

	if (igtype==3):
		Ngl = 34
		nre = Ngl-1
		Eg = numpy.zeros(Ngl); gsdpa = numpy.zeros(Ngl)
		Eg = engrp3()
            
	if (igtype==4):
		Ngl = 239
		nre = Ngl-1
		Eg = numpy.zeros(Ngl); gsdpa = numpy.zeros(Ngl)
		Eg = engrp4()
            
	if (igtype==5): 
		Ngl = 199
		nre = Ngl-1
		Eg = numpy.zeros(Ngl); gsdpa = numpy.zeros(Ngl)
		Eg = engrp5()
            
	if (igtype==6):
		Ngl = 710
		nre = Ngl-1
		Eg = numpy.zeros(Ngl); gsdpa = numpy.zeros(Ngl)
		Eg = engrp6()
            
	if (igtype==7):
		Ngl = 641
		nre = Ngl-1
		Eg = numpy.zeros(Ngl); gsdpa = numpy.zeros(Ngl)
		Eg = engrp7()
            
	if (igtype==8):
		Ngl = 101
		nre = Ngl-1
		Eg = numpy.zeros(Ngl); gsdpa = numpy.zeros(Ngl)
		Eg = engrp8()
            
	if (igtype==9):
		Ngl = 48
		nre = Ngl-1
		Eg = numpy.zeros(Ngl); gsdpa = numpy.zeros(Ngl)
		Eg = engrp9()
            
	if (igtype==10):
		Ngl = 101			# DLC-2 group structure
		nre = Ngl-1
		Eg = numpy.zeros(Ngl); gsdpa = numpy.zeros(Ngl)
		Eg = engrp10()
			
	if (igtype==11):
		Ngl = 229			# 229 group structure
		nre = Ngl-1
		Eg = numpy.zeros(Ngl); gsdpa = numpy.zeros(Ngl)
		Eg = engrp11()
            
	if (igtype==12):
		Ngl = 229			# 229 group structure
		nre = Ngl-1
		Eg = numpy.zeros(Ngl); gsdpa = numpy.zeros(Ngl)
		Eg = engrp12()

	ifg = 0
	for i in range (Ngl-1):
		if (Eg[i] <= E[0] and E[0] <= Eg[i+1]):
			ifg = i
			break

	for i in range (ifg, Ngl-1):
		Eg1 = Eg[i]
		Eg2 = Eg[i+1]
		Nsect = 100
		h = (Eg2-Eg1)/Nsect
		s1 = 0; f1 = 0; f2 = 0
		s2 = 0; g1 = 0; g2 = 0
		for j in range (Nsect):
			t = Eg1 + j*h
			if (t < 0.1):
				L = 1
			if (0.1 <= t and t < 820.3e+03):
				L = 2
			if (t >= 820.3e+03):
				L = 3
			if (insp == 0):
				flux1 = spectrum (t,L)
			if (insp == 1):
				flux1 = numpy.interp(t, Eg, fi)	#srchintrp3 (Eg,fi,Ngl,t)
			
			sigma1 = numpy.interp(t, E, sdpa)
			f1 = sigma1*flux1
			g1 = flux1
				
			if (t+h < 0.1):
				L = 1
			if (0.1 <= t+h and t+h < 820.3e+03):
				L = 2
			if (t+h >= 820.3e+03):
				L = 3
			if (insp == 0):
				flux2 = spectrum (t+h,L)
			if (insp == 1):
				flux2 = numpy.interp(t+h, Eg, fi)	#srchintrp3 (Eg,fi,Ngl,t+h)
			
			sigma2 = numpy.interp(t+h, E, sdpa)
			f2 = sigma2*flux2
			g2 = flux2

			s1 = s1 + (h/2)*(f1+f2)
			s2 = s2 + (h/2)*(g1+g2)

		if (s1 != 0 and s2 != 0):
			gsdpa[i] = s1/s2
	
	return (gsdpa)

# GROUPMULTI function completes here

#==================================================================
	
	## The subroutines engrp1, engrp2, etc. are different types of 
	## energy group structures.
	
def engrp1():
	Ngl = 176
	Eg = [1.000E-05,1.000E-01,4.140E-01,5.320E-01,6.830E-01,8.760E-01,
	1.130E+00,1.450E+00,1.860E+00,2.380E+00,3.060E+00,3.930E+00,
	5.040E+00,6.480E+00,8.320E+00,1.070E+01,1.370E+01,1.760E+01,
	2.260E+01,2.900E+01,3.730E+01,4.790E+01,6.140E+01,7.890E+01,
	1.010E+02,1.300E+02,1.670E+02,2.140E+02,2.750E+02,3.540E+02,
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
	1.490E+07,1.570E+07,1.650E+07,1.690E+07,1.730E+07,1.960E+07]
		
	return(Eg)

def engrp2():
	Ngl = 27
	Eg = [1.0E-05,0.0253E+00,0.4642E+00,1.0E+00,2.1544E+00,
	4.6416E+00,1.0E+01,2.1544E+01,4.6416E+01,1.0E+02,2.1544E+02,
	4.6416E+02,1.0E+03,2.1544E+03,4.6416E+03,1.0E+04,2.1544E+04,
	4.6416E+04,1.0E+05,2.0E+05,0.4E+06,0.8E+06,1.4E+06,2.5E+06,
	4.0E+06,6.5E+06,1.05E+07]
	return(Eg)

def engrp3():
	Ngl = 34
	Eg = [1.000e-05,1.000e-01,5.400e-01,4.000e+00,8.315e+00,
	1.371e+01,2.260e+01,4.017e+01,6.790e+01,9.166e+01,1.486e+02,
	3.043e+02,4.540e+02,7.485e+02,1.230e+03,2.030e+03,3.355e+03,
	5.531e+03,9.119e+03,1.503e+04,2.479e+04,4.087e+04,6.738e+04,
	1.111e+05,1.832e+05,3.020e+05,4.979e+05,8.209e+05,1.353e+06,
	2.231e+06,3.679e+06,6.065e+06,1.000e+07,1.964e+07]
	return(Eg)

def engrp4():
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
	
	return(Eg)

def engrp5 ():
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
	
	return(Eg)

def engrp6():
	Ngl = 710
	Eg = [1.0E-05,1.047129E-05,1.096478E-05,1.148154E-05,1.202264E-05,
	1.258925E-05,1.318257E-05,1.380384E-05,1.445440E-05,
	1.513561E-05,1.584893E-05,1.659587E-05,1.737801E-05,
	1.819701E-05,1.905461E-05,1.995262E-05,2.089296E-05,
	2.187762E-05,2.290868E-05,2.398833E-05,2.511886E-05,
	2.630268E-05,2.754229E-05,2.884032E-05,3.019952E-05,
	3.162278E-05,3.311311E-05,3.467369E-05,3.630781E-05,
	3.801894E-05,3.981072E-05,4.168694E-05,4.365158E-05,
	4.570882E-05,4.786301E-05,5.011872E-05,5.248075E-05,
	5.495409E-05,5.754399E-05,6.025596E-05,6.309573E-05,
	6.606934E-05,6.918310E-05,7.244360E-05,7.585776E-05,
	7.943282E-05,8.317638E-05,8.709636E-05,9.120108E-05,
	9.549926E-05,1.000000E-04,1.047129E-04,1.096478E-04,
	1.148154E-04,1.202264E-04,1.258925E-04,1.318257E-04,
	1.380384E-04,1.445440E-04,1.513561E-04,1.584893E-04,
	1.659587E-04,1.737801E-04,1.819701E-04,1.905461E-04,
	1.995262E-04,2.089296E-04,2.187762E-04,2.290868E-04,
	2.398833E-04,2.511886E-04,2.630268E-04,2.754229E-04,
	2.884032E-04,3.019952E-04,3.162278E-04,3.311311E-04,
	3.467369E-04,3.630781E-04,3.801894E-04,3.981072E-04,
	4.168694E-04,4.365158E-04,4.570882E-04,4.786301E-04,
	5.011872E-04,5.248075E-04,5.495409E-04,5.754399E-04,
	6.025596E-04,6.309573E-04,6.606934E-04,6.918310E-04,
	7.244360E-04,7.585776E-04,7.943282E-04,8.317638E-04,
	8.709636E-04,9.120108E-04,9.549926E-04,1.000000E-03,
	1.047129E-03,1.096478E-03,1.148154E-03,1.202264E-03,
	1.258925E-03,1.318257E-03,1.380384E-03,1.445440E-03,
	1.513561E-03,1.584893E-03,1.659587E-03,1.737801E-03,
	1.819701E-03,1.905461E-03,1.995262E-03,2.089296E-03,
	2.187762E-03,2.290868E-03,2.398833E-03,2.511886E-03,
	2.630268E-03,2.754229E-03,2.884032E-03,3.019952E-03,
	3.162278E-03,3.311311E-03,3.467369E-03,3.630781E-03,
	3.801894E-03,3.981072E-03,4.168694E-03,4.365158E-03,
	4.570882E-03,4.786301E-03,5.011872E-03,5.248075E-03,
	5.495409E-03,5.754399E-03,6.025596E-03,6.309573E-03,
	6.606934E-03,6.918310E-03,7.244360E-03,7.585776E-03,
	7.943282E-03,8.317638E-03,8.709636E-03,9.120108E-03,
	9.549926E-03,1.000000E-02,1.047129E-02,1.096478E-02,
	1.148154E-02,1.202264E-02,1.258925E-02,1.318257E-02,
	1.380384E-02,1.445440E-02,1.513561E-02,1.584893E-02,
	1.659587E-02,1.737801E-02,1.819701E-02,1.905461E-02,
	1.995262E-02,2.089296E-02,2.187762E-02,2.290868E-02,
	2.398833E-02,2.511886E-02,2.630268E-02,2.754229E-02,
	2.884032E-02,3.019952E-02,3.162278E-02,3.311311E-02,
	3.467369E-02,3.630781E-02,3.801894E-02,3.981072E-02,
	4.168694E-02,4.365158E-02,4.570882E-02,4.786301E-02,
	5.011872E-02,5.248075E-02,5.495409E-02,5.754399E-02,
	6.025596E-02,6.309573E-02,6.606934E-02,6.918310E-02,
	7.244360E-02,7.585776E-02,7.943282E-02,8.317638E-02,
	8.709636E-02,9.120108E-02,9.549926E-02,1.000000E-01,
	1.047129E-01,1.096478E-01,1.148154E-01,1.202264E-01,
	1.258925E-01,1.318257E-01,1.380384E-01,1.445440E-01,
	1.513561E-01,1.584893E-01,1.659587E-01,1.737801E-01,
	1.819701E-01,1.905461E-01,1.995262E-01,2.089296E-01,
	2.187762E-01,2.290868E-01,2.398833E-01,2.511886E-01,
	2.630268E-01,2.754229E-01,2.884032E-01,3.019952E-01,
	3.162278E-01,3.311311E-01,3.467369E-01,3.630781E-01,
	3.801894E-01,3.981072E-01,4.168694E-01,4.365158E-01,
	4.570882E-01,4.786301E-01,5.011872E-01,5.248075E-01,
	5.495409E-01,5.754399E-01,6.025596E-01,6.309573E-01,
	6.606934E-01,6.918310E-01,7.244360E-01,7.585776E-01,
	7.943282E-01,8.317638E-01,8.709636E-01,9.120108E-01,
	9.549926E-01,1.000000E+00,1.047129E+00,1.096478E+00,
	1.148154E+00,1.202264E+00,1.258925E+00,1.318257E+00,
	1.380384E+00,1.445440E+00,1.513561E+00,1.584893E+00,
	1.659587E+00,1.737801E+00,1.819701E+00,1.905461E+00,
	1.995262E+00,2.089296E+00,2.187762E+00,2.290868E+00,
	2.398833E+00,2.511886E+00,2.630268E+00,2.754229E+00,
	2.884032E+00,3.019952E+00,3.162278E+00,3.311311E+00,
	3.467369E+00,3.630781E+00,3.801894E+00,3.981072E+00,
	4.168694E+00,4.365158E+00,4.570882E+00,4.786301E+00,
	5.011872E+00,5.248075E+00,5.495409E+00,5.754399E+00,
	6.025596E+00,6.309573E+00,6.606934E+00,6.918310E+00,
	7.244360E+00,7.585776E+00,7.943282E+00,8.317638E+00,
	8.709636E+00,9.120108E+00,9.549926E+00,1.000000E+01,
	1.047129E+01,1.096478E+01,1.148154E+01,1.202264E+01,
	1.258925E+01,1.318257E+01,1.380384E+01,1.445440E+01,
	1.513561E+01,1.584893E+01,1.659587E+01,1.737801E+01,
	1.819701E+01,1.905461E+01,1.995262E+01,2.089296E+01,
	2.187762E+01,2.290868E+01,2.398833E+01,2.511886E+01,
	2.630268E+01,2.754229E+01,2.884032E+01,3.019952E+01,
	3.162278E+01,3.311311E+01,3.467369E+01,3.630781E+01,
	3.801894E+01,3.981072E+01,4.168694E+01,4.365158E+01,
	4.570882E+01,4.786301E+01,5.011872E+01,5.248075E+01,
	5.495409E+01,5.754399E+01,6.025596E+01,6.309573E+01,
	6.606934E+01,6.918310E+01,7.244360E+01,7.585776E+01,
	7.943282E+01,8.317638E+01,8.709636E+01,9.120108E+01,
	9.549926E+01,1.000000E+02,1.047129E+02,1.096478E+02,
	1.148154E+02,1.202264E+02,1.258925E+02,1.318257E+02,
	1.380384E+02,1.445440E+02,1.513561E+02,1.584893E+02,
	1.659587E+02,1.737801E+02,1.819701E+02,1.905461E+02,
	1.995262E+02,2.089296E+02,2.187762E+02,2.290868E+02,
	2.398833E+02,2.511886E+02,2.630268E+02,2.754229E+02,
	2.884032E+02,3.019952E+02,3.162278E+02,3.311311E+02,
	3.467369E+02,3.630781E+02,3.801894E+02,3.981072E+02,
	4.168694E+02,4.365158E+02,4.570882E+02,4.786301E+02,
	5.011872E+02,5.248075E+02,5.495409E+02,5.754399E+02,
	6.025596E+02,6.309573E+02,6.606934E+02,6.918310E+02,
	7.244360E+02,7.585776E+02,7.943282E+02,8.317638E+02,
	8.709636E+02,9.120108E+02,9.549926E+02,1.000000E+03,
	1.047129E+03,1.096478E+03,1.148154E+03,1.202264E+03,
	1.258925E+03,1.318257E+03,1.380384E+03,1.445440E+03,
	1.513561E+03,1.584893E+03,1.659587E+03,1.737801E+03,
	1.819701E+03,1.905461E+03,1.995262E+03,2.089296E+03,
	2.187762E+03,2.290868E+03,2.398833E+03,2.511886E+03,
	2.630268E+03,2.754229E+03,2.884032E+03,3.019952E+03,
	3.162278E+03,3.311311E+03,3.467369E+03,3.630781E+03,
	3.801894E+03,3.981072E+03,4.168694E+03,4.365158E+03,
	4.570882E+03,4.786301E+03,5.011872E+03,5.248075E+03,
	5.495409E+03,5.754399E+03,6.025596E+03,6.309573E+03,
	6.606934E+03,6.918310E+03,7.244360E+03,7.585776E+03,
	7.943282E+03,8.317638E+03,8.709636E+03,9.120108E+03,
	9.549926E+03,1.000000E+04,1.047129E+04,1.096478E+04,
	1.148154E+04,1.202264E+04,1.258925E+04,1.318257E+04,
	1.380384E+04,1.445440E+04,1.513561E+04,1.584893E+04,
	1.659587E+04,1.737801E+04,1.819701E+04,1.905461E+04,
	1.995262E+04,2.089296E+04,2.187762E+04,2.290868E+04,
	2.398833E+04,2.511886E+04,2.630268E+04,2.754229E+04,
	2.884032E+04,3.019952E+04,3.162278E+04,3.311311E+04,
	3.467369E+04,3.630781E+04,3.801894E+04,3.981072E+04,
	4.168694E+04,4.365158E+04,4.570882E+04,4.786301E+04,
	5.011872E+04,5.248075E+04,5.495409E+04,5.754399E+04,
	6.025596E+04,6.309573E+04,6.606934E+04,6.918310E+04,
	7.244360E+04,7.585776E+04,7.943282E+04,8.317638E+04,
	8.709636E+04,9.120108E+04,9.549926E+04,1.000000E+05,
	1.047129E+05,1.096478E+05,1.148154E+05,1.202264E+05,
	1.258925E+05,1.318257E+05,1.380384E+05,1.445440E+05,
	1.513561E+05,1.584893E+05,1.659587E+05,1.737801E+05,
	1.819701E+05,1.905461E+05,1.995262E+05,2.089296E+05,
	2.187762E+05,2.290868E+05,2.398833E+05,2.511886E+05,
	2.630268E+05,2.754229E+05,2.884032E+05,3.019952E+05,
	3.162278E+05,3.311311E+05,3.467369E+05,3.630781E+05,
	3.801894E+05,3.981072E+05,4.168694E+05,4.365158E+05,
	4.570882E+05,4.786301E+05,5.011872E+05,5.248075E+05,
	5.495409E+05,5.754399E+05,6.025596E+05,6.309573E+05,
	6.606934E+05,6.918310E+05,7.244360E+05,7.585776E+05,
	7.943282E+05,8.317638E+05,8.709636E+05,9.120108E+05,
	9.549926E+05,1.000000E+06,1.047129E+06,1.096478E+06,
	1.148154E+06,1.202264E+06,1.258925E+06,1.318257E+06,
	1.380384E+06,1.445440E+06,1.513561E+06,1.584893E+06,
	1.659587E+06,1.737801E+06,1.819701E+06,1.905461E+06,
	1.995262E+06,2.089296E+06,2.187762E+06,2.290868E+06,
	2.398833E+06,2.511886E+06,2.630268E+06,2.754229E+06,
	2.884032E+06,3.019952E+06,3.162278E+06,3.311311E+06,
	3.467369E+06,3.630781E+06,3.801894E+06,3.981072E+06,
	4.168694E+06,4.365158E+06,4.570882E+06,4.786301E+06,
	5.011872E+06,5.248075E+06,5.495409E+06,5.754399E+06,
	6.025596E+06,6.309573E+06,6.606934E+06,6.918310E+06,
	7.244360E+06,7.585776E+06,7.943282E+06,8.317638E+06,
	8.709636E+06,9.120108E+06,9.549926E+06,1.00E+07,
	1.02E+07,1.04E+07,1.06E+07,1.08E+07,1.10E+07,1.12E+07,
	1.14E+07,1.16E+07,1.18E+07,1.20E+07,1.22E+07,1.24E+07,
	1.26E+07,1.28E+07,1.30E+07,1.32E+07,1.34E+07,1.36E+07,
	1.38E+07,1.40E+07,1.42E+07,1.44E+07,1.46E+07,1.48E+07,
	1.50E+07,1.52E+07,1.54E+07,1.56E+07,1.58E+07,1.60E+07,
	1.62E+07,1.64E+07,1.66E+07,1.68E+07,1.70E+07,1.72E+07,
	1.74E+07,1.76E+07,1.78E+07,1.80E+07,1.82E+07,1.84E+07,
	1.86E+07,1.88E+07,1.90E+07,1.92E+07,1.94E+07,1.96E+07,
	1.98E+07,2.00E+07,2.10E+07,2.20E+07,2.3E+07,2.4E+07,
	2.5E+07,2.6E+07,2.7E+07,2.8E+07,2.9E+07,3.0E+07,3.2E+07,
	3.4E+07,3.6E+07,3.8E+07,4.0E+07,4.2E+07,4.4E+07,4.6E+07,
	4.8E+07,5.0E+07,5.2E+07,5.4E+07,5.6E+07,5.8E+07,6.0E+07,
	6.5E+07,7.0E+07,7.5E+07,8.0E+07,9.0E+07,1.0E+08,1.1E+08,
	1.2E+08,1.3E+08,1.4E+08,1.5E+08,1.6E+08,1.8E+08,2.0E+08,
	2.4E+08,2.8E+08,3.2E+08,3.6E+08,4.0E+08,4.4E+08,4.8E+08,
	5.2E+08,5.6E+08,6.0E+08,6.4E+08,6.8E+08,7.2E+08,7.6E+08,
	8.0E+08,8.4E+08,8.8E+08,9.2E+08,9.6E+08,1.0E+09]
	
	return(Eg)

def engrp7():
	Ngl = 641
	Eg = [1.000E-05,1.050E-04,1.100E-04,1.150E-04,1.200E-04,
	1.280E-04,1.350E-04,1.430E-04,1.500E-04,1.600E-04,1.700E-04,
	1.800E-04,1.900E-04,2.000E-04,2.100E-04,2.200E-04,2.300E-04,
	2.400E-04,2.550E-04,2.700E-04,2.800E-04,3.000E-04,3.200E-04,
	3.400E-04,3.600E-04,3.800E-04,4.000E-04,4.250E-04,4.500E-04,
	4.750E-04,5.000E-04,5.250E-04,5.500E-04,5.750E-04,6.000E-04,
	6.300E-04,6.600E-04,6.900E-04,7.200E-04,7.600E-04,8.000E-04,
	8.400E-04,8.800E-04,9.200E-04,9.600E-04,1.000E-03,1.050E-03,
	1.100E-03,1.150E-03,1.200E-03,1.280E-03,1.350E-03,1.430E-03,
	1.500E-03,1.600E-03,1.700E-03,1.800E-03,1.900E-03,2.000E-03,
	2.100E-03,2.200E-03,2.300E-03,2.400E-03,2.550E-03,2.700E-03,
	2.800E-03,3.000E-03,3.200E-03,3.400E-03,3.600E-03,3.800E-03,
	4.000E-03,4.250E-03,4.500E-03,4.750E-03,5.000E-03,5.250E-03,
	5.500E-03,5.750E-03,6.000E-03,6.300E-03,6.600E-03,6.900E-03,
	7.200E-03,7.600E-03,8.000E-03,8.400E-03,8.800E-03,9.200E-03,
	9.600E-03,1.000E-02,1.050E-02,1.100E-02,1.150E-02,1.200E-02,
	1.280E-02,1.310E-02,1.430E-02,1.500E-02,1.600E-02,1.700E-02,
	1.800E-02,1.900E-02,2.000E-02,2.100E-02,2.200E-02,2.300E-02,
	2.400E-02,2.550E-02,2.700E-02,2.800E-02,3.000E-02,3.200E-02,
	3.400E-02,3.600E-02,3.800E-02,4.000E-02,4.250E-02,4.500E-02,
	4.750E-02,5.000E-02,5.250E-02,5.500E-02,5.750E-02,6.000E-02,
	6.300E-02,6.600E-02,6.900E-02,7.200E-02,7.600E-02,8.000E-02,
	8.400E-02,8.800E-02,9.200E-02,9.600E-02,1.000E-01,1.050E-01,
	1.100E-01,1.150E-01,1.200E-01,1.280E-01,1.350E-01,1.430E-01,
	1.500E-01,1.600E-01,1.700E-01,1.800E-01,1.900E-01,2.000E-01,
	2.100E-01,2.200E-01,2.300E-01,2.400E-01,2.550E-01,2.700E-01,
	2.800E-01,3.000E-01,3.200E-01,3.400E-01,3.600E-01,3.800E-01,
	4.000E-01,4.250E-01,4.500E-01,4.750E-01,5.000E-01,5.250E-01,
	5.500E-01,5.750E-01,6.000E-01,6.300E-01,6.600E-01,6.900E-01,
	7.200E-01,7.600E-01,8.000E-01,8.400E-01,8.800E-01,9.200E-01,
	9.600E-01,1.000E+00,1.050E+00,1.100E+00,1.150E+00,1.200E+00,
	1.280E+00,1.350E+00,1.430E+00,1.500E+00,1.600E+00,1.700E+00,
	1.800E+00,1.900E+00,2.000E+00,2.100E+00,2.200E+00,2.300E+00,
	2.400E+00,2.550E+00,2.700E+00,2.800E+00,3.000E+00,3.200E+00,
	3.400E+00,3.600E+00,3.800E+00,4.000E+00,4.250E+00,4.500E+00,
	4.750E+00,5.000E+00,5.250E+00,5.500E+00,5.750E+00,6.000E+00,
	6.300E+00,6.600E+00,6.900E+00,7.200E+00,7.600E+00,8.000E+00,
	8.400E+00,8.800E+00,9.200E+00,9.600E+00,1.000E+01,1.050E+01,
	1.100E+01,1.150E+01,1.200E+01,1.280E+01,1.350E+01,1.430E+01,
	1.500E+01,1.600E+01,1.700E+01,1.800E+01,1.900E+01,2.000E+01,
	2.100E+01,2.200E+01,2.300E+01,2.400E+01,2.550E+01,2.700E+01,
	2.800E+01,3.000E+01,3.200E+01,3.400E+01,3.600E+01,3.800E+01,
	4.000E+01,4.250E+01,4.500E+01,4.750E+01,5.000E+01,5.250E+01,
	5.500E+01,5.750E+01,6.000E+01,6.300E+01,6.600E+01,6.900E+01,
	7.200E+01,7.600E+01,8.000E+01,8.400E+01,8.800E+01,9.200E+01,
	9.600E+01,1.000E+02,1.050E+02,1.100E+02,1.150E+02,1.200E+02,
	1.280E+02,1.350E+02,1.430E+02,1.500E+02,1.600E+02,1.700E+02,
	1.800E+02,1.900E+02,2.000E+02,2.100E+02,2.200E+02,2.300E+02,
	2.400E+02,2.550E+02,2.700E+02,2.800E+02,3.000E+02,3.200E+02,
	3.400E+02,3.600E+02,3.800E+02,4.000E+02,4.250E+02,4.500E+02,
	4.750E+02,5.000E+02,5.250E+02,5.500E+02,5.750E+02,6.000E+02,
	6.300E+02,6.600E+02,6.900E+02,7.200E+02,7.600E+02,8.000E+02,
	8.400E+02,8.800E+02,9.200E+02,9.600E+02,1.000E+03,1.050E+03,
	1.100E+03,1.150E+03,1.200E+03,1.280E+03,1.350E+03,1.430E+03,
	1.500E+03,1.600E+03,1.700E+03,1.800E+03,1.900E+03,2.000E+03,
	2.010E+03,2.200E+03,2.300E+03,2.400E+03,2.550E+03,2.700E+03,
	2.800E+03,3.000E+03,3.200E+03,3.400E+03,3.600E+03,3.800E+03,
	4.000E+03,4.250E+03,4.500E+03,4.750E+03,5.000E+03,5.250E+03,
	5.500E+03,5.750E+03,6.000E+03,6.300E+03,6.600E+03,6.900E+03,
	7.200E+03,7.600E+03,8.000E+03,8.400E+03,8.800E+03,9.200E+03,
	9.600E+03,1.000E+04,1.050E+04,1.100E+04,1.150E+04,1.200E+04,
	1.280E+04,1.350E+04,1.430E+04,1.500E+04,1.600E+04,1.700E+04,
	1.800E+04,1.900E+04,2.000E+04,2.100E+04,2.200E+04,2.300E+04,
	2.400E+04,2.550E+04,2.700E+04,2.800E+04,3.000E+04,3.200E+04,
	3.400E+04,3.600E+04,3.800E+04,4.000E+04,4.250E+04,4.500E+04,
	4.750E+04,5.000E+04,5.250E+04,5.500E+04,5.750E+04,6.000E+04,
	6.300E+04,6.600E+04,6.900E+04,7.200E+04,7.600E+04,8.000E+04,
	8.400E+04,8.800E+04,9.200E+04,9.600E+04,1.000E+05,1.050E+05,
	1.100E+05,1.150E+05,1.200E+05,1.275E+05,1.350E+05,1.425E+05,
	1.500E+05,1.600E+05,1.700E+05,1.800E+05,1.900E+05,2.000E+05,
	2.100E+05,2.200E+05,2.300E+05,2.400E+05,2.550E+05,2.700E+05,
	2.800E+05,3.000E+05,3.200E+05,3.400E+05,3.600E+05,3.800E+05,
	4.000E+05,4.250E+05,4.500E+05,4.750E+05,5.000E+05,5.250E+05,
	5.500E+05,5.750E+05,6.000E+05,6.300E+05,6.600E+05,6.900E+05,
	7.200E+05,7.600E+05,8.000E+05,8.400E+05,8.800E+05,9.200E+05,
	9.600E+05,1.000E+06,1.100E+06,1.200E+06,1.300E+06,1.400E+06,
	1.500E+06,1.600E+06,1.700E+06,1.800E+06,1.900E+06,2.000E+06,
	2.100E+06,2.200E+06,2.300E+06,2.400E+06,2.500E+06,2.600E+06,
	2.700E+06,2.800E+06,2.900E+06,3.000E+06,3.100E+06,3.200E+06,
	3.300E+06,3.400E+06,3.500E+06,3.600E+06,3.700E+06,3.800E+06,
	3.900E+06,4.000E+06,4.100E+06,4.200E+06,4.300E+06,4.400E+06,
	4.500E+06,4.600E+06,4.700E+06,4.800E+06,4.900E+06,5.000E+06,
	5.100E+06,5.200E+06,5.300E+06,5.400E+06,5.500E+06,5.600E+06,
	5.700E+06,5.800E+06,5.900E+06,6.000E+06,6.100E+06,6.200E+06,
	6.300E+06,6.400E+06,6.500E+06,6.600E+06,6.700E+06,6.800E+06,
	6.900E+06,7.000E+06,7.100E+06,7.200E+06,7.300E+06,7.400E+06,
	7.500E+06,7.600E+06,7.700E+06,7.800E+06,7.900E+06,8.000E+06,
	8.100E+06,8.200E+06,8.300E+06,8.400E+06,8.500E+06,8.600E+06,
	8.700E+06,8.800E+06,8.900E+06,9.000E+06,9.100E+06,9.200E+06,
	9.300E+06,9.400E+06,9.500E+06,9.600E+06,9.700E+06,9.800E+06,
	9.900E+06,1.000E+07,1.010E+07,1.020E+07,1.030E+07,1.040E+07,
	1.050E+07,1.060E+07,1.070E+07,1.080E+07,1.090E+07,1.100E+07,
	1.110E+07,1.120E+07,1.130E+07,1.140E+07,1.150E+07,1.160E+07,
	1.170E+07,1.180E+07,1.190E+07,1.200E+07,1.210E+07,1.220E+07,
	1.230E+07,1.240E+07,1.250E+07,1.260E+07,1.270E+07,1.280E+07,
	1.290E+07,1.300E+07,1.310E+07,1.320E+07,1.330E+07,1.340E+07,
	1.350E+07,1.360E+07,1.370E+07,1.380E+07,1.390E+07,1.400E+07,
	1.410E+07,1.420E+07,1.430E+07,1.440E+07,1.450E+07,1.460E+07,
	1.470E+07,1.480E+07,1.490E+07,1.500E+07,1.510E+07,1.520E+07,
	1.530E+07,1.540E+07,1.550E+07,1.560E+07,1.570E+07,1.580E+07,
	1.590E+07,1.600E+07,1.610E+07,1.620E+07,1.630E+07,1.640E+07,
	1.650E+07,1.660E+07,1.670E+07,1.680E+07,1.690E+07,1.700E+07,
	1.710E+07,1.720E+07,1.730E+07,1.740E+07,1.750E+07,1.760E+07,
	1.770E+07,1.780E+07,1.790E+07,1.800E+07,1.810E+07,1.820E+07,
	1.830E+07,1.840E+07,1.850E+07,1.860E+07,1.870E+07,1.880E+07,
	1.890E+07,1.900E+07,1.910E+07,1.920E+07,1.930E+07,1.940E+07,
	1.950E+07,1.960E+07,1.970E+07,1.980E+07,1.990E+07,2.000E+07]

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

def engrp9():
	Ngl = 48
	Eg = [1.00E-05,1.24E-02,3.06E-02,4.28E-02,5.69E-02,8.20E-02,
		1.12E-01,1.46E-01,1.84E-01,2.71E-01,3.58E-01,5.03E-01,
		6.25E-01,7.82E-01,9.10E-01,9.71E-01,1.01E+00,1.07E+00,
		1.13E+00,1.17E+00,1.24E+00,1.46E+00,1.86E+00,2.38E+00,
		3.93E+00,4.45E+00,5.04E+00,5.72E+00,6.48E+00,7.34E+00,
		8.32E+00,1.21E+01,1.37E+01,2.90E+01,4.79E+01,7.89E+01,
		1.30E+02,2.03E+03,9.12E+03,6.74E+04,1.83E+05,4.98E+05,
		8.21E+05,1.35E+06,2.23E+06,3.68E+06,6.07E+06,2.00E+07]

	return(Eg)
	
def engrp10():	# DLC-2
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
	
	return(Eg)

def engrp11():	# EXTENDED FROM 198 GROUP STRUCTURE
	# BY SELF UP TO 200 MeV
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

	return(Eg)

def engrp12():	# EXTENDED FROM 198 GROUP STRUCTURE
# BY SELF UP TO 150 MeV
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

	return(Eg)
