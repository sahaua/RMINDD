## >> Code: RMINDD.py, module file: EngdepU.py
## >> Perform: Recoil energy distributions of nuclei from basic ENDF-6 files
## >> Author: Uttiyoarnab Saha
## >> Version and Date: 1.0 and 25/03/2021
## >> Last modified: 25/03/2021, Kolkata
## >> Update: 25/03/2021
## >> Major changes:
##
## ======================================================================================================

import numpy
import math
import UtilsU

#=======Linear interpolation=======*		
        
def crstd(x,x1,x2,y1,y2):
	if ((x2-x1) != 0):
		y = y1 + ((y2-y1)*(x-x1)/(x2-x1))
	return(y)

# =======Interpolate cross section to unique common energy=======*
	
	# When the basic / dpa / heating cross sections from partial reactions
	# have to be computed corresponding to the unique energy points from
	# from their individual energy arrays, linear interpolation is performed
	# to find the cross sections corresponding to unique energy points.
	
def trptuqce (E,s1,Etu):
	s2 = numpy.interp(Etu, E, s1)
	return(s2)

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

def line_type3_info(filehandle,numdata,numvariables):
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

#=======Integrate over displacement model:: 'al' coeffcients =======*
		
	# Average damage energy for discrete level reactions from Legendre
	# expansion coefficients representation of secondary particle angular 
	# distribution.
	
def discrec1(alc1,NLa,E,Q,bta,Z1,A1,Z2,A2,Ed,mdisp,bad,cad):
	nquad = 64
	Pl1 = [0]*NLa
	xabc = [0]*nquad; wg = [0]*nquad
	A = A2
	U = (A+1)*abs(Q)/A
	if (Q > 0):
		U = 0
	R = math.sqrt((A*(A+1-bta)/bta)*(1-U/E))
	R1 = bta*R/(A+1-bta)

	(xabc, wg) = GQ()	# Gauss-quadrature values and weights

	s = 0
	for i in range(64):
		Pl1[0] = 1
		Pl1[1] = xabc[i]
		T = E*(A+1-bta)*(1+R1*R1-2*R1*xabc[i])/(A+1)**2
		if (NLa > 3):
			for j in range(2, NLa):
				Pl1[j] = (((2 * (j-1)+1) * xabc[i] * Pl1[j-1])-((j-1) * Pl1[j-2])) / j
			p1 = 0
			for l in range(NLa):
				p1 = p1 + (((2 * l)+1) * Pl1[l] * alc1[l])/2
			s = s + (wg[i] * p1 * ADisM(Z1,A1,Z2,A2,T,Ed,mdisp,bad,cad))

		if (NLa == 3):
			s = s + (wg[i] * 0.5 * ADisM(Z1,A1,Z2,A2,T,Ed,mdisp,bad,cad))

	return(s)

#=======Integrate over displacement model:: tabulated fmuE =======*
	
	# Average damage energy for discrete level reactions from tabulated
	# mu vs. f(mu,E) representation of secondary particle angular 
	# distribution.
	
def discrec2(fpr,E,Q,bta,Z1,A1,Z2,A2,Ed,mdisp,bad,cad):
	nquad = 64
	xabc = [0]*nquad; wg = [0]*nquad
	A = A2
	U = (A+1)*abs(Q)/A
	if (Q > 0):
		U = 0
	R = math.sqrt((A*(A+1-bta)/bta)*(1-U/E))
	R1 = bta*R/(A+1-bta)
	
	(xabc, wg) = GQ()	# Gauss-quadrature values and weights

	s = 0
	for i in range(nquad):
		T = E * (A+1-bta) * (1+R1*R1-2*R1*xabc[i]) / (A+1)**2
		s = s + (wg[i] * fpr[i] * ADisM(Z1,A1,Z2,A2,T,Ed,mdisp,bad,cad))
		
	return(s)
		
#=======Energy partition and displacement models=======*
		
	# Atom-displacement damage models: NRT and arc-dpa to find damage
	# energy for a given recoil energy in a given target material.
	
def ADisM(Z1,A1,Z2,A2,En,Ed,mdisp,bad,cad):
		
# z2 = 6			!14	! Z !				!z1, A1 = recoil nuclei
# A2 = 11.898 		!28.085	!			!z2, A2 = lattice nuclei
# Ed = 31.0!25.0	!42.0!24.0!20.0 !35.0 !25.0! 	!31.0 
		
	twothd = 0.666666667
	threeq = 0.75
	sixth = 0.166666667
	onep5 = 1.5
	c1 = 30.724
	c2 = 0.0793
	c3 = 3.4008
	c4 = 0.40244
	dont = 40.0
	el = c1*Z1*Z2*(math.sqrt(Z1**twothd+Z2**twothd))*(A1+A2)/A2
	rel = 1/el
	fl1 = c2*Z1**twothd*(math.sqrt(Z2))*(A1+A2)**onep5
	fl2 = ((Z1**twothd+Z2**twothd)**threeq*A1**onep5*(math.sqrt(A2)))
	fl = fl1/fl2
	ep = En*rel
	dam = En / (1 + fl * (c3*ep**sixth+c4*ep**threeq+ep))
		   
	df = dam
		
	if (Ed == 0):
		if (Z2==5):
			Ed = 40.0
		if (Z2==6):
			Ed = 31.0
		if (Z2==7):
			Ed = 40.0
		if (Z2==14):
			Ed = 25.0
		if (Z2==22):
			Ed = 40.0
		if (Z2==23):
			Ed = 40.0
		if (Z2==24):
			Ed = 40.0
		if (Z2==25):
			Ed = 40.0
		if (Z2==26):
			Ed = 40.0
		if (Z2==28):
			Ed = 40.0
		if (Z2==29):
			Ed = 40.0		# Cu
		if (Z2==39):
			Ed = 36.0		# Y
		if (Z2==40):
			Ed = 40.0		# Zr
		if (Z2==42):
			Ed = 60.0		# Mo
		if (Z2==46):
			Ed = 40.0		# Pd
		if (Z2==73):
			Ed = 90.0		# Ta
		if (Z2==74):
			Ed = 70.0		# W
		if (Z2==78):
			Ed = 40.0		# Pt
	if (Z1 == 0):
		df = 0.0
	if (En < Ed):
		df = 0.0
	if (Ed <= En and En < 2*Ed/0.8):
		df = 2*Ed/0.8
	if (En >= 2*Ed/0.8):
		if (mdisp == 1):
			df = dam
		if (mdisp == 2):
			if (bad == 0):
				if (Z2 == 5):
					bad = -1.0
				if (Z2 == 6):
					bad = -1.0
				if (Z2 == 7):
					bad = -1.0
				if (Z2 == 14):
					bad = -1.0
				if (Z2 == 22):
					bad = -1.0
				if (Z2 == 23):
					bad = -1.0
				if (Z2 == 24):
					bad = -1.0
				if (Z2 == 25):
					bad = -1.0
				if (Z2 == 26):
					bad = -0.568
				if (Z2 == 28):
					bad = -1.007
				if (Z2 == 29):
					bad = -0.68		# Cu
				if (Z2 == 39):
					bad = -1.0		# Y
				if (Z2 == 40):
					bad = -1.0		# Zr
				if (Z2 == 42):
					bad = -1.0		# Mo
				if (Z2 == 46):
					bad = -0.877		# Pd
				if (Z2 == 73):
					bad = -1.0		# Ta
				if (Z2 == 74):
					bad = -0.564		# W
				if (Z2 == 78):
					bad = -1.122		# Pt
			
			if (cad == 0):
				if (Z2 == 5):
					cad = 0.58
				if (Z2 == 6):
					cad = 0.71
				if (Z2 == 7):
					cad = 0.5 	# not mentioned in Konobeyev et al. paper Nucl.Eng.Tech. 3 (2017) 169-175
				if (Z2 == 14):
					cad = 0.5
				if (Z2 == 22):
					cad = 0.83
				if (Z2 == 23):
					cad = 0.51
				if (Z2 == 24):
					cad = 0.37
				if (Z2 == 25):
					cad = 0.33
				if (Z2 == 26):
					cad = 0.286
				if (Z2 == 28):
					cad = 0.227
				if (Z2 == 29):
					cad = 0.16
				if (Z2 == 39):
					cad = 0.5
				if (Z2 == 40):
					cad = 0.7
				if (Z2 == 42):
					cad = 0.46
				if (Z2 == 46):
					cad = 0.15
				if (Z2 == 73):
					cad = 0.72
				if (Z2 == 74):
					cad = 0.119
				if (Z2 == 78):
					cad = 0.11
			
			zi = (1-cad)*((0.8*dam)/(2*Ed))**bad + cad
			df = zi*dam
	
	return(df)
	
#=======Integrate over secondary energy transfer kernel :: 'al' coeffcients =======*
	
	# Average energy of recoil nucleus for heating from Legendre
	# expansion coefficient representation.
	
def discrec3(alc1,NLa,E,Q,bta,A2):
	nquad = 64
	Pl1 = [0]*NLa
	xabc = [0]*nquad; wg = [0]*nquad
	A = A2
	U = (A+1)*abs(Q)/A
	if (Q > 0):
		U = 0
	R = math.sqrt((A*(A+1-bta)/bta)*(1-U/E))
	R1 = bta*R/(A+1-bta)

	(xabc, wg) = GQ()	# Gauss-quadrature values and weights

	s = 0
	for i in range(64):
		Pl1[0] = 1
		Pl1[1] = xabc[i]
		T = E * (A+1-bta) * (1+R1*R1-2*R1*xabc[i]) / (A+1)**2
		if (NLa > 3):
			for j in range(2,NLa):
				Pl1[j] = (((2*(j-1)+1)*xabc[i]*Pl1[j-1])-((j-1)*Pl1[j-2])) / j
			p1 = 0
			for l in range(NLa):
				p1 = p1 + (((2*l)+1)*Pl1[l] * alc1[l])/2
			s = s + wg[i]*p1*T
		
		if (NLa == 3):
			s = s + (wg[i]*0.5*T)

	return(s)

#=================================================================

	# Average energy of charged particle for heating from Legendre 
	# expansion coefficient representation.
	
def discrec5(alc1,NLa,E,Q,bta,A2):
	nquad = 64
	Pl1 = [0]*NLa
	xabc = [0]*nquad; wg = [0]*nquad
	A = A2
	U = (A+1)*abs(Q)/A
	if (Q > 0):
		U = 0
	R = math.sqrt((A*(A+1-bta)/bta)*(1-U/E))
	R1 = bta*R/(A+1-bta)

	(xabc, wg) = GQ()	# Gauss-quadrature values and weights

	s = 0
	for i in range(64):
		Pl1[0] = 1
		Pl1[1] = xabc[i]
		T = E * bta * (1+R*R+2*R*xabc[i]) / (A+1)**2
		if (NLa > 3):
			for j in range(2,NLa):
				Pl1[j] = (((2*(j-1)+1)*xabc[i]*Pl1[j-1])-((j-1)*Pl1[j-2])) / j
			p1 = 0
			for l in range(NLa):
				p1 = p1 + (((2 * l)+1) * Pl1[l] * alc1[l]) / 2
			s = s + wg[i]*p1*T

		if (NLa == 3):
			s = s + (wg[i]*0.5*T)

	return(s)

#=======Integrate over secondary energy transfer kernel :: tabulated fmuE =======*

	# Average energy of recoil nucleus for discrete level reactions
	# from tabulated mu vs. f(mu,E) representation of secondary 
	# particle angular distribution.
 	
def discrec4(fpr,E,Q,bta,A2):
	nquad = 64	
	xabc = [0]*nquad; wg = [0]*nquad
	A = A2
	U = (A+1)*abs(Q)/A
	if (Q > 0):
		U = 0
	R = math.sqrt((A*(A+1-bta)/bta)*(1-U/E))
	R1 = bta*R/(A+1-bta)

	(xabc, wg) = GQ()	# Gauss-quadrature values and weights

	s = 0
	for i in range(nquad):
		T = E * (A+1-bta) * (1+R1*R1-2*R1*xabc[i]) / (A+1)**2
		s = s + wg[i]*fpr[i]*T
	return(s)

# ---------------------------------------
	
	# Average energy of charged particle for discrete level reactions
	# from tabulated mu vs. f(mu,E) representation of secondary 
	# particle angular distribution.
	
def discrec6(fpr,E,Q,bta,A2):
	nquad = 64
	xabc = [0]*nquad; wg = [0]*nquad
	A = A2
	U = (A+1)*abs(Q)/A
	if (Q > 0):
		U = 0
	R = math.sqrt((A*(A+1-bta)/bta)*(1-U/E))
	R1 = bta*R/(A+1-bta)

	(xabc, wg) = GQ()	# Gauss-quadrature values and weights

	s = 0
	for i in range(nquad):
		T = E * bta * (1+R*R+2*R*xabc[i]) / (A+1)**2
		s = s + wg[i]*fpr[i]*T
	return(s)

#=======Continuum (n,particle) when recoil data are present========*
		
	# Average damage energy from tabulated energy distributions in 
	# File 6.

def Tinteg1(z1,A1,z2,A2,Ed,ter6,tf6,n,mdisp,bad,cad):
	s = 0
	erl = 0
	fl = 0
	p1 = 0
	for i in range(n):
		eru = ter6[i]
		fu = tf6[i]
		p2 = ADisM(z1,A1,z2,A2,eru,Ed,mdisp,bad,cad)
		p = fu*p2
		ti = p2*fl + p1
		d = (eru-erl)*ti/2
		s = s + d
		erl = eru
		fl = fu
		p1 = p

	return(s)

#=======Calculate heating from continuum particle and recoil data=======*

	# Average recoil energy from tabulated energy distributions in 
	# File 6.

def Tinteg1heat(A2,ter6,tf6,n):
	s = 0
	erl = 0
	fl = 0
	p1 = 0
	for i in range(n):
		eru = ter6[i]
		fu = tf6[i]
		p2 = eru
		p = fu*p2
		ti = p2*fl + p1
		d = (eru-erl)*ti/2
		s = s + d
		erl = eru
		fl = fu
		p1 = p

	return(s)

#======= Integration =======*  		
    
	# Average damage energy from tabulated energy distributions

def conint(Z1,A1,Z2,A2,En,Ed,Ep1,fEEp,nfe,Nrc,mdisp,bad,cad):
	s = 0
	erl = 0
	fl = 0
	p1 = 0
	for i in range(nfe):
		eru = Ep1[i]
		et = (En+eru)/A1
		fu = fEEp[i]
		if (Nrc > 0):
			p2 = ADisM(Z1,A1,Z2,A2,eru,Ed,mdisp,bad,cad)
		if (Nrc == 0):
			p2 = eru # for using emitted neutron data
		p = fu*p2
		ti = p2*fl + p1
		d = (eru-erl)*ti/2
		s = s + d
		erl = eru
		fl = fu
		p1 = p

	return(s)

	# Average recoil energy from tabulated energy distributions	

def conintheat(Z1,A1,Z2,A2,En,Ep1,fEEp,nfe,Nrc):	#alc2d,NLa,xabc,wg,
	s = 0
	erl = 0
	fl = 0
	p1 = 0
	for i in range(nfe):
		eru = Ep1[i]
		et = (En+eru)/A1
		fu = fEEp[i]
		if (Nrc >= 0):
			p2 = eru
		p = fu*p2
		ti = p2*fl + p1
		d = (eru-erl)*ti/2
		s = s + d
		erl = eru
		fl = fu
		p1 = p

	return(s)

def angav(A,Z):
	nquad = 64
	wg = [0]*nquad; xabc = [0]*nquad; xabl = [0]*nquad
	(xabc,wg) = GQ()
	for i in range(64):
		xabl[i] = (1+(A*xabc[i]))/math.sqrt(1+2*A*xabc[i]+A*A)
	s = 0
	for i in range(64):
		s = s + (wg[i] * xabl[i])
	
	return(s)

	# Average damage energy from File 5 emitted neutron energy 
	# distributions
	
def conint1(Z1,A1,Z2,A2,En,Ed,Ep1,fEEp,nfe,mdisp,bad,cad):
	s = 0
	erl = 0
	fl = 0
	p1 = 0
	for i in range(nfe):
		eru = Ep1[i]
		fu = fEEp[i]
		p2 = dscrs3(Z1,A1,Z2,A2,En,Ed,eru,mdisp,bad,cad)
		p = fu*p2
		ti = p2*fl + p1
		d = (eru-erl)*ti/2
		s = s + d
		erl = eru
		fl = fu
		p1 = p

	return(s)

	# Average recoil energy from File 5 emitted neutron energy 
	# distributions

def conint1heat(Z1,A1,Z2,A2,En,Ep1,fEEp,nfe):
	s = 0
	erl = 0
	fl = 0
	p1 = 0
	for i in range(nfe):
		eru = Ep1[i]
		fu = fEEp[i]
		p2 = dscrs3heat(Z1,A1,Z2,A2,En,eru)
		p = fu*p2
		ti = p2*fl + p1
		d = (eru-erl)*ti/2
		s = s + d
		erl = eru
		fl = fu
		p1 = p

	return(s)
#===============================================	
	# Calculation of average daamge recoil energy from average neutron energy

def dscrs3(Z1,A1,Z2,A2,En,Ed,f,mdisp,bad,cad):
	nquad = 64
	wg = [0]*nquad; xabc = [0]*nquad; xabl = [0]*nquad
	(xabc,wg) = GQ()
	for i in range(nquad):
 		xabl[i] = (1+(A2*xabc[i]))/math.sqrt(1+2*A2*xabc[i]+A2*A2)
	s = 0
	for i in range(nquad):
		T = (En-2*math.sqrt(En*f)*xabl[i]+f)/A2
		s = s + (wg[i]*ADisM(Z1,A1,Z2,A2,T,Ed,mdisp,bad,cad))
	return(s)	
#===============================================
def dscrs3heat(Z1,A1,Z2,A2,En,f):
	nquad = 64
	wg = [0]*nquad; xabc = [0]*nquad; xabl = [0]*nquad
	(xabc,wg) = GQ()
	for i in range(nquad):
		xabl[i] = (1+(A2*xabc[i]))/math.sqrt(1+2*A2*xabc[i]+A2*A2)       
	s = 0
	for i in range(nquad):
		T = (En-2*math.sqrt(En*f)*xabl[i]+f)/A2
		s = s + (wg[i]*T)	
	return(s)
#===============================================

	# Calculation of average damage and recoil energy from evaporation 
	# spectrum data given for the emitted neutrons

def Tinteg3 (A,theta,Em1max,Z,Ed,En,mdisp,bad,cad):
	qx1 = [-0.86114,-0.33998,0.33998,0.86114]
	qx2 = [0.34785,0.65215,0.65215,0.34785]
	nq = 4
	IE = theta**2 * (1 - ((1 + Em1max/theta)*math.exp(-Em1max/theta)))
	ll = 0.0
	ul = Em1max
	Nsect = 100
	h = (ul-ll)/Nsect
	s = 0.0
	for j in range(Nsect):
		r = ll + (j-1)*h
		f11 = 0.
		for i in range(nq):
			T = (En-2*qx1[i]*math.sqrt(En*r)+r)/A
			f11 = f11 + qx2[i]*r*math.exp(-r/theta)*ADisM(Z,A,Z,A,T,Ed,mdisp,bad,cad)/2
		r1 = r+h
		f22 = 0.0
		for i in range(nq):
			T = (En-2*qx1[i]*math.sqrt(En*r1)+r1)/A
			f22 = f22 + qx2[i]*(r1)*math.exp(-(r1)/theta)*ADisM(Z,A,Z,A,T,Ed,mdisp,bad,cad)/2

		s = s + (h/2) * (f11 + f22)
	return(s/IE)

#===================================================
		
def Tinteg3heat(A,theta,Em1max,Z,En):
	qx1 = [-0.86114,-0.33998,0.33998,0.86114]
	qx2 = [0.34785,0.65215,0.65215,0.34785]
	nq = 4
	IE = theta**2 * (1 - ((1 + Em1max/theta)*math.exp(-Em1max/theta)))
	ll = 0.0
	ul = Em1max
	Nsect = 100
	h = (ul-ll)/Nsect
	s = 0.0
	for j in range(Nsect):
		r = ll + (j-1)*h
		f11 = 0.0
		for i in range(nq):
			T = (En-2*qx1[i]*math.sqrt(En*r)+r)/A
			f11 = f11 + qx2[i]*r*math.exp(-r/theta)*T/2

		r1 = r+h
		f22 = 0.0
		for i in range(nq):
			T = (En-2*qx1[i]*math.sqrt(En*r1)+r1)/A
			f22 = f22 + qx2[i]*(r1)*math.exp(-(r1)/theta)*T/2

		s = s + (h/2) * (f11 + f22)

	return(s/IE)

#=======Continuum (n,particle) when recoil data are absent========*
 	
	# Calculate average damage energy using one-particle recoil 
	# approximation	
	
def Tinteg(z1,A1,z2,A2,Ed,E,Q,bta,mdisp,bad,cad,n3,f1,f2,f3):
	nquad = 64
	xabc = [0]*nquad; wg = [0]*nquad
	(xabc,wg) = GQ()	
	s = 0
	for i in range(nquad):
		td = n3*(f1-(f2*xabc[i])+f3)
		s = s + (wg[i]*ADisM(z1,A1,z2,A2,td,Ed,mdisp,bad,cad))

# 		if (z2 == 5) s = s/1.5
# 		if (z2 == 6 .or. z2 == 14) s = s/2
# 		if (z2 == 22) s = s/2

# 		A = A2
# 		U = (A+1)*abs(Q)/A
# 		if (Q>0) U = 0
# 		R = sqrt((A*(A+1-bta)/bta)*(1-U/E))
# 		R1 = bta*R/(A+1-bta)
# 		s = 0
# 		!if (Q>0 .and. E<Q) then
# 		!	s = 0
# 		!else
# 		!s = 0
# 		do i = 1, 64
# 			T = E*(A+1-bta)*(1+R1*R1-2*R1*xabc(i))/(A+1)**2
# 			s = s + (wg(i)*0.5*T)!ADisM(Z1,A1,Z2,A2,T,Ed,mdisp,bad,cad))
# 		end do
# 		!end if
	return(s)

#=======Interpolation Schemes=======*
        
# Definitions of different interpolation schemes
def TERPOLIN (intflg,x,x1,x2,y1,y2):
	y = y1
	if (intflg==1 or intflg==11 or intflg==21):
		y = y1
	if ((intflg==2 or intflg==12 or intflg==22) and (x2-x1) != 0):
		y = y1+(y2-y1)*(x-x1)/(x2-x1)
	if (intflg==3 or intflg==13 or intflg==23):
		y = y1+(y2-y1)*(log(x/x1))/log(x2/x1)
	if ((intflg==4 or intflg==14 or intflg==24) and (x2-x1) != 0):
		y = y1*exp((x-x1)*log(y2/y1)/(x2-x1))
	if (intflg==5 or intflg==15 or intflg==25):
		y = y1*exp(log(x/x1)*log(y2/y1)/log(x2/x1))
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

#=======Interpolate yield data in required energy=======*

	# To interpolate the yields of particles from the File 6 energy
	# points to the energy points in File 3.
	
def yldterpolin (Nrule,irule,n1,E1,Y1,n2):
	E1 = [0]*n1; Y1 = [0]*n1
	E2 = [0]*n2; Y2 = [0]*n2
	Nrule = [0]*20; irule = [0]*20
	for i in range(n2):
		for j in range(n1):
			if (E2[i] == E1[j]):
				Y2[i] = Y1[j]
				break
			if (E1[j]<E2[i] and E2[i]<E1[j+1]):
				diff1 = E2[i] - E1[j]
				diff2 = E2[i] - E2[i-1]
				for k in range(20):
					if (j <= Nrule[k]):
						intflg = irule(k)
						break
				if (diff1 <= diff2):
					Y2[i] = TERPOLIN(intflg,E2[i],E1[j],E1[j+1],Y1[j],Y1[j+1])
				if (diff2 < diff1):
					Y2[i] = TERPOLIN(intflg,E2[i],E2[i-1],E1[j+1],Y2[i-1],Y1[j+1])
				break
	return(E2,Y2)

#=================================================	
	# Calculate average recoil energy using one-particle recoil 
	# approximation to estimate heating
	
def Tintegheat1(A2,E,Q,bta,n3,f1,f2,f3):
	nquad = 64
	xabc = [0]*nquad; wg = [0]*nquad
	(xabc,wg) = GQ()	
	s = 0
	for i in range(nquad):
		td = n3*(f1-(f2*xabc[i])+f3)
		s = s + (wg[i]*td)

# 		if (z2 == 5) s = s/1.5
# 		if (z2 == 6 .or. z2 == 14) s = s/2
# 		if (z2 == 22) s = s/2
		
# 		A = A2
# 		U = (A+1)*abs(Q)/A
# 		if (Q>0) U = 0
# 		R = sqrt((A*(A+1-bta)/bta)*(1-U/E))
# 		R1 = bta*R/(A+1-bta)
# 		s = 0
#		if (Q>0 .and. E<Q) then
#			s = 0
#		else
#		s = 0
# 		do i = 1, 64
# 			T = E*(A+1-bta)*(1+R1*R1-2*R1*xabc(i))/(A+1)**2
# 			s = s + (wg(i)*T)
# 		end do
#		end if
	return(s)	#*2

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
		
	ifile = open('tape01', 'r')
	
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

	ifile.close()

	return (NBTp,INTrp,Nyld,Eyld,Yld)

#=======Other Nonelastic=======*
	
	## Calculation of neutron dpa and heating cross sections due to
	## charged particle out reactions of neutrons.
	
def n_CPO (ofile_outRMINDD,MTi,lpr,mdisp,Ed,bad,cad,NPt,Etu):
	sdpat = numpy.zeros(NPt); snhtt = numpy.zeros(NPt)
	signcpol = numpy.zeros(NPt); tot_energy_products1 = numpy.zeros(NPt); num_of_displ1 = numpy.zeros(NPt)
	# for energy-balance heating
	snhtt_EB = numpy.zeros(NPt)
	tot_energy_n_photons1 = numpy.zeros(NPt)
	# ------
	beta = numpy.zeros(281); n1 = numpy.zeros(281); ze = numpy.zeros(281)
	cbe = numpy.zeros(281)
	MTdthl = numpy.zeros(250)
	xabc = numpy.zeros(64); wg = numpy.zeros(64)

	# for file6
	fiso = [0]*5; LAW = [0]*5; Nyld = [0]*5
	LG = [0]*5; NE6 = [0]*5
	NEP = numpy.zeros((5,200)); NL = numpy.zeros((5,200)); NBT = numpy.zeros((5,20))
	INTr = numpy.zeros((5,20)); ND6g = [0]*400
	Nyld = [0]*5
	Eint6 = numpy.zeros((5,400))
	yi = numpy.zeros((5,400))
	En = numpy.zeros((5,400))
	Enp = numpy.zeros((5,400,1000))
	f = numpy.zeros((5,400,1000))		  	# dimension increased from 1000 to 10000 for W 180
	al6 = numpy.zeros((5,400,65))
	cdata6 = numpy.zeros((5,400,201)); fdata6 = numpy.zeros((5,400,201))
	Ball = [0]*10000				# dimension increased from 1000 to 10000 for W 180
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
#--------------------------------------------------------------------

	print('', file = ofile_outRMINDD)
	print('CPO reaction MT=',MTi,' dpa/heating cross section', file = ofile_outRMINDD)
	print('-------------------------------------------------', file = ofile_outRMINDD)
#--------------------------------------------------------------------

	iflpresent = 0 		# flag for the presence of cross sections MF=3. 

	ifile = open ('tape02', 'r')

	## extraction of cross sections
	while True:
		line = ifile.readline()
		if (line == ''):
			break
		data = eachlineinfo(line)
		MAT = int(data[6]); MF = int(data[7]); MT  = int(data[8])
		if ( MT == 0 and MF != 0):
			line = ifile.readline()
			data = eachlineinfo(line)
			iflspace = 0
			for element in data:
				if (element == ''):
					iflspace = 1
					break
			if (iflspace == 0):
				ZAv = float(data[0]); AWRv = float(data[1]); L0 = int(data[2]);\
				L1 = int(data[3]); NKv = int(data[4]); L2 = int(data[5]);\
				MAT = int(data[6]); MF = int(data[7]); MT = int(data[8])

		if (MAT != -1):
			if (MF == 3):
				if (MT == MTi):
					iflpresent = 1
					ifdthl = 0
					for i in range (250):
						if (MT == int(MTdthl[i])):
							ifdthl = 1
							break
					ZA = ZAv
					AWR = AWRv
					line = ifile.readline()
					data = eachlineinfo(line)
					QM = float(data[0]); QI =  float(data[1]); LR = int(data[3])
					NR = int(data[4]); NP = int(data[5])
					ifile.readline()
					if (NP > NPt):
						Eall = [0]*NPt; sall = [0]*NPt; sdall = [0]*NPt
						shall = [0]*NPt; Eall1 = [0]*NP; sall1 = [0]*NP
						num_of_displ = [0]*NPt; tot_energy_products = [0]*NPt
						shall_EB = [0]*NPt
						tot_energy_n_photons = [0]*NPt
						dn4_neutrons = [0]*NPt; dn5_photons = [0]*NPt
						(Eall1, sall1) = line_type3_info(ifile,NP,2)
					if (NP <= NPt):
						Eall = [0]*NP; sall = [0]*NP; sdall = [0]*NP
						shall = [0]*NP; num_of_displ = [0]*NP; tot_energy_products = [0]*NP
						shall_EB = [0]*NP
						tot_energy_n_photons = [0]*NP
						dn4_neutrons = [0]*NP; dn5_photons = [0]*NP
						(Eall, sall) = line_type3_info(ifile,NP,2)
					if (ifdthl == 1):
						alfull = [[0]*65]*NP; alc = [0]*65
		else:
			break

	ifile.close()

	print('', file = ofile_outRMINDD)
	print('Number of energy points given is ',NP, file = ofile_outRMINDD)

	Z = ZA//1000
	A = AWR

## only if MF3 cross sections are available then do the following

	if (iflpresent == 1):

## extraction of Legendre polynomial coefficients and Probability
## from File 4

		ifl4 = 0
		ifspad4 = 0 #flag for secondary particle angular data
		ifspad4al = 0 # flag for secondary particle angular data in 'al' coefficients
		ifspad4muf = 0 # flag for secondary particle angular data in 'mu,f' form

		if (MTi >= int(MTdthl[0])):
			ifile = open ('tape01', 'r')
			while True:
				line = ifile.readline()
				if (line == ''):
					break
				data = eachlineinfo(line)
				MAT = int(data[6]); MF = int(data[7]); MT  = int(data[8])
				if (MF==3 and MT==0 ):
					line = ifile.readline()
					(ZA,AWR,l1,LTT,NK,l2,MAT,MF,MT) = line_type1_info(line)
				if (MAT != -1):
					if (MF == 4):
						if (MT == MTi):
							ifl4 = 1; ifspad4 = 1
							ifile.readline()
							if (LTT == 3 or LTT == 1):
								ifspad4al = 1
								line = ifile.readline()
								data = eachlineinfo(line)
								NE1 = int(data[5])
								## Legendre polynomial coefficients
								ifile.readline()
								al4 = [[0]*65]*NE1; EL = [0]*NE1
								for i in range (NE1):
									line = ifile.readline()
									data = eachlineinfo(line)
									T = float(data[0]); EL[i] = float(data[1]); NL4 = int(data[4])
									NL4 = NL4 + 1
									al4[i][0] = 1
									al4[i][1:-1] = 0.0
									temporary = [0]*NL4
									temporary = line_type3_info (ifile,NL4-1,1)
									for j, value in enumerate (temporary, 1):
										al4[i][j] = value

							if (LTT == 3 or LTT == 2):
								ifspad4muf = 1
								NE2 = int(ifile.readline().split()[5])
								## Probability
								ifile.readline()
								Enf = [0]*NE2; NPr = [0]*NE2; cdata4 = numpy.zeros((NE2,201))
								fdata4 = numpy.zeros((NE2,201)); ftotal= numpy.zeros((NE2,64))
								fmuE = numpy.zeros((NP,64)); fpr = [0]*64
								for i in range (NE2):
									line = ifile.readline()
									data = eachlineinfo(line)
									T = float(data[0]); Enf[i] = float(data[1]) 
									NPr[i] = int(ifile.readline().split()[0])
									temporary1 = [0]*NPr[i]
									temporary2 = [0]*NPr[i]
									(temporary1,temporary2) = line_type3_info(ifile,NPr[i],2)
									for j, (value1, value2) in enumerate(zip(temporary1, temporary2), 0):
										cdata4[i][j] = value1
										fdata4[i][j] = value2
				else:
					break
			ifile.close()

		# Data from File 6

		ifl6 = 0
		iflsp = 0 #flag to know presence of secondary particle data in continuum form in MF = 6
		iflr = 0  #flag to know presence of recoil data in MF = 6
		iflawiso = 0 #flag for isotropic distribution
		ifspad6 = 0 #flag for secondary particle angular data
		iffdlcd = 0 #flag for fdata linear with cdata (LAW=2)
		iflgfdlcd = 0 #flag for log(fdata) linear with cdata (LAW=2)
		ifllab = 0 	# secondary energy and angle in lab. system
		iflcom = 0 # secondary energy and angle in c.o.m. system
		iflcomlab = 0 # secondary energy and angle in c.o.m (A<=4), lab. (A>4)

		ifile = open ('tape01', 'r')

		while True:
			line = ifile.readline()
			if (line == ''):
				break
			data = eachlineinfo(line)
			MAT = int(data[6]); MF = int(data[7]); MT  = int(data[8])
			if (MT == 0):
				line = ifile.readline()
				data = eachlineinfo(line)
				iflspace = 0
				for element in data:
					if (element == ''):
						iflspace = 1
						break
				if (iflspace == 0):
					(ZAv,AWRv,l1,LCT,NKv,l2,MAT,MF,MT) = line_type1_info(line)
			if (MAT != -1):
				if (MF == 6):
					if (MT == MTi):
						ifl6 = 1
						NK = NKv
						if (LCT == 1):
							ifllab = 1
						if (LCT == 2):
							iflcom = 1
						if (LCT == 3):
							iflcomlab = 1
						for NSS in range (NK):
							line = ifile.readline()
							data = eachlineinfo(line)
							(ZAP[NSS],AWP[NSS],LIP,LAW[NSS],NR6,Nyld[NSS],MAT,MF,MT) = line_type1_info(line)

							if (float(AWP[NSS]) > int(beta[lpr])):
								irs = NSS
								iflr = 1
							temporary1 = [0]*NR6
							temporary2 = [0]*NR6
							(temporary1,temporary2) = line_type3_info(ifile,NR6,2)
							for N, (value1, value2) in enumerate(zip(temporary1, temporary2), 0):
								NBT[NSS][N] = value1
								INTr[NSS][N] = value2
							temporary1 = [0]*int(Nyld[NSS])
							temporary2 = [0]*int(Nyld[NSS])
							(temporary1,temporary2) = line_type3_info(ifile,int(Nyld[NSS]),2)
							for N, (value1, value2) in enumerate(zip(temporary1, temporary2), 0):
								Eint6[NSS][N] = value1
								yi[NSS][N] = value2

							if (LAW[NSS] == 3):
								iflawiso = 1
							if (LAW[NSS] == 2):
								if(int(math.ceil(AWP[NSS])) == beta[lpr]):
									ifspad6 = 1
									isps = NSS
								line = ifile.readline()
								data = eachlineinfo(line)
								c1 = float(data[0]); c2 = float(data[1]); l3 = int(data[2])
								l4 = int(data[3]); NR6 = int(data[4]); NE6[NSS] = int(data[5])
								MAT = int(data[6]); MF = int(data[7]); MT = int(data[8])
								temporary1 = [0]*NR6
								temporary2 = [0]*NR6
								(temporary1,temporary2) = line_type3_info(ifile,NR6,2)
								for N, (value1, value2) in enumerate(zip(temporary1, temporary2), 0):
									NBT[NSS][N] = value1
									INTr[NSS][N] = value2

								for i in range (int(NE6[NSS])):
									line = ifile.readline()
									data = eachlineinfo(line)
									c1 = float(data[0]); En[NSS][i] = float(data[1]); LG[NSS] = int(data[2])
									l2 = int(data[3]); NW = int(data[4]); NL[NSS][i] = int(data[5])
									MAT = int(data[6]); MF = int(data[7]); MT = int(data[8])

									if (int(LG[NSS]) == 0):
										NL[NSS][i] = NL[NSS][i] + 1
										al6[NSS][i][0] = 1
										al6[NSS][i][1:-1] = 0
										temporary = [0]*int(NL[NSS][i])
										temporary = line_type3_info (ifile,int(NL[NSS][i])-1,1)
										for j, value in enumerate (temporary, 1):
											al6[NSS][i][j] = value

									if (int(LG[NSS]) > 0):
										if (int(LG[NSS]) == 12):
											iffdlcd = 1
										if (int(LG[NSS]) == 14):
											iflgfdlcd = 1
										temporary1 = [0]*int(NL[NSS][i])
										temporary2 = [0]*int(NL[NSS][i])
										(temporary1,temporary2) = line_type3_info(ifile,int(NL[NSS][i]),2)
										for j, (value1, value2) in enumerate(zip(temporary1, temporary2), 0):
											cdata6[NSS][i][j] = value1
											fdata6[NSS][i][j] = value2

							if (int(LAW[NSS]) == 1):
								line = ifile.readline()
								data = eachlineinfo(line)
								c1 = float(data[0]); c2 = float(data[1]); LG[NSS] = int(data[2])
								LP = int(data[3]); NR6 = int(data[4]); NE6[NSS] = int(data[5])
								MAT = int(data[6]); MF = int(data[7]); MT = int(data[8])
								temporary1 = [0]*NR6
								temporary2 = [0]*NR6
								(temporary1,temporary2) = line_type3_info(ifile,NR6,2)
								for N, (value1, value2) in enumerate(zip(temporary1, temporary2), 0):
									NBT[NSS][N] = value1
									INTr[NSS][N] = value2
								for i in range (int(NE6[NSS])):
									line = ifile.readline()
									data = eachlineinfo(line)
									c1 = float(data[0]); En[NSS][i] = float(data[1]); ND = int(data[2])
									NA = int(data[3]) 
									NW = int(data[4])
									NEP[NSS][i] = int(data[5])
									MAT = int(data[6]); MF = int(data[7]); MT = int(data[8])

									if (NSS == NK):
										ND6g[i] = ND
									if (NA != 0):
										if (int(LG[NSS]) == 2 or int(LG[NSS]) == 1):
											fiso[NSS] = 1
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
												K1 = K1 + NA + 2
												K2 = K2 + NA + 2
									if (int(LG[NSS]) == 1 and NA == 0):
										fiso[NSS] = 1
										temporary1 = [0]*int(NEP[NSS][i])
										temporary2 = [0]*int(NEP[NSS][i])
										(temporary1,temporary2) = line_type3_info(ifile,int(NEP[NSS][i]),2)
										for j, (value1, value2) in enumerate(zip(temporary1, temporary2), 0):
											Enp[NSS][i][j] = value1
											f[NSS][i][j] = value2
			else:
				break

		ifile.close()

		ifllargeNP = 0
		if (NP > NPt):
			for i in range (NPt):
				Eall[i] = Etu[i]
				for j in range (NP):
					if (Eall[i] == Eall1[j]):
						sall[i] = sall1[j]
						break
					if (Eall1[j] < Eall[i] and Eall[i] < Eall1[j+1]):
						sall[i] = crstd(Eall[i],Eall1[j],Eall1[j+1],sall1[j],sall1[j+1])
						break
			ifllargeNP = 1
			NP = NPt

	# n1 = (A+1-A')/(A+1)
	# cbe = Coulomb barrier energy

		n1[lpr] = (A+1-beta[lpr])/(A+1)
		cbe[lpr] = (1.029e+6 * ze[lpr]*Z) / (beta[lpr]**(1/3) + A**(1/3))

		n2 = A/(A+1)
		n3 = 1.0/(A+1)

		if (LR == 0 or QI == 0):
			QI = QM
			if (QI < 0):
				U = abs(QI/n2)

		nquad = 64
		xabc = [0]*nquad; wg = [0]*nquad
		(xabc, wg) = GQ()

		if (ifl4 == 1):
			print('', file = ofile_outRMINDD)
			print('Emitted charged particle data are given in File 4', file = ofile_outRMINDD)

			if (LTT==3 or LTT==2):

	## LOG- LINEAR INTERPOLATION BETWEEN MU AND F(MU,E) 

				for i in range (NE2):
					for j in range (64):
						for k in range (int(NPr[i])):
							if (xabc[j] == cdata4[i,k]):
								ftotal[i][j] = fdata4[i][k]
								break
							if (xabc[j]>cdata4[i,k] and xabc[j]<cdata4[i][k+1]):
								x = xabc[j]
								x1 = cdata4[i][k]
								x2 = cdata4[i][k+1]
								y1 = fdata4[i][k]
								y2 = fdata4[i][k+1]
								ftotal[j][k] = y1*math.exp((x-x1)*math.log(y2/y1)/(x2-x1))   
								break

	## TO GET THE F(MU,E) FOR THE FULL ENERGY RANGE

				for i in range (NP):
					for j in range (64):
						fmuE[i][j] = 0.5

				for i in range (NP):
					for j in range (NE2):
						if (Eall[i] == Enf[j]):
							for k in range(64):
								fmuE[i][k] = ftotal[j][k]
							break
						if (Eall[i]>Enf[j] and Eall[i]<Enf[j+1]):
							x = Eall[i]
							x1 = Enf[j]
							x2 = Enf[j+1]
							for k in range (64):
								y1 = ftotal[j][k]
								y2 = ftotal[j+1][k]
								fmuE[i][k] = crstd(x,x1,x2,y1,y2)
							break

			if (LTT == 0):
				for i in range (NP):
					alfull[i][0] = 1
					for k in range (1, 65):
						alfull[i][k] = 0

			if (LTT==3 or LTT==1):
				for i in range (NP):
					alfull[i][0] = 1
					for k in range (1, 65):
						alfull[i][k] = 0

				for i in range (NP):
					if (sall[i] != 0):
						for j in range (NE1-1):
							if (Eall[i] == EL[j]):
								for k in range (1, 65):
									alfull[i][k] = al4[j][k]
								break
							if (EL[j] < Eall[i] and Eall[i] <= EL[j+1]):
								diff1 = Eall[i] - EL[j]
								diff2 = Eall[i] - Eall[i-1]
								if (diff1 <= diff2):
									for k in range (1, 65):
										x = Eall[i]
										x1 = EL[j]
										x2 = EL[j+1]
										y1 = al4[j][k]
										y2 = al4[j+1][k]
										alfull[i][k] = crstd(x,x1,x2,y1,y2)
								if (diff2 < diff1):
									for k in range (1, 65):
										x = Eall[i]
										x1 = Eall[i-1]
										x2 = EL[j+1]
										y1 = alfull[i-1][k]
										y2 = al4[j+1][k]
										alfull[i][k] = crstd(x,x1,x2,y1,y2)
								for k in range (1, 65):
									if (alfull[i][k] == 0):
										alfull[i][k] = alfull[i-1][k]
								break

	# discrec1-- average damage energy from Legendre Polynomial 
	# expansion Coeffiecients
	# discrec3-- average energy of recoil nucleus from " " " "
	# discrec5-- average energy of charged particle from " " " " 
	# discrec2-- average damage energy from tabulated mu vs. f(mu,E)
	# discrec4-- average energy of recoil nucleus from " " " "
	# discrec6-- average energy of charged particle from " " " " 

			NLa = 65
			if (LTT==0):
				NLa = 3
			for i in range (NP):
				if (sall[i] != 0):		#.and. Eall(i)>=abs(QI)
					if((LTT==3 and Eall[i]<=EL[NE1]) or LTT==1 or LTT==0):
						for k in range (NLa):
							alc[k] = alfull[i][k] 
						dn1 = abs(discrec1(alc,NLa,Eall[i],QI,beta[lpr],Z-ze[lpr], \
										A+1-beta[lpr],Z,A,Ed,mdisp,bad,cad))

						dn2 = abs(discrec3(alc,NLa,Eall[i],QI,beta[lpr],A))	
						dn3 = abs(discrec5(alc,NLa,Eall[i],QI,beta[lpr],A))

					if((LTT == 3 and Eall[i]>Enf[0]) or LTT == 2):
						for j in range (64):
							fpr[j] = fmuE[i][j]
						dn1 = abs(discrec2(fpr,Eall[i],QI,beta[lpr],Z-ze[lpr], \
									A+1-beta[lpr],Z,A,Ed,mdisp,bad,cad))

						dn2 = abs(discrec4(fpr,Eall[i],QI,beta[lpr],A))
						dn3 = abs(discrec6(fpr,Eall[i],QI,beta[lpr],A))

					if (sall[i] == 0):
						sdall[i] = 0
						shall[i] = 0
					if (sall[i] != 0):
						dn1 = dn1*0.8/(2*Ed)
						num_of_displ[i] = dn1
						tot_energy_products[i] = dn2 + dn3
						sdall[i] = sall[i]*dn1
						shall[i] = sall[i]*(dn2+dn3)
					if (sdall[i] < 1e-12):
						sdall[i] = 0
					if (shall[i] < 1e-12):
						shall[i] = 0


		if (ifl6 == 1):
			if (iflawiso == 1):
				for i in range (NP):
					alfull[i][0] = 1
					for k in range (1, 66):
						alfull[i][k] = 0

			if (ifspad6 == 1):
				print('', file = ofile_outRMINDD)
				print('Emitted charged particle data are given in File 6', file = ofile_outRMINDD)

				if (LG[isps] == 0):
					for i in range (NP):
						alfull[i][0] = 1
						for k in range (1, 65):
							alfull[i][k] = 0
					for i in range (NP):
						if (sall[i] != 0):
							for j in range (int(NE6[isps])-1):
								if (Eall[i] ==  En[isps][j]):
									for k in range (1, 65):
										alfull[i][k] = al6[isps][j][k]
									break
								if (En[isps][j] < Eall[i] and Eall[i] <= En[isps][j+1]):
									diff1 = Eall[i] - En[isps][j]
									diff2 = Eall[i] - Eall[i-1]
									if (diff1 <= diff2):
										for k in range (1, 65):
											x = Eall[i]
											x1 = En[isps][j]
											x2 = En[isps][j+1]
											y1 = al6[isps][j][k]
											y2 = al6[isps][j+1][k]
											alfull[i][k] = crstd(x,x1,x2,y1,y2)
									if (diff2 < diff1):
										for k in range (1, 65):
											x = Eall[i]
											x1 = Eall[i-1]
											x2 = En[isps][j+1]
											y1 = alfull[i-1][k]
											y2 = al6[isps][j+1][k]
											alfull[i][k] = crstd(x,x1,x2,y1,y2)
									for k in range (1, 65):
										if (alfull[i][k] == 0):
											alfull[i][k] = alfull[i-1][k]
									break

				if (LG[isps]>0):

	 ## INTERPOLATION BETWEEN MU AND F(MU,E) 

					for i in range (int(NE6[isps])):
						for j in range (64):
							for k in range (int(NL[isps][i])):
								if (xabc[j] == cdata6[isps][i][k]):
									ftotal[i][j] = fdata6[isps][i][k]
									break
								if(xabc[j]>cdata6[isps][i][k] and xabc[j]<cdata6[isps][i][k+1]):
									x = xabc[j]
									x1 = cdata6[isps][i][k]
									x2 = cdata6[isps][i][k+1]
									y1 = fdata6[isps][i][k]
									y2 = fdata6[isps][i][k+1]
									if (LG[isps] == 14):
										ftotal[i][j] = y1*math.exp((x-x1)*math.log(y2/y1)/(x2-x1))
									if(LG[isps] == 12):
										ftotal[i][j] = crstd(x,x1,x2,y1,y2)
									break

		## TO GET THE F(MU,E) FOR THE FULL ENERGY RANGE

					for i in range (NP):
						for j in range (64):
							fmuE[j][i] = 0.5
					for i in range (NP):
						for j in range (int(NE6[isps])-1):
							if (Eall[i] == En[isps][j]):
								for k in range (64):
									fmuE[i][k] = ftotal[j][k]
								break
							if (Eall[i] > En[isps][j] and Eall[i] <= En[isps][j+1]):
								x = Eall[i]
								x1 = En[isps][j]
								x2 = En[isps][j+1]
								for k in range (64):
									y1 = ftotal[j][k]
									y2 = ftotal[j+1][k]
									fmuE[i][k] = crstd(x,x1,x2,y1,y2)
								break

				NLa = 65
				if (iflawiso == 1):
					NLa = 3
				for i in range (NP):
					if (sall[i] != 0):	 #  .and. Eall(i)>=abs(QI)
						if (LG[isps] == 0 or iflawiso == 1):
							for k in range (NLa):
								alc[k] = alfull[i][k] 
							dn1 = abs(discrec1(alc,NLa,Eall[i],QI,beta[lpr],Z-ze[lpr],\
											A+1-beta[lpr],Z,A,Ed,mdisp,bad,cad))
							dn2 = abs(discrec3(alc,NLa,Eall[i],QI,beta[lpr],A))
							dn3 = abs(discrec5(alc,NLa,Eall[i],QI,beta[lpr],A))
						if (LG[isps] > 0):			
							for j in range (64):
								fpr[j] = fmuE[i][j]
							dn1 = abs(discrec2(fpr,Eall[i],QI,beta[lpr],Z-ze[lpr],\
											A+1-beta[lpr],Z,A,Ed,mdisp,bad,cad))	
							dn2 = abs(discrec4(fpr,Eall[i],QI,beta[lpr],A))
							dn3 = abs(discrec6(fpr,Eall[i],QI,beta[lpr],A))

						if (sall[i] == 0):
							sdall[i] = 0
							shall[i] = 0
						if (sall[i] != 0):
							dn1 = dn1*0.8/(2*Ed)
							num_of_displ[i] = dn1
							tot_energy_products[i] = dn2 + dn3
							sdall[i] = sall[i]*dn1
							shall[i] = sall[i]*(dn2+dn3)
						if (sdall[i]<1e-12):
							sdall[i] = 0
						if (shall[i]<1e-12):
							shall[i] = 0

	# When the energy distribution of the recoil nucleus is present

			if (iflr==1 and ifspad6==0):
				print('', file = ofile_outRMINDD)
				print('Recoil nucleus and charged particle energy', file = ofile_outRMINDD)
				print('distribution data are given in File 6', file = ofile_outRMINDD)

				kdim = int(numpy.amax(NEP))
				ter6 = [0]*kdim; tf6 = [0]*kdim; ntm = [0]*NP

				for i in range (NP):
					for j in range (int(NE6[irs])):
						if (En[irs][j] == Eall[i]):
							for k in range (int(NEP[irs][j])-1):
								ter6[k] = Enp[irs][j][k]
								tf6[k] = f[irs][j][k]
							ntm[i] = NEP[irs][j]
							break
						if(En[irs][j] < Eall[i] and Eall[i] <= En[irs][j+1]):
							for k in range (int(NEP[irs][j])):	#TENDL-2017 Ni59 --> NEP(irs,j+1)
								x = Eall[i]
								x1 = En[irs][j]	
								x2 = En[irs][j+1]
								y1 = Enp[irs][j][k]
								y2 = Enp[irs][j+1][k]
								y11 = f[irs][j][k]
								y22 = f[irs][j+1][k]
								ter6[k] = crstd(x,x1,x2,y1,y2)
								tf6[k] = crstd(x,x1,x2,y11,y22)

							ntm[i] = int(NEP[irs][j]) #TENDL-2017 Ni59 --> NEP(irs,j+1)
							break

					#do i = 1, NP
					#	ter6 = 0
					#	tf6 = 0
					for k in range (int(ntm[i])):
						#ter6(k) = Er6(i,k)
						#tf6(k) = fr6(i,k)
						if (ter6[k] == 1.0):
							ter6[k] = 0
						if (tf6[k] >= 1.0):
							tf6[k] = 0
					dn1 = 0
					dn2 = 0
					if (sall[i] != 0):	 #  .and. Eall(i)>=abs(QI)
						dn1 = Tinteg1(Z-ze[lpr],A+1-beta[lpr],Z,A,Ed,ter6,tf6,int(ntm[i]),mdisp,bad,cad)
						dn2 = Tinteg1heat(A,ter6,tf6,int(ntm[i]))

					dn1 = abs(dn1)
					dn2 = abs(dn2)

					if (sall[i] == 0):
						sdall[i] = 0
						shall[i] = 0
					if (sall[i] != 0):
						dn1 = dn1*0.8/(2*Ed)
						num_of_displ[i] = dn1
						tot_energy_products[i] = dn2
						sdall[i] = sall[i]*dn1
						shall[i] = sall[i]*dn2
					if (sdall[i] < 1e-12):
						sdall[i] = 0
					if (shall[i] < 1e-12):
						shall[i] = 0

			## THESE PART (BELOW) IS FOR THE ENERGY DEPOSITED BY THE 
			## LIGHT CHARGED PARTICLE DATA PRODUCED IN NEUTRON REACTION, IT
			## IS DONE ONLY IF BOTH CONTINUUM RECOIL AND PARTICLE SECTIONS
			## ARE GIVEN. IT CORRESPONDS TO ONLY HEATING CALCULATIONS, NOT
			## IN THE DISPLACEMENT CROSS SECTIONS.
			## THE ENERGY DEPOSITED BY THE CHARGED PARTICLES WILL BE ADDED 
			## TO THAT DEPOSITED BY THE RECOIL NUCLEUS.

				ipstart = 0

				if (MTi==11 or MTi==22 or MTi==23 or MTi==24 or MTi==25 or
				MTi==28 or MTi==29 or MTi==30 or MTi==32 or MTi==33 or
				MTi==34 or MTi==35 or MTi==36 or MTi==41 or MTi==42 or MTi==44 or MTi==45):

					ipstart = 1

				for iparticle in range (ipstart, irs):
					ntm = [0]*NP
					for i in range (NP):
						ter6 = [0]*kdim
						tf6 = [0]*kdim
						for j in range (1, int(NE6[iparticle])-1):
							if (En[iparticle][j] == Eall[i]):
								for k in range (int(NEP[iparticle][j])):
									ter6[k] = Enp[iparticle][j][k]
									tf6[k] = f[iparticle][j][k]
								ntm[i] = int(NEP[iparticle][j])
								break

							if (En[iparticle][j] < Eall[i] and Eall[i] <= En[iparticle][j+1]):
								for k in range (int(NEP[iparticle][j])):
									x = Eall[i]
									x1 = En[iparticle][j]
									x2 = En[iparticle][j+1]
									y1 = Enp[iparticle][j][k]
									y2 = Enp[iparticle][j+1][k]
									y11 = f[iparticle][j][k]
									y22 = f[iparticle][j+1][k]
									ter6[k] = crstd(x,x1,x2,y1,y2)
									tf6[k] = crstd(x,x1,x2,y11,y22)

								ntm[i] = int(NEP[iparticle][j])
								break

						for k in range (int(ntm[i])):
							#if(iflcomlab == 1):
							#	ter6[k] = ter6[k] + (AWP[iparticle]/(A+1)**2)*Eall[i]
							#if(iflcom==1):
							#	ter6[k] = ter6[k] + (AWP[iparticle]/(A+1)**2)*Eall[i]
							if (ter6[k] == 1.0):
								ter6[k] = 0
							if (tf6[k]>=1.0):
								tf6[k] = 0
						dn3 = 0
						if (sall[i] != 0):		#  .and. Eall(i)>=abs(QI)
							dn3 = yi[iparticle][0]*Tinteg1heat(A,ter6,tf6,ntm[i])

						dn3 = abs(dn3)

						if (sall[i] == 0):
							shall[i] = 0
						if (sall[i] != 0):
							tot_energy_products[i] = tot_energy_products[i] + dn3
							shall[i] = shall[i] + sall[i]*dn3
						if (shall[i] < 1e-12):
							shall[i] = 0

			if (iflr == 0 and ifspad6 == 0):	#.or. ifl6==1

				iflnrcyp = 0 
				ipstart = 0

				if (MTi==11 or MTi==22 or MTi==23 or MTi==24 or MTi==25 or
				MTi==28 or MTi==29 or MTi==30 or MTi==32 or MTi==33 or
				MTi==34 or MTi==35 or MTi==36 or MTi==41 or MTi==42 or
				MTi==44 or MTi==45):

					ipstart = 1

			# Subsection before gamma energy distribution - NKprtl
				if (NK > 1):
					NKprtl = NK-2
				if (NK == 1):
					NKprtl = NK-1

			# neutron ZAP=1, AWP=1; charged particles ZAP>1, AWP<=4
				for NSS in range (NK):
					if (int(ZAP[NSS]) > 1 and float(AWP[NSS]) <= beta[lpr]):
						iflnrcyp = 1
						break
				for i in range (NP):
					Estar = n1[lpr]*Eall[i]
					availE = QI + (n2*Eall[i])
					if (availE < cbe[lpr]):
						Ea = availE
					if (cbe[lpr] < availE):
						Ea  = cbe[lpr]
					f1 = Estar
					f2 = 2*math.sqrt(beta[lpr]*Estar*Ea)
					f3 = beta[lpr]*Ea
					dn1 = 0
					dn2 = 0
					dn3 = 0
					if (sall[i] != 0): # .and. Eall(i)>=abs(QI)
						dn1 = abs(Tinteg(Z-ze[lpr],A+1-beta[lpr],Z,A,Ed,Eall[i],QI,beta[lpr],mdisp,bad,cad,n3,f1,f2,f3)) #,
						dn2 = abs(Tintegheat1(A,Eall[i],QI,beta[lpr],n3,f1,f2,f3))

						if (iflnrcyp == 0):
							dn3 = abs(Ea)

					if (sall[i] == 0):
						sdall[i] = 0
						shall[i] = 0
					if (sall[i] != 0):
						dn1 = dn1*0.8/(2*Ed)
						num_of_displ[i] = dn1
						sdall[i] = sall[i]*dn1
						if ((dn2+dn3)<availE):
							tot_energy_products[i] = dn2 + dn3
							shall[i] = sall[i]*(dn2+dn3)
						if ((dn2+dn3)>availE):
							tot_energy_products[i] = availE
							shall[i] = sall[i]*availE
					if (sdall[i] < 1e-12):
						sdall[i] = 0
					if (shall[i] < 1e-12):
						shall[i] = 0

				if (iflnrcyp == 1):
					print('', file = ofile_outRMINDD)
					print('Energy distribution data is given for emitted', file = ofile_outRMINDD)
					print('charged particle, but not for recoil nucleus', file = ofile_outRMINDD)

					kdim = int(numpy.amax(NEP))

					for iparticle in range (ipstart, NKprtl+1):
						ntm = [0]*NP
						for i in range (NP):
							ter6 = [0]*kdim
							tf6 = [0]*kdim
							for j in range (int(NE6[iparticle])-1):
								if (En[iparticle][j] == Eall[i]):
									for k in range (int(NEP[iparticle][j])):
										ter6[k] = Enp[iparticle][j][k]
										tf6[k] = f[iparticle][j][k]
									ntm[i] = NEP[iparticle][j]
									break
								if (En[iparticle][j] < Eall[i] and Eall[i] <= En[iparticle][j+1]):
									for k in range (int(NEP[iparticle][j])):
										x = Eall[i]
										x1 = En[iparticle][j]
										x2 = En[iparticle][j+1]
										y1 = Enp[iparticle][j][k]
										y2 = Enp[iparticle][j+1][k]
										y11 = f[iparticle][j][k]
										y22 = f[iparticle][j+1][k]
										ter6[k] = crstd(x,x1,x2,y1,y2)
										tf6[k] = crstd(x,x1,x2,y11,y22)
									ntm[i] = int(NEP[iparticle][j])
									break

							for k in range (int(ntm[i])):
								#if (iflcomlab == 1):
								#	ter6[k] = ter6[k] + (AWP[iparticle]/(A+1)**2)*Eall[i]
								#if (iflcom==1):
								#	ter6[k] = ter6[k] + (AWP[iparticle]/(A+1)**2)*Eall[i]
								if (ter6[k] == 1.0):
									ter6[k] = 0
								if (tf6[k] >= 1.0):
									tf6[k] = 0
							dn3 = 0
							if (sall[i] != 0): #  .and. Eall(i)>=abs(QI)
								dn3 = yi[iparticle][1]*Tinteg1heat(A,ter6,tf6,ntm[i])

							dn3 = abs(dn3)

							if (sall[i] == 0):
								shall[i] = 0
							if (sall[i] != 0):
								tot_energy_products[i] = tot_energy_products[i] + dn3
								shall[i] = shall[i] + sall[i]*dn3
							if (shall[i] <1e-12):
								shall[i] = 0

				# if for iflnrcyp = 1 (particle data there, recoil not)

		# THE ENERGY CARRIED AWAY BY NEUTRONS AND PHOTONS. 
		# IT CORRESPONDS TO ONLY HEATING CALCULATIONS BY ENERGY BALANCE METHOD.

			for NSS in range(NK):
				tnYldg = numpy.zeros(NP)
				tnYldg = TERPOL (NBT[NSS][:], INTr[NSS][:], Nyld[NSS], Eint6[NSS][:], yi[NSS][:], NP, Eall)
				ntm = [0]*NP; ntmd = [0]*NP
				if ( (AWP[NSS] == 0 and ZAP[NSS] == 0) or (AWP[NSS] == 1.0 and ZAP[NSS] == 1.0) ):
					print('', file = ofile_outRMINDD)
					print('Secondary neutrons and photons energy', file = ofile_outRMINDD)
					print('distribution data are given in MF 6 MT 5', file = ofile_outRMINDD)
	
					# for neutrons
					if ( AWP[NSS] == 1.0 and ZAP[NSS] == 1.0 ):
						print("Neutrons (Z and A):: ", int(ZAP[NSS]), AWP[NSS])
					# for photons ----
					if ( AWP[NSS] == 0.0 and ZAP[NSS] == 0.0 ):
						print("Photons (Z and A):: ", int(ZAP[NSS]), AWP[NSS])
					for i in range (NP):
						ter6 = [0]*kdim; tf6 = [0]*kdim
						for j in range (int(NE6[NSS])):
							if (En[NSS][j] == Eall[i]):
								for k in range (int(NEP[NSS][j])):
									ter6[k] = Enp[NSS][j][k]
									tf6[k] = f[NSS][j][k]
								ntm[i] = int(NEP[NSS][j])
								break
	
							if (En[NSS][j] < Eall[i] and Eall[i] < En[NSS][j+1]):
								for k in range (int(NEP[NSS][j])):
									x = Eall[i]
									x1 = En[NSS][j]
									x2 = En[NSS][j+1]
									y1 = Enp[NSS][j][k]
									y2 = Enp[NSS][j+1][k]
									y11 = f[NSS][j][k]
									y22 = f[NSS][j+1][k]
									if ((x2 - x1) != 0):
										ter6[k] = TERPOLIN(2,x,x1,x2,y1,y2)
										tf6[k] = TERPOLIN(2,x,x1,x2,y11,y22)
								ntm[i] = int(NEP[NSS][j])
								break
						for k in range (int(ntm[i])):
						#	if (iflcomlab == 1):
						#		ter6[k] = ter6[k] + (AWP[NSS]/(AWP[NSS]+1)**2)*Eall[i]
						#	if (iflcom == 1):
						#		ter6[k] = ter6[k] + (AWP[NSS]/(AWP[NSS]+1)**2)*Eall[i]
							if (ter6[k] == 1.0):
								ter6[k] = 0
							if (tf6[k] >= 1.0):
								tf6[k] = 0
	
						if (sall[i] != 0 ):
							if ( AWP[NSS] == 1.0 and ZAP[NSS] == 1.0 ):
								dn4_neutrons[i] = abs(tnYldg[i]*Tinteg1heat(AWP[NSS],ter6,tf6,int(ntm[i])))
	
							if ( AWP[NSS] == 0.0 and ZAP[NSS] == 0.0 ):
								dn5_photons[i] = abs(tnYldg[i]*Tinteg1heat(AWP[NSS],ter6,tf6,int(ntm[i])))
	
							tot_energy_n_photons[i] = dn4_neutrons[i] + dn5_photons[i]
							shall_EB[i] = sall[i]*(Eall[i] - tot_energy_n_photons[i])		

		# if for the File 6 data are present

		if (ifl4 == 0 and ifl6 == 0):
			print('', file = ofile_outRMINDD)
			print('Reaction product distribution data are not given', file = ofile_outRMINDD)
			print('in File 4 and File 6', file = ofile_outRMINDD)

			for i in range (NP):
				Estar = n1[lpr]*Eall[i]
				availE = QI + (n2*Eall[i])	# if (QI<0)
    			#if (QI>0) availE = n2*Eall(i)
				if (availE < cbe[lpr]):
					Ea = availE
				if (cbe[lpr] < availE):
					Ea  = cbe[lpr]
				f1 = Estar
				f2 = 2*math.sqrt(beta[lpr]*Estar*Ea)
				f3 = beta[lpr]*Ea
				dn1 = 0
				dn2 = 0
				dn3 = 0
				if (sall[i] != 0): # .and. Eall(i)>=abs(QI)
					dn1 = abs(Tinteg(Z-ze[lpr],A+1-beta[lpr],Z,A,Ed,Eall[i],QI,beta[lpr],mdisp,bad,cad,n3,f1,f2,f3)) #,
					dn2 = abs(Tintegheat1(A,Eall[i],QI,beta[lpr],n3,f1,f2,f3))
					dn3 = abs(Ea)
				if (sall[i] == 0):
					sdall[i] = 0
					shall[i] = 0
				if (sall[i] != 0):
					dn1 = dn1*0.8/(2*Ed)
					num_of_displ[i] = dn1
					sdall[i] = sall[i]*dn1
					if ((dn2+dn3) < availE):
						tot_energy_products[i] = dn2+dn3
						shall[i] = sall[i]*(dn2+dn3)
					if ((dn2+dn3) > availE):
						tot_energy_products[i] = availE
						shall[i] = sall[i]*availE
				if (sdall[i] < 1e-12):
					sdall[i] = 0
				if (shall[i] < 1e-12):
					shall[i] = 0

		# if neither File 4 nor File 6 is present.

		if (600 <= MTi and MTi <= 849):
			if (ifspad4 == 1 or ifspad6 == 1):
				ifdpd = 1		

		if (ifllargeNP == 0):
			Etu = numpy.asarray(Etu)
			Eall = numpy.asarray(Eall)
			sdall = numpy.asarray(sdall)
			shall = numpy.asarray(shall)
			shall_EB = numpy.asarray(shall_EB)
			sall = numpy.asarray(sall)
			num_of_displ = numpy.asarray(num_of_displ)
			tot_energy_products = numpy.asarray(tot_energy_products)
			tot_energy_n_photons = numpy.asarray(tot_energy_n_photons)
			sdpat = trptuqce(Eall,sdall,Etu)
			snhtt = trptuqce(Eall,shall,Etu)
			snhtt_EB = trptuqce(Eall,shall_EB,Etu)
			signcpol = trptuqce(Eall,sall,Etu)
			num_of_displ1 = trptuqce(Eall,num_of_displ,Etu)
			tot_energy_products1 = trptuqce(Eall,tot_energy_products,Etu)
			tot_energy_n_photons1 = trptuqce(Eall,tot_energy_n_photons,Etu)

		if (ifllargeNP == 1):
			sdpat = sdall
			snhtt = shall
			snhtt_EB = shall_EB
			signcpol = sall
			num_of_displ1 = num_of_displ
			tot_energy_products1 = tot_energy_products
			tot_energy_n_photons1 = tot_energy_n_photons

	#--------------
	# **** the above is done only if MF3 for that MT is present ****

	return (signcpol,num_of_displ1,tot_energy_products1,tot_energy_n_photons1,sdpat,snhtt,snhtt_EB,ifdpd,iflpresent)

#=======Contributions from (n, xn) reactions=======*

	# Controls the computation of dpa and heating cross sections due to 
	# the (n,2n), (n,3n) and (n,4n) reactions. Contributions from these
	# three reactions are added and kept in an arbitrarily assigned MT
	# number 1601. This number is also used in the file name giving the
	# total from (n,xn) reactions.

def CONTROL_nxn (ofile_outRMINDD,NPt,Etu,mdisp,Ed,bad,cad):
	dpaxn1 = [0]*NPt; snhtxn1 = [0]*NPt; signxntot = [0]*NPt
	num_of_displnxntot = [0]*NPt; tot_energy_productsnxntot = [0]*NPt
		# snhtt = Heating cross section from (n,(i)n)				
		# snhtxn1=total Heating cross section from (n,xn) 
	MTnum = [0]*3
	MTnum[0] = 16 # (n,2n) 
	MTnum[1] = 17 # (n,3n) 
	MTnum[2] = 37 # (n,4n)

	print( 'n, xn .....')

	for i in range(3):
		iflMTpr = 0
		MTfind = MTnum[i]
		iflMTpr = FindMT (MTfind)
		if (iflMTpr == 1):
			(signxnl, num_disp, tot_en_products, sdpat, snhtt, iflpresent) = n_xn (ofile_outRMINDD,MTnum[i],NPt,Etu,mdisp,Ed,bad,cad)
			if (iflpresent == 1):
				MTc = MTnum[i]
				print ( 'MT = ', MTc)
				for isum in range (NPt):
					signxntot[isum] = signxntot[isum] + signxnl[isum]
					num_of_displnxntot[isum] = num_of_displnxntot[isum] + num_disp[isum]
					tot_energy_productsnxntot[isum] = tot_energy_productsnxntot[isum] + tot_en_products[isum]
					dpaxn1[isum] = dpaxn1[isum] + sdpat[isum]
					snhtxn1[isum] = snhtxn1[isum] + snhtt[isum]

	# 1601 = id for (n,xn) output file, 1 = DPA
	printtofile (NPt,Etu,signxntot, num_of_displnxntot, dpaxn1,1601,0,1)
	printtofile (NPt,Etu,signxntot, tot_energy_productsnxntot, snhtxn1,1601,0,2)

#==============================================

#=======Contributions from thresholds giving out charged particles=======*

	# Controls the computation of dpa and heating cross sections due to
	# the (n,CPO) reactions. Contributions from these reactions are added
	# and kept in an arbitrarily assigned MT number 3001. This number is 
	# also used in the file name giving the total from (n,CPO) reactions.

def CONTROL_nCPO (ofile_outRMINDD,NPt,Etu,mdisp,Ed,bad,cad):

	tmdpaMT103 = [0]*NPt; tmdpaMT104 = [0]*NPt
	tmdpaMT105 = [0]*NPt; tmdpaMT106 = [0]*NPt; tmdpaMT107 = [0]*NPt
	dpaMT103 = [0]*NPt; dpaMT104 = [0]*NPt; dpaMT105 = [0]*NPt
	dpaMT106 = [0]*NPt; dpaMT107 = [0]*NPt; dpa3 = [0]*NPt; dpaMT5 = [0]*NPt
	tmnhtMT103 = [0]*NPt; tmnhtMT104 = [0]*NPt; tmnhtMT105 = [0]*NPt
	tmnhtMT106 = [0]*NPt; tmnhtMT107 = [0]*NPt; snhtMT103 = [0]*NPt
	snhtMT104 = [0]*NPt; snhtMT105 = [0]*NPt; snhtMT106 = [0]*NPt
	snhtMT107 = [0]*NPt; snht3 = [0]*NPt; snhtMT5 = [0]*NPt
	tmnhtMT103_EB = [0]*NPt; tmnhtMT104_EB = [0]*NPt; tmnhtMT105_EB = [0]*NPt
	tmnhtMT106_EB = [0]*NPt; tmnhtMT107_EB = [0]*NPt; snhtMT103_EB = [0]*NPt
	snhtMT104_EB = [0]*NPt; snhtMT105_EB = [0]*NPt; snhtMT106_EB = [0]*NPt
	snhtMT107_EB = [0]*NPt; snht3_EB = [0]*NPt; snhtMT5_EB = [0]*NPt

	signcpo3 = [0]*NPt; num_of_displ3 = [0]*NPt; tot_energy_products3 = [0]*NPt; tot_energy_n_photons3 = [0]*NPt
	signcpoMT103 = [0]*NPt; num_of_displMT103 = [0]*NPt; tot_energy_productsMT103 = [0]*NPt; tot_energy_n_photonsMT103 = [0]*NPt
	signcpoMT104 = [0]*NPt; num_of_displMT104 = [0]*NPt; tot_energy_productsMT104 = [0]*NPt; tot_energy_n_photonsMT104 = [0]*NPt
	signcpoMT105 = [0]*NPt; num_of_displMT105 = [0]*NPt; tot_energy_productsMT105 = [0]*NPt; tot_energy_n_photonsMT105 = [0]*NPt
	signcpoMT106 = [0]*NPt; num_of_displMT106 = [0]*NPt; tot_energy_productsMT106 = [0]*NPt; tot_energy_n_photonsMT106 = [0]*NPt
	signcpoMT107 = [0]*NPt; num_of_displMT107 = [0]*NPt; tot_energy_productsMT107 = [0]*NPt; tot_energy_n_photonsMT107 = [0]*NPt

	# dapnth, snhtt = individual MTnum(i) contribution
	# dpa3, snht3 = total dpa,heating for (n,CPO)
	# dpaMT5 = dpa contribution from MT=5

	MTnum = numpy.zeros(281)	# 31 + 250

	# MT103: 600 to 649 
	# MT104: 650 to 699
	# MT105: 700 to 749 
	# MT106: 750 to 799 
	# MT107: 800 to 849 

	MTtpnum = numpy.zeros(5)	# for MT=203,204,205,206,207 in JENDL and the likes

	Eyld = numpy.zeros((5,200)); Yld = numpy.zeros((5,200))
		# Energy and Yield of each CP from file 6 (p,d,t,3He,4He),200=maximum energy points
	Nyld = numpy.zeros(5); ifltfr5 = numpy.zeros(5)		# Read from file 6
	NBTp = numpy.zeros((5,20)); INTrp = numpy.zeros((5,20))				# Read from file 6
	NBTpp = [0]*20; NBTpd = [0]*20; NBTpt = [0]*20; NBTp3He = [0]*20
	NBTpal = [0]*20; INTrpp = [0]*20; INTrpd = [0]*20; INTrpt = [0]*20
	INTrp3He = [0]*20; INTrpal = [0]*20

	beta = [0]*5; ze = [0]*5; cbe = [0]*5; n1 = [0]*5
		# beta = mass of emitted particle/mass of neutron, ze = charge
		# ze = 1 for p,d,t; 2 for 3He and al
		# cbe = Coulomb Barrier Energy, pg 129 njoy16 manual
		# n1 = E*/E in pg 126 njoy16 manual

		# n2=A/A+1, n3 = 1.0/A+1

	print('n, CPO .....', file = ofile_outRMINDD)
	print('n, CPO .....')

	MTnum[0]=11;MTnum[1]=22;MTnum[2]=23;MTnum[3]=24
	MTnum[4]=25;MTnum[5]=28;MTnum[6]=29;MTnum[7]=30
	MTnum[8]=32;MTnum[9]=33;MTnum[10]=34;MTnum[11]=35
	MTnum[12]=36;MTnum[13]=41;MTnum[14]=42;MTnum[15]=44
	MTnum[16]=45;MTnum[17]=103;MTnum[18]=104;MTnum[19]=105
	MTnum[20]=106;MTnum[21]=107;MTnum[22]=108;MTnum[23]=109
	MTnum[24]=111;MTnum[25]=112;MTnum[26]=113;MTnum[27]=114
	MTnum[28]=115;MTnum[29]=116;MTnum[30]=117

	for i in range (31, 281):
		MTnum[i] = i+569

	MTtpnum[0]=203;MTtpnum[1]=204;MTtpnum[2]=205;MTtpnum[3]=206
	MTtpnum[4]=207

	# ifdisc103, etc. -- for discrete data representation
	# ifdisdata103, etc. -- for presence of data for discrete levels
	# iffldd103, etc. -- for presence of full discrete + continuum data

	# ifdisc = if discrete level data are present i.e 600+ for MT 103

	ifdisc103=0;ifdisc104=0;ifdisc105=0;ifdisc106=0;ifdisc107=0 
	ifdisdata103=0;ifdisdata104=0;ifdisdata105=0;ifdisdata106=0
	ifdisdata107=0

	iffldd103=0;iffldd104=0;iffldd105=0;iffldd106=0;iffldd107=0

	ifl103pr = 0 # flags for presence of 
	ifl104pr = 0 # MT=103, 104, .., 107
	ifl105pr = 0
	ifl106pr = 0
	ifl107pr = 0

	# flag to indicate whether do for individual particle MTs or not
	iflMTtppr = 0 

	# IF MT = 203, 204, 205, 206, 207 ARE PRESENT
	# i.e. if (n,xp), (n,xd), (n,xt), (n,x3He), (n,xa) are present
	# it is asked at during execution whether to use these MTs
	# for calculation		

	for j in range(5):
		iflMTpr = 0
		iflMTpr = FindMT(MTtpnum[j])
		if (iflMTpr==1):
			print('', file = ofile_outRMINDD)
			print(':: MESSAGE :: CPO reaction cross section', file = ofile_outRMINDD)
			print('-------------------------------------------------', file = ofile_outRMINDD)
			print('Total charged particle production MTs are found.', file = ofile_outRMINDD)
			print('', file = ofile_outRMINDD)
			print('Total charged particle production MTs are found.')
			print('Calculate from these reaction cross sections?')
			print('..... Enter: 0 / 1 for No / Yes .....')
			iflMTtppr = input()
			if(iflMTtppr == 1):
				print('calculating from these MTs')
			if(iflMTtppr == 1):
				print('User chose to consider these MTs', file = ofile_outRMINDD)
			if(iflMTtppr == 0):
				print('not considering these MTs')
			if(iflMTtppr == 0):
				print('User chose not to consider the total particle', file = ofile_outRMINDD)
				print('production MTs. So the individual (n,CPO) MTs', file = ofile_outRMINDD)
				print('are considered', file = ofile_outRMINDD)
			break

	if (iflMTtppr == 1):
		for j in range(5):
			iflMTpr = 0
			MTfind = MTtpnum[j]
			if (MTfind == 203):
				lpr = 18
			if (MTfind == 204):
				lpr = 19
			if (MTfind == 205):
				lpr = 20
			if (MTfind == 206):
				lpr = 21
			if (MTfind == 207):
				lpr = 22
			iflMTpr = FindMT(MTfind)
			if (iflMTpr == 1):
				(signcpol,num_of_displ1,tot_energy_products1, \
				tot_energy_n_photons1,dpanth, snhtt, snhtt_EB, ifdpd, iflpresent) = \
											n_CPO (ofile_outRMINDD,MTtpnum[j],lpr,mdisp,Ed,bad,cad,NPt,Etu)
				if (iflpresent == 1):
					MTc = MTtpnum[j]
					print('MT = ', MTc)
					if (MTc == 203):
						dpaMT103 = dpanth
						snhtMT103 = snhtt
						snhtMT103 = snhtt_EB
						signcpoMT103 = signcpol
						num_of_displMT103 = num_of_displ1
						tot_energy_productsMT103 = tot_energy_products1
						tot_energy_n_photonsMT103 = tot_energy_n_photons1
						printtofile(NPt,Etu,signcpol,num_of_displ1,dpaMT103,203,0,1)
						printtofile(NPt,Etu,signcpol,tot_energy_products1,snhtMT103,203,0,2)
						printtofile(NPt,Etu,signcpol,tot_energy_products1,snhtMT103_EB,203,0,3)
					if (MTc == 204):
						dpaMT104 = dpanth
						snhtMT104 = snhtt
						snhtMT104 = snhtt_EB
						signcpoMT104 = signcpol
						num_of_displMT104 = num_of_displ1
						tot_energy_productsMT104 = tot_energy_products1
						tot_energy_n_photonsMT104 = tot_energy_n_photons1
						printtofile(NPt,Etu,signcpol,num_of_displ1,dpaMT104,204,0,1)
						printtofile(NPt,Etu,signcpol,tot_energy_products1,snhtMT104,204,0,2)
						printtofile(NPt,Etu,signcpol,tot_energy_products1,snhtMT104_EB,204,0,3)
					if (MTc == 205):
						dpaMT105 = dpanth
						snhtMT105 = snhtt
						snhtMT105 = snhtt_EB
						signcpoMT105 = signcpol
						num_of_displMT105 = num_of_displ1
						tot_energy_productsMT105 = tot_energy_products1
						tot_energy_n_photonsMT105 = tot_energy_n_photons1
						printtofile(NPt,Etu,signcpol,num_of_displ1,dpaMT105,205,0,1)
						printtofile(NPt,Etu,signcpol,tot_energy_products1,snhtMT105,205,0,2)
						printtofile(NPt,Etu,signcpol,tot_energy_products1,snhtMT105_EB,205,0,3)
					if (MTc == 206):
						dpaMT106 = dpanth
						snhtMT106 = snhtt
						snhtMT106 = snhtt_EB
						signcpoMT106 = signcpol
						num_of_displMT106 = num_of_displ1
						tot_energy_productsMT106 = tot_energy_products1
						tot_energy_n_photonsMT106 = tot_energy_n_photons1
						printtofile(NPt,Etu,signcpol,num_of_displ1,dpaMT106,206,0,1)
						printtofile(NPt,Etu,signcpol,tot_energy_products1,snhtMT106,206,0,2)
						printtofile(NPt,Etu,signcpol,tot_energy_products1,snhtMT106_EB,206,0,3)
					if (MTc == 207):
						dpaMT107 = dpanth
						snhtMT107 = snhtt
						snhtMT107 = snhtt_EB
						signcpoMT107 = signcpol
						num_of_displMT107 = num_of_displ1
						tot_energy_productsMT107 = tot_energy_products1
						tot_energy_n_photonsMT107 = tot_energy_n_photons1
						printtofile(NPt,Etu,signcpol,num_of_displ1,dpaMT107,207,0,1)
						printtofile(NPt,Etu,signcpol,tot_energy_products1,snhtMT107,207,0,2)
						printtofile(NPt,Etu,signcpol,tot_energy_products1,snhtMT107_EB,207,0,3)
				# for the iflpresent
 			# for the iflMTpr
 		# for loop
	# for the iflMTtppr
	#--------------------------------

	# if MT=203, 204, etc are not present
	# if MT=203, 204, etc are present, but not considering
	# if individual particle MTs are considered

	if (iflMTtppr == 0):
		for j in range (281):
			iflMTpr = 0
			MTfind = MTnum[j]	
			iflMTpr = FindMT(MTfind)
			if (iflMTpr == 1):
				if (MTfind == 103):
					ifl103pr = 1
				if (MTfind == 104):
					ifl104pr = 1
				if (MTfind == 105):
					ifl105pr = 1
				if (MTfind == 106):
					ifl106pr = 1
				if (MTfind == 107):
					ifl107pr = 1

				(signcpol,num_of_displ1,tot_energy_products1, \
				tot_energy_n_photons1,dpanth, snhtt, snhtt_EB, ifdpd, iflpresent) = \
											n_CPO (ofile_outRMINDD,int(MTnum[j]),j,mdisp,Ed,bad,cad,NPt,Etu)
				if (iflpresent == 1):
					MTc = int(MTnum[j])
					print('MT = ', MTc)
	# dpa and heat from MT=103,..,107 stored temporarily and used
	# if disc. + cont. data are not there / only disc. data is there
	# but cont. data is not given
				if (MTc == 103):
					tmdpaMT103 = dpanth
					tmnhtMT103 = snhtt
					tmnhtMT103_EB = snhtt_EB
					tmsigncpoMT103 = signcpol
					tmnum_ofdisplMT103 = num_of_displ1
					tmtot_energy_productsMT103 = tot_energy_products1
					tmtot_energy_n_photonsMT103 = tot_energy_n_photons1
				if (MTc == 104):
					tmdpaMT104 = dpanth
					tmnhtMT104 = snhtt
					tmnhtMT104_EB = snhtt_EB
					tmsigncpoMT104 = signcpol
					tmnum_ofdisplMT104 = num_of_displ1
					tmtot_energy_productsMT104 = tot_energy_products1
					tmtot_energy_n_photonsMT104 = tot_energy_n_photons1
				if (MTc == 105):
					tmdpaMT105 = dpanth
					tmnhtMT105 = snhtt
					tmnhtMT105_EB = snhtt_EB
					tmsigncpoMT105 = signcpol
					tmnum_ofdisplMT105 = num_of_displ1
					tmtot_energy_productsMT105 = tot_energy_products1
					tmtot_energy_n_photonsMT105 = tot_energy_n_photons1
				if (MTc == 106):
					tmdpaMT106 = dpanth
					tmnhtMT106 = snhtt
					tmnhtMT106_EB = snhtt_EB
					tmsigncpoMT106 = signcpol
					tmnum_ofdisplMT106 = num_of_displ1
					tmtot_energy_productsMT106 = tot_energy_products1
					tmtot_energy_n_photonsMT106 = tot_energy_n_photons1
				if (MTc == 107):
					tmdpaMT107 = dpanth
					tmnhtMT107 = snhtt
					tmnhtMT107_EB = snhtt_EB
					tmsigncpoMT107 = signcpol
					tmnum_ofdisplMT107 = num_of_displ1
					tmtot_energy_productsMT107 = tot_energy_products1
					tmtot_energy_n_photonsMT107 = tot_energy_n_photons1

				if (MTc != 103 and MTc != 104 and MTc != 105 and \
					MTc != 106 and MTc != 107):

					if (MTc < 600):

						printtofile(NPt,Etu,signcpol,num_of_displ1,dpanth,MTc,0,1)
						printtofile(NPt,Etu,signcpol,tot_energy_products1,snhtt,MTc,0,2)
						printtofile(NPt,Etu,signcpol,tot_energy_n_photons1,snhtt_EB,MTc,0,3)

						for isum in range (NPt):
							dpa3[isum] = dpa3[isum] + dpanth[isum]
							snht3[isum] = snht3[isum] + snhtt[isum]
							snht3_EB[isum] = snht3_EB[isum] + snhtt_EB[isum]
							signcpo3[isum] = signcpo3[isum] + signcpol[isum]
							num_of_displ3[isum] = num_of_displ3[isum] + num_of_displ1[isum]
							tot_energy_products3[isum] = tot_energy_products3[isum] + tot_energy_products1[isum]
							tot_energy_n_photons3[isum] = tot_energy_n_photons3[isum] + tot_energy_n_photons1[isum]

				if (600<=MTc and MTc<=649):
					ifdisdata103 = 1
					if (ifdpd == 1):
						ifdisc103 = 1
						for isum in range (NPt):
							dpaMT103[isum] = dpaMT103[isum] + dpanth[isum]
							snhtMT103[isum] = snhtMT103[isum] + snhtt[isum]
							snhtMT103_EB[isum] = snhtMT103_EB[isum] + snhtt_EB[isum]
							signcpoMT103[isum] = signcpoMT103[isum] + signcpol[isum]
							num_of_displMT103[isum] = num_of_displMT103[isum] + num_of_displ1[isum]
							tot_energy_productsMT103[isum] = tot_energy_productsMT103[isum] + tot_energy_products1[isum]
							tot_energy_n_photonsMT103[isum] = tot_energy_n_photonsMT103[isum] + tot_energy_n_photons1[isum]
						printtofile(NPt,Etu,signcpol,num_of_displ1,dpanth,MTc,0,1)
						printtofile(NPt,Etu,signcpol,tot_energy_products1,snhtt,MTc,0,2)
						printtofile(NPt,Etu,signcpol,tot_energy_n_photons1,snhtt_EB,MTc,0,3)

					if (ifdpd == 0 and ifdisc103 == 1):
						for isum in range(NPt):
							dpaMT103[isum] = dpaMT103[isum] + dpanth[isum]
							snhtMT103[isum] = snhtMT103[isum] + snhtt[isum]
							snhtMT103_EB[isum] = snhtMT103_EB[isum] + snhtt_EB[isum]
							signcpoMT103[isum] = signcpoMT103[isum] + signcpol[isum]
							num_of_displMT103[isum] = num_of_displMT103[isum] + num_of_displ1[isum]
							tot_energy_productsMT103[isum] = tot_energy_productsMT103[isum] + tot_energy_products1[isum]
							tot_energy_n_photonsMT103[isum] = tot_energy_n_photonsMT103[isum] + tot_energy_n_photons1[isum]
						printtofile(NPt,Etu,signcpol,num_of_displ1,dpanth,MTc,0,1)
						printtofile(NPt,Etu,signcpol,tot_energy_products1,snhtt,MTc,0,2)
						printtofile(NPt,Etu,signcpol,tot_energy_n_photons1,snhtt_EB,MTc,0,3)
					if (ifdpd == 0 and ifdisc103 == 0):
						dpaMT103 = dpanth 
						snhtMT103 = snhtt
						snhtMT103_EB[isum] = snhtt_EB[isum]
						signcpoMT103 = signcpol
						num_of_displMT103 = num_of_displ1
						tot_energy_productsMT103 = tot_energy_products1
						tot_energy_n_photonsMT103[isum] = tot_energy_n_photons1[isum]
						printtofile(NPt,Etu,signcpol,num_of_displ1,dpanth,MTc,0,1)
						printtofile(NPt,Etu,signcpol,tot_energy_products1,snhtt,MTc,0,2)
						printtofile(NPt,Etu,signcpol,tot_energy_n_photons1,snhtt_EB,MTc,0,3)
					if (MTc == 649):
						iffldd103 = 1
						print('..... taking disc.+cont. (n,p)')
						printtofile(NPt,Etu,signcpol,num_of_displ1,dpanth,MTc,0,1)
						printtofile(NPt,Etu,signcpol,tot_energy_products1,snhtt,MTc,0,2)
						printtofile(NPt,Etu,signcpol,tot_energy_n_photons1,snhtt_EB,MTc,0,3)
						printtofile(NPt,Etu,signcpoMT103,num_of_displMT103,dpaMT103,103,0,1)
						printtofile(NPt,Etu,signcpoMT103,tot_energy_productsMT103,snhtMT103,103,0,2)
						printtofile(NPt,Etu,signcpoMT103,tot_energy_n_photonsMT103,snhtMT103_EB,103,0,3)
				#-------------------

				if (650 <= MTc and MTc <= 699):
					ifdisdata104 = 1
					if (ifdpd == 1):
						ifdisc104 = 1
						for isum in range (NPt):
							dpaMT104[isum] = dpaMT104[isum] + dpanth[isum]
							snhtMT104[isum] = snhtMT104[isum] + snhtt[isum]
							snhtMT104_EB[isum] = snhtMT104_EB[isum] + snhtt_EB[isum]
							signcpoMT104[isum] = signcpoMT104[isum] + signcpol[isum]
							num_of_displMT104[isum] = num_of_displMT104[isum] + num_of_displ1[isum]
							tot_energy_productsMT104[isum] = tot_energy_productsMT104[isum] + tot_energy_products1[isum]
							tot_energy_n_photonsMT104[isum] = tot_energy_n_photonsMT104[isum] + tot_energy_n_photons1[isum]
						printtofile(NPt,Etu,signcpol,num_of_displ1,dpanth,MTc,0,1)
						printtofile(NPt,Etu,signcpol,tot_energy_products1,snhtt,MTc,0,2)
						printtofile(NPt,Etu,signcpol,tot_energy_n_photons1,snhtt_EB,MTc,0,3)
					if (ifdpd == 0 and ifdisc104 == 1):
						for isum in range (NPt):
							dpaMT104[isum] = dpaMT104[isum] + dpanth[isum]
							snhtMT104[isum] = snhtMT104[isum] + snhtt[isum]
							snhtMT104_EB[isum] = snhtMT104_EB[isum] + snhtt_EB[isum]
							signcpoMT104[isum] = signcpoMT104[isum] + signcpol[isum]
							num_of_displMT104[isum] = num_of_displMT104[isum] + num_of_displ1[isum]
							tot_energy_productsMT104[isum] = tot_energy_productsMT104[isum] + tot_energy_products1[isum]
							tot_energy_n_photonsMT104[isum] = tot_energy_n_photonsMT104[isum] + tot_energy_n_photons1[isum]
						printtofile(NPt,Etu,signcpol,num_of_displ1,dpanth,MTc,0,1)
						printtofile(NPt,Etu,signcpol,tot_energy_products1,snhtt,MTc,0,2)
						printtofile(NPt,Etu,signcpol,tot_energy_n_photons1,snhtt_EB,MTc,0,3)
					if (ifdpd == 0 and ifdisc104 == 0):
						dpaMT104 = dpanth
						snhtMT104 = snhtt
						snhtMT104_EB = snhtt_EB
						signcpoMT104 = signcpol
						num_of_displMT104 = num_of_displ1
						tot_energy_productsMT104 = tot_energy_products1
						tot_energy_n_photonsMT104 = tot_energy_n_photons1
						printtofile(NPt,Etu,signcpol,num_of_displ1,dpanth,MTc,0,1)
						printtofile(NPt,Etu,signcpol,tot_energy_products1,snhtt,MTc,0,2)
						printtofile(NPt,Etu,signcpol,tot_energy_n_photons1,snhtt_EB,MTc,0,3)
					if (MTc == 699):
						iffldd104 = 1
						print('..... taking disc.+cont. (n,d)')
						printtofile(NPt,Etu,signcpol,num_of_displ1,dpanth,MTc,0,1)
						printtofile(NPt,Etu,signcpol,tot_energy_products1,snhtt,MTc,0,2)
						printtofile(NPt,Etu,signcpol,tot_energy_n_photons1,snhtt_EB,MTc,0,3)
						printtofile(NPt,Etu,signcpoMT104,num_of_displMT104,dpaMT104,104,0,1)
						printtofile(NPt,Etu,signcpoMT104,tot_energy_productsMT104,snhtMT104,104,0,2)
						printtofile(NPt,Etu,signcpoMT104,tot_energy_n_photonsMT104,snhtMT104_EB,104,0,3)
				#-------------------

				if (700 <= MTc and MTc <= 749):
					ifdisdata105 = 1
					if (ifdpd == 1):
						ifdisc105 = 1
						for isum in range(NPt):
							dpaMT105[isum] = dpaMT105[isum] + dpanth[isum]
							snhtMT105[isum] = snhtMT105[isum] + snhtt[isum]
							snhtMT105_EB[isum] = snhtMT105_EB[isum] + snhtt_EB[isum]
							signcpoMT105[isum] = signcpoMT105[isum] + signcpol[isum]
							num_of_displMT105[isum] = num_of_displMT105[isum] + num_of_displ1[isum]
							tot_energy_productsMT105[isum] = tot_energy_productsMT105[isum] + tot_energy_products1[isum]
							tot_energy_n_photonsMT105[isum] = tot_energy_n_photonsMT105[isum] + tot_energy_n_photons1[isum]
						printtofile(NPt,Etu,signcpol,num_of_displ1,dpanth,MTc,0,1)
						printtofile(NPt,Etu,signcpol,tot_energy_products1,snhtt,MTc,0,2)
						printtofile(NPt,Etu,signcpol,tot_energy_n_photons1,snhtt_EB,MTc,0,3)
					if (ifdpd == 0 and ifdisc105==1): 		
						for isum in range(NPt):
							dpaMT105[isum] = dpaMT105[isum] + dpanth[isum]
							snhtMT105[isum] = snhtMT105[isum] + snhtt[isum]
							snhtMT105_EB[isum] = snhtMT105_EB[isum] + snhtt_EB[isum]
							signcpoMT105[isum] = signcpoMT105[isum] + signcpol[isum]
							num_of_displMT105[isum] = num_of_displMT105[isum] + num_of_displ1[isum]
							tot_energy_productsMT105[isum] = tot_energy_productsMT105[isum] + tot_energy_products1[isum]
							tot_energy_n_photonsMT105[isum] = tot_energy_n_photonsMT105[isum] + tot_energy_n_photons1[isum]
						printtofile(NPt,Etu,signcpol,num_of_displ1,dpanth,MTc,0,1)
						printtofile(NPt,Etu,signcpol,tot_energy_products1,snhtt,MTc,0,2)
						printtofile(NPt,Etu,signcpol,tot_energy_n_photons1,snhtt_EB,MTc,0,3)
					if (ifdpd == 0 and ifdisc105 == 0):
						dpaMT105 = dpanth
						snhtMT105 = snhtt
						snhtMT105_EB = snhtt_EB
						signcpoMT105 = signcpol
						num_of_displMT105 = num_of_displ1
						tot_energy_productsMT105 = tot_energy_products1
						tot_energy_n_photonsMT105 = tot_energy_n_photons1
						printtofile(NPt,Etu,signcpol,num_of_displ1,dpanth,MTc,0,1)
						printtofile(NPt,Etu,signcpol,tot_energy_products1,snhtt,MTc,0,2)
						printtofile(NPt,Etu,signcpol,tot_energy_n_photons1,snhtt_EB,MTc,0,3)
					if (MTc == 749):
						iffldd105 = 1
						print('..... taking disc.+cont. (n,t)')
						printtofile(NPt,Etu,signcpol,num_of_displ1,dpanth,MTc,0,1)
						printtofile(NPt,Etu,signcpol,tot_energy_products1,snhtt,MTc,0,2)
						printtofile(NPt,Etu,signcpol,tot_energy_n_photons1,snhtt_EB,MTc,0,3)
						printtofile(NPt,Etu,signcpoMT105,num_of_displMT105,dpaMT105,105,0,1)
						printtofile(NPt,Etu,signcpoMT105,tot_energy_productsMT105,snhtMT105,105,0,2)
						printtofile(NPt,Etu,signcpoMT105,tot_energy_n_photonsMT105,snhtMT105_EB,105,0,3)
				#-------------------

				if (750 <= MTc and MTc <= 799):
					ifdisdata106 = 1
					if (ifdpd == 1):
						ifdisc106 = 1
						for isum in range(NPt):
							dpaMT106[isum] = dpaMT106[isum] + dpanth[isum]
							snhtMT106[isum] = snhtMT106[isum] + snhtt[isum]
							snhtMT106_EB[isum] = snhtMT106_EB[isum] + snhtt_EB[isum]
							signcpoMT106[isum] = signcpoMT106[isum] + signcpol[isum]
							num_of_displMT106[isum] = num_of_displMT106[isum] + num_of_displ1[isum]
							tot_energy_productsMT106[isum] = tot_energy_productsMT106[isum] + tot_energy_products1[isum]
							tot_energy_n_photonsMT106[isum] = tot_energy_n_photonsMT106[isum] + tot_energy_n_photons1[isum]
						printtofile(NPt,Etu,signcpol,num_of_displ1,dpanth,MTc,0,1)
						printtofile(NPt,Etu,signcpol,tot_energy_products1,snhtt,MTc,0,2)
						printtofile(NPt,Etu,signcpol,tot_energy_n_photons1,snhtt_EB,MTc,0,3)
					if (ifdpd == 0 and ifdisc106 == 1):
						for isum in range (NPt):
							dpaMT106[isum] = dpaMT106[isum] + dpanth[isum]
							snhtMT106[isum] = snhtMT106[isum] + snhtt[isum]
							snhtMT106_EB[isum] = snhtMT106_EB[isum] + snhtt_EB[isum]
							signcpoMT106[isum] = signcpoMT106[isum] + signcpol[isum]
							num_of_displMT106[isum] = num_of_displMT106[isum] + num_of_displ1[isum]
							tot_energy_productsMT106[isum] = tot_energy_productsMT106[isum] + tot_energy_products1[isum]
							tot_energy_n_photonsMT106[isum] = tot_energy_n_photonsMT106[isum] + tot_energy_n_photons1[isum]
						printtofile(NPt,Etu,signcpol,num_of_displ1,dpanth,MTc,0,1)
						printtofile(NPt,Etu,signcpol,tot_energy_products1,snhtt,MTc,0,2)
						printtofile(NPt,Etu,signcpol,tot_energy_n_photons1,snhtt_EB,MTc,0,3)
					if (ifdpd == 0 and ifdisc106 == 0):
						dpaMT106 = dpanth
						snhtMT106 = snhtt
						snhtMT106_EB = snhtt_EB
						signcpoMT106 = signcpol
						num_of_displMT106 = num_of_displ1
						tot_energy_productsMT106 = tot_energy_products1
						tot_energy_n_photonsMT106 = tot_energy_n_photons1
						printtofile(NPt,Etu,signcpol,num_of_displ1,dpanth,MTc,0,1)
						printtofile(NPt,Etu,signcpol,tot_energy_products1,snhtt,MTc,0,2)
						printtofile(NPt,Etu,signcpol,tot_energy_n_photons1,snhtt_EB,MTc,0,3)
					if (MTc == 799):
						iffldd106 = 1
						print('..... taking disc.+cont. (n,3He)')
						printtofile(NPt,Etu,signcpol,num_of_displ1,dpanth,MTc,0,1)
						printtofile(NPt,Etu,signcpol,tot_energy_products1,snhtt,MTc,0,2)
						printtofile(NPt,Etu,signcpol,tot_energy_n_photons1,snhtt_EB,MTc,0,3)
						printtofile(NPt,Etu,signcpoMT106,num_of_displMT106,dpaMT106,106,0,1)
						printtofile(NPt,Etu,signcpoMT106,tot_energy_productsMT106,snhtMT106,106,0,2)
						printtofile(NPt,Etu,signcpoMT106,tot_energy_n_photonsMT106,snhtMT106_EB,106,0,3)
				# -------------------
				
				if (800 <= MTc and MTc <= 849):
					ifdisdata107 = 1
					if (ifdpd == 1):
						ifdisc107 = 1
						for isum in range(NPt):
							dpaMT107[isum] = dpaMT107[isum] + dpanth[isum]
							snhtMT107[isum] = snhtMT107[isum] + snhtt[isum]
							snhtMT107_EB[isum] = snhtMT107_EB[isum] + snhtt_EB[isum]
							signcpoMT107[isum] = signcpoMT107[isum] + signcpol[isum]
							num_of_displMT107[isum] = num_of_displMT107[isum] + num_of_displ1[isum]
							tot_energy_productsMT107[isum] = tot_energy_productsMT107[isum] + tot_energy_products1[isum]
							tot_energy_n_photonsMT107[isum] = tot_energy_n_photonsMT107[isum] + tot_energy_n_photons1[isum]
						printtofile(NPt,Etu,signcpol,num_of_displ1,dpanth,MTc,0,1)
						printtofile(NPt,Etu,signcpol,tot_energy_products1,snhtt,MTc,0,2)
						printtofile(NPt,Etu,signcpol,tot_energy_n_photons1,snhtt_EB,MTc,0,3)
					if (ifdpd == 0 and ifdisc107 == 1):
						for isum in range(NPt):
							dpaMT107[isum] = dpaMT107[isum] + dpanth[isum]
							snhtMT107[isum] = snhtMT107[isum] + snhtt[isum]
							snhtMT107_EB[isum] = snhtMT107_EB[isum] + snhtt_EB[isum]
							signcpoMT107[isum] = signcpoMT107[isum] + signcpol[isum]
							num_of_displMT107[isum] = num_of_displMT107[isum] + num_of_displ1[isum]
							tot_energy_productsMT107[isum] = tot_energy_productsMT107[isum] + tot_energy_products1[isum]
							tot_energy_n_photonsMT107[isum] = tot_energy_n_photonsMT107[isum] + tot_energy_n_photons1[isum]
						printtofile(NPt,Etu,signcpol,num_of_displ1,dpanth,MTc,0,1)
						printtofile(NPt,Etu,signcpol,tot_energy_products1,snhtt,MTc,0,2)
						printtofile(NPt,Etu,signcpol,tot_energy_n_photons1,snhtt_EB,MTc,0,3)
					if (ifdpd == 0 and ifdisc107 == 0):
						dpaMT107 = dpanth
						snhtMT107 = snhtt
						snhtMT107_EB = snhtt_EB
						signcpoMT107 = signcpol
						num_of_displMT107 = num_of_displ1
						tot_energy_productsMT107 = tot_energy_products1
						tot_energy_n_photonsMT107 = tot_energy_n_photons1
						printtofile(NPt,Etu,signcpol,num_of_displ1,dpanth,MTc,0,1)
						printtofile(NPt,Etu,signcpol,tot_energy_products1,snhtt,MTc,0,2)
						printtofile(NPt,Etu,signcpol,tot_energy_n_photons1,snhtt_EB,MTc,0,3)
					if (MTc == 849):
						iffldd107 = 1
						print('..... taking disc.+cont. (n,a)')
						printtofile(NPt,Etu,signcpol,num_of_displ1,dpanth,MTc,0,1)
						printtofile(NPt,Etu,signcpol,tot_energy_products1,snhtt,MTc,0,2)
						printtofile(NPt,Etu,signcpol,tot_energy_n_photons1,snhtt_EB,MTc,0,3)
						printtofile(NPt,Etu,signcpoMT107,num_of_displMT107,dpaMT107,107,0,1)
						printtofile(NPt,Etu,signcpoMT107,tot_energy_productsMT107,snhtMT107,107,0,2)
						printtofile(NPt,Etu,signcpoMT107,tot_energy_n_photonsMT107,snhtMT107_EB,107,0,3)
					# -------------------
				# for the iflpresent
			# for the iflMTpr
		#
		#-------------------
	
		if ((ifdisdata103 == 0 or iffldd103 == 0) and ifl103pr==1):
			dpaMT103 = tmdpaMT103
			snhtMT103 = tmnhtMT103
			snhtMT103_EB = tmnhtMT103_EB
			signcpoMT103 = tmsigncpoMT103
			num_of_displMT103 = tmnum_ofdisplMT103
			tot_energy_productsMT103 = tmtot_energy_productsMT103
			tot_energy_n_photonsMT103 = tmtot_energy_n_photonsMT103
			if (ifdisdata103==1):
				print('..... taking only MT=103')
			printtofile(NPt,Etu,signcpoMT103,num_of_displMT103,dpaMT103,103,0,1)
			printtofile(NPt,Etu,signcpoMT103,tot_energy_productsMT103,snhtMT103,103,0,2)
			printtofile(NPt,Etu,signcpoMT103,tot_energy_n_photonsMT103,snhtMT103_EB,103,0,3)
		# -------------------
		if ((ifdisdata104 == 0 or iffldd104 == 0) and ifl104pr == 1):
			dpaMT104 = tmdpaMT104
			snhtMT104 = tmnhtMT104
			snhtMT104_EB = tmnhtMT104_EB
			signcpoMT104 = tmsigncpoMT104
			num_of_displMT104 = tmnum_ofdisplMT104
			tot_energy_productsMT104 = tmtot_energy_productsMT104
			tot_energy_n_photonsMT104 = tmtot_energy_n_photonsMT104
			if (ifdisdata104 == 1):
				print('..... taking only MT=104')
			printtofile(NPt,Etu,signcpoMT104,num_of_displMT104,dpaMT104,104,0,1)
			printtofile(NPt,Etu,signcpoMT104,tot_energy_productsMT104,snhtMT104,104,0,2)
			printtofile(NPt,Etu,signcpoMT104,tot_energy_n_photonsMT104,snhtMT104_EB,104,0,3)
		# -------------------
		if ((ifdisdata105 == 0 or iffldd105 == 0) and ifl105pr == 1):
			dpaMT105 = tmdpaMT105
			snhtMT105 = tmnhtMT105
			snhtMT105_EB = tmnhtMT105_EB
			signcpoMT105 = tmsigncpoMT105
			num_of_displMT105 = tmnum_ofdisplMT105
			tot_energy_productsMT105 = tmtot_energy_productsMT105
			tot_energy_n_photonsMT105 = tmtot_energy_n_photonsMT105
			if (ifdisdata105 == 1):
				print('..... taking only MT=105')
			printtofile(NPt,Etu,signcpoMT105,num_of_displMT105,dpaMT105,105,0,1)
			printtofile(NPt,Etu,signcpoMT105,tot_energy_productsMT105,snhtMT105,105,0,2)
			printtofile(NPt,Etu,signcpoMT105,tot_energy_n_photonsMT105,snhtMT105_EB,105,0,3)
		# -------------------
		if ((ifdisdata106 == 0 or iffldd106 == 0) and ifl106pr == 1):
			dpaMT106 = tmdpaMT106
			snhtMT106 = tmnhtMT106
			snhtMT106_EB = tmnhtMT106_EB
			signcpoMT106 = tmsigncpoMT106
			num_of_displMT106 = tmnum_ofdisplMT106
			tot_energy_productsMT106 = tmtot_energy_productsMT106
			tot_energy_n_photonsMT106 = tmtot_energy_n_photonsMT106
			if (ifdisdata106 == 1):
				print('..... taking only MT=106')
			printtofile(NPt,Etu,signcpoMT106,num_of_displMT106,dpaMT106,106,0,1)
			printtofile(NPt,Etu,signcpoMT106,tot_energy_productsMT106,snhtMT106,106,0,2)
			printtofile(NPt,Etu,signcpoMT106,tot_energy_n_photonsMT106,snhtMT106_EB,106,0,3)
		# -------------------
		if ((ifdisdata107 == 0 or iffldd107 == 0) and ifl107pr == 1):
			dpaMT107 = tmdpaMT107
			snhtMT107 = tmnhtMT107
			snhtMT107_EB = tmnhtMT107_EB
			signcpoMT107 = tmsigncpoMT107
			num_of_displMT107 = tmnum_ofdisplMT107
			tot_energy_productsMT107 = tmtot_energy_productsMT107
			tot_energy_n_photonsMT107 = tmtot_energy_n_photonsMT107
			if (ifdisdata107 == 1):
				print('..... taking only MT=107')
			printtofile(NPt,Etu,signcpoMT107,num_of_displMT107,dpaMT107,107,0,1)
			printtofile(NPt,Etu,signcpoMT107,tot_energy_productsMT107,snhtMT107,107,0,2)
			printtofile(NPt,Etu,signcpoMT107,tot_energy_n_photonsMT107,snhtMT107_EB,107,0,3)
		# --------------------

		# If required take Contributions coming from MT = 5

		# Sometimes the cross sections corresponding to MT = 103, .., 107
		# cannot be obtained because of incomplete disc.+cont. data and
		# absence of MTs = 103, .., 107. Then MT = 5 cross sections along
		# with the yields of charged particles from File 6 are considered.
		# This may give higher predictions of dpa and heating cross sections.
		
		#===========
		# ::NOTE:: n and gamma for heating_EB is not written for this part (May be done later, if needed)
		#===========

		if ((ifdisdata103==1 and iffldd103==0 and ifl103pr==0) or\
			(ifdisdata104==1 and iffldd104==0 and ifl104pr==0) or\
			(ifdisdata105==1 and iffldd105==0 and ifl105pr==0) or\
			(ifdisdata106==1 and iffldd106==0 and ifl106pr==0) or\
			(ifdisdata107==1 and iffldd107==0 and ifl107pr==0)):

			print('', file = ofile_outRMINDD)
			print(':: MESSAGE :: CPO reaction cross section', file = ofile_outRMINDD)
			print('-------------------------------------------------', file = ofile_outRMINDD)
			print('Contributions are taken from MT = 5, since any of', file = ofile_outRMINDD)
			print('the following reactions: (n,p), (n,d), (n,t), ', file = ofile_outRMINDD)
			print('(n,3He), (n,a) cross section(s) is/are incomplete', file = ofile_outRMINDD)

			dpaMT5 = 0
			snhtMT5 = 0
			print( 'MT = 5')

			ifile = open ("tape02", 'r')

			# extraction of cross sections
			NP1 = 0
			while (True):
				line = ifile.readline()
				data = eachlineinfo(line)
				MAT = int(data[6]); MF = int(data[7]); MT  = int(data[8])
				if (MT == 0 and MF != 0):
					line = ifile.readline()
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

				if (MAT != -1):
					if (MF == 3):
						if (MT == 5):
							ZA = ZAv
							AWR = AWRv
							line = ifile.readline()
							data = eachlineinfo(line)
							QM = float(data[0]); QI =  float(data[1]); LR = int(data[3])
							NR = int(data[4]); NP1 = int(data[5])
							E5 = [0]*NP1; sig5 = [0]*NP1
							sdMT5 = [0]*NP1; shMT5 = [0]*NP1; shMT5_EB = [0]*NP1
							Yldp5 = [0]*NP1; Yldd5 = [0]*NP1; Yldtr5 = [0]*NP1; Yld3He5 = [0]*NP1
							YldHe5 = [0]*NP1; Yldarr = [0]*NP1
							num_of_displMT5 = [0]*NP1; tot_energy_productsMT5 = [0]*NP1
							tot_energy_n_photonsMT5 = [0]*NP1
							ifile500.readline()
							i = 0
							while (i<NP1):
								line = ifile.readline()
								data = eachlineinfo(line)
								for j in range(0,6,2):
									if (data[j] != ''):
										E5[i] = float(data[j])
										sig5[i] = float(data[j+1])
										i += 1
									else:
										i += 1
										break
				else:
					break

			ifile.close()

			Z = ZA//1000
			A = AWR
			# get energy-yield data by calling file6 MT5
			(NBTp,INTrp,Nyld,Eyld,Yld) = gtYMf6Mt5()

			if (ifdisdata103 == 1 and iffldd103 == 0 and ifl103pr == 0):
				Eyldp = [0]*int(Nyld[0]); Yldp = [0]*int(Nyld[0])
				Eyldp = Eyld[0][0:Nyld[0]]
				Yldp = Yld[0][0:Nyld[0]]
				NBTpp = NBTp[0][:]
				INTrpp = INTrp[0][:]
				(E5,Yldp5) = yldterpolin(NBTpp,INTrpp,int(Nyld[0]),Eyldp,Yldp,NP1)
				ze[0] = 1; beta[0] = 1
				ifltfr5[0] = 1
				print( '.....for incomplete MT = 103')
				print( '.....contributions taken from MT = 5')
				print('', file = ofile_outRMINDD)
				print('.....for incomplete MT = 103', file = ofile_outRMINDD) 

			if (ifdisdata104 == 1 and iffldd104 == 0 and ifl104pr == 0):
				Eyldd = [0]*int(Nyld[1]); Yldd = [0]*int(Nyld[1])
				Eyldd = Eyld[1][0:Nyld[1]]
				Yldd = Yld[1][0:Nyld[1]]
				NBTpd = NBTp[1][:]
				INTrpd = INTrp[1][:]
				(E5,Yldd5) = yldterpolin(NBTpd,INTrpd,int(Nyld[1]),Eyldd,Yldd,NP1)
				ze[1] = 1; beta[1] = 2
				ifltfr5[1] = 1
				print( '.....for incomplete MT = 104')
				print( '.....contributions taken from MT = 5')
				print('', file = ofile_outRMINDD)
				print('.....for incomplete MT = 104', file = ofile_outRMINDD)

			if (ifdisdata105 == 1 and iffldd105 == 0 and ifl105pr == 0):
				Eyldtr = [0]*int(Nyld[2]); Yldtr = [0]*int(Nyld[2])
				Eyldtr = Eyld[2][0:Nyld[2]]
				Yldtr = Yld[2][0:Nyld[2]]
				NBTpt = NBTp[2][:]
				INTrpt = INTrp[2][:]
				(E5,Yldtr5) = yldterpolin(NBTpt,INTrpt,int(Nyld[2]),Eyldtr,Yldtr,NP1)
				ze[2] = 1; beta[2] = 3
				ifltfr5[2] = 1
				print( '.....for incomplete MT = 105')
				print( '.....contributions taken from MT = 5')
				print('', file = ofile_outRMINDD)
				print('.....for incomplete MT = 105', file = ofile_outRMINDD)

			if (ifdisdata106 == 1 and iffldd106 == 0 and ifl106pr == 0):
				Eyld3He = [0]*int(Nyld[3]); Yld3He = [0]*int(Nyld[3])
				Eyld3He = Eyld[3][0:Nyld[3]]
				Yld3He = Yld[3][0:Nyld[3]]
				NBTp3He = NBTp[3][:]
				INTrp3He = INTrp[3][:]
				(E5,Yld3He5) = yldterpolin(NBTp3He,INTrp3He,int(Nyld[3]),Eyld3He,Yld3He,NP1)
				ze[3] = 2; beta[3] = 3
				ifltfr5[3] = 1
				print( '.....for incomplete MT = 106')
				print( '.....contributions taken from MT = 5')
				print('', file = ofile_outRMINDD)
				print('.....for incomplete MT = 106', file = ofile_outRMINDD)

			if (ifdisdata107 == 1 and iffldd107 == 0 and ifl107pr == 0):
				EyldHe = [0]*int(Nyld[4]); YldHe = [0]*int(Nyld[4])
				EyldHe = Eyld[4][0:Nyld[4]]
				YldHe = Yld[4][0:Nyld[4]]
				NBTpal = NBTp[4][:]
				INTrpal = INTrp[4][:]
				(E5,YldHe5) = yldterpolin(NBTpal,INTrpal,int(Nyld[4]),EyldHe,YldHe,NP1)
				ze[4] = 2; beta[4] = 4
				ifltfr5[4] = 1
				print( '.....for incomplete MT = 107')
				print( '.....contributions taken from MT = 5')
				print('', file = ofile_outRMINDD)
				print('.....for incomplete MT = 107', file = ofile_outRMINDD)

			n2 = A/(A+1)
			n3 = 1.0/(A+1)

			if (LR == 0 or QI == 0):
				QI = QM
				U = abs(QI/n2)

			for  i in range(NP1):
				for inp in range(5):
					if (ifltfr5[inp] == 1):
						if (inp == 1):
							Yldarr = Yldp5
						if (inp == 2):
							Yldarr = Yldd5
						if (inp == 3):
							Yldarr = Yldtr5
						if (inp == 4):
							Yldarr = Yld3He5
						if (inp == 5):
							Yldarr = YldHe5
						n1[inp] = (A+1-beta[inp])/(A+1)
						cbe[inp] = (1.029e+6*ze[inp]*Z)/(beta[inp]**(1/3)+A**(1/3))
						Estar = n1[inp]*E5[i]
						availE = QI + (n2*E5[i])	# if (QI<0)
						#if (QI>0) availE = n2*E5(i)
						if (availE < cbe[inp]):
							Ea = availE
						if (cbe[inp] < availE):
							Ea = cbe[inp]
						f1 = Estar
						f2 = 2*math.sqrt(beta[inp]*Estar*Ea)
						f3 = beta[inp]*Ea 
						dn1 = 0
						dn2 = 0
						dn3 = 0
						if (sig5[i] != 0): # .and. E5(i)>=abs(QI)
							dn1 = abs(Tinteg(Z-ze[inp],A+1-beta[inp],Z,A,Ed,E5[i],QI,\
											beta[inp],mdisp,bad,cad,n3,f1,f2,f3))	
							dn2 = abs(Tintegheat1(A,E5[i],QI,beta[inp],n3,f1,f2,f3))
							dn3 = abs(Ea)
						if (sig5[i] == 0):
							sdMT5[i] = 0
							shMT5[i] = 0
						if (sig5[i] != 0):
							dn1 = dn1*0.8/(2*Ed)
							num_of_displMT5[i] = num_of_displMT5[i] + dn1
							sdMT5[i] = sdMT5[i] + Yldarr[i]*sig5[i]*dn1
							if((dn2+dn3)<availE):
								tot_energy_productsMT5[i] = tot_energy_productsMT5[i] + (dn2+dn3)
								shMT5[i] = shMT5[i] + Yldarr[i]*sig5[i]*(dn2 + dn3)
							if ((dn2+dn3)>availE):
								tot_energy_productsMT5[i] = tot_energy_productsMT5[i] + availE
								shMT5[i] = shMT5[i] + Yldarr[i]*sig5[i]*availE

						if (sdMT5[i]<1e-12):
							sdMT5[i] = 0
						if (shMT5[i]<1e-12):
							shMT5[i] = 0

			sig5 = numpy.asarray(sig5)
			sdMT5 = numpy.asarray(sdMT5)
			shMT5 = numpy.asarray(shMT5)
			num_of_displMT5 = numpy.asarray(num_of_displMT5)
			tot_energy_productsMT5 = numpy.asarray(tot_energy_productsMT5)
			sigetMT5 = trptuqce(E5, sig5, Etu)
			dpaMT5 = trptuqce(E5, sdMT5, Etu)
			snhtMT5 = trptuqce(E5, shMT5, Etu)
			num_of_displ3 = trptuqce(E5, num_of_displMT5, Etu)
			tot_energy_products3 = trptuqce(E5, tot_energy_productsMT5, Etu)

			printtofile(NPt,Etu,sigetMT5,num_of_displ3,dpaMT5,5,0,1)
			printtofile(NPt,Etu,sigetMT5,tot_energy_products3,snhtMT5,5,0,2)

		# if required, then take contributions from MT = 5

		# -------------------

 	# for the iflMTtppr

	# --------------------------------

	for isum in range(NPt):
		dpa3[isum] = dpa3[isum] + dpaMT103[isum] + dpaMT104[isum] + \
		dpaMT105[isum] + dpaMT106[isum] + dpaMT107[isum] + dpaMT5[isum]

		snht3[isum] = snht3[isum] + snhtMT103[isum] + snhtMT104[isum] + \
		snhtMT105[isum] + snhtMT106[isum] + snhtMT107[isum] + snhtMT5[isum]

		snht3_EB[isum] = snht3_EB[isum] + snhtMT103_EB[isum] + snhtMT104_EB[isum] + \
		snhtMT105_EB[isum] + snhtMT106_EB[isum] + snhtMT107_EB[isum]

		signcpo3[isum] = signcpo3[isum] + signcpoMT103[isum] + signcpoMT104[isum] + \
		signcpoMT105[isum] + signcpoMT106[isum] + signcpoMT107[isum]
		
		num_of_displ3[isum] = num_of_displ3[isum] + num_of_displMT103[isum] + \
		num_of_displMT104[isum] + num_of_displMT105[isum] + num_of_displMT106[isum] + \
		num_of_displMT107[isum]
		
		tot_energy_products3[isum] = tot_energy_products3[isum] + tot_energy_productsMT103[isum] + \
		tot_energy_productsMT104[isum] + tot_energy_productsMT105[isum] + tot_energy_productsMT106[isum] + \
		tot_energy_productsMT107[isum]

		tot_energy_n_photons3[isum] = tot_energy_n_photons3[isum] + tot_energy_n_photonsMT103[isum] + \
		tot_energy_n_photonsMT104[isum] + tot_energy_n_photonsMT105[isum] + tot_energy_n_photonsMT106[isum] + \
		tot_energy_n_photonsMT107[isum]

	# THIS PART IS ADDED TO REMOVE SUDDEN NON-ZERO VALUES BELOW THE 
	# STARTING ENERGY POINT
	# --------------

	for i in reversed(range (len(Etu))):
		if (dpa3[i] == 0 and (dpa3[i-1] > 0 and dpa3[i+1] > 0)):
			for j in range(i):
				dpa3[j] = 0
			break

	for i in reversed(range (len(Etu))):
		if (snht3[i] == 0 and (snht3[i-1] > 0 and snht3[i+1] > 0)):
			for j in range(i):
				snht3[j] = 0
			break
		if (snht3_EB[i] == 0 and (snht3_EB[i-1] > 0 and snht3_EB[i+1] > 0)):
			for j in range(i):
				snht3_EB[j] = 0
			break
	# --------------

	printtofile(NPt,Etu,signcpo3,num_of_displ3,dpa3,3001,0,1)
	printtofile(NPt,Etu,signcpo3,tot_energy_products3,snht3,3001,0,2)
	printtofile(NPt,Etu,signcpo3,tot_energy_n_photons3,snht3_EB,3001,0,3)

#=======Average Eg2 from Continuum and discrete (File 6)=======*	

	# To find Egamma square from continuum of File 15 and 
	# from discrete plus continuum of File 6

	# Used to calculate the average value of squared Eg,
	# needed in the estimation of recoil energy from radiative capture reaction

def avegsqc(Ep1,fEEp,nfe,iplaw,ndisc,Yd,iFile):
	nquad = 64
	xabc = [0]*nquad; wg = [0]*nquad
	(xabc, wg) = GQ()
	s1 = 0
	s2 = 0
	s3 = 0
	s4 = 0
	if (ndisc != 0 and iFile == 6):
		for i in range(ndisc):
			s1 = s1 + Ep1[i]*Ep1[i]*fEEp[i]
			s2 = s2 + fEEp[i]
	eru = Ep1[ndisc+1]
	fu = fEEp[ndisc+1]
	for i in range(ndisc+2, nfe):
		erl = eru
		eru = Ep1[i]
		fl = fu
		fu = fEEp[i]
		de = eru-erl
		if (de != 0):
			for j in range(nquad):
				eg = erl+(1+xabc[j])*de/2 
				fg = TERPOLIN(iplaw,eg,erl,eru,fl,fu)
				p2 = eg*eg
				s3 = s3 + wg[j]*p2*fg*de
				s4 = s4 + wg[j]*fg*de

	# matching with recoil energy from NJOY21 
	# for W180, W182 TENDL-2015, ENDF/B-VII.1
	if (iFile == 6):
		value_to_return = (s1+s3)/(s2+s4)
	if (iFile==15):
		value_to_return = s3*Yd/s4

	return(value_to_return)

#=======Radiative Capture=======*

	## Calculation of neutron dpa and heating cross sections due to
	## radiative capture of neutron reaction.

def RADIATIVE_CAPTURE (ofile_outRMINDD,NPt,Etu,mdisp,Ed,bad,cad):
	siget = numpy.zeros(NPt); sdpa = numpy.zeros(NPt); snht = numpy.zeros(NPt)
	s1 = numpy.zeros(NPt); s2 = numpy.zeros(NPt)
	num_of_displ = numpy.zeros(NPt); tot_energy_products = numpy.zeros(NPt)

	# file 12
	NBT = numpy.zeros(20); INTr = numpy.zeros(20)
	# File 12 -- only the LO = 1 (MULTIPLICITIES) option is included
	# File 12 -- LO = 2 (TRANSITION PROBABILITY ARRAYS) option is not included.

	# file 6
	NPg6 = numpy.zeros(2000); ND6 = numpy.zeros(2000); NBT6 = numpy.zeros(20);INTr6 = numpy.zeros(20)
	Eg6 = numpy.zeros((2000,1000)); g6 = numpy.zeros((2000,1000))
	eu = numpy.zeros(1000); gu = numpy.zeros(1000)
	E6 = numpy.zeros(2000); Y6 = numpy.zeros(2000); En6 = numpy.zeros(2000)

	# file 15
	NBT15a = numpy.zeros(20); INTr15a = numpy.zeros(20); NPg15 = numpy.zeros(200)
	Eg15 = numpy.zeros((200,200)); g15 = numpy.zeros((200,200))
	E15 = numpy.zeros(100); Y15 = numpy.zeros(100); En15 = numpy.zeros(200)

	print("n, g .....")

	print('', file = ofile_outRMINDD)
	print(' Radiative capture dpa/heating cross section', file = ofile_outRMINDD)
	print('------------------------------------------------', file = ofile_outRMINDD)

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
			data = eachlineinfo(line)
			iflspace = 0
			for element in data:
				if (element == ''):
					iflspace = 1
					break
			if (iflspace == 0):
				ZAv = float(data[0]); AWRv = float(data[1]); L0 = int(data[2]);\
				LCT = int(data[3]); NKv = int(data[4]); L2 = int(data[5]);\
				MAT = int(data[6]); MF = int(data[7]); MT = int(data[8])

		if (MAT != -1):
			if (MF == 6):
				if (MT == 102):
					ifl6 = 1
					NK6 = NKv
					line = ifile.readline()
					data = eachlineinfo(line)
					C1 = float(data[0]); C2 =  float(data[1]); LIP = int(data[2])
					LAW = int(data[3]); NR = int(data[4]); NP6 = int(data[5]); MAT = int(data[6])
					MF = int(data[7]); MT = int(data[8])
					N = 0
					while (N < NR):
						line = ifile.readline()
						data = eachlineinfo(line)
						for i in range(0,6,2):
							if (data[i] != ''):
								NBT[N] = int(data[i])
								INTr[N] = int(data[i+1])
								N += 1
							else:
								N += 1
								break
					N = 0
					while (N < NP6):
						line = ifile.readline()
						data = eachlineinfo(line)
						for i in range(0,6,2):
							if (data[i] != ''):
								E6[N] = float(data[i])
								Y6[N] = float(data[i+1])
								N += 1
							else:
								N += 1
								break
					Y6tot = [0]*NPt
	
	# Interpolating to find the total yields corresponding to energy
	# points in unique energy array
	
					for i in range(NPt):
						for j in range(NP6):
							if (Etu[i] == E6[j]):
								Y6tot[i] = Y6[j]
							if (E6[j]<Etu[i] and Etu[i]<E6[j+1]):
								for k in range(NR):
									if (j <= NBT[k]):
										intflg = INTr[k]
								Y6tot[i] = TERPOLIN(intflg,Etu[i],E6[j],E6[j+1],Y6[j],Y6[j+1])
					
					line = ifile.readline()
					data = eachlineinfo(line)
					C1 = float(data[0]); C2 =  float(data[1]); LANG = int(data[2])
					LEP = int(data[3]); NR6 = int(data[4]); NE6 = int(data[5]); MAT = int(data[6])
					MF = int(data[7]); MT = int(data[8])
					N = 0
					while (N < NR6):
						line = ifile.readline()
						data = eachlineinfo(line)
						for i in range(0,6,2):
							if (data[i] != ''):
								NBT6[N] = int(data[i])
								INTr6[N] = int(data[i+1])
								N += 1
							else:
								N += 1
								break
					for i in range(NE6):
						line = ifile.readline()
						data = eachlineinfo(line)
						C1 = float(data[0]); En6[i] =  float(data[1]); ND6[i] = int(data[2])
						NA = int(data[3]); NW = int(data[4]); NPg6[i] = int(data[5]); MAT = int(data[6])
						MF = int(data[7]); MT = int(data[8])
						N = 0
						while (N < NPg6[i]):
							line = ifile.readline()
							data = eachlineinfo(line)
							for j in range(0,6,2):
								if (data[j] != ''):
									Eg6[i][N] = float(data[j])
									g6[i][N] = float(data[j+1])
									N += 1
								else:
									N += 1
									break
		else:
			break
	ifile.close()

	ifl12 = 0
	ifile = open("tape01", 'r')
	while True:
		line = ifile.readline()
		if (line == ''):
			break
		data = eachlineinfo(line)
		MAT = int(data[6]); MF = int(data[7]); MT  = int(data[8])
		if (MF <= 12 and MT == 0):
			if (MAT != -1):
				line = ifile.readline()
				data = eachlineinfo(line)
				iflspace = 0
				for element in data:
					if (element == ''):
						iflspace = 1
						break
				if (iflspace == 0):
					ZAv = float(data[0]); AWRv = float(data[1]); Lo = int(data[2]);\
					L1 = int(data[3]); NKv = int(data[4]); L2 = int(data[5]);\
					MAT = int(data[6]); MF = int(data[7]); MT = int(data[8])

		if (MF == 0):
			if (MAT != -1):
				line = ifile.readline()
				data = eachlineinfo(line)
				iflspace = 0
				for element in data:
					if (element == ''):
						iflspace = 1
						break
				if (iflspace == 0):
					ZAv = float(data[0]); AWRv = float(data[1]); Lo = int(data[2]);\
					L1 = int(data[3]); NKv = int(data[4]); L2 = int(data[5]);\
					MAT = int(data[6]); MF = int(data[7]); MT = int(data[8])

		if (MAT != -1):
			if (MF == 12):
				if (MT == 102):
					ifl12 = 1
					NK = NKv
					AWR = AWRv
					ZA = ZAv
					NPn12 = [0]*NK; En12 = numpy.zeros((NK,200)); Yk12 = numpy.zeros((NK,200))
					EGk = [0]*NK; ESk = [0]*NK; LP = [0]*NK; NBTk12 = numpy.zeros((NK,20))
					INTrk12 = numpy.zeros((NK,20)); NRk12 = [0]*NK
					line = ifile.readline()
					data = eachlineinfo(line)
					C1 = float(data[0]); C2 =  float(data[1]); L1 = int(data[2])
					L2 = int(data[3]); NR = int(data[4]); NP12 = int(data[5]); MAT = int(data[6])
					MF = int(data[7]); MT = int(data[8])
					E12 = [0]*NP12; Y12 = [0]*NP12
					N = 0
					while (N < NR):
						line = ifile.readline()
						data = eachlineinfo(line)
						for i in range(0,6,2):
							if (data[i] != ''):
								NBT[N] = int(data[i])
								INTr[N] = int(data[i+1])
								N += 1
							else:
								N += 1
								break
					N = 0
					while (N < NP12):
						line = ifile.readline()
						data = eachlineinfo(line)
						for i in range(0,6,2):
							if (data[i] != ''):
								E12[N] = float(data[i])
								Y12[N] = float(data[i+1])
								N += 1
							else:
								N += 1
								break
					Y12tot = [0]*NPt
	
	# Interpolating to find the total yields corresponding to energy
	# points in unique energy array
	
					for i in range(NPt):
						for j in range(NP12):
							if (Etu[i] == E12[j]):
								Y12tot[i] = Y12[j]
								break
							if (E12[j]<Etu[i] and Etu[i]<E12[j+1]):
								for k in range(NR):
									if (j<=NBT[k]):
										intflg = INTr[k]
										break
								Y12tot[i] = TERPOLIN(intflg,Etu[i],E12[j],E12[j+1],Y12[j],Y12[j+1])
								break
					
					if (NK == 1):
						for N in range(NP12):
							Yk12[NK][N] = Y12[N]
					if (NK > 1):
						for i in range(NK):
							line = ifile.readline()
							data = eachlineinfo(line)
							EGk[i] = float(data[0]); ESk[i] =  float(data[1]); LP[i] = int(data[2])
							LF = int(data[3]); NRk12[i] = int(data[4])
							NPn12[i] = int(data[5]); MAT = int(data[6])
							MF = int(data[7]); MT = int(data[8])
							N = 0
							while (N < NRk12[i]):
								line = ifile.readline()
								data = eachlineinfo(line)
								for j in range(0,6,2):
									if (data[j] != ''):
										NBTk12[i][N] = int(data[j])
										INTrk12[i][N] = int(data[j+1])
										N += 1
									else:
										N += 1
										break
							N = 0
							while (N < NPn12[i]):
								line = ifile.readline()
								data = eachlineinfo(line)
								for j in range(0,6,2):
									if (data[j] != ''):
										En12[i][N] = float(data[j])
										Yk12[i][N] = float(data[j+1])
										N += 1
									else:
										N += 1
										break
		
	# Yield of discrete and continuum g for all neutron energy 
	# (unique energy array) by interpolation

					if (ifl12 == 1):
						Yk12tot = numpy.zeros((NK,NPt))
						if (NK == 1):
							for j in range(NPt):
								Yk12tot[NK][j] = Y12tot[j]
						if (NK > 1):
							for i in range(NK):
								for j in range(NPt):
									for k in range(NPn12[i]):
										if (Etu[j] == En12[i][k]):
											Yk12tot[i][j] = Yk12[i][k]
											break
										if (En12[i][k]<Etu[j] and Etu[j]<En12[i][k+1]):
											for l in range(NRk12[i]):
												if (k <= NBTk12[i][l]):
													intflg = INTrk12[i][l]
													break
											Yk12tot[i][j] = TERPOLIN(intflg,Etu[j],En12[i][k],En12[i][k+1],Yk12[i][k],Yk12[i][k+1])
											break
		else:
			break
	ifile.close()

	ifl15 = 0
	ifile = open ("tape01", 'r')
	if (ifl6 == 0):
		while True:
			line = ifile.readline()
			if (line == ''):
				break
			data = eachlineinfo(line)
			MAT = int(data[6]); MF = int(data[7]); MT  = int(data[8])
			if (MT == 0):
				line = ifile.readline()
				(ZAv,AWRv,L0,L1,NCv,L2,MAT,MF,MT) = line_type1_info(line)
	
			if (MAT != -1):
				if (MF == 15):
					if (MT == 102):
						ifl15 = 1
						NC = NCv
						line = ifile.readline()
						(C1,C2,L1,LF,NR,NP15,MAT,MF,MT) = line_type2_info(line)
						N = 0
						while (N < NR):
							line = ifile.readline()
							data = eachlineinfo(line)
							for j in range(0,6,2):
								if (data[j] != ''):
									NBT[N] = int(data[j])
									INTr[N] = int(data[j+1])
									N += 1
								else:
									N += 1
									break
						N = 0
						while (N < NP15):
							line = ifile.readline()
							data = eachlineinfo(line)
							for j in range(0,6,2):
								if (data[j] != ''):
									E15[N] = float(data[j])
									Y15[N] = float(data[j+1])
									N += 1
								else:
									N += 1
									break
						line = ifile.readline()
						(C1,C2,L1,L2,NR15a,NE15,MAT,MF,MT) = line_type2_info(line)
						N = 0
						while (N < NR15a):
							line = ifile.readline()
							data = eachlineinfo(line)
							for j in range(0,6,2):
								if (data[j] != ''):
									NBT15a[N] = int(data[j])
									INTr15a[N] = int(data[j+1])
									N += 1
								else:
									N += 1
									break
						for i in range(NE15):
							line = ifile.readline()
							(C1,En15[i],L1,L2,NR,NPg15[i],MAT,MF,MT) = line_type2_info(line)
							N = 0
							while (N < NR):
								line = ifile.readline()
								data = eachlineinfo(line)
								for j in range(0,6,2):
									if (data[j] != ''):
										NBT[N] = int(data[j])
										INTr[N] = int(data[j+1])
										N += 1
									else:
										N += 1
										break
							N = 0
							while (N < NPg15[i]):
								line = ifile.readline()
								data = eachlineinfo(line)
								for j in range(0,6,2):
									if (data[j] != ''):
										Eg15[i][N] = float(data[j])
										g15[i][N] = float(data[j+1])
										N += 1
									else:
										N += 1
										break
			else:
				break
		ifile.close()

# ---- Reading the (n,g) cross sections ----

	ifile = open ("tape02", 'r')
	ifile.readline()
	line = ifile.readline()
	(ZA,AWR,L0,L1,L2,L3,MAT,MF,MT) = line_type1_info(line)
	while True:
		line = ifile.readline()
		if (line == ''):
			break
		data = eachlineinfo(line)
		MAT = int(data[6]); MF = int(data[7]); MT = int(data[8])
		if (MAT != -1):
			if (MF == 3):
				if (MT == 102):
					line = ifile.readline()
					data = eachlineinfo(line)
					QM = float(data[0]); QI =  float(data[1]); LR = int(data[3])
					NR = int(data[4]); NP = int(data[5])
					ifile.readline()
					E = [0]*NP; sig = [0]*NP
					i = 0 
					while (i < NP):
						line = ifile.readline()
						data = eachlineinfo(line)
						for j in range(0,6,2):
							if (data[j] != ''):
								E[i] = float(data[j])
								sig[i] = float(data[j+1])
								i += 1
							else:
								i += 1
								break
		else:
			break
	ifile.close()

	print('', file = ofile_outRMINDD)
	print('Number of energy points given is ',NP, file = ofile_outRMINDD)

#-------------------------------------------------------------
	# Finding basic cross sections in unique energy array because 
	# there can repetition of energy points in MT = 102 
	
	k = 0
	for i in range(NPt):
		for j in range(k,NP):
			if (Etu[i] == E[j]):
				siget[i] = sig[j]
				break
			if (Etu[i] > E[j] and Etu[i] <= E[j+1]):
				if (E[j] == E[j+1]):
					siget[i] = sig[j+1]
					k = j+1
				else:
					siget[i] = crstd(Etu[i],E[j],E[j+1],sig[j],sig[j+1])
					k = j
				break
	
	if (ifl6 == 1):
		print('', file = ofile_outRMINDD)
		print('Emitted photon data are given in File 6', file = ofile_outRMINDD)
		print('', file = ofile_outRMINDD)
		print('Neutron energy/ discrete gamma/ continuum gamma', file = ofile_outRMINDD)
		for i in range(0,NE6,5):
			print('', file = ofile_outRMINDD)
			print(En6[i],' ',ND6[i],' ',NPg6[i]-ND6[i], file = ofile_outRMINDD)
	
	if (ifl12 == 1):
		print('', file = ofile_outRMINDD)
		print('Emitted photon data are given in File 12', file = ofile_outRMINDD)

	if (ifl15 == 1):
		print('', file = ofile_outRMINDD)
		print('Emitted photon data are given in File 15 ', file = ofile_outRMINDD)
		
	# ENERGY OF THE EMITTED DISCRETE PHOTON FROM FILE 12
	
	Eg2Av = [0]*NPt
	if ((ifl12 == 1 and NK > 1) or (ifl12==1 and NK==1 and ifl15==0)):
		EGkp = numpy.zeros((NK,NPt))
		NKd = NK
		if (ifl15 == 1):
			NKd = NK-1
		for i in range(NKd):
			for j in range(NPt):
				if (LP[i] == 2):
					EGkp[i][j] = EGk[i] + (AWR*Etu[j]/(AWR+1))
				else:
					EGkp[i][j] = EGk[i]
		print('', file = ofile_outRMINDD)
		print('Total number of emitted gamma sections', file = ofile_outRMINDD)
		print(NKd,' discrete and',NK-NKd,' continuum', file = ofile_outRMINDD)

	# For each neutron energy, Total yield over all NK contributions
	# must be normalized to Y12tot. Total yields are given only
	# when NK > 1

	if (ifl12 == 1 and NK > 1):
		for j in range(NPt):
			sY12 = 0
			for i in range(NK):
				sY12 = sY12 + Yk12tot[i][j]
			for i in range(NK):
				Yk12tot[i][j] = Yk12tot[i][j] * Y12tot[j]/sY12

 	# AVERAGE OF (Egamma)SQUARE FROM FILE 6 AND FILE 15

	if (ifl15 == 1):
		Eg15t = numpy.zeros((NPt,1000)); g15t = numpy.zeros((NPt,1000)); ntm = [0]*NPt
		for i in range(NPt):
			if (Yk12tot[NK][i] != 0):
				for j in range(NE15):
					if(Etu[i] == En15[j]):
						for k in range(NPg15[j]):
							Eg15t[i][k] = Eg15[j][k]
							g15t[i][k] = g15[j][k]
						ntm[i] = NPg15[j]
						break
					if (En15[j]<Etu[i] and Etu[i]<En15[j+1]):
						for k1 in range (NR15a):
							if (j<=NBT15a[k1]):
								iplaw = INTr15a[k1]
								break

						#if (Etu(i)<1.0e+6) iplaw = 1
						diff1 = Etu[i] - En15[j]		
						diff2 = Etu[i] - Etu[i-1]
						if (diff1 <= diff2):
							for k in range(NPg15[j]):
								x = Etu[i]
								x1 = En15[j]
								x2 = En15[j+1]
								y1 = Eg15[j][k]
								y2 = Eg15[j+1][k]
								y11 = g15[j][k]
								y22 = g15[j+1][k]
								Eg15t[i][k] = TERPOLIN(iplaw,x,x1,x2,y1,y2)
								g15t[i][k] = TERPOLIN(iplaw,x,x1,x2,y11,y22)
						if (diff2 < diff1):
							for k in range(NPg15[j]):
								x = Etu[i]
								x1 = Etu[i-1]
								x2 = En15[j+1]
								y1 = Eg15t[i-1][k]
								y2 = Eg15[j+1][k]
								y11 = g15t[i-1][k]
								y22 = g15[j+1][k]
								Eg15t[i][k] = TERPOLIN(iplaw,x,x1,x2,y1,y2)
								g15t[i][k] = TERPOLIN(iplaw,x,x1,x2,y11,y22)
						
						ntm[i] = NPg15[j]
						break
		
		for i in range(NPt):
			if (Yk12tot[NK][i] != 0):
				eu = 0
				gu = 0
				for j in range(ntm[i]):
					eu[j] = Eg15t[i][j]
					gu[j] = g15t[i][j]
				s1[i] = avegsqc(eu,gu,ntm[i],1,0,Yk12tot[NK][i],15)


	if (ifl6 == 1):
		Eg6t = numpy.zeros((NPt,1000)); g6t = numpy.zeros((NPt,1000)); ntm = [0]*NPt; ntmd = [0]*NPt
		for i in range(NPt):
			for j in range(NE6):
				if(Etu[i] == En6[j]):
					for k in range(int(NPg6[j])):
						Eg6t[i][k] = Eg6[j][k]
						g6t[i][k] = g6[j][k]
					ntm[i] = NPg6[j]
					ntmd[i] = ND6[j]
					break

			# histogram interpolations are used to validate 
			# recoil energy (e.g. in W182 from ENDF/B-VII.1)
				
				if (En6[j]<Etu[i] and Etu[i]<En6[j+1]):
					for k1 in range(NR6):
						if (j<=int(NBT6[k1])):
							iplaw = INTr6[k1]
							break
					# if (Etu(i)<1.0e+6) iplaw = 1
					diff1 = Etu[i] - En6[j]		
					diff2 = Etu[i] - Etu[i-1]
					if (diff1 <= diff2):
						for k in range(int(NPg6[j])):
							x = Etu[i]
							x1 = En6[j]
							x2 = En6[j+1]
							y1 = Eg6[j][k]
							y2 = Eg6[j+1][k]
							y11 = g6[j][k]
							y22 = g6[j+1][k]
							Eg6t[i][k] = TERPOLIN(iplaw,x,x1,x2,y1,y2)
							g6t[i][k] = TERPOLIN(iplaw,x,x1,x2,y11,y22)
					
					if (diff2 < diff1):
						for k in range(int(NPg6[j])):
							x = Etu[i]
							x1 = Etu[i-1]
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
		
		for i in range(NPt):
			ncont = int(ntm[i] - ntmd[i])
			ndisc = int(ntmd[i])
			for j in range(int(ntm[i])):
				eu[j] = Eg6t[i][j]
				gu[j] = g6t[i][j]
			s1[i] = avegsqc(eu,gu,int(ntm[i]),LEP,ndisc,Y6tot[i],6)

	if ((ifl12 == 1 and NK > 1) or (ifl12 == 1 and NK == 1 and ifl15 == 0)):
		NKd = NK
		if (ifl15 == 1):
			NKd = NK-1
		for i in range(NPt):
			s2[i] = 0
			for j in range(NKd):
				s2[i] = s2[i] + (Yk12tot[j][i]*EGkp[j][i]*EGkp[j][i])
		
	# Total Egamma square from discrete and continuum
	
	for i in range(NPt):
		Eg2Av[i] = s1[i] + s2[i]
 
	Z = int(ZA/1000)
	A1 = AWR + 1
	A1f = 1/A1
	if (LR == 0):
		QM = abs(QI)

 	# POINT DPA/heating CROSS SECTION

	emc2 = 939.512e6 
	tm = emc2*A1
	rtm = 1/tm

	for i in range(NPt):
		Er = Eg2Av[i]*rtm/2
		dn1 = ADisM(Z,AWR+1,Z,AWR,Er,Ed,mdisp,bad,cad)*0.8/(2*Ed)
		dn2 = Er + Etu[i]*A1f
		num_of_displ[i] = dn1
		tot_energy_products[i] = dn2
		if (ifl12 == 0 and ifl6==0 and ifl15==0):
			Ea = Etu[i]*A1f
			Ea1 =  QM  			# AWR*Ea +		
			Er = Ea1*Ea1*rtm/2				
			dn1 =  ADisM(Z,AWR+1,Z,AWR,Er,Ed,mdisp,bad,cad)*0.8/(2*Ed)
			dn2 = Er
			num_of_displ[i] = dn1
			tot_energy_products[i] = dn2

		if (siget[i] == 0):
			sdpa[i] = 0
		if (siget[i] == 0):
			snht[i] = 0
		if (siget[i] != 0):
			sdpa[i] = abs(siget[i]*dn1)
			snht[i] = abs(siget[i]*dn2)

	printtofile(NPt,Etu,siget,num_of_displ,sdpa,102,0,1)
	printtofile(NPt,Etu,siget,tot_energy_products,snht,102,0,2)

#=========== ELASTIC INTERACTION ===========*
		
    # Calculation of neutron dpa and heating cross sections due to the elastic
	# scattering interactions.
	
def ELASTIC_SCATTERING(ofile_outRMINDD,NPt,Etu,mdisp,Ed,bad,cad):

	siget = [0]*NPt; sdpat = [0]*NPt; snht = [0]*NPt
	alfull = numpy.zeros((NPt,65)); fmuE = numpy.zeros((NPt,64))
	num_of_displ = numpy.zeros(NPt); tot_energy_products = numpy.zeros(NPt)
	nquad = 64
	xabc = [0]*nquad; wg = [0]*nquad; fpr = [0]*nquad
# ------------------------------------------------------------------		
	print( "n, n .....")
	print('', file = ofile_outRMINDD)
	print(' Elastic scattering dpa/heating cross section', file = ofile_outRMINDD)
	print('------------------------------------------------', file = ofile_outRMINDD)

# ------------------------------------------------------------------

	# Extraction of Elastic Cross Sections
	ifile = open("tape02", 'r')
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
					LR = int(ifile.readline().split()[3])
					(E, sig) = line_type3_info(ifile,NP,2)
		else:
			break
	ifile.close()

	print('', file = ofile_outRMINDD)
	print('Number of energy points given is ',NP, file = ofile_outRMINDD)

 	# Extraction of Legendre Polynomial Coefficients and 
	# Tabulated Probability

	ifile = open("tape01", 'r')
	while True:
		line = ifile.readline()
		if (line == ''):
			break
		data = eachlineinfo(line)
		MAT = int(data[6]); MF = int(data[7]); MT  = int(data[8])
		if (MF == 3 and MT == 0):
			ifile.readline()
			line = ifile.readline()
			data = eachlineinfo(line)
			iflspace = 0
			for element in data:
				if (element == ''):
					iflspace = 1
					break
			if (iflspace == 0):
				ZA = float(data[0]); AWR = float(data[1]); l1 = int(data[2]);\
				LTT = int(data[3]); NK = int(data[4]); l2 = int(data[5]);\
				MAT = int(data[6]); MF = int(data[7]); MT = int(data[8])
		if (MAT	!= -1):
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
						al = numpy.zeros((NE1,65)); EL = [0]*NE1; alc = [0]*65
						for i in range(NE1):
							line = ifile.readline()
							data = eachlineinfo(line)
							T = float(data[0]); EL[i] = float(data[1]); NL = int(data[4])
							NL = NL+1
							al[i][0] = 1
							temporary = [0]*NL
							temporary = line_type3_info(ifile,NL-1,1)
							for j, value in enumerate(temporary, 1):
								al[i][j] = value
					if (LTT == 3 or LTT == 2):
						line = ifile.readline()
						data = eachlineinfo(line)
						NE2 = int(data[5])
						# Tabulated Probability
						ifile.readline()
						Enf = [0]*NE2; NPr = [0]*NE2; cdata = numpy.zeros((NE2,201))
						fdata = numpy.zeros((NE2,201)); ftotal = numpy.zeros((NE2,64))
						for i in range(NE2):
							line = ifile.readline()
							data = eachlineinfo(line)
							T = float(data[0]); Enf[i] = float(data[1])
							line = ifile.readline()
							data = eachlineinfo(line)
							NPr[i] = int(data[0])
							temporary1 = [0]*201
							temporary2 = [0]*201
							(temporary1,temporary2) = line_type3_info(ifile,NPr[i],2)
							for j, (value1, value2) in enumerate(zip(temporary1, temporary2), 0):
								cdata[i][j] = value1
								fdata[i][j] = value2
		else:
			break
	ifile.close()

	print('', file = ofile_outRMINDD)
	if(LTT == 3):
		print('Legendre coefficients and tabulated probability', file = ofile_outRMINDD)
		print('data representations of angular distribution', file = ofile_outRMINDD)
	if(LTT == 1):
		print('Legendre coefficients representation of angular', file = ofile_outRMINDD)
		print('distribution', file = ofile_outRMINDD)
	if(LTT == 2):
		print('Tabulated probability data representation of', file = ofile_outRMINDD)
		print('angular distribution', file = ofile_outRMINDD)
	if(LTT == 0):
		print('All energy scattering are isotropic', file = ofile_outRMINDD)

#-------------------------------------------------------------------------
	A = AWR
	Z = int(ZA/1000)
	y = 4*A/((A+1)**2)
	nquad = 64
	(xabc, wg) = GQ()	# Gauss-quadrature values and weights
#-------------------------------------------------------------------------		
	# Finding basic cross sections in unique energy array because 
	# there can repetition of energy points in MT=2 
	
	k = 0
	for i in range(NPt):
		for j in range(k, NP):
			if (Etu[i] == E[j]):
				siget[i] = sig[j]
				break
			if (Etu[i] > E[j] and Etu[i] <= E[j+1]):
				if (E[j] == E[j+1]):
					siget[i] = sig[j+1]
					k = j+1
				else:
					siget[i] = crstd(Etu[i],E[j],E[j+1],sig[j],sig[j+1])
					k = j
				break

	# LOG- LINEAR INTERPOLATION BETWEEN MU AND F(MU,E) 

	for i in range(NE2):
		for j in range(64):
			for k in range(NPr[i]):
				if (xabc[j] == cdata[i][k]):
					ftotal[i][j] = fdata[i][k]
					break
				if(xabc[j]>cdata[i][k] and xabc[j]<cdata[i][k+1]):
					x = xabc[j]
					x1 = cdata[i][k]
					x2 = cdata[i][k+1]
					y1 = fdata[i][k]  
					y2 = fdata[i][k+1]        
					ftotal[i][j] = y1*math.exp((x-x1)*math.log(y2/y1)/(x2-x1))
					break
 		
	# TO GET THE F(MU,E) FOR THE FULL ENERGY RANGE
		
	for i in range(NPt):
		for j in range(64):
			fmuE[i][j] = 0.5
		
	for i in range(NPt):
		for j in range(NE2):
			if (Etu[i] == Enf[j]):
				for k in range(64):
					fmuE[i][k] = ftotal[j][k]
				exit
			if (Etu[i] > Enf[j] and Etu[i] < Enf[j+1]):
				x = Etu[i]
				x1 = Enf[j]
				x2 = Enf[j+1]
				for k in range(64):
					y1 = ftotal[j][k]
					y2 = ftotal[j+1][k]
					fmuE[i][k] = (y1 + ((y2-y1)*(x-x1)/(x2-x1)))
				exit

	# TO GET THE AL(E) FOR THE FULL ENERGY RANGE

	if (LTT == 0):		# isotropic angular distribution
		for i in range(NPt):
			alfull[i][0] = 1
			for k in range(1, 65):
				alfull[i][k] = 0
	if (LTT != 0):
		for i in range(NPt):
			alfull[i][0] = 1
			for k in range(1, 65):	# the number 65 has to be changed to 60 everywhere in the following for Ti48 JENDL-4.0
				alfull[i][k] = 0
	
		for i in range(NPt):
			for j in range(NE1-1):
				if (Etu[i] == EL[j]):
					for k in range(1, 65):
						alfull[i][k] = al[j][k]
					break
				if (Etu[i] > EL[j] and Etu[i] <= EL[j+1]):
					diff1 = Etu[i] - EL[j]
					diff2 = Etu[i] - Etu[i-1]
					if (diff1 <= diff2):
						for k in range(1, 65):
							x = Etu[i]
							x1 = EL[j]
							x2 = EL[j+1]
							y1 = al[j,k]
							y2 = al[j+1,k]
							alfull[i][k] = y1+((y2-y1)*(x-x1)/(x2-x1))

					if (diff2 < diff1):
						for k in range(1, 65):
							x = Etu[i]
							x1 = Etu[i-1]
							x2 = EL[j+1]
							y1 = alfull[i-1][k]
							y2 = al[j+1][k]
							alfull[i][k] = y1+((y2-y1)*(x-x1)/(x2-x1))

					for k in range(1, 65):
						if (alfull[i][k] == 0):
							alfull[i][k] = alfull[i-1][k]
					break
	# if for LTT != 0

	NLa = 65
	if (LTT == 0):
		NLa = 3

	# Average damage energy is collected in dn1 
	# (for dpa cross sections)
	# Average recoil energy is collected in dn2
	# (for heating cross sections)

	for i in range(NPt):
		if((LTT == 3 and Etu[i]<=EL[NE1-1]) or LTT==1 or LTT==0):
			for k in range(NLa):
				if (abs(alfull[i][k])>1):
					alfull[i][k] = 0			# TO AVOID VERY LARGE NUMBERS
				if (abs(alfull[i][k])<1.0E-12):
					alfull[i][k] = 0	# TO AVOID VERY SMALL NUMBERS
				if (abs(alfull[i][k])>0):
					alc[k] = alfull[i][k]	# TO AVOID NaN
			dn1 = abs(discrec1(alc,NLa,Etu[i],0.0,1.0,Z,A,Z,A,Ed,mdisp,bad,cad))*0.8/(2*Ed)
			dn2 = abs(discrec3(alc,NLa,Etu[i],0.0,1.0,A))
			num_of_displ[i] = dn1
			tot_energy_products[i] = dn2
			
		if((LTT == 3 and Etu[i] > Enf[0]) or LTT == 2):
			for j in range(64):
				fpr[j] = fmuE[i][j]

			dn1 = abs(discrec2(fpr,Etu[i],0.0,1.0,Z,A,Z,A,Ed,mdisp,bad,cad))*0.8/(2*Ed)

			dn2 = abs(discrec4(fpr,Etu[i],0.0,1.0,A))
			num_of_displ[i] = dn1
			tot_energy_products[i] = dn2

		sdpat[i] = siget[i]*dn1
		snht[i] = siget[i]*dn2

	printtofile (NPt,Etu,siget,num_of_displ,sdpat,2,0,1) 
	printtofile (NPt,Etu,siget,tot_energy_products,snht,2,0,2)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# NOT USED	
# Function to estimate the number of displacements using
# K-P model, not applied anywhere.
def AKPel(A,Z,En,Ed):
	Ec = A*1000
	if (En < Ed):
		Nd = 0.0
	if (En >= Ed and En <= (2*Ed)):
		Nd = 1.0
	if (En > (2.0*Ed) and En <= Ec):
		Nd = En/(2.0*Ed)
	if (En > Ec):
		Nd = Ec/(2.0*Ed)
	AKPel = Nd
	return(AKPel)
		
#=======INELASTIC=======*
	
	# Calculation of neutron dpa and heating cross sections due to the inelastic
	# scattering interactions.
	
def INELASTIC_SCATTERING(ofile_outRMINDD,NPt,Etu,mdisp,Ed,bad,cad):

	sdpa = [0]*NPt; sdpatemp = [0]*NPt; snhttemp = [0]*NPt
	snht = [0]*NPt; alc = [0]*65; alfull = numpy.zeros((NPt,65))
	num_of_displ = numpy.zeros(NPt); tot_energy_products = numpy.zeros(NPt)
	siginel_temp = [0]*NPt; tot_siginel4 = numpy.zeros(NPt)
	# for file4
	En4 = [0]*2000; al4 = numpy.zeros((2000,65))
	# for file5
	NP5 = [0]*3; NE5 = [0]*3; LF = [0]*3
	NF5 = numpy.zeros((3,100))
	Eint = numpy.zeros((3,50)); En5 = numpy.zeros((3,100)); tht = numpy.zeros((3,100))
	Enp5 = numpy.zeros((3,100,500)); f5 = numpy.zeros((3,100,500))
	# for file6 
	fiso = [0]*3; LAW = [0]*3; NP6 = [0]*3; LG = [0]*3; NE6 = [0]*3
	NEP = numpy.zeros((3,100)); NBT = [0]*50; INTr = [0]*50
	Eint6 = numpy.zeros((3,50)); yi = numpy.zeros((3,50)); En = numpy.zeros((3,100))
	Enp = numpy.zeros((3,100,1000)); f = numpy.zeros((3,100,1000)); al6 = numpy.zeros((2000,65))

#--------------------------------------------------------------------
	print( "n, n' .....")

	print('', file = ofile_outRMINDD)
	print(' Inelastic scattering dpa/heating cross section', file = ofile_outRMINDD)
	print('------------------------------------------------', file = ofile_outRMINDD)

	ifile = open("tape01", 'r')
	ifile.readline()
	line = ifile.readline()
	(ZA,AWR,L0,L1,L2,L3,MAT,MF,MT) = line_type1_info(line)
	ifile.close()

	Z = int(ZA/1000)
	A = AWR
	y = A/((A+1)**2)
	n1 = 1/(A+1)
	n2 = A/(A+1)
#--------------------------------------------------------------------
	mta = 50
	for l in range(41):
		iflMTpr = 0
		mta = mta + 1
		MTfind = mta
		iflMTpr = FindMT(MTfind)
		if (iflMTpr == 1):
			iflpr = 0
			ifile = open("tape02", 'r')
			while True:
				line = ifile.readline()
				if (line == ''):
					break
				data = eachlineinfo(line)
				MAT = int(data[6]); MF = int(data[7]); MT = int(data[8])
				if (MAT != -1):
					if (MF == 3):
						if (MT == mta):
							iflpr = 1
							line = ifile.readline()
							data = eachlineinfo(line)
							QM = float(data[0]); QI =  float(data[1])
							LR = int(data[3]); NR = int(data[4]); NP = int(data[5])
							E = [0]*NP; sig = [0]*NP; sdpal = [0]*NP; snhtl = [0]*NP
							ifile.readline()
							(E, sig) = line_type3_info(ifile,NP,2)
				else:
					break
			ifile.close()

			print('', file = ofile_outRMINDD)
			print('Level MT=',mta,' has ',NP,' cross section points', file = ofile_outRMINDD)

			if (iflpr == 1):	# do the following only if MT is present in MF=3
				if4d = 0
				if4c = 0
				ifile = open ("tape01", 'r')
				while True:
					line = ifile.readline()
					if (line == ''):
						break
					data = eachlineinfo(line)
					MAT = int(data[6]); MF = int(data[7]); MT = int(data[8])
					if (MT == 0):
						line = ifile.readline()
						(ZA,AWR,L1,LTTv,L2,L3,MAT,MF,MT) = line_type1_info(line)
					if (MAT != -1):
						if (MF == 4):
							if (MT == mta):
								if (MT < 91):
									if4d = 1
								if (MT == 91):
									if4c = 1
								LTT = LTTv
								line = ifile.readline()
								(ZA,AWR,LI,LCT,L2,L3,MAT,MF,MT) = line_type1_info(line)
								if (LTT == 1  and  LI == 0):
									line = ifile.readline()
									(C1,C2,L1,L2,NR,NE4,MAT,MF,MT) = line_type2_info(line)
									(NBT, INTr) = line_type3_info(ifile,NR,2)
									for i in range(NE4):
										line = ifile.readline()
										data = eachlineinfo(line)
										c1 = float(data[0]); En4[i] = float(data[1]); NL4 = int(data[4])
										NL4 = NL4 + 1
										al4[i][0] = 1
										al4[i][1:-1] = 0.0
										temporary = [0]*NL4
										temporary = line_type3_info (ifile,NL4-1,1)
										for j, value in enumerate(temporary, 1):
											al4[i][j] = value
					else:
						break
				ifile.close()

				ift5c = 0
				if5c = 0
				if (mta == 91 and if4c == 0):
					ifile = open ("tape01", 'r')
					while True:
						line = ifile.readline()
						if (line == ''):
							break
						data = eachlineinfo(line)
						MAT = int(data[6]); MF = int(data[7]); MT = int(data[8])
						if (MT == 0):
							line = ifile.readline()
							(ZA,AWR,L1,L2,NKv,L3,MAT,MF,MT) = line_type1_info(line)
						if (MAT != -1):
							if (MF == 5):
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
											#read(700,12) (NBT(N), INTr(N), N=1, NR)
											#read(700,13) (Eint(l,NSS,N), p(l,NSS,N), N = 1, NP5(l,NSS))
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
						else:
							break
					ifile.close()

				iflawiso = 0 		# for isotropic emitted particle distribution
				if6d = 0
				if6c = 0

				if (if4d == 0 or if4c == 0 or if5c == 0):
					ifile = open ("tape01", 'r')
					while True:
						line = ifile.readline()
						if (line == ''):
							break
						data = eachlineinfo(line)
						MAT = int(data[6]); MF = int(data[7]); MT = int(data[8])
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
								if (MT == mta):
									if (MT < 91):
										if6d = 1
									if (MT == 91):
										if6c = 1
									NK = NKv
									for NSS in range(NK):
										line = ifile.readline()
										data = eachlineinfo(line)
										(ZAP,AWP,LIP,LAW[NSS],NR,NP6[NSS],MAT,MF,MT) = line_type1_info(line)
										(NBT, INTr) = line_type3_info(ifile,NR,2)
										temporary1 = [0]*NP6[NSS]
										temporary2 = [0]*NP6[NSS]
										(temporary1,temporary2) = line_type3_info(ifile,NP6[NSS],2)
										for j, (value1, value2) in enumerate(zip(temporary1, temporary2), 0):
											Eint6[NSS][j] = value1
											yi[NSS][j] = value2
										if (LAW[NSS] == 3):
											iflawiso = 1
										if (LAW[NSS] == 2):
											line = ifile.readline()
											(c1,c2,l3,l4,NR,NE6[NSS],MAT,MF,MT) = line_type2_info(line)
											(NBT, INTr) = line_type3_info(ifile,NR,2)
											for i in range (NE6[NSS]):
												line = ifile.readline()
												(c1,En[NSS][i],LG[NSS],l2,NW,NL6,MAT,MF,MT) = line_type2_info(line)
												NL6 = NL6 + 1
												al6[i][0] = 1
												al6[i][1:-1] = 0.0
												temporary = [0]*NL6
												temporary = line_type3_info (ifile,NL6-1,1)
												for j, value in enumerate(temporary, 1):
													al6[i][j] = value
										if (LAW[NSS] == 1):
											line = ifile.readline()
											(c1,c2,LG[NSS],LP,NR,NE6[NSS],MAT,MF,MT) = line_type2_info(line)
											(NBT, INTr) = line_type3_info(ifile,NR,2)
											for i in range (NE6[NSS]):
												line = ifile.readline()
												(c1,En[NSS][i],ND,NA,NW,NEP[NSS][i],MAT,MF,MT) = line_type2_info(line)
												if (NA != 0):
													if (LG[NSS] == 2 or LG[NSS] == 1):
														fiso[NSS] = 1
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
													fiso[NSS] = 1
													temporary1 = [0]*int(NEP[NSS][i])
													temporary2 = [0]*int(NEP[NSS][i])
													(temporary1,temporary2) = line_type3_info(ifile,int(NEP[NSS][i]),2)
													for j, (value1, value2) in enumerate(zip(temporary1, temporary2), 0):
														Enp[NSS][i][j] = value1
														f[NSS][i][j] = value2
						else:
							break
					ifile.close()

				if (mta == 91 and if4c == 0):
					kdim = int(numpy.max(NEP))
					if (if5c == 1):
						kdim = int(numpy.max(NF5))
					NFi = [0]*NP; Ep = numpy.zeros((kdim,NP)); fEEp1 = numpy.zeros((kdim,NP))
					Ep1 = [0]*kdim; fEEp = [0]*kdim

		#allocate (NF(NE6(41,1)),Ep(2000,5000),fEEp1(2000,5000)) !! dimension is increased from 2000 to 5000 for Mo94 ENDF/B-VII.1
		#allocate (NFi(5000))
#-----------------------------------------------------------------------------------*		

	# When angular distributions for discrete levels 
	# are given in File 6

	# TO GET THE AL(E) FOR THE FULL ENERGY RANGE
				if (if6d == 1):
					print ('', file = ofile_outRMINDD)
					print ('Angular distribution data of discrete level', file = ofile_outRMINDD)
					print ('scattering are given in File 6', file = ofile_outRMINDD)

					if (iflawiso == 1):
						for i in range(NP):
							alfull[i][0] = 1

					if (iflawiso != 1):
						for i in range(NP):
							alfull[i][0] = 1
						for i in range(NP):
							for j in range(NE6[1]-1):
								if (E[i] ==  En[1][j]):
									for k in range (1, 65):
										alfull[i][k] = al6[j][k]
									break
								if (E[i] > En[1][j] and E[i] <= En[1][j+1]):
									diff1 = E[i] - En[1][j]
									diff2 = E[i] - E[i-1]
									if (diff1 <= diff2):
										for k in range (1, 65):
											x = E[i]
											x1 = En[1][j]
											x2 = En[1][j+1]
											y1 = al6[j][k]
											y2 = al6[j+1][k]
											alfull[i][k] = y1+((y2-y1)*(x-x1)/(x2-x1))
									if (diff2 < diff1):
										for k in range (1, 65):
											x = E[i]
											x1 = E[i-1]
											x2 = En[1][j+1]
											y1 = alfull[i-1][k]
											y2 = al6[j+1][k]
											alfull[i][k] = y1+((y2-y1)*(x-x1)/(x2-x1))
									for k in range (1, 65):
										if (alfull[i][k] == 0):
											alfull[i][k] = alfull[i-1][k]
									break

	# When angular distributions for discrete levels and continuum
	# are given in File 4

				if (if4d == 1 or if4c == 1):
					print ('', file = ofile_outRMINDD)
					print ('Angular distribution data of discrete level', file = ofile_outRMINDD)
					print ('scattering are given in File 4', file = ofile_outRMINDD)
					if (if4c == 1):
						print ('', file = ofile_outRMINDD)
						print ('Angular distribution data for continuum', file = ofile_outRMINDD)
						print (' scattering are given in File 4', file = ofile_outRMINDD)
					if (LTT == 0):
						for i in range(NP):
							alfull[i][0] = 1

					if (LTT != 0):
						for i in range (NP):
							alfull[i][0] = 1
						for i in range (NP):
							for j in range (NE4-1):
								if (E[i] ==  En4[j]):
									for k in range (1, 65):
										alfull[i][k] = al4[j][k]
									break
								if (E[i] > En4[j] and E[i] <= En4[j+1]):
									diff1 = E[i] - En4[j]
									diff2 = E[i] - E[i-1]
									if (diff1 <= diff2):
										for k in range (1, 65):
											x = E[i]
											x1 = En4[j]
											x2 = En4[j+1]
											y1 = al4[j][k]
											y2 = al4[j+1][k]
											alfull[i][k] = y1+((y2-y1)*(x-x1)/(x2-x1))
									if (diff2 < diff1):
										for k in range (1, 65):
											x = E[i]
											x1 = E[i-1]
											x2 = En4[j+1]	
											y1 = alfull[i-1][k]
											y2 = al4[j+1][k]
											alfull[i][k] = y1+((y2-y1)*(x-x1)/(x2-x1))
									for k in range (1, 65):
										if (alfull[i][k] == 0):
											alfull[i][k] = alfull[i-1][k]
									break

 	# Point DPA/heating Cross Section -- discrete

				NLa = 65
				if (if4d == 1 and LTT == 0):
					NLa = 3
				if (if6d == 1 and iflawiso == 1):
					NLa = 3
				for i in range (NP):
					Q = QI
					U = (A+1)*(-Q)/A
					if (U < E[i]):
						for k in range (NLa):
							alc[k] = alfull[i][k]
						if (sig[i] == 0):
							dn1 = 0
						if (sig[i] == 0):
							dn2 = 0
						if (sig[i] != 0):
							bta = 1.0
							dn1 = abs(discrec1(alc,NLa,E[i],Q,bta,Z,A,Z,A, Ed,mdisp,bad,cad))*0.8/(2*Ed)
							dn2 = abs(discrec3(alc,NLa,E[i],Q,bta,A))
							num_of_displ[i] = dn1
							tot_energy_products[i] = dn2

						sdpal[i] = sig[i]*dn1
						snhtl[i] = sig[i]*dn2

	# When secondary energy distributions for continuum are given
	# in File 6.

				if (mta == 91  and if4c == 0):
					if (if6c == 1):
						if (NK == 3):
							Nrc = 1		# since arrays and lists are indexed from 0 - (n-1)
						if (NK == 2):
							Nrc = 0
						if (NK == 1):
							Nrc = 0
						if (Nrc == 1):
							print ('', file = ofile_outRMINDD)
							print ('Energy distribution of recoil is given in File 6', file = ofile_outRMINDD)
						if (Nrc == 0):
							print ('', file = ofile_outRMINDD)
							print ('Energy distribution of recoil is not given in', file = ofile_outRMINDD) 
							print ('File 6, calculating using neutron energy' , file = ofile_outRMINDD)
							print ('distribution data', file = ofile_outRMINDD)

	# Interpolation of secondary energies and their distributions 
	# corresponding to the incident energies in File 3.

						for i in range (NP):
							for j in range (int(NE6[Nrc])):
								if (E[i] ==  En[Nrc][j]):
									for k in range (int(NEP[Nrc][j])):
										Ep[k][i] = Enp[Nrc][j][k]
										fEEp1[k][i] = f[Nrc][j][k]
									NFi[i] = NEP[Nrc][j]
									break
								if (E[i] > En[Nrc][j] and E[i] < En[Nrc][j+1]):
									diff1 = E[i] - En[Nrc][j]
									diff2 = E[i] - E[i-1]
									if (diff1 <= diff2):
										for k in range (int(NEP[Nrc][j])):
											x = E[i]
											x1 = En[Nrc][j]	
											x2 = En[Nrc][j+1]
											y1 = Enp[Nrc][j][k]
											y2 = Enp[Nrc][j+1][k]
											y11 = f[Nrc][j][k]
											y22 = f[Nrc][j+1][k]
											Ep[k][i] = y1 + ((y2-y1)*(x-x1)/(x2-x1))
											fEEp1[k][i] = y11 + ((y22-y11)*(x-x1)/(x2-x1))
									if (diff2 < diff1):	    			# "j.ne.1" condition added for Mn55 in ENDF/B-VII.1
										for k in range (int(NEP[Nrc][j])):
											x = E[i]
											x1 = E[i-1]	
											x2 = En[Nrc][j+1]
											y1 = Ep[k][i-1]
											y2 = Enp[Nrc][j+1][k]
											y11 = fEEp1[k][i-1]
											y22 = f[Nrc][j+1][k]
											Ep[k][i] = y1 + ((y2-y1)*(x-x1)/(x2-x1))
											fEEp1[k][i] = y11 + ((y22-y11)*(x-x1)/(x2-x1))
									NFi[i] = NEP[Nrc][j]
									break

	# If energy distribution of emitted neutron in continuum are given 
	# in File 5. Only a few options are coded, since such instances are
	# found rarely.

					if (if5c == 1):
						print('', file = ofile_outRMINDD)
						print('Emitted neutron distribution for continuum', file = ofile_outRMINDD)
						print('reaction is given in File 5', file = ofile_outRMINDD)

						if (LF[0] != 9):
							for i in range (NP):
								for j in range (int(NE5[0])):
									if (E[i] ==  En5[0][j]):
										for k in range (int(NF5[0][j])):
											Ep[k][i] = Enp5[0][j][k]
											fEEp1[k][i] = f5[0][j][k]
										NFi[i] = NF5[0][j]
										break
									if (E[i] > En5[0][j] and E[i] < En5[0][j+1]):
										diff1 = E[i] - En5[0][j]
										diff2 = E[i] - E[i-1]
										if (diff1 <= diff2):
											for k in range (int(NF5[0][j])):
												x = E[i]
												x1 = En5[0][j]
												x2 = En5[0][j+1]
												y1 = Enp5[0][j][k]
												y2 = Enp5[0][j+1][k]
												y11 = f5[0][j][k]
												y22 = f5[0][j+1][k]
												Ep[k][i] = y1 + ((y2-y1)*(x-x1)/(x2-x1))
												fEEp1[k][i] = y11 + ((y22-y11)*(x-x1)/(x2-x1))
										if (diff2 < diff1):
											for k in range (int(NF5[0][j])):
												x = E[i]
												x1 = E[i-1]	
												x2 = En5[0][j+1]
												y1 = Ep[k][i-1]
												y2 = Enp5[0][j+1][k]
												y11 = fEEp1[k][i-1]
												y22 = f5[0][j+1][k]
												Ep[k][i] = y1 + ((y2-y1)*(x-x1)/(x2-x1))
												fEEp1[k][i] = y11 + ((y22-y11)*(x-x1)/(x2-x1))
										NFi[i] = NF5[0][j]
										break

						if (LF[0] == 9):
							theta = tht[0][0]

	# Point DPA/heating Cross Section -- continuum

					for i in range (NP):
						Ex = E[i]
						Q = QI
						for k in range (int(NFi[i])):
							Ep1[k] = Ep[k][i]
							fEEp[k] = fEEp1[k][i]
							if (fEEp[k] >= 1.0):
								fEEp[k] = 0
							if (fiso[Nrc] == 1 or ift5c == 1) :
								fEEp[k] = fEEp1[k][i]
								if (fEEp[k] >= 1.0):
									fEEp[k] = 0
						if (sig[i] == 0):
							dn1 = 0
							dn2 = 0

	# If recoil nucleus energy distributions are given, then the
	# average damage energies and recoil energies are computed directly
	# from these data using functions conint and conintheat respectively.
	# If recoil data are not present, then emitted neutron energy 
	# distributions are used and functions dscrs3 and dscrs3heat are used 
	# along with conint and conintheat.

						if (if6c == 1 and sig[i] != 0):
							dn1 = abs(conint(Z,A,Z,A,Ex,Ed,Ep1,fEEp,int(NFi[i]),Nrc, mdisp,bad,cad))	#alc2d,65,xabc,wg,
							dn2 = abs(conintheat(Z,A,Z,A,Ex,Ep1,fEEp,int(NFi[i]),Nrc))
							if (Nrc == 0):
								dnx = dn1
								dn1 = abs(dscrs3(Z,A,Z,A,Ex,Ed,dnx,mdisp,bad,cad))
								dnx = dn2
								dn2 = abs(dscrs3heat(Z,A,Z,A,Ex,dnx))

	# For tabulated emitted neutron energy distribution data, functions
	# conint1 and conint1heat are used to find average damage energy and
	# recoil energy. For evaporation spectrum data of emitted neutrons,
	# functions Tinteg3 and Tinteg3heat are used.

						if (ift5c == 1 and sig[i] != 0):
							dn1 = abs(conint1(Z,A,Z,A,Ex,Ed,Ep1,fEEp,NFi[i],mdisp,bad,cad))
							dn2 = abs(conint1heat(Z,A,Z,A,Ex,Ep1,fEEp,NFi[i]))
			
						if (if5c == 1 and sig[i] != 0):
							if (ift5c != 1):
								th = Q/n2
								Em1max = Ex - abs(th)
								dn1 = abs(Tinteg3(A,theta,Em1max,Z,Ed,Ex,mdisp,bad,cad))
								dn2 = abs(Tinteg3heat(A,theta,Em1max,Z,Ex))

						dn1 = dn1*0.8/(2*Ed)
						num_of_displ[i] = dn1
						tot_energy_products[i] = dn2
						sdpal[i] = sig[i]*dn1
						snhtl[i] = sig[i]*dn2

#--------------------------------------------------------------------
	# adding all discrete and continuum dpa cross sections
	# into the sdpa array from sdpal at each l

				sdpatemp = trptuqce (E,sdpal,Etu)
				for i in range (NPt):
					sdpa[i] = sdpa[i] + sdpatemp[i]
				snhttemp = trptuqce (E,snhtl,Etu)
				for i in range (NPt):
					snht[i] = snht[i] + snhttemp[i]
				siginel_temp = trptuqce (E,sig,Etu)
				for i in range (NPt):
					tot_siginel4[i] = tot_siginel4[i] + siginel_temp[i]
#--------------------------------------------------------------------	

		if (mta == 91):
			break

	printtofile (NPt,Etu,tot_siginel4,num_of_displ,sdpa,4,0,1)
	printtofile (NPt,Etu,tot_siginel4,tot_energy_products,snht,4,0,2)

#=============================================

#=======(n,xn)=======*

	# Calculation of neutron dpa and heating cross sections due to
	# (n, 2n), (n, 3n) and (n, 4n) reactions.

def n_xn (ofile_outRMINDD,MTi,NPt,Etu,mdisp,Ed,bad,cad):

	sdpat = [0]*NPt; snhtt = [0]*NPt; signxnl = [0]*NPt; alc = [0]*65
	num_of_displ = [0]*NPt; tot_energy_products = [0]*NPt

	# for file5
	NP5 = [0]*3; NE5 = [0]*3; LF = [0]*3
	NF5 = numpy.zeros((3,100))
	Eint = numpy.zeros((3,50)); p = numpy.zeros((3,50)); En5 = numpy.zeros((3,100))
	tht = numpy.zeros((3,100)); Enp5 = numpy.zeros((3,100,200)); f5 = numpy.zeros((3,100,200))

	# for file6
	LAW = [0]*3; NP6 = [0]*3; LG = [0]*3; NE6 = [0]*3; fiso = [0]*3
	NEP = numpy.zeros((3,100)); NBT = [0]*50; INTr = [0]*50
	Eint6 = numpy.zeros((3,50)); yi = numpy.zeros((3,50)); En = numpy.zeros((3,100))
	Enp = numpy.zeros((3,100,200)); f = numpy.zeros((3,100,200))

#-----------------------------------------------------------------------------------*
	iflpresent = 0

	print ('', file = ofile_outRMINDD)
	print ('  (n, xn) MT=',MTi,' dpa/heating cross section', file = ofile_outRMINDD)
	print ('------------------------------------------------', file = ofile_outRMINDD)

	ifile = open ("tape01", 'r')
	ifile.readline() 
	line = ifile.readline()
	(ZA,AWR,L0,L1,L2,L3,MAT,MF,MT) = line_type1_info(line)
	ifile.close()

	Z = int(ZA/1000)
	A = AWR

	ifile = open ("tape02", 'r')
	while True:
		line = ifile.readline()
		if (line == ''):
			break
		data = eachlineinfo(line)
		MAT = int(data[6]); MF = int(data[7]); MT = int(data[8])
		if (MAT != -1):
			if (MF == 3):
				if (MT == MTi):
					iflpresent = 1
					line = ifile.readline()
					data = eachlineinfo(line)
					QM = float(data[0]); QI =  float(data[1])
					LR = int(data[3]); NR = int(data[4]); NP = int(data[5])
					ifile.readline()
					E = [0]*NP; sig = [0]*NP; sdpa = [0]*NP; snht = [0]*NP
					(E, sig) = line_type3_info(ifile,NP,2)
		else:
			break
	ifile.close()

	print('', file = ofile_outRMINDD)
	print('Number of energy points given is ', NP, file = ofile_outRMINDD)

	if (iflpresent == 1):
		ift5 = 0
		if5 = 0
		ifile = open ("tape01", 'r')
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
					if (MT == MTi):
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

		if6 = 0
		if (if5 == 0):
			ifile = open ("tape01", 'r')
			while True:
				line = ifile.readline()
				if (line == ''):
					break
				data = eachlineinfo(line)
				MAT = int(data[6]); MF = int(data[7]); MT = int(data[8])	
				if (MT == 0):
					if (MF == 6):
						line = ifile.readline()
						(ZAv,AWRv,l1,LCT,NKv,l2,MAT,MF,MT) = line_type1_info(line)
					if (MF < 6):
						ifile.readline()
						line = ifile.readline()
						(ZAv,AWRv,l1,LCT,NKv,l2,MAT,MF,MT) = line_type1_info(line)
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
								(NBT, INTr) = line_type3_info(ifile,NR,2)
								temporary1 = [0]*NP6[NSS]
								temporary2 = [0]*NP6[NSS]
								(temporary1,temporary2) = line_type3_info(ifile,NP6[NSS],2)
								for N, (value1, value2) in enumerate(zip(temporary1, temporary2), 0):
									Eint6[NSS][N] = value1
									yi[NSS][N] = value2
		
								if (LAW[NSS] == 2):
									line = ifile.readline()
									(c1,c2,l3,l4,NR,NE6[NSS],MAT,MF,MT) = line_type2_info(line)
									al6 = numpy.zeros((NE6[NSS],65))
									(NBT, INTr) = line_type3_info(ifile,NR,2)
									for i in range (NE6[NSS]):
										line = ifile.readline()
										(c1,En[NSS][i],LG[NSS],l2,NW,NL6,MAT,MF,MT) = line_type2_info(line)
										NL6 = NL6 + 1
										al6[i][0] = 1
										temporary = [0]*NL6
										temporary = line_type3_info (ifile,NL6,1)
										for j, value in enumerate(temporary, 1):
											al6[i][j] = value

								if (LAW[NSS] == 1):
									line = ifile.readline()
									(c1,c2,LG[NSS],LP,NR,NE6[NSS],MAT,MF,MT) = line_type2_info(line)
									temporary1 = [0]*NR
									temporary2 = [0]*NR
									(temporary1,temporary2) = line_type3_info(ifile,NR,2)
									for N, (value1, value2) in enumerate(zip(temporary1, temporary2), 0):
										NBT[N] = value1
										INTr[N] = value2
										# 		read(601,*)
										# 		read(601,*)
									for i in range (NE6[NSS]):
										line = ifile.readline()
										(c1,En[NSS][i],ND,NA,NW,NEP[NSS][i],MAT,MF,MT) = line_type2_info(line)
										if (NA != 0):
											if (LG[NSS] == 2 or LG[NSS] == 1):
												fiso[NSS] = 1
												Ball = [0]*NW
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
											fiso[NSS] = 1
											temporary1 = [0]*int(NEP[NSS][i])
											temporary2 = [0]*int(NEP[NSS][i])
											(temporary1,temporary2) = line_type3_info(ifile,int(NEP[NSS][i]),2)
											for j, (value1, value2) in enumerate(zip(temporary1, temporary2), 0):
												Enp[NSS][i][j] = value1
												f[NSS][i][j] = value2
				else:
					break
		ifile.close()

		kdim = int(numpy.max(NEP))
		if (if5 == 1):
			kdim = int(numpy.max(NF5))
		NFi = [0]*NP; Ep = numpy.zeros((kdim,NP)); fEEp1 = numpy.zeros((kdim,NP)); Ep1 = [0]*kdim; fEEp = [0]*kdim
#--------------------------------------------------------------------
	# If 3 sub-sections are there in File 6 (NK=3), then taking the 
	# data in the recoil sub-section, i.e. Nrc=2. Otherwise, emitted
	# neutron energy distribution data are used with Nrc=1.
	
		if (NK == 3):
			Nrc = 1		# since arrays and lists are indexed from 0 - (n-1)
		if (NK == 2):
			Nrc = 0
		if (NK == 1):
			Nrc = 0

		if (if6 == 1):
			if (Nrc == 1):
				print ('', file = ofile_outRMINDD)
				print ('Energy distribution of recoil nucleus is given in File 6', file = ofile_outRMINDD)
			if (Nrc == 0):
				print ('', file = ofile_outRMINDD)
				print ('Energy distribution of recoil nucleus is', file = ofile_outRMINDD)
				print ('not given in File 6, calculations performed', file = ofile_outRMINDD)
				print ('using emitted neutron distribution data', file = ofile_outRMINDD)

			for i in range (NP):
				for j in range (int(NE6[Nrc])-1):
					if (E[i] ==  En[Nrc][j]):
						for k in range (int(NEP[Nrc][j])):
							Ep[k][i] = Enp[Nrc][j][k]
							fEEp1[k][i] = f[Nrc][j][k]
						NFi[i] = NEP[Nrc][j]
						break
					if (E[i] > En[Nrc][j] and E[i] <= En[Nrc][j+1]):
						diff1 = E[i] - En[Nrc][j]
						diff2 = E[i] - E[i-1]
						if (diff1 <= diff2):
							for k in range (int(NEP[Nrc][j])):
								x = E[i]
								x1 = En[Nrc][j]
								x2 = En[Nrc][j+1]
								y1 = Enp[Nrc][j][k]
								y2 = Enp[Nrc][j+1][k]
								y11 = f[Nrc][j][k]
								y22 = f[Nrc][j+1][k]
								Ep[k][i] = y1 + ((y2-y1)*(x-x1) / (x2-x1))
								fEEp1[k][i] = y11 + ((y22-y11)*(x-x1) / (x2-x1))
						if (diff2 < diff1):
							for k in range (int(NEP[Nrc][j])):
								x = E[i]
								x1 = E[i-1]
								x2 = En[Nrc][j+1]
								y1 = Ep[k][i-1]
								y2 = Enp[Nrc][j+1][k]
								y11 = fEEp1[k][i-1]
								y22 = f[Nrc][j+1][k]
								Ep[k][i] = y1 + ((y2-y1)*(x-x1) / (x2-x1))
								fEEp1[k][i] = y11 + ((y22-y11)*(x-x1) / (x2-x1))
						NFi[i] = NEP[Nrc][j]
						break

#-----------------------------------------------------------------------
	# Only the tabulated data and evaporation spectrum representation
	# of energy distribution of emitted neutrons are coded (as in case
	# of inelastic scattering C - ENDF/B-VII.1). Other cases
	# have not occured.

		if (if5 == 1):
			print ('', file = ofile_outRMINDD)
			print ('Emitted neutron energy distribution data', file = ofile_outRMINDD)
			print ('are given in File 5', file = ofile_outRMINDD)

			if (LF[1] != 9):
				for i in range (NP):
					for j in range (NE5[0]-1):
						if (E[i] ==  En5[0][j]):
							for k in range (NF5[0][j]):
								Ep[k][i] = Enp5[0][j][k]
								fEEp1[k][i] = f5[0][j][k]
							NFi[i] = NF5[0][j]
							break
						if (E[i] > En5[0][j] and E[i] <= En5[0][j+1]):
							diff1 = E[i] - En5[0][j]
							diff2 = E[i] - E[i-1]
							if (diff1 <= diff2):
								for k in range (NF5[0][j]):
									x = E[i]
									x1 = En5[0][j]
									x2 = En5[0][j+1]
									y1 = Enp5[0][j][k]
									y2 = Enp5[0][j+1][k]
									y11 = f5[0][j][k]
									y22 = f5[0][j+1][k]
									Ep[k][i] = y1 + ((y2-y1)*(x-x1) / (x2-x1))
									fEEp1[k][i] = y11 + ((y22-y11)*(x-x1) / (x2-x1))
							if (diff2 < diff1):
								for k in range (NF5[0][j]):
									x = E[i]
									x1 = E[i-1]	
									x2 = En5[0][j+1]
									y1 = Ep[k][i-1]
									y2 = Enp5[0][j+1][k]
									y11 = fEEp1[k][i-1]
									y22 = f5[0][j+1][k]
									Ep[k][i] = y1 + ((y2-y1)*(x-x1) / (x2-x1))
									fEEp1[k][i] = y11 + ((y22-y11)*(x-x1) / (x2-x1))
							NFi[i] = NF5[0][j]
							break

			if (LF[0] == 9):
				theta = tht[0][0]
#-------------------------------------------------------------------------------

		if (MTi == 16):
			A2 = A-1		# recoil nuclei (n,2n) reaction
		if (MTi == 17):
			A2 = A-2		# recoil nuclei (n,3n) reaction
		if (MTi == 37):
			A2 = A-3		# recoil nuclei (n,4n) reaction

		for i in range (NP):
			Ex = E[i]
			if (NFi[i] > kdim):
				NFi[i] = kdim
			for k in range (int(NFi[i])):
				Ep1[k] = Ep[k][i]
				if (fiso[Nrc] ==1 or (if5 == 1 and ift5 == 0)):
					fEEp[k] = fEEp1[k][i]
					if (fEEp[k] >= 1.0):
						fEEp[k] = 0

			if (if6 == 1):
				dn1 = abs(conint(Z,A2,Z,A,Ex,Ed,Ep1,fEEp,int(NFi[i]), Nrc,mdisp,bad,cad))
				dn2 = abs(conintheat(Z,A2,Z,A,Ex,Ep1,fEEp,int(NFi[i]),Nrc))
				if (Nrc == 0):
					dnx = dn1
					dn1 = abs(dscrs3(Z,A2,Z,A,Ex,Ed,dnx,mdisp,bad,cad))
					dnx = dn2
					dn2 = abs(dscrs3heat(Z,A2,Z,A,Ex,dnx))

			if (if5 == 1 and ift5 == 0):
				dn1 = abs(conint1(Z,A2,Z,A,Ex,Ed,Ep1,fEEp,int(NFi[i]), mdisp,bad,cad))
				dn2 = abs(conint1heat(Z,A2,Z,A,Ex,Ep1,fEEp,int(NFi[i])))

			if (if5 == 1 and sig[i] != 0):
				if (ift5 == 1):
					th = U			# QM*(A+1)/A
					Em1max = Ex - abs(th)
					dn1 = abs(Tinteg3(A,theta,Em1max,Z,Ed,Ex,mdisp,bad,cad))
					dn2 = abs(Tinteg3heat(A,theta,Em1max,Z,Ex))

			if (sig[i] == 0):
				dn1 = 0
				dn2 = 0
			dn1 = dn1 * 0.8/(2*Ed)
			num_of_displ[i] = dn1
			tot_energy_products[i] = dn2 
			sdpa[i] = sig[i] * dn1
			snht[i] = sig[i] * dn2

		sdpat = trptuqce (E,sdpa,Etu)
		signxnl = trptuqce (E,sig,Etu)
		printtofile (NPt,Etu,signxnl,num_of_displ,sdpat,MTi,0,1)
		snhtt = trptuqce (E,snht,Etu)
		printtofile (NPt,Etu,signxnl,tot_energy_products,snhtt,MTi,0,2)

	return(signxnl, num_of_displ, tot_energy_products, sdpat, snhtt, iflpresent)
#============================================

#=======(n, anything)=======*
	# Calculation of neutron dpa and heating cross sections due to
	# all inexplicitly given neutron reactions.

def anytnMF6MT5 (ofile_outRMINDD,NPt,Etu,mdisp,Ed,bad,cad):
	sdpat = [0]*NPt; snhtt = [0]*NPt; snhtt_EB = [0]*NPt
	num_of_displ = [0]*NPt; tot_energy_products = [0]*NPt; tot_energy_n_photons = [0]*NPt
	siget = [0]*NPt

	print("n, anything .....")
	print('', file = ofile_outRMINDD)
	print('n, anything reaction MT=',5,' dpa/heating cross section', file = ofile_outRMINDD)
	print('-------------------------------------------------', file = ofile_outRMINDD)

	iflpresent = 0 		# flag for the presence of cross sections MF=3. 
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
					Z = ZA//1000
					AWR = float(AWR)
					A = AWR
					line = ifile.readline()
					data = eachlineinfo(line)
					QM = float(data[0]); QI =  float(data[1])
					LR = int(data[3]); NR = int(data[4]); NP = int(data[5])
					Eall = [0]*NP; sall = [0]*NP; sdall = [0]*NP; shall = [0]*NP; shall_EB = [0]*NP
					num_of_displ1 = [0]*NP; tot_energy_products1 = [0]*NP; tot_energy_n_photons1 = [0]*NP
					dn4_neutrons = [0]*NP; dn5_photons = [0]*NP
					ifile.readline()
					(Eall, sall) = line_type3_info(ifile,NP,2)
		else:
			break
	ifile.close()

	print('', file = ofile_outRMINDD)
	print('Number of energy points given is ', NP, file = ofile_outRMINDD)

## only if MF3 cross sections are available then do the following

	if (iflpresent == 1):

		ifllab = 0 	# secondary energy and angle in lab. system
		iflcom = 0 # secondary energy and angle in c.o.m. system
		iflcomlab = 0 # secondary energy and angle in c.o.m (A<=4), lab. (A>4)

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
						ND = numpy.zeros((NK,400))
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
						if (LCT == 1):
							ifllab = 1
						if (LCT == 2):
							iflcom = 1
						if (LCT == 3):
							iflcomlab = 1
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
									(c1,En[NSS][i],ND[NSS][i],NA,NW,NEP[NSS][i],MAT,MF,MT) = line_type2_info(line)
									if (NA != 0):
										if (LG[NSS] == 1 or LG[NSS] == 2):
											Ball = [0]*NW
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

		ifllargeNP = 0
		if (NP > NPt):
			for i in range (NPt):
				Eall[i] = Etu[i]
				for j in range (NP):
					if (Eall[i] == Eall1[j]):
						sall[i] = sall1[j]
						break
					if (Eall1[j] < Eall[i] and Eall[i] < Eall1[j+1]):
						sall[i] = crstd(Eall[i],Eall1[j],Eall1[j+1],sall1[j],sall1[j+1])
						break
			ifllargeNP = 1
			NP = NPt

		# for declaring dimensions below
		kdim = int(numpy.max(NEP))

		for NSS in range (NK):
			tnYldg = numpy.zeros(NP)
			Zvalue = int(ZAP[NSS]/1000)
			Avalue = AWP[NSS]
			tnYldg = TERPOL(NBT[NSS][:],INTr[NSS][:],Nyld[NSS], Eint6[NSS][:],Yi[NSS][:],NP,Eall)
			ntm = [0]*NP; ntmd = [0]*NP
			if (AWP[NSS] != 0 and ZAP[NSS] != 0 and AWP[NSS] != 1.0 and ZAP[NSS] != 1.0):
				print('', file = ofile_outRMINDD)
				print('Recoil nucleus and charged particle energy', file = ofile_outRMINDD)
				print('distribution data are given in MF 6 MT 5', file = ofile_outRMINDD)

				# this part is for heavy recoil nuclei
				if (Avalue >= 4 and Zvalue > 2):
					print(Zvalue, UtilsU.elementFromZValue(Zvalue), int(ZAP[NSS])%1000,'AWP = ', AWP[NSS])
					for i in range (NP):
						ter6 = [0]*kdim; tf6 = [0]*kdim
						for j in range (NE6[NSS]):
							if (En[NSS][j] == Eall[i]):
								for k in range (int(NEP[NSS][j])):
									ter6[k] = Enp[NSS][j][k]
									tf6[k] = f[NSS][j][k]
								ntm[i] = NEP[NSS][j]
								ntmd[i] = ND[NSS][j]
								break
							if (En[NSS][j] < Eall[i] and Eall[i] < En[NSS][j+1]):
								for k in range (int(NEP[NSS][j])):		#TENDL-2017 Ni59 --> NEP(NSS,j+1)
									x = Eall[i]
									x1 = En[NSS][j]	
									x2 = En[NSS][j+1]
									y1 = Enp[NSS][j][k]
									y2 = Enp[NSS][j+1][k]
									y11 = f[NSS][j][k]
									y22 = f[NSS][j+1][k]
									ter6[k] = TERPOLIN(2,x,x1,x2,y1,y2)
									tf6[k] = TERPOLIN(2,x,x1,x2,y11,y22)
								ntm[i] = NEP[NSS][j] 		 	#TENDL-2017 Ni59 --> NEP(irs,j+1)
								ntmd[i] = ND[NSS][j]
								break

						for k in range (int(ntm[i])):
							if (ter6[k] == 1.0):
								ter6[k] = 0
							if (tf6[k] >= 1.0):
								tf6[k] = 0

						dn1 = 0
						dn2 = 0
						if (sall[i] != 0): 	#  .and. Eall(i)>=abs(QI)
							dn1 = Tinteg1 (Zvalue,Avalue,Z,A,Ed,ter6,tf6,int(ntm[i]),mdisp,bad,cad)
							dn2 = Tinteg1heat (Avalue,ter6,tf6,int(ntm[i]))

						dn1 = abs(dn1)
						dn2 = abs(dn2)

						if (sall[i] == 0):
							sdthis = 0
							shthis = 0
						if (sall[i] != 0):
							dn1 = dn1*0.8/(2*Ed)
							sdthis = tnYldg[i]*sall[i]*dn1
							shthis = tnYldg[i]*sall[i]*dn2
						sdall[i] = sdall[i] + sdthis
						shall[i] = shall[i] + shthis
						num_of_displ1[i] = num_of_displ1[i] + dn1
						tot_energy_products1[i] = tot_energy_products1[i] + dn2

		# THESE PART (BELOW) IS FOR THE ENERGY DEPOSITED BY THE 
		# LIGHT CHARGED PARTICLE DATA PRODUCED IN NEUTRON REACTION, IT
		# IS DONE ONLY IF BOTH CONTINUUM RECOIL AND PARTICLE SECTIONS
		# ARE GIVEN. IT CORRESPONDS TO ONLY HEATING CALCULATIONS, NOT
		# IN THE DISPLACEMENT CROSS SECTIONS.
		# THE ENERGY DEPOSITED BY THE CHARGED PARTICLES WILL BE ADDED 
		# TO THAT DEPOSITED BY THE RECOIL NUCLEUS.

				# this part is only for p, d, t, 3He, 4He (Z<=2 and A<=4, but A!=1)
				if ( Zvalue <= 2 and Avalue <= 4 ):	
					print("Light CP (Z and A):: ", Zvalue, Avalue)
					for i in range (NP):
						ter6 = [0]*kdim; tf6 = [0]*kdim
						for j in range (int(NE6[NSS])):
							if (En[NSS][j] == Eall[i]):
								for k in range (int(NEP[NSS][j])):
									ter6[k] = Enp[NSS][j][k]
									tf6[k] = f[NSS][j][k]
								ntm[i] = NEP[NSS][j]
								ntmd[i] = ND[NSS][j]
								break

							if (En[NSS][j] < Eall[i] and Eall[i] < En[NSS][j+1]):
								for k in range (int(NEP[NSS][j])):
									x = Eall[i]
									x1 = En[NSS][j]
									x2 = En[NSS][j+1]
									y1 = Enp[NSS][j][k]
									y2 = Enp[NSS][j+1][k]
									y11 = f[NSS][j][k]
									y22 = f[NSS][j+1][k]
									ter6[k] = TERPOLIN(2,x,x1,x2,y1,y2)
									tf6[k] = TERPOLIN(2,x,x1,x2,y11,y22)
								ntm[i] = NEP[NSS][j]
								ntmd[i] = ND[NSS][j]
								break

						for k in range (int(ntm[i])):
							#if (iflcomlab == 1):
							#	ter6[k] = ter6[k] + (AWP[NSS]/(AWP[NSS]+1)**2)*Eall[i]
							#if (iflcom == 1):
							#	ter6[k] = ter6[k] + (AWP[NSS]/(AWP[NSS]+1)**2)*Eall[i]
							if (ter6[k] == 1.0):
								ter6[k] = 0
							if (tf6[k] >= 1.0):
								tf6[k] = 0
						dn3 = 0
						if (sall[i] != 0): 		#  .and. Eall(i)>=abs(QI)
							dn3 = tnYldg[i]*Tinteg1heat(AWP[NSS],ter6,tf6,int(ntm[i]))

						dn3 = abs(dn3)
						if (sall[i] == 0):
							shthis = 0
						if (sall[i] != 0):
							shthis = sall[i]*dn3
						shall[i] = shall[i] + shthis
						tot_energy_products1[i] = tot_energy_products1[i] + dn3

		# THE ENERGY CARRIED AWAY BY NEUTRONS AND PHOTONS. 
		# IT CORRESPONDS TO ONLY HEATING CALCULATIONS BY ENERGY BALANCE METHOD.

			if ( (AWP[NSS] == 0 and ZAP[NSS] == 0) or (AWP[NSS] == 1.0 and ZAP[NSS] == 1.0) ):
				print('', file = ofile_outRMINDD)
				print('Secondary neutrons and photons energy', file = ofile_outRMINDD)
				print('distribution data are given in MF 6 MT 5', file = ofile_outRMINDD)

				# for neutrons
				if ( AWP[NSS] == 1.0 and ZAP[NSS] == 1.0 ):
					print("Neutrons (Z and A):: ", int(ZAP[NSS]), AWP[NSS])
				# for photons ----
				if ( AWP[NSS] == 0.0 and ZAP[NSS] == 0.0 ):
					print("Photons (Z and A):: ", int(ZAP[NSS]), AWP[NSS])
				for i in range (NP):
					ter6 = [0]*kdim; tf6 = [0]*kdim
					for j in range (int(NE6[NSS])):
						if (En[NSS][j] == Eall[i]):
							for k in range (int(NEP[NSS][j])):
								ter6[k] = Enp[NSS][j][k]
								tf6[k] = f[NSS][j][k]
							ntm[i] = NEP[NSS][j]
							ntmd[i] = ND[NSS][j]
							break

						if (En[NSS][j] < Eall[i] and Eall[i] < En[NSS][j+1]):
							for k in range (int(NEP[NSS][j])):
								x = Eall[i]
								x1 = En[NSS][j]
								x2 = En[NSS][j+1]
								y1 = Enp[NSS][j][k]
								y2 = Enp[NSS][j+1][k]
								y11 = f[NSS][j][k]
								y22 = f[NSS][j+1][k]
								ter6[k] = TERPOLIN(2,x,x1,x2,y1,y2)
								tf6[k] = TERPOLIN(2,x,x1,x2,y11,y22)
							ntm[i] = NEP[NSS][j]
							ntmd[i] = ND[NSS][j]
							break
					for k in range (int(ntm[i])):
					#	if (iflcomlab == 1):
					#		ter6[k] = ter6[k] + (AWP[NSS]/(AWP[NSS]+1)**2)*Eall[i]
					#	if (iflcom == 1):
					#		ter6[k] = ter6[k] + (AWP[NSS]/(AWP[NSS]+1)**2)*Eall[i]
						if (ter6[k] == 1.0):
							ter6[k] = 0
						if (tf6[k] >= 1.0):
							tf6[k] = 0

					if (sall[i] != 0 ):
						if ( AWP[NSS] == 1.0 and ZAP[NSS] == 1.0 ):
							dn4_neutrons[i] = abs(tnYldg[i]*Tinteg1heat(AWP[NSS],ter6,tf6,int(ntm[i])))

						if ( AWP[NSS] == 0.0 and ZAP[NSS] == 0.0 ):
							dn5_photons[i] = abs(tnYldg[i]*Tinteg1heat(AWP[NSS],ter6,tf6,int(ntm[i])))

						tot_energy_n_photons1[i] = dn4_neutrons[i] + dn5_photons[i]
						shall_EB[i] = sall[i]*(Eall[i] - tot_energy_n_photons1[i])

		if (ifllargeNP == 0):
			Etu = numpy.asarray(Etu)
			Eall = numpy.asarray(Eall)
			sdall = numpy.asarray(sdall)
			shall = numpy.asarray(shall)
			shall_EB = numpy.asarray(shall_EB)
			sall = numpy.asarray(sall)
			num_of_displ1 = numpy.asarray(num_of_displ1)
			tot_energy_products1 = numpy.asarray(tot_energy_products1)
			tot_energy_n_photons1 = numpy.asarray(tot_energy_n_photons1)
			sdpat = trptuqce(Eall,sdall,Etu)
			snhtt = trptuqce(Eall,shall,Etu)
			snhtt_EB = trptuqce(Eall,shall_EB,Etu)
			siget = trptuqce(Eall,sall,Etu)
			num_of_displ = trptuqce(Eall,num_of_displ1,Etu)
			tot_energy_products = trptuqce(Eall,tot_energy_products1,Etu)
			tot_energy_n_photons = trptuqce(Eall,tot_energy_n_photons1,Etu)

		if (ifllargeNP == 1):
			sdpat = sdall
			snhtt = shall
			snhtt_EB = shall_EB
			num_of_displ = num_of_displ1
			tot_energy_products = tot_energy_products1
			tot_energy_n_photons = tot_energy_n_photons1

		# THIS PART IS ADDED TO REMOVE SUDDEN NON-ZERO VALUES BELOW THE 
		# STARTING ENERGY POINT (WHEN DOING W180 FROM ENDF/B-VII.1)
		#-------------
	
		#for i in range (NPt):
		#	if (Etu[i] >= 2.0e+7):
		#		nbstart = i - 1
		#		break

		for i in reversed(range (len(Etu))):
			if (sdpat[i] == 0 and (sdpat[i-1] > 0 and sdpat[i+1] > 0)):
				for j in range(i):
					sdpat[j] = 0
				break

		for i in reversed(range (len(Etu))):
			if (snhtt[i] == 0 and (snhtt[i-1] > 0 and snhtt[i+1] > 0)):
				for j in range(i):
					snhtt[j] = 0
				break
		# --------------

		#for i = nbstart, 1, -1
		#		if (sdpat(i)==0) then
		#			do j = 1, i
		#				sdpat(j)=0
		#			end do
		#			exit
		#		end if
		#	end do
		#	do i = nbstart, 1, -1
		#		if (snhtt(i)==0) then
		#			do j = 1, i
		#				snhtt(j)=0
		#			end do
		#			exit
		#		end if
		#	end do

		printtofile (NPt,Etu,siget,num_of_displ,sdpat,5001,0,1)
		printtofile (NPt,Etu,siget,tot_energy_products,snhtt,5001,0,2)
		printtofile (NPt,Etu,siget,tot_energy_n_photons,snhtt_EB,5001,0,3)
	
		#--------------
	# **** the above is done only if MF3 for that MT is present ****

	# anytnMF6MT5 (n, anything) completes here
#====================================================

#=======Make the required unique common energy and call reactions=======*
	
	# It creates an unique energy array out of the MF3 MT1 energy points
	# given in the pre-processed ENDF-6 file which is used as the common 
	# energy array for partial and total dpa and heating cross section 
	# computation. This energy array is called unique because any repetition
	# (may occur at the dense resonances regions) in the pre-processed energy 
	# points are found out and removed, so that each energy point is present
	# only once. It calls other reaction-specific subroutines and their required
	# multigrouping according to the inputs given.
	
def uqce (ofile_outRMINDD,insp,nra,nreac,mdisp,Ed,bad,cad,mgyn,igtype):

	# Et=Energy array in MT=1, Etu=unique of Et
	# extraction of total energy points

	ifile = open ('tape02', 'r') # tape02=output of (reconr + broadr) PENDF

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
					QM = float(data[0]); QI =  float(data[1]); NR = int(data[4]); NPt = int(data[5])
					Et = [0]*NPt
					LR = int(ifile.readline().split()[1])
					i = 0
					while (i < NPt):
						line = ifile.readline()
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
	ifile.close()

	print('', file = ofile_outRMINDD)
	print(NPt,' Total cross sections energy points', file = ofile_outRMINDD)
	
	# make unique common energy
	
	Etu = numpy.array(Et)
	Etu = numpy.unique(Etu)
	NPt = len(Etu)
	
	print('', file = ofile_outRMINDD)
	print(NPt,' Unique total cross sections energy points', file = ofile_outRMINDD)

	# call reactions
	
	# binary variables required in case  nra(i)=6
		
	irct1y=0; irct2y=0; irct3y=0; irct4y=0; irct5y=0; irct6y=0

	for i in range(nreac):

		irct = nra[i]
		if (irct == 1):
			RADIATIVE_CAPTURE (ofile_outRMINDD,NPt,Etu,mdisp,Ed,bad,cad) # MT=102
			if (mgyn==1):
				groupmulti(insp,102,igtype,1) 	# 102=id for output file, 1=DPA
				groupmulti(insp,102,igtype,2) 	# 102=id for output file, 2=Heating
			irct1y=1
		if (irct==2):
			ELASTIC_SCATTERING (ofile_outRMINDD,NPt,Etu,mdisp,Ed,bad,cad) 	# MT=2
			if (mgyn==1): 
				groupmulti(insp,2,igtype,1)
				groupmulti(insp,2,igtype,2)
			irct2y=1
		if (irct==3):
			INELASTIC_SCATTERING (ofile_outRMINDD,NPt,Etu,mdisp,Ed,bad,cad) 		# MT=51 to 91
			if (mgyn==1):
				groupmulti(insp,4,igtype,1)
				groupmulti(insp,4,igtype,2)
			irct3y=1
		if (irct==4):							 # (n,2n), (n,3n) and (n,4n)
			CONTROL_nxn (ofile_outRMINDD,NPt,Etu,mdisp,Ed,bad,cad) 
			if (mgyn==1):
				groupmulti(insp,1601,igtype,1)
				groupmulti(insp,1601,igtype,2)
			irct4y=1
		if (irct==5):
			CONTROL_nCPO(ofile_outRMINDD,NPt,Etu,mdisp,Ed,bad,cad) 		# (n,CPO)
			if (mgyn==1):
				groupmulti(insp,3001,igtype,1)
				groupmulti(insp,3001,igtype,2)
			irct5y=1
		if (irct==6):
			anytnMF6MT5(ofile_outRMINDD,NPt,Etu,mdisp,Ed,bad,cad)		# (n, anything)
			if (mgyn==1):
				groupmulti(insp,5001,igtype,1)
				groupmulti(insp,5001,igtype,2)
				groupmulti(insp,5001,igtype,3)
			irct6y=1
		if (irct==7): 				# Total DPA and Heating due to incident neutron
			if(irct1y==0):
				RADIATIVE_CAPTURE (ofile_outRMINDD,NPt,Etu,mdisp,Ed,bad,cad)
			if(irct2y==0):
				ELASTIC_SCATTERING (ofile_outRMINDD,NPt,Etu,mdisp,Ed,bad,cad)
			if(irct3y==0):
				INELASTIC_SCATTERING (ofile_outRMINDD,NPt,Etu,mdisp,Ed,bad,cad)
			if(irct4y==0):
				CONTROL_nxn(ofile_outRMINDD,NPt,Etu,mdisp,Ed,bad,cad)
			if(irct5y==0):
				CONTROL_nCPO(ofile_outRMINDD,NPt,Etu,mdisp,Ed,bad,cad)
			if(irct6y==0):
				anytnMF6MT5(ofile_outRMINDD,NPt,Etu,mdisp,Ed,bad,cad)
			total(NPt,Etu)				# Add contributions to DPA,Heating from all partial reactions
			if (mgyn==1):
				groupmulti(insp,1,igtype,1) 		# 1= id for output file, 1 = DPA
				groupmulti(insp,1,igtype,2) 		# 1= id for output file, 2 = Heating
				groupmulti(insp,1,igtype,3) 		# 1= id for output file, 3 = Heating_Energy_balance

# ==============================================================

# =======Total=======*
# Calculation of total neutron dpa and heating cross sections
# by adding together the partial contributions from individual
# reactions.

def total (NPt, Etu):
	print("n, all .....")

	for iprd in range(1, 4):
		tot_sig = [0]*NPt
		tot_dhval = [0]*NPt
		tot_dhsig = [0]*NPt
		for nrct in range (1, 7):
			s1 = [0]*NPt
			if (nrct == 1 and iprd == 1):
				ifile = open("ndpa2.txt", 'r')
			if (nrct == 2 and iprd == 1):
				ifile = open("ndpa4.txt", 'r')
			if (nrct == 3 and iprd == 1):
				ifile = open("ndpa1601.txt", 'r')
			if (nrct == 4 and iprd == 1):
				ifile = open("ndpa102.txt", 'r') 
			if (nrct == 5 and iprd == 1):
				ifile = open("ndpa3001.txt", 'r')
			if (nrct == 6 and iprd == 1):
				ifile = open("ndpa5001.txt", 'r') 

			if (nrct == 1 and iprd == 2):
				ifile = open("nheat2.txt", 'r')
			if (nrct == 2 and iprd == 2):
				ifile = open("nheat4.txt", 'r')
			if (nrct == 3 and iprd == 2):
				ifile = open("nheat1601.txt", 'r')
			if (nrct == 4 and iprd == 2):
				ifile = open("nheat102.txt", 'r') 
			if (nrct == 5 and iprd == 2):
				ifile = open("nheat3001.txt", 'r') 
			if (nrct == 6 and iprd == 2):
				ifile = open("nheat5001.txt", 'r')

			if (nrct == 1 and iprd == 3):
				ifile = open("nheat2.txt", 'r')
			if (nrct == 2 and iprd == 3):
				ifile = open("nheat4.txt", 'r')
			if (nrct == 3 and iprd == 3):
				ifile = open("nheat1601.txt", 'r')
			if (nrct == 4 and iprd == 3):
				ifile = open("nheat102.txt", 'r') 
			if (nrct == 5 and iprd == 3):
				ifile = open("nheat_EB3001.txt", 'r') 
			if (nrct == 6 and iprd == 3):
				ifile = open("nheat_EB5001.txt", 'r')

			npmt = int(ifile.readline().split()[0])
			if (npmt != 0):
				for i in range(NPt):
					data = ifile.readline().split()
					s1[i] = float(data[1])
					tot_sig[i] = tot_sig[i] + s1[i]
					s1[i] = float(data[2])
					tot_dhval[i] = tot_dhval[i] + s1[i]
					s1[i] = float(data[3])
					tot_dhsig[i] = tot_dhsig[i] + s1[i]
			ifile.close()
		printtofile (NPt,Etu,tot_sig,tot_dhval,tot_dhsig,1,0,iprd)


# ===================================================================

# =======Print the dpa / heating cross sections to files=======*
		
	# The point and multigrouped dpa / heating cross sections are
	# written in files in the energy -- cross sections two-column format.
	# The files are named according to point / grouped, dpa / heat and
	# the MT number of the partial contribution. The total point dpa 
	# and heating cross sections are written only up to 20 MeV since in 
	# general information above 20 MeV is available only for (n,n) partial.
		
def printtofile (NPt,Etu,siget,disp_heat_value,sdpat,MTtp,iflgrouped,ifldh):
	MTname = str(MTtp)

	# NPt = Total points in Etu, sdpat=cross section, 
	# MTtp=id for output file, iflgrouped=0, 
	# point data: 1, grouped data, ifldh, 1 = DPA: 2 = Heating 
	# Create o/p files with name as: [ndpa(grouped)+id for output file],
	# [nheat(grouped)+ id for outputfile].
		
	if (ifldh == 1):
		if (iflgrouped == 0):
			dpafilename_point = 'ndpa' + MTname + '.txt'
			ofile100 = open(dpafilename_point, 'w')
		if (iflgrouped == 1):
			dpafilename_grouped = 'ndpagrouped' + MTname + '.txt'
			ofile100 = open(dpafilename_grouped, 'w')

	if (ifldh==2):
		if (iflgrouped == 0):
			heatfilename_point = 'nheat' + MTname + '.txt'
			ofile100 = open(heatfilename_point, 'w')
		if (iflgrouped == 1):
			heatfilename_grouped = 'nheatgrouped' + MTname + '.txt'
			ofile100 = open(heatfilename_grouped, 'w')

	if (ifldh==3):
		if (iflgrouped == 0):
			heatfilename_point = 'nheat_EB' + MTname + '.txt'
			ofile100 = open(heatfilename_point, 'w')
		if (iflgrouped == 1):
			heatfilename_grouped = 'nheatgrouped_EB' + MTname + '.txt'
			ofile100 = open(heatfilename_grouped, 'w')

	NPtlast = NPt
	
	print(NPtlast, file = ofile100)
	if (NPt != 0):
		for i in range(NPtlast):
			print('%13.6e%13.6e%14.6e%14.6e' %(Etu[i], siget[i], disp_heat_value[i], sdpat[i]), file = ofile100)

	ofile100.close()

#=======Get MFs and MTs from File 1 of raw data=======*
 	
	## The directory of Files (MFs) and the corresponding Sections (MTs)
	## given in the evaluation are read from File 1 and returned.
	
def file1(): 
		# maximum of NXC = 350 (ENDF-102), 
		# but deviates for Mn55 ENDF/B-VII.1, so changed to 1000 
	ifile = open('tape01', 'r')
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
	ifile.close()

	return(nfiles,MFs,MTs)

#=======To find for availble MT=======*
	
	# The presence of a particular reaction cross section (MT in MF3) is
	# searched from the information about evaluation given in File 1 of raw
	# ENDF-6 file.
	
def FindMT(MTfind):		# MTfind = raw ENDF MT, iflMTpr = output from FindMT
	MFs = numpy.zeros(1000); MTs = numpy.zeros(1000)
	(nfiles,MFs,MTs) = file1()
		
	iflMTpr = 0
	for i in range(nfiles):
		if (MFs[i]==3 and MTs[i]==MTfind):
			iflMTpr = 1
			break
	return(iflMTpr)

#=======Multigroup=======*

	## Multigroup, according to requirement, the point dpa and heating
	## cross sections into the chosen neutron energy group structure.
	
def groupmulti(insp,MTg,igtype,ifldh):
	MTgname =  str(MTg)
		
	print("Group .....")
		
	if (insp == 1):
		ifile = open('NeutronSpectrum.txt', 'r')
		nre = int(ifile.readline().split()[-1])
		Ngl = nre + 1
		fi = numpy.zeros(Ngl)
		for i in reversed(range(nre)):
			fi[i] = float(ifile.readline().split()[0])
		fi[-1] = fi[-2]
		ifile.close()

	if (ifldh == 1):
		filename = 'ndpa' + MTgname + '.txt'
		ifile = open(filename, 'r')
	if (ifldh == 2):
		filename = 'nheat' + MTgname + '.txt'
		ifile = open(filename, 'r')
	if (ifldh == 3):
		filename = 'nheat_EB' + MTgname + '.txt'
		ifile = open(filename, 'r')

	NP = int(ifile.readline().split()[-1])
	if (NP != 0):					#do the following only if cross sections are present
		E = numpy.zeros(NP); siget = numpy.zeros(NP); disp_heat_val = numpy.zeros(NP); sdpa = numpy.zeros(NP)
		for i in range(NP):
			line = ifile.readline(); data = line.split()
			E[i] = float(data[0]); siget[i] = float(data[1])
			disp_heat_val[i] = abs(float(data[2])); sdpa[i] = abs(float(data[3]))
		ifile.close()

		if (igtype==0):
			ifile = open('Energy-GroupLimits.txt', 'r')
			Ngl = int(ifile.readline().split()[-1])
			Eg = numpy.zeros(Ngl); gsiget = numpy.zeros(Ngl)
			gdhval = numpy.zeros(Ngl); gsdpa = numpy.zeros(Ngl)
			for i in reversed(range(Ngl)):
				Eg[i] = float(ifile.readline().split()[0])
			ifile.close()

		if (igtype==1):
			Ngl = 176
			nre = Ngl-1
			Eg = numpy.zeros(Ngl); gsiget = numpy.zeros(Ngl)
			gdhval = numpy.zeros(Ngl); gsdpa = numpy.zeros(Ngl)
			Eg = engrp1()

		if (igtype==2):
			Ngl = 27
			nre = Ngl-1
			Eg = numpy.zeros(Ngl); gsiget = numpy.zeros(Ngl)
			gdhval = numpy.zeros(Ngl); gsdpa = numpy.zeros(Ngl)
			Eg = engrp2()

		if (igtype==3):
			Ngl = 34
			nre = Ngl-1
			Eg = numpy.zeros(Ngl); gsiget = numpy.zeros(Ngl)
			gdhval = numpy.zeros(Ngl); gsdpa = numpy.zeros(Ngl)
			Eg = engrp3()

		if (igtype==4):
			Ngl = 239
			nre = Ngl-1
			Eg = numpy.zeros(Ngl); gsiget = numpy.zeros(Ngl)
			gdhval = numpy.zeros(Ngl); gsdpa = numpy.zeros(Ngl)
			Eg = engrp4()

		if (igtype==5): 
			Ngl = 199
			nre = Ngl-1
			Eg = numpy.zeros(Ngl); gsiget = numpy.zeros(Ngl)
			gdhval = numpy.zeros(Ngl); gsdpa = numpy.zeros(Ngl)
			Eg = engrp5()

		if (igtype==6):
			Ngl = 710
			nre = Ngl-1
			Eg = numpy.zeros(Ngl); gsiget = numpy.zeros(Ngl)
			gdhval = numpy.zeros(Ngl); gsdpa = numpy.zeros(Ngl)
			Eg = engrp6()

		if (igtype==7):
			Ngl = 641
			nre = Ngl-1
			Eg = numpy.zeros(Ngl); gsiget = numpy.zeros(Ngl)
			gdhval = numpy.zeros(Ngl); gsdpa = numpy.zeros(Ngl)
			Eg = engrp7()

		if (igtype==8):
			Ngl = 101
			nre = Ngl-1
			Eg = numpy.zeros(Ngl); gsiget = numpy.zeros(Ngl)
			gdhval = numpy.zeros(Ngl); gsdpa = numpy.zeros(Ngl)
			Eg = engrp8()

		if (igtype==9):
			Ngl = 48
			nre = Ngl-1
			Eg = numpy.zeros(Ngl); gsiget = numpy.zeros(Ngl)
			gdhval = numpy.zeros(Ngl); gsdpa = numpy.zeros(Ngl)
			Eg = engrp9()

		if (igtype==10):
			Ngl = 101			# DLC-2 group structure
			nre = Ngl-1
			Eg = numpy.zeros(Ngl); gsiget = numpy.zeros(Ngl)
			gdhval = numpy.zeros(Ngl); gsdpa = numpy.zeros(Ngl)
			Eg = engrp10()

		if (igtype==11):
			Ngl = 229			# 229 group structure
			nre = Ngl-1
			Eg = numpy.zeros(Ngl); gsiget = numpy.zeros(Ngl)
			gdhval = numpy.zeros(Ngl); gsdpa = numpy.zeros(Ngl)
			Eg = engrp11()

		if (igtype==12):
			Ngl = 229			# 229 group structure
			nre = Ngl-1
			Eg = numpy.zeros(Ngl); gsiget = numpy.zeros(Ngl)
			gdhval = numpy.zeros(Ngl); gsdpa = numpy.zeros(Ngl)
			Eg = engrp12()

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

				bcs1 = numpy.interp(t+h, E, siget) * flux1
				dhval1 = numpy.interp(t+h, E, disp_heat_val) * flux1
				dhcs1 = numpy.interp(t+h, E, sdpa) * flux1
				deno1 = flux1

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

				bcs2 = numpy.interp(t+h, E, siget) * flux2
				dhval2 = numpy.interp(t+h, E, disp_heat_val) * flux2
				dhcs2 = numpy.interp(t+h, E, sdpa) * flux2
				deno2 = flux2

				bcs = bcs + (h/2)*(bcs1 + bcs2)
				dhval = dhval + (h/2)*(dhval1 + dhval2)
				dhcs = dhcs + (h/2)*(dhcs1 + dhcs2)
				denominator = denominator + (h/2)*(deno1 + deno2)

			if (bcs != 0 and denominator != 0):
				gsiget[i] = bcs/denominator
			if (dhval != 0 and denominator != 0):
				gdhval[i] = dhval/denominator
			if (dhcs != 0 and denominator != 0):
				gsdpa[i] = dhcs/denominator

		printtofile(Ngl,Eg,gsiget,gdhval,gsdpa,MTg,1,ifldh)

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
