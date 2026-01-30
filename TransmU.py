!	>> Code: RMINDD.py, module file: TransmU.py
!	>> Perform: Activation, Gas production and Transmutation calculations 
!	>> 			from basic ENDF-6 files
!	>> Author: Dr. Uttiyoarnab Saha
!	>> Version and Date: 1.0 and 14/03/2021
!	>> Last modified: 16/03/2021, Kolkata
!	>> Update: 16/03/2021
!	>> Major changes:
!
!======================================================================================================

import numpy
import math
import os

def TransmCal (eliso,igtyp,insp,iconcyn):
	
	!! Calculation of gas production and total activation cross section
	!! due to charged particle production reactions and (n,xn) and (n,g)
	!! reactions of neutron.

tnEyld = numpy.zeros((200,200)); tnYld = numpy.zeros((200,200))
tnNBTp = numpy.zeros((200,20)); tnINTrp = numpy.zeros((200,20))
tnNyld = [0]*200; tnZp = [0]*200

Eyld = numpy.zeros((5,200)); Yld = numpy.zeros((5,200))
Nyld = [0]*5; iflMTtppr = [0]*5
NBTp = numpy.zeros((5,20)); INTrp = numpy.zeros((5,20))
NBTpp = [0]*20; NBTpd = [0]*20; NBTpt = [0]*20; NBTp3He = [0]*20; NBTpal = [0]*20
INTrpp = [0]*20; INTrpd = [0]*20; INTrpt = [0]*20; INTrp3He = [0]*20; INTrpal = [0]*20

!! Assuming maximum 200 transmuted isotopes from MF6 MT5 and 
!! remaining are explicit transmutation reactions
Zvaltrack = [0]*230; Avaltrack = [0]*230; cntrack = [0]*230
	
ofile400 = open("TransmGas-group.txt", "a")
ofile401 = open("TransmNucl-group.txt", "a")
ofile402 = open("TransmGas-point.txt", "w")
ofile403 = open("TransmNucl-MF5-point.txt", "w")
ofile_outRMINDD = open("Output_RadEMC-TransmU.txt", "a")
ofile405 = open("TransmNucl-Net-group.txt", "a")

		
print('E/(n,xp)/(n,xd)/(n,xt)/(n,x3He)/(n,xa)/(n,act.)', file = ofile400)
print('Energy(eV) - Cross section(barns)', file = ofile401)
print('E/(n,xp)/(n,xd)/(n,xt)/(n,x3He)/(n,xa)/(n,act.)', file = ofile402)
print('Energy(eV) - Cross section(barns)', file = ofile403)
print('Energy (eV) / Cross-section (barns)', file = ofile405)

nrab = 37
MTnum = [0]*42; iflag = [0]*42

MTnum = [5,11,16,17,22,23,24,25,28,29,30,32,33,34,35,36,37,41, \
	42,44,45,102,103,104,105,106,107,108,109,110,111,112,113,114, \
	115,116,117,203,204,205,206,207]


	!! Read and Store the Values of Z and A of the target nucleus

ifile = open ('tape02', 'r')
while True:
	read(500,1) MAT, MF, MT, NS
 	if (MT==0 .and. MF==0) then
		read(500,11)ZAv,AWRv,L0,L1,NKv,L2,MAT,MF,MT,NS
		Ztarget = int(ZAv)/1000
		Atarget =  mod(int(ZAv),1000)
		exit
 	end if
end do
close(unit=500)


open (unit = 500, file = "tape02")

do 
	read (500,5)  MAT, MF, MT
 	if (MAT.ne.-1) then
 		if (MF.eq.3) then
			if (MT.eq.1) then
				read(500,6) QM, QI, NR, NP
				allocate (Et(NP), E(NP))
				read(500,7) LR
				read(500,8) (Et(i), s, i = 1,NP)
			end if
 		end if
	else
 		exit
 	end if
end do
		
close(unit=500)
		
write(404,*) NP, ' Total cross section energy points'
		
write(*,*) 'Unique total energy array .....'
k = 1
k1 = 1
do i = 1, NP
 	if (k1 <= NP) then
 		E(k) = Et(k1)
 		k = k+1
 	end if
	cnt = 1
  	do j = k1+1, NP
 		if (Et(k1) == Et(j)) then
 			cnt = cnt+1
 		end if
 	end do
 	k1 = k1+cnt
end do
		
NP = k-1
		
write(404,*) NP, ' Unique energy points'
		
do i = 1, NP
	if (E(i)>=20.0e+06) then
		ifull = i
		exit
	end if
end do
NP = ifull
		
write(404,*) NP, ' Energy points up to 20 MeV'
				
write(*,*) 'Reading Data .....'
	
	!! Gas producion cross sections MT = 203, ...., 207	
	!! may be given sometimes; nreac = 38 to 42 => these MTs
	
NPt = size(E)

allocate(sig5tot(NPt),sigparttot(NPt),sppoint(NPt),sdpoint(NPt), &
		strpoint(NPt),s3Hepoint(NPt),sapoint(NPt),sigtpoint(NPt))
		
do nreac = 38, 42
	open (unit = 500, file = "tape02")

	iflag(nreac) = 0
	iflMTtppr(nreac-37) = 0
	
	MTi = MTnum(nreac)
		
	!! extraction of cross sections
	NP1 = 0
 	do 
		read(500,1) MAT, MF, MT, NS
		if ( MT==0 .and. MF.ne.0) then	! .and. NS==99999
			read(500,11)ZAv,AWRv,L0,L1,NKv,L2,MAT,MF,MT,NS
		end if
		if (MAT.ne.-1) then
			if (MF.eq.3) then
				if (MT.eq.MTi) then
					iflag(nreac) = 1
					iflMTtppr(nreac-37) = 1
					write(*,*) 'MT = ', MTi, '.....'
					write(404,*) 'Gas production MT found: ', MTi
					ZA = ZAv
					AWR = AWRv
					read(500,2) QM, QI, LR, NR, NP1
					allocate(E1(NP1),sig1(NP1))
					read(500,*) 
					read(500,3) (E1(j1), sig1(j1), j1 = 1, NP1)
				end if
			end if
		else
			exit
		end if
 	end do

	if (iflag(nreac)==1) then
		
		if(MTi==203) then
			allocate (E203(NP1),sig203(NP1))
			E203 = E1
			sig203 = 0
			sig203 = sig1
 		end if
		if(MTi==204) then
			allocate (E204(NP1),sig204(NP1))
			E204 = E1
			sig204 = 0
			sig204 = sig1
 		end if
		if(MTi==205) then
			allocate (E205(NP1),sig205(NP1))
			E205 = E1
			sig205 = 0
			sig205 = sig1
 		end if
		if(MTi==206) then
			allocate (E206(NP1),sig206(NP1))
			E206 = E1
			sig206 = 0
			sig206 = sig1
 		end if
		if(MTi==207) then
			allocate (E207(NP1),sig207(NP1))
			E207 = E1
			sig207 = 0
			sig207 = sig1
 		end if
		
		E1 = 0
		sig1 = 0
		deallocate(E1,sig1)
		
	end if
		
	close(unit=500)
end do
		
iflp103 = 0
iflp104 = 0
iflp105 = 0
iflp106 = 0
iflp107 = 0
	
	!! Calculation of gas production and activation from individual
	!! CPO and other reactions; nreac = 1 to nrab(=37) => these MTs.
			
do nreac = 1, nrab
	open (unit = 500, file = "tape02")

	iflag(nreac) = 0
	
	MTi = MTnum(nreac)
	
	!! extraction of cross sections
	NP1 = 0
 	do 
 		read(500,1) MAT, MF, MT, NS
 		if (MT==0 .and. MF.ne.0) then	! .and. NS==99999
			read(500,11)ZAv,AWRv,L0,L1,NKv,L2,MAT,MF,MT,NS
 		end if
 		if (MAT.ne.-1) then
			if (MF.eq.3) then
				if (MT.eq.MTi) then
					iflag(nreac) = 1
					if (MT==103) iflp103 = 1
					if (MT==104) iflp104 = 1
					if (MT==105) iflp105 = 1
					if (MT==106) iflp106 = 1
					if (MT==107) iflp107 = 1
					write(*,*) 'MT = ', MTi, '.....'
					ZA = ZAv
					AWR = AWRv
					read(500,2) QM, QI, LR, NR, NP1
					allocate(E1(NP1),sig1(NP1))
					read(500,*) 
					read(500,3) (E1(j), sig1(j), j = 1, NP1)
				end if
			end if
		else
			exit
 		end if
 	end do

	if (iflag(nreac)==1) then
		
		if(MTi==5) then
			write(404,*)
			write(404,*) '(n, anything) MT = 5 given'
			allocate (E5(NP1),sig5(NP1))
			E5 = E1
			sig5 = 0
			sig5 = sig1
 		end if
		if(MTi==11) then
			write(404,*)
			write(404,*) 'Activation and Gas production MT = 11 given'
			allocate (E11(NP1),sig11(NP1))
			E11 = E1
			sig11 = 0
			sig11 = sig1
 		end if
		if(MTi==16) then
			write(404,*)
			write(404,*) 'Activation MT = 16 given'
			allocate (E16(NP1),sig16(NP1))
			E16 = E1
			sig16 = 0
			sig16 = sig1
 		end if
		if(MTi==17) then
			write(404,*)
			write(404,*) 'Activation MT = 17 given'
			allocate (E17(NP1),sig17(NP1))
			E17 = E1
			sig17 = 0
			sig17 = sig1
 		end if
		if(MTi==22) then
			write(404,*)
			write(404,*) 'Activation and Gas production MT = 22 given'
			allocate (E22(NP1),sig22(NP1))
			E22 = E1
			sig22 = 0
			sig22 = sig1
 		end if
		if(MTi==23) then
			write(404,*)
			write(404,*) 'Activation and Gas production MT = 23 given'
			allocate (E23(NP1),sig23(NP1))
			E23 = E1
			sig23 = 0
			sig23 = sig1
 		end if
		if(MTi==24) then
			write(404,*)
			write(404,*) 'Activation and Gas production MT = 24 given'
			allocate (E24(NP1),sig24(NP1))
			E24 = E1
			sig24 = 0
			sig24 = sig1
 		end if
		if(MTi==25) then
			write(404,*)
			write(404,*) 'Activation and Gas production MT = 25 given'
			allocate (E25(NP1),sig25(NP1))
			E25 = E1
			sig25 = 0
			sig25 = sig1
 		end if
		if(MTi==28) then
			write(404,*)
			write(404,*) 'Activation and Gas production MT = 28 given'
			allocate (E28(NP1),sig28(NP1))
			E28 = E1
			sig28 = 0
			sig28 = sig1
 		end if
		if(MTi==29) then
			write(404,*)
			write(404,*) 'Activation and Gas production MT = 29 given'
			allocate (E29(NP1),sig29(NP1))
			E29 = E1
			sig29 = 0
			sig29 = sig1
 		end if
		if(MTi==30) then
			write(404,*)
			write(404,*) 'Activation and Gas production MT = 30 given'
			allocate (E30(NP1),sig30(NP1))
			E30 = E1
			sig30 = 0
			sig30 = sig1
 		end if
		if(MTi==32) then
			write(404,*)
			write(404,*) 'Activation and Gas production MT = 32 given'
			allocate (E32(NP1),sig32(NP1))
			E32 = E1
			sig32 = 0
			sig32 = sig1
 		end if
		if(MTi==33) then
			write(404,*)
			write(404,*) 'Activation and Gas production MT = 33 given'
			allocate (E33(NP1),sig33(NP1))
			E33 = E1
			sig33 = 0
			sig33 = sig1
 		end if
		if(MTi==34) then
			write(404,*)
			write(404,*) 'Activation and Gas production MT = 34 given'
			allocate (E34(NP1),sig34(NP1))
			E34 = E1
			sig34 = 0
			sig34 = sig1
 		end if
		if(MTi==35) then
			write(404,*)
			write(404,*) 'Activation and Gas production MT = 35 given'
			allocate (E35(NP1),sig35(NP1))
			E35 = E1
			sig35 = 0
			sig35 = sig1
 		end if
		if(MTi==36) then
			write(404,*)
			write(404,*) 'Activation and Gas production MT = 36 given'
			allocate (E36(NP1),sig36(NP1))
			E36 = E1
			sig36 = 0
			sig36 = sig1
 		end if
		if(MTi==37) then
			write(404,*)
			write(404,*) 'Activation MT = 37 given'
			allocate (E37(NP1),sig37(NP1))
			E37 = E1
			sig37 = 0
			sig37 = sig1
 		end if
		if(MTi==41) then
			write(404,*)
			write(404,*) 'Activation and Gas production MT = 41 given'
			allocate (E41(NP1),sig41(NP1))
			E41 = E1
			sig41 = 0
			sig41 = sig1
 		end if
		if(MTi==42) then
			write(404,*)
			write(404,*) 'Activation and Gas production MT = 42 given'
			allocate (E42(NP1),sig42(NP1))
			E42 = E1
			sig42 = 0
			sig42 = sig1
 		end if
		if(MTi==44) then
			write(404,*)
			write(404,*) 'Activation and Gas production MT = 44 given'
			allocate (E44(NP1),sig44(NP1))
			E44 = E1
			sig44 = 0
			sig44 = sig1
 		end if
		if(MTi==45) then
			write(404,*)
			write(404,*) 'Activation and Gas production MT = 45 given'
			allocate (E45(NP1),sig45(NP1))
			E45 = E1
			sig45 = 0
			sig45 = sig1
 		end if
		if(MTi==102) then
			write(404,*)
			write(404,*) 'Activation MT = 102 given'
			allocate (E102(NP1),sig102(NP1))
			E102 = E1
			sig102 = 0
			sig102 = sig1
 		end if
		if(MTi==103) then
			write(404,*)
			write(404,*) 'Activation and Gas production MT = 103 given'
			allocate (E103(NP1),sig103(NP1))
			E103 = E1
			sig103 = 0
			sig103 = sig1
 		end if
		if(MTi==104) then
			write(404,*)
			write(404,*) 'Activation and Gas production MT = 104 given'
			allocate (E104(NP1),sig104(NP1))
			E104 = E1
			sig104 = 0
			sig104 = sig1
 		end if
		if(MTi==105) then
			write(404,*)
			write(404,*) 'Activation and Gas production MT = 105 given'
			allocate (E105(NP1),sig105(NP1))
			E105 = E1
			sig105 = 0
			sig105 = sig1
 		end if
		if(MTi==106) then
			write(404,*)
			write(404,*) 'Activation and Gas production MT = 106 given'
			allocate (E106(NP1),sig106(NP1))
			E106 = E1
			sig106 = 0
			sig106 = sig1
 		end if
		if(MTi==107) then
			write(404,*)
			write(404,*) 'Activation and Gas production MT = 107 given'
			allocate (E107(NP1),sig107(NP1))
			E107 = E1
			sig107 = 0
			sig107 = sig1
 		end if
		if(MTi==108) then
			write(404,*)
			write(404,*) 'Activation and Gas production MT = 108 given'
			allocate (E108(NP1),sig108(NP1))
			E108 = E1
			sig108 = 0
			sig108 = sig1
 		end if
		if(MTi==109) then
			write(404,*)
			write(404,*) 'Activation and Gas production MT = 109 given'
			allocate (E109(NP1),sig109(NP1))
			E109 = E1
			sig109 = 0
			sig109 = sig1
 		end if
		if(MTi==110) then
			write(404,*)
			write(404,*) 'Activation and Gas production MT = 110 given'
			allocate (E110(NP1),sig110(NP1))
			E110 = E1
			sig110 = 0
			sig110 = sig1
 		end if
		if(MTi==111) then
			write(404,*)
			write(404,*) 'Activation and Gas production MT = 111 given'
			allocate (E111(NP1),sig111(NP1))
			E111 = E1
			sig111 = 0
			sig111 = sig1
 		end if
		if(MTi==112) then
			write(404,*)
			write(404,*) 'Activation and Gas production MT = 112 given'
			allocate (E112(NP1),sig112(NP1))
			E112 = E1
			sig112 = 0
			sig112 = sig1
 		end if
 		if(MTi==113) then
			write(404,*)
			write(404,*) 'Activation and Gas production MT = 113 given'
			allocate (E113(NP1),sig113(NP1))
			E113 = E1
			sig113 = 0
			sig113 = sig1
 		end if
		if(MTi==114) then
			write(404,*)
			write(404,*) 'Activation and Gas production MT = 114 given'
			allocate (E114(NP1),sig114(NP1))
			E114 = E1
			sig114 = 0
			sig114 = sig1
 		end if
		if(MTi==115) then
			write(404,*)
			write(404,*) 'Activation and Gas production MT = 115 given'
			allocate (E115(NP1),sig115(NP1))
			E115 = E1
			sig115 = 0
			sig115 = sig1
 		end if
		if(MTi==116) then
			write(404,*)
			write(404,*) 'Activation and Gas production MT = 116 given'
			allocate (E116(NP1),sig116(NP1))
			E116 = E1
			sig116 = 0
			sig116 = sig1
 		end if
		if(MTi==117) then
			write(404,*)
			write(404,*) 'Activation and Gas production MT = 117 given'
			allocate (E117(NP1),sig117(NP1))
			E117 = E1
			sig117 = 0
			sig117 = sig1
 		end if
			
		E1 = 0
		sig1 = 0
		deallocate(E1,sig1)
		
	end if
		
	close(unit=500)
end do
		
	!! The calculation from MT = 600 to 849 is performed when 
	!! MT = 103, ....., 107 are not given. These are performed in
	!! the subroutine "adddiscnth".
		
iflMTpr = 0
		
if (iflp103==0) then
	allocate (E103(NPt),sig103(NPt))
    call adddiscnth(iflMTpr,103,E,NPt,sig103)
	if (iflMTpr==1) then
		write(404,*)
		write(404,*) 'Discrete + Continuum (n, p) data represented'
		write(*,*) "From discrete (n,p) ....."
		write(*,*) 'MT = ', 103, '.....'
		E103 = E
		iflag(23)=1
	end if
end if
		
iflMTpr = 0
		
if (iflp104==0) then
	allocate (E104(NPt),sig104(NPt))
    call adddiscnth(iflMTpr,104,E,NPt,sig104)
	if (iflMTpr==1) then
		write(404,*)
		write(404,*) 'Discrete + Continuum (n, d) data represented'
		write(*,*) "From discrete (n,d) ....."
		write(*,*) 'MT = ', 104, '.....'
		E104 = E
		iflag(24)=1
	end if
end if
		
iflMTpr = 0
		
if (iflp105==0) then
	allocate (E105(NPt),sig105(NPt))
    call adddiscnth(iflMTpr,105,E,NPt,sig105)
	if (iflMTpr==1) then
		write(404,*)
		write(404,*) 'Discrete + Continuum (n, t) data represented'
		write(*,*) "From discrete (n,t) ....."
		write(*,*) 'MT = ', 105, '.....'
		E105 = E
		iflag(25)=1
	end if
end if
		
iflMTpr = 0
		
if (iflp106==0) then
	allocate (E106(NPt),sig106(NPt))
    call adddiscnth(iflMTpr,106,E,NPt,sig106)
	if (iflMTpr==1) then
		write(404,*)
		write(404,*) 'Discrete + Continuum (n, 3He) data represented'
		write(*,*) "From discrete (n,3He) ....."
		write(*,*) 'MT = ', 106, '.....'
		E106 = E
		iflag(26)=1
	end if
end if
		
iflMTpr = 0
		
if (iflp107==0) then
	allocate (E107(NPt),sig107(NPt))
    call adddiscnth(iflMTpr,107,E,NPt,sig107)
	if (iflMTpr==1) then
		write(404,*)
		write(404,*) 'Discrete + Continuum (n, a) data represented'
		write(*,*) "From discrete (n,a) ....."
		write(*,*) 'MT = ', 107, '.....'
		E107 = E
		iflag(27)=1
	end if
end if
				
	!! The yields for the production of charged particles are collected
	!! from MF = 6, MT = 5 in order to compute the respective gas production 
	!! contributions from cross sections in MF = 3, MT = 5.
		
call gtYMf6Mt5 (NBTp,INTrp,Nyld,Eyld,Yld,tnNBTp,tnINTrp,tnNyld,tnEyld, &
				tnYld,tnZp,itnt)
		
allocate(Eyldp(Nyld(1)),Eyldd(Nyld(2)),Eyldtr(Nyld(3)), &
		Eyld3He(Nyld(4)),EyldHe(Nyld(5)), &
		Yldp(Nyld(1)),Yldd(Nyld(2)),Yldtr(Nyld(3)), &
		Yld3He(Nyld(4)),YldHe(Nyld(5)))
	
Eyldp = Eyld(1,1:Nyld(1))
Eyldd = Eyld(2,1:Nyld(2))
Eyldtr = Eyld(3,1:Nyld(3))
Eyld3He = Eyld(4,1:Nyld(4))
EyldHe = Eyld(5,1:Nyld(5))

Yldp = Yld(1,1:Nyld(1))
Yldd = Yld(2,1:Nyld(2))
Yldtr = Yld(3,1:Nyld(3))
Yld3He = Yld(4,1:Nyld(4))
YldHe = Yld(5,1:Nyld(5))

NBTpp = NBTp(1,1:20)
NBTpd = NBTp(2,1:20)
NBTpt = NBTp(3,1:20)
NBTp3He = NBTp(4,1:20)
NBTpal = NBTp(5,1:20)

INTrpp = INTrp(1,1:20)
INTrpd = INTrp(2,1:20)
INTrpt = INTrp(3,1:20)
INTrp3He = INTrp(4,1:20)
INTrpal = INTrp(5,1:20)

write(404,*)
write(404,*) 'Multigroup activation cross sections asked in:'
write(404,*)

if (igtyp==0) write(404,*) '0 -- User-defined energy groups'
if (igtyp==1) write(404,*) '1 -- VITAMIN-J 175 groups'
if (igtyp==2) write(404,*) '2 -- 238 groups'
if (igtyp==3) write(404,*) '3 -- 198 groups'
if (igtyp==4) write(404,*) '4 -- 33 groups'
if (igtyp==5) write(404,*) '5 -- 26 groups'
if (igtyp==6) write(404,*) '6 -- DLC-2 100 groups'
if (igtyp==7) write(404,*) '7 -- 198 groups to 229 groups (200 MeV)'
if (igtyp==8) write(404,*) '8 -- 198 groups to 229 groups (150 MeV)'
if (igtyp==9) write(404,*) '9 -- 616 groups (from DEMO-HCPB-FW)'

		
write(404,*)
write(404,*)'Activation cross sections can be found in file:'
write(404,*)'GasProDatatest'
write(404,*)
		
	!! Multigrouping the required cross sections in the given energy
	!! group structure and then adding into separate collections for 
	!! particular gas production species and activation.
	
if (igtyp==0) then
	open(unit=406,file='Energy-GroupLimits.txt')
	read(406,23) Ngl
end if
if (igtyp==1) Ngl=176
if (igtyp==2) Ngl=239
if (igtyp==3) Ngl=199
if (igtyp==4) Ngl=34
if (igtyp==5) Ngl=27
if (igtyp==6) Ngl=101
if (igtyp==7) Ngl=229
if (igtyp==8) Ngl=229
if (igtyp==9) Ngl=617
	
write(*,*) 'Grouping: ', Ngl-1, 'groups .....'
	
allocate(Eg(Ngl),gsig(Ngl),gsp(Ngl),gsd(Ngl),gstr(Ngl), &
		gs3He(Ngl),gsa(Ngl),gsigt(Ngl),gsig5(Ngl),Yldpg(Ngl), &
		Ylddg(Ngl),Yldtrg(Ngl),Yld3Heg(Ngl),YldHeg(Ngl),tnYldg(Ngl), &
		sigtrack(Ngl,230))
		
sigtrack = 0

if (igtyp==0) then
	read(406,24) (Eg(i), i = Ngl,1,-1)
	close(unit=406)
end if	
if (igtyp==1) call engrp1(Eg,Ngl)
if (igtyp==2) call engrp2(Eg,Ngl)
if (igtyp==3) call engrp3(Eg,Ngl)
if (igtyp==4) call engrp4(Eg,Ngl)
if (igtyp==5) call engrp5(Eg,Ngl)
if (igtyp==6) call engrp6(Eg,Ngl)
if (igtyp==7) call engrp7(Eg,Ngl)
if (igtyp==8) call engrp8(Eg,Ngl)
if (igtyp==9) call engrp9(Eg,Ngl)
	
	Yldpg=0
	Ylddg=0
	Yldtrg=0
	Yld3Heg=0
	YldHeg=0
		
    call terpol(NBTpp,INTrpp,Nyld(1),Eyldp,Yldp,Ngl,Eg,Yldpg) 
    call terpol(NBTpd,INTrpd,Nyld(2),Eyldd,Yldd,Ngl,Eg,Ylddg)
    call terpol(NBTpt,INTrpt,Nyld(3),Eyldtr,Yldtr,Ngl,Eg,Yldtrg)
    call terpol(NBTp3He,INTrp3He,Nyld(4),Eyld3He,Yld3He,Ngl,Eg, &
	  Yld3Heg)
    call terpol(NBTpal,INTrpal,Nyld(5),EyldHe,YldHe,Ngl,Eg, &
	  YldHeg)
		
	gsig5 = 0
	gsigt = 0
	gsp = 0
	gsd = 0
	gstr = 0
	gs3He = 0
	gsa = 0
	sigtpoint = 0
	sig5tot = 0
	sppoint = 0
	sdpoint = 0
	strpoint = 0
	s3Hepoint = 0
	sapoint = 0
        
	do nreac = 1, 42
		
		sigparttot = 0
		
		MTi = MTnum(nreac)

		gsig = 0
		
		if (iflag(nreac)==1) then
		
			if(MTi==5) then
				NPtg = size(E5)
				allocate(E1p(NPtg),sig(NPtg))
				E1p = E5
				sig = sig5
				call terpolAPR(2,NPtg,E5,sig5,NPt,E,sig5tot)
			end if
			if(MTi==11) then
				NPtg = size(E11)
				allocate(E1p(NPtg),sig(NPtg))
				E1p = E11
				sig = sig11
				call terpolAPR(2,NPtg,E11,sig11,NPt,E,sigparttot)
			end if
			if(MTi==16) then
				NPtg = size(E16)
				allocate(E1p(NPtg),sig(NPtg))
				E1p = E16
				sig = sig16
				call terpolAPR(2,NPtg,E16,sig16,NPt,E,sigparttot)
			end if
			if(MTi==17) then
				NPtg = size(E17)
				allocate(E1p(NPtg),sig(NPtg))
				E1p = E17
				sig = sig17
				call terpolAPR(2,NPtg,E17,sig17,NPt,E,sigparttot)
			end if
			if(MTi==22) then
				NPtg = size(E22)
				allocate(E1p(NPtg),sig(NPtg))
				E1p = E22
				sig = sig22
				call terpolAPR(2,NPtg,E22,sig22,NPt,E,sigparttot)
			end if
			if(MTi==23) then
				NPtg = size(E23)
				allocate(E1p(NPtg),sig(NPtg))
				E1p = E23
				sig = sig23
				call terpolAPR(2,NPtg,E23,sig23,NPt,E,sigparttot)
			end if
			if(MTi==24) then
				NPtg = size(E24)
				allocate(E1p(NPtg),sig(NPtg))
				E1p = E24
				sig = sig24
				call terpolAPR(2,NPtg,E24,sig24,NPt,E,sigparttot)
			end if
			if(MTi==25) then
				NPtg = size(E25)
				allocate(E1p(NPtg),sig(NPtg))
				E1p = E25
				sig = sig25
				call terpolAPR(2,NPtg,E25,sig25,NPt,E,sigparttot)
			end if
			if(MTi==28) then
				NPtg = size(E28)
				allocate(E1p(NPtg),sig(NPtg))
				E1p = E28
				sig = sig28
				call terpolAPR(2,NPtg,E28,sig28,NPt,E,sigparttot)
			end if
			if(MTi==29) then
				NPtg = size(E29)
				allocate(E1p(NPtg),sig(NPtg))
				E1p = E29
				sig = sig29
				call terpolAPR(2,NPtg,E29,sig29,NPt,E,sigparttot)
			end if
			if(MTi==30) then
				NPtg = size(E30)
				allocate(E1p(NPtg),sig(NPtg))
				E1p = E30
				sig = sig30
				call terpolAPR(2,NPtg,E30,sig30,NPt,E,sigparttot)
			end if
			if(MTi==32) then
				NPtg = size(E32)
				allocate(E1p(NPtg),sig(NPtg))
				E1p = E32
				sig = sig32
				call terpolAPR(2,NPtg,E32,sig32,NPt,E,sigparttot)
			end if
			if(MTi==33) then
				NPtg = size(E33)
				allocate(E1p(NPtg),sig(NPtg))
				E1p = E33
				sig = sig33
				call terpolAPR(2,NPtg,E33,sig33,NPt,E,sigparttot)
			end if
			if(MTi==34) then
				NPtg = size(E34)
				allocate(E1p(NPtg),sig(NPtg))
				E1p = E34
				sig = sig34
				call terpolAPR(2,NPtg,E34,sig34,NPt,E,sigparttot)
			end if
			if(MTi==35) then
				NPtg = size(E35)
				allocate(E1p(NPtg),sig(NPtg))
				E1p = E35
				sig = sig35
				call terpolAPR(2,NPtg,E35,sig35,NPt,E,sigparttot)
			end if
			if(MTi==36) then
				NPtg = size(E36)
				allocate(E1p(NPtg),sig(NPtg))
				E1p = E36
				sig = sig36
				call terpolAPR(2,NPtg,E36,sig36,NPt,E,sigparttot)
			end if
			if(MTi==37) then
				NPtg = size(E37)
				allocate(E1p(NPtg),sig(NPtg))
				E1p = E37
				sig = sig37
				call terpolAPR(2,NPtg,E37,sig37,NPt,E,sigparttot)
			end if
			if(MTi==41) then
				NPtg = size(E41)
				allocate(E1p(NPtg),sig(NPtg))
				E1p = E41
				sig = sig41
				call terpolAPR(2,NPtg,E41,sig41,NPt,E,sigparttot)
			end if
			if(MTi==42) then
				NPtg = size(E42)
				allocate(E1p(NPtg),sig(NPtg))
				E1p = E42
				sig = sig42
				call terpolAPR(2,NPtg,E42,sig42,NPt,E,sigparttot)
			end if
			if(MTi==44) then
				NPtg = size(E44)
				allocate(E1p(NPtg),sig(NPtg))
				E1p = E44
				sig = sig44
				call terpolAPR(2,NPtg,E44,sig44,NPt,E,sigparttot)
			end if
			if(MTi==45) then
				NPtg = size(E45)
				allocate(E1p(NPtg),sig(NPtg))
				E1p = E45
				sig = sig45
				call terpolAPR(2,NPtg,E45,sig45,NPt,E,sigparttot)
			end if
			if(MTi==102) then
				NPtg = size(E102)
				allocate(E1p(NPtg),sig(NPtg))
				E1p = E102
				sig = sig102
				call terpolAPR(2,NPtg,E102,sig102,NPt,E,sigparttot)
			end if
			if(MTi==103) then
				NPtg = size(E103)
				allocate(E1p(NPtg),sig(NPtg))
				E1p = E103
				sig = sig103
				call terpolAPR(2,NPtg,E103,sig103,NPt,E,sigparttot)
			end if
			if(MTi==104) then
				NPtg = size(E104)
				allocate(E1p(NPtg),sig(NPtg))
				E1p = E104
				sig = sig104
				call terpolAPR(2,NPtg,E104,sig104,NPt,E,sigparttot)
			end if
			if(MTi==105) then
				NPtg = size(E105)
				allocate(E1p(NPtg),sig(NPtg))
				E1p = E105
				sig = sig105
				call terpolAPR(2,NPtg,E105,sig105,NPt,E,sigparttot)
			end if
			if(MTi==106) then
				NPtg = size(E106)
				allocate(E1p(NPtg),sig(NPtg))
				E1p = E106
				sig = sig106
				call terpolAPR(2,NPtg,E106,sig106,NPt,E,sigparttot)
			end if
			if(MTi==107) then
				NPtg = size(E107)
				allocate(E1p(NPtg),sig(NPtg))
				E1p = E107
				sig = sig107
				call terpolAPR(2,NPtg,E107,sig107,NPt,E,sigparttot)
			end if
			if(MTi==108) then
				NPtg = size(E108)
				allocate(E1p(NPtg),sig(NPtg))
				E1p = E108
				sig = sig108
				call terpolAPR(2,NPtg,E108,sig108,NPt,E,sigparttot)
			end if
			if(MTi==109) then
				NPtg = size(E109)
				allocate(E1p(NPtg),sig(NPtg))
				E1p = E109
				sig = sig109
				call terpolAPR(2,NPtg,E109,sig109,NPt,E,sigparttot)
			end if
			if(MTi==110) then
				NPtg = size(E110)
				allocate(E1p(NPtg),sig(NPtg))
				E1p = E110
				sig = sig110
				call terpolAPR(2,NPtg,E110,sig110,NPt,E,sigparttot)
			end if
			if(MTi==111) then
				NPtg = size(E111)
				allocate(E1p(NPtg),sig(NPtg))
				E1p = E111
				sig = sig111
				call terpolAPR(2,NPtg,E111,sig111,NPt,E,sigparttot)
			end if
			if(MTi==112) then
				NPtg = size(E112)
				allocate(E1p(NPtg),sig(NPtg))
				E1p = E112
				sig = sig112
				call terpolAPR(2,NPtg,E112,sig112,NPt,E,sigparttot)
			end if
			if(MTi==113) then
				NPtg = size(E113)
				allocate(E1p(NPtg),sig(NPtg))
				E1p = E113
				sig = sig113
				call terpolAPR(2,NPtg,E113,sig113,NPt,E,sigparttot)
			end if
			if(MTi==114) then
				NPtg = size(E114)
				allocate(E1p(NPtg),sig(NPtg))
				E1p = E114
				sig = sig114
				call terpolAPR(2,NPtg,E114,sig114,NPt,E,sigparttot)
			end if
			if(MTi==115) then
				NPtg = size(E115)
				allocate(E1p(NPtg),sig(NPtg))
				E1p = E115
				sig = sig115
				call terpolAPR(2,NPtg,E115,sig115,NPt,E,sigparttot)
			end if
			if(MTi==116) then
				NPtg = size(E116)
				allocate(E1p(NPtg),sig(NPtg))
				E1p = E116
				sig = sig116
				call terpolAPR(2,NPtg,E116,sig116,NPt,E,sigparttot)
			end if
			if(MTi==117) then
				NPtg = size(E117)
				allocate(E1p(NPtg),sig(NPtg))
				E1p = E117
				sig = sig117
				call terpolAPR(2,NPtg,E117,sig117,NPt,E,sigparttot)
			end if
			if(MTi==203) then
				NPtg = size(E203)
				allocate(E1p(NPtg),sig(NPtg))
				E1p = E203
				sig = sig203
				call terpolAPR(2,NPtg,E203,sig203,NPt,E,sigparttot)
			end if
			if(MTi==204) then
				NPtg = size(E204)
				allocate(E1p(NPtg),sig(NPtg))
				E1p = E204
				sig = sig204
				call terpolAPR(2,NPtg,E204,sig204,NPt,E,sigparttot)
			end if
			if(MTi==205) then
				NPtg = size(E205)
				allocate(E1p(NPtg),sig(NPtg))
				E1p = E205
				sig = sig205
				call terpolAPR(2,NPtg,E205,sig205,NPt,E,sigparttot)
			end if
			if(MTi==206) then
				NPtg = size(E206)
				allocate(E1p(NPtg),sig(NPtg))
				E1p = E206
				sig = sig206
				call terpolAPR(2,NPtg,E206,sig206,NPt,E,sigparttot)
			end if
			if(MTi==207) then
				NPtg = size(E207)
				allocate(E1p(NPtg),sig(NPtg))
				E1p = E207
				sig = sig207
				call terpolAPR(2,NPtg,E207,sig207,NPt,E,sigparttot)
			end if
			
		!! Cross section gets multigrouped based on n-spectrum here .....
			
			call groupmulti (insp,E1p,sig,NPtg,Eg,Ngl,gsig) 
		
			if (MTi==5) then
				gsig5 = gsig
				do i = 1, Ngl
					gsp(i) = gsp(i)+Yldpg(i)*gsig(i)
					gsd(i) = gsd(i)+Ylddg(i)*gsig(i)
					gstr(i) = gstr(i)+Yldtrg(i)*gsig(i)
					gs3He(i) = gs3He(i)+Yld3Heg(i)*gsig(i)
					gsa(i) = gsa(i)+YldHeg(i)*gsig(i)
					gsigt(i) = gsigt(i) + gsig(i)
				end do
				deallocate(Yldpg,Ylddg,Yldtrg,Yld3Heg,YldHeg)
				allocate(Yldpg(NPt),Ylddg(NPt),Yldtrg(NPt), &
				Yld3Heg(NPt),YldHeg(NPt))
				Yldpg=0
				Ylddg=0
				Yldtrg=0
				Yld3Heg=0
				YldHeg=0
				
			call terpol(NBTpp,INTrpp,Nyld(1),Eyldp,Yldp,NPt,E,Yldpg) 
			call terpol(NBTpd,INTrpd,Nyld(2),Eyldd,Yldd,NPt,E,Ylddg)
			call terpol(NBTpt,INTrpt,Nyld(3),Eyldtr,Yldtr,NPt,E,Yldtrg)
			call terpol(NBTp3He,INTrp3He,Nyld(4),Eyld3He,Yld3He,NPt,E, &
				Yld3Heg)
			call terpol(NBTpal,INTrpal,Nyld(5),EyldHe,YldHe,NPt,E, &
				YldHeg)
				
				do i = 1, NPt
					sppoint(i) = sppoint(i)+Yldpg(i)*sig5tot(i)
					sdpoint(i) = sdpoint(i)+Ylddg(i)*sig5tot(i)
					strpoint(i) = strpoint(i)+Yldtrg(i)*sig5tot(i)
					s3Hepoint(i) = s3Hepoint(i)+Yld3Heg(i)*sig5tot(i)
					sapoint(i) = sapoint(i)+YldHeg(i)*sig5tot(i)
					sigtpoint(i) = sigtpoint(i) + sig5tot(i)
				end do
			end if
		
			if (iflMTtppr(1) == 0) then
			if (MTi==28.or.MTi==41.or.MTi==42.or.MTi==44.or.MTi==45.or. &
			MTi==103.or.MTi==111.or.MTi==112.or.MTi==115.or.MTi==116)then 
				do i = 1, Ngl
					gsp(i) = gsp(i)+gsig(i)
					gsigt(i) = gsigt(i) + gsig(i)
				end do
				do i = 1, NPt
					sppoint(i) = sppoint(i) + sigparttot(i)
					sigtpoint(i) = sigtpoint(i) + sigparttot(i)
				end do
			end if
			end if 
			
			if (MTi==203) then
				do i = 1, Ngl
					gsp(i) = gsig(i)
					gsigt(i) = gsigt(i) + gsig(i)
				end do
				do i = 1, NPt
					sppoint(i) = sigparttot(i)
					sigtpoint(i) = sigtpoint(i) + sigparttot(i)
				end do
			end if
			
			if (iflMTtppr(2) == 0) then
			if (MTi==11.or.MTi==32.or.MTi==35.or.MTi==104.or. &
			MTi==114.or.MTi==115.or.MTi==117)then 
				do i = 1, Ngl
					gsd(i) = gsd(i)+gsig(i)
					gsigt(i) = gsigt(i) + gsig(i)
				end do
				do i = 1, NPt
					sdpoint(i) = sdpoint(i) + sigparttot(i)
					sigtpoint(i) = sigtpoint(i) + sigparttot(i)
				end do
			end if
			end if
			
			if (MTi==204) then
				do i = 1, Ngl
					gsd(i) = gsig(i)
					gsigt(i) = gsigt(i) + gsig(i)
				end do
				do i = 1, NPt
					sdpoint(i) = sigparttot(i)
					sigtpoint(i) = sigtpoint(i) + sigparttot(i)
				end do
			end if
			
			if (iflMTtppr(3) == 0) then
			if (MTi==33.or.MTi==36.or.MTi==105.or.MTi==113.or.MTi==116)then 
				do i = 1, Ngl
					gstr(i) = gstr(i)+gsig(i)
					gsigt(i) = gsigt(i) + gsig(i)
				end do
				do i = 1, NPt
					strpoint(i) = strpoint(i) + sigparttot(i)
					sigtpoint(i) = sigtpoint(i) + sigparttot(i)
				end do
			end if
			end if
			
			if (MTi==205) then
				do i = 1, Ngl
					gstr(i) = gsig(i)
					gsigt(i) = gsigt(i) + gsig(i)
				end do
				do i = 1, NPt
					strpoint(i) = sigparttot(i)
					sigtpoint(i) = sigtpoint(i) + sigparttot(i)
				end do
			end if
			
			if (iflMTtppr(4) == 0) then
			if (MTi==34.or.MTi==106)then 
				do i = 1, Ngl
					gs3He(i) = gs3He(i)+gsig(i)
					gsigt(i) = gsigt(i) + gsig(i)
				end do
				do i = 1, NPt
					s3Hepoint(i) = s3Hepoint(i) + sigparttot(i)
					sigtpoint(i) = sigtpoint(i) + sigparttot(i)
				end do
			end if
			end if
			
			if (MTi==206) then
				do i = 1, Ngl
					gs3He(i) = gsig(i)
					gsigt(i) = gsigt(i) + gsig(i)
				end do
				do i = 1, NPt
					s3Hepoint(i) = sigparttot(i)
					sigtpoint(i) = sigtpoint(i) + sigparttot(i)
				end do
			end if
			
			if (iflMTtppr(5) == 0) then
			if (MTi==22.or.MTi==23.or.MTi==24.or.MTi==25.or.MTi==29.or. &
			MTi==30.or.MTi==35.or.MTi==36.or.MTi==45.or.MTi==107.or. &
			MTi==108.or.MTi==109.or.MTi==112.or.MTi==113.or.MTi==114.or. &
			MTi==117) then 
				do i = 1, Ngl
					gsa(i) = gsa(i)+gsig(i)
					gsigt(i) = gsigt(i) + gsig(i)
				end do
				do i = 1, NPt
					sapoint(i) = sapoint(i) + sigparttot(i)
					sigtpoint(i) = sigtpoint(i) + sigparttot(i)
				end do
			end if
			end if
			
			if (MTi==207) then
				do i = 1, Ngl
					gsa(i) = gsig(i)
					gsigt(i) = gsigt(i) + gsig(i)
				end do
				do i = 1, NPt
					sapoint(i) = sigparttot(i)
					sigtpoint(i) = sigtpoint(i) + sigparttot(i)
				end do
			end if
			
		!! For the n,g in Ni58 and Mo98
			if (eliso=='Ni58'.or.eliso=='Mo98') then
			open(unit=4000,file="ng_Ni58&Mo98.txt",position='append')
			if (MTi==102) then
				write(4000,*) eliso
				write(4000,12) (gsig(i), i = 1, Ngl)
			end if
			close(unit=4000)
			end if
			
			if (MTi==16 .or. MTi==17 .or. MTi==37 .or. MTi==102) then
			do i = 1, Ngl
				gsigt(i) = gsigt(i) + gsig(i)
			end do
			do i = 1, NPt
				sigtpoint(i) = sigtpoint(i) + sigparttot(i) 
			end do
			end if
			
			write(*,*) 'MT = ',MTi,'..... done'
			
			deallocate(E1p,sig)
		
			if (MTi.ne.5) then
			if (MTi==11.or.MTi==33.or.MTi==42) then
				Zval = Ztarget-1 
				Aval = Atarget-3
				Zvaltrack(1) = Zval
				Avaltrack(1) = Aval
				cntrack(1) = cntrack(1)+1
				do k = 1, Ngl
					sigtrack(k,1) = sigtrack(k,1) + gsig(k)
				end do
			end if
			if (MTi==16) then
				Zval = Ztarget
				Aval = Atarget-1
				Zvaltrack(2) = Zval
				Avaltrack(2) = Aval
				cntrack(2) = cntrack(2)+1
				do k = 1, Ngl
					sigtrack(k,2) = sigtrack(k,2) + gsig(k)
				end do
			end if
			if (MTi==17) then
				Zval = Ztarget
				Aval = Atarget-2
				Zvaltrack(3) = Zval
				Avaltrack(3) = Aval
				cntrack(3) = cntrack(3)+1
				do k = 1, Ngl
					sigtrack(k,3) = sigtrack(k,3) + gsig(k)
				end do
			end if
			if (MTi==22) then
				Zval = Ztarget-2
				Aval = Atarget-4
				Zvaltrack(4) = Zval
				Avaltrack(4) = Aval
				cntrack(4) = cntrack(4)+1
				do k = 1, Ngl
					sigtrack(k,4) = sigtrack(k,4) + gsig(k)
				end do
			end if
			if (MTi==23) then
				Zval = Ztarget-6
				Aval = Atarget-12
				Zvaltrack(5) = Zval
				Avaltrack(5) = Aval
				cntrack(5) = cntrack(5)+1
				do k = 1, Ngl
					sigtrack(k,5) = sigtrack(k,5) + gsig(k)
				end do
			end if
			if (MTi==24) then
				Zval = Ztarget-2
				Aval = Atarget-5
				Zvaltrack(6) = Zval
				Avaltrack(6) = Aval
				cntrack(6) = cntrack(6)+1
				do k = 1, Ngl
					sigtrack(k,6) = sigtrack(k,6) + gsig(k)
				end do
			end if
			if (MTi==25) then
				Zval = Ztarget-2
				Aval = Atarget-6
				Zvaltrack(7) = Zval
				Avaltrack(7) = Aval
				cntrack(7) = cntrack(7)+1
				do k = 1, Ngl
					sigtrack(k,7) = sigtrack(k,7) + gsig(k)
				end do
			end if
			if (MTi==28.or.MTi==104) then
				Zval = Ztarget-1
				Aval = Atarget-1
				Zvaltrack(8) = Zval
				Avaltrack(8) = Aval
				cntrack(8) = cntrack(8)+1
				do k = 1, Ngl
					sigtrack(k,8) = sigtrack(k,8) + gsig(k)
				end do
			end if
			if (MTi==29) then 
				Zval = Ztarget-4
				Aval = Atarget-8
				Zvaltrack(9) = Zval
				Avaltrack(9) = Aval
				cntrack(9) = cntrack(9)+1
				do k = 1, Ngl
					sigtrack(k,9) = sigtrack(k,9) + gsig(k)
				end do
			end if
			if (MTi==30) then 
				Zval = Ztarget-4
				Aval = Atarget-9
				Zvaltrack(10) = Zval
				Avaltrack(10) = Aval
				cntrack(10) = cntrack(10)+1
				do k = 1, Ngl
					sigtrack(k,10) = sigtrack(k,10) + gsig(k)
				end do
			end if
			if (MTi==32.or.MTi==41.or.MTi==105) then
				Zval = Ztarget-1
				Aval = Atarget-2
				Zvaltrack(11) = Zval
				Avaltrack(11) = Aval
				cntrack(11) = cntrack(11)+1
				do k = 1, Ngl
					sigtrack(k,11) = sigtrack(k,11) + gsig(k)
				end do
			end if
			if (MTi==34.or.MTi==107.or.MTi==116) then
				Zval = Ztarget-2
				Aval = Atarget-3
				Zvaltrack(12) = Zval
				Avaltrack(12) = Aval
				cntrack(12) = cntrack(12)+1
				do k = 1, Ngl
					sigtrack(k,12) = sigtrack(k,12) + gsig(k)
				end do
			end if
			if (MTi==35.or.MTi==113) then
				Zval = Ztarget-5
				Aval = Atarget-10
				Zvaltrack(13) = Zval
				Avaltrack(13) = Aval
				cntrack(13) = cntrack(13)+1
				do k = 1, Ngl
					sigtrack(k,13) = sigtrack(k,13) + gsig(k)
				end do
			end if
			if (MTi==36) then
				Zval = Ztarget-5
				Aval = Atarget-11
				Zvaltrack(14) = Zval
				Avaltrack(14) = Aval
				cntrack(14) = cntrack(14)+1
				do k = 1, Ngl
					sigtrack(k,14) = sigtrack(k,14) + gsig(k)
				end do
			end if
			if (MTi==37) then
				Zval = Ztarget
				Aval = Atarget-3
				Zvaltrack(15) = Zval
				Avaltrack(15) = Aval
				cntrack(15) = cntrack(15)+1
				do k = 1, Ngl
					sigtrack(k,15) = sigtrack(k,15) + gsig(k)
				end do
			end if
			if (MTi==44.or.MTi==106.or.MTi==115) then
				Zval = Ztarget-2
				Aval = Atarget-2
				Zvaltrack(16) = Zval
				Avaltrack(16) = Aval
				cntrack(16) = cntrack(16)+1
				do k = 1, Ngl
					sigtrack(k,16) = sigtrack(k,16) + gsig(k)
				end do
			end if
			if (MTi==45.or.MTi==117) then
				Zval = Ztarget-3
				Aval = Atarget-5
				Zvaltrack(17) = Zval
				Avaltrack(17) = Aval
				cntrack(17) = cntrack(17)+1
				do k = 1, Ngl
					sigtrack(k,17) = sigtrack(k,17) + gsig(k)
				end do
			end if
			if (MTi==102) then
				Zval = Ztarget
				Aval = Atarget+1
				Zvaltrack(18) = Zval
				Avaltrack(18) = Aval
				cntrack(18) = cntrack(18)+1
				do k = 1, Ngl
					sigtrack(k,18) = sigtrack(k,18) + gsig(k)
				end do 
			end if
			if (MTi==103) then
				Zval = Ztarget-1
				Aval = Atarget
				Zvaltrack(19) = Zval
				Avaltrack(19) = Aval
				cntrack(19) = cntrack(19)+1
				do k = 1, Ngl
					sigtrack(k,19) = sigtrack(k,19) + gsig(k)
				end do
			end if
			if (MTi==108) then
				Zval = Ztarget-4
				Aval = Atarget-7
				Zvaltrack(20) = Zval
				Avaltrack(20) = Aval
				cntrack(20) = cntrack(20)+1
				do k = 1, Ngl
					sigtrack(k,20) = sigtrack(k,20) + gsig(k)
				end do
			end if
			if (MTi==109) then
				Zval = Ztarget-6
				Aval = Atarget-11
				Zvaltrack(21) = Zval
				Avaltrack(21) = Aval
				cntrack(21) = cntrack(21)+1
				do k = 1, Ngl
					sigtrack(k,21) = sigtrack(k,21) + gsig(k)
				end do
			end if
			if (MTi==111) then
				Zval = Ztarget-2
				Aval = Atarget-1
				Zvaltrack(22) = Zval
				Avaltrack(22) = Aval
				cntrack(22) = cntrack(22)+1
				do k = 1, Ngl
					sigtrack(k,22) = sigtrack(k,22) + gsig(k)
				end do
			end if
			if (MTi==112) then
				Zval = Ztarget-3
				Aval = Atarget-4
				Zvaltrack(23) = Zval
				Avaltrack(23) = Aval
				cntrack(23) = cntrack(23)+1
				do k = 1, Ngl
					sigtrack(k,23) = sigtrack(k,23) + gsig(k)
				end do
			end if
			if (MTi==114) then
				Zval = Ztarget-5
				Aval = Atarget-9
				Zvaltrack(24) = Zval
				Avaltrack(24) = Aval
				cntrack(24) = cntrack(24)+1
				do k = 1, Ngl
					sigtrack(k,24) = sigtrack(k,24) + gsig(k)
				end do
			end if
				write(401,21) MTi, eliso, Zval, Aval
				write(401,20)(Eg(k), gsig(k), k = 1, Ngl)
				write(401,*)'-----------------'
			end if
		end if
	end do
		
		write(400,*) eliso, Ngl
		write(400,4) (Eg(k),gsp(k),gsd(k),gstr(k),gs3He(k),gsa(k), &
		gsigt(k),k = 1, Ngl)
		write(400,*) '-------------------'
		
!! Print the Point Cross Sections (Energies corresponding
!! to the MF=3 MT=1) of Production of Light Charged Particles

		write(402,*) eliso, NPt
		
		write(402,4) (E(k),sppoint(k),sdpoint(k),strpoint(k), &
		s3Hepoint(k),sapoint(k),sigtpoint(k), k = 1, NPt)
		write(402,*) '-------------------'
	
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!! Find and Print the Multigrouped cross sections of Transmuted Nuclides
	
	do itn = 1, itnt
		tnYldg = 0
		itnZ = tnZp(itn)/1000
		itnA =  mod(tnZp(itn),1000)
		
		call terpol(tnNBTp(itn,1:20),tnINTrp(itn,1:20),tnNyld(itn), &
				tnEyld(itn,1:200),tnYld(itn,1:200),Ngl,Eg,tnYldg)
				
		iflsametrnsN = 0
		do i = 1, 230
			if (Zvaltrack(i)==itnZ .and. Avaltrack(i)==itnA) then
				cntrack(i) = cntrack(i) + 1
				iflsametrnsN = 1
				do k = 1, Ngl
					sigtrack(k,i) = sigtrack(k,i) + tnYldg(k)*gsig5(k)
				end do
				exit
			end if
		end do
		if (iflsametrnsN == 0) then
		do i = 1, 230
			if (Zvaltrack(i)==0 .and. Avaltrack(i)==0) then
				Zvaltrack(i) = itnZ
				Avaltrack(i) = itnA
				cntrack(i) = cntrack(i) + 1
				do k = 1, Ngl
					sigtrack(k,i) = sigtrack(k,i) + tnYldg(k)*gsig5(k)
				end do
				exit
			end if
		end do
		end if
				
		write(401,21)itn,eliso,itnZ,itnA
		write(401,20)(Eg(k), tnYldg(k)*gsig5(k), k = 1, Ngl)
		write(401,*)'-----------------'
		
	end do
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
write(401,*) 'List of Transmuted Nuclides Produced and Number of Times Each is Found'
write(401,*) 'Sl. No.Z/A/Counts'
!! Find Net Number of Transmuted Nuclides
do i = 230, 1, -1
	if (Zvaltrack(i).ne.0.and.Avaltrack(i).ne.0.and.cntrack(i).ne.0) then
		itrnstot = i
		exit
	end if
end do 
write(401,22) (i,Zvaltrack(i), Avaltrack(i), cntrack(i), i = 1, itrnstot+1)
write(401,*) '------------------'

do i = 1, itrnstot
	write(405,21) i, eliso, Zvaltrack(i), Avaltrack(i)
	write(405,20) (Eg(k), sigtrack(k,i), k = 1, Ngl)
	write(405,*) '------------------'
end do

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	! To calculate the concentrations of transmuted species at 
	! different times of neutron irradiation
if (iconcyn == 1) then
write(*,*) 'Calculating concentrations of transmuted species .....'
call transpNeuspec(eliso,Ngl,itrnstot,Zvaltrack,Avaltrack,sigtrack, &
gsp,gsd,gstr,gs3He,gsa,gsigt)
write(*,*) '..... done'
end if
	
	deallocate(Eg,gsig,gsp,gsd,gstr,gs3He,gsa,gsigt,gsig5,tnYldg)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!! Find and Print the Point Cross Sections (Energies corresponding
!! to the MF=3 MT=1) of Transmuted Nuclides

allocate(tnYldg(NPt))

do itn = 1, itnt
	tnYldg = 0
	itnZ = tnZp(itn)/1000
	itnA =  mod(tnZp(itn),1000)
	call terpol(tnNBTp(itn,1:20),tnINTrp(itn,1:20),tnNyld(itn), &
			tnEyld(itn,1:200),tnYld(itn,1:200),NPt,E,tnYldg)
	write(403,*)itn,eliso,itnZ,itnA
	do k = 1, NPt
		if (tnYldg(k).ne.0) then
			kstart = k-1
			exit
		end if
	end do 
	write(403,20)(E(k), tnYldg(k)*sig5tot(k), k = kstart, NPt)
	write(403,*)'------------------'
end do

deallocate(tnYldg)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1 	format(66x,I4,I2,I3,I5)
2	format(2e11.0,11x,3i11)
3	format(6e11.0)
5	format(66x,i4,i2,i3)
6	format(2e11.0,22x,2i11)
7	format(11x,i11)
8	format(6e11.0)
11 	format(2e11.0,4I11,I4,I2,I3,I5) 
4 	format(1p7e13.6)
12 	format(1p1e13.6)
20 	format(1p2e13.6)
21	format(i4,a8,i5,i6)
22 	format(4i5)
23	format(i11)
24	format(e11.0)

close (unit=400)
close (unit=401)
close (unit=402)
close (unit=403)
close (unit=404)
close (unit=405)

end subroutine TransmCal
		
!=======================================================================
	
	!! To energy multigroup cross sections depending on neutron spectrum
	
subroutine groupmulti (insp,E1p,sig,NPtg,Eg,Ngl,gsig)
real,dimension(:):: E1p(NPtg), sig(NPtg), Eg(Ngl), gsig(Ngl)
real,allocatable,dimension(:):: fi

if (insp == 1) then
	allocate(fi(Ngl))
	fi = 0
	open(unit=407,file='NeutronSpectrum.txt')
	read(407,23) nre
	read(407,24) (fi(i), i = nre,1,-1)
	close(unit=407)
end if

gsig = 0

do i = 1, Ngl-1
	if (Eg(i)<=E1p(1).and.E1p(1)<=Eg(i+1)) then
		ifg = i
		exit
	end if
end do
k1=1;
do i = ifg, Ngl-1
	Eg1 = Eg(i)
	Eg2 = Eg(i+1)
	Nsect = 10000
	h = (Eg2-Eg1)/Nsect
	s1 = 0; f1 = 0; f2 = 0
	s2 = 0; g1 = 0; g2 = 0
	!if (igrps==3) h = (Eg2-Eg1)/1000
	do j = 1, Nsect
		t = Eg1 + (j-1)*h
		do k = k1,NPtg
			if (insp == 0) then
				if (E1p(k)<0.1) then
				L = 1
				end if
				if (0.1<=E1p(k).and.E1p(k)<820.3e+03) then
				L = 2
				end if
				if (E1p(k)>=820.3e+03) then
				L = 3
				end if
			end if
			if (E1p(k)==t)then
				if (insp == 0) flux1 = spectrum (t,L)
				if (insp == 1) flux1 = srchintrp3 (Eg,fi,Ngl,t)
				f1=sig(k)*flux1
				g1 = flux1
				k1=k
				exit
			end if
			if (E1p(k)<t .and. t<E1p(k+1))then
				if (insp == 0) flux1 = spectrum (t,L)
				if (insp == 1) flux1 = srchintrp3 (Eg,fi,Ngl,t)
				y=crstd(t,E1p(k),E1p(k+1),sig(k),sig(k+1))
				f1=y*flux1
				g1 = flux1
				k1=k
				exit
			end if
		end do
		do k = k1,NPtg
			if (insp == 0) then
				if (E1p(k)<0.1) then
				L = 1
				end if
				if (0.1<=E1p(k).and.E1p(k)<820.3e+03) then
				L = 2
				end if
				if (E1p(k)>=820.3e+03) then
				L = 3
				end if
			end if
			if (E1p(k)==t+h)then
				if (insp == 0) flux2 = spectrum (t+h,L)
				if (insp == 1) flux2 = srchintrp3 (Eg,fi,Ngl,t+h)
				f2 = sig(k)*flux2
				g2 = flux2
				k1=k
				exit
			end if
			if (E1p(k)<t+h .and. t+h<E1p(k+1))then
				if (insp == 0) flux2 = spectrum (t+h,L)
				if (insp == 1) flux2 = srchintrp3 (Eg,fi,Ngl,t+h)
				y=crstd(t+h,E1p(k),E1p(k+1),sig(k),sig(k+1))
				f2=y*flux2
				g2 = flux2
				k1=k
				exit
			end if
		end do
		s1=s1+(h/2)*(f1+f2)
		s2=s2+(h/2)*(g1+g2)
	end do
	if (s1.ne.0 .and. s2.ne.0) gsig(i) = s1/s2
end do

23	format(i11)
24	format(e11.0)

end subroutine groupmulti

!=======================================================================
	
	!! To calculate concentration (appm) of transmuted species for a 
	!! particular incident neutron spectrum

subroutine transpNeuspec(eliso,Ngl,itrnstot,Ztrans,Atrans,sigtrans, &
gsp,gsd,gstr,gs3He,gsa,gsigt)

character(8):: eliso
real,dimension(Ngl,itrnstot):: sigtrans
real,dimension(itrnstot):: transsigone 
real,dimension(Ngl):: gsp,gsd,gstr,gs3He,gsa,gsigt,fiinc
integer,dimension(itrnstot):: Ztrans, Atrans
real,allocatable,dimension(:):: times,concp,concd,conctr,conc3He,conca
real,allocatable,dimension(:,:):: conctrans
integer,allocatable,dimension(:):: itime

open(unit=400, file="NeutronSpectrum.txt")
open(unit=401, file="TranspNeuspec-input.txt")
open(unit=402, file="Concentrations-TransmutedSpecies.txt")

data btcm2/1e-24/
data fmln/1e+06/

write(*,*) 'Reading neutron spectrum .....'
read(400,1) nre
if (nre.ne.(Ngl-1)) then
	write(*,*) 'Given spectrum is not as per requirement.'
	write(*,*) 'Aborting calculation of concentrations.'
	goto 100
end if

fiinc = 0
read(400,2) (fiinc(i), i = nre,1,-1)

close(unit=400)

fisum = sum(fiinc)

write(*,*) 'Reading isotopic mass and irradiation times .....'
read(401,*)
read(401,*) amass
read(401,*)
read(401,*) ntimes

allocate(itime(ntimes),times(ntimes),concp(ntimes),concd(ntimes), &
conctr(ntimes),conc3He(ntimes),conca(ntimes),conctrans(itrnstot,ntimes))

concp = 0
concd = 0
conctr = 0
conc3He = 0
conca = 0
conctrans = 0

do i = 1, ntimes
	read(401,*) itime(i)
	times(i) = itime(i)*86400 	! itime is in days, times is in seconds
end do
close(unit=401)

	!! Finding one-group cross sections
psigone = 0
dsigone = 0
tsigone = 0
tHesigone = 0
alsigone = 0
absigone = 0
transsigone = 0

do i = 1, Ngl
	psigone = psigone + gsp(i)*fiinc(i)
	dsigone = dsigone + gsd(i)*fiinc(i)
	tsigone = tsigone + gstr(i)*fiinc(i)
	tHesigone = tHesigone + gs3He(i)*fiinc(i)
	alsigone = alsigone + gsa(i)*fiinc(i)
	absigone = absigone + gsigt(i)*fiinc(i)
	do j = 1, itrnstot
		transsigone(j) = transsigone(j) + sigtrans(i,j)*fiinc(i)
	end do
end do

!! These are total reaction rates : (sigma*fi) = sum(sigma_g * fi_g)
psigone1 = psigone*btcm2
dsigone1 = dsigone*btcm2
tsigone1 = tsigone*btcm2
tHesigone1 = tHesigone*btcm2
alsigone1 = alsigone*btcm2
absigone1 = absigone*btcm2

!! These are neutron spectrum-averaged one group cross sections in barns
psigone = psigone/fisum
dsigone = dsigone/fisum
tsigone = tsigone/fisum
tHesigone = tHesigone/fisum
alsigone = alsigone/fisum
absigone = absigone/fisum



do j = 1, itrnstot
	transsigone(j) = transsigone(j)*btcm2
end do

do i = 1, ntimes
	factintg = trspintegr(absigone1,times(i))
	concp(i) = psigone1*factintg*(fmln*amass)
	concd(i) = dsigone1*factintg*(fmln*amass)
	conctr(i) = tsigone1*factintg*(fmln*amass)
	conc3He(i) = tHesigone1*factintg*(fmln*amass)
	conca(i) = alsigone1*factintg*(fmln*amass)
	do j = 1, itrnstot
	if(transsigone(j).ne.0) conctrans(j,i) = transsigone(j)*factintg*(fmln*amass)
	end do
end do

write(*,*) 'Writing .....'
write(402,*) 'Calculation for:: ', eliso
write(402,*) 'Total Flux of the Incident Particle: ', fisum
write(402,*)
write(402,*) 'Neutron Spectrum-Averaged Onegroup Cross sections (barns)'
write(402,*)'Hydrogen/ Deuterium/ Tritium/ 3Helium/ Helium/ Activation'
write(402,3) psigone,dsigone,tsigone,tHesigone,alsigone,absigone 
write(402,*)
write(402,4) (Ztrans(i), Atrans(i), transsigone(i)/(btcm2*fisum), i = 1, itrnstot)
write(402,*) '--------------------------'
write(402,*)
write(402,*) 'Concentrations of Transmuted species (appm)'
write(402,*)'Time(d)/Hydrogen/ Deuterium/ Tritium/ 3Helium/ Helium'
write(402,5) (itime(i),concp(i),concd(i),conctr(i),conc3He(i),conca(i), &
i = 1, ntimes)
do i = 1, itrnstot
	if (transsigone(i).ne.0) then
		write(402,6) Ztrans(i), Atrans(i)
		write(402,7) (itime(j), conctrans(i,j), j = 1, ntimes)
	end if
end do

close(unit=402)

1	format(i11)
2	format(e11.0)
3	format(6e13.6)
4 	format(2i5,e13.6)
5 	format(i8,5e13.6)
6 	format(2i5)
7 	format(i8,e13.6)

100 end subroutine transpNeuspec

!================================================================

	!! Integrate for concentration of transmuted species
	
function trspintegr(sig,time)
	trspintegr = 0
	if (sig.ne.0) trspintegr = (1-exp(-sig*time))/sig
end function trspintegr

!================================================================

	!! To interpolate into required energy points using the specific
	!! scheme in the range.
	
subroutine terpol (NRin,Intrin,n1,E1,Y1,n2,E2,Y2)
real,dimension(n1):: E1,Y1
real,dimension(n2):: E2,Y2
integer, dimension(20):: NRin,Intrin
		
do i = 1, n2
	do j = 1, n1
		if (E2(i)==E1(j)) then
			Y2(i) = Y1(j)
			exit
		end if
		if (E1(j)<E2(i).and.E2(i)<E1(j+1)) then
			diff1 = E2(i) - E1(j)		
 			diff2 = E2(i) - E2(i-1)
				
			do k = 1, 20
				if (j<=NRin(k)) then
					intflg = Intrin(k)
					exit
				end if
			end do
				
 			if (diff1 <= diff2) then
				Y2(i)=terpolin(intflg,E2(i),E1(j),E1(j+1),Y1(j),Y1(j+1))
 			end if
 			if (diff2 < diff1) then
			Y2(i)=terpolin(intflg,E2(i),E2(i-1),E1(j+1),Y2(i-1),Y1(j+1))		
			end if
			exit
 		end if
	end do
end do
		
end subroutine terpol
!===============================================================
	!! Schemes of interpolation
	
function terpolin (intflg,x,x1,x2,y1,y2)
	is = intflg
	if (is==1 .or. is==11 .or. is==21) then
		y = y1
	end if
	if (is==2 .or. is==12 .or. is==22) then
		y = y1+(y2-y1)*(x-x1)/(x2-x1)
	end if
	if (is==3 .or. is==13 .or. is==23) then
		y = y1+(y2-y1)*(log(x/x1))/log(x2/x1)
	end if
	if (is==4 .or. is==14 .or. is==24) then
		y = y1*exp((x-x1)*log(y2/y1)/(x2-x1))
	end if
	if (is==5 .or. is==15 .or. is==25) then
		y = y1*exp(log(x/x1)*log(y2/y1)/log(x2/x1))
	end if
	terpolin = y
end function terpolin

!===============================================================

	!! Interpolation of Arrays via Particular Rule 
	
subroutine terpolAPR (iprule,n1,x1,y1,n2,x2,y2)
real,dimension(n1):: x1,y1
real,dimension(n2):: x2,y2
do i = 1, n2
	do j = 1, n1
		if (x2(i)==x1(j)) then
			y2(i) = y1(j)
			exit
		end if
		if (x1(j)<x2(i).and.x2(i)<x1(j+1)) then
			diff1 = x2(i) - x1(j)		
 			diff2 = x2(i) - x2(i-1)
				
 			if (diff1 <= diff2) then
				y2(i)=terpolin(iprule,x2(i),x1(j),x1(j+1),y1(j),y1(j+1))
 			end if
 			if (diff2 < diff1) then
			y2(i)=terpolin(iprule,x2(i),x2(i-1),x1(j+1),y2(i-1),y1(j+1))		
			end if
			exit
 		end if
	end do
end do
		
end subroutine terpolAPR

!====================================================================
	
	!! To collect the yields of various transmutation light charged 
	!! particles and heavy nuclides
		
subroutine gtYMf6Mt5 (NBTp,INTrp,Nyld,Eyld,Yld,tnNBTp,tnINTrp, &
        tnNyld,tnEyld,tnYld,tnZp,itnt)
real,dimension(5,200):: Eyld, Yld
real,dimension(200,200):: tnEyld, tnYld
integer,dimension(5):: iZp, iAp, Nyld
integer,dimension(200):: tnZp, tnNyld
double precision Bvall
integer, dimension(5,20):: NBTp, INTrp
integer,dimension(200,20):: tnNBTp, tnINTrp
		
ifl6mt5 = 0
Nyld = 0
Eyld = 0 
Yld = 0

iZp(1) = 1001
iZp(2) = 1002
iZp(3) = 1003
iZp(4) = 2003
iZp(5) = 2004

iAp(1) = 1
iAp(2) = 2
iAp(3) = 3
iAp(4) = 3
iAp(5) = 4

open(unit=100,file='tape01')
do
   	read(100,10) MAT, MF, MT, NS
	if (MAT==-1) then
		ifl6mt5 = 0
		exit
	end if
  	if (MT==0) then	
		if (MF==6) then
			read(100,11)ZAv,AWRv,l1,LCT,NKv,l2, MAT, MF, MT, NS
		end if
		if (MF<6) then
			read(100,*)
			read(100,11)ZAv,AWRv,l1,LCT,NKv,l2, MAT, MF, MT, NS
  		end if
	end if
				
  	if (MAT.ne.-1) then
   		if (MF.eq.6) then
			if (MT.eq.5) then 
				ifl6mt5 = 1
		
				NK = NKv
		
		!allocate(tnEyld(NK-5,200),tnYld(NK-5,200),tnZp(NK-5), &
		!tnAp(NK-5),tnNyld(NK-5),tnNBTp(NK-5,20), tnINTrp(NK-5,20))
		
				itn = 0
				do NSS = 1, NK
					read(100,11)ZAP,AWP,LIP,LAW,NR6,NP6,MAT,MF,MT,NS	
					iflgp = 0
					if(int(ZAP)==iZp(1).or.int(ZAP)==iZp(2) &
		.or.int(ZAP)==iZp(3).or.int(ZAP)==iZp(4).or.int(ZAP)==iZp(5) &
		.or.int(ZAP)>iZp(5)) iflgp = 1
		
			if (int(ZAP)==iZp(1).and.int(ceiling(AWP))==iAp(1)) ip = 1
			if (int(ZAP)==iZp(2).and.int(ceiling(AWP))==iAp(2)) ip = 2
			if (int(ZAP)==iZp(3).and.int(ceiling(AWP))==iAp(3)) ip = 3
			if (int(ZAP)==iZp(4).and.int(ceiling(AWP))==iAp(4)) ip = 4
			if (int(ZAP)==iZp(5).and.int(ceiling(AWP))==iAp(5)) ip = 5
			if (int(ZAP)>iZp(5).and.int(ceiling(AWP))>iAp(5)) then
				ip = 6
				itn = itn + 1
			end if
		
					if (iflgp == 0) then
						read(100,12) (NBT, INTr, N=1, NR6)
						read(100,13) (E, Y, N = 1, NP6)
					end if
		
					if (iflgp == 1) then 
						if (ip < 6) then
							Nyld(ip) = NP6
						read(100,12) (NBTp(ip,N), INTrp(ip,N), N=1, NR6)
						read(100,13) (Eyld(ip,N), Yld(ip,N), N = 1, NP6)
						end if
						if (ip >= 6) then
						tnZp(itn) = int(ZAP)
						tnNyld(itn) = NP6
						read(100,12)(tnNBTp(itn,N),tnINTrp(itn,N),N=1,NR6)
						read(100,13)(tnEyld(itn,N),tnYld(itn,N),N=1,NP6)
						end if
					end if
					if (LAW.ne.0) then	! for calculating from JEFF-3.3
						read(100,11)c1,c2,l3,l4,NR6,NE6,MAT,MF,MT,NS
						read(100,12) (NBT, INTr, N=1, NR6)
						do i = 1, NE6
							read(100,11)c1,En,ND,NA,NW,NEP,MAT,MF,MT,NS
							READ(100,13) (Bvall, J1 = 1, NW)
						end do
					end if
		
				end do
				itnt = itn
			end if
		end if
	else 
		exit
	end if
end do

10 		format(66x,I4,I2,I3,I5)
11 		format(2e11.0,4I11,I4,I2,I3,I5)
12 		format(6i11)
13      FORMAT(6E11.0)
		
		close(unit=100)
		end subroutine gtYMf6Mt5
!============================================================
	
	!! To find the presence of a particular MT in MF = 3
		
    subroutine FindMT(MTfind,iflMTpr)
		integer,dimension(1000):: MFs,MTs
		call file1 (nfiles,MFs,MTs)
		
		iflMTpr = 0
		do i = 1, nfiles
			if (MFs(i)==3 .and. MTs(i)==MTfind) then
				iflMTpr = 1
				exit
			end if
		end do
	end subroutine FindMT
!========================================
	!! To collect the directory of reactions evaluated and given in
	!! the directory of File 1 in tape01.
		
    subroutine file1(nfiles,MFs,MTs)
		integer,dimension(1000):: MFs,MTs
		! max. NXC = 350 (according to ENDF-102 formats), 
		! but deviates for Mn55 ENDF/B-VII.1, so increased to 1000
		open(unit=100,file='tape01')
		read(100,*)
		read(100,1)ZA,AWR,LRP,LFI,NLIB,NMOD,MAT,MF,MT
		read(100,1)ELIS,STA,LIS,LISO,num,NFOR,MAT,MF,MT
		read(100,1)AWI,EMAX,LREL,num,NSUB,NVER,MAT,MF,MT
		read(100,1)TEMP,c2,LDRV,num,NWD,NXC,MAT,MF,MT
		do i = 1,NWD 
		read(100,*)
		end do
		do i = 1, NXC
		read(100,1)blnk,blnk,MFs(i),MTs(i),NCn,MODn,MAT,MF,MT
		end do
		nfiles = NXC
1 		format(2e11.0,4I11,I4,I2,I3)
		
		close(unit=100)

	end subroutine file1
!==================================================================

	!! Calculation of total (n,p), (n,d), (n,tr), (n,3He) and (n,a)
	!! cross sections from MT = 600 to 849.
		
    subroutine adddiscnth(iflMTpr,MTnth,E,NP,signth)
		real,dimension(NP):: E, signth
		real,allocatable,dimension(:) :: E1, sig1
		iflMTpr = 0
		if (MTnth==103) then 
			MTi = 600
			MTimax = 649
		end if
		if (MTnth==104) then 
			MTi = 650
			MTimax = 699
		end if
		if (MTnth==105) then 
			MTi = 700
			MTimax = 749
		end if
		if (MTnth==106) then 
			MTi = 750
			MTimax = 799
		end if
		if (MTnth==107) then 
			MTi = 800
			MTimax = 849
		end if

		signth = 0
		open (unit = 500, file = "tape02")

		do
		NP1 = 0
		iflp = 0
		if (MTi<=MTimax) then
 		read(500,1) MAT, MF, MT, NS
 		if ( MT==0 .and. MF.ne.0) then	! .and. NS==99999
  		read(500,11)ZAv,AWRv,L0,L1,NKv,L2,MAT,MF,MT,NS
 		end if
 		if (MAT.ne.-1) then
 		if (MF.eq.3) then
 		if (MT.eq.MTi .or. MT.eq.MTimax) then
			if (MT==MTimax) MTi = MT
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
 		else
 		exit
 		end if
		if (iflp==1) then
		do i = 1, NP
			if (E(i)>=E1(1)) then
				istart = i
				exit
			end if
		end do
		do i = istart, NP
 			do j = 1, NP1
 				if (E(i) == E1(j)) then
				   signth(i) = signth(i) + sig1(j)
				   exit
 				end if
 				if (E(i)>E1(j) .and. E(i)<E1(j+1)) then
		zf = crstd(E(i),E1(j),E1(j+1),sig1(j),sig1(j+1))
 		signth(i) = signth(i) + zf
				exit
 				end if
			end do
		end do
		deallocate(E1,sig1)
		end if

		else
		exit
		end if
		
		end do
		
		close(unit=500)
		
1 		format(66x,I4,I2,I3,I5)
2		format(2e11.0,11x,3i11)
3		format(6e11.0)
11 		format(2e11.0,4I11,I4,I2,I3,I5) 
		
	end subroutine adddiscnth
!==================================================================

	!! The weighting spectrum for multigrouping of cross sections
		
        function spectrum(En,L)
 		real En
 		select case(L)
 		case (1)
 		fi = En * exp(-En/0.0253)		
 		case (2)
 		fi = 1/En
 		case (3)
 		fi = sqrt(En) * exp(-En/(1.415e+6))
 		case (4)
 		fi = 1
 		case default
 		fi = 0
 		end select
 		spectrum = fi
 		end function spectrum
		
        function srchintrp3 (E,fi,Ngl,En)
		real, dimension(Ngl) :: E, fi
		do i = 1, Ngl
			if (E(i) == En) then
				cr1 = fi(i)
				exit
			end if
			if (i<=Ngl-1) then
			if (E(i) < En .and. En < E(i+1)) then
				x = En
				x1 = E(i)
				x2 = E(i+1)
				y1 = fi(i)
				y2 = fi(i+1)
				cr1 = crstd(x,x1,x2,y1,y2)
				exit
			end if
			end if
		end do
		srchintrp3 = cr1
		end function srchintrp3
		
		real function crstd(x,x1,x2,y1,y2)
		real x,y,x1,x2,y1,y2
		y = y1+((y2-y1)*(x-x1)/(x2-x1))
		crstd = y
		end function crstd
!=======================================================================

		!! Energy Group Structures
		
        subroutine engrp1 (Eg,Ngl)
		real, dimension (:):: Eg(176)
		Ngl = 176
		Eg = (/1.000E-05,1.000E-01,4.140E-01,5.320E-01,6.830E-01, &
	8.760E-01,1.130E+00,1.450E+00,1.860E+00,2.380E+00, &	
	3.060E+00,3.930E+00,5.040E+00,6.480E+00,8.320E+00, &	
	1.070E+01,1.370E+01,1.760E+01,2.260E+01,2.900E+01, &	
	3.730E+01,4.790E+01,6.140E+01,7.890E+01,1.010E+02, &	
	1.300E+02,1.670E+02,2.140E+02,2.750E+02,3.540E+02, &	
	4.540E+02,5.830E+02,7.490E+02,9.610E+02,1.230E+03, &	
	1.580E+03,2.030E+03,2.250E+03,2.490E+03,2.610E+03, &	
	2.750E+03,3.040E+03,3.350E+03,3.710E+03,4.310E+03, &	
	5.530E+03,7.100E+03,9.120E+03,1.060E+04,1.170E+04, &	
	1.500E+04,1.930E+04,2.190E+04,2.360E+04,2.420E+04, &	
	2.480E+04,2.610E+04,2.700E+04,2.850E+04,3.180E+04, &	
	3.430E+04,4.090E+04,4.630E+04,5.250E+04,5.660E+04, &	
	6.740E+04,7.200E+04,7.950E+04,8.250E+04,8.650E+04, &	
	9.800E+04,1.110E+05,1.170E+05,1.230E+05,1.290E+05, &	
	1.360E+05,1.430E+05,1.500E+05,1.580E+05,1.660E+05, &	
	1.740E+05,1.830E+05,1.930E+05,2.020E+05,2.130E+05, &	
	2.240E+05,2.350E+05,2.470E+05,2.730E+05,2.870E+05, &	
	2.950E+05,2.970E+05,2.990E+05,3.020E+05,3.340E+05, &	
	3.690E+05,3.880E+05,4.080E+05,4.500E+05,4.980E+05, &	
	5.230E+05,5.500E+05,5.780E+05,6.080E+05,6.390E+05, &	
	6.720E+05,7.070E+05,7.430E+05,7.810E+05,8.210E+05, &	
	8.630E+05,9.070E+05,9.620E+05,1.000E+06,1.110E+06, &	
	1.160E+06,1.220E+06,1.290E+06,1.350E+06,1.420E+06, &	
	1.500E+06,1.570E+06,1.650E+06,1.740E+06,1.830E+06, &	
	1.920E+06,2.020E+06,2.120E+06,2.230E+06,2.310E+06, &	
	2.350E+06,2.370E+06,2.390E+06,2.470E+06,2.590E+06, &	
	2.730E+06,2.870E+06,3.010E+06,3.170E+06,3.330E+06, &	
	3.680E+06,4.070E+06,4.490E+06,4.720E+06,4.970E+06, &	
	5.220E+06,5.490E+06,5.770E+06,6.070E+06,6.380E+06, &	
	6.590E+06,6.700E+06,7.050E+06,7.410E+06,7.790E+06, &	
	8.190E+06,8.610E+06,9.050E+06,9.510E+06,1.000E+07, &	
	1.050E+07,1.110E+07,1.160E+07,1.220E+07,1.250E+07, &	
	1.280E+07,1.350E+07,1.380E+07,1.420E+07,1.460E+07, &	
	1.490E+07,1.570E+07,1.650E+07,1.690E+07,1.730E+07, &	
	1.960E+07/)	
			end subroutine engrp1

        subroutine engrp2(Eg,Ngl)
		real, dimension (:):: Eg(239)
		Ngl = 239
		Eg=(/0.1000E-04,0.1000E-03,0.5000E-03,0.7500E-03,0.1000E-02, &
	0.1200E-02,0.1500E-02,0.2000E-02,0.2500E-02,0.3000E-02, & 	
	0.4000E-02,0.5000E-02,0.7500E-02,0.1000E-01,0.2530E-01, & 	
	0.3000E-01,0.4000E-01,0.5000E-01,0.6000E-01,0.7000E-01, & 	
	0.8000E-01,0.9000E-01,0.1000E+00,0.1250E+00,0.1500E+00, & 	
	0.1750E+00,0.2000E+00,0.2250E+00,0.2500E+00,0.2750E+00, & 	
	0.3000E+00,0.3250E+00,0.3500E+00,0.3750E+00,0.4000E+00, & 	
	0.4500E+00,0.5000E+00,0.5500E+00,0.6000E+00,0.6250E+00, & 	
	0.6500E+00,0.7000E+00,0.7500E+00,0.8000E+00,0.8500E+00, & 	
	0.9000E+00,0.9250E+00,0.9500E+00,0.9750E+00,0.1000E+01, & 	
	0.1010E+01,0.1020E+01,0.1030E+01,0.1040E+01,0.1050E+01, & 	
	0.1060E+01,0.1070E+01,0.1080E+01,0.1090E+01,0.1100E+01, & 	
	0.1110E+01,0.1120E+01,0.1130E+01,0.1140E+01,0.1150E+01, & 	
	0.1175E+01,0.1200E+01,0.1225E+01,0.1250E+01,0.1300E+01, & 	
	0.1350E+01,0.1400E+01,0.1450E+01,0.1500E+01,0.1590E+01, & 	
	0.1680E+01,0.1770E+01,0.1860E+01,0.1940E+01,0.2000E+01, & 	
	0.2120E+01,0.2210E+01,0.2300E+01,0.2380E+01,0.2470E+01, & 	
	0.2570E+01,0.2670E+01,0.2770E+01,0.2870E+01,0.2970E+01, & 	
	0.3000E+01,0.3050E+01,0.3150E+01,0.3500E+01,0.3730E+01, & 	
	0.4000E+01,0.4750E+01,0.5000E+01,0.5400E+01,0.6000E+01, & 	
	0.6250E+01,0.6500E+01,0.6750E+01,0.7000E+01,0.7150E+01, & 	
	0.8100E+01,0.9100E+01,0.1000E+02,0.1150E+02,0.1190E+02, & 	
	0.1290E+02,0.1375E+02,0.1440E+02,0.1510E+02,0.1600E+02, & 	
	0.1700E+02,0.1850E+02,0.1900E+02,0.2000E+02,0.2100E+02, & 	
	0.2250E+02,0.2500E+02,0.2750E+02,0.3000E+02,0.3125E+02, & 	
	0.3175E+02,0.3325E+02,0.3375E+02,0.3460E+02,0.3550E+02, & 	
	0.3700E+02,0.3800E+02,0.3910E+02,0.3960E+02,0.4100E+02, & 	
	0.4240E+02,0.4400E+02,0.4520E+02,0.4700E+02,0.4830E+02, & 	
	0.4920E+02,0.5060E+02,0.5200E+02,0.5340E+02,0.5900E+02, & 	
	0.6100E+02,0.6500E+02,0.6750E+02,0.7200E+02,0.7600E+02, & 	
	0.8000E+02,0.8200E+02,0.9000E+02,0.1000E+03,0.1080E+03, & 	
	0.1150E+03,0.1190E+03,0.1220E+03,0.1860E+03,0.1925E+03, & 	
	0.2075E+03,0.2100E+03,0.2400E+03,0.2850E+03,0.3050E+03, & 	
	0.5500E+03,0.6700E+03,0.6830E+03,0.9500E+03,0.1150E+04, & 	
	0.1500E+04,0.1550E+04,0.1800E+04,0.2200E+04,0.2290E+04, & 	
	0.2580E+04,0.3000E+04,0.3740E+04,0.3900E+04,0.6000E+04, & 	
	0.8030E+04,0.9500E+04,0.1300E+05,0.1700E+05,0.2500E+05, & 	
	0.3000E+05,0.4500E+05,0.5000E+05,0.5200E+05,0.6000E+05, & 	
	0.7300E+05,0.7500E+05,0.8200E+05,0.8500E+05,0.1000E+06, & 	
	0.1283E+06,0.1500E+06,0.2000E+06,0.2700E+06,0.3300E+06, & 	
	0.4000E+06,0.4200E+06,0.4400E+06,0.4700E+06,0.4995E+06, & 	
	0.5500E+06,0.5730E+06,0.6000E+06,0.6700E+06,0.6790E+06, & 	
	0.7500E+06,0.8200E+06,0.8611E+06,0.8750E+06,0.9000E+06, & 	
	0.9200E+06,0.1010E+07,0.1100E+07,0.1200E+07,0.1250E+07, & 	
	0.1317E+07,0.1356E+07,0.1400E+07,0.1500E+07,0.1850E+07, & 	
	0.2354E+07,0.2479E+07,0.3000E+07,0.4304E+07,0.4800E+07, & 	
	0.6434E+07,0.8187E+07,0.1000E+08,0.1284E+08,0.1384E+08, & 	
	0.1455E+08,0.1568E+08,0.1733E+08,0.2000E+08/) 	
	
		end subroutine engrp2

		subroutine engrp3(Eg,Ngl)
		real,dimension (:):: Eg(199)
		Ngl = 199
		Eg=(/1.0000E-05,5.0000E-04,2.0000E-03,5.0000E-03,1.0000E-02, &
		1.4500E-02,2.1000E-02,3.0000E-02,4.0000E-02,5.0000E-02, &
		7.0000E-02,1.0000E-01,1.2500E-01,1.5000E-01,1.8400E-01, &
		2.2500E-01,2.7500E-01,3.2500E-01,3.6680E-01,4.1399E-01, &
		5.0000E-01,5.3158E-01,6.2506E-01,6.8256E-01,8.0000E-01, &
		8.7643E-01,1.0000E+00,1.0400E+00,1.0800E+00,1.1253E+00, &
		1.3000E+00,1.4450E+00,1.8554E+00,2.3824E+00,3.0590E+00, &
		3.9279E+00,5.0435E+00,6.4760E+00,8.3153E+00,1.0677E+01, &
		1.3710E+01,1.7604E+01,2.2603E+01,2.9023E+01,3.7266E+01, &
		4.7851E+01,6.1442E+01,7.8893E+01,1.0130E+02,1.3007E+02, &
		1.6702E+02,2.1445E+02,2.7536E+02,3.5357E+02,4.5400E+02, &
		5.8295E+02,7.4852E+02,9.6112E+02,1.2341E+03,1.5846E+03, &
		2.0347E+03,2.2487E+03,2.4852E+03,2.6126E+03,2.7465E+03, &
		3.0354E+03,3.3546E+03,3.7074E+03,4.3074E+03,5.5308E+03, &
		7.1017E+03,9.1188E+03,1.0595E+04,1.1709E+04,1.5034E+04, &
		1.9305E+04,2.1875E+04,2.3579E+04,2.4176E+04,2.4788E+04, &
		2.6058E+04,2.7000E+04,2.8501E+04,3.1828E+04,3.4307E+04, &
		4.0868E+04,4.6309E+04,5.2475E+04,5.6562E+04,6.7379E+04, &
		7.1998E+04,7.9499E+04,8.2503E+04,8.6517E+04,9.8037E+04, &
		1.1109E+05,1.1679E+05,1.2277E+05,1.2907E+05,1.3569E+05, &
		1.4264E+05,1.4996E+05,1.5764E+05,1.6573E+05,1.7422E+05, &
		1.8316E+05,1.9255E+05,2.0242E+05,2.1280E+05,2.2371E+05, &
		2.3518E+05,2.4724E+05,2.7324E+05,2.8725E+05,2.9452E+05, &
		2.9721E+05,2.9849E+05,3.0197E+05,3.3373E+05,3.6883E+05, &
		3.8774E+05,4.0762E+05,4.5049E+05,4.9787E+05,5.2340E+05, &
		5.5023E+05,5.7844E+05,6.0810E+05,6.3928E+05,6.7206E+05, &
		7.0651E+05,7.4274E+05,7.8082E+05,8.2085E+05,8.6294E+05, &
		9.0718E+05,9.6164E+05,1.0026E+06,1.1080E+06,1.1648E+06, &
		1.2246E+06,1.2874E+06,1.3534E+06,1.4227E+06,1.4957E+06, &
		1.5724E+06,1.6530E+06,1.7377E+06,1.8268E+06,1.9205E+06, &
		2.0190E+06,2.1225E+06,2.2313E+06,2.3069E+06,2.3457E+06, &
		2.3653E+06,2.3852E+06,2.4660E+06,2.5924E+06,2.7253E+06, &
		2.8651E+06,3.0119E+06,3.1664E+06,3.3287E+06,3.6788E+06, &
		4.0657E+06,4.4933E+06,4.7237E+06,4.9659E+06,5.2205E+06, &
		5.4881E+06,5.7695E+06,6.0653E+06,6.3763E+06,6.5924E+06, &
		6.7032E+06,7.0469E+06,7.4082E+06,7.7880E+06,8.1873E+06, &
		8.6071E+06,9.0484E+06,9.5123E+06,1.0000E+07,1.0513E+07, &
		1.1052E+07,1.1618E+07,1.2214E+07,1.2523E+07,1.2840E+07, &
		1.3499E+07,1.3840E+07,1.4191E+07,1.4550E+07,1.4918E+07, &
		1.5683E+07,1.6487E+07,1.6905E+07,1.9640E+07/)

		end subroutine engrp3
		
		subroutine engrp4(Eg,Ngl)
		real, dimension (:):: Eg(34)
		Ngl = 34
		Eg=(/1.000e-05,1.000e-01,5.400e-01,4.000e+00,8.315e+00, &
		1.371e+01,2.260e+01,4.017e+01,6.790e+01,9.166e+01,1.486e+02, &
		3.043e+02,4.540e+02,7.485e+02,1.230e+03,2.030e+03,3.355e+03, &
		5.531e+03,9.119e+03,1.503e+04,2.479e+04,4.087e+04,6.738e+04, &
		1.111e+05,1.832e+05,3.020e+05,4.979e+05,8.209e+05,1.353e+06, &
		2.231e+06,3.679e+06,6.065e+06,1.000e+07,1.964e+07/)
	
		end subroutine engrp4
		
		subroutine engrp5 (Eg,Ngl)
		real, dimension (:):: Eg(27)
		Ngl = 27
		Eg = (/1.0E-05,0.0253E+00,0.4642E+00,1.0E+00,2.1544E+00, &
		4.6416E+00,1.0E+01,2.1544E+01,4.6416E+01,1.0E+02,2.1544E+02, &
		4.6416E+02,1.0E+03,2.1544E+03,4.6416E+03,1.0E+04,2.1544E+04, &
		4.6416E+04,1.0E+05,2.0E+05,0.4E+06,0.8E+06,1.4E+06,2.5E+06,  &
		4.0E+06,6.5E+06,1.05E+07/)
		
		end subroutine engrp5
		
		subroutine engrp6 (Eg,Ngl)
		real, dimension (:):: Eg(101) ! DLC-2
		Ngl = 101
		Eg = (/1.0000E-05,4.1399E-01,5.3158E-01,6.8256E-01,8.7642E-01, &
		1.1254E+00,1.4450E+00,1.8554E+00,2.3824E+00,3.0590E+00, &
		3.9279E+00,5.0435E+00,6.4760E+00,8.3153E+00,1.0677E+01, &
		1.3710E+01,1.7603E+01,2.2603E+01,2.9023E+01,3.7267E+01, &
		4.7851E+01,6.1442E+01,7.8893E+01,1.0130E+02,1.3007E+02, &
		1.6702E+02,2.1445E+02,2.7536E+02,3.5358E+02,4.5400E+02, &
		5.8295E+02,7.4852E+02,9.6112E+02,1.2341E+03,1.5846E+03, &
		2.0347E+03,2.6126E+03,3.3546E+03,4.3074E+03,5.5308E+03, &
		7.1017E+03,9.1188E+03,1.1709E+04,1.5034E+04,1.9305E+04, &
		2.4788E+04,3.1828E+04,4.0868E+04,5.2475E+04,6.7379E+04, &
		8.6517E+04,1.1109E+05,1.2277E+05,1.3569E+05,1.4996E+05, &
		1.6573E+05,1.8316E+05,2.0242E+05,2.2371E+05,2.4724E+05, &
		2.7324E+05,3.0197E+05,3.3373E+05,3.6883E+05,4.0762E+05, &
		4.5049E+05,4.9787E+05,5.5023E+05,6.0810E+05,6.7206E+05, &
		7.4274E+05,8.2085E+05,9.0718E+05,1.0026E+06,1.1080E+06, &
		1.2246E+06,1.3534E+06,1.4957E+06,1.6530E+06,1.8268E+06, &
		2.0190E+06,2.2313E+06,2.4660E+06,2.7253E+06,3.0119E+06, &
		3.3287E+06,3.6788E+06,4.0657E+06,4.4933E+06,4.9659E+06, &
		5.4881E+06,6.0653E+06,6.7032E+06,7.4082E+06,8.1873E+06, &
		9.0484E+06,1.0000E+07,1.1052E+07,1.2214E+07,1.3500E+07, &
		1.5000E+07/)
		
		end subroutine engrp6
		
		subroutine engrp7 (Eg,Ngl)	! EXTENDED FROM 198 GROUP STRUCTURE
		real, dimension (:):: Eg(229) ! BY SELF UP TO 200 MeV
		Ngl = 229
		Eg=(/1.0000E-05,5.0000E-04,2.0000E-03,5.0000E-03,1.0000E-02, &
		1.4500E-02,2.1000E-02,3.0000E-02,4.0000E-02,5.0000E-02, &
		7.0000E-02,1.0000E-01,1.2500E-01,1.5000E-01,1.8400E-01, &
		2.2500E-01,2.7500E-01,3.2500E-01,3.6680E-01,4.1399E-01, &
		5.0000E-01,5.3158E-01,6.2506E-01,6.8256E-01,8.0000E-01, &
		8.7643E-01,1.0000E+00,1.0400E+00,1.0800E+00,1.1253E+00, &
		1.3000E+00,1.4450E+00,1.8554E+00,2.3824E+00,3.0590E+00, &
		3.9279E+00,5.0435E+00,6.4760E+00,8.3153E+00,1.0677E+01, &
		1.3710E+01,1.7604E+01,2.2603E+01,2.9023E+01,3.7266E+01, &
		4.7851E+01,6.1442E+01,7.8893E+01,1.0130E+02,1.3007E+02, &
		1.6702E+02,2.1445E+02,2.7536E+02,3.5357E+02,4.5400E+02, &
		5.8295E+02,7.4852E+02,9.6112E+02,1.2341E+03,1.5846E+03, &
		2.0347E+03,2.2487E+03,2.4852E+03,2.6126E+03,2.7465E+03, &
		3.0354E+03,3.3546E+03,3.7074E+03,4.3074E+03,5.5308E+03, &
		7.1017E+03,9.1188E+03,1.0595E+04,1.1709E+04,1.5034E+04, &
		1.9305E+04,2.1875E+04,2.3579E+04,2.4176E+04,2.4788E+04, &
		2.6058E+04,2.7000E+04,2.8501E+04,3.1828E+04,3.4307E+04, &
		4.0868E+04,4.6309E+04,5.2475E+04,5.6562E+04,6.7379E+04, &
		7.1998E+04,7.9499E+04,8.2503E+04,8.6517E+04,9.8037E+04, &
		1.1109E+05,1.1679E+05,1.2277E+05,1.2907E+05,1.3569E+05, &
		1.4264E+05,1.4996E+05,1.5764E+05,1.6573E+05,1.7422E+05, &
		1.8316E+05,1.9255E+05,2.0242E+05,2.1280E+05,2.2371E+05, &
		2.3518E+05,2.4724E+05,2.7324E+05,2.8725E+05,2.9452E+05, &
		2.9721E+05,2.9849E+05,3.0197E+05,3.3373E+05,3.6883E+05, &
		3.8774E+05,4.0762E+05,4.5049E+05,4.9787E+05,5.2340E+05, &
		5.5023E+05,5.7844E+05,6.0810E+05,6.3928E+05,6.7206E+05, &
		7.0651E+05,7.4274E+05,7.8082E+05,8.2085E+05,8.6294E+05, &
		9.0718E+05,9.6164E+05,1.0026E+06,1.1080E+06,1.1648E+06, &
		1.2246E+06,1.2874E+06,1.3534E+06,1.4227E+06,1.4957E+06, &
		1.5724E+06,1.6530E+06,1.7377E+06,1.8268E+06,1.9205E+06, &
		2.0190E+06,2.1225E+06,2.2313E+06,2.3069E+06,2.3457E+06, &
		2.3653E+06,2.3852E+06,2.4660E+06,2.5924E+06,2.7253E+06, &
		2.8651E+06,3.0119E+06,3.1664E+06,3.3287E+06,3.6788E+06, &
		4.0657E+06,4.4933E+06,4.7237E+06,4.9659E+06,5.2205E+06, &
		5.4881E+06,5.7695E+06,6.0653E+06,6.3763E+06,6.5924E+06, &
		6.7032E+06,7.0469E+06,7.4082E+06,7.7880E+06,8.1873E+06, &
		8.6071E+06,9.0484E+06,9.5123E+06,1.0000E+07,1.0513E+07, &
		1.1052E+07,1.1618E+07,1.2214E+07,1.2523E+07,1.2840E+07, &
		1.3499E+07,1.3840E+07,1.4191E+07,1.4550E+07,1.4918E+07, &
		1.5683E+07,1.6487E+07,1.6905E+07,1.9640E+07, &
		2.10E+07,2.40E+07,2.60E+07,2.80E+07,3.00E+07,3.20E+07, &
		3.50E+07,3.80E+07,4.00E+07,4.50E+07,5.00E+07,5.50E+07, &
		6.00E+07,7.00E+07,8.00E+07,9.00E+07,1.00E+08,1.10E+08, &
		1.20E+08,1.30E+08,1.40E+08,1.50E+08,1.60E+08,1.70E+08, &
		1.80E+08,1.90E+08,1.95E+08,1.98E+08,1.99E+08,2.00E+08/)
		end subroutine engrp7
		
		subroutine engrp8 (Eg,Ngl)	! EXTENDED FROM 198 GROUP STRUCTURE
		real, dimension (:):: Eg(229) ! BY SELF UP TO 150 MeV
		Ngl = 229
		Eg=(/1.0000E-05,5.0000E-04,2.0000E-03,5.0000E-03,1.0000E-02, &
		1.4500E-02,2.1000E-02,3.0000E-02,4.0000E-02,5.0000E-02, &
		7.0000E-02,1.0000E-01,1.2500E-01,1.5000E-01,1.8400E-01, &
		2.2500E-01,2.7500E-01,3.2500E-01,3.6680E-01,4.1399E-01, &
		5.0000E-01,5.3158E-01,6.2506E-01,6.8256E-01,8.0000E-01, &
		8.7643E-01,1.0000E+00,1.0400E+00,1.0800E+00,1.1253E+00, &
		1.3000E+00,1.4450E+00,1.8554E+00,2.3824E+00,3.0590E+00, &
		3.9279E+00,5.0435E+00,6.4760E+00,8.3153E+00,1.0677E+01, &
		1.3710E+01,1.7604E+01,2.2603E+01,2.9023E+01,3.7266E+01, &
		4.7851E+01,6.1442E+01,7.8893E+01,1.0130E+02,1.3007E+02, &
		1.6702E+02,2.1445E+02,2.7536E+02,3.5357E+02,4.5400E+02, &
		5.8295E+02,7.4852E+02,9.6112E+02,1.2341E+03,1.5846E+03, &
		2.0347E+03,2.2487E+03,2.4852E+03,2.6126E+03,2.7465E+03, &
		3.0354E+03,3.3546E+03,3.7074E+03,4.3074E+03,5.5308E+03, &
		7.1017E+03,9.1188E+03,1.0595E+04,1.1709E+04,1.5034E+04, &
		1.9305E+04,2.1875E+04,2.3579E+04,2.4176E+04,2.4788E+04, &
		2.6058E+04,2.7000E+04,2.8501E+04,3.1828E+04,3.4307E+04, &
		4.0868E+04,4.6309E+04,5.2475E+04,5.6562E+04,6.7379E+04, &
		7.1998E+04,7.9499E+04,8.2503E+04,8.6517E+04,9.8037E+04, &
		1.1109E+05,1.1679E+05,1.2277E+05,1.2907E+05,1.3569E+05, &
		1.4264E+05,1.4996E+05,1.5764E+05,1.6573E+05,1.7422E+05, &
		1.8316E+05,1.9255E+05,2.0242E+05,2.1280E+05,2.2371E+05, &
		2.3518E+05,2.4724E+05,2.7324E+05,2.8725E+05,2.9452E+05, &
		2.9721E+05,2.9849E+05,3.0197E+05,3.3373E+05,3.6883E+05, &
		3.8774E+05,4.0762E+05,4.5049E+05,4.9787E+05,5.2340E+05, &
		5.5023E+05,5.7844E+05,6.0810E+05,6.3928E+05,6.7206E+05, &
		7.0651E+05,7.4274E+05,7.8082E+05,8.2085E+05,8.6294E+05, &
		9.0718E+05,9.6164E+05,1.0026E+06,1.1080E+06,1.1648E+06, &
		1.2246E+06,1.2874E+06,1.3534E+06,1.4227E+06,1.4957E+06, &
		1.5724E+06,1.6530E+06,1.7377E+06,1.8268E+06,1.9205E+06, &
		2.0190E+06,2.1225E+06,2.2313E+06,2.3069E+06,2.3457E+06, &
		2.3653E+06,2.3852E+06,2.4660E+06,2.5924E+06,2.7253E+06, &
		2.8651E+06,3.0119E+06,3.1664E+06,3.3287E+06,3.6788E+06, &
		4.0657E+06,4.4933E+06,4.7237E+06,4.9659E+06,5.2205E+06, &
		5.4881E+06,5.7695E+06,6.0653E+06,6.3763E+06,6.5924E+06, &
		6.7032E+06,7.0469E+06,7.4082E+06,7.7880E+06,8.1873E+06, &
		8.6071E+06,9.0484E+06,9.5123E+06,1.0000E+07,1.0513E+07, &
		1.1052E+07,1.1618E+07,1.2214E+07,1.2523E+07,1.2840E+07, &
		1.3499E+07,1.3840E+07,1.4191E+07,1.4550E+07,1.4918E+07, &
		1.5683E+07,1.6487E+07,1.6905E+07,1.9640E+07, &
		2.10E+07,2.40E+07,2.60E+07,2.80E+07,3.00E+07,3.20E+07, &
		3.50E+07,3.80E+07,4.00E+07,4.50E+07,5.00E+07,5.50E+07, &
		6.00E+07,6.50E+07,7.00E+07,7.50E+07,8.00E+07,8.50E+07, &
		9.00E+07,9.50E+07,1.00E+08,1.05E+08,1.10E+08,1.20E+08, &
		1.25E+08,1.30E+08,1.35E+08,1.40E+08,1.45E+08,1.50E+08/)
		end subroutine engrp8
		
		subroutine engrp9 (Eg,Ngl)	! 616 energy groups DEMO-HCPB-FW
		real, dimension (:):: Eg(617)
		Ngl = 617
		Eg=(/1.00E-05,1.05E-05,1.10E-05,1.15E-05,1.20E-05,1.26E-05, &
		1.32E-05,1.38E-05,1.45E-05,1.51E-05,1.58E-05,1.66E-05,1.74E-05, &
		1.82E-05,1.91E-05,2.00E-05,2.09E-05,2.19E-05,2.29E-05,2.40E-05, &
		2.51E-05,2.63E-05,2.75E-05,2.88E-05,3.02E-05,3.16E-05,3.31E-05, &
		3.47E-05,3.63E-05,3.80E-05,3.98E-05,4.17E-05,4.37E-05,4.57E-05, &
		4.79E-05,5.01E-05,5.25E-05,5.50E-05,5.75E-05,6.03E-05,6.31E-05, &
		6.61E-05,6.92E-05,7.24E-05,7.59E-05,7.94E-05,8.32E-05,8.71E-05, &
		9.12E-05,9.55E-05,1.00E-04,1.05E-04,1.10E-04,1.15E-04,1.20E-04, &
		1.26E-04,1.32E-04,1.38E-04,1.45E-04,1.51E-04,1.58E-04,1.66E-04, &
		1.74E-04,1.82E-04,1.91E-04,2.00E-04,2.09E-04,2.19E-04,2.29E-04, &
		2.40E-04,2.51E-04,2.63E-04,2.75E-04,2.88E-04,3.02E-04,3.16E-04, &
		3.31E-04,3.47E-04,3.63E-04,3.80E-04,3.98E-04,4.17E-04,4.37E-04, &
		4.57E-04,4.79E-04,5.01E-04,5.25E-04,5.50E-04,5.75E-04,6.03E-04, &
		6.31E-04,6.61E-04,6.92E-04,7.24E-04,7.59E-04,7.94E-04,8.32E-04, &
		8.71E-04,9.12E-04,9.55E-04,1.00E-03,1.05E-03,1.10E-03,1.15E-03, &
		1.20E-03,1.26E-03,1.32E-03,1.38E-03,1.45E-03,1.51E-03,1.58E-03, &
		1.66E-03,1.74E-03,1.82E-03,1.91E-03,2.00E-03,2.09E-03,2.19E-03, &
		2.29E-03,2.40E-03,2.51E-03,2.63E-03,2.75E-03,2.88E-03,3.02E-03, &
		3.16E-03,3.31E-03,3.47E-03,3.63E-03,3.80E-03,3.98E-03,4.17E-03, &
		4.37E-03,4.57E-03,4.79E-03,5.01E-03,5.25E-03,5.50E-03,5.75E-03, &
		6.03E-03,6.31E-03,6.61E-03,6.92E-03,7.24E-03,7.59E-03,7.94E-03, &
		8.32E-03,8.71E-03,9.12E-03,9.55E-03,1.00E-02,1.05E-02,1.10E-02, &
		1.15E-02,1.20E-02,1.26E-02,1.32E-02,1.38E-02,1.45E-02,1.51E-02, &
		1.58E-02,1.66E-02,1.74E-02,1.82E-02,1.91E-02,2.00E-02,2.09E-02, &
		2.19E-02,2.29E-02,2.40E-02,2.51E-02,2.63E-02,2.75E-02,2.88E-02, &
		3.02E-02,3.16E-02,3.31E-02,3.47E-02,3.63E-02,3.80E-02,3.98E-02, &
		4.17E-02,4.37E-02,4.57E-02,4.79E-02,5.01E-02,5.25E-02,5.50E-02, &
		5.75E-02,6.03E-02,6.31E-02,6.61E-02,6.92E-02,7.24E-02,7.59E-02, &
		7.94E-02,8.32E-02,8.71E-02,9.12E-02,9.55E-02,1.00E-01,1.05E-01, &
		1.10E-01,1.15E-01,1.20E-01,1.26E-01,1.32E-01,1.38E-01,1.45E-01, &
		1.51E-01,1.58E-01,1.66E-01,1.74E-01,1.82E-01,1.91E-01,2.00E-01, &
		2.09E-01,2.19E-01,2.29E-01,2.40E-01,2.51E-01,2.63E-01,2.75E-01, &
		2.88E-01,3.02E-01,3.16E-01,3.31E-01,3.47E-01,3.63E-01,3.80E-01, &
		3.98E-01,4.17E-01,4.37E-01,4.57E-01,4.79E-01,5.01E-01,5.25E-01, &
		5.50E-01,5.75E-01,6.03E-01,6.31E-01,6.61E-01,6.92E-01,7.24E-01, &
		7.59E-01,7.94E-01,8.32E-01,8.71E-01,9.12E-01,9.55E-01,1.00E+00, &
		1.05E+00,1.10E+00,1.15E+00,1.20E+00,1.26E+00,1.32E+00,1.38E+00, &
		1.45E+00,1.51E+00,1.58E+00,1.66E+00,1.74E+00,1.82E+00,1.91E+00, &
		2.00E+00,2.09E+00,2.19E+00,2.29E+00,2.40E+00,2.51E+00,2.63E+00, &
		2.75E+00,2.88E+00,3.02E+00,3.16E+00,3.31E+00,3.47E+00,3.63E+00, &
		3.80E+00,3.98E+00,4.17E+00,4.37E+00,4.57E+00,4.79E+00,5.01E+00, &
		5.25E+00,5.50E+00,5.75E+00,6.03E+00,6.31E+00,6.61E+00,6.92E+00, &
		7.24E+00,7.59E+00,7.94E+00,8.32E+00,8.71E+00,9.12E+00,9.55E+00, &
		1.00E+01,1.05E+01,1.10E+01,1.15E+01,1.20E+01,1.26E+01,1.32E+01, &
		1.38E+01,1.45E+01,1.51E+01,1.58E+01,1.66E+01,1.74E+01,1.82E+01, &
		1.91E+01,2.00E+01,2.09E+01,2.19E+01,2.29E+01,2.40E+01,2.51E+01, &
		2.63E+01,2.75E+01,2.88E+01,3.02E+01,3.16E+01,3.31E+01,3.47E+01, &
		3.63E+01,3.80E+01,3.98E+01,4.17E+01,4.37E+01,4.57E+01,4.79E+01, &
		5.01E+01,5.25E+01,5.50E+01,5.75E+01,6.03E+01,6.31E+01,6.61E+01, &
		6.92E+01,7.24E+01,7.59E+01,7.94E+01,8.32E+01,8.71E+01,9.12E+01, &
		9.55E+01,1.00E+02,1.05E+02,1.10E+02,1.15E+02,1.20E+02,1.26E+02, &
		1.32E+02,1.38E+02,1.45E+02,1.51E+02,1.58E+02,1.66E+02,1.74E+02, &
		1.82E+02,1.91E+02,2.00E+02,2.09E+02,2.19E+02,2.29E+02,2.40E+02, &
		2.51E+02,2.63E+02,2.75E+02,2.88E+02,3.02E+02,3.16E+02,3.31E+02, &
		3.47E+02,3.63E+02,3.80E+02,3.98E+02,4.17E+02,4.37E+02,4.57E+02, &
		4.79E+02,5.01E+02,5.25E+02,5.50E+02,5.75E+02,6.03E+02,6.31E+02, &
		6.61E+02,6.92E+02,7.24E+02,7.59E+02,7.94E+02,8.32E+02,8.71E+02, &
		9.12E+02,9.55E+02,1.00E+03,1.05E+03,1.10E+03,1.15E+03,1.20E+03, &
		1.26E+03,1.32E+03,1.38E+03,1.45E+03,1.51E+03,1.58E+03,1.66E+03, &
		1.74E+03,1.82E+03,1.91E+03,2.00E+03,2.09E+03,2.19E+03,2.29E+03, &
		2.40E+03,2.51E+03,2.63E+03,2.75E+03,2.88E+03,3.02E+03,3.16E+03, &
		3.31E+03,3.47E+03,3.63E+03,3.80E+03,3.98E+03,4.17E+03,4.37E+03, &
		4.57E+03,4.79E+03,5.01E+03,5.25E+03,5.50E+03,5.75E+03,6.03E+03, &
		6.31E+03,6.61E+03,6.92E+03,7.24E+03,7.59E+03,7.94E+03,8.32E+03, &
		8.71E+03,9.12E+03,9.55E+03,1.00E+04,1.05E+04,1.10E+04,1.15E+04, &
		1.20E+04,1.26E+04,1.32E+04,1.38E+04,1.45E+04,1.51E+04,1.58E+04, &
		1.66E+04,1.74E+04,1.82E+04,1.91E+04,2.00E+04,2.09E+04,2.19E+04, &
		2.29E+04,2.40E+04,2.51E+04,2.63E+04,2.75E+04,2.88E+04,3.02E+04, &
		3.16E+04,3.31E+04,3.47E+04,3.63E+04,3.80E+04,3.98E+04,4.17E+04, &
		4.37E+04,4.57E+04,4.79E+04,5.01E+04,5.25E+04,5.50E+04,5.75E+04, &
		6.03E+04,6.31E+04,6.61E+04,6.92E+04,7.24E+04,7.59E+04,7.94E+04, &
		8.32E+04,8.71E+04,9.12E+04,9.55E+04,1.00E+05,1.05E+05,1.10E+05, &
		1.15E+05,1.20E+05,1.26E+05,1.32E+05,1.38E+05,1.45E+05,1.51E+05, &
		1.58E+05,1.66E+05,1.74E+05,1.82E+05,1.91E+05,2.00E+05,2.09E+05, &
		2.19E+05,2.29E+05,2.40E+05,2.51E+05,2.63E+05,2.75E+05,2.88E+05, &
		3.02E+05,3.16E+05,3.31E+05,3.47E+05,3.63E+05,3.80E+05,3.98E+05, &
		4.17E+05,4.37E+05,4.57E+05,4.79E+05,5.01E+05,5.25E+05,5.50E+05, &
		5.75E+05,6.03E+05,6.31E+05,6.61E+05,6.92E+05,7.24E+05,7.59E+05, &
		7.94E+05,8.32E+05,8.71E+05,9.12E+05,9.55E+05,1.00E+06,1.05E+06, &
		1.10E+06,1.15E+06,1.20E+06,1.26E+06,1.32E+06,1.38E+06,1.45E+06, &
		1.51E+06,1.58E+06,1.66E+06,1.74E+06,1.82E+06,1.91E+06,2.00E+06, &
		2.09E+06,2.19E+06,2.29E+06,2.40E+06,2.51E+06,2.63E+06,2.75E+06, &
		2.88E+06,3.02E+06,3.16E+06,3.31E+06,3.47E+06,3.63E+06,3.80E+06, &
		3.98E+06,4.17E+06,4.37E+06,4.57E+06,4.79E+06,5.01E+06,5.25E+06, &
		5.50E+06,5.75E+06,6.03E+06,6.31E+06,6.61E+06,6.92E+06,7.24E+06, &
		7.59E+06,7.94E+06,8.32E+06,8.71E+06,9.12E+06,9.55E+06,1.00E+07, &
		1.05E+07,1.10E+07,1.15E+07,1.20E+07,1.26E+07,1.32E+07,1.38E+07, &
		1.45E+07,1.51E+07,1.58E+07,1.66E+07,1.74E+07,1.82E+07,1.91E+07, &
		2.00E+07,2.00E+07/)
		end subroutine engrp9
		
end module TransmU
