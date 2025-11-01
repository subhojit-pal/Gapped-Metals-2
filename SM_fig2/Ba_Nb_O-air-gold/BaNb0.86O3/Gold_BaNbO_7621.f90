PROGRAM Casimir

!skip definitions and go down to the program!

!some of these parameters are no longer in use		
DOUBLE PRECISION :: x1,x2,wp,x,G 	
DOUBLE PRECISION :: vikt,viktb,sum,sumb,sum3W,sum3IN
DOUBLE PRECISION :: dq,q,D1,D2,D3,D2b,D3b
INTEGER :: L,k,dloop,ii,jj,jjj

INTEGER, PARAMETER :: kmax=7000 !steps in q integration vary to test accuracy 
	
DOUBLE PRECISION :: bsum,bsum2,bsum3,qmax
DOUBLE PRECISION :: qsum,qsum2,qsum3,D,qweight,qbweight
DOUBLE PRECISION :: weight,bweight
DOUBLE PRECISION :: qe, pi 
DOUBLE PRECISION :: eps1,eps2,eps3,eps4
DOUBLE PRECISION :: gamma1,gamma2,gamma3,gamma4 
DOUBLE PRECISION :: D234,D241,D232,D134,D141,D132
DOUBLE PRECISION :: film1,film2,bdist 
DOUBLE PRECISION :: Kb,T,qmin,c,s,h,Ng
DOUBLE PRECISION :: dist,D21,D43,D12,D20,D10,Gdisp,Gdisp1,Gdisp2
DOUBLE PRECISION :: wp1,wp2,wp3,DTM12,DTM23,DTE12,DTE23	
DOUBLE PRECISION :: B1,B2, wIR,wUV,CIR,CUV

!Parameters for dielectric function of gold
DOUBLE PRECISION :: AualphaModal_1,AufreqModal_1,AualphaModal_2,AufreqModal_2,AualphaModal_3,AufreqModal_3
DOUBLE PRECISION :: AualphaModal_4,AufreqModal_4,AualphaModal_5,AufreqModal_5,AualphaModal_6,AufreqModal_6
DOUBLE PRECISION :: AualphaModal_7,AufreqModal_7,AualphaModal_8,AufreqModal_8,AualphaModal_9,AufreqModal_9
DOUBLE PRECISION :: AualphaModal_10,AufreqModal_10,AualphaModal_11,AufreqModal_11
DOUBLE PRECISION :: AualphaModal_12,AufreqModal_12 
DOUBLE PRECISION :: AualphaModal_13,AufreqModal_13,AualphaModal_14,AufreqModal_14,AualphaModal_15
DOUBLE PRECISION :: AufreqModal_15,AualphaModal_16,AufreqModal_16
DOUBLE PRECISION :: epsAu

!parameters for dielectric function of silica
DOUBLE PRECISION :: eps0,epsINF,Wlo,wto,epsII,GGG1, GGG2, GGG3,SiO2epsp3
DOUBLE PRECISION :: SiO2alphaModal_1,SiO2alphaModal_2,SiO2alphaModal_3,SiO2alphaModal_4
DOUBLE PRECISION :: SiO2freqModal_1,SiO2freqModal_2,SiO2freqModal_3,SiO2freqModal_4

!ignore these parameters or check which are not used and remove those
DOUBLE PRECISION :: e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11
DOUBLE PRECISION :: f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11
DOUBLE PRECISION :: g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11
DOUBLE PRECISION :: ee1,ee2,ee3,ee4,ee5,ee6,ee7,ee8,ee9
DOUBLE PRECISION :: ff1,ff2,ff3,ff4,ff5,ff6,ff7,ff8,ff9
DOUBLE PRECISION :: gg1,gg2,gg3,gg4,gg5,gg6,gg7,gg8,gg9
DOUBLE PRECISION :: W1,W2,W3,W4,W5,A1,A2,A3,A4,A5
DOUBLE PRECISION :: cc1,cc2,tau1,tau2,temp,omega
DOUBLE PRECISION :: ome1,ome2,ome3,ome4,ome5,ome6,ome7,ome8
DOUBLE PRECISION :: ome9,ome10,ome11,ome12,ome13,ome14,ome15,ome16
DOUBLE PRECISION :: c1,c2,c3,c4,c5,c6,c7,c8,c9
DOUBLE PRECISION :: c10,c11,c12,c13,c14,c15,c16
DOUBLE PRECISION :: amma1,amma2,amma3,amma4,amma5 
DOUBLE PRECISION :: amma6,amma7,amma8
DOUBLE PRECISION :: amma9,amma10,amma11,amma12 
DOUBLE PRECISION :: amma13,amma14,amma15,amma16
DOUBLE PRECISION :: epswaterES,phi,epsLaTE
DOUBLE PRECISION :: alphaModal_1,freqModal_1,alphaModal_2,freqModal_2,alphaModal_3,freqModal_3
DOUBLE PRECISION :: alphaModal_4,freqModal_4,alphaModal_5,freqModal_5,alphaModal_6,freqModal_6
DOUBLE PRECISION :: alphaModal_7,freqModal_7,alphaModal_8,freqModal_8,alphaModal_9,freqModal_9
DOUBLE PRECISION :: alphaModal_10,freqModal_10,alphaModal_11,freqModal_11,alphaModal_12,freqModal_12
DOUBLE PRECISION :: alphaModal_13,freqModal_13
 

!!! Please focus on the lines below!!
!!! You need to define pi as a real or double precision in your own program and so on
!!! similar with integers like dloop
 

pi=3.14159265359D0
h=1.054571800D-34  !as noted by the numbers this is what is usually called "hbar" 
Kb=1.38064852D-23  !Boltzmann constant
c=2.99792458D8     !velocity of light 
T=300.0D0          !temperature in K
qe=1.60217662D-19 !unit charge in SI units




OPEN(UNIT=42, FILE = "Gold_BaNbO_7621")
dist=0D-10    !thickness air gap


DO dloop=1,800

dist=dist+1.D-10  ! calculate at different air gap thicknesses in SI units



bsum=0.D0
bsum2=0.D0
bsum3=0.D0


 
OPEN(UNIT=34, FILE = "Gold_BaNbO_7621_wp_eps1") !GM
OPEN(UNIT=35, FILE = "Gold_BaNbO_7621_wp_eps2") !Gold
DO i=0,2000 !3998 !2000      		

wp=(2.D0*pi*(i)*Kb*T/h) !Matsubara frequencies

AualphaModal_1 =36852.7364D0
AufreqModal_1 =(qe/h)*0.05D0
AualphaModal_2 =11.2904D0
AufreqModal_2 =(qe/h)*0.0546D0
AualphaModal_3 =0.7657D0
AufreqModal_3 =(qe/h)*2.4785D0
AualphaModal_4 =2.5141D0
AufreqModal_4 =(qe/h)*4.2668D0
AualphaModal_5 =1.6004D0
AufreqModal_5 =(qe/h)*9.4797D0
AualphaModal_6 =0.9337D0
AufreqModal_6 =(qe/h)*21.4693D0
AualphaModal_7 =0.2958D0
AufreqModal_7 =(qe/h)*46.3635D0
AualphaModal_8 =0.0217D0
AufreqModal_8 =(qe/h)*105.3024D0
AualphaModal_9 =0.0008D0
AufreqModal_9 =(qe/h)*162.2798D0
AualphaModal_10 =0.0D0
AufreqModal_10 =(qe/h)*275.9329D0
AualphaModal_11 =0.0001D0
AufreqModal_11 =(qe/h)*504.1258D0
AualphaModal_12 =0.0D0
AufreqModal_12 =(qe/h)*980.2226D0
AualphaModal_13 =0.0D0
AufreqModal_13 =(qe/h)*1951.4302D0
AualphaModal_14 =0.0D0
AufreqModal_14 =(qe/h)*3902.8604D0
AualphaModal_15 =0.0D0
AufreqModal_15 =(qe/h)*7805.7208D0
AualphaModal_16 =0.0D0
AufreqModal_16 =(qe/h)*15611.4415D0

epsAu=1.D0+(AualphaModal_1/(1.D0+(wp/AufreqModal_1)**2.D0))+&
&(AualphaModal_2/(1.D0+(wp/AufreqModal_2)**2.D0))+&
&(AualphaModal_3/(1.D0+(wp/AufreqModal_3)**2.D0))+&
&(AualphaModal_4/(1.D0+(wp/AufreqModal_4)**2.D0))+&
&(AualphaModal_5/(1.D0+(wp/AufreqModal_5)**2.D0))+&
&(AualphaModal_6/(1.D0+(wp/AufreqModal_6)**2.D0))+&
&(AualphaModal_7/(1.D0+(wp/AufreqModal_7)**2.D0))+&
&(AualphaModal_8/(1.D0+(wp/AufreqModal_8)**2.D0))+&
&(AualphaModal_9/(1.D0+(wp/AufreqModal_9)**2.D0))+&
&(AualphaModal_10/(1.D0+(wp/AufreqModal_10)**2.D0))+&
&(AualphaModal_11/(1.D0+(wp/AufreqModal_11)**2.D0))+&
&(AualphaModal_12/(1.D0+(wp/AufreqModal_12)**2.D0))+&
&(AualphaModal_13/(1.D0+(wp/AufreqModal_13)**2.D0))+&
&(AualphaModal_14/(1.D0+(wp/AufreqModal_14)**2.D0))+&
&(AualphaModal_15/(1.D0+(wp/AufreqModal_15)**2.D0))+&
&(AualphaModal_16/(1.D0+(wp/AufreqModal_16)**2.D0))



eps2=epsAu !final parameterised dielectric function for Gold 




!Parameterised dielectric function for a specific gapped metal surface
!!Ba7Nb6O21_VNb, composition_7_6_21, kpoints_20_20_20

alphaModal_1 =0.6139D0
freqModal_1 =0.0215D0*(qe/h)
alphaModal_2 =4.5368D0
freqModal_2 =0.0438D0*(qe/h)
alphaModal_3 =7.3998D0
freqModal_3 =0.0872D0*(qe/h)
alphaModal_4 =80.6752D0
freqModal_4 =0.1967D0*(qe/h)
alphaModal_5 =0.0D0
freqModal_5 =0.2227D0*(qe/h)
alphaModal_6 =6.027D0
freqModal_6 =0.6448D0*(qe/h)
alphaModal_7 =0.6325D0
freqModal_7 =2.5296D0*(qe/h)
alphaModal_8 =1.2681D0
freqModal_8 =5.0772D0*(qe/h)
alphaModal_9 =1.7721D0
freqModal_9 =8.4405D0*(qe/h)
alphaModal_10 =0.5311D0
freqModal_10 =16.0637D0*(qe/h)
alphaModal_11 =0.3452D0
freqModal_11 =27.1476D0*(qe/h)
alphaModal_12 =0.0074D0
freqModal_12 =47.1511D0*(qe/h)
alphaModal_13 =0.0114D0
freqModal_13 =79.6918D0*(qe/h)


 	
epsLaTE=1.D0+(alphaModal_1/(1.D0+(wp/freqModal_1)**2.D0))+&
&(alphaModal_2/(1.D0+(wp/freqModal_2)**2.D0))+&
&(alphaModal_3/(1.D0+(wp/freqModal_3)**2.D0))+&
&(alphaModal_4/(1.D0+(wp/freqModal_4)**2.D0))+&
&(alphaModal_5/(1.D0+(wp/freqModal_5)**2.D0))+&
&(alphaModal_6/(1.D0+(wp/freqModal_6)**2.D0))+&
&(alphaModal_7/(1.D0+(wp/freqModal_7)**2.D0))+&
&(alphaModal_8/(1.D0+(wp/freqModal_8)**2.D0))+&
&(alphaModal_9/(1.D0+(wp/freqModal_9)**2.D0))+&
&(alphaModal_10/(1.D0+(wp/freqModal_10)**2.D0))+&
&(alphaModal_11/(1.D0+(wp/freqModal_11)**2.D0))+&
&(alphaModal_12/(1.D0+(wp/freqModal_12)**2.D0))+&
&(alphaModal_13/(1.D0+(wp/freqModal_13)**2.D0))


!gapped metal 
eps1=epsLaTE !dielectric function for gapped metal surface
eps4=eps1


eps3=1.D0 !dielectric function of air


! eps1-eps4(bdist=0 in this case)-eps3(dist)-eps2
!gapped Metal-gapped Metal (zero thickness so not included)  - air- SiO2 


qsum=0.D0
qsum2=0.D0
qsum3=0.D0



qmax=(25D0)/(dist)
qmin=0     

  
dq= (qmax-qmin)/DFLOAT(kmax)

 

DO k=1,kmax  !q integration from qmin to qmax using simpson's rule

q=qmin+k*dq

 

	

! weight in simpsons formula
	
  	
IF (k/2*2/=k) THEN 
  		
qweight=4.D0 
		
qbweight=0.D0
  	
END IF
  
  	
IF (k/2*2==k) THEN 
  		 
qweight=2.D0
		
IF ((k)/4*4/=k) THEN
		
qbweight=4.D0
		
END IF
		
IF ((k)/4*4==k) THEN
		 
qbweight=2.D0
		
END IF
	
   	
END IF


	
IF ((k==0) .OR. (k==kmax) ) THEN 
  		
qweight=1.D0 
		
qbweight=1.D0
  	
END IF





!eps1-eps3(dist)-eps2
!gapped Metal-gapped Metal  - vapor- SiO2  

gamma1=DSQRT(q*q+eps1*wp*wp/(c*c) )	!retarded wave vectors in media 1
gamma2=DSQRT(q*q+eps2*wp*wp/(c*c) )
gamma3=DSQRT(q*q+eps3*wp*wp/(c*c) )
gamma4=DSQRT(q*q+eps4*wp*wp/(c*c) )

D234=-(eps4*gamma3-eps3*gamma4)/(eps4*gamma3+eps3*gamma4) !the 2 indicated that it is TM refl coeff
D241=-(eps1*gamma4-eps4*gamma1)/(eps1*gamma4+eps4*gamma1)
D232=-(eps2*gamma3-eps3*gamma2)/(eps2*gamma3+eps3*gamma2)


D134=-(gamma3-gamma4)/(gamma3+gamma4)			!!the 1 indicated that it is TE refl coeff
D141=-(gamma4-gamma1)/(gamma4+gamma1)
D132=-(gamma3-gamma2)/(gamma3+gamma2)



film1=D234*D232         !product of TM reflection coefficient
film2=D134*D132		!product of TE reflection coefficient


!(see the expression for Lifshitz energy contributions from TM and TE gives below)

Gdisp1=1.D0-DEXP(-2.D0*dist*gamma3)*film1	!retardation TM thin films 
Gdisp2=1.D0-DEXP(-2.D0*dist*gamma3)*film2	!retardation TE thin films

s=(Kb*T/(2.D0*pi))*(LOG(Gdisp1)+LOG(Gdisp2) )  
!Expression for free energy/unit area except for constants and sum+integration 
 
!use the derivative of free energy, below, if you want to calculate pressure rather than free energy
!B1=DEXP(-2.D0*dist*gamma3)*film1   !not used for free energy. If you want to calculate force use it
!B2= DEXP(-2.D0*dist*gamma3)*film2 !not used for free energy. If you want to calculate force use it
!s=-(Kb*T/(pi))*gamma3*((B1/Gdisp1)+(B2/Gdisp2)) !pressure



IF (i==0) THEN
	s=s/2.D0	!original sum is from minus to plus infinity so when going over to sum
			!from zero to infinity the zero frequency term should be multiplied by 0.5
END IF
 

qsum=qsum+s*dq*q*qweight/3.D0           !q integration
!qsum2=qsum2+s*2.D0*dq*q*qbweight/3.D0 	!q integration with different steplength to test accuracy
					!check Simpsons rule


END DO		!end of q loop


!qsum3=(qsum+(qsum-qsum2)/15)   !a way to improve the accuracy often just using qsum is ok
qsum3=qsum


bsum=bsum+qsum3 		!summand in frequency summation
 

!WRITE(40,*) eps2,'  ',eps1

WRITE(34,*) wp,' , ',eps1
WRITE(35,*) wp,' , ',eps2

END DO		!End of i loop (Matsubara frequency summation)

close(34)
close(35)
 


WRITE(*,*) dist*1D9,'   ',bsum   !print on screen dist in nm and free energy in J/m^2
WRITE(42,*) dist*1D9,' ,  ',bsum !print to file dist in nm and free energy in J/m^2
 



 !Always test the dielectric function data to make sure they are generated correctly.
!save data for dieelectric function to data file and then plot often

 
 

 
END DO !dloop

close(42)




END PROGRAM Casimir

