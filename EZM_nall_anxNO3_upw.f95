!11/22/23: cleaned up to send to Konstantin Weber and Mick

PROGRAM EZM

IMPLICIT NONE

INTEGER,PARAMETER :: Hp=2000, Lp=1e7, & 
    ndays=1e5, &
	dy=1e5,dz=10, &  
	ny=Lp/dy, nz=Hp/dz, &
	Aij=dy*dz !area of each grid box
REAL*8,PARAMETER :: dt=0.05D0, &
	!ndays=0.05, & !for fractions of days
	H=2000D0, & !!!!!change these scales,too!!!!!
	euz=25D0, & !60-110 euphotic layer depth from Paulmier, 2009 (ML depth is deeper!no?!)
	mlz=20D0, & !40m mixed layer depth
	o2sat=0.2121D0, & !mol/m3 from calc_oxsat(25+273,35) in matlab. WOCE clim-avg surf T at 10S, E. Pac.
	Kg=3D-5, & !m/s
	Ws=10D0, & !m/day !3 from Anderson 2007
    Ws2=0D0, & !m/day !sinking speed for second pool (called "DOM") -- Ws2=0 for actual DOM
	kappazmin=1D-5, & !m2/s -spin up, then change to 1e-5, value for most of the deep ocean (higher at top and bottom)
	!
	!BACTERIA METABOLISMS: 
	!
	!Bo: Aerobic Heterotroph
	yd_bo=0.14D0, & !mol cells/mol Detritus - aerobic bacteria yield
	yo_bo=yd_bo/20*4/(1-yd_bo)*1D0, & !mol cells/mol O2 -where cells have 1mol N 
	enh4_bo=(1/yd_bo-1)*1D0, & !production of ammonia per mol cells produced
	!
	!Bdno3: Anaerobic Denitrifer Step 1: NO3 to NO2
	yd_bdno3=0.11D0, & !mol cells/mol Detritus for full denitr (NO3 to N2)	
	!yno3_bdno3=ynd/29.2*5/(1-ynd)*1D0, & !for full denitr (NO3 to N2)	
	yno3_bdno3=yd_bdno3/20*2/(1-yd_bdno3)*1D0, & !Check this!! (NO3 to NO2)	
	enh4_bdno3=(1D0/yd_bdno3-1)*1D0, & !ammonium excretion
	eno2_bdno3=1D0/yno3_bdno3*1D0, & !NO2 excretion
	!
	!Bdno2: Anaerobic Denitrifer Step 2 AND 3: NO2 to N2 
	yd_bdno2=0.11D0, & !mol cells/mol Detritus for full denitr (NO3 to N2)
	yno2_bdno2=yd_bdno2/20*3/(1-yd_bdno2)*1D0, &	
	enh4_bdno2=(1D0/yd_bdno2-1)*1D0, &
	!en2o_bdno2=1/yno2_bdno2*1D0, & !to N2o means divide by 2 here, but, we're accounting for total N
	!
	!Bana: Anammox -yield matches Strous 1998, not theory: using f=0.1- Check:
	!Strous 1998, balanced
    !ynh4_bana=1D0/223D0, &!1/(1+20/3/fauto)*1D0, & !=0.0102*1D0 !mol cells/mol NH4 
    ynh4_bana=1D0/153.7778D0, &!1/(1+20/3/fauto)*1D0, & !=0.0102*1D0 !mol cells/mol NH4 
	yno2_bana=1D0/215.5556D0, & !fauto/20*3/(1-fauto)*1D0, & !=0.0114*1D0, & !mol cells/mol oxygen
	en2_bana=326.6667D0, &!N2 produced (1/ynh4_bana+1/yno2_bana-1)*1D0, & 
	eno3_bana=41.6667D0, &!N2 produced (1/ynh4_bana+1/yno2_bana-1)*1D0, & 
	!
	!Bnh4: Ammonia oxidizer
	ynh4_bnh4=1D0/112D0, &!1/(1+20/6/fauto)*1D0, & !mol cells/mol NH4 !0.0187*1D0, & !mol cells/mol NH4 
	yo_bnh4=1D0/162D0, &!fauto/20*4/(1-fauto)*1D0, & !0.0140*1D0, & !mol cells/mol oxygen
	eno2_bnh4=(1D0/ynh4_bnh4-1)*1D0, & 
	!
	!Bno2: Nitrite oxidizer
	yno2_bno2=1D0/334D0, &!1/(1+20/2/fauto)*1D0, & !0.0043*1D0, & !mol cells/mol NO2
	yo_bno2=1D0/162D0, &!fauto/20*4/(1-fauto)*1D0, & !0.0094*1D0, & !mol cells/mol oxygen
	eno3_bno2=(1D0/yno2_bno2-1)*1D0, & 
	!
	!P: phytoplankton
	!eo_p=8, & !number for O2 output as average of all N 
	umaxp2=3D0, & !max growth rate for p
	umaxp3=0.515D0, & !max growth rate for p
	kI=10D0, &
	!kBio=0.005D3, & !biomass attenuation coefficient
	kBio=0.04D3, & !biomass attenuation coefficient
	Iinmax=1000D0, & !then calc Iin from Iinmax -- use euz to calculate I at each depth:
	!New: (Iinmax from nitr was 1400)
    chl2cmax=0.2D0, & !mg Chl/mmol C from Dutk 2015 
    chl2cmin=0D0, & !.02D0, & !min chl2c 
    phimax=40D0, & !mmol C/mol photons/Ein  -quantum yield (mol C/mol photons)
    a_chlp2=0.01D0, & !diatom
    a_chlp3=0.04D0, & !ll pro m2/mg chl a -absorption parameter (in paper as m2/mg chla) avg over all wavelengths 
    a_chlD=0.0D3, & !.04D3, & !for light attenutation: plus DOM
    convI=2.77D18/6.02D23*86400D0, & !from W/m2 to Ein/m2/d 2.77e18[quanta/W/s]*1/6.02e23[Einstein/quanta]. from MBARI: http://www3.mbari.org/bog/nopp/par.html
	!
	RredO=467D0/4D0/16D0, &!mol O2/mol N to balance (=7.3). 9.4D0 for avg of  Anderson 1995. or, 467D0/4D0/16D0, &!from my stoich, i think (11/14/17) 1D0/(yo_bo*(1/yo_bo-1)), & !implied ratio of export production to balance surface oxygen (see 3/27/15 notes)
	!
	!UPTAKE parameterizations
	po_coef=2329100D0, & !m3/mol biomass N/day - multiply this by O2 (in mol/m3) and by oxygen yield (mol biomass N/mol O2) to get growth rate
    pnh4_max= 50.8D0, & !16.64D0, & !??Litchman (2-fold lower) mol NO3 uptake/mol cellular N/day 
    pnox_max= 50.8D0, & !8.32D0, & !Litchman (4-fold lower affinity) mol NO3 uptake/mol cellular N/day 
	knh4= 0.133D-3, & !0.0472D-3, &!(2fold lower for NH4)0.043D-3, & !mol/m3 - DIN uptake half-sat; Litchman
	knox= 0.133D-3, & !0.0944D-3, &!0.16D-3, & !mol/m3 - DIN uptake half-sat; Litchman
	!Effective ksats (Verdy)
    knh4_effp2= 0.164D-3, & !mol/m3 Litchman scaling, with aV^b with a for knh4 from Ward 2014 
    knox_effp2= 0.327D-3, & !mol/m3 Litchman
    knh4_effp3= 0.0018D-3, & !mol/m3 Litchman scaling, with aV^b with a for knh4 from Ward 2014 
    knox_effp3= 0.0036D-3, & !mol/m3 Litchman
    amminhib=4.6*1e3*1D0, & !m3/mol NH4 (from 1/uM nh4 in follows 2007)
	pd_max= 0.5D0, & !1/day - normalized max detritus uptake to match nitrate max uptake
	kd= 0.01D-3, & !org matter uptake half sat; used to be 5
	!D/POM processing:
	alpha=2D0, &
    mortf=0.5D0, & !fraction of mortality that goes to DOM vs POM
	!
	!OTHER
	mp=1D-2, &
	mb=1D-2, & !linear mortality
	mquad = 0.1D3, &!.1D3, & !1D3, &
	!
	!Zooplankton:
	mz=0.5D3, & !for quadratic mortality, so that mz*Z*Z = 1*0.1*0.1 and is like linear mz=0.1/day.
	g=2D3, &
	gam=0.2D0, & !growth yield for zoo w.r.t. B (10% efficiency..)
    kO2 = 1D-3, & !mol/m3 O2 limitation for zoo3
    !Temperature, as in GUD/Nitrify 1D model:
    TempAeArr = -4D3, &
    TemprefArr = 293.15D0, &
    Tkel = 273.15D0, &
    TempCoeffArr = 0.8D0

INTEGER :: t,mlboxes,i,j,ic,jc,nt,startSS,ind,recordDaily,dommodel,dailycycle
REAL*8 :: ym(ny) = (/(i,i=0+dy/2,Lp-dy/2, dy)/) !middle of boxes
REAL*8 :: y(ny+1) = (/(i,i=0,Lp, dy)/) !edges of boxes
REAL*8 :: zm(nz) = (/(j,j=0+dz/2,Hp-dz/2, dz)/)
REAL*8 :: z(nz+1) = (/(j,j=0,Hp, dz)/)
REAL*8 :: koverh, distn, adv, diff, cputime,dayint,dayint2,Iin !!!12/4- can take ozoo out after implementing MYRK!
REAL*8,DIMENSION(:),ALLOCATABLE :: time, burial
REAL*8,DIMENSION(:,:),ALLOCATABLE :: sumall
REAL*8,DIMENSION(3) :: specs
REAL*8,DIMENSION(nz+4) :: olost
REAL*8,DIMENSION(nz+1,ny+1) :: psi
REAL*8,DIMENSION(nz,ny+1) :: v,Ky
REAL*8,DIMENSION(nz+1,ny) :: w,wd,Kz,wd2
REAL*8,DIMENSION(nz+4,ny+4) :: eqmask,inmask,Iz, & !lam
	u_bo,u_bo_pa,po,pd,pdom,pnh4,pnh4b,pno2,pno2b,pno3,pno3b,pn2o,Temp,TempFun,Q10r,O2limZ, &  !these are replaced each time, don't need k and A,B,C terms.
	u_bdno3,u_bdno2,u_bdn2o,u_bana,u_bnh4,u_bno2, &
	u_p1,u_p2,u_p3,u_p3nolim,Biot,Chlt,nlimtot,limnh4,limno2,limno3,inhibnh4, & !phytopl vars
    limnh4b,limno2b,limno3b,limnh4p3,limno2p3,limno3p3, &
	nh4,no2,no3,n2o,ntot,o,d,dom,zoo,ozooall,olosszoo, & !
	bo,bo_pa,bdno3,bdno2,bdn2o,bnh4,bno2,bana,bt,btsq, & !8 types
	zoo2,zoo3,p1,p2,p3,pt,ptsq, & !3 p types and separate p grazer
	knh4A,knh4B,knh4C,knh4D,nh4A,nh4B,nh4C, &
	kno2A,kno2B,kno2C,kno2D,no2A,no2B,no2C, &
	kno3A,kno3B,kno3C,kno3D,no3A,no3B,no3C, &
	kn2oA,kn2oB,kn2oC,kn2oD,n2oA,n2oB,n2oC, &
	koA,koB,koC,koD,oA,oB,oC, &
	kboA,kboB,kboC,kboD,boA,boB,boC, &
	kbo_paA,kbo_paB,kbo_paC,kbo_paD,bo_paA,bo_paB,bo_paC, &
	kbnh4A,kbnh4B,kbnh4C,kbnh4D,bnh4A,bnh4B,bnh4C, &
	kbno2A,kbno2B,kbno2C,kbno2D,bno2A,bno2B,bno2C, &
	kbdno3A,kbdno3B,kbdno3C,kbdno3D,bdno3A,bdno3B,bdno3C, &
	kbdno2A,kbdno2B,kbdno2C,kbdno2D,bdno2A,bdno2B,bdno2C, &
	kbdn2oA,kbdn2oB,kbdn2oC,kbdn2oD,bdn2oA,bdn2oB,bdn2oC, &
	kbanaA,kbanaB,kbanaC,kbanaD,banaA,banaB,banaC, &
	kdA,kdB,kdC,kdD,dA,dB,dC,kdomA,kdomB,kdomC,kdomD,domA,domB,domC, &
	kzooA,kzooB,kzooC,kzooD,zooA,zooB,zooC, &
	kzoo2A,kzoo2B,kzoo2C,kzoo2D,zoo2A,zoo2B,zoo2C, &
	kzoo3A,kzoo3B,kzoo3C,kzoo3D,zoo3A,zoo3B,zoo3C, &
	kp1A,kp1B,kp1C,kp1D,p1A,p1B,p1C, &
	kp2A,kp2B,kp2C,kp2D,p2A,p2B,p2C, &
	kp3A,kp3B,kp3C,kp3D,p3A,p3B,p3C, &
    no3uptakeP,no2emitP, &
    PC1,PCmax1,PC2,PCmax2,PC3,PCmax3,a_Ip2,a_Ip3,chl2c,chl2c_p1,chl2c_p2,chl2c_p3, &!for Chl:C model
    domf !fraction of OM pool that is dom (vs particulate sinking d)

startSS=1 !SS right now is 1e5+18000 days - updated 8/29
dayint=100*1D0 !for recording each timestep for movies
dayint2=100*1D0 !for recording sums at each timestep for time series only, resulting in 10,000 points.
recordDaily=0 !for recording each dt for movies/2D
dommodel=1 !1 to include dom
dailycycle=1 !daily light cycle

specs(1)=ndays
specs(2)=dayint
specs(3)=dt

kdA(:,:)=0D0
kdB(:,:)=0D0
kdC(:,:)=0D0
kdD(:,:)=0D0
!
kdomA(:,:)=0D0
kdomB(:,:)=0D0
kdomC(:,:)=0D0
kdomD(:,:)=0D0
!
knh4A(:,:)=0D0
knh4B(:,:)=0D0
knh4C(:,:)=0D0
knh4D(:,:)=0D0
kno2A(:,:)=0D0
kno2B(:,:)=0D0
kno2C(:,:)=0D0
kno2D(:,:)=0D0
kno3A(:,:)=0D0
kno3B(:,:)=0D0
kno3C(:,:)=0D0
kno3D(:,:)=0D0
kn2oA(:,:)=0D0
kn2oB(:,:)=0D0
kn2oC(:,:)=0D0
kn2oD(:,:)=0D0
!
kboA(:,:)=0D0
kboB(:,:)=0D0
kboC(:,:)=0D0
kboD(:,:)=0D0
!
kbo_paA(:,:)=0D0
kbo_paB(:,:)=0D0
kbo_paC(:,:)=0D0
kbo_paD(:,:)=0D0
!
kbdno2A(:,:)=0D0
kbdno2B(:,:)=0D0
kbdno2C(:,:)=0D0
kbdno2D(:,:)=0D0
kbdno3A(:,:)=0D0
kbdno3B(:,:)=0D0
kbdno3C(:,:)=0D0
kbdno3D(:,:)=0D0
kbdn2oA(:,:)=0D0
kbdn2oB(:,:)=0D0
kbdn2oC(:,:)=0D0
kbdn2oD(:,:)=0D0
!
kbnh4A(:,:)=0D0
kbnh4B(:,:)=0D0
kbnh4C(:,:)=0D0
kbnh4D(:,:)=0D0
kbno2A(:,:)=0D0
kbno2B(:,:)=0D0
kbno2C(:,:)=0D0
kbno2D(:,:)=0D0
!
kbanaA(:,:)=0D0
kbanaB(:,:)=0D0
kbanaC(:,:)=0D0
kbanaD(:,:)=0D0
!
koA(:,:)=0D0
koB(:,:)=0D0
koC(:,:)=0D0
koD(:,:)=0D0
!
kzooA(:,:)=0D0
kzooB(:,:)=0D0
kzooC(:,:)=0D0
kzooD(:,:)=0D0
!
kzoo2A(:,:)=0D0
kzoo2B(:,:)=0D0
kzoo2C(:,:)=0D0
kzoo2D(:,:)=0D0
!
kzoo3A(:,:)=0D0
kzoo3B(:,:)=0D0
kzoo3C(:,:)=0D0
kzoo3D(:,:)=0D0
!
kp1A(:,:)=0D0
kp1B(:,:)=0D0
kp1C(:,:)=0D0
kp1D(:,:)=0D0
!
kp2A(:,:)=0D0
kp2B(:,:)=0D0
kp2C(:,:)=0D0
kp2D(:,:)=0D0
!
kp3A(:,:)=0D0
kp3B(:,:)=0D0
kp3C(:,:)=0D0
kp3D(:,:)=0D0


print*,'Run for total n of days:';print*,ndays
print*,'ny is:'; print*,ny
print*,'nz is:'; print*,nz

nt=ndays/dt

print*,'Number of days:'
print*,ndays
print*,'Number of timesteps:'
print*,nt

mlboxes=100/dz*1D0 !discrete n of boxes in the mixed layer, close to 100m total sum
koverh=Kg/100/mlboxes *60*60*24*1D0 !gas transfer coefficient for each of the n boxes comprising the ml

ALLOCATE(time(nt))
ALLOCATE(burial(nt+1))
ind=ndays/dayint2
ALLOCATE(sumall(ind,22))


!import velocity fields:

!fields u and w 4000 are from FLOW folder (run for 100 days at 4000m depth and checked for stability)
!fields u and w 2000 are from a 100 day run at 2000m depth, with **stronger flow!!** (realistic f for -10 degrees)

OPEN(UNIT=3,FILE='u2000_10m.txt',ACCESS='SEQUENTIAL',BLANK='ZERO') 
	DO I=1,ny+1
		read(3,*) (v(J,I),J=1,nz)
	END DO
CLOSE(3)

OPEN(UNIT=3,FILE='w2000_10m.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
	DO I=1,ny
		read(3,*) (w(J,I),J=1,nz+1)
	END DO
CLOSE(3)

!sinking speed for first pool: "d"
wd(:,:)=Ws
wd(1,:)=0D0
wd(nz+1,:)=0D0
wd=wd+w; !vertical velocity combination for detritus

!sinking speed for second pool: "dom"
wd2(:,:)=Ws2
wd2(1,:)=0D0
wd2(nz+1,:)=0D0
wd2=wd2+w; !vertical velocity combination for detritus


Temp(:,:)=0D0 !not sure if it matters but just in case

do j=1,nz
	!! took out July 2016: lam(j+2,3:ny+2)=lammax*exp(-zm(j)/euz) !lambda decreasing exp with depth
	!Kz(j,:)=(1e-1*exp(-zm(j)/euz)+kappazmin)*3600*24D0 !depth dependent mixing
	
	!had 5e-2 for bottom boundary peak value in paper
	Kz(j,:)=(1e-2*exp(-zm(j)/mlz)+kappazmin+1e-2*exp((zm(j)-H)/100))*3600*24D0 !larger at bottom boundary, too
	!Kz(j,:)=kappazmin*3600*24D0
	
				!n(j+2,3:ny+2)=0.03*(1-exp(-zm(j)/200)) ! n increases with depth

		Temp(j+2,3:ny+2)=11*exp(-zm(j)/100)+11*exp(-zm(j)/700)+2D0 !temperature to match WOCE data

end do

Q10r=exp(0.6931/10*(Temp-24)) !Q10 relationship
!from Nitr model (GUD-like temp):
TempFun = TempCoeffArr*exp(TempAeArr*(1D0/(Temp+Tkel)-1D0/TemprefArr))

Ky(:,:)=1e3*3600*24D0 !m2/s horizontal mixing
Ky(:,1)=0D0
Ky(:,ny+1)=0D0 !no flux out the sides

Kz(1,:)=0D0
Kz(nz+1,:)=0D0

eqmask(:,:)=0D0
eqmask(3:mlboxes+2,3:ny+2)=1D0; !mask for air-sea equilibration

inmask(:,:)=0D0
inmask(3:nz+2,3:ny+2)=1D0

!import previous steady state as IC:

if (startSS.eq.1) then 

		OPEN(UNIT=3,FILE='nh4_fSS_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
		DO I=1,ny+4
		read(3,*) (nh4(J,I),J=1,nz+4)
		END DO
		CLOSE(3)
		
		OPEN(UNIT=3,FILE='no2_fSS_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
		DO I=1,ny+4
		read(3,*) (no2(J,I),J=1,nz+4)
		END DO
		CLOSE(3)
		
		OPEN(UNIT=3,FILE='no3_fSS_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
		DO I=1,ny+4
		read(3,*) (no3(J,I),J=1,nz+4)
		END DO
		CLOSE(3)
		
		OPEN(UNIT=3,FILE='n2o_fSS_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
		DO I=1,ny+4
		read(3,*) (n2o(J,I),J=1,nz+4)
		END DO
		CLOSE(3)
		
		OPEN(UNIT=3,FILE='d_fSS_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
		DO I=1,ny+4
		read(3,*) (d(J,I),J=1,nz+4)
		END DO
		CLOSE(3)
		
		OPEN(UNIT=3,FILE='dom_fSS_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
		DO I=1,ny+4
		read(3,*) (dom(J,I),J=1,nz+4)
		END DO
		CLOSE(3)

		OPEN(UNIT=3,FILE='o_fSS_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
		DO I=1,ny+4
		read(3,*) (o(J,I),J=1,nz+4)
		END DO
		CLOSE(3)

		OPEN(UNIT=3,FILE='bo_fSS_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
		DO I=1,ny+4
		read(3,*) (bo(J,I),J=1,nz+4)
		END DO
		CLOSE(3)
		
		OPEN(UNIT=3,FILE='bo_pa_fSS_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
		DO I=1,ny+4
		read(3,*) (bo_pa(J,I),J=1,nz+4)
		END DO
		CLOSE(3)
		
		OPEN(UNIT=3,FILE='bdno3_fSS_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
		DO I=1,ny+4
		read(3,*) (bdno3(J,I),J=1,nz+4)
		END DO
		CLOSE(3)
		
		OPEN(UNIT=3,FILE='bdno2_fSS_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
		DO I=1,ny+4
		read(3,*) (bdno2(J,I),J=1,nz+4)
		END DO
		CLOSE(3)
		
		OPEN(UNIT=3,FILE='bdn2o_fSS_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
		DO I=1,ny+4
		read(3,*) (bdn2o(J,I),J=1,nz+4)
		END DO
		CLOSE(3)
		
		OPEN(UNIT=3,FILE='bana_fSS_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
		DO I=1,ny+4
		read(3,*) (bana(J,I),J=1,nz+4)
		END DO
		CLOSE(3)
		
		OPEN(UNIT=3,FILE='bnh4_fSS_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
		DO I=1,ny+4
		read(3,*) (bnh4(J,I),J=1,nz+4)
		END DO
		CLOSE(3)
		
		OPEN(UNIT=3,FILE='bno2_fSS_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
		DO I=1,ny+4
		read(3,*) (bno2(J,I),J=1,nz+4)
		END DO
		CLOSE(3)

		OPEN(UNIT=3,FILE='zoo_fSS_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
		DO I=1,ny+4
		read(3,*) (zoo(J,I),J=1,nz+4)
		END DO
		CLOSE(3)
		
		OPEN(UNIT=3,FILE='zoo2_fSS_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
		DO I=1,ny+4
		read(3,*) (zoo2(J,I),J=1,nz+4)
		END DO
		CLOSE(3)
		
		OPEN(UNIT=3,FILE='zoo3_fSS_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
		DO I=1,ny+4
		read(3,*) (zoo3(J,I),J=1,nz+4)
		END DO
		CLOSE(3)
		
		OPEN(UNIT=3,FILE='p1_fSS_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
		DO I=1,ny+4
		read(3,*) (p1(J,I),J=1,nz+4)
		END DO
		CLOSE(3)
	
        OPEN(UNIT=3,FILE='p2_fSS_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
		DO I=1,ny+4
		read(3,*) (p2(J,I),J=1,nz+4)
		END DO
		CLOSE(3)
		
		OPEN(UNIT=3,FILE='p3_fSS_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
		DO I=1,ny+4
		read(3,*) (p3(J,I),J=1,nz+4)
		END DO
		CLOSE(3)

else
		!!Initial Conditions:
		nh4(:,:)=0D0 !inmask*1e-3*1D0
		no2(:,:)=0D0 !inmask*1e-3*1D0
		no3(:,:)=inmask*30*1e-3*1D0 
		n2o(:,:)=0D0 !start with none bc there's no sink...
		
		bo(:,:)=inmask*1e-4*1D0
		bo_pa(:,:)=inmask*1e-4*1D0		
		bnh4(:,:)=inmask*1e-4*1D0
		bno2(:,:)=inmask*1e-4*1D0
		bdno3(:,:)=inmask*1e-4*1D0
		bdno2(:,:)=inmask*1e-4*1D0
		bdn2o(:,:)=inmask*1e-4*1D0
		bana(:,:)=inmask*1e-4*1D0
		p1(:,:)=inmask*1e-4*1D0
		p2(:,:)=inmask*1e-4*1D0
		p3(:,:)=inmask*1e-4*1D0
		
		o(:,:)=inmask*150*1024/1e6*1D0 !mol/m3 crude estimate 
		
		zoo(:,:)=inmask*1e-4*1D0
		zoo2(:,:)=inmask*1e-4*1D0
		zoo3(:,:)=inmask*1e-4*1D0

		
		do j=1,nz
			!ndecay profile fit:
			no3(j+2,3:ny+2)=0.03D0*(1-exp(-zm(j)/200)) ! n increases with depth
		end do

		d(:,:)=0D0 !start with none
		dom(:,:)=0D0
		nh4(:,:)=0D0
		no2(:,:)=0D0

!start with SS oxygen, if not starting at SS:
		OPEN(UNIT=3,FILE='o_fSS_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
		DO I=1,ny+4
		read(3,*) (o(J,I),J=1,nz+4)
		END DO
		CLOSE(3)
		
end if

!no N2O:
bdn2o(:,:)=0D0

!only one P:
p1(:,:)=0D0

if (dommodel.eq.0) then 
	dom(:,:)=0D0
	bo_pa(:,:)=0D0
end if
!but i am also not using pa bacteria now, even with dommodel on, so:
bo_pa(:,:)=0D0

OPEN(UNIT=5,FILE='time_record.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',status='REPLACE')
CLOSE(5)
OPEN(UNIT=5,FILE='oxygenproblems.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',STATUS='REPLACE')
CLOSE(5)

if (recordDaily.eq.1) then

OPEN(UNIT=5,FILE='Iin.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',status='REPLACE')
CLOSE(5)

OPEN(UNIT=5,FILE='Iz.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',status='REPLACE')
CLOSE(5)

OPEN(UNIT=5,FILE='time_all.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',status='REPLACE')
CLOSE(5)

OPEN(UNIT=5,FILE='nh4_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',STATUS='REPLACE')
CLOSE(5)	
OPEN(UNIT=5,FILE='no2_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',STATUS='REPLACE')
CLOSE(5)
OPEN(UNIT=5,FILE='no3_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',STATUS='REPLACE')
CLOSE(5)
OPEN(UNIT=5,FILE='d_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',status='REPLACE')
CLOSE(5)
OPEN(UNIT=5,FILE='dom_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',status='REPLACE')
CLOSE(5)
OPEN(UNIT=5,FILE='o_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',status='REPLACE')
CLOSE(5)	
OPEN(UNIT=5,FILE='bo_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',status='REPLACE')
CLOSE(5)	
OPEN(UNIT=5,FILE='bo_pa_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',status='REPLACE')
CLOSE(5)	
OPEN(UNIT=5,FILE='bnh4_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',status='REPLACE')
CLOSE(5)
OPEN(UNIT=5,FILE='bno2_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',status='REPLACE')
CLOSE(5)
OPEN(UNIT=5,FILE='ubo_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',status='REPLACE')
CLOSE(5)
OPEN(UNIT=5,FILE='ubnh4_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',status='REPLACE')
CLOSE(5)
OPEN(UNIT=5,FILE='ubno2_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',status='REPLACE')
CLOSE(5)
OPEN(UNIT=5,FILE='up_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',status='REPLACE')
CLOSE(5)
OPEN(UNIT=5,FILE='zoo_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',status='REPLACE')
CLOSE(5)
OPEN(UNIT=5,FILE='zoo2_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',status='REPLACE')
CLOSE(5)
OPEN(UNIT=5,FILE='zoo3_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',status='REPLACE')
CLOSE(5)
OPEN(UNIT=5,FILE='p1_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',status='REPLACE')
CLOSE(5)
OPEN(UNIT=5,FILE='p2_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',status='REPLACE')
CLOSE(5)
OPEN(UNIT=5,FILE='p3_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',status='REPLACE')
CLOSE(5)


end if 

print *,'Starting time loop:'
do t=1,nt


!LIGHT:
if (dailycycle.eq.1) then !Incoming light daily cycle:
    Iin=Iinmax/2D0*(cos(t*dt*2D0*3.1416D0)+1D0)
else !No daily cycle:
    Iin=Iinmax/2D0
end if

!total biomass:
!Biot=bo+bo_pa+bnh4+bno2+p1+p2+p3+bana+bdno3+bdno2 !just pp and b: zoo and dom? ask steph
!chl impact on light:
Chlt=(p1*chl2c_p1+p2*chl2c_p2+p3*chl2c_p3)*6.6D0 !molN/m3 *mgChl/mmolC *5molC/molN = gChl

do i=1,ny
	do j=1,nz

		!Iz(j+2,i+2)=Iin*exp(-zm(j)*(1/euz+sum(Biot(3:j+2,i+2)*kBio))) !sums at each depth like dutk.
        !updated:
        Iz(j+2,i+2)=Iin*exp(-zm(j)*(1/euz+sum(Chlt(3:j+2,i+2)*a_chlD))) !this one sums the k at each depth

		
	end do
end do


i=dayint2*1000
j=t*dt*1000

if ((MOD(j,i).eq.0).OR.(t.eq.1)) then

	!trace the integral:
	ind=t*dt/dayint2
	
	!print*,ind
	
	time(ind)=(t-1)*dt !start at zero- weird just bc i can't just easily write out the vector

	ntot=nh4+no2+no3+n2o

	sumall(ind,1)=sum(o)*Aij !mol/m3 times volume (with dx=1,dy=1)
	sumall(ind,2)=sum(d)*Aij
	sumall(ind,3)=sum(dom)*Aij
	sumall(ind,4)=sum(bo)*Aij
	sumall(ind,5)=sum(bo_pa)*Aij
	sumall(ind,6)=sum(bdno3)*Aij
	sumall(ind,7)=sum(bdno2)*Aij
	sumall(ind,8)=sum(bdn2o)*Aij
	sumall(ind,9)=sum(bana)*Aij
	sumall(ind,10)=sum(bnh4)*Aij
	sumall(ind,11)=sum(bno2)*Aij
	sumall(ind,12)=sum(zoo)*Aij
	sumall(ind,13)=sum(ntot)*Aij
	sumall(ind,14)=sum(nh4)*Aij
	sumall(ind,15)=sum(no2)*Aij
	sumall(ind,16)=sum(no3)*Aij
	sumall(ind,17)=sum(n2o)*Aij
	sumall(ind,18)=sum(p1)*Aij
	sumall(ind,19)=sum(p2)*Aij
	sumall(ind,20)=sum(p3)*Aij
	sumall(ind,21)=sum(zoo2)*Aij
	sumall(ind,22)=sum(zoo3)*Aij

end if

call MYRK(nh4,no3,no2,n2o,d,dom,o, &
				zoo,zoo2,zoo3, &
				p1,p2,p3, &
				bo,bo_pa,bdno3,bdno2,bdn2o, &
				bnh4,bno2,bana, & 
				knh4A,kno3A,kno2A,kn2oA,kdA,kdomA,koA, &
				kzooA,kzoo2A,kzoo3A, &
				kp1A,kp2A,kp3A, &
				kboA,kbo_paA,kbdno3A,kbdno2A,kbdn2oA, &
				kbnh4A,kbno2A,kbanaA, &
				nh4A,no3A,no2A,n2oA,dA,domA,oA, &
				zooA,zoo2A,zoo3A, &
				p1A,p2A,p3A, &
				boA,bo_paA,bdno3A,bdno2A,bdn2oA, &
				bnh4A,bno2A,banaA) !, & 
				!po,pd,pdom,pnh4,pno2,pno3,pn2o,u_bo,u_bo_pa,u_bdno3,u_bdno2,u_bdn2o, &
				!u_bnh4,u_bno2,u_bana,distn)				
		

call MYRK(nh4A,no3A,no2A,n2oA,dA,domA,oA,zooA,zoo2A,zoo3A,p1A,p2A,p3A, &
				boA,bo_paA,bdno3A,bdno2A,bdn2oA,bnh4A,bno2A,banaA, & 
				knh4B,kno3B,kno2B,kn2oB,kdB,kdomB,koB,kzooB,kzoo2B,kzoo3B,kp1B,kp2B,kp3B, &
				kboB,kbo_paB,kbdno3B,kbdno2B,kbdn2oB,kbnh4B,kbno2B,kbanaB, &
				nh4B,no3B,no2B,n2oB,dB,domB,oB,zooB,zoo2B,zoo3B,p1B,p2B,p3B, &
				boB,bo_paB,bdno3B,bdno2B,bdn2oB,bnh4B,bno2B,banaB) !, & 
				
call MYRK(nh4B,no3B,no2B,n2oB,dB,domB,oB,zooB,zoo2B,zoo3B,p1B,p2B,p3B, &
				boB,bo_paB,bdno3B,bdno2B,bdn2oB,bnh4B,bno2B,banaB, & 
				knh4C,kno3C,kno2C,kn2oC,kdC,kdomC,koC,kzooC,kzoo2C,kzoo3C,kp1C,kp2C,kp3C, &
				kboC,kbo_paC,kbdno3C,kbdno2C,kbdn2oC,kbnh4C,kbno2C,kbanaC, &
				nh4C,no3C,no2C,n2oC,dC,domC,oC,zooC,zoo2C,zoo3C,p1C,p2C,p3C, &
				boC,bo_paC,bdno3C,bdno2C,bdn2oC,bnh4C,bno2C,banaC) !, & 


call MYRK(nh4C,no3C,no2C,n2oC,dC,domC,oC,zooC,zoo2C,zoo3C,p1C,p2C,p3C, &
				boC,bo_paC,bdno3C,bdno2C,bdn2oC,bnh4C,bno2C,banaC, & 
				knh4D,kno3D,kno2D,kn2oD,kdD,kdomD,koD,kzooD,kzoo2D,kzoo3D,kp1D,kp2D,kp3D, &
				kboD,kbo_paD,kbdno3D,kbdno2D,kbdn2oD,kbnh4D,kbno2D,kbanaD, &
				nh4A,no3A,no2A,n2oA,dA,domA,oA,zooA,zoo2A,zoo3A,p1A,p2A,p3A, & !a placeholder- unneccessary calc		
				boA,bo_paA,bdno3A,bdno2A,bdn2oA,bnh4A,bno2A,banaA) !, & 
								
nh4 = nh4 + dt/6D0*(knh4A + 2D0*knh4B + 2D0*knh4C + knh4D);
no2 = no2 + dt/6D0*(kno2A + 2D0*kno2B + 2D0*kno2C + kno2D);
no3 = no3 + dt/6D0*(kno3A + 2D0*kno3B + 2D0*kno3C + kno3D);
n2o = n2o + dt/6D0*(kn2oA + 2D0*kn2oB + 2D0*kn2oC + kn2oD);
bo = bo + dt/6D0*(kboA + 2D0*kboB + 2D0*kboC + kboD);
bo_pa = bo_pa + dt/6D0*(kbo_paA + 2D0*kbo_paB + 2D0*kbo_paC + kbo_paD);
bdno3 = bdno3 + dt/6D0*(kbdno3A + 2D0*kbdno3B + 2D0*kbdno3C + kbdno3D);
bdno2 = bdno2 + dt/6D0*(kbdno2A + 2D0*kbdno2B + 2D0*kbdno2C + kbdno2D);
bnh4 = bnh4 + dt/6D0*(kbnh4A + 2D0*kbnh4B + 2D0*kbnh4C + kbnh4D);
bno2 = bno2 + dt/6D0*(kbno2A + 2D0*kbno2B + 2D0*kbno2C + kbno2D);
bana = bana + dt/6D0*(kbanaA + 2D0*kbanaB + 2D0*kbanaC + kbanaD);
d = d + dt/6D0*(kdA + 2D0*kdB + 2D0*kdC + kdD);
dom = dom + dt/6D0*(kdomA + 2D0*kdomB + 2D0*kdomC + kdomD);
o = o + dt/6D0*(koA + 2D0*koB + 2D0*koC + koD);	
zoo = zoo + dt/6D0*(kzooA + 2D0*kzooB + 2D0*kzooC + kzooD);			
!p1 = p1 + dt/6D0*(kp1A + 2*kp1B + 2*kp1C + kp1D);	
p2 = p2 + dt/6D0*(kp2A + 2D0*kp2B + 2D0*kp2C + kp2D);	
p3 = p3 + dt/6D0*(kp3A + 2D0*kp3B + 2D0*kp3C + kp3D);	
zoo2 = zoo2 + dt/6D0*(kzoo2A + 2D0*kzoo2B + 2D0*kzoo2C + kzoo2D);		
zoo3 = zoo3 + dt/6D0*(kzoo3A + 2D0*kzoo3B + 2D0*kzoo3C + kzoo3D);		

if (MOD(t*dt,1.00).eq.0) then
	print*,t*dt		
	OPEN(UNIT=7,FILE='time_record.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',POSITION='APPEND')
	WRITE(7,*) (t*dt)
	CLOSE(7)	
end if

if (recordDaily.eq.1) then

!append at every time step:
print*,t*dt		
	OPEN(UNIT=7,FILE='Iin.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',POSITION='APPEND')
	WRITE(7,*) (Iin)
	CLOSE(7)

OPEN(UNIT=5,FILE='Iz.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',POSITION='APPEND')
WRITE(5,*) (Iz)
CLOSE(5)	

print*,t*dt		
	OPEN(UNIT=7,FILE='time_all.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',POSITION='APPEND')
	WRITE(7,*) (t*dt)
	CLOSE(7)	

OPEN(UNIT=5,FILE='nh4_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',POSITION='APPEND')
WRITE(5,*) (nh4)
CLOSE(5)	
OPEN(UNIT=5,FILE='no2_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',POSITION='APPEND')
WRITE(5,*) (no2)
CLOSE(5)
OPEN(UNIT=5,FILE='no3_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',POSITION='APPEND')
WRITE(5,*) (no3)
CLOSE(5)
OPEN(UNIT=5,FILE='d_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',POSITION='APPEND')
WRITE(5,*) (d)
CLOSE(5)
OPEN(UNIT=5,FILE='dom_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',POSITION='APPEND')
WRITE(5,*) (dom)
CLOSE(5)
OPEN(UNIT=5,FILE='o_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',POSITION='APPEND')
WRITE(5,*) (o)
CLOSE(5)	
OPEN(UNIT=5,FILE='bo_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',POSITION='APPEND')
WRITE(5,*) (bo)
CLOSE(5)	
OPEN(UNIT=5,FILE='bo_pa_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',POSITION='APPEND')
WRITE(5,*) (bo_pa)
CLOSE(5)	
OPEN(UNIT=5,FILE='bnh4_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',POSITION='APPEND')
WRITE(5,*) (bnh4)
CLOSE(5)
OPEN(UNIT=5,FILE='bno2_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',POSITION='APPEND')
WRITE(5,*) (bno2)
CLOSE(5)
OPEN(UNIT=5,FILE='ubo_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',POSITION='APPEND')
WRITE(5,*) (u_bo)
CLOSE(5)
OPEN(UNIT=5,FILE='ubnh4_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',POSITION='APPEND')
WRITE(5,*) (u_bnh4)
CLOSE(5)
OPEN(UNIT=5,FILE='ubno2_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',POSITION='APPEND')
WRITE(5,*) (u_bno2)
CLOSE(5)
OPEN(UNIT=5,FILE='up1_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',POSITION='APPEND')
WRITE(5,*) (u_p1)
CLOSE(5)
OPEN(UNIT=5,FILE='up2_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',POSITION='APPEND')
WRITE(5,*) (u_p2)
CLOSE(5)
OPEN(UNIT=5,FILE='up3_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',POSITION='APPEND')
WRITE(5,*) (u_p3)
CLOSE(5)
OPEN(UNIT=5,FILE='zoo_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',POSITION='APPEND')
WRITE(5,*) (zoo)
CLOSE(5)
OPEN(UNIT=5,FILE='zoo2_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',POSITION='APPEND')
WRITE(5,*) (zoo2)
CLOSE(5)
OPEN(UNIT=5,FILE='zoo3_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',POSITION='APPEND')
WRITE(5,*) (zoo3)
CLOSE(5)
OPEN(UNIT=5,FILE='p1_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',POSITION='APPEND')
WRITE(5,*) (p1)
CLOSE(5)	
OPEN(UNIT=5,FILE='p2_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',POSITION='APPEND')
WRITE(5,*) (p2)
CLOSE(5)	
OPEN(UNIT=5,FILE='p3_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',POSITION='APPEND')
WRITE(5,*) (p3)
CLOSE(5)	

end if
	

if ((MOD(t*dt,dayint).eq.0).OR.(t*dt.eq.ndays).OR.(t.eq.1)) then

	
OPEN(UNIT=5,FILE='specs_f_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
DO I=1,3
WRITE(5,*) (specs(I))
END DO
CLOSE(5)

OPEN(UNIT=5,FILE='wd_f_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
DO I=1,ny
WRITE(5,*) (wd(J,I),J=1,nz+1)
END DO
CLOSE(5)

OPEN(UNIT=5,FILE='kz_f_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
DO I=1,ny
WRITE(5,*) (kz(J,I),J=1,nz+1)
END DO
CLOSE(5)

OPEN(UNIT=5,FILE='Iz_f_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
DO I=1,ny+4
WRITE(5,*) (Iz(J,I),J=1,nz+4)
END DO
CLOSE(5)

OPEN(UNIT=5,FILE='nh4_f_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
DO I=1,ny+4
WRITE(5,*) (nh4(J,I),J=1,nz+4)
END DO
CLOSE(5)

OPEN(UNIT=5,FILE='no3_f_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
DO I=1,ny+4
WRITE(5,*) (no3(J,I),J=1,nz+4)
END DO
CLOSE(5)

OPEN(UNIT=5,FILE='no2_f_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
DO I=1,ny+4
WRITE(5,*) (no2(J,I),J=1,nz+4)
END DO
CLOSE(5)

OPEN(UNIT=5,FILE='n2o_f_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
DO I=1,ny+4
WRITE(5,*) (n2o(J,I),J=1,nz+4)
END DO
CLOSE(5)

OPEN(UNIT=5,FILE='d_f_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
DO I=1,ny+4
WRITE(5,*) (d(J,I),J=1,nz+4)
END DO
CLOSE(5)

OPEN(UNIT=5,FILE='dom_f_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
DO I=1,ny+4
WRITE(5,*) (dom(J,I),J=1,nz+4)
END DO
CLOSE(5)

OPEN(UNIT=5,FILE='o_f_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
DO I=1,ny+4
WRITE(5,*) (o(J,I),J=1,nz+4)
END DO
CLOSE(5)

OPEN(UNIT=5,FILE='bo_f_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
DO I=1,ny+4
WRITE(5,*) (bo(J,I),J=1,nz+4)
END DO
CLOSE(5)

OPEN(UNIT=5,FILE='bo_pa_f_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
DO I=1,ny+4
WRITE(5,*) (bo_pa(J,I),J=1,nz+4)
END DO
CLOSE(5)

OPEN(UNIT=5,FILE='bdno3_f_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
DO I=1,ny+4
WRITE(5,*) (bdno3(J,I),J=1,nz+4)
END DO
CLOSE(5)

OPEN(UNIT=5,FILE='bdno2_f_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
DO I=1,ny+4
WRITE(5,*) (bdno2(J,I),J=1,nz+4)
END DO
CLOSE(5)

OPEN(UNIT=5,FILE='bdn2o_f_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
DO I=1,ny+4
WRITE(5,*) (bdn2o(J,I),J=1,nz+4)
END DO
CLOSE(5)

OPEN(UNIT=5,FILE='bnh4_f_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
DO I=1,ny+4
WRITE(5,*) (bnh4(J,I),J=1,nz+4)
END DO
CLOSE(5)

OPEN(UNIT=5,FILE='bno2_f_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
DO I=1,ny+4
WRITE(5,*) (bno2(J,I),J=1,nz+4)
END DO
CLOSE(5)

OPEN(UNIT=5,FILE='bana_f_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
DO I=1,ny+4
WRITE(5,*) (bana(J,I),J=1,nz+4)
END DO
CLOSE(5)

OPEN(UNIT=5,FILE='zoo_f_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
DO I=1,ny+4
WRITE(5,*) (zoo(J,I),J=1,nz+4)
END DO
CLOSE(5)

OPEN(UNIT=5,FILE='zoo2_f_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
DO I=1,ny+4
WRITE(5,*) (zoo2(J,I),J=1,nz+4)
END DO
CLOSE(5)

OPEN(UNIT=5,FILE='zoo3_f_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
DO I=1,ny+4
WRITE(5,*) (zoo3(J,I),J=1,nz+4)
END DO
CLOSE(5)

OPEN(UNIT=5,FILE='p1_f_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
DO I=1,ny+4
WRITE(5,*) (p1(J,I),J=1,nz+4)
END DO
CLOSE(5)

OPEN(UNIT=5,FILE='p2_f_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
DO I=1,ny+4
WRITE(5,*) (p2(J,I),J=1,nz+4)
END DO
CLOSE(5)

OPEN(UNIT=5,FILE='p3_f_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
DO I=1,ny+4
WRITE(5,*) (p3(J,I),J=1,nz+4)
END DO
CLOSE(5)

OPEN(UNIT=5,FILE='time_f_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
DO I=1,ind
WRITE(5,*) (time(I))
END DO
CLOSE(5)

OPEN(UNIT=5,FILE='y_f_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
DO I=1,ny
WRITE(5,*) (ym(I))
END DO
CLOSE(5)

OPEN(UNIT=5,FILE='z_f_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
DO I=1,nz
WRITE(5,*) (zm(I))
END DO
CLOSE(5)

OPEN(UNIT=5,FILE='sumall_f_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
DO I=1,ind
WRITE(5,*) (sumall(I,J),J=1,22)
END DO
CLOSE(5)


!!Add oxygen budget terms:

!OPEN(UNIT=5,FILE='exportO_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
!DO I=1,ny+4
!!WRITE(5,*) (1/yoe*lam(J,I)*n(J,I),J=1,nz+4)
!WRITE(5,*) (RredO*lam(J,I)*n(J,I),J=1,nz+4)
!END DO
!CLOSE(5)

!OPEN(UNIT=5,FILE='bactOcons_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
!DO I=1,ny+4
!WRITE(5,*) (-(u_bo(J,I)*bo(J,I)+u_bo_pa(J,I)*bo_pa(J,I))/yoe,J=1,nz+4)
!END DO
!CLOSE(5)

!for this first one, olosszoo is actually this term
OPEN(UNIT=5,FILE='zoo1_Ocons_f_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
DO I=1,ny+4
WRITE(5,*) (-olosszoo(J,I),J=1,nz+4)
END DO
CLOSE(5)

!now, use olosszoo as a placeholder to record other zoo O2 consumptions
olosszoo=RredO*(1D0-gam)*g*p2*zoo2*O2limZ*TempFun
OPEN(UNIT=5,FILE='zoo2_Ocons_f_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
DO I=1,ny+4
WRITE(5,*) (-olosszoo(J,I),J=1,nz+4)
END DO
CLOSE(5)

olosszoo=RredO*(1D0-gam)*g*p3*zoo3*O2limZ*TempFun
OPEN(UNIT=5,FILE='zoo3_Ocons_f_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
DO I=1,ny+4
WRITE(5,*) (-olosszoo(J,I),J=1,nz+4)
END DO
CLOSE(5)


OPEN(UNIT=5,FILE='airsea_f_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
DO I=1,ny+4
WRITE(5,*) (koverh*(o2sat-o(J,I))*eqmask(J,I),J=1,nz+4)
END DO
CLOSE(5)

do i=1,ny
	do j=1,nz;	ic=i+2; jc=j+2	
		call myupwind(o,v,w,j,i,dy,dz,Aij,nz,ny,adv); !advection scheme outputs advective fluxes         
		call mydiff(o,Ky,Kz,j,i,dy,dz,Aij,nz,ny,diff)
		koD(jc,ic)=-adv+diff
	end do
end do

OPEN(UNIT=5,FILE='advdiff_f_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
DO I=1,ny+4
WRITE(5,*) (koD(J,I),J=1,nz+4)
END DO
CLOSE(5)

OPEN(UNIT=5,FILE='ddtO_f_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
DO I=1,ny+4
WRITE(5,*) (koA(J,I),J=1,nz+4)
END DO
CLOSE(5)

OPEN(UNIT=5,FILE='u_bo_f_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
DO I=1,ny+4
WRITE(5,*) (u_bo(J,I),J=1,nz+4)
END DO
CLOSE(5)

OPEN(UNIT=5,FILE='u_bo_pa_f_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
DO I=1,ny+4
WRITE(5,*) (u_bo_pa(J,I),J=1,nz+4)
END DO
CLOSE(5)

OPEN(UNIT=5,FILE='u_bdno3_f_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
DO I=1,ny+4
WRITE(5,*) (u_bdno3(J,I),J=1,nz+4)
END DO
CLOSE(5)

OPEN(UNIT=5,FILE='u_bdno2_f_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
DO I=1,ny+4
WRITE(5,*) (u_bdno2(J,I),J=1,nz+4)
END DO
CLOSE(5)

OPEN(UNIT=5,FILE='u_bana_f_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
DO I=1,ny+4
WRITE(5,*) (u_bana(J,I),J=1,nz+4)
END DO
CLOSE(5)

OPEN(UNIT=5,FILE='u_bnh4_f_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
DO I=1,ny+4
WRITE(5,*) (u_bnh4(J,I),J=1,nz+4)
END DO
CLOSE(5)

OPEN(UNIT=5,FILE='u_bno2_f_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
DO I=1,ny+4
WRITE(5,*) (u_bno2(J,I),J=1,nz+4)
END DO
CLOSE(5)

OPEN(UNIT=5,FILE='u_p3_f_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
DO I=1,ny+4
WRITE(5,*) (u_p3(J,I),J=1,nz+4)
END DO
CLOSE(5)

OPEN(UNIT=5,FILE='u_p2_f_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
DO I=1,ny+4
WRITE(5,*) (u_p2(J,I),J=1,nz+4)
END DO
CLOSE(5)

OPEN(UNIT=5,FILE='chl2c_p1_f_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
DO I=1,ny+4
WRITE(5,*) (chl2c_p1(J,I),J=1,nz+4)
END DO
CLOSE(5)

OPEN(UNIT=5,FILE='chl2c_p2_f_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
DO I=1,ny+4
WRITE(5,*) (chl2c_p2(J,I),J=1,nz+4)
END DO
CLOSE(5)

OPEN(UNIT=5,FILE='chl2c_p3_f_2d.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
DO I=1,ny+4
WRITE(5,*) (chl2c_p3(J,I),J=1,nz+4)
END DO
CLOSE(5)

    end if  !end the mod- recording after X days
	
	

end do !end time loop


print*,'Total CPU time in seconds:'
call CPU_TIME(cputime)
print*,cputime


contains
SUBROUTINE MYQUICK(C,v,w,j,i,dy,dz,Aij,nz,ny,adv)
implicit none
REAL*8 :: C(nz+4,ny+4),v(nz,ny+1),w(nz+1,ny),adv
REAL*8 :: vp1,vn1,vp,vn,wp1,wn1,wp,wn,Dx1,Dx,Dxn1,Dy1,Dy2,Dyn1,Fr,Fl,Fu,Fd
INTEGER :: j,i,Aij,nz,ny,dy,dz
INTEGER :: jc,ic

jc=j+2
ic=i+2
		!at right face:
        vp1=(v(j,i+1)+abs(v(j,i+1)))/2D0;
        vn1=(v(j,i+1)-abs(v(j,i+1)))/2D0;
        !at left face:
        vp=(v(j,i)+abs(v(j,i)))/2D0;
        vn=(v(j,i)-abs(v(j,i)))/2D0;
        !at top face:
        wp1=(w(j+1,i)+abs(w(j+1,i)))/2D0;
        wn1=(w(j+1,i)-abs(w(j+1,i)))/2D0;
        !at bottom face:
        wp=(w(j,i)+abs(w(j,i)))/2D0;
        wn=(w(j,i)-abs(w(j,i)))/2D0;
        
       !GET DxC's and DyC's
                            
        Dx1=C(jc,ic+2)-2D0*C(jc,ic+1)+C(jc,ic);
        Dx=C(jc,ic+1)-2D0*C(jc,ic)+C(jc,ic-1);
        Dxn1=C(jc,ic)-2D0*C(jc,ic-1)+C(jc,ic-2);
               
        Dy1=C(jc+2,ic)-2D0*C(jc+1,ic)+C(jc,ic);
        Dy2=C(jc+1,ic)-2D0*C(jc,ic)+C(jc-1,ic);
        Dyn1=C(jc,ic)-2D0*C(jc-1,ic)+C(jc-2,ic);
        
       !CALC FLUXES at each face:
        Fr=v(j,i+1)/2D0*(C(jc,ic)+C(jc,ic+1)) - vp1/8D0*Dx - vn1/8D0*Dx1;
        Fl=v(j,i)/2D0*(C(jc,ic-1)+C(jc,ic)) - vp/8D0*Dxn1 - vn/8D0*Dx;
        Fu=w(j+1,i)/2D0*(C(jc,ic)+C(jc+1,ic)) - wp1/8D0*Dy2 - wn1/8D0*Dy1;
        Fd=w(j,i)/2D0*(C(jc-1,ic)+C(jc,ic)) - wp/8D0*Dyn1 - wn/8D0*Dy2;
                
        Fr=Fr*dz;
        Fl=Fl*dz;
        Fu=Fu*dy;
        Fd=Fd*dy;
        
        adv=(Fr-Fl+Fu-Fd)/Aij;  
        
END SUBROUTINE MYQUICK

SUBROUTINE MYUPWIND(C,v,w,j,i,dy,dz,Aij,nz,ny,adv)
implicit none
REAL*8 :: C(nz+4,ny+4),v(nz,ny+1),w(nz+1,ny),adv
REAL*8 :: vp1,vn1,vp,vn,wp1,wn1,wp,wn,Dx1,Dx,Dxn1,Dy1,Dy2,Dyn1,Fr,Fl,Fu,Fd
INTEGER :: j,i,Aij,nz,ny,dy,dz
INTEGER :: jc,ic

jc=j+2
ic=i+2
		!at right face:
        vp1=(v(j,i+1)+abs(v(j,i+1)))/2D0;
        vn1=(v(j,i+1)-abs(v(j,i+1)))/2D0;
        !at left face:
        vp=(v(j,i)+abs(v(j,i)))/2D0;
        vn=(v(j,i)-abs(v(j,i)))/2D0;
        !at top face:
        wp1=(w(j+1,i)+abs(w(j+1,i)))/2D0;
        wn1=(w(j+1,i)-abs(w(j+1,i)))/2D0;
        !at bottom face:
        wp=(w(j,i)+abs(w(j,i)))/2D0;
        wn=(w(j,i)-abs(w(j,i)))/2D0;
       
        Fr=vp1*C(jc,ic)+vn1*C(jc,ic+1);
        Fl=vp*C(jc,ic-1)+vn*C(jc,ic);
        Fu=wp1*C(jc,ic)+wn1*C(jc+1,ic);
        Fd=wp*C(jc-1,ic)+wn*C(jc,ic);
               
        Fr=Fr*dz;
        Fl=Fl*dz;
        Fu=Fu*dy;
        Fd=Fd*dy;
        
        adv=(Fr-Fl+Fu-Fd)/Aij;

END SUBROUTINE MYUPWIND

SUBROUTINE MYDIFF(C,Ky,Kz,j,i,dy,dz,Aij,nz,ny,diff)
implicit none
REAL*8 :: C(nz+4,ny+4),Ky(nz,ny+1),Kz(nz+1,ny),diff
REAL*8 :: Fr,Fl,Fu,Fd
INTEGER :: j,i,Aij,nz,ny,dy,dz
INTEGER :: jc,ic

jc=j+2
ic=i+2
        
        Fr=Ky(j,i+1)*(C(jc,ic+1)-C(jc,ic))/dy;
        Fl=Ky(j,i)*(C(jc,ic)-C(jc,ic-1))/dy;
        Fu=Kz(j+1,i)*(C(jc+1,ic)-C(jc,ic))/dz;
        Fd=Kz(j,i)*(C(jc,ic)-C(jc-1,ic))/dz;
               
        Fr=Fr*dz; !this cancels out dividing by dz below (Aij) to get just d/dy(Flux in y)
        Fl=Fl*dz;
        Fu=Fu*dy;
        Fd=Fd*dy;
        
        diff=(Fr-Fl+Fu-Fd)/Aij;
        
END SUBROUTINE MYDIFF

SUBROUTINE ZOOO(o,ozoo,j,i,olost) !!written for equal area boxes!!
implicit none
REAL*8 :: o(nz+4,ny+4),ozoo,olost(nz+4)
REAL*8 :: ot
INTEGER :: j,i
INTEGER :: jc,ic,box,boxm1,hbox,nboxt,mino

nboxt=201 !201 !total number of boxes to average zoo consumption over. odd for even spacing.
mino=10e-3*1D0 !minimum oxygen conc (mol/m3) at which zooplankton breathe.

jc=j+2
ic=i+2		

ot=0
olost(:)=0

!if (ozoo.le.mino) then !!only do this if total zoo consumption is <=8uM for this timestep- should be. because if not, need to make sure there's enough O2 in each box that O2 is being taken out of. 8 is the greatest that would be (if all in one box). but actually this competes with B O2 use, too! hmm... need O2 limit on zoo growth?!? one solution: calculate derivative already and make sure that this new consumption isn't more than that change. (if ozoo.le.kdO)

hbox=(nboxt-1)/2 !bottom box here

	!1. get total sum of oxygen over the nboxes
	do box=1,nboxt
		boxm1=box-1
		if ((jc-hbox+boxm1).ge.3.and.(jc-hbox+boxm1).le.(nz+2).and.o(jc-hbox+boxm1,ic).gt.mino) then
			ot=ot+o(jc-hbox+boxm1,ic)
            !print*,ot
		end if
	end do
	
	!2. get weighted consumption: amt of o2 used by zoo, ozoo, weighted by o in each box
	do box=1,nboxt
		boxm1=box-1
		if ((jc-hbox+boxm1).ge.3.and.(jc-hbox+boxm1).le.(nz+2).and.o(jc-hbox+boxm1,ic).gt.mino) then
			olost(jc-hbox+boxm1)=o(jc-hbox+boxm1,ic)/ot*ozoo !weighted consumption of o for this box
		end if
	end do

	if (ot.eq.0) then
	
	!append files:
	OPEN(UNIT=1,FILE='oxygenproblems.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',POSITION='APPEND')
	!DO I=1,ny+4
	WRITE(1,*) 'equals 0! at: (j,i)'
	!END DO
	CLOSE(1)

	!print*,'equals 0! at: (j,i)'
	!print*,j
	!print*,i
	
	end if

!else if (ozoo.gt.8e-3*1D0) then
!	print*,'ozoo is greater than 8uM at: (j,i)'
!	print*,j
!	print*,i
!end if
        
END SUBROUTINE ZOOO


SUBROUTINE MYRK(nh4_one,no3_one,no2_one,n2o_one,d_one,dom_one,o_one, &
				zoo_one,zoo2_one,zoo3_one, &
				p1_one,p2_one,p3_one, &
				bo_one,bo_pa_one,bdno3_one,bdno2_one,bdn2o_one, &
				bnh4_one,bno2_one,bana_one, & 
				knh4_two,kno3_two,kno2_two,kn2o_two,kd_two,kdom_two,ko_two, &
				kzoo_two,kzoo2_two,kzoo3_two, &
				kp1_two,kp2_two,kp3_two, &
				kbo_two,kbo_pa_two,kbdno3_two,kbdno2_two,kbdn2o_two, &
				kbnh4_two,kbno2_two,kbana_two, &
				nh4_two,no3_two,no2_two,n2o_two,d_two,dom_two,o_two, &
				zoo_two,zoo2_two,zoo3_two, &
                p1_two,p2_two,p3_two, &
				bo_two,bo_pa_two,bdno3_two,bdno2_two,bdn2o_two, &
				bnh4_two,bno2_two,bana_two) !, & 
				!po,pd,pdom,pnh4,pno2,pno3,pn2o,u_bo,u_bo_pa,u_bdno3,u_bdno2,u_bdn2o, &
				!u_bnh4,u_bno2,u_bana,distn)
implicit none
REAL*8, dimension(nz+4,ny+4), intent(in) :: nh4_one,no3_one,no2_one,n2o_one, &
				d_one,dom_one,o_one,zoo_one,zoo2_one,zoo3_one,p1_one,p2_one,p3_one, &
				bo_one,bo_pa_one,bdno3_one,bdno2_one,bdn2o_one,bnh4_one,bno2_one,bana_one
REAL*8, dimension(nz+4,ny+4), intent(out) :: knh4_two,kno3_two,kno2_two,kn2o_two, &
				kd_two,kdom_two,ko_two,kzoo_two,kzoo2_two,kzoo3_two, &
                kp1_two,kp2_two,kp3_two, &
				kbo_two,kbo_pa_two,kbdno3_two,kbdno2_two,kbdn2o_two, &
				kbnh4_two,kbno2_two,kbana_two, &
				nh4_two,no3_two,no2_two,n2o_two,d_two,dom_two,o_two,zoo_two,zoo2_two,zoo3_two, &
                p1_two,p2_two,p3_two, &
				bo_two,bo_pa_two,bdno3_two,bdno2_two,bdn2o_two,bnh4_two,bno2_two,bana_two
REAL*8 :: ozoo !,distn -tinker with grazing routine later
				
!REAL*8, dimension(nz+4,ny+4) :: bt,po,pd,pdom,pnh4,pno2,pno3,pn2o,u_bo,u_bo_pa,u_bdno3,u_bdno2,u_bdn2o,u_bnh4,u_bno2,u_bana,ozooall


bt=bo_one+bo_pa_one+bdno3_one+bdno2_one+bdn2o_one+bnh4_one+bno2_one+bana_one
!linear grazing:
ozooall=RredO*(1-gam)*g*bt*zoo_one*TempFun
!explicity grazing:
!g=gmax*bt/(bt+kB)*TempFun
!ozooall=RredO*(1-gam)*g*zoo_one
olosszoo(:,:)=0D0
do i=1,ny
    do j=1,nz ; ic=i+2; jc=j+2;       
        call myupwind(nh4_one,v,w,j,i,dy,dz,Aij,nz,ny,adv); !%advection scheme outputs advective fluxes             
        call mydiff(nh4_one,Ky,Kz,j,i,dy,dz,Aij,nz,ny,diff)
        knh4_two(jc,ic)=-adv+diff
        call myupwind(no3_one,v,w,j,i,dy,dz,Aij,nz,ny,adv); !%advection scheme outputs advective fluxes             
        call mydiff(no3_one,Ky,Kz,j,i,dy,dz,Aij,nz,ny,diff)
        kno3_two(jc,ic)=-adv+diff
        call myupwind(no2_one,v,w,j,i,dy,dz,Aij,nz,ny,adv); !%advection scheme outputs advective fluxes             
        call mydiff(no2_one,Ky,Kz,j,i,dy,dz,Aij,nz,ny,diff)
        kno2_two(jc,ic)=-adv+diff
        !call myupwind(n2o_one,v,w,j,i,dy,dz,Aij,nz,ny,adv); !%advection scheme outputs advective fluxes             
        !call mydiff(n2o_one,Ky,Kz,j,i,dy,dz,Aij,nz,ny,diff)
        !kn2o_two(jc,ic)=-adv+diff
        !
        call myupwind(bo_one,v,w,j,i,dy,dz,Aij,nz,ny,adv); !%advection scheme outputs advective fluxes      
        call mydiff(bo_one,Ky,Kz,j,i,dy,dz,Aij,nz,ny,diff)
        kbo_two(jc,ic)=-adv+diff
        !
        call myupwind(bdno3_one,v,w,j,i,dy,dz,Aij,nz,ny,adv); !%advection scheme outputs advective fluxes      
        call mydiff(bdno3_one,Ky,Kz,j,i,dy,dz,Aij,nz,ny,diff)
        kbdno3_two(jc,ic)=-adv+diff
        call myupwind(bdno2_one,v,w,j,i,dy,dz,Aij,nz,ny,adv); !%advection scheme outputs advective fluxes      
        call mydiff(bdno2_one,Ky,Kz,j,i,dy,dz,Aij,nz,ny,diff)
        kbdno2_two(jc,ic)=-adv+diff
        !call myupwind(bdn2o_one,v,w,j,i,dy,dz,Aij,nz,ny,adv); !%advection scheme outputs advective fluxes      
        !call mydiff(bdn2o_one,Ky,Kz,j,i,dy,dz,Aij,nz,ny,diff)
        !kbdn2o_two(jc,ic)=-adv+diff
        !!
        call myupwind(bnh4_one,v,w,j,i,dy,dz,Aij,nz,ny,adv); !%advection scheme outputs advective fluxes      
        call mydiff(bnh4_one,Ky,Kz,j,i,dy,dz,Aij,nz,ny,diff)
        kbnh4_two(jc,ic)=-adv+diff
        call myupwind(bno2_one,v,w,j,i,dy,dz,Aij,nz,ny,adv); !%advection scheme outputs advective fluxes      
        call mydiff(bno2_one,Ky,Kz,j,i,dy,dz,Aij,nz,ny,diff)
        kbno2_two(jc,ic)=-adv+diff
        call myupwind(bana_one,v,w,j,i,dy,dz,Aij,nz,ny,adv); !%advection scheme outputs advective fluxes      
        call mydiff(bana_one,Ky,Kz,j,i,dy,dz,Aij,nz,ny,diff)
        kbana_two(jc,ic)=-adv+diff
        !
        call myupwind(d_one,v,wd,j,i,dy,dz,Aij,nz,ny,adv); !%!advection scheme outputs advective fluxes      
        call mydiff(d_one,Ky,Kz,j,i,dy,dz,Aij,nz,ny,diff)
        kd_two(jc,ic)=-adv+diff
        !
        call myupwind(o_one,v,w,j,i,dy,dz,Aij,nz,ny,adv); !%advection scheme outputs advective fluxes         
        call mydiff(o_one,Ky,Kz,j,i,dy,dz,Aij,nz,ny,diff)
        ko_two(jc,ic)=-adv+diff
        !
        call myupwind(zoo_one,v,w,j,i,dy,dz,Aij,nz,ny,adv); !%advection scheme outputs advective fluxes 
        call mydiff(zoo_one,Ky,Kz,j,i,dy,dz,Aij,nz,ny,diff)
        kzoo_two(jc,ic)=-adv+diff       
        ! 
        	!distributed zooplankton oxygen consumption (for B grazers only)
		ozoo=ozooall(jc,ic)
        	call ZOOO(o_one,ozoo,j,i,olost) !olost is out
        	olosszoo(:,ic)=olosszoo(:,ic) + olost !all positive quantities
        !
        call myupwind(zoo2_one,v,w,j,i,dy,dz,Aij,nz,ny,adv); !%advection scheme outputs advective fluxes 
        call mydiff(zoo2_one,Ky,Kz,j,i,dy,dz,Aij,nz,ny,diff)
        kzoo2_two(jc,ic)=-adv+diff 
        !
        call myupwind(zoo3_one,v,w,j,i,dy,dz,Aij,nz,ny,adv); !%advection scheme outputs advective fluxes 
        call mydiff(zoo3_one,Ky,Kz,j,i,dy,dz,Aij,nz,ny,diff)
        kzoo3_two(jc,ic)=-adv+diff 
        !
        !call myupwind(p1_one,v,w,j,i,dy,dz,Aij,nz,ny,adv); !%advection scheme outputs advective fluxes 
        !call mydiff(p1_one,Ky,Kz,j,i,dy,dz,Aij,nz,ny,diff)
        !kp1_two(jc,ic)=-adv+diff 
        call myupwind(p2_one,v,w,j,i,dy,dz,Aij,nz,ny,adv); !%advection scheme outputs advective fluxes 
        call mydiff(p2_one,Ky,Kz,j,i,dy,dz,Aij,nz,ny,diff)
        kp2_two(jc,ic)=-adv+diff 
        call myupwind(p3_one,v,w,j,i,dy,dz,Aij,nz,ny,adv); !%advection scheme outputs advective fluxes 
        call mydiff(p3_one,Ky,Kz,j,i,dy,dz,Aij,nz,ny,diff)
        kp3_two(jc,ic)=-adv+diff 
        
        if (dommodel.eq.1) then 
        	call myupwind(bo_pa_one,v,wd,j,i,dy,dz,Aij,nz,ny,adv); !%advection scheme outputs advective fluxes      
        	call mydiff(bo_pa_one,Ky,Kz,j,i,dy,dz,Aij,nz,ny,diff)
        	kbo_pa_two(jc,ic)=-adv+diff
        	!
        	call myupwind(dom_one,v,wd2,j,i,dy,dz,Aij,nz,ny,adv); !%!advection scheme outputs advective fluxes      
        	call mydiff(dom_one,Ky,Kz,j,i,dy,dz,Aij,nz,ny,diff)
        	kdom_two(jc,ic)=-adv+diff    	
        end if
        	
    end do
end do

!Uptake for bacteria, not pp:
po=po_coef*o_one !not really uptake- already normalized into 1/day by diving by Qmin (above)
!pd=pd_lin*d !linear uptake normalized to growth rates
pd=pd_max*(d_one/(d_one+kd)) 
pdom=pd_max*(dom_one/(dom_one+kd))
!pn=pn_max*(n_one/(n_one+kn))
!pn2o=pn_max*(n2o_one/(n2o_one+kn))*TempFun !check that this doesn't cause "double trouble" for the hets- though they should be limited by DOM, not N
limnh4b=(nh4_one/(nh4_one+knh4))!*TempFun
limno2b=(no2_one/(no2_one+knox))!*TempFun
limno3b=(no3_one/(no3_one+knox))!*TempFun

pnh4b=pnh4_max*limnh4b !check that this doesn't cause "double trouble" for the hets- though they are limited by DOM, not N
pno2b=pnox_max*limno2b !check that this doesn't cause "double trouble" for the hets- though they are limited by DOM, not N
pno3b=pnox_max*limno3b

!bacteria growth rates:
if (dommodel.eq.0) then 
	!to save changing the below pdom limit to using d instead of dom, just:
	pdom = pd
end if

!new way: having uptake of all hets as a f of BOTH: call this pdom to fit into old format. 
if (dommodel.eq.1) then 
    pdom = pd_max*((d_one + dom_one)/(d_one+dom_one+kd))
    domf = dom_one/(dom_one + d_one + 1D-38) 
end if
    
u_bo=min(po*yo_bo,pdom*yd_bo)*TempFun !free-living, like the rest of the hets
u_bo_pa=min(po*yo_bo,pd*yd_bo)*TempFun !particle assoc u controlled by POM
u_bdno3=min(pno3b*yno3_bdno3,pdom*yd_bdno3)*TempFun
u_bdno2=min(pno2b*yno2_bdno2,pdom*yd_bdno2)*TempFun
!u_bdn2o=min(pn2o*yn2o_bdn2o,pdom*yd_bdn2o)
u_bana=min(pno2b*yno2_bana,pnh4b*ynh4_bana)*TempFun
u_bnh4=min(pnh4b*ynh4_bnh4,po*yo_bnh4)*TempFun
u_bno2=min(pno2b*yno2_bno2,po*yo_bno2)*TempFun

!phytopl uptake and growth rate:
inhibnh4 = exp(-amminhib*nh4_one) !from GUD
!for 2:
limnh4=(nh4_one/(nh4_one+knh4_effp2))!*TempFun
limno2=(no2_one/(no2_one+knox_effp2))*inhibnh4!*TempFun
limno3=(no3_one/(no3_one+knox_effp2))*inhibnh4!*TempFun
!u_p2=umaxp2*TempFun*min(Iz/(Iz+kI),limnh4+limno2+limno3)!*exp(-Iz*.01D0) !multiply by 0.9 for penalty?
PCmax2 = umaxp2*TempFun*min(1D0,limnh4+limno2+limno3)

!again now for 3:
limnh4p3=(nh4_one/(nh4_one+knh4_effp3))!*TempFun
limno2p3=(no2_one/(no2_one+knox_effp3))*inhibnh4!*TempFun
limno3p3=(no3_one/(no3_one+knox_effp3))*inhibnh4!*TempFun
nlimtot=limnh4p3+limno2p3+limno3p3
PCmax3 = umaxp3*TempFun*min(1D0,nlimtot)

!Geider chlorophyll:c based growth rates:
a_Ip2 = phimax*a_chlp2*Iz*convI !mmol C/mol Ein * m2/mg chla * Ein/m2/d = mmol C/mg chla/d
a_Ip3 = phimax*a_chlp3*Iz*convI !mmol C/mol Ein * m2/mg chla * Ein/m2/d = mmol C/mg chla/d
chl2c_p1 = 0D0!max(chl2cmin, min(chl2cmax, chl2cmax/(1D0+chl2cmax*a_I/2D0/PCmax1)))
chl2c_p2 = max(chl2cmin, min(chl2cmax, chl2cmax/(1D0+chl2cmax*a_Ip2/2D0/PCmax2)))*inmask
chl2c_p3 = max(chl2cmin, min(chl2cmax, chl2cmax/(1D0+chl2cmax*a_Ip3/2D0/PCmax3)))*inmask
!chl2c = chl2c_p3
!PC1 = PCmax1*(1D0 - exp(-a_I*chl2c_p1/PCmax1))
PC2 = PCmax2*(1D0 - exp(-a_Ip2*chl2c_p2/PCmax2))
PC2(1:2,:)=0D0
PC2(nz+3:nz+4,:)=0D0
PC2(:,1:2)=0D0
PC2(:,ny+3:ny+4)=0D0
PC3 = PCmax3*(1D0 - exp(-a_Ip3*chl2c_p3/PCmax3))
PC3(1:2,:)=0D0
PC3(nz+3:nz+4,:)=0D0
PC3(:,1:2)=0D0
PC3(:,ny+3:ny+4)=0D0

!u_p1=PC1
u_p2=PC2
u_p3=PC3


!NO2 excretion:
!uptake of NO3 does not have f factor:
        no3uptakeP=u_p3nolim*p3_one*limno3/(nlimtot+1D-38)
    !amount of NO2 emitted is then: 
        no2emitP=0D0!(u_p3nolim-u_p3)*p3_one*limno3/(nlimtot+1D-38)

!totals for grazing and mortality
!above for grazing: bt=bo_one+bo_pa_one+bdno3_one+bdno2_one+bdn2o_one+bnh4_one+bno2_one+bana_one
pt=p1_one+p2_one+p3_one
btsq=bo_one*bo_one+bo_pa_one*bo_pa_one+bdno3_one*bdno3_one+bdno2_one*bdno2_one &
        +bdn2o_one*bdn2o_one+bnh4_one*bnh4_one+bno2_one*bno2_one+bana_one*bana_one
ptsq=p1_one*p1_one+p2_one*p2_one+p3_one*p3_one
!g_p=gmax_p*pt/(pt+kP)*TempFun

!no bdn2o:
!distn= ( 1/yn2o_bdn2o*sum(u_bdn2o*bdn2o) + (1/ynh4_bana + 1/yno2_bana -1)*sum(u_bana*bana) &
!       - sum(koverh*(0D0-n2o_one)*eqmask) )/nz/ny !nfix represents loss of n2o to the atmosphere 
!Until Jul 2019:
!distn= ( 1D0/yno2_bdno2*sum(u_bdno2*bdno2_one) + (1D0/ynh4_bana + 1D0/yno2_bana -1D0)*sum(u_bana*bana_one) )/nz/ny
!Jul 2019 new anx stoich:
distn= ( 1D0/yno2_bdno2*sum(u_bdno2*bdno2_one) + en2_bana*sum(u_bana*bana_one) )/nz/ny

!!Oxygen limitation for pro grazer zoo3 (and also zoo2) 
O2limZ = o_one/(o_one + kO2) 

knh4_two = knh4_two &
			+ enh4_bo*(u_bo*bo_one + u_bo_pa*bo_pa_one) & !het B
			+ (1D0-gam)*g*(bt)*zoo_one*TempFun  & !zoo
			+ (1D0-gam)*g*p2_one*zoo2_one*O2limZ*TempFun  & !zoo2
			+ (1D0-gam)*g*p3_one*zoo3_one*O2limZ*TempFun  & !zoo3
			+ enh4_bdno3*u_bdno3*bdno3_one & !source: 3 anaerobic hets
			+ enh4_bdno2*u_bdno2*bdno2_one &
			!+ enh4_bdn2o*u_bdn2o*bdn2o_one  & 
			- 1D0/ynh4_bnh4*u_bnh4*bnh4_one & !NH4 oxidizer consumption
			- 1D0/ynh4_bana*u_bana*bana_one & !anammox consumption
			!- u_p1*p1_one & !phytoplankton
            !- u_p2*p2_one*limnh4/(limnh4+limno2+1D-38) &
            - u_p2*p2_one*limnh4/(limnh4+limno2+limno3+1D-38) &
            - u_p3*p3_one*limnh4p3/(limnh4p3+limno2p3+limno3p3+1D-38)
            !distn? but NO3 seems better!

kno2_two = kno2_two &
			+ eno2_bnh4*u_bnh4*bnh4_one & !source: NH4 oxidizer
			+ eno2_bdno3*u_bdno3*bdno3_one & !source: het NO3 denitrifier
			- 1D0/yno2_bno2*u_bno2*bno2_one & !sink: aerobic NO2 oxidizer,
			- 1D0/yno2_bdno2*u_bdno2*bdno2_one & !sink: het NO2 reducer
			- 1D0/yno2_bana*u_bana*bana_one &  !sink: anammox
			!- u_p2*p2_one*limno2/(limnh4+limno2+1D-38) & !phytoplankton!
            - u_p2*p2_one*limno2/(limnh4+limno2+limno3+1D-38) & !or no lim?
            - u_p3*p3_one*limno2p3/(limnh4p3+limno2p3+limno3p3+1D-38)
            !+ no2emitP !phytoplankton leak
		

kno3_two = kno3_two &
			+ eno3_bno2*u_bno2*bno2_one  & !source: aerobic NO2 oxidizer
	   		- 1D0/yno3_bdno3*u_bdno3*bdno3_one & !sink: het NO3 denitrifier
            + en2_bana*u_bana*bana_one & !source: anammox    
	   		- u_p2*p2_one*limno3/(limnh4+limno2+limno3+1D-38) & !no3uptakeP
	   		- u_p3*p3_one*limno3p3/(limnh4p3+limno2p3+limno3p3+1D-38) &
	   		+  distn*inmask  !representing far away N fixation, assume it enters in laterally...
	   		
	   		
!kn2o_two = kn2o_two + en2o_bdno2*u_bdno2*bdno2_one - 1/yn2o_bdn2o*u_bdn2o*bdn2o_one &
!            + koverh*(0D0-n2o_one)*eqmask !koverh*(satconc-conc) source: het NO2 denitrifier, sink: het N2O denitrifier	and air-sea flux   		
				
kbo_two= kbo_two +  bo_one*(u_bo - mb*TempFun - mquad*bo_one*TempFun - g*zoo_one*TempFun)!/(bt+1D-38)) !grazing should have a different Q10!!! 
kbo_pa_two= kbo_pa_two +  bo_pa_one*(u_bo_pa - mb*TempFun - mquad*bo_pa_one*TempFun - g*zoo_one*TempFun)!/(bt+1D-38))
kbdno3_two= kbdno3_two +  bdno3_one*(u_bdno3 - mb*TempFun - mquad*bdno3_one*TempFun - g*zoo_one*TempFun)!/(bt+1D-38)) 
kbdno2_two= kbdno2_two +  bdno2_one*(u_bdno2 - mb*TempFun - mquad*bdno2_one*TempFun - g*zoo_one*TempFun)!/(bt+1D-38)) 
kbdn2o_two= kbdn2o_two +  bdn2o_one*(u_bdn2o - mb*TempFun - mquad*bdn2o_one*TempFun - g*zoo_one*TempFun)!/(bt+1D-38)) 
kbnh4_two= kbnh4_two +  bnh4_one*(u_bnh4 - mb*TempFun - mquad*bnh4_one*TempFun - g*zoo_one*TempFun)!/(bt+1D-38)) 
kbno2_two= kbno2_two +  bno2_one*(u_bno2 - mb*TempFun - mquad*bno2_one*TempFun - g*zoo_one*TempFun)!/(bt+1D-38)) 
kbana_two= kbana_two +  bana_one*(u_bana - mb*TempFun - mquad*bana_one*TempFun - g*zoo_one*TempFun)!/(bt+1D-38)) 

!kp1_two = kp1_two + p1_one*(u_p1 - mp - g*zoo2_one)
kp2_two = kp2_two + p2_one*(u_p2 - mp*TempFun - mquad*p2_one*TempFun - g*zoo2_one*O2limZ*TempFun)!/(pt+1D-38))
kp3_two = kp3_two + p3_one*(u_p3 - mp*TempFun - mquad*p3_one*TempFun - g*zoo3_one*O2limZ*TempFun)!/(pt+1D-38))

kzoo_two=kzoo_two + gam*g*(bt)*zoo_one*TempFun - mz*zoo_one*zoo_one*TempFun         
!kzoo2_two=kzoo2_two + gam*g*(p1_one+p2_one+p3_one)*zoo2_one - mz*zoo2_one*zoo2_one 
kzoo2_two=kzoo2_two + gam*g*p2_one*zoo2_one*O2limZ*TempFun - mz*zoo2_one*zoo2_one*TempFun 
kzoo3_two=kzoo3_two + gam*g*p3_one*zoo3_one*O2limZ*TempFun - mz*zoo3_one*zoo3_one*TempFun 

ko_two= ko_two &
		+ RredO*(u_p1*p1_one + u_p2*p2_one + u_p3*p3_one) & !source: PP
		+ koverh*(o2sat-o_one)*eqmask &  !source: air-sea
		- 1D0/yo_bo*u_bo*bo_one & !sink: 2 aerobic hets
		- 1D0/yo_bo*u_bo_pa*bo_pa_one &
		- 1D0/yo_bnh4*u_bnh4*bnh4_one & !sink: 2 aerobic N oxidizers
		- 1D0/yo_bno2*u_bno2*bno2_one &
		- olosszoo & !zoo1 (deep B eaters)
		- RredO*(1D0-gam)*g*p2_one*zoo2_one*O2limZ*TempFun & !zoo2 (zoo is already accounted for above)
		- RredO*(1D0-gam)*g*p3_one*zoo3_one*O2limZ*TempFun !zoo2 (zoo is already accounted for above)
       ! - RredO*(1-gam)*g*pt*zoo_one*TempFun !check that this is the same as doing all O2 above

if (dommodel.eq.1) then !only PA hets eat d

	kd_two= kd_two &
			+ (1D0-mortf)*mb*bt*TempFun &
			+ (1D0-mortf)*mz*zoo_one*zoo_one*TempFun & !Q10 for grazing!
			+ (1D0-mortf)*mz*zoo2_one*zoo2_one*TempFun & !Q10 for grazing!
			+ (1D0-mortf)*mz*zoo3_one*zoo3_one*TempFun & !Q10 for grazing!
			+ (1D0-mortf)*mp*pt*TempFun & !Q10 here too? ask steph
			+ (1D0-mortf)*mquad*ptsq*TempFun & !Q10 here too? ask steph
			+ (1D0-mortf)*mquad*btsq*TempFun & !Q10 here too? ask steph
            - 1D0/yd_bo*u_bo*bo_one*(1D0 - domf) & !sink: all heterotrophs
            - 1D0/yd_bdno3*u_bdno3*bdno3_one*(1D0 - domf) &
            - 1D0/yd_bdno2*u_bdno2*bdno2_one*(1D0 - domf) 
			!- alpha*u_bo_pa*bo_pa_one/yd_bo  !sink: only the PA bacteria!
		
	kdom_two= kdom_two &
			!+ (alpha-1)*1/yd_bo*u_bo_pa*bo_pa_one & !source: the PA bacteria
			+ mortf*mb*bt*TempFun &
			+ mortf*mz*zoo_one*zoo_one*TempFun & !Q10 for grazing!
			+ mortf*mz*zoo2_one*zoo2_one*TempFun & !Q10 for grazing!
			+ mortf*mz*zoo3_one*zoo3_one*TempFun & !Q10 for grazing!
			+ mortf*mp*pt*TempFun & !Q10 here too? ask steph
			+ mortf*mquad*ptsq*TempFun & !Q10 here too? ask steph
			+ mortf*mquad*btsq*TempFun & !Q10 here too? ask steph
            - 1D0/yd_bo*u_bo*bo_one*domf & !sink: all heterotrophs
            - 1D0/yd_bdno3*u_bdno3*bdno3_one*domf &
            - 1D0/yd_bdno2*u_bdno2*bdno2_one*domf 
			!- 1/yd_bdn2o*u_bdn2o*bdn2o_one 
			
else !all d and dom together in d: 
	
	kd_two= kd_two &
			+ mb*bt*TempFun &
			+ mz*zoo_one*zoo_one*TempFun & !Q10 for grazing!
			+ mz*zoo2_one*zoo2_one*TempFun & !Q10 for grazing!
			+ mz*zoo3_one*zoo3_one*TempFun & !Q10 for grazing!
			+ mp*pt*TempFun & !Q10 here too? ask steph
			+ mquad*ptsq*TempFun & !Q10 here too? ask steph
			+ mquad*btsq*TempFun & !Q10 here too? ask steph
			- 1D0/yd_bo*u_bo*bo_one & !sink: all 4 heterotrophs
			- 1D0/yd_bdno3*u_bdno3*bdno3_one &
			- 1D0/yd_bdno2*u_bdno2*bdno2_one 
			!- 1D0/yd_bdn2o*u_bdn2o*bdn2o_one 
			
end if


!TIMESTEP:
!i assume that bc i checked this out earlier, when i wrote this, and made sure numbers were exact with the non-RK model, that the following works, even though i didn't specify bo, d, dom in the subroutine:
nh4_two = nh4 + dt/2D0*knh4_two; 
no2_two = no2 + dt/2D0*kno2_two; 
no3_two = no3 + dt/2D0*kno3_two; 
n2o_two = n2o + dt/2D0*kn2o_two; 
bo_two = bo + dt/2D0*kbo_two; 
bo_pa_two = bo_pa + dt/2D0*kbo_pa_two; 
bdno3_two = bdno3 + dt/2D0*kbdno3_two; 
bdno2_two = bdno2 + dt/2D0*kbdno2_two; 
!bdn2o_two = bdn2o + dt/2*kbdn2o_two; 
bana_two = bana + dt/2D0*kbana_two; 
bnh4_two = bnh4 + dt/2D0*kbnh4_two; 
bno2_two = bno2 + dt/2D0*kbno2_two; 
d_two = d + dt/2D0*kd_two; 
dom_two = dom + dt/2D0*kdom_two; 
o_two = o + dt/2D0*ko_two; 
zoo_two = zoo + dt/2D0*kzoo_two; 
zoo2_two = zoo2 + dt/2D0*kzoo2_two; 
zoo3_two = zoo3 + dt/2D0*kzoo3_two; 
!p1_two = p1 + dt/2*kp1_two; 
p2_two = p2 + dt/2D0*kp2_two; 
p3_two = p3 + dt/2D0*kp3_two; 

        
END SUBROUTINE MYRK

END PROGRAM EZM
