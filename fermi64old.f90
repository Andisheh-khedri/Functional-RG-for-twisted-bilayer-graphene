!-----------------------------------------------------------------------------------------------------
! *****************--compared to fermi63.f90 trying to employ mpi to distribute jobs
!-----------------------------------------------------------------------------------------------------
module commons
INTEGER, PARAMETER::nspin=2, norbital=2, ntype=2, nband=nspin*norbital*ntype&  !
                   ,npatch=6*2,Ntot=(nband*((nband)**3))*(npatch**3)&
                   ,nprotpatch=6,nprot=npatch/nprotpatch, Ntotal=((npatch**3)*(nband**4))&
                   ,idpermute=0,IdFSrotsym=1,Idantisym=1,Idantisymextra=0,ISOC=0,Idbandsym=1
DOUBLE PRECISION,PARAMETER::t=1.0d+00,u=4.0d+00*t,beta=1.0d+00/((10.0d+00**(-6.0d+00))*t)&
                            ,pi=atan2(0.0d+00,-t),Cfield=0.0d+00*t&
                            ,tprime=(0.0d+00)*t,tdprime=(0.0d+00)*t&
                            ,Jcoup=-0.0d+00*u,uprime=0.0d+00*u&!u-1.0d+00*Jcoup &!3.0d+00*t&!
                            ,Jprimespin=0.0d+00*Jcoup,Jprimeorb=0.0d+00*Jcoup&  !0.0d+00,Jprime=0.0d+00!
                            ,lamSO=0.0d+00*t,Bfieldx=((sqrt(2.0d+00))/2.0d+00)*t*0.0d+00,Bfieldy=Bfieldx,Bfieldz=0.0d+00&
                            ,Lconst=1.0d+00,BZconst=4.0d+00*pi/(3.0d+00*Lconst*sqrt(3.0d+00))&
                            ,Vint=0.0d+00*t
CHARACTER*3,PARAMETER::patchflag='sim' !options: 'sim' or 'com'
CHARACTER*3,PARAMETER::angleflag='eqd' ! options: 'eqd' or 'sop'
INTEGER::icall
COMPLEX*16, PARAMETER::ic=cmplx(0.0d+00,1.0d+00),VSOx12=lamSO,VSOx21=-lamSO,VSOy12=-ic*LamSO,VSOy21=-ic*LamSO

! Technical parameters:

DOUBLE PRECISION,PARAMETER::patcherr1=t*(10.0d+00**(-15.0d+00)),patcherr2=t*(10.0d+00**(-14.0d+00))&
                           ,zeroHex=t*(10.0d+00**(-12.0d+00)),zeroHexedge=t*(10.0d+00**(-12.0d+00))&
                           ,edgeboarder=10.0d+00**(-10.0d+00)

!-----------------------------------------------------------------------------------------------------
!------------------------------------------------------------------
end module
!------------------------------------------------------------------
!------------------------------------------------------------------
program HFRG

use commons

implicit none

common eF,maxvertex,imat,start_time,irout,p4int,thetap,KFSpmat,imax!,Yvertex,imax
!------------------------------------------------------------------
!  Parameter declaration:
!------------------------------------------------------------------
Integer,parameter::nangle=3001,nrot=4
Double precision::density,eF,angle,zero,scalecommon,k1,kzerop,kmax&
                 ,kFS,resol,maxvertex,groupvelx,groupvely,base&
                 ,kpx,kpy,angle1,angle2,angle3,KFS1,KFS2,KFS3,anglek,anglered
DOUBLE PRECISION,DIMENSION(1)::Yang,YDang
Integer::eFint,angleint,k1int,nroots,NEQ,n1,n2,n3,n4,p1,p2,p3,i&
         ,icount,irout,k,p4v,kxint,kyint,n1spin,n2spin,n3spin,n4spin&
         ,n1orb,n2orb,n3orb,n4orb,n1orbminus,l1,l2,nroots1,nroots2,kint&
         ,n1type,n2type,n3type,n4type,n1p,n2p,n3p,n4p,imax,p3max,kk

DOUBLE PRECISION,dimension(1:nband,1:npatch)::thetap,KFSpmat
INTEGER,DIMENSION(1:nband,1:npatch)::nrootspmat

!------------------------------------------------------------------
! parameters for the ode-solver:
!------------------------------------------------------------------

INTEGER:: openstatus,Flag, ISTATE, IOPT, MF,MU, ML, LRW, LIW
DOUBLE PRECISION::y,yout, RTOL, ATOL,JAC

DOUBLE PRECISION,allocatable,dimension(:)::RWORK
INTEGER,allocatable,dimension(:)::IWORK
!------------------------------------------------------------------
! External functiions
!------------------------------------------------------------------
DOUBLE PRECISION:: fermi,dispersion,groupvel,numgroupvel
! their corresponding variables:
DOUBLE PRECISION::mom,theta
INTEGER,DIMENSION(1:Ntotal):: n1mat,n2mat,n3mat,n4mat,p1mat,p2mat,p3mat
INTEGER,DIMENSION(1:nband,1:nband,1:nband,1:nband,1:npatch,1:npatch,1:npatch)::jmat
INTEGER,DIMENSION(5:8,5:8,5:8,5:8,1:npatch,1:npatch,1:npatch)::imat
INTEGER,DIMENSION(1:nband,1:nband,1:nband,1:npatch,1:npatch,1:npatch)::p4int,nfold,p4nf
DOUBLE PRECISION::kdelta
!------------------------------------------------------------------
! parameters for zero_rc library
!------------------------------------------------------------------
DOUBLE PRECISION::zeroer,zeromin,zeromax,zeroarg,zerovalue,zerofun
INTEGER*4::zerostatus
!------------------------------------------------------------------
! Unnecessary (perhaps!):
!------------------------------------------------------------------
INTEGER,parameter::ITOL=1,nrowpd=2,NEQ1=1
DOUBLE PRECISION,DIMENSION(nrowpd,NEQ1)::pd1
!------------------------------------------------------------------
! Actual flow parameters:
!------------------------------------------------------------------
DOUBLE PRECISION::x,xout
!------------------------------------------------------------------
DOUBLE PRECISION,DIMENSION(1:2*Ntotal)::Yvertexcommon,Yvertex,YDvertex
DOUBLE PRECISION,DIMENSION(1:2*Ntot)::Yvertexred,YDvertexred
DOUBLE PRECISION,DIMENSION(1:Ntot)::Yvertexre,YDvertexre
INTEGER::nchcommon,icommon
!------------------------------------------------------------------
INTEGER,DIMENSION(1:npatch,1:npatch,1:npatch)::nrep,p4extra,nmean
INTEGER,DIMENSION(1:npatch,1:npatch,1:npatch,1:npatch)::p3mean
INTEGER,DIMENSION(1:npatch,1:npatch)::nmis
INTEGER,DIMENSION(1:nband)::nspinmat,norbmat,ntypemat
!------------------------------------------------------------------
! spin to band transformation:
!------------------------------------------------------------------
DOUBLE PRECISION::VSOx,VSOy,VSOz,ESOp,ESOn
DOUBLE PRECISION,DIMENSION(1:nband)::Ematp4,Emateisp4
INTEGER::s1lab,s2lab,s3lab,s4lab,iprime,s1plab,s2plab,s3plab,s4plab,idprime
COMPLEX*16::VSO12,VSO21
DOUBLE PRECISION,DIMENSION(1:nband,1:nband,1:nband,1:npatch,1:npatch,1:npatch,1:2)::p4vec
COMPLEX*16::summat
COMPLEX,DIMENSION(1:Ntotal)::Yvertexcomplex,YvertexAB
DOUBLE PRECISION,DIMENSION(1:Ntotal)::YvertexABinitial
INTEGER,DIMENSION(1:npatch,1:nband)::i1dmat
INTEGER,DIMENSION(1:npatch,1:npatch,1:nband,1:nband)::i2dmat
INTEGER,DIMENSION(1:npatch,1:npatch,1:npatch,1:nband,1:nband,1:nband)::i3dmat
INTEGER,DIMENSION(1:nband,1:nband,1:nband,1:nband,1:npatch,1:npatch,1:npatch)::i4dmat
!------------------------------------------------------------------
! parameters for the Eispack library:
!------------------------------------------------------------------
INTEGER*4::Mdim,EVerr
REAL*8,ALLOCATABLE,DIMENSION(:,:)::Mmatre,Mmatim,MEVecre,MEVecim
REAL*8,ALLOCATABLE,DIMENSION(:)::MEValre,MEValim
LOGICAL::EVflag
!------------------------------------------------------------------
! checking!
!------------------------------------------------------------------
DOUBLE PRECISION::Hexdisp,kpxr,kpyr,kpxri,kpyri,KFSr,angler&
                  ,dispbefore,KFSmax,kpxrnew,kpyrnew
INTEGER::BZcount,BZid,nkxint,nkyint,neiint,loopcount,nneimax&
        ,neicount,numsign,nrootsmat,j
DOUBLE PRECISION,DIMENSION(1:6,1:2)::BZV 
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::BZVnei,BZVneibe
DOUBLE PRECISION,DIMENSION(1:20)::KFSext,KFSextbe,KFSmat
!------------------------------------------------------------------
! Going from type-basis to diag-basis 
!------------------------------------------------------------------
COMPLEX*16,DIMENSION(1:2,1:2,1:npatch)::Umattype,Mmattype
COMPLEX*16,DIMENSION(1:npatch)::fun
DOUBLE PRECISION,DIMENSION(1:npatch)::phase
COMPLEX*16,DIMENSION(1:2,1:2)::Umatp4
DOUBLE PRECISION::Yverdummy,YverdummyV,Yverdummyu
COMPLEX*16::funcomplex
INTEGER,DIMENSION(1:ntype,1:norbital,1:nspin)::nbandmat
!INTEGER,DIMENSION(1:norbital,1:nspin)::nbandmat
INTEGER::n1prime,n2prime,n3prime,n4prime
COMPLEX,DIMENSION(1:ntype,1:ntype,1:ntype,1:ntype)::sumcheck,sumcheckback
INTEGER,DIMENSION(1:nband,1:nband,1:nband,1:nband,1:npatch,1:npatch,1:npatch)::iredmat
INTEGER,DIMENSION(1:Ntotal)::imattoiredmat
INTEGER,DIMENSION(1:Ntot)::iredmattoimat
!------------------------------------------------------------------
DOUBLE PRECISION::thetamin,thetamax,thetabefore,KFSthetamax,energy,kxmax,kymax
DOUBLE PRECISION,DIMENSION(1:npatch)::anglematmin,anglematmax,anglemat
INTEGER::idfold
character*6::nbands
character*4::phis
!------------------------------------------------------------------
!------------------------------------------------------------------
real :: start_time, stop_time
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
call cpu_time(start_time)
!-----------------------------------------------------------------------
print*,'BZconst=',BZconst
!print*,'finding the error:',LRWc,20+16*Ntot
!-----------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------
open(unit=1,file="Aug1testa1r1b1np12free.txt",status="new",action="write",position="rewind",IOSTAT=openstatus)
!-----------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------
!open(unit=2,file="Aug1test2np12.txt",status="new",action="write",position="rewind",IOSTAT=openstatus)
!-----------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------
!open(unit=3,file="Aug1test3np12.txt",status="new",action="write",position="rewind",IOSTAT=openstatus)
!-----------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------
!open(unit=41,file="PHI2",status="old",action="read",position="rewind",IOSTAT=openstatus)
!-----------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------
!open(unit=3,file="testmatsum57_500_i2.txt",status="new",action="write",position="rewind",IOSTAT=openstatus)
!----------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------
!open(unit=4,file="Aug1test4np12.txt",status="new",action="write",position="rewind",IOSTAT=openstatus)
!-----------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------
!open(unit=5,file="Aug1test5np12.txt",status="new",action="write",position="rewind",IOSTAT=openstatus)
!-----------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------
print*,'int-check:',anint(0.0d+00),anint(-0.49d+00),anint(-0.50d+00),anint(-0.99d+00),anint(-1.0d+00),anint(-1.1d+00)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Technical parameters:
!----------------------------------------------------------------------
kzerop=t*((10.0d+00)**(-5.0d+00))
kmax=(sqrt(2.0d+00))*pi
!----------------------------------------------------------------------
!----------------------------------------------------------------------
JAC=0.0d+00
print*,((nband**4)*(npatch**3)),Ntot
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!if(pi<0)then
!----------------------------------------------------------------------
! Plotting the BZ for Hexagonal lattice!
!----------------------------------------------------------------------
Do kyint=1,2
Do kxint=1,101

kpx=1.0d+00*BZconst*(-1.0d+00+2.0d+00*((kxint-1)/100.0d+00))   

if(kpx .LE. BZconst .and. kpx> (BZconst/2.0d+00))then
kpy=(sqrt(3.0d+00))*(BZconst-kpx)
elseif(kpx .LE. (BZconst/2.0d+00) .and. kpx>(-BZconst/2.0d+00))then
kpy=((sqrt(3.0d+00))/2.0d+00)*BZconst
elseif(kpx .LE. (-BZconst/2.0d+00) .and. kpx>(-BZconst))then
kpy=(sqrt(3.0d+00))*(kpx+BZconst)
endif

if(kyint==1)then
kpy=kpy
elseif(kyint==2)then
kpy=-kpy
endif

CALL Folding(kpx,kpy,kpx,kpy)

CALL Hexagonal(kpx,kpy,BZid)

if(BZid==0)then
print*,'Outside of BZ!'
endif

!write(5,*)kpx/BZconst,kpy/BZconst

end do !kxint
end do !kyint
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
Do eFint=1,1

!----------------------------------------------------------------------
! Fermi energy:
!----------------------------------------------------------------------
eF=-t*1.0d+00
!----------------------------------------------------------------------
!----------------------------------------------------------------------
nrootspmat(1:nband,1:npatch)=0
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Plotting dispersion:
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! plotting dispersion:

Do l1=1,5
Do n1=1,3
!----------------------------------------------------------------------
Do kint=1,101
!----------------------------------------------------------------------
if(n1==1)then      !Gamma -> M

kpx=0.0d+00
kpy=((sqrt(3.0d+00))/2.0d+00)*BZconst*(0.0d+00+1.0d+00*((kint-1)/100.0d+00)) 
base=0.0d+00+kpy

elseif(n1==2)then  !M -> K

kpy=((sqrt(3.0d+00))/2.0d+00)*BZconst
kpx=(0.5d+00)*BZconst*(0.0d+00+1.0d+00*((kint-1)/100.0d+00)) 
base=((sqrt(3.0d+00))/2.0d+00)*BZconst+kpx

elseif(n1==3)then

kpx=(0.5d+00)*BZconst*(1.0d+00-1.0d+00*((kint-1)/100.0d+00)) 
kpy=(sqrt(3.0d+00))*kpx
base=((sqrt(3.0d+00))/2.0d+00)*BZconst+(BZconst/2.0d+00)+(BZconst-KFS)
endif
!----------------------------------------------------------------------
KFS=sqrt((kpx**2)+(kpy**2))
angle=atan2(kpy,kpx)
angle=(-SIGN(1.0d+00,angle)+1.0d+00)*(pi)+angle

!write(2,*)kpx/BZconst,kpy/BZconst

if(l1==1)then
!write(2,*)kpx/BZconst,kpy/BZconst,(base)/BZconst,dispersion(l1,KFS,angle,0.0d+00)
elseif(l1==5)then
!write(3,*)kpx/BZconst,kpy/BZconst,(base)/BZconst,dispersion(l1,KFS,angle,0.0d+00)
endif


end do !kint
end do !n1
end do !l1
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Discretization the mom-dep:
!----------------------------------------------------------------------
Do k=1,npatch

angle=(2.0d+00*pi/npatch)*(k-0.5d+00)!(k-0.5d+00+(10.0d+00**(-3.0d+00)))

!----------------------------------------------------------------------
Do l1=1,nband
!----------------------------------------------------------------------

CALL roots_dispersion(l1,angle,eF,nroots,KFS)


if(nroots<1)then

CALL BZ_edge(angle,KFS)

nroots=1

endif

kpx=KFS*cos(angle)
kpy=KFS*sin(angle)

CALL Hexagonal(kpx,kpy,BZid)

if(BZid==0 .and. nroots==1)then

CALL Folding(kpx,kpy,kpx,kpy)

KFS=sqrt((kpx**2)+(kpy**2))

anglek=atan2(kpy,kpx)
anglek=(-SIGN(1.0d+00,anglek)+1.0d+00)*(pi)+anglek

else

anglek=angle

endif

kpx=KFS*cos(anglek)
kpy=KFS*sin(anglek)

CALL Hexagonal(kpx,kpy,BZid)

if(BZid==0)then
print*,'Careful! OUT!!'
endif

!----------------------------------------------------------------------
! finding the patch index associated to (angle,KFS):
!----------------------------------------------------------------------

CALL patch_lab_Ge(patchflag,anglek,KFS,p1,anglered)

!----------------------------------------------------------------------
if(l1==1 .or. l1==5)then
!if(k==5)then
print*,'roots-step4:',l1,k,angle/(2.0d+00*pi),k,nroots,KFS,p1,BZid
!endif
endif
!----------------------------------------------------------------------
if(p1<1 .or. p1>npatch)then
print*,'patch-lab-not.working-properly'
endif
!----------------------------------------------------------------------

KFS=sqrt(((kpx)**2)+((kpy)**2))

thetap(l1,p1)=anglek
KFSpmat(l1,p1)=KFS
nrootspmat(l1,p1)=nroots

!----------------------------------------------------------------------
!----------------------------------------------------------------------
end do !l1
!----------------------------------------------------------------------
end do !k
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
Do l1=1,nband
Do k=1,npatch

anglek=thetap(l1,k)
KFS=KFSpmat(l1,k)

kpx=KFS*cos(anglek)
kpy=KFS*sin(anglek)

!write(4,*)l1,k,kpx/BZconst,kpy/BZconst,dispersion(l1,KFS,anglek,eF),thetap(l1,k)/(2.0d+00*pi)

end do !k
end do !l1
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!if (pi<0.0d+00)then
!----------------------------------------------------------------------
Do k=1,npatch

Umattype(1,1,k)=(1.0d+00/sqrt(2.0d+00))*cmplx(-1.0d+00,0.0d+00)
Umattype(1,2,k)=(1.0d+00/sqrt(2.0d+00))*cmplx(1.0d+00,0.0d+00)
Umattype(2,1,k)=Umattype(1,2,k)
Umattype(2,2,k)=Umattype(1,2,k)

end do !k
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Vec-add on the Fermi Surface!
!----------------------------------------------------------------------
Do n1=1,nband

Do p1=1,npatch

angle1=thetap(n1,p1)
KFS1=KFSpmat(n1,p1)

Do n2=1,nband

Do p2=1,npatch

angle2=thetap(n2,p2)
KFS2=KFSpmat(n2,p2)

Do n3=1,nband

Do p3=1,npatch

angle3=thetap(n3,p3)
KFS3=KFSpmat(n3,p3)

kpx=KFS1*(cos(angle1))+KFS2*(cos(angle2))-KFS3*(cos(angle3))
kpy=KFS1*(sin(angle1))+KFS2*(sin(angle2))-KFS3*(sin(angle3))


CALL Hexagonal(kpx,kpy,BZid)

nfold(n1,n2,n3,p1,p2,p3)=BZid

if(BZid==0)then

CALL Folding(kpx,kpy,kpxr,kpyr)

CALL New_Folding(kpx,kpy,kpxrnew,kpyrnew,idfold)

else

kpxr=kpx
kpyr=kpy

kpxrnew=kpx
kpyrnew=kpy

endif

CALL Hexagonal(kpxr,kpyr,BZid)

if(BZid==0)then
print*,'careful p4 not on the FS!',kpxr/BZconst,kpyr/BZconst,p1,p2,p3,n1,n2,n3
end if

anglek=atan2(kpyr,kpxr)
anglek=(-SIGN(1.0d+00,anglek)+1.0d+00)*(pi)+anglek


p4vec(n1,n2,n3,p1,p2,p3,1)=kpxr!new
p4vec(n1,n2,n3,p1,p2,p3,2)=kpyr!new


if(sqrt((kpxr)**2+(kpyr)**2).LE. 10.0d+00**(-4.0d+00))then

p4int(n1,n2,n3,p1,p2,p3)=0

else

!------------------------------------------------------------------
! Finding the corresponding patch index:
!------------------------------------------------------------------
!Do n4=1,nband

CALL patch_lab_Ge(patchflag,anglek,sqrt(((kpxr)**2)+((kpyr)**2)),k,anglered)
p4int(n1,n2,n3,p1,p2,p3)=k

if(k<1 .or. k>npatch)then
print*,'patch-lab-Ge-not-working-properly-here is p4'
endif

!end do !n4
!------------------------------------------------------------------

endif  !p4 small

!------------------------------------------------------------------
!------------------------------------------------------------------

end do ! p3
end do ! n3
end do ! p2
end do ! n2
end do ! p1
end do ! n1
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
Yverdummyu=cmplx(0.0d+00,0.0d+00)
Yverdummyv=cmplx(0.0d+00,0.0d+00)
Yverdummy=cmplx(0.0d+00,0.0d+00)
!-----------------------------------------------------------------------
n1=0

Do n1type=1,ntype
Do n1orb=1,norbital
Do n1spin=1,nspin
n1=n1+1

norbmat(n1)=n1orb
nspinmat(n1)=n1spin
ntypemat(n1)=n1type

nbandmat(n1type,n1orb,n1spin)=n1

!nbandmat(n1orb,n1spin)=n1

end do !n1spin
end do !n1orb
end do !n1type

!-----------------------------------------------------------------------
NEQ=Ntot

Yvertexcomplex(1:Ntotal)=cmplx(0.0d+00,0.0d+00)
YvertexABinitial(1:Ntotal)=0.0d+00
YvertexAB(1:Ntotal)=cmplx(0.0d+00,0.0d+00)

Yvertex(1:2*Ntotal)=0.0d+00
Yvertexre(1:NEQ)=0.0d+00
!-----------------------------------------------------------------------

i=0
k=0

Do n1=1,nband
Do n2=1,nband
Do n3=1,nband
Do n4=1,nband

Do p1=1,npatch
Do p2=1,npatch
Do p3=1,npatch


i=i+1


jmat(n1,n2,n3,n4,p1,p2,p3)=i

n1mat(i)=n1
n2mat(i)=n2
n3mat(i)=n3
n4mat(i)=n4

p1mat(i)=p1
p2mat(i)=p2
p3mat(i)=p3

!if(n1>4 .and. n2>4 .and. n3>4 .and. n4>4)then

!k=k+1

!imat(n1,n2,n3,n4,p1,p2,p3)=k

!endif


end do
end do
end do

end do
end do
end do
end do

k=0

Do n1=5,8!1,nband
Do n2=5,8!1,nband
Do n3=5,8!1,nband
Do n4=5,8!1,nband

Do p1=1,npatch
Do p2=1,npatch
Do p3=1,npatch


k=k+1


imat(n1,n2,n3,n4,p1,p2,p3)=k


end do
end do
end do

end do
end do
end do
end do

!------------------------------------------------------------------
!------------------------------------------------------------------
! Initial Conditions:
!------------------------------------------------------------------
!------------------------------------------------------------------
Do p1=1,npatch
Do p2=1,npatch
Do p3=1,npatch

Do n1orb=1,norbital
Do n2orb=1,norbital
Do n3orb=1,norbital
Do n4orb=1,norbital

Do n1spin=1,nspin
Do n2spin=1,nspin
Do n3spin=1,nspin
Do n4spin=1,nspin

Do n1type=1,ntype
Do n2type=1,ntype
Do n3type=1,ntype
Do n4type=1,ntype


n1=nbandmat(n1type,n1orb,n1spin)
n2=nbandmat(n2type,n2orb,n2spin)
n3=nbandmat(n3type,n3orb,n3spin)
n4=nbandmat(n4type,n4orb,n4spin)

!n1=n1spin
!n2=n2spin
!n3=n3spin
!n4=n4spin

!n1=nbandmat(n1orb,n1spin)
!n2=nbandmat(n2orb,n2spin)
!n3=nbandmat(n3orb,n3spin)
!n4=nbandmat(n4orb,n4spin)


i=jmat(n1,n2,n3,n4,p1,p2,p3)


!if(p1==1 .and. p2==1 .and. p3==1)then


kpx=p4vec(n1,n2,n3,p1,p2,p3,1)
Kpy=p4vec(n1,n2,n3,p1,p2,p3,2)

CALL Unitary(kpx,kpy,Umatp4)


Do n1p=1,ntype
Do n2p=1,ntype
Do n3p=1,ntype
Do n4p=1,ntype


!----------------------------------------------------------------------
! Local interaction:
!----------------------------------------------------------------------

n1prime=nbandmat(n1p,n1orb,n1spin)
n2prime=nbandmat(n2p,n2orb,n2spin)
n3prime=nbandmat(n3p,n3orb,n3spin)
n4prime=nbandmat(n4p,n4orb,n4spin)

iprime=jmat(n1prime,n2prime,n3prime,n4prime,p1,p2,p3)



YvertexABinitial(iprime)=((kdelta(n1p,n2p))*(kdelta(n3p,n4p))*(kdelta(n1p,n3p)))*&
            ((u)*(kdelta(n1orb,n2orb)*kdelta(n3orb,n4orb)*kdelta(n1orb,n3orb))&
              *(((kdelta(n1spin,n3spin))*(kdelta(n2spin,n4spin)))-((kdelta(n1spin,n4spin))*(kdelta(n2spin,n3spin))))&
            +(u)*(kdelta(n1spin,n2spin)*kdelta(n3spin,n4spin)*kdelta(n1spin,n3spin))&                          !(uprime-Jcoup)
              *(((kdelta(n1orb,n3orb))*(kdelta(n2orb,n4orb)))-((kdelta(n1orb,n4orb))*(kdelta(n2orb,n3orb))))&
            +(u)*(1.0d+00-kdelta(n2orb,n1orb))*(1.0d+00-kdelta(n1spin,n2spin))&
                     *((kdelta(n1spin,n3spin))*(kdelta(n2spin,n4spin))*(kdelta(n1orb,n3orb))*(kdelta(n2orb,n4orb))&
                      -(kdelta(n2spin,n3spin))*(kdelta(n1spin,n4spin))*(kdelta(n2orb,n3orb))*(kdelta(n1orb,n4orb))))


!YvertexABinitial(iprime)=((kdelta(n1p,n2p))*(kdelta(n3p,n4p))*(kdelta(n1p,n3p)))*&
!          (u)*(((kdelta(n1spin,n3spin))*(kdelta(n2spin,n4spin)))-((kdelta(n1spin,n4spin))*(kdelta(n2spin,n3spin))))

!----------------------------------------------------------------------

Yverdummyu=cmplx(YvertexABinitial(iprime),0.0d+00)

Yverdummy=Yverdummyu


Yvertexcomplex(i)=Yvertexcomplex(i)+Yverdummy*(CONJG(Umattype(n1p,n1type,p1)))*(CONJG(Umattype(n2p,n2type,p2)))&
                 *Umattype(n3p,n3type,p3)*Umatp4(n4p,n4type)


end do !n4p
end do !n3p
end do !n2p
end do !n1p

!----------------------------------------------------------------------
!endif !(p1=p2=p3=1)
!----------------------------------------------------------------------

end do !t4
end do !t3
end do !t2
end do !t1

end do !s4
end do !s3
end do !s2
end do !s1

end do !o4
end do !o3
end do !o2
end do !o1

end do !p3
end do !p2
end do !p1
!------------------------------------------------------------------
!------------------------------------------------------------------
!----------------------------------------------------------------------
! Solving the flow equations:
!----------------------------------------------------------------------
MF=10
ML=1
MU=1
!NEQ=Ntot 
!----------------------------------------------------------------------
!Yvertexre(1:NEQ)=real(Yvertexcomplex(1:NEQ))

NEQ=((nband/2)**4)*(npatch**3)

Do n1=5,8
Do n2=5,8
Do n3=5,8
Do n4=5,8

Do p1=1,npatch
Do p2=1,npatch
Do p3=1,npatch

i=imat(n1,n2,n3,n4,p1,p2,p3)
j=jmat(n1,n2,n3,n4,p1,p2,p3)


Yvertexre(i)=real(Yvertexcomplex(j))


end do
end do
end do

end do
end do
end do
end do

!----------------------------------------------------------------------
! Dimension Declaration of WORK:
!----------------------------------------------------------------------
If(MF==10)then
LRW=20+16*NEQ
LIW=20
elseif(MF==21 .or. MF==22)then
LRW=22+9*NEQ+NEQ**2 
LIW=20+NEQ
elseif(MF==24 .or. MF==25)then
LRW=22+10*NEQ+(2*ML+MU)*NEQ
LIW=20+NEQ
endif

allocate(RWORK(LRW),IWORK(LIW))
!------------------------------------------------------------------
! Initialization of the solver:
!------------------------------------------------------------------
Flag=1
!------------------------------------------------------------------
IOPT=1            
RWORK(5:10)=0.0d+00 
IWORK(5:10)=0
IWORK(6)=20000 
resol=1.20d+00!1.05d+00!1.05d+00!
!------------------------------------------------------------------
maxvertex=Yvertexre(imat(5,6,5,6,1,1,1))
print*,'maxvertex_outside=',maxvertex
!------------------------------------------------------------------
print*,'final xout=',10.00d+00*t/((resol)**(700.0d+00))
!------------------------------------------------------------------
!icall=0

Do i=1,700

irout=i

if(abs(maxvertex)<14.0d+00*t)then
!------------------------------------------------------------------
RTOL=t*10.0d+00**(-3.0d+00)   !converged:10**(-5)
ATOL=t*10.0d+00**(-3.0d+00)
ISTATE=1
x=10.00d+00*t/((resol)**(i-1))
xout=10.00d+00*t/((resol)**(i))
!------------------------------------------------------------------
icall=0
!------------------------------------------------------------------

CALL Fflow (NEQ, x , Yvertexre, YDvertexre)

!Yvertexre(1:NEQ)=YDvertexre(1:NEQ)

CALL DLSODE (Fflow, NEQ, Yvertexre, x, xout, ITOL, RTOL, ATOL, Flag,ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF)


print*,"Fcur",i,ISTATE,10.00d+00*t/((resol)**(i-1)),xout/t,maxvertex

!write(3,*)i,10.00d+00*t/((resol)**(i-1)),xout/t,maxvertex,p3max

print*,'max-vertex belongs to:',imax,n1mat(imax),n2mat(imax),n3mat(imax),n4mat(imax)&
      ,p1mat(imax),p2mat(imax),p3mat(imax)
!------------------------------------------------------------------
endif
!------------------------------------------------------------------
end do !i
!------------------------------------------------------------------
deallocate(RWORK,IWORK)
!------------------------------------------------------------------
!------------------------------------------------------------------
!if(xout<0)then

Do p1=1,npatch
Do p2=1,npatch
Do p3=1,npatch

Do n1=5,8!1,nband
Do n2=5,8!1,nband
Do n3=5,8!1,nband
Do n4=5,8!1,nband


i=imat(n1,n2,n3,n4,p1,p2,p3)

k=p4int(n1,n2,n3,p1,p2,p3)

write(1,*)Yvertexre(i),k,nfold(n1,n2,n3,p1,p2,p3)

end do !n4
end do !n3
end do !n2
end do !n1

end do !p3
end do !p2
end do !p1

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!end if ! sup
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
end do ! ieF: Fermi-Energy
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!end if ! sup
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
call cpu_time(stop_time)
print*,"total time",stop_time - start_time
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
end program HFRG
!------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------
!******* Flow equations***********************
!------------------------------------------------------------------------------------------------------------
subroutine Fflow(NEQ, x, Yvertexre, YDvertexre)

use commons

implicit none

common eF,maxvertex,imat,start_time,irout,p4int,thetap,KFSpmat,imax!,Yvertex,imax

!-------------------------------------------------------------------------------------------------------
! Parameter declaration:
!-------------------------------------------------------------------------------------------------------
INTEGER::k,i,j1,j2,NEQ,n1,n2,n3,n4,p1,p2,p3,p4,kp,l1,l2,kp2,nrootskp2,j3,j4,nch,j,irot,irout&
         ,nroots1,nroots2,nroots3,nrootsk,nrootskp,sint,icommon,nchcommon,imin,imax,kr,kv,kstart,kfinal&
         ,p2max,p3max,n2max,n4max,iasym,iasym2,kint,kk
DOUBLE PRECISION::eF,y,scalecommon,groupvel,Matsumpart1,Matsumpart2,step,step1,step2&
                  ,angle1,angle2,angle3,KFS1,KFS2,KFS3,KFSk,KFSkp,stepsize&
                  ,anglek,anglekp,energy,Matsum,anglekp2,KFSkp2,radkp,radkp2&
                  ,anglekpred,kpx,kpy,maxvertex,kp2x,kp2y,anglekp2red&
                  ,x,minvertex,r3,dis
DOUBLE PRECISION,DIMENSION(1:2)::sgn
!-------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------
! External functions:
!-------------------------------------------------------------------------------------------------------
DOUBLE PRECISION::dispersion,fermi
INTEGER,DIMENSION(1:Ntot):: n1mat,n2mat,n3mat,n4mat,p1mat,p2mat,p3mat
INTEGER,DIMENSION(1:Ntot):: nrep
!INTEger,DIMENSION(1:nband,1:nband,1:nband,1:nband,1:npatch,1:npatch,1:npatch)::imat
INTEGER,DIMENSION(5:8,5:8,5:8,5:8,1:npatch,1:npatch,1:npatch)::imat
INTEger,DIMENSION(1:nband,1:nband,1:nband,1:npatch,1:npatch,1:npatch)::p4int
!-------------------------------------------------------------------------------------------------------
DOUBLE PRECISION::Yang,YDang
DOUBLE PRECISION,DIMENSION(1:NEQ)::Taupp,Tauphc,Tauphd
!-------------------------------------------------------------------------------------------------------
DOUBLE PRECISION,DIMENSION(1:2,1:nband,1:npatch)::KFSmat,gvelFSmat
INTEGER,DIMENSION(1:2,1:nband,1:npatch)::nrootsmat
DOUBLE PRECISION::Matsumpp,Matsumph

!DOUBLE PRECISION, DIMENSION(1:nband,1:nband,1:nband,1:nband,1:npatch,1:npatch,1:npatch)::Matsumpp1,Matsumpp2&
!                                                                                        ,Matsumph1,Matsumph2

DOUBLE PRECISION, DIMENSION(4:5,4:5,4:5,4:5,1:npatch,1:npatch,1:npatch)::Matsumpp1,Matsumpp2&
                                                                                        ,Matsumph1,Matsumph2

INTEGER,DIMENSION(1:2,1:npatch,1:npatch,1:npatch)::kpmat

DOUBLE PRECISION,DIMENSION(1:nband,1:npatch)::thetap
DOUBLE PRECISION,DIMENSION(1:nband,1:npatch)::KFSpmat
INTEGER,DIMENSION(1:nband,1:npatch)::nrootspmat


DOUBLE PRECISION,DIMENSION(1:2,1:nband,1:npatch)::thetakmat

INTEGER,DIMENSION(1:2,1:2,1:npatch,1:npatch,1:npatch)::kpFSmat

INTEGER::BZid,n1a,n2a,l1a,l2a,n1b,n2b,l1b,l2b,n1i,n2i,l1i,l2i

INTEGER,DIMENSION(1:nband)::ns

INTEGER,DIMENSION(1:nband,1:nband,1:nband,1:nband,1:npatch,1:npatch,1:npatch)::iredmat
DOUBLE PRECISION,DIMENSION(1:2*Ntotal)::Yvertexcommon,Yvertex,YDvertex
INTEGER,DIMENSION(1:Ntot)::iredmattoimat

DOUBLE PRECISION,DIMENSION(1:NEQ)::Yvertexre,YDvertexre

DOUBLE PRECISION::KFSthetamax
INTEGER,DIMENSION(1:2,1:nband,1:npatch)::kintmat
INTEGER,DIMENSION(1:nband)::nref,nmin,nmax
INTEGER::nref1,nref2,n1m,n1p,n2m,n2p,l1m,l1p,l2m,l2p&
         ,n1r,n2r,l1r,l2r,n3r,nrootstotal

REAL:: start_time, stop_time
!------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------
scalecommon=x
icall=icall+1
stepsize=10.0d+00**(-6.0d+00)

!------------------------------------------------------------------
kpmat(1:2,1:npatch,1:npatch,1:npatch)=0
kpFSmat(1:2,1:2,1:npatch,1:npatch,1:npatch)=0

Matsumpp=0.0d+00
Matsumph=0.0d+00

Matsumpp1(4:5,4:5,4:5,4:5,1:npatch,1:npatch,1:npatch)=0.0d+00!(1:nband,1:nband,1:nband,1:nband,1:npatch,1:npatch,1:npatch)=0.0d+00
Matsumpp2(4:5,4:5,4:5,4:5,1:npatch,1:npatch,1:npatch)=0.0d+00!(1:nband,1:nband,1:nband,1:nband,1:npatch,1:npatch,1:npatch)=0.0d+00

Matsumph1(4:5,4:5,4:5,4:5,1:npatch,1:npatch,1:npatch)=0.0d+00!(1:nband,1:nband,1:nband,1:nband,1:npatch,1:npatch,1:npatch)=0.0d+00
Matsumph2(4:5,4:5,4:5,4:5,1:npatch,1:npatch,1:npatch)=0.0d+00!(1:nband,1:nband,1:nband,1:nband,1:npatch,1:npatch,1:npatch)=0.0d+00


gvelFSmat(1:2,1:nband,1:npatch)=0.0d+00
!KFSmat(1:2,1:npatch)=0.0d+00

nrootsmat(1:2,1:nband,1:npatch)=0

nrep(1:NEQ)=0

sgn(1)=-1.0d+00
sgn(2)=1.0d+00


if(Idbandsym==1)then

nref(1:(nband/2))=(nband/2)
nref(((nband/2)+1):nband)=(nband/2)+1

nmin((nband/2))=1
nmax((nband/2))=(nband/2)
nmin((nband/2)+1)=(nband/2)+1
nmax((nband/2)+1)=nband
nref1=5!(nband/2)
nref2=5!(nband/2)+1

elseif(Idbandsym==0)then

nref1=1
nref2=nband

endif

YDvertexre(1:NEQ)=0.0d+00
!------------------------------------------------------------------
!------------------------------------------------------------------
nrootstotal=0
!------------------------------------------------------------------
Do k=1,npatch

Do sint=1,2

energy=eF+(sgn(sint)*x)

Do l1=1,nband

anglek=thetap(l1,k)
 
!------------------------------------------------------------------

CALL roots_dispersion(l1,anglek,energy,nrootsk,KFSk)

if(nrootsk>0 )then

!------------------------------------------------------------------ 

kpx=KFSk*cos(anglek)
kpy=KFSk*sin(anglek)

CALL Hexagonal(kpx,kpy,BZid)

if(BZid==0)then

CALL Folding(kpx,kpy,kpx,kpy)

KFSk=sqrt((kpx**2)+(kpy**2))

anglek=atan2(kpy,kpx)
anglek=(-SIGN(1.0d+00,anglek)+1.0d+00)*(pi)+anglek

endif

!------------------------------------------------------------------ !add-end

if(KFSk .LE. 0.0001d+00)then
nrootsk=0
endif

endif
!------------------------------------------------------------------
!------------------------------------------------------------------
if(nrootsk>0)then

CALL patch_lab_Ge(patchflag,anglek,KFSk,kint,anglekpred)

thetakmat(sint,l1,kint)=anglek
KFSmat(sint,l1,kint)=KFSk
nrootsmat(sint,l1,kint)=nrootsk
gvelFSmat(sint,l1,kint)=abs(groupvel(1,l1,KFSK,anglek))

!------------------------------------------------------------------
!------------------------------------------------------------------

kintmat(sint,l1,k)=kint

endif
!------------------------------------------------------------------

nrootstotal=nrootstotal+nrootsk

end do !l1

end do !sint
end do !k

print*,'nrootstotal=',nrootstotal
!------------------------------------------------------------------
if(nrootstotal>0)then
!------------------------------------------------------------------

Do sint=1,2

energy=eF+(sgn(sint)*x)

Do l1=nref1,nref2!1,nband

Do k=1,npatch

anglek=thetakmat(sint,l1,k)

if (nrootsmat(sint,l1,k)>0)then

Do n1=nref1,nref2!1,nband

Do p1=1,npatch

angle1=thetap(n1,p1)
KFS1=KFSpmat(n1,p1)

!p3=p1

Do n2=nref1,nref2!1,nband

Do p2=1,npatch

angle2=thetap(n2,p2)
KFS2=KFSpmat(n2,p2)

Do l2=nref1,nref2!1,nband

!------------------------------------------------------------------ 
!------------------------------------------------------------------ 

kpx=KFSpmat(n1,p1)*cos(thetap(n1,p1))+KFSpmat(n2,p2)*cos(thetap(n2,p2))-KFSmat(sint,l1,k)*cos(thetakmat(sint,l1,k))
kpy=KFSpmat(n1,p1)*sin(thetap(n1,p1))+KFSpmat(n2,p2)*sin(thetap(n2,p2))-KFSmat(sint,l1,k)*sin(thetakmat(sint,l1,k))

anglekp=atan2(kpy,kpx)
anglekp=(-SIGN(1.0d+00,anglekp)+1.0d+00)*(pi)+anglekp

radkp=sqrt(((kpx)**2)+((kpy)**2))

!------------------------------------------------------------------
! Regulator:
!------------------------------------------------------------------
if(abs(dispersion(l2,radkp,anglekp,eF))>scalecommon+stepsize)then
step=1.0d+00
elseif(abs(dispersion(l2,radkp,anglekp,eF))<scalecommon-stepsize)then
step=0.0d+00
else
step=0.5d+00
endif
!------------------------------------------------------------------ 
!------------------------------------------------------------------
!------------------------------------------------------------------
! 1.pp channel
!------------------------------------------------------------------

!if(abs(dispersion(l2,radkp,anglekp,eF)+energy-eF)< (10.0d+00**(-12.0d+00)))then

!Matsumpp=(beta/4.0d+00)*(((tanh(beta*(energy-eF)/2.0d+00))**2)-1.0d+00)


!else


!Matsumpp=((fermi(energy-eF)-fermi(-dispersion(l2,radkp,anglekp,eF)))&
!                                      /(energy-eF+dispersion(l2,radkp,anglekp,eF)))

!endif


!------------------------------------------------------------------
!new:

if(energy-eF < 0.0d+00 .and. dispersion(l2,radkp,anglekp,eF)>0.0d+00)then
Matsumpp=0.0d+00
elseif(energy-eF .GE. 0.0d+00 .and. dispersion(l2,radkp,anglekp,eF) .LE. 0.0d+00)then
Matsumpp=0.0d+00
elseif(energy-eF < 0.0d+00)then
Matsumpp=1.0d+00/(dispersion(l1,KFSmat(sint,l1,k),thetakmat(sint,l1,k),eF)+dispersion(l2,radkp,anglekp,eF))
elseif(energy-eF .GE. 0.0d+00)then
Matsumpp=-1.0d+00/(dispersion(l1,KFSmat(sint,l1,k),thetakmat(sint,l1,k),eF)+dispersion(l2,radkp,anglekp,eF))
endif


if(abs(dispersion(l1,KFSmat(sint,l1,k),thetakmat(sint,l1,k),eF)+dispersion(l2,radkp,anglekp,eF))<10.0d+00**(-14.0d+00))then
!print*,'careful! devided by zero!'
Matsumpp=(beta/4.0d+00)*(((tanh(beta*(energy-eF)/2.0d+00))**2)-1.0d+00)
endif
!------------------------------------------------------------------

Matsumpp=step*Matsumpp*(((KFSmat(sint,l1,k))/(gvelFSmat(sint,l1,k))))

!------------------------------------------------------------------

if(sint==1)then
Matsumpp1(n1,n2,l1,l2,p1,p2,k)=Matsumpp
elseif(sint==2)then
Matsumpp2(n1,n2,l1,l2,p1,p2,k)=Matsumpp
endif

!------------------------------------------------------------------
!------------------------------------------------------------------
! 2.ph channel
!------------------------------------------------------------------

kpx=KFSpmat(n1,p1)*cos(thetap(n1,p1))-KFSpmat(n2,p2)*cos(thetap(n2,p2))+KFSmat(sint,l1,k)*cos(thetakmat(sint,l1,k))
kpy=KFSpmat(n1,p1)*sin(thetap(n1,p1))-KFSpmat(n2,p2)*sin(thetap(n2,p2))+KFSmat(sint,l1,k)*sin(thetakmat(sint,l1,k))

anglekp=atan2(kpy,kpx)
anglekp=(-SIGN(1.0d+00,anglekp)+1.0d+00)*(pi)+anglekp

radkp=sqrt(((kpx)**2)+((kpy)**2))

!------------------------------------------------------------------
! Regulator:
!------------------------------------------------------------------
if(abs(dispersion(l2,radkp,anglekp,eF))>scalecommon+stepsize)then
step=1.0d+00
elseif(abs(dispersion(l2,radkp,anglekp,eF))<scalecommon-stepsize)then
step=0.0d+00
else
step=0.5d+00
endif
!------------------------------------------------------------------ 
!------------------------------------------------------------------

!if(abs(dispersion(l2,radkp,anglekp,eF)-energy+eF)< (10.0d+00**(-12.0d+00)))then

!Matsumph=(beta/4.0d+00)*(((tanh(beta*(energy-eF)/2.0d+00))**2)-1.0d+00)

!else

!Matsumph=((fermi(energy-eF))-(fermi(dispersion(l2,radkp,anglekp,eF))))/(energy-eF-dispersion(l2,radkp,anglekp,eF))

!end if

!------------------------------------------------------------------
!new:

if(energy-eF < 0.0d+00 .and. dispersion(l2,radkp,anglekp,eF)<0.0d+00)then
Matsumph=0.0d+00
elseif(energy-eF .GE. 0.0d+00 .and. dispersion(l2,radkp,anglekp,eF) .GE. 0.0d+00)then
Matsumph=0.0d+00
elseif(energy-eF < 0.0d+00)then
Matsumph=1.0d+00/(dispersion(l1,KFSmat(sint,l1,k),thetakmat(sint,l1,k),eF)-dispersion(l2,radkp,anglekp,eF))
elseif(energy-eF .GE. 0.0d+00)then
Matsumph=-1.0d+00/(dispersion(l1,KFSmat(sint,l1,k),thetakmat(sint,l1,k),eF)-dispersion(l2,radkp,anglekp,eF))
endif

if(abs(dispersion(l1,KFSmat(sint,l1,k),thetakmat(sint,l1,k),eF)-dispersion(l2,radkp,anglekp,eF))<10.0d+00**(-14.0d+00))then

!print*,'careful! devided by zero!'

Matsumph=(beta/4.0d+00)*(((tanh(beta*(energy-eF)/2.0d+00))**2)-1.0d+00)

!print*,n1,n2,l1,l2,p1,p2,k,Matsumph,((tanh(beta*(energy-eF)/2.0d+00))**2),(beta/4.0d+00),energy-eF

endif

!------------------------------------------------------------------

Matsumph=step*Matsumph*((KFSmat(sint,l1,k))/(gvelFSmat(sint,l1,k)))

!------------------------------------------------------------------

if(sint==1)then
Matsumph1(n1,n2,l1,l2,p1,p2,k)=Matsumph
elseif(sint==2)then
Matsumph2(n1,n2,l1,l2,p1,p2,k)=Matsumph
endif
!------------------------------------------------------------------
!------------------------------------------------------------------ 
!------------------------------------------------------------------ 
!------------------------------------------------------------------
!------------------------------------------------------------------ 
end do !l2
!------------------------------------------------------------------
!------------------------------------------------------------------
end do !p2
end do !n2
end do !p1
end do !n1
!------------------------------------------------------------------
end if !k on FS
end do !k
end do !l1

end do !sint
!------------------------------------------------------------------
!------------------------------------------------------------------
!------------------------------------------------------------------

Do n1=5,8!1,nband

!------------------------------------------------------------------
! Employing anti-sym prop of vertex functions: 
!------------------------------------------------------------------
if(Idantisym==0)then
n2max=nband
elseif(Idantisym==1)then
n2max=n1
endif
!------------------------------------------------------------------

Do n2=5,n2max!1,n2max

!------------------------------------------------------------------
!------------------------------------------------------------------
Do n4=5,8!1,nband!n4max

Do n3=5,8!1,nband

Do p1=1,npatch

Do p2=1,npatch

!------------------------------------------------------------------
if(idFSrotsym==1 .and. p1>nprot .and. p2>nprot)then
p3max=nprot
else
p3max=npatch
endif
!------------------------------------------------------------------

Do p3=1,p3max

i=imat(n1,n2,n3,n4,p1,p2,p3)

if(nrep(i)==0)then


p4=p4int(n1,n2,n3,p1,p2,p3)

if(p4==0)then

YDvertexre(i)=0.0d+00
nrep(i)=1

elseif(p4>0)then

!------------------------------------------------------------------
!------------------------------------------------------------------
!------------------------------------------------------------------
!------------------------------------------------------------------
if(i==1)then
maxvertex=Yvertexre(i)
imax=1
elseif(abs(Yvertexre(i))>abs(maxvertex))then
maxvertex=Yvertexre(i)
imax=i
endif

if(i==1)then
minvertex=Yvertexre(i)
imin=1
elseif(abs(Yvertexre(i))<abs(minvertex))then
minvertex=Yvertexre(i)
imin=i
endif

!------------------------------------------------------------------
YDvertexre(i)=0.0d+00
!------------------------------------------------------------------
nchcommon=nch
!------------------------------------------------------------------
!------------------------------------------------------------------
! Right-hand side of the flow eqs:
!------------------------------------------------------------------
Taupp(i)=0.0d+00
Tauphc(i)=0.0d+00
Tauphd(i)=0.0d+00
!------------------------------------------------------------------
!------------------------------------------------------------------
Do sint=1,2

Do k=1,npatch

!------------------------------------------------------------------
!------------------------------------------------------------------
! Evaluation of the vertex:
!------------------------------------------------------------------
Do l1=5,8!1,nband

if (nrootsmat(sint,l1,k)>0)then

Do l2=5,8!1,nband 
!------------------------------------------------------------------
if(Idbandsym==1)then
n1r=nref(n1)
n2r=nref(n2)
l1r=nref(l1)
l2r=nref(l2)
n3r=nref(n3)
elseif(Idbandsym==0)then
n1r=n1
n2r=n2
n3r=n3
l1r=l1
l2r=l2
endif
!------------------------------------------------------------------
! (p-p channel)
!------------------------------------------------------------------
!------------------------------------------------------------------

p4=p4int(n1,n2,n3,p1,p2,p3)

j1=imat(n1,n2,l1,l2,p1,p2,k)  !term1
j2=imat(n3,n4,l1,l2,p3,p4,k)

if(sint==1)then

Taupp(i)=Taupp(i)-(1.0d+00/((1.0d+00)*(npatch)))*Yvertexre(j1)*Yvertexre(j2)&
                 *Matsumpp1(n1r,n2r,l1r,l2r,p1,p2,k)

elseif(sint==2)then

Taupp(i)=Taupp(i)-(1.0d+00/((1.0d+00)*(npatch)))*Yvertexre(j1)*Yvertexre(j2)&
                 *Matsumpp2(n1r,n2r,l1r,l2r,p1,p2,k)

endif
!------------------------------------------------------------------
!end if !kp on FS
!------------------------------------------------------------------
! (p-h-c channel)
!------------------------------------------------------------------

p4=p4int(n1,n2,n3,p1,p2,p3)

j1=imat(n3,l1,n1,l2,p3,k,p1)  !term4
j2=imat(n2,l1,n4,l2,p2,k,p4)

j3=imat(n1,l1,n3,l2,p1,k,p3)  !term2
j4=imat(n4,l1,n2,l2,p4,k,p2)


if(sint==1)then

Tauphc(i)=Tauphc(i)-1.0d+00*(1.0d+00/((1.0d+00)*(npatch)))&
            *(((Yvertexre(j1)*Yvertexre(j2))*(Matsumph1(n3r,n1r,l1r,l2r,p3,p1,k)))&
             +((Yvertexre(j3)*Yvertexre(j4))*(Matsumph1(n1r,n3r,l1r,l2r,p1,p3,k)))) 

elseif(sint==2)then 

Tauphc(i)=Tauphc(i)-1.0d+00*(1.0d+00/((1.0d+00)*(npatch)))&
            *(((Yvertexre(j1)*Yvertexre(j2))*(Matsumph2(n3r,n1r,l1r,l2r,p3,p1,k)))&
             +((Yvertexre(j3)*Yvertexre(j4))*(Matsumph2(n1r,n3r,l1r,l2r,p1,p3,k)))) 

endif
!------------------------------------------------------------------
!------------------------------------------------------------------
!(p-h-d channel)
!------------------------------------------------------------------

p4=p4int(n1,n2,n3,p1,p2,p3)           

j1=imat(n4,l1,n1,l2,p4,k,p1)  !term3
j2=imat(n2,l1,n3,l2,p2,k,p3)

j3=imat(n1,l1,n4,l2,p1,k,p4)  !term5
j4=imat(n3,l1,n2,l2,p3,k,p2)


if(sint==1)then

Tauphd(i)=Tauphd(i)+1.0d+00*(1.0d+00/((1.0d+00)*(npatch)))&
            *(((Yvertexre(j1)*Yvertexre(j2))*(Matsumph1(n2r,n3r,l1r,l2r,p2,p3,k)))&
             +((Yvertexre(j3)*Yvertexre(j4))*(Matsumph1(n3r,n2r,l1r,l2r,p3,p2,k))))

elseif(sint==2)then

Tauphd(i)=Tauphd(i)+1.0d+00*(1.0d+00/((1.0d+00)*(npatch)))&
            *(((Yvertexre(j1)*Yvertexre(j2))*(Matsumph2(n2r,n3r,l1r,l2r,p2,p3,k)))&
             +((Yvertexre(j3)*Yvertexre(j4))*(Matsumph2(n3r,n2r,l1r,l2r,p3,p2,k))))

endif
!------------------------------------------------------------------
end do !l2

endif ! k on the const-energy surface

end do !l1
!------------------------------------------------------------------
!------------------------------------------------------------------

end do !k

end do !sint
!----------------------------------------------------------------------


YDvertexre(i)=(Taupp(i)+Tauphc(i)+Tauphd(i))/(2.0d+00*pi)

nrep(i)=1

!------------------------------------------------------------------
! vertex recostruction employing anti-symmetric properties:
!------------------------------------------------------------------
if(Idantisym==1)then

iasym=imat(n2,n1,n3,n4,p2,p1,p3)

if(nrep(iasym)==0 .and. n1>n2)then

YDvertexre(iasym)=-YDvertexre(i)                   ! real part reconstruction 
!YDvertex(iasym+(NEQ/2))=-YDvertex(i+(NEQ/2))   ! imaginary part reconstruction 

nrep(iasym)=1

endif

endif

!----------------------------------------------------------------------
!----------------------------------------------------------------------
! vertex recostruction due to the rot-sym of Fermi surface!
!----------------------------------------------------------------------

Do j1=1,nprotpatch-1

if(p1+j1*nprot<(npatch+1) .and. p2+j1*nprot<(npatch+1) .and. p3+j1*nprot<(npatch+1) .and. idFSrotsym==1)then

irot=imat(n1,n2,n3,n4,p1+j1*nprot,p2+j1*nprot,p3+j1*nprot)


if(nrep(irot)==0)then

YDvertexre(irot)=YDvertexre(i)                   ! real part reconstruction
!YDvertex(irot+(NEQ/2))=YDvertex(i+(NEQ/2))   ! imaginary part reconstruction

nrep(irot)=1
endif

if(Idantisym==1)then

iasym=imat(n2,n1,n3,n4,p2+j1*nprot,p1+j1*nprot,p3+j1*nprot)


if(nrep(iasym)==0)then

YDvertexre(iasym)=-YDvertexre(i)                   ! real part reconstruction
!YDvertex(iasym+(NEQ/2))=-YDvertex(i+(NEQ/2))   ! imaginary part reconstruction
nrep(iasym)=1

endif

endif


end if ! p1,p2,p3 inside BZ!

end do !j1
!----------------------------------------------------------------------
!----------------------------------------------------------------------

end if !k4 very small

endif !if (nrep(i)==0)


if(n4>n3 .and. nrep(i)==0)then
print*,'not covered!!!!!!'
endif

end do! p3

!----------------------------------------------------------------------
!----------------------------------------------------------------------

end do! p2
end do! p1

end do! n4
end do! n3
!----------------------------------------------------------------------
!----------------------------------------------------------------------
end do! n2
end do! n1
!----------------------------------------------------------------------
end if !nrootstotal>0
!----------------------------------------------------------------------
!----------------------------------------------------------------------
print*,'icall=',icall,'scale/t=',scalecommon/t,'maxV=',maxvertex,'minV=',minvertex
print*,'max-vertex belongs to:',imax
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
end subroutine Fflow
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
double precision function fermi(w)

use commons

implicit none

common eF

double precision,intent(in)::w
DOUBLE PRECISION::eF,scalecommon

fermi=(1.0d+00-tanh(beta*(w)/2.0d+00))/2.0d+00

end function fermi
!----------------------------------------------------------------------
!----------------------------------------------------------------------
DOUBLE PRECISION function dispersion(blab,mom,theta,energy)

use commons

implicit none

common eF

DOUBLE PRECISION,INTENT(IN)::mom,theta,energy
INTEGER, INTENT(IN)::blab
DOUBLE PRECISION::eF,scalecommon,kx,ky
!----------------------------------------------------------------------
! Spin-orbit coupling:
!----------------------------------------------------------------------
DOUBLE PRECISION, DIMENSION(1:nband)::Emat
COMPLEX*16::VSO12,VSO21
DOUBLE PRECISION::VSOx,VSOy,VSOz,ESOp,ESOn
!----------------------------------------------------------------------
!------------------------------------------------------------------
! parameters for the Eispack library:
!------------------------------------------------------------------
INTEGER*4::Mdim,EVerr
REAL*8,ALLOCATABLE,DIMENSION(:,:)::Mmatre,Mmatim,MEVecre,MEVecim
REAL*8,ALLOCATABLE,DIMENSION(:)::MEValre,MEValim
LOGICAL::EVflag
!------------------------------------------------------------------
INTEGER::n1,n2
DOUBLE PRECISION::Bfun,Dfun,Efun
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! for Hexagonal lattice:
!----------------------------------------------------------------------

kx=mom*cos(theta)
ky=mom*sin(theta)

dispersion=0.0d+00

if(blab .LE. (nband/2.0d+00))then        

dispersion=1.0*t*sqrt(3.0d+00+2.0d+00*(cos(Lconst*kx*(sqrt(3.0d+00))))&
                +4.0d+00*(Cos(Lconst*kx*(sqrt(3.0d+00))/2.0d+00))*(cos(Lconst*1.5d+00*ky)))-energy

else 

dispersion=-1.0d+00*t*sqrt(3.0d+00+2.0d+00*(cos(Lconst*kx*(sqrt(3.0d+00))))&
                +4.0d+00*(Cos(Lconst*kx*(sqrt(3.0d+00))/2.0d+00))*(cos(Lconst*1.5d+00*ky)))-energy


endif
!----------------------------------------------------------------------
!----------------------------------------------------------------------
end function dispersion
!----------------------------------------------------------------------
!----------------------------------------------------------------------
DOUBLE PRECISION function groupvel(Id,blab,mom,theta)

use commons

implicit none

common eF

DOUBLE PRECISION,INTENT(IN)::mom,theta
INTEGER,INTENT(IN)::id,blab
DOUBLE PRECISION::eF,scalecommon,grouprad,groupangle,kx,ky &
                 ,groupradSO,groupangleSO&
                 ,deltakx,deltaky,groupvelkx,groupvelky,dispersion&
                 ,anglekxp,anglekxn,anglekyp,anglekyn&
                 ,kradkxp,kradkxn,kradkyp,kradkyn&
                 ,ek0,ekxp,ekxn,ekyp,ekyn,angle,dis,r3
INTEGER::sgn

!----------------------------------------------------------------------
if(blab .LE. nband/2.0d+00)then
sgn=1
else
sgn=-1
endif
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!--------------------------------------
! Analytical:
!--------------------------------------
kx=mom*(cos(theta))
ky=mom*(sin(theta))


!CALL Folding(kx,ky,kx,ky)

!angle=atan2(ky,kx)
!angle=(-SIGN(1.0d+00,angle)+1.0d+00)*(pi)+angle


!groupvelkx=(-sgn)*(t/2.0d+00)*(1.0d+00/(sqrt(3.0d+00+(2.0d+00*(cos(Lconst*(sqrt(3.0d+00))*kx)))&
!                  +(4.0d+00*(cos(((sqrt(3.0d+00))/2.0d+00)*Lconst*kx))*(cos(1.5d+00*Lconst*ky))))))&
!           *((sqrt(3.0d+00))*Lconst*2.0d+00)*((sin((sqrt(3.0d+00))*Lconst*kx))&
!           +((sin((sqrt(3.0d+00)/2.0d+00)*Lconst*kx))*(cos(1.5d+00*Lconst*ky))))

!groupvelky=(-sgn)*(t/2.0d+00)*(1.0d+00/(sqrt(3.0d+00+(2.0d+00*(cos(Lconst*(sqrt(3.0d+00))*kx)))&
!                  +(4.0d+00*(cos(((sqrt(3.0d+00))/2.0d+00)*Lconst*kx))*(cos(1.5d+00*Lconst*ky))))))&
!           *(6.0d+00*Lconst)*(cos(((sqrt(3.0d+00))/2.0d+00)*Lconst*kx))*(sin(1.5d+00*Lconst*ky))


groupvelkx=(-sgn)*(t/2.0d+00)*(1.0d+00/(t*sqrt(3.0d+00+2.0d+00*(cos(Lconst*kx*(sqrt(3.0d+00))))&
                +4.0d+00*(cos(Lconst*kx*(sqrt(3.0d+00))/2.0d+00))*(cos(Lconst*1.5d+00*ky)))))&
           *((2.0d+00*(sqrt(3.0d+00))*Lconst*(-sin(kx*Lconst*(sqrt(3.0d+00)))))&
            +(2.0d+00*(sqrt(3.0d+00))*Lconst*(-sin(kx*Lconst*(sqrt(3.0d+00))/2.0d+00))*(cos(1.5d+00*ky*Lconst))))


groupvelky=(-sgn)*(t/2.0d+00)*(1.0d+00/(t*sqrt(3.0d+00+2.0d+00*(cos(Lconst*kx*(sqrt(3.0d+00))))&
                +4.0d+00*(cos(Lconst*kx*(sqrt(3.0d+00))/2.0d+00))*(cos(Lconst*1.5d+00*ky)))))&
           *(6.0d+00*Lconst*(cos(kx*Lconst*(sqrt(3.0d+00))/2.0d+00))*(-sin(Lconst*1.5d+00*ky)))

!--------------------------------------

!grouprad=(groupvelkx)*(cos(theta))+(groupvelky)*(sin(theta))

!--------------------------------------

grouprad=(-sgn)*(t)*(Lconst)*(1.0d+00/sqrt(3.0d+00+2.0d+00*(cos(Lconst*kx*(sqrt(3.0d+00))))&
                +4.0d+00*(Cos(Lconst*kx*(sqrt(3.0d+00))/2.0d+00))*(cos(Lconst*1.5d+00*ky))))&
         *((((sqrt(3.0d+00))*cos(theta))*(sin((sqrt(3.0d+00))*kx*Lconst)&
                                        +(sin((sqrt(3.0d+00))*kx*Lconst/2.0d+00))*(cos(1.5d+00*ky*Lconst))))&
          +(3.0d+00*(sin(theta))*(cos((sqrt(3.0d+00))*kx*Lconst/2.0d+00))*(sin(1.5d+00*ky*Lconst))))


!--------------------------------------

kx=mom*(cos(theta))
ky=mom*(sin(theta))

r3=sqrt(3.0d+00)
dis=sqrt(1.0d+00+4.0d+00*((cos(kx*Lconst*r3/2.0d+00))**2.0d+00)&
           +4.0d+00*(cos(kx*Lconst*r3/2.0d+00))*(cos(ky*Lconst*3.0d+00/2.0d+00)))


grouprad=(sgn)*((-0.5d+00)*t*(Lconst)*(1.0d+00/(dis))&
                         *((2.0d+00*(sin(kx*Lconst*r3))*(r3)*(cos(theta)))&
                        +(4.0d+00*(sin(kx*Lconst*r3/2.0d+00))*(cos(ky*Lconst*3.0d+00/2.0d+00))*(r3/2.0d+00)*(cos(theta)))&
                        +(4.0d+00*(cos(kx*Lconst*r3/2.0d+00))*(sin(ky*Lconst*3.0d+00/2.0d+00))*(3.0d+00/2.0d+00)*(sin(theta)))))


!--------------------------------------

groupangle=(groupvelky)*(cos(theta))-(groupvelkx)*(sin(theta))

!print*, 'difference-inside-route:',abs(sqrt(((grouprad)**2)+((groupangle)**2)))-abs(grouprad),blab,mom,theta/(2.0d+00*pi)

if(id==1)then

groupvel=grouprad

elseif(id==-1)then

groupvel=groupangle

elseif(id==0)then

groupvel=sqrt(((grouprad)**2)+((groupangle)**2))

endif

!----------------------------------------------------------------------
!----------------------------------------------------------------------
end function groupvel
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Subroutine to find the roots of the dipersion:
!----------------------------------------------------------------------
subroutine roots_dispersion(blab,angle,energy,nroots,KFS)

use commons

implicit none

common eF

!------------------------------------------------------------------ 
! Parameter declaration:
!------------------------------------------------------------------ 
DOUBLE PRECISION,INTENT(IN)::angle,energy
DOUBLE PRECISION,INTENT(OUT)::KFS
INTEGER,INTENT(OUT)::nroots
INTEGER,INTENT(IN)::blab

INTEGER::icommon,nchcommon,kint1,blabb
DOUBLE PRECISION::eF,scalecommon,kmax,kzerop,k1,groupvelroot&
                  ,kx,ky,anglered

!------------------------------------------------------------------ 
! parameters for zero_rc library
!------------------------------------------------------------------
DOUBLE PRECISION::zeromin,zeromax,zeroarg,zerovalue,zerofun&
                  ,zeromingr,zeromaxgr
INTEGER*4::zerostatus
!------------------------------------------------------------------
!------------------------------------------------------------------
! External functiions
!------------------------------------------------------------------
DOUBLE PRECISION:: fermi,dispersion,groupvel
! their corresponding variables:
DOUBLE PRECISION::mom,theta
!------------------------------------------------------------------
!------------------------------------------------------------------
INTEGER::icountersp,BZid,numsign,kint,nrootsmat
DOUBLE PRECISION::gvsign,kxr,kyr,KFSmin,KFSmax,kxbound,kybound,zeroer
DOUBLE PRECISION,DIMENSION(1:20)::KFSext,KFSmat
!------------------------------------------------------------------
!------------------------------------------------------------------
!------------------------------------------------------------------
kzerop=t*((10.0d+00)**(-8.0d+00))

zeroer = t*(10.0d+00**(-15.0d+00))
!----------------------------------------------------------------------            
!Initialization:
!----------------------------------------------------------------------

blabb=blab! important checking!!!!

CALL BZ_edge(angle,KFSmax)

!----------------------------------------------------------------------
! Initialization:
!----------------------------------------------------------------------

zeromin = kzerop
zeromax = KFSmax

zerostatus=0
nroots=-1
KFS=0.0d+00

icountersp=0

!----------------------------------------------------------------------
if(((dispersion(blabb,zeromin,angle,energy))*dispersion(blabb,zeromax,angle,energy))<0.0d+00)then

Do

call zero_rc ( zeromin, zeromax, zeroer, zeroarg, zerostatus, zerovalue )

if ( zerostatus == 0 ) then
exit 
end if

zerovalue=dispersion(blabb,zeroarg,angle,energy)


end do

kFS=zeroarg
nroots=1

elseif(abs(dispersion(blabb,zeromin,angle,energy))< zeroer)then

kFS=zeromin
nroots=1

elseif(abs(dispersion(blabb,zeromax,angle,energy))< zeroer)then

kFS=zeromax
nroots=1

elseif(((dispersion(blabb,zeromin,angle,energy))*(dispersion(blabb,zeromax,angle,energy)))>0.0d+00)then

nroots=0


endif

!----------------------------------------------------------------------

if(nroots>0 .and. abs(dispersion(blabb,kFS,angle,energy))>10.0d+00*zeroer)then
print*,'Note: roots of dispersion has not been found properly!',dispersion(blabb,kFS,angle,energy),angle/pi,eF
endif

if(nroots .LE. 0)then
KFS=0.0d+00
end if

!print*, 'Inside:', dispersion(blabb,KFS,angle,energy),blabb,KFS,angle

!-----------------------------------------------------------------------
end subroutine roots_dispersion
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!******* mapping integers to arrays:***********************
!-----------------------------------------------------------------------
subroutine onedtomultd(ndirection,i,n1,n2,n3,n4,p1,p2,p3)

use commons

implicit none

INTEGER::i,n1,n2,n3,n4,p1,p2,p3,k,n1v,n2v,n3v,n4v,p1v,p2v,p3v,ndirection
!------------------------------------------------------------------
k=0

Do n1v=1,nband
Do n2v=1,nband
Do n3v=1,nband
Do n4v=1,nband

Do p1v=1,npatch
Do p2v=1,npatch
Do p3v=1,npatch

if(ndirection==1)then  !1d to mult-d 

if(k<i)then
k=k+1

n1=n1v
n2=n2v
n3=n3v
n4=n4v

p1=p1v
p2=p2v
p3=p3v

end if

elseif(ndirection==-1)then  ! mult-d to 1d

k=k+1

if(n1v==n1 .and. n2v==n2 .and. n3v==n3 .and. n4v==n4 .and. p1v==p1 .and. p2v==p2 .and. p3v==p3)then

i=k

endif

endif !ndirection

end do !p3
end do !p2
end do !p1

end do !n4
end do !n3
end do !n2
end do !n1
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
end subroutine onedtomultd
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!*******Addition***********************
!-----------------------------------------------------------------------
subroutine Vec_add(angle1,KFS1,angle2,KFS2,angle3,KFS3,energy,anglep,radp)

use commons

implicit none

common eF

INTEGER::nroots1,nroots2,nroots3,kpint
DOUBLE PRECISION::eF,y,scalecommon
DOUBLE PRECISION::energy,angle1,angle2,angle3,KFS1,KFS2,KFS3,anglep&
                  ,kxp,kyp,radp
!----------------------------------------------------------------------
!----------------------------------------------------------------------

kxp=KFS1*cos(angle1)+KFS2*cos(angle2)-KFS3*cos(angle3)
kyp=KFS1*sin(angle1)+KFS2*sin(angle2)-KFS3*sin(angle3)

radp=sqrt(((kxp)**2)+((kyp)**2))
anglep=atan2(kyp,kxp)
anglep=(-SIGN(1.0d+00,anglep)+1.0d+00)*(pi)+anglep

!----------------------------------------------------------------------
!----------------------------------------------------------------------
end subroutine Vec_add
!----------------------------------------------------------------------
! K-delta
!----------------------------------------------------------------------
DOUBLE PRECISION function kdelta(n1,n2)

use commons

implicit none

Integer,intent(in)::n1,n2

if(n1==n2)then
kdelta=1.0d+00
else
kdelta=0.0d+00
endif

end function kdelta
!----------------------------------------------------------------------
!-----------------------------------------------------------------------
!*******patch index***********************
!-----------------------------------------------------------------------
subroutine patch_lab(flag,blab,theta,mom,nplab,anglek)

use commons

implicit none


CHARACTER*3, INTENT (IN)::flag
DOUBLE PRECISION, INTENT(IN)::theta,mom
INTEGER, INTENT(OUT)::nplab
DOUBLE PRECISION, INTENT(OUT)::anglek
INTEGER, INTENT (IN)::blab
DOUBLE PRECISION::kpx,kpy,kxprime,kyprime,base
!----------------------------------------------------------------------
!----------------------------------------------------------------------
kpy=mom*sin(theta)
kpx=mom*cos(theta)

!------------------------------------------------------------------
! Finding the corresponding patch index:
!------------------------------------------------------------------

anglek=atan2(kpy,kpx)
anglek=(-SIGN(1.0d+00,anglek)+1.0d+00)*(pi)+anglek

nplab=min(npatch,(1+int((anglek*npatch)/(2.0d+00*pi))))

!------------------------------------------------------------------
! 'complex' patching scheme:
!------------------------------------------------------------------
if (flag=='com')then 

if(kpx .GE. 0.0d+00 .and. kpy .GE. 0.0d+00)then

if(kpy>pi-kpx)then

kyprime=pi-kpx
kxprime=pi-kpy
base=0.0d+00

anglek=base+atan2(kyprime,kxprime)
anglek=(-SIGN(1.0d+00,anglek)+1.0d+00)*(pi)+anglek
nplab=min(npatch/4,(1+int((anglek*npatch)/(2.0d+00*pi))))

endif

elseif(kpx<0.0d+00 .and. kpy .GE. 0.0d+00)then

if(kpy>pi+kpx)then

kyprime=pi-kpy
kxprime=pi-abs(kpx) 
base=pi/2.0d+00

anglek=base+atan2(kyprime,kxprime)
anglek=(-SIGN(1.0d+00,anglek)+1.0d+00)*(pi)+anglek
nplab=min(npatch/2,(1+int((anglek*npatch)/(2.0d+00*pi))))

endif

elseif(kpx<0.0d+00 .and. kpy<0.0d+00)then

if(kpy < -pi-kpx)then

kyprime=pi-abs(kpx)
kxprime=pi-abs(kpy)
base=pi

anglek=base+atan2(kyprime,kxprime)
anglek=(-SIGN(1.0d+00,anglek)+1.0d+00)*(pi)+anglek
nplab=min(3*npatch/4,(1+int((anglek*npatch)/(2.0d+00*pi))))

endif

elseif(kpx .GE. 0.0d+00 .and. kpy<0.0d+00)then

if(kpy < -pi+kpx)then

kyprime=pi-abs(kpy)
kxprime=pi-kpx
base=1.5d+00*pi

anglek=base+atan2(kyprime,kxprime)
anglek=(-SIGN(1.0d+00,anglek)+1.0d+00)*(pi)+anglek
nplab=min(npatch,(1+int((anglek*npatch)/(2.0d+00*pi))))

endif

endif 

endif

!----------------------------------------------------------------------
!----------------------------------------------------------------------
end subroutine patch_lab
!----------------------------------------------------------------------
! Matrix elements of the spin-orbit Unitary matrix!
!----------------------------------------------------------------------
SUBROUTINE SPIN_ORBIT(kx,ky,Emat,Umat,Umatdag,HSO0)
!----------------------------------------------------------------------
use commons
implicit none
common eF
!----------------------------------------------------------------------
DOUBLE PRECISION, INTENT(IN)::kx,ky
COMPLEX*16, DIMENSION(1:nband,1:nband), INTENT(OUT)::Umat,Umatdag,HSO0
DOUBLE PRECISION, DIMENSION(1:nband), INTENT(OUT)::Emat
COMPLEX*16::VSO12,VSO21
DOUBLE PRECISION::VSOx,VSOy,VSOz,ESOp,ESOn,disp&
                  ,krad,kangle,disp1,disp2,eF,dispersion
!----------------------------------------------------------------------
!------------------------------------------------------------------
! parameters for the Eispack library:
!------------------------------------------------------------------
INTEGER*4::Mdim,EVerr
REAL*8,ALLOCATABLE,DIMENSION(:,:)::Mmatre,Mmatim,MEVecre,MEVecim
REAL*8,ALLOCATABLE,DIMENSION(:)::MEValre,MEValim
LOGICAL::EVflag
!------------------------------------------------------------------
INTEGER::n1,n2
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Type of the spin-orbi coupling:
!----------------------------------------------------------------------
! Rashba:
!----------------------------------------------------------------------
!VSO12=lamSO*((ic)*((VSOx12+(CONJG(VSOx21)))*(sin(kx)))+((VSOy12+(CONJG(VSOy21)))*(sin(ky))))
!VSO21=lamSO*((ic)*((VSOx21+(CONJG(VSOx12)))*(sin(kx)))+((VSOy21+(CONJG(VSOy12)))*(sin(ky))))

VSOx=LamSO*(sin(ky))!real(VSO12)
VSOy=-LamSO*(sin(kx))!-aimag(VSO12)
VSOz=Bfieldz/2.0d+00!0.0d+00

ESOp=sqrt(((VSOx)**2)+((VSOy)**2)+((VSOz)**2))
ESOn=-sqrt(((VSOx)**2)+((VSOy)**2)+((VSOz)**2))

Emat(1)=ESOn!ESOp
Emat(2)=ESOp!ESOn
!----------------------------------------------------------------------
!----------------------------------------------------------------------

Umat(1,2)=cmplx(1.0d+00/(sqrt(2.0d+00)),0.0d+00)
Umat(2,2)=(1.0d+00/(sqrt(2.0d+00)))*(VSOx+ic*(VSOy))/(sqrt(((VSOx)**2)+((VSOy)**2)))

Umat(1,1)=cmplx(1.0d+00/(sqrt(2.0d+00)),0.0d+00)
Umat(2,1)=(-1.0d+00/(sqrt(2.0d+00)))*(VSOx+ic*(VSOy))/(sqrt(((VSOx)**2)+((VSOy)**2)))


Umatdag(1,1)=CONJG(Umat(1,1))
Umatdag(1,2)=CONJG(Umat(2,1))

Umatdag(2,1)=CONJG(Umat(1,2))
Umatdag(2,2)=CONJG(Umat(2,2))

!----------------------------------------------------------------------
!----------------------------------------------------------------------
HSO0(1,1)=cmplx(-2.0d+00*t*(cos(kx)+cos(ky)),0.0d+00)
HSO0(2,2)=cmplx(-2.0d+00*t*(cos(kx)+cos(ky)),0.0d+00)
HSO0(1,2)=VSOx-ic*VSOy
HSO0(2,1)=CONJG(HSO0(1,2))
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Finding the elements of Umat using Eispack library:
!----------------------------------------------------------------------
Mdim=2
ALLOCATE(Mmatre(1:Mdim,1:Mdim),Mmatim(1:Mdim,1:Mdim),MEVecre(1:Mdim,1:Mdim),MEVecim(1:Mdim,1:Mdim)&
        ,MEValre(1:Mdim),MEValim(1:Mdim))


!disp=-2.0d+00*t*(cos(kx)+cos(ky))

!--------------------------------------------------
krad=sqrt((kx**2)+(ky**2))

kangle=atan2(ky,kx)
kangle=(-SIGN(1.0d+00,kangle)+1.0d+00)*(pi)+kangle

!disp1=dispersion(1,krad,kangle,eF) ! maybe better to define base dispersion
!disp2=dispersion(2,krad,kangle,eF)

disp=-2.0d+00*t*(cos(krad*(cos(kangle)))+cos(krad*(sin(kangle))))-eF&
           -4.0d+00*tprime*(cos(krad*(cos(kangle))))*(cos(krad*(sin(kangle))))

!--------------------------------------------------

Mmatre(1,1)=disp+VSOz
Mmatim(1,1)=0.0d+00

Mmatre(1,2)=VSOx+(Bfieldx)
Mmatim(1,2)=-VSOy-(Bfieldy)

Mmatre(2,1)=VSOx+(Bfieldx)
Mmatim(2,1)=VSOy+(Bfieldy)

Mmatre(2,2)=disp-VSOz
Mmatim(2,2)=0.0d+00

EVflag=.true.

CALL ch(Mdim,Mmatre,Mmatim,MEValre,EVflag,MEVecre,MEVecim,EVerr)

!CALL cg_lr (Mdim,Mmatre,Mmatim,MEValre,MEValim,EVflag,MEVecre,MEVecim,EVerr)


Do n1=1,Mdim
Do n2=1,Mdim

Umat(n1,n2)=cmplx(MEVecre(n1,n2),MEVecim(n1,n2))

end do  !n2

Emat(n1)=MEvalre(n1)-disp

end do  !n1



DEALLOCATE(Mmatre,Mmatim,MEVecre,MEVecim,MEValre,MEValim)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
END SUBROUTINE SPIN_ORBIT
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
SUBROUTINE Hexagonal(kx,ky,id)
!----------------------------------------------------------------------
use commons
implicit none
common eF
!----------------------------------------------------------------------
DOUBLE PRECISION, INTENT(IN)::kx,ky
INTEGER, INTENT(OUT)::id
DOUBLE PRECISION::angle,eF,KFSmax,kxmax,kymax!,zeroHex,zeroHexedge
!----------------------------------------------------------------------
!zeroHex = t*(10.0d+00**(-12.0d+00))
!zeroHexedge= t*(10.0d+00**(-12.0d+00))

angle=atan2(ky,kx)
angle=(-SIGN(1.0d+00,angle)+1.0d+00)*(pi)+angle

if(abs((2.0d+00*pi/npatch)*int(angle/(2.0d+00*pi/npatch))-angle)<patcherr1)then

angle=angle+patcherr1

if(angle>2.0d+00*pi)then
angle=angle-2.0d+00*pi
endif

endif

CALL BZ_edge(angle,KFSmax)

kxmax=KFSmax*cos(angle)
kymax=KFSmax*sin(angle)

!----------------------------------------------------------------------
id=0
!----------------------------------------------------------------------
if(abs(kx) .LE. BZconst .or. abs(kx-zeroHex) .LE. BZconst .or. abs(kx+zeroHex) .LE. BZconst)then

if(angle < (pi/3.0d+00))then
if (ky .LE.(sqrt(3.0d+00))*(BZconst-kx) .or. ky-zeroHex .LE.(sqrt(3.0d+00))*(BZconst-kx) &
       .or. ky+zeroHex .LE.(sqrt(3.0d+00))*(BZconst-kx))then
id=1
endif
elseif(angle .GE.(pi/3.0d+00) .and. angle < (2.0d+00*pi/3.0d+00))then
if(ky .LE. ((sqrt(3.0d+00))/2.0d+00)*BZconst .or. ky-zeroHex .LE. ((sqrt(3.0d+00))/2.0d+00)*BZconst &
      .or. ky+zeroHex .LE. ((sqrt(3.0d+00))/2.0d+00)*BZconst)then
id=1
endif
elseif(angle .GE.(2.0d+00*pi/3.0d+00) .and. angle < (pi))then
if(ky .LE. ((sqrt(3.0d+00))*(kx+BZconst)) .or. ky-zeroHex .LE. ((sqrt(3.0d+00))*(kx+BZconst))&
   .or. ky+zeroHex .LE. ((sqrt(3.0d+00))*(kx+BZconst)))then
id=1
endif
elseif(angle .GE.(pi) .and. angle < (4.0d+00*pi/3.0d+00))then
if(ky .GE. (-(sqrt(3.0d+00))*(kx+BZconst)) .or. ky-zeroHex .GE. (-(sqrt(3.0d+00))*(kx+BZconst))&
      .or. ky+zeroHex .GE. (-(sqrt(3.0d+00))*(kx+BZconst)))then
id=1
endif
elseif(angle .GE.(4.0d+00*pi/3.0d+00) .and. angle < (5.0d+00*pi/3.0d+00))then
if(ky .GE. -((sqrt(3.0d+00))/2.0d+00)*BZconst .or. ky-zeroHex .GE. -((sqrt(3.0d+00))/2.0d+00)*BZconst &
      .or. ky+zeroHex .GE. -((sqrt(3.0d+00))/2.0d+00)*BZconst)then
id=1
endif
elseif(angle .GE.(5.0d+00*pi/3.0d+00) .and. angle .LE. (2.0d+00*pi))then
if(ky .GE. ((sqrt(3.0d+00))*(kx-BZconst)) .or. ky-zeroHex .GE. ((sqrt(3.0d+00))*(kx-BZconst))&
      .or. ky+zeroHex .GE. ((sqrt(3.0d+00))*(kx-BZconst)))then
id=1
endif
endif


else

id=0

endif
!----------------------------------------------------------------------
if(abs(kx-kxmax)<zeroHexedge .and. abs(ky-kymax)<zeroHexedge)then

id=2

endif
!----------------------------------------------------------------------

END SUBROUTINE Hexagonal
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! A subroutine for folding for Hexagonal lattice:
!----------------------------------------------------------------------
!----------------------------------------------------------------------
SUBROUTINE Folding(kx,ky,kxr,kyr)
!----------------------------------------------------------------------
use commons
implicit none
common eF
!----------------------------------------------------------------------
DOUBLE PRECISION,INTENT(IN)::kx,ky
DOUBLE PRECISION,INTENT(OUT)::kxr,kyr

DOUBLE PRECISION,DIMENSION(1:6,1:2)::BZV
DOUBLE PRECISION::KFS,angle,eF,kpxri,kpyri,KFSmax,kxmax,kymax
INTEGER::BZcount,BZid,nkxint,nkyint,neiint,nneimax,neicount,n1,n2,n3 
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::BZVnei,BZVneibe
!----------------------------------------------------------------------
angle=atan2(ky,kx)
angle=(-SIGN(1.0d+00,angle)+1.0d+00)*(pi)+angle

KFS=sqrt(((kx)**2)+((ky)**2))

!----------------------------------------------------------------------
! Restriction to the BZ:
!----------------------------------------------------------------------

BZV(1,1)=BZconst*(1.5d+00)
BZV(1,2)=BZconst*(sqrt(3.0d+00))/2.0d+00

BZV(2,1)=0.0d+00
BZV(2,2)=BZconst*(sqrt(3.0d+00))

BZV(3,1)=-BZV(1,1)
BZV(3,2)=BZV(1,2)

BZV(4,1:2)=-BZV(1,1:2)

BZV(5,1:2)=-BZV(2,1:2)

BZV(6,1)=BZV(1,1)
BZV(6,2)=-BZV(1,2)


!print*,'Check-Folding-begin:',kx/BZconst,ky/BZconst


CALL Hexagonal(kx,ky,BZid)

if(BZid==2)then

CALL BZ_edge(angle,KFSmax)

kxmax=KFSmax*cos(angle)
kymax=KFSmax*sin(angle)

kxr=kxmax
kyr=kymax

elseif(BZid==1)then

kxr=kx
kyr=ky

else


BZid=0
neiint=1
BZcount=1

Do

if(BZcount==1)then

nneimax=6*(BZcount)

if(neiint==1)then

Allocate(BZVnei(1:nneimax,1:2))

BZVnei(1:nneimax,1:2)=BZV(1:nneimax,1:2)

endif

endif

kpxri=kx-BZVnei(neiint,1)
kpyri=ky-BZVnei(neiint,2)

CALL Hexagonal(kpxri,kpyri,BZid)

!print*,'Check-Folding:',BZcount,neiint,kpxri/BZconst,kpyri/BZconst,BZid

!write(52,*)BZcount,neiint,kpxri/BZconst,kpyri/BZconst,BZid

if(BZid==1)then

kxr=kpxri
kyr=kpyri

!print*,kxint,kyint,'BZcount=',BZcount,'neiint=',neiint,kpx/pi,kpy/pi

Deallocate(BZVnei)

exit
endif


if(neiint<nneimax)then
neiint=neiint+1
else

Allocate(BZVneibe(1:nneimax,1:2))
BZVneibe(1:nneimax,1:2)=BZVnei(1:nneimax,1:2)

Deallocate(BZVnei)

BZcount=BZcount+1

nneimax=6*(BZcount)

Allocate(BZVnei(1:nneimax,1:2))

neicount=1

Do n1=1,6

if(n1==6)then
n2=1
else
n2=n1+1
endif

Do n3=0,BZcount-1

BZVnei(neicount,1:2)=(BZcount-n3)*BZV(n1,1:2)+(n3)*BZV(n2,1:2)       

neicount=neicount+1

end do !n3

end do !n1

Deallocate(BZVneibe)

neiint=1

endif

if(BZcount>2)then
print*,'too far away from origin!',kx/BZconst,ky/BZconst
Deallocate(BZVnei)
exit
endif


end do

endif

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------

END SUBROUTINE Folding
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! A Subroutine to find intervals of chaning sign!
!----------------------------------------------------------------------
!----------------------------------------------------------------------
SUBROUTINE Extermum(l1,energy,angle,numsign,KFSext)
!----------------------------------------------------------------------
use commons
implicit none
common eF
!----------------------------------------------------------------------
INTEGER, INTENT(IN)::l1
DOUBLE PRECISION, INTENT(IN)::angle,energy
INTEGER, INTENT(OUT)::numsign
DOUBLE PRECISION,DIMENSION(1:20),INTENT(OUT)::KFSext
!----------------------------------------------------------------------
INTEGER:: kint,BZid
DOUBLE PRECISION::eF,dispersion,dispbefore,KFS,kx,ky,KFSmax,KFSedge
!----------------------------------------------------------------------
! initialization:

numsign=0
KFSext(1:20)=0.0d+00

!CALL  BZ_edge(angle,KFSedge)

KFSedge=BZconst

Do kint=1,501

KFS=0.0d+00+((KFSedge)*(kint-1.0d+00)/(500.0d+00))

kx=KFS*(cos(angle))
ky=KFS*(sin(angle))

CALL Hexagonal(kx,ky,BZid)

if(BZid>0)then

if(kint==1)then
dispbefore=dispersion(l1,KFS,angle,eF)
endif

if(dispersion(l1,KFS,angle,eF).NE. SIGN(dispersion(l1,KFS,angle,eF),dispbefore))then

dispbefore=dispersion(l1,KFS,angle,eF)
numsign=numsign+1

KFSext(numsign)=KFS

endif  


!write(54,*)angle/pi,KFS/BZconst,numsign,KFSmax,dispersion(l1,KFS,angle,eF)

endif !BZid

if(numsign>20)then
print*,'Error!:', 'numsign-max should be increased!'
endif

!----------------------------------------------------------------------
end do !kint
!----------------------------------------------------------------------

END SUBROUTINE Extermum
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! A subroutine for folding for Hexagonal lattice:
!----------------------------------------------------------------------
!----------------------------------------------------------------------
SUBROUTINE aFolding(kx,ky,kxr,kyr)
!----------------------------------------------------------------------
use commons
implicit none
common eF
!----------------------------------------------------------------------
DOUBLE PRECISION,INTENT(IN)::kx,ky
DOUBLE PRECISION,INTENT(OUT)::kxr,kyr

DOUBLE PRECISION,DIMENSION(1:6,1:2)::BZV
DOUBLE PRECISION::KFS,angle,eF,kxri,kyri,sgn,KFSr,kxd,kyd,deltax,deltay&
                  ,bvecrad
INTEGER::BZcount,BZid,nkxint,nkyint,neiint,nneimax,neicount,n1,n2,n3,n4 
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::BZVnei,BZVneibe


DOUBLE PRECISION,DIMENSION(1:2,1:2)::bvec
!----------------------------------------------------------------------
angle=atan2(ky,kx)
angle=(-SIGN(1.0d+00,angle)+1.0d+00)*(pi)+angle

KFS=sqrt(((kx)**2)+((ky)**2))

kxd=kx
kyd=ky


kxri=kx
kyri=ky   
!----------------------------------------------------------------------
! Restriction to the BZ:
!----------------------------------------------------------------------

bvec(1,1)=BZconst*1.5d+00
bvec(1,2)=BZconst*((sqrt(3.0d+00))/2.0d+00)

bvec(2,1)=-bvec(1,1)
bvec(2,2)=bvec(1,2)

bvecrad=sqrt(((bvec(1,1))**2)+((bvec(1,2))**2))

n3=1+int(KFS/bvecrad)

n4=0

print*,'check radius:',bvecrad/BZconst,KFS/bvecrad,n3

BZcount=0
n1=1
n2=1

print*,'inside-routine-begin',kx/BZconst,ky/BZconst,KFS/BZconst

Do

if(BZcount>10)then
print*,'careful too far away from Gamma point!'
exit
endif


if(n2==1)then
sgn=1.0d+00
elseif(n2==2)then
sgn=-1.0d+00
endif

kxri=kxd+(n3)*(sgn)*bvec(n1,1)
kyri=kyd+(n3)*(sgn)*bvec(n1,2)

CALL Hexagonal(kxri,kyri,BZid)

if(BZid==1)then
kxr=kxri
kyr=kyri
exit

else

KFSr=sqrt(kxri**2+kyri**2)

print*,'checking-inside routine:',BZcount,kxri/BZconst,kyri/BZconst,BZid,KFSr/BZconst,n1,n2,n4

if(n4==1 .and. KFSr.LE. n3*bvecrad)then
kxd=kxri
kyd=kyri
endif

!KFS=sqrt(kxd**2+kyd**2)
if(n1<2)then
n1=n1+1

print*,n1,n2
elseif(n2<2)then
n2=n2+1
n1=1
elseif(n1==2 .and. n2==2)then
n4=1
n1=1
n2=1
endif

!endif

BZcount=BZcount+1


endif

if(BZid==0)then
kxr=kxd
kyr=kyd
endif

end do
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------

END SUBROUTINE aFolding
!----------------------------------------------------------------------
!-----------------------------------------------------------------------
!*******patch index***********************
!-----------------------------------------------------------------------
subroutine patch_lab_Ge(flag,theta,mom,nplab,anglek)

use commons

implicit none


CHARACTER*3, INTENT (IN)::flag
DOUBLE PRECISION, INTENT(IN)::theta,mom
INTEGER, INTENT(OUT)::nplab
DOUBLE PRECISION, INTENT(OUT)::anglek
!INTEGER, INTENT (IN)::blab
DOUBLE PRECISION::kx,ky,kxprime,kyprime,base,KFSmax,kxmax,kymax!,patcherr1,patcherr2
INTEGER::nplab0,nplabm,nplabp,BZid
!----------------------------------------------------------------------
!patcherr1=10.0d+00**(-15.0d+00)
!patcherr2=10.0d+00**(-14.0d+00)
!----------------------------------------------------------------------
!------------------------------------------------------------------
! Finding the corresponding patch index:
!------------------------------------------------------------------
ky=mom*sin(theta)
kx=mom*cos(theta)

CALL Hexagonal(kx,ky,BZid)

!print*,'inside-patch-lab-Ge-step1:',theta/(2.0d+00*pi),BZid

if(BZid==2)then

CALL BZ_edge(theta,KFSmax)

kxmax=KFSmax*cos(theta)
kymax=KFSmax*sin(theta)

ky=kymax
kx=kxmax

anglek=theta

elseif(BZid==0)then

CALL Folding(kx,ky,kx,ky)

anglek=atan2(ky,kx)
anglek=(-SIGN(1.0d+00,anglek)+1.0d+00)*(pi)+anglek

elseif(BZid==1)then

anglek=theta

end if

!----------------------------------------------

if(anglek==2.0d+00*pi)then
anglek=0.0d+00
elseif(anglek<0.0d+00 .or. anglek>2.0d+00*pi)then
print*, 'from patch_Lab_Ge: UPS:angle ouside [0.2*pi]'
endif

!----------------------------------------------

nplab0=min(npatch,(1+int(((anglek+patcherr1)*npatch)/(2.0d+00*pi))))

!----------------------------------------------
ky=(mom-patcherr2)*sin(theta)
kx=(mom-patcherr2)*cos(theta)

CALL Hexagonal(kx,ky,BZid)

!print*,'inside-patch-lab-Ge-step1:',theta/(2.0d+00*pi),BZid

if(BZid==2)then

CALL BZ_edge(theta,KFSmax)

kxmax=KFSmax*cos(theta)
kymax=KFSmax*sin(theta)

ky=kymax
kx=kxmax

anglek=theta

elseif(BZid==0)then

CALL Folding(kx,ky,kx,ky)

anglek=atan2(ky,kx)
anglek=(-SIGN(1.0d+00,anglek)+1.0d+00)*(pi)+anglek

elseif(BZid==1)then

anglek=theta

end if

!----------------------------------------------

if(anglek==2.0d+00*pi)then
anglek=0.0d+00
elseif(anglek<0.0d+00 .or. anglek>2.0d+00*pi)then
!print*, 'from patch_Lab_Ge: UPS:angle ouside [0.2*pi]'
endif

!----------------------------------------------

nplabm=min(npatch,(1+int(((anglek+patcherr1)*npatch)/(2.0d+00*pi))))

!----------------------------------------------
ky=(mom+patcherr2)*sin(theta)
kx=(mom+patcherr2)*cos(theta)

CALL Hexagonal(kx,ky,BZid)

!print*,'inside-patch-lab-Ge-step1:',theta/(2.0d+00*pi),BZid

if(BZid==2)then

CALL BZ_edge(theta,KFSmax)

kxmax=KFSmax*cos(theta)
kymax=KFSmax*sin(theta)

ky=kymax
kx=kxmax

anglek=theta

elseif(BZid==0)then

CALL Folding(kx,ky,kx,ky)

anglek=atan2(ky,kx)
anglek=(-SIGN(1.0d+00,anglek)+1.0d+00)*(pi)+anglek

elseif(BZid==1)then

anglek=theta

end if

!----------------------------------------------

if(anglek==2.0d+00*pi)then
anglek=0.0d+00
elseif(anglek<0.0d+00 .or. anglek>2.0d+00*pi)then
!print*, 'from patch_Lab_Ge: UPS:angle ouside [0.2*pi]'
endif

!----------------------------------------------

nplabp=min(npatch,(1+int(((anglek+patcherr1)*npatch)/(2.0d+00*pi))))

!----------------------------------------------
if(nplabm .ne. nplab0 .or. nplabp .ne. nplab0)then

print*,'ambiguity',nplab0,nplabm,nplabp

endif
!----------------------------------------------

nplab=min(nplab0,nplabp,nplabm)

!print*,'inside-patch-lab-Ge-step2:',nplab,nplabm,nplabp,anglek/(2.0d+00*pi),BZid,(anglek*npatch)/(2.0d+00*pi)

!------------------------------------------------------------------
! 'complex' patching scheme:
!------------------------------------------------------------------
if (flag=='com')then 

if(theta .LE. (pi/6.0d+00) .and. kx> (3.0d+00*BZconst/4.0d+00))then
kyprime=ky
kxprime=3.0d+00*(BZconst-kx)
anglek=atan2(kyprime,kxprime)

anglek=(-SIGN(1.0d+00,anglek)+1.0d+00)*(pi)+anglek
nplab=min(npatch/12,(1+int((anglek*npatch)/(2.0d+00*pi))))

elseif(theta> (pi/6.0d+00) .and. theta .LE. (pi/3.0d+00) .and. &
       ky .GE. ((-(sqrt(3.0d+00))/3.0d+00)*(kx-(1.5d+00*BZconst))) )then

kxprime=((sqrt(3.0d+00))/2.0d+00)*BZconst-ky
kyprime=(-(sqrt(3.0d+00))/3.0d+00)*ky+BZconst-kx

anglek=atan2(kyprime,kxprime)

anglek=(pi/3.0d+00)-atan2((tan((pi/3.0d+00)-anglek))/3.0d+00,1.0d+00)

anglek=(-SIGN(1.0d+00,anglek)+1.0d+00)*(pi)+anglek
nplab=min(npatch/6,(1+int((anglek*npatch)/(2.0d+00*pi))))

elseif(theta > (pi/3.0d+00) .and. theta .LE. (pi/2.0d+00) .and. &
       ky .GE. ((-(sqrt(3.0d+00))/3.0d+00)*(kx-(1.5d+00*BZconst))) )then

kyprime=((sqrt(3.0d+00))*BZconst/2.0d+00)-ky
kxprime=(BZconst/2.0d+00)-kx

anglek=(pi/3.0d+00)-atan2(kyprime,kxprime)

anglek=(pi/3.0d+00)+atan2(((tan(anglek))/3.0d+00),1.0d+00)
anglek=(-SIGN(1.0d+00,anglek)+1.0d+00)*(pi)+anglek
nplab=min(npatch/4,(1+int((anglek*npatch)/(2.0d+00*pi))))


elseif(theta> (pi/2.0d+00) .and. theta .LE. (2.0d+00*pi/3.0d+00) .and. &
       ky .GE. (((sqrt(3.0d+00))/3.0d+00)*(kx+(1.5d+00*BZconst))) )then

kyprime=((sqrt(3.0d+00))*BZconst/2.0d+00)-ky
kxprime=(BZconst/2.0d+00)-abs(kx)

anglek=atan2(kyprime,kxprime)

anglek=(2.0d+00*pi/3.0d+00)-atan2(((tan((pi/3.0d+00)-anglek))/3.0d+00),1.0d+00)
anglek=(-SIGN(1.0d+00,anglek)+1.0d+00)*(pi)+anglek
nplab=min(npatch/3,(1+int((anglek*npatch)/(2.0d+00*pi))))


elseif(theta> (2.0d+00*pi/3.0d+00) .and. theta .LE. (5.0d+00*pi/6.0d+00) .and. &
       ky .GE. (((sqrt(3.0d+00))/3.0d+00)*(kx+(1.5d+00*BZconst))) )then

kxprime=((sqrt(3.0d+00))*BZconst/2.0d+00)-ky
kyprime=abs((ky/(sqrt(3.0d+00)))-BZconst)-abs(kx)

anglek=(pi/3.0d+00)-atan2(kyprime,kxprime)

anglek=(2.0d+00*pi/3.0d+00)+atan2(((tan(anglek))/3.0d+00),1.0d+00)
anglek=(-SIGN(1.0d+00,anglek)+1.0d+00)*(pi)+anglek
nplab=min(5*npatch/12,(1+int((anglek*npatch)/(2.0d+00*pi))))


elseif(theta> (5.0d+00*pi/6.0d+00) .and. theta .LE. (pi) .and. &
       abs(kx) .GE. (3.0d+00*BZconst/4.0d+00) )then

kxprime=3.0d+00*((BZconst)-abs(kx))
kyprime=ky

!anglek=(pi/3.0d+00)-atan2(kyprime,kxprime)

anglek=pi-atan2(kyprime,kxprime)
anglek=(-SIGN(1.0d+00,anglek)+1.0d+00)*(pi)+anglek
nplab=min(6*npatch/12,(1+int((anglek*npatch)/(2.0d+00*pi))))


!----------------------------------------------------------------------
! Lower-half
!----------------------------------------------------------------------
elseif(theta> (pi) .and. theta .LE. (13.0d+00*pi/12.0d+00) .and. &
       abs(kx) .GE. (3.0d+00*BZconst/4.0d+00) )then

kxprime=3.0d+00*((BZconst)-abs(kx))
kyprime=abs(ky)

anglek=(pi)+atan2(kyprime,kxprime)
anglek=(-SIGN(1.0d+00,anglek)+1.0d+00)*(pi)+anglek
nplab=min(7*npatch/12,(1+int((anglek*npatch)/(2.0d+00*pi))))


elseif(theta> (7.0d+00*pi/6.0d+00) .and. theta .LE. (8.0d+00*pi/6.0d+00) .and. &
       abs(ky) .GE. abs((-(sqrt(3.0d+00))/3.0d+00)*(kx+(1.5d+00*BZconst))) )then

kxprime=(BZconst*(sqrt(3.0d+00))/2.0d+00)-abs(ky)
kyprime=abs((((-sqrt(3.0d+00))/3.0d+00)*ky)-BZconst)-abs(kx)

anglek=atan2(kyprime,kxprime)
anglek=(4.0d+00*pi/3.0d+00)-atan2(((tan((pi/3.0d+00)-anglek))/3.0d+00),1.0d+00)

anglek=(-SIGN(1.0d+00,anglek)+1.0d+00)*(pi)+anglek
nplab=min(8*npatch/12,(1+int((anglek*npatch)/(2.0d+00*pi))))

elseif(theta> (8.0d+00*pi/6.0d+00) .and. theta .LE. (9.0d+00*pi/6.0d+00) .and. &   !checked!
       abs(ky) .GE. abs((-(sqrt(3.0d+00))/3.0d+00)*(kx+(1.5d+00*BZconst))) )then

kyprime=(BZconst*(sqrt(3.0d+00))/2.0d+00)-abs(ky)
kxprime=(BZconst/2.0d+00)-abs(kx)

anglek=atan2(kyprime,kxprime)
anglek=(4.0d+00*pi/3.0d+00)+atan2(((tan((pi/3.0d+00)-anglek))/3.0d+00),1.0d+00)

anglek=(-SIGN(1.0d+00,anglek)+1.0d+00)*(pi)+anglek
nplab=min(9*npatch/12,(1+int((anglek*npatch)/(2.0d+00*pi))))

elseif(theta> (9.0d+00*pi/6.0d+00) .and. theta .LE. (10.0d+00*pi/6.0d+00) .and. &    !checked
       abs(ky) .GE. abs(((sqrt(3.0d+00))/3.0d+00)*(kx-(1.5d+00*BZconst))) )then

kyprime=(BZconst*(sqrt(3.0d+00))/2.0d+00)-abs(ky)
kxprime=(BZconst/2.0d+00)-abs(kx)

anglek=atan2(kyprime,kxprime)
anglek=(5.0d+00*pi/3.0d+00)-atan2(((tan((pi/3.0d+00)-anglek))/3.0d+00),1.0d+00)

anglek=(-SIGN(1.0d+00,anglek)+1.0d+00)*(pi)+anglek
nplab=min(10*npatch/12,(1+int((anglek*npatch)/(2.0d+00*pi))))


elseif(theta> (10.0d+00*pi/6.0d+00) .and. theta .LE. (11.0d+00*pi/6.0d+00) .and. &  !checked!
       abs(ky) .GE. abs(((sqrt(3.0d+00))/3.0d+00)*(kx-(1.5d+00*BZconst))) )then

kxprime=(BZconst*(sqrt(3.0d+00))/2.0d+00)-abs(ky)
kyprime=abs((ky/(sqrt(3.0d+00)))+BZconst)-abs(kx)

anglek=atan2(kyprime,kxprime)
anglek=(5.0d+00*pi/3.0d+00)+atan2(((tan((pi/3.0d+00)-anglek))/3.0d+00),1.0d+00)

anglek=(-SIGN(1.0d+00,anglek)+1.0d+00)*(pi)+anglek
nplab=min(11*npatch/12,(1+int((anglek*npatch)/(2.0d+00*pi))))

elseif(theta> (11.0d+00*pi/6.0d+00) .and. theta .LE. (12.0d+00*pi/6.0d+00) .and. &   !checked!
       abs(kx) .GE. (3.0d+00*BZconst/4.0d+00) )then

kyprime=abs(ky)
kxprime=(BZconst)-abs(kx)

anglek=atan2(kyprime,kxprime)
anglek=(2.0d+00*pi)-atan2(((tan(anglek))/3.0d+00),1.0d+00)

anglek=(-SIGN(1.0d+00,anglek)+1.0d+00)*(pi)+anglek
nplab=min(12*npatch/12,(1+int((anglek*npatch)/(2.0d+00*pi))))

endif

endif
!----------------------------------------------------------------------
!----------------------------------------------------------------------
end subroutine patch_lab_Ge
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Matrix elements of the spin-orbit Unitary matrix!
!----------------------------------------------------------------------
SUBROUTINE Unitary(kx,ky,Umat)
!----------------------------------------------------------------------
use commons
implicit none
common eF
!----------------------------------------------------------------------
DOUBLE PRECISION, INTENT(IN)::kx,ky
COMPLEX*16, DIMENSION(1:2,1:2), INTENT(OUT)::Umat
COMPLEX*16::VSO12,VSO21
DOUBLE PRECISION::VSOx,VSOy,VSOz,ESOp,ESOn,disp&
                  ,krad,kangle,disp1,disp2,eF,dispersion
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! parameters for the Eispack library:
!----------------------------------------------------------------------
INTEGER*4::Mdim,EVerr
REAL*8,ALLOCATABLE,DIMENSION(:,:)::Mmatre,Mmatim,MEVecre,MEVecim
REAL*8,ALLOCATABLE,DIMENSION(:)::MEValre,MEValim
LOGICAL::EVflag
!----------------------------------------------------------------------
INTEGER::n1,n2
COMPLEX*16,DIMENSION(1:2,1:2)::Mmattype
COMPLEX*16::fun
DOUBLE PRECISION::phase
!----------------------------------------------------------------------

!fun=cos(ky*Lconst)-ic*sin(ky*Lconst)&
!       +(cos(ky*Lconst/2.0d+00)+ic*sin(ky*Lconst/2.0d+00))*2.0d+00*cos(kx*Lconst*(sqrt(3.0d+00))/2.0d+00)

!phase=atan2(aimag(fun),real(fun)) 

!Mmattype(1,1)=cmplx(-eF,0.0d+00)
!Mmattype(2,2)=cmplx(-eF,0.0d+00)

!Mmattype(1,2)=fun
!Mmattype(2,1)=CONJG(fun)

!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Using Eispack library for the Eigenvalue problem:
!----------------------------------------------------------------------
!Mdim=2
!ALLOCATE(Mmatre(1:Mdim,1:Mdim),Mmatim(1:Mdim,1:Mdim),MEVecre(1:Mdim,1:Mdim),MEVecim(1:Mdim,1:Mdim)&
!        ,MEValre(1:Mdim),MEValim(1:Mdim))

!----------------------------------------------------------------------
!Mmatre(1,1)=real(Mmattype(1,1))
!Mmatim(1,1)=aimag(Mmattype(1,1))

!Mmatre(1,2)=real(Mmattype(1,2))
!Mmatim(1,2)=aimag(Mmattype(1,2))

!Mmatre(2,1)=real(Mmattype(2,1))
!Mmatim(2,1)=aimag(Mmattype(2,1))

!Mmatre(2,2)=real(Mmattype(2,2))
!Mmatim(2,2)=aimag(Mmattype(2,2))


!EVflag=.true.

!CALL ch(Mdim,Mmatre,Mmatim,MEValre,EVflag,MEVecre,MEVecim,EVerr)


!Do n1=1,Mdim
!Do n2=1,Mdim

!Umat(n1,n2)=cmplx(MEVecre(n1,n2),MEVecim(n1,n2))

!end do  !n2
!end do  !n1


Umat(1,1)=(1.0d+00/sqrt(2.0d+00))*cmplx(-1.0d+00,0.0d+00)
Umat(1,2)=(1.0d+00/sqrt(2.0d+00))*cmplx(1.0d+00,0.0d+00)
Umat(2,1)=Umat(1,2)
Umat(2,2)=Umat(1,2)

!Umat(1,1)=cmplx(1.0d+00,0.0d+00)
!Umat(1,2)=cmplx(0.0d+00,0.0d+00)
!Umat(2,1)=Umat(1,2)
!Umat(2,2)=Umat(1,1)

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
END SUBROUTINE Unitary
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
SUBROUTINE BZ_edge(angle,KFSmax)
!----------------------------------------------------------------------
use commons
implicit none
common eF
!----------------------------------------------------------------------
DOUBLE PRECISION, INTENT(IN)::angle
DOUBLE PRECISION, INTENT(OUT)::KFSmax
DOUBLE PRECISION::eF,BZconstsh!,edgeboarder
!----------------------------------------------------------------------
!----------------------------------------------------------------------

!edgeboarder=10.0d+00**(-10.0d+00)

BZconstsh=BZconst-edgeboarder
if(angle>2.0d+00*pi .or. angle <0.0d+00)then

print*,'from-BZ-edge: Ups!angle out of [0,2*pi]'

endif


if(angle < (pi/3.0d+00))then

KFSmax=(sqrt(3.0d+00))*BZconstsh/(sin(angle)+(sqrt(3.0d+00))*cos(angle))

elseif(angle .GE.(pi/3.0d+00) .and. angle < (2.0d+00*pi/3.0d+00))then

KFSmax=((sqrt(3.0d+00))/2.0d+00)*BZconstsh/sin(angle)

elseif(angle .GE.(2.0d+00*pi/3.0d+00) .and. angle < (pi))then

KFSmax=(sqrt(3.0d+00))*BZconstsh/(sin(angle)-(sqrt(3.0d+00))*cos(angle))

elseif(angle .GE.(pi) .and. angle < (4.0d+00*pi/3.0d+00))then

KFSmax=-(sqrt(3.0d+00))*BZconstsh/(sin(angle)+(sqrt(3.0d+00))*cos(angle))

elseif(angle .GE.(4.0d+00*pi/3.0d+00) .and. angle < (5.0d+00*pi/3.0d+00))then

KFSmax=-((sqrt(3.0d+00))/2.0d+00)*BZconstsh/sin(angle)

elseif(angle .GE.(5.0d+00*pi/3.0d+00) .and. angle .LE. (2.0d+00*pi))then

KFSmax=-(sqrt(3.0d+00))*BZconstsh/(sin(angle)-(sqrt(3.0d+00))*cos(angle))

endif


END SUBROUTINE BZ_edge
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
SUBROUTINE New_Folding(kx,ky,kxe,kye,id)
!----------------------------------------------------------------------
use commons
implicit none
common eF
!----------------------------------------------------------------------
DOUBLE PRECISION, INTENT(IN)::kx,ky
DOUBLE PRECISION, INTENT(OUT)::kxe,kye
INTEGER, INTENT(OUT)::id
DOUBLE PRECISION::z,ax,ay,bx,by,eF,x1,x2,angle,a&
                 ,kminx,kminy,KFSmax,ktryx,ktryy,ktryabs,kminabs,keabs,maxabsk,phired
INTEGER::ix1,ix2,BZid,s1,s2,icount
DOUBLE PRECISION,DIMENSION(1:3)::signval
!----------------------------------------------------------------------
!----------------------------------------------------------------------
signval(1)=-1.0d+00
signval(2)=0.0d+00
signval(3)=1.0d+00

a=BZconst*2.0d+00

icount=0

z=(sqrt(3.0d+00))*a/2.0d+00

ax=a*3.0d+00/4.0d+00
ay=a*(sqrt(3.0d+00))/4.0d+00

bx=0.0d+00
by=a*(sqrt(3.0d+00))/2.0d+00

x1=kx/ax
x2=(ky-ay*(kx/ax))/by

ix1=anint(x1)
ix2=anint(x2)

kxe=kx-(ix1*ax+ix2*bx)
kye=ky-(ix1*ay+ix2*by)

!print*,'Inside new Folding:',x1,x2,ix1,ix2,kxe/BZconst,kye/BZconst
!print*,'ax,by',ax/BZconst,by/BZconst

angle=atan2(kye,kxe)
angle=(-SIGN(1.0d+00,angle)+1.0d+00)*(pi)+angle

phired=angle-(pi/3.0d+00)*int(angle/(pi/3.0d+00))
maxabsk=(a/2.0d+00)*(tan(pi/3.0d+00))/(sin(phired)+(tan(pi/3.0d+00))*(cos(phired)))

kminx=kxe
kminy=kye
kminabs=sqrt((kminx**2)+(kminy**2))

id=1

Do

CALL Hexagonal(kxe,kye,BZid)
if(BZid==2 .or. BZid==1)then
exit
endif


if(maxabsk .GE. kminabs)then    !BZid .ne. 0)then
exit
endif
	Do s1=1,3
	Do s2=1,3

	ktryx=kxe+(signval(s1))*ax+(signval(s2))*bx
	ktryy=kye+(signval(s1))*ay+(signval(s2))*by
	ktryabs=sqrt((ktryx**2)+(ktryy**2))

		if(ktryabs<kminabs)then

		kminx=ktryx
		kminy=ktryy
		kminabs=ktryabs

		kxe=kminx
		kye=kminy
		angle=atan2(kye,kxe)
		angle=(-SIGN(1.0d+00,angle)+1.0d+00)*(pi)+angle

                phired=angle-(pi/3.0d+00)*int(angle/(pi/3.0d+00))
                maxabsk=(a/2.0d+00)*(tan(pi/3.0d+00))/(sin(phired)+(tan(pi/3.0d+00))*(cos(phired)))

		keabs=sqrt((kxe**2)+(kye**2))

                !CALL BZ_edge(angle,KFSmax)

			if(abs(keabs-maxabsk) <10.0d+00**(-10.0d+00))then
			kxe=(keabs-10.0d+00**(-12.0d+00))*cos(angle)
			kye=(keabs-10.0d+00**(-12.0d+00))*sin(angle)
			keabs=sqrt((kxe**2)+(kye**2))

			kminabs=keabs

			end if

		endif

icount=icount+1


	end do !s2
	end do !s1

if(icount>100)then
print*,'New-folding: Ups: too far from oorigin'
id=0
exit
endif

end do

END SUBROUTINE New_Folding
!----------------------------------------------------------------------
