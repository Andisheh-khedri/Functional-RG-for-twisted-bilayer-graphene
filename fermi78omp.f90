!-------------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------------
module commons
!-------------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------------

DOUBLE PRECISION,PARAMETER::t=1.0d+00                      !The nearest-neighbor hopping amplitude
DOUBLE PRECISION,PARAMETER::u=4.0d+00*t                    !The strength of local Coulomb interaction

DOUBLE PRECISION,PARAMETER::chemical_potential=-1.0d+00*t  !The chemical potential

INTEGER, PARAMETER::npatch=6*2                             !number of patches (multiples of 6 for Hexagonal lattice)


INTEGER, PARAMETER::Ndis=500                               !Number of kx,ky points for the linear interpolation of dispersion! (Ndis*Ndis regular grid)


DOUBLE PRECISION,PARAMETER::beta=1.0d+00&
                            /((10.0d+00**(-3.0d+00))*t)    !Inverse temperature (i.e., broadening of delta functions appearing due to Fermi-Dirac statistics)
                                                           !for analytic fu model (beta=10.0d+00**(-6.0d+00))  

DOUBLE PRECISION,PARAMETER::error_roots_disp=&
                            10.0d+00**(-4.0d+00)           !Error for finding the Fermi Surfaces (should be smaller than 1/beta)
                                                           !Always order(s) of magnitute less than beta

INTEGER, PARAMETER::num_ODE_steps=400                      !Numbr of ODE-Solver steps depends on 'ODE_steps_size'

DOUBLE PRECISION, PARAMETER::ODE_steps_size=1.05d+00       !ODE-Solver step-size 

DOUBLE PRECISION, PARAMETER::max_vertex_required=14.0d+00  !maximum vertex after which the program stops!


INTEGER, PARAMETER::nspin=2                                !number of spin degrees of freedom
INTEGER, PARAMETER::norbital=2                             !number orbital degrees of freedom
INTEGER, PARAMETER::ntype=2                                !number of atoms per unit cell (A/B) / number of bands per spin per orbit
INTEGER, PARAMETER::nband=nspin*norbital*ntype             !number of bands


DOUBLE PRECISION,PARAMETER::Lconst=1.0d+00                 !Lattice constant
DOUBLE PRECISION,PARAMETER::pi=atan2(0.0d+00,-t)           !pi
DOUBLE PRECISION,PARAMETER::BZconst=4.0d+00*pi&
                           /(3.0d+00*Lconst*sqrt(3.0d+00)) !Reciprocical lattice constrant


DOUBLE PRECISION, PARAMETER::deltak_gv=1.0d+00*pi &
                             /(10.0d+00**(4.0d+00))        !the inftesimal kx/ky to obtain group velocity


INTEGER, PARAMETER::IdFSrotsym=1                           !If '=1', the rotational symmetry will be employed!
INTEGER, PARAMETER::Idantisym=1                            !If '=1', the anti-symmetry properties of vertex will be employed! 
INTEGER, PARAMETER::Idbandsym=1                            !If '=1', the SU(2)-symmetry of bands will be employed!

INTEGER, PARAMETER::Ntot=(((nband)**4))*(npatch**3)        !total dimension of vertex functions

INTEGER, PARAMETER::nprotpatch=6                           !integer to specify rotational symmetry 2*pi/nprotpatch

INTEGER, PARAMETER::nprot=npatch/nprotpatch                !number of patches per sector 2*pi/nprotpatch


INTEGER, PARAMETER::Ntotal=2*((npatch**3)*(nband**4))      !number of differential equations to be solved


CHARACTER*4,PARAMETER::vertex_flag='real'                  !'real' if the initial vertex is real and 'comp' is it is complex
 
CHARACTER*4,PARAMETER::run_flag='real'                     !'real' reads the dft data and interpolates on a reular grid
                                                           !'test' reads from file the interpolated data

!-------------------------------------------------------------------------------------------------------------------------------------------
CHARACTER*3,PARAMETER::patchflag='sim'                     !options: 'sim' for simple patching and 'com' for complicated patching

CHARACTER*3,PARAMETER::disp_flag='ana'                     ! options: 'ana' for the analytical dispersiob and 'num' for numerical dispersion 
                                                           ! and 'try' to test multi Fermi Surface scheme

INTEGER::icall                                             ! number of times the differential equation solver is being called

!---------------------------------------
! Technical parameters-DFT-related:
!---------------------------------------
INTEGER*4,PARAMETER::Nnumdisp=127141                       !Number of kpoints in DFT data


INTEGER, PARAMETER::num_theta_ired=100                     !Number of angles in the irreducible wedge (0,pi/3) to extract information about the bands                   
INTEGER, PARAMETER::num_kvec_thf=50                        !Number of k-vectors to find the extermums of dispersion for a given angle                   


!---------------------------------------
! Technical parameters:
!---------------------------------------

COMPLEX*16, PARAMETER::ic=cmplx(0.0d+00,1.0d+00)

DOUBLE PRECISION,PARAMETER::patcherr1=t*(10.0d+00**(-15.0d+00))
DOUBLE PRECISION,PARAMETER::patcherr2=t*(10.0d+00**(-14.0d+00))
DOUBLE PRECISION,PARAMETER::zeroHex=t*(10.0d+00**(-12.0d+00))
DOUBLE PRECISION,PARAMETER::zeroHexedge=t*(10.0d+00**(-12.0d+00))
DOUBLE PRECISION,PARAMETER::edgeboarder=10.0d+00**(-10.0d+00)         

INTEGER, PARAMETER::idpermute=0,Idantisymextra=0,ISOC=0 
!-----------------------------------------------------------------------------------------------------
!------------------------------------------------------------------
end module
!------------------------------------------------------------------
!------------------------------------------------------------------
program HFRG

use commons

implicit none

common eF,kxsamp,kysamp,dissamp,maxvertex,p4int,thetap,KFSpmat&
       ,imax,divflag,bands_info_mat,Kext_patch_mat,n_angles_mat_a,angles_mat_a,nFS_mat_a

!------------------------------------------------------------------
!  Parameter declaration:
!------------------------------------------------------------------
Integer,parameter::nangle=3001,nrot=4
Double precision::density,eF,angle,zero,scalecommon,k1,kzerop,kmax&
                 ,kFS,resol,maxvertex,groupvelx,groupvely,base&
                 ,kpx,kpy,angle1,angle2,angle3,KFS1,KFS2,KFS3,anglek,anglered
Integer::eFint,angleint,k1int,nroots,NEQ,n1,n2,n3,n4,p1,p2,p3,i&
         ,icount,irout,k,p4v,kxint,kyint,n1spin,n2spin,n3spin,n4spin&
         ,n1orb,n2orb,n3orb,n4orb,n1orbminus,l1,l2,nroots1,nroots2,kint&
         ,n1type,n2type,n3type,n4type,n1p,n2p,n3p,n4p,imax,p3max,kk,divflag,xcount

DOUBLE PRECISION,dimension(1:nband,1:npatch)::thetap,KFSpmat,thetapn,KFSpmatn
INTEGER,DIMENSION(1:nband,1:npatch)::nrootspmat,nrootspmatn

!------------------------------------------------------------------

INTEGER:: openstatus,Flag, ISTATE, IOPT, MF,MU, ML, LRW, LIW
DOUBLE PRECISION::y,yout, RTOL, ATOL,JAC

DOUBLE PRECISION,allocatable,dimension(:)::RWORK
INTEGER,allocatable,dimension(:)::IWORK
!------------------------------------------------------------------
! External functions
!------------------------------------------------------------------
DOUBLE PRECISION:: fermi,dispersion,groupvel,numgroupvel
! their corresponding variables:
DOUBLE PRECISION::mom,theta

INTEGER,DIMENSION(1:nband,1:nband,1:nband,1:nband,1:npatch,1:npatch,1:npatch)::imat
INTEGER,allocatable,DIMENSION(:):: n1mat,n2mat,n3mat,n4mat,p1mat,p2mat,p3mat

INTEGER,DIMENSION(1:nband,1:nband,1:nband,1:nband,1:npatch,1:npatch,1:npatch)::p4int
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
DOUBLE PRECISION,allocatable,DIMENSION(:)::Yvertexre,YDvertexre,Yvertexrebefore&
                                          ,Yvertex,YDvertex,Yvertexbefore
INTEGER::nchcommon,icommon
!------------------------------------------------------------------
INTEGER,DIMENSION(1:npatch,1:npatch,1:npatch)::nrep,p4extra,nmean
INTEGER,DIMENSION(1:npatch,1:npatch,1:npatch,1:npatch)::p3mean
INTEGER,DIMENSION(1:npatch,1:npatch)::nmis
INTEGER,DIMENSION(1:nband)::nspinmat,norbmat,ntypemat
!------------------------------------------------------------------
! spin to band transformation:
!------------------------------------------------------------------
DOUBLE PRECISION,DIMENSION(1:nband)::Ematp4,Emateisp4
INTEGER::s1lab,s2lab,s3lab,s4lab,iprime,s1plab,s2plab,s3plab,s4plab,idprime
COMPLEX*16::VSO12,VSO21
COMPLEX*16::summat

DOUBLE PRECISION,allocatable,DIMENSION(:,:,:,:,:,:,:)::p4vec

COMPLEX,allocatable,DIMENSION(:)::Yvertexcomplex
DOUBLE PRECISION,allocatable,DIMENSION(:)::YvertexABinitial

INTEGER,DIMENSION(1:npatch,1:nband)::i1dmat
INTEGER,DIMENSION(1:npatch,1:npatch,1:nband,1:nband)::i2dmat
INTEGER,DIMENSION(1:npatch,1:npatch,1:npatch,1:nband,1:nband,1:nband)::i3dmat
INTEGER,DIMENSION(1:nband,1:nband,1:nband,1:nband,1:npatch,1:npatch,1:npatch)::i4dmat
!------------------------------------------------------------------
!------------------------------------------------------------------
! checking!
!------------------------------------------------------------------
DOUBLE PRECISION::Hexdisp,kpxr,kpyr,kpxri,kpyri,KFSr,angler&
                  ,dispbefore,KFSmax,kpxrnew,kpyrnew
INTEGER::BZcount,BZid,nkxint,nkyint,neiint,loopcount,nneimax&
        ,neicount,numsign,nrootsmat
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
INTEGER::n1prime,n2prime,n3prime,n4prime
COMPLEX,DIMENSION(1:ntype,1:ntype,1:ntype,1:ntype)::sumcheck,sumcheckback
INTEGER,DIMENSION(1:nband,1:nband,1:nband,1:nband,1:npatch,1:npatch,1:npatch)::iredmat
!------------------------------------------------------------------
DOUBLE PRECISION::thetamin,thetamax,thetabefore,KFSthetamax,energy,kxmax,kymax,dispersionana,discheck
DOUBLE PRECISION,DIMENSION(1:npatch)::anglematmin,anglematmax,anglemat
INTEGER::idfold
character*6::nbands
character*4::phis


DOUBLE PRECISION,DIMENSION(1:Ndis)::kxsamp,kysamp
DOUBLE PRECISION,DIMENSION(1:nband,1:Ndis,1:Ndis)::dissamp

INTEGER::nsign,flag_ex,unitm
DOUBLE PRECISION::ex_min,ex_max
DOUBLE PRECISION,DIMENSION(1:10)::KFSext_ex

DOUBLE PRECISION,DIMENSION(1:nband,1:npatch)::kmin_mat,kmax_mat

DOUBLE PRECISION,DIMENSION(1:Nnumdisp)::kxdftmat,kydftmat,dispdftmat

INTEGER*4,parameter::N_interp=Ndis*Ndis
REAL*8,DIMENSION(1:N_interp)::kpxmat_interp,kpymat_interp,dispmat_interp

DOUBLE PRECISION::p_interp,dispcheck,kpxmin,dispersion_theta,e_min,e_max

INTEGER::i_min,i_max,icounter,n_angles,nFS,j,ns
DOUBLE PRECISION::angle_min,angle_max,energy1,energy2,kred
DOUBLE PRECISION,DIMENSION(1:num_theta_ired)::tmat
!------------------------------------------------------------------
!------------------------------------------------------------------
real :: start_time, stop_time,stack_size_check
integer :: count_i,count_f,count_start_flow,count_end_flow
real::count_rate
!-----------------------------------------------------------------------
DOUBLE PRECISION,DIMENSION(1:nband,1:2,1:npatch,1:2)::kmin_mat_a,kmax_mat_a
INTEGER,DIMENSION(1:nband,1:2,1:npatch)::nFS_mat_a
DOUBLE PRECISION,DIMENSION(1:nband,1:npatch)::kred_mat_a


DOUBLE PRECISION,DIMENSION(1:nband)::angle_min_mat,angle_max_mat,e_min_mat,e_max_mat
INTEGER,DIMENSION(1:nband)::nFS_max_mat

INTEGER::i_phase
CHARACTER*3::patch_flag

DOUBLE PRECISION,DIMENSION(1:5,1:nband)::bands_info_mat
DOUBLE PRECISION,DIMENSION(1:2,1:nband,1:2,1:npatch,1:2)::Kext_patch_mat
INTEGER,DIMENSION(1:nband,1:2)::n_angles_mat_a
DOUBLE PRECISION,DIMENSION(1:nband,1:2,1:npatch)::angles_mat_a

INTEGER::n1max,n2max,n3max,n4max
!-----------------------------------------------------------------------
call cpu_time(start_time)
CALL SYSTEM_CLOCK(count_i, count_rate)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------
open(unit=2,file="maxvertex_flow_num_ana.txt",status="new",action="write",position="rewind",IOSTAT=openstatus)
!-----------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------
open(unit=31,file="band1_real.txt",status="new",action="write",position="rewind",IOSTAT=openstatus)
open(unit=32,file="band2_real.txt",status="new",action="write",position="rewind",IOSTAT=openstatus)
open(unit=33,file="band3_real.txt",status="new",action="write",position="rewind",IOSTAT=openstatus)
open(unit=34,file="band4_real.txt",status="new",action="write",position="rewind",IOSTAT=openstatus)
!-----------------------------------------------------------------------------------------------------
open(unit=61,file="band1_interp.txt",status="old",action="read",position="rewind",IOSTAT=openstatus)
open(unit=62,file="band2_interp.txt",status="old",action="read",position="rewind",IOSTAT=openstatus)
open(unit=63,file="band3_interp.txt",status="old",action="read",position="rewind",IOSTAT=openstatus)
open(unit=64,file="band4_interp.txt",status="old",action="read",position="rewind",IOSTAT=openstatus)
!-----------------------------------------------------------------------------------------------------
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
! Allocation of memory:
!----------------------------------------------------------------------
allocate(n1mat(1:Ntot),n2mat(1:Ntot),n3mat(1:Ntot),n4mat(1:Ntot)&
        ,p1mat(1:Ntot),p2mat(1:Ntot),p3mat(1:Ntot))
allocate(p4vec(1:nband,1:nband,1:nband,1:npatch,1:npatch,1:npatch,1:2))

allocate(Yvertexcomplex(1:Ntot),YvertexABinitial(1:Ntot))
!----------------------------------------------------------------------
if(vertex_flag=='real')then
NEQ=Ntot
elseif(vertex_flag=='comp')then
NEQ=Ntotal
endif
!----------------------------------------------------------------------
allocate(Yvertex(1:NEQ),YDvertex(1:NEQ),Yvertexbefore(1:NEQ))
!----------------------------------------------------------------------
Yvertex(1:NEQ)=0.0d+00
!----------------------------------------------------------------------
! Initialization of variables:
!---------------------------------------------------------------------
kxsamp(1:Ndis)=0.0d+00
kysamp(1:Ndis)=0.0d+00
dissamp(1:nband,1:Ndis,1:Ndis)=0.0d+00

kmin_mat(1:nband,1:npatch)=0.0d+00
kmax_mat(1:nband,1:npatch)=BZconst
nrootspmat(1:nband,1:npatch)=0
!----------------------------------------------------------------------
Yverdummyu=cmplx(0.0d+00,0.0d+00)
Yverdummyv=cmplx(0.0d+00,0.0d+00)
Yverdummy=cmplx(0.0d+00,0.0d+00)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
Yvertexcomplex(1:Ntot)=cmplx(0.0d+00,0.0d+00)
Yvertex(1:NEQ)=0.0d+00
!----------------------------------------------------------------------
! Initialization for the solver:
!----------------------------------------------------------------------
MF=10
ML=1
MU=1
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Dimension Declaration of WORK:
!----------------------------------------------------------------------
If(MF==10)then
LRW=20+20*NEQ
print*,'LRW is=',LRW
LIW=20
elseif(MF==21 .or. MF==22)then
LRW=22+9*NEQ+NEQ**2 
LIW=20+NEQ
elseif(MF==24 .or. MF==25)then
LRW=22+10*NEQ+(2*ML+MU)*NEQ
LIW=20+NEQ
endif

!------------------------------------------------------------------
allocate(RWORK(LRW),IWORK(LIW))
!------------------------------------------------------------------
!------------------------------------------------------------------
Flag=1
!------------------------------------------------------------------
IOPT=1            
RWORK(5:10)=0.0d+00 
IWORK(5:10)=0
IWORK(6)=20000 
resol=ODE_steps_size
!------------------------------------------------------------------
!------------------------------------------------------------------
divflag=0
CALL SYSTEM_CLOCK(count_start_flow, count_rate)
!----------------------------------------------------------------------
eF=chemical_potential
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Step 0: Generating Data (grid+dispersion)
!----------------------------------------------------------------------
!CALL GENERATE_DATA_GIVEN_GRID(127141)    !!CALL GENERATE_DATA(200)
!-------------------------------------------------------------
!-------------------------------------------------------------
!CALL Sample_points(Ndis,kxsamp,kysamp,kpxmat_interp,kpymat_interp)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!Do l1=1,nband/2
!Do j=1,Ndis*Ndis


!mom=sqrt(((kpxmat_interp(j))**2)+((kpymat_interp(j))**2))

!angle=atan2(kpymat_interp(j),kpxmat_interp(j))
!angle=(-SIGN(1.0d+00,angle)+1.0d+00)*(pi)+angle

!i=60+l1

!read(i,*)dispcheck

!i=30+l1

!write(i,*)kpxmat_interp(j)/BZconst,kpymat_interp(j)/BZconst,dispcheck,dispersionana(l1*2,mom,angle,0.0d+00)

!end do
!end do
!----------------------------------------------------------------------
!----------------------------------------------------------------------
if(pi<0.0d+00)then
!----------------------------------------------------------------------
if(disp_flag=='num')then
!----------------------------------------------------------------------
! Step 1:
!----------------------------------------------------------------------
! Reading DFT data and interpolating for reqular grid:
!----------------------------------------------------------------------
! (a) defining a regular grid:
!-------------------------------------------------------------
CALL Sample_points(Ndis,kxsamp,kysamp,kpxmat_interp,kpymat_interp)
!-------------------------------------------------------------
! (b) intepolation for the reqular grid:
!-------------------------------------------------------------
CALL Transform_regular_grid(run_flag,Nnumdisp,Ndis,dissamp)
!----------------------------------------------------------------------
! (c) comparing interpolated dipersion with DFT data
!----------------------------------------------------------------------
!CALL Comp_interp_DFT(Nnumdisp)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Step 2: Obtaining information about the bands (directly from DFT data)
!----------------------------------------------------------------------
CALL Bands_Info(run_flag,num_theta_ired,angle_min_mat,angle_max_mat,e_min_mat,e_max_mat,nFS_max_mat)

bands_info_mat(1,1:nband)=1.0d+00*nFS_max_mat(1:nband)
bands_info_mat(2,1:nband)=angle_min_mat(1:nband)
bands_info_mat(3,1:nband)=angle_max_mat(1:nband)
bands_info_mat(4,1:nband)=e_min_mat(1:nband)
bands_info_mat(5,1:nband)=e_max_mat(1:nband)

!bands_info_mat(1,1:nband)=1.0d+00

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Step 3: Obtaining information about the Fermi Surfaces for a given
!         chemical potential. This can help to set up the 
!         patching scheme!
!----------------------------------------------------------------------

CALL patch_assist(num_theta_ired,num_kvec_thf,n_angles_mat_a,angles_mat_a,nFS_mat_a,kmin_mat_a,kmax_mat_a,kred_mat_a)

Kext_patch_mat(1,1:nband,1:2,1:npatch,1:2)=kmin_mat_a(1:nband,1:2,1:npatch,1:2)
Kext_patch_mat(2,1:nband,1:2,1:npatch,1:2)=kmax_mat_a(1:nband,1:2,1:npatch,1:2)

!----------------------------------------------------------------------
else

bands_info_mat(1,1:nband)=1.0d+00

endif
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Step 4: Introducing patch-scheme:
!----------------------------------------------------------------------
!----------------------------------------------------------------------
CALL patch_generator(eF,thetap,KFSpmat,nrootspmat)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Step 5: Finding patching index associated to
!         vec(p4)=vec(p1)+vec(p2)-vec(p3), assuming that vec(p_1,2,3)
!         lies on the Fermi Surface:
!----------------------------------------------------------------------
CALL p4_generator(Npatch,thetap,KFSpmat,eF,bands_info_mat,p4vec,p4int)
!----------------------------------------------------------------------
print*,'p4_generator done!'
!-----------------------------------------------------------------------
!----------------------------------------------------------------------
! Step 6: Introducing the initial vertex tensor:
!----------------------------------------------------------------------
! (a) Transforming (n1,n2,n3,n4,p1,2,p3) array to one dimensional
!----------------------------------------------------------------------
CALL imat_generator(nband,npatch,imat,n1mat,n2mat,n3mat,n4mat)
!----------------------------------------------------------------------
print*,'imat_generator done!'
!----------------------------------------------------------------------
! (b) the initial vertex:
!-----------------------------------------------------------------------
CALL Initial_vertex('local',thetap,KFSpmat,Yvertexcomplex)
!----------------------------------------------------------------------
print*,'initial vertex done!'
!----------------------------------------------------------------------
!--------------------------------------------------------
Yvertex(1:Ntot)=real(Yvertexcomplex(1:Ntot))
!--------------------------------------------------------
if(vertex_flag=='comp')then
Yvertex((Ntot+1):(2*Ntot))=aimag(Yvertexcomplex(1:Ntot))
endif
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Step 7: Solving the Flow equations:
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
maxvertex=Yvertex(imat(1,2,1,2,1,1,1))
!------------------------------------------------------------------
Do i=1,num_ODE_steps

irout=i

if(abs(maxvertex)<max_vertex_required)then
!------------------------------------------------------------------
RTOL=t*10.0d+00**(-3.0d+00)   
ATOL=t*10.0d+00**(-3.0d+00)
ISTATE=1
x=10.00d+00*t/((resol)**(i-1))
xout=10.00d+00*t/((resol)**(i))
!------------------------------------------------------------------
icall=0
!------------------------------------------------------------------

CALL Fflow (NEQ, x , Yvertex, YDvertex)

CALL DLSODE (Fflow, NEQ, Yvertex, x, xout, ITOL, RTOL, ATOL, Flag,ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF)

!------------------------------------------------------------------
if(divflag==0)then
Yvertexbefore(1:NEQ)=Yvertex(1:NEQ)
elseif(divflag==1)then
Yvertex(1:NEQ)=Yvertexbefore(1:NEQ)
end if
!------------------------------------------------------------------

print*,"Fcur",i,ISTATE,10.00d+00*t/((resol)**(i-1)),xout/t,maxvertex

print*,'max_vertex belongs to:','i=',imax,'n1=',n1mat(imax),'n2=',n2mat(imax),'n3=',n3mat(imax),'n4=',n4mat(imax)

write(2,*)xout/t,maxvertex

!------------------------------------------------------------------
endif
!------------------------------------------------------------------
end do !i
!------------------------------------------------------------------
!------------------------------------------------------------------
CALL SYSTEM_CLOCK(count_end_flow, count_rate)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Step 8: printing out final vertex in the txt-file
!-----------------------------------------------------------------------
CALL write_vertex(NEQ,Yvertex,n1mat(imax),n2mat(imax),n3mat(imax),n4mat(imax))
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
endif !pi<0
!-----------------------------------------------------------------------
!----------------------------------------------------------------------
! Deallocation of memory:
!----------------------------------------------------------------------
deallocate(RWORK,IWORK)
deallocate(n1mat,n2mat,n3mat,n4mat,p1mat,p2mat,p3mat,p4vec)
deallocate(Yvertexcomplex,YvertexABinitial)
deallocate(Yvertex,YDvertex,Yvertexbefore)
!-----------------------------------------------------------------------
call cpu_time(stop_time)
CALL SYSTEM_CLOCK(count_f, count_rate)
print*,"total time",stop_time - start_time,real(count_f-count_i)/count_rate,real(count_end_flow-count_start_flow)/count_rate
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
end program HFRG
!------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------
!******* Flow equations***********************
!------------------------------------------------------------------------------------------------------------
subroutine Fflow(NEQ, x, Yvertex, YDvertex)

USE OMP_LIB

use commons

implicit none

common eF,kxsamp,kysamp,dissamp,maxvertex,p4int,thetap,KFSpmat&
       ,imax,divflag,bands_info_mat,Kext_patch_mat,n_angles_mat_a,angles_mat_a,nFS_mat_a

!-------------------------------------------------------------------------------------------------------
! Parameter declaration:
!-------------------------------------------------------------------------------------------------------
INTEGER::k,i,j1,j2,NEQ,n1,n2,n3,n4,p1,p2,p3,p4,kp,l1,l2,kp2,nrootskp2,j3,j4,nch,j,irot,irout&
         ,nroots1,nroots2,nroots3,nrootsk,nrootskp,sint,icommon,nchcommon,imin,imax,kr,kv,kstart,kfinal&
         ,p2max,p3max,n2max,n4max,iasym,iasym2,kint,kk,divflag,xcount
DOUBLE PRECISION::eF,y,scalecommon,groupvel,Matsumpart1,Matsumpart2,step,step1,step2&
                  ,angle1,angle2,angle3,KFS1,KFS2,KFS3,KFSk,KFSkp,stepsize&
                  ,anglek,anglekp,energy,Matsum,anglekp2,KFSkp2,radkp,radkp2&
                  ,anglekpred,kpx,kpy,maxvertex,kp2x,kp2y,anglekp2red&
                  ,x,minvertex,r3,dis,maxvertex_n
DOUBLE PRECISION,DIMENSION(1:2)::sgn
!-------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------
! External functions:
!-------------------------------------------------------------------------------------------------------
DOUBLE PRECISION::dispersion,fermi
INTEGER,DIMENSION(1:Ntot):: n1mat,n2mat,n3mat,n4mat,p1mat,p2mat,p3mat

!INTEGER,DIMENSION(1:Ntot):: nrep
INTEGER,ALLOCATABLE,DIMENSION(:)::nrep

!INTEger,DIMENSION(1:nband,1:nband,1:nband,1:nband,1:npatch,1:npatch,1:npatch)::imat
integer::imat
INTEger,DIMENSION(1:nband,1:nband,1:nband,1:nband,1:npatch,1:npatch,1:npatch)::p4int
!-------------------------------------------------------------------------------------------------------
DOUBLE PRECISION::Yang,YDang
COMPLEX*16,DIMENSION(1:NEQ)::Taupp,Tauphc,Tauphd
!-------------------------------------------------------------------------------------------------------
DOUBLE PRECISION,DIMENSION(1:2,1:nband,1:npatch)::KFSmat,gvelFSmat
INTEGER,DIMENSION(1:2,1:nband,1:npatch)::nrootsmat
DOUBLE PRECISION::Matsumpp,Matsumph

DOUBLE PRECISION, ALLOCATABLE,DIMENSION(:,:,:)::Matsumpp1,Matsumpp2,Matsumph1,Matsumph2&
                                                        ,Matsumph3,Matsumph4

INTEGER,DIMENSION(1:2,1:npatch,1:npatch,1:npatch)::kpmat

DOUBLE PRECISION,DIMENSION(1:nband,1:npatch)::thetap
DOUBLE PRECISION,DIMENSION(1:nband,1:npatch)::KFSpmat
INTEGER,DIMENSION(1:nband,1:npatch)::nrootspmat


DOUBLE PRECISION,DIMENSION(1:2,1:nband,1:npatch)::thetakmat

INTEGER,DIMENSION(1:2,1:2,1:npatch,1:npatch,1:npatch)::kpFSmat

INTEGER::BZid,n1a,n2a,l1a,l2a,n1b,n2b,l1b,l2b,n1i,n2i,l1i,l2i

INTEGER,DIMENSION(1:nband)::ns

DOUBLE PRECISION,DIMENSION(1:NEQ)::Yvertex,YDvertex

DOUBLE PRECISION::KFSthetamax,dispersioncall,dispersioncalll1
INTEGER,DIMENSION(1:2,1:nband,1:npatch)::kintmat
INTEGER,DIMENSION(1:nband)::nref,nmin,nmax
INTEGER::nref1,nref2,n1m,n1p,n2m,n2p,l1m,l1p,l2m,l2p&
         ,n1r,n2r,l1r,l2r,n3r,nrootstotal,p1p,p2p,p3p,n3p,n4p,nrefgap

INTEGER,DIMENSION(1:npatch)::xcount_mat

DOUBLE PRECISION,DIMENSION(1:nband,1:nband,1:npatch,1:nband,1:nband,1:npatch,1:npatch)::Yphcj1,Yphcj3,Yppj1 &
                                                                                      ,Yphcj1_im,Yphcj3_im,Yppj1_im

DOUBLE PRECISION,DIMENSION(1:nband,1:nband,1:npatch,1:npatch,1:nband,1:nband,1:npatch)::Yphcj2,Yphcj4,Yppj2 &
                                                                                       ,Yphcj2_im,Yphcj4_im,Yppj2_im

INTEGER,DIMENSION(1:nband,1:nband,1:npatch,1:nband,1:nband,1:npatch,1:npatch)::iphmat,ippmat,p4ppmat,p4phmat

INTEGER::p3red,p1red,p2red


INTEGER,DIMENSION(1:nband)::nbar_mat

DOUBLE PRECISION,DIMENSION(1:Ndis)::kxsamp,kysamp
DOUBLE PRECISION,DIMENSION(1:nband,1:Ndis,1:Ndis)::dissamp

DOUBLE PRECISION,DIMENSION(1:nband,1:npatch)::kmin_mat,kmax_mat
!------------------------------------------------------------------
real:: start_time, stop_time
integer :: count_i_ph,count_f_ph,count_i_pp,count_f_pp
real::count_rate_ph,count_rate_pp
!-----------------------------------------------------------------------
DOUBLE PRECISION,DIMENSION(1:5,1:nband)::bands_info_mat
DOUBLE PRECISION,DIMENSION(1:2,1:nband,1:2,1:npatch,1:2)::Kext_patch_mat
INTEGER,DIMENSION(1:nband,1:2)::n_angles_mat_a
DOUBLE PRECISION,DIMENSION(1:nband,1:2,1:npatch)::angles_mat_a
INTEGER,DIMENSION(1:nband,1:2,1:npatch)::nFS_mat_a
!-----------------------------------------------------------------------
!--------------------------------
! parameter for OMP:
!--------------------------------
INTEGER:: NTHREADS, TID,num_threads
!----------------------------------------------------
!----------------------------------------------------
! Allocating memory:
!----------------------------------------------------
ALLOCATE(nrep(1:Ntot))

ALLOCATE(Matsumpp1(1:nband,1:nband,1:npatch)&
        ,Matsumpp2(1:nband,1:nband,1:npatch)&
        ,Matsumph1(1:nband,1:nband,1:npatch)&
        ,Matsumph2(1:nband,1:nband,1:npatch)&
        ,Matsumph3(1:nband,1:nband,1:npatch)&
        ,Matsumph4(1:nband,1:nband,1:npatch))
!----------------------------------------------------
!-------------------------------------------------------------------------------------------------------
scalecommon=x
icall=icall+1
stepsize=10.0d+00**(-6.0d+00)
xcount=0
!------------------------------------------------------------------
kpmat(1:2,1:npatch,1:npatch,1:npatch)=0
kpFSmat(1:2,1:2,1:npatch,1:npatch,1:npatch)=0

Matsumpp=0.0d+00
Matsumph=0.0d+00

Matsumpp1(1:nband,1:nband,1:npatch)=0.0d+00
Matsumpp2(1:nband,1:nband,1:npatch)=0.0d+00

Matsumph1(1:nband,1:nband,1:npatch)=0.0d+00
Matsumph2(1:nband,1:nband,1:npatch)=0.0d+00
Matsumph3(1:nband,1:nband,1:npatch)=0.0d+00
Matsumph4(1:nband,1:nband,1:npatch)=0.0d+00


gvelFSmat(1:2,1:nband,1:npatch)=0.0d+00
!KFSmat(1:2,1:npatch)=0.0d+00

nrootsmat(1:2,1:nband,1:npatch)=0

nrep(1:Ntot)=0

sgn(1)=-1.0d+00
sgn(2)=1.0d+00

!-------------------------------------------------
!print*,'nbar_mat def:'
Do n1=1,nband

if(n1==2*int(n1/2))then
nbar_mat(n1)=n1-1
else
nbar_mat(n1)=n1+1
endif

!print*,n1,nbar_mat(n1)

end do
!-------------------------------------------------

if(Idbandsym==1)then

nref(1)=1
nref(2)=1

nref(3)=3
nref(4)=3

nref(5)=5
nref(6)=5

nref(7)=7
nref(8)=7

nrefgap=2

elseif(Idbandsym==0)then

Do n1=1,nband

nref(n1)=n1

end do

nrefgap=1

endif
!------------------------------------------------------------------
!print*,'I am in Fflow...initializing YD...'
!------------------------------------------------------------------
YDvertex(1:NEQ)=0.0d+00
!------------------------------------------------------------------
!print*,'I am in Fflow...finding KFS...'
!------------------------------------------------------------------
nrootstotal=0
!------------------------------------------------------------------
Do l1=1,nband

Do sint=1,2

energy=eF+(sgn(sint)*x)

if(disp_flag=='num')then

!-------------------------------------------------------------------

kint=0

if(bands_info_mat(1,l1)<2 .or. energy<bands_info_mat(4,l1) .or. energy>bands_info_mat(5,l1))then

i=1

else

i=2

endif

Do k=1,n_angles_mat_a(l1,i)/6

Do j=1,nFS_mat_a(l1,i,k)

kint=kint+1

anglek=angles_mat_a(l1,i,k)

kmin_mat(l1,k)=kext_patch_mat(1,l1,i,k,j)
kmax_mat(l1,k)=kext_patch_mat(2,l1,i,k,j) 
!------------------------------------------------------------------

CALL roots_dispersion_theta2(l1,anglek,energy,kmin_mat(l1,k),kmax_mat(l1,k),nrootsk,KFSk)

!print*,'nroots-inside flow',nrootsk,anglek/(2.0d+00*pi),KFSk/BZconst

!------------------------------------------------------------------
!------------------------------------------------------------------
if(nrootsk>0)then

Do j1=0,5

thetakmat(sint,l1,kint+j1*(npatch/6))=anglek+j1*(pi/3.0d+00)
KFSmat(sint,l1,kint+j1*(npatch/6))=KFSk
nrootsmat(sint,l1,kint+j1*(npatch/6))=nrootsk
gvelFSmat(sint,l1,kint+j1*(npatch/6))=abs(groupvel(1,l1,KFSK,anglek))

!print*,'KFS inside-num:',anglek+j1*(pi/3.0d+00),KFSk,nrootsk,abs(groupvel(1,l1,KFSK,anglek))

end do !j1

!------------------------------------------------------------------
endif
!------------------------------------------------------------------

nrootstotal=nrootstotal+nrootsk

end do !j
end do !k

!-------------------------------------------------------------------
else

Do k=1,npatch

anglek=thetap(l1,k)


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

thetakmat(sint,l1,k)=anglek
KFSmat(sint,l1,k)=KFSk
nrootsmat(sint,l1,k)=nrootsk
gvelFSmat(sint,l1,k)=abs(groupvel(1,l1,KFSK,anglek))

!------------------------------------------------------------------
!print*,'KFS inside:',thetakmat(sint,l1,k),KFSmat(sint,l1,k),nrootsmat(sint,l1,k),gvelFSmat(sint,l1,k)
!------------------------------------------------------------------

kintmat(sint,l1,k)=kint

endif
!------------------------------------------------------------------

nrootstotal=nrootstotal+nrootsk


end do !k

endif

end do !sint
end do !l1
!------------------------------------------------------------------
if(nrootstotal>0)then
!------------------------------------------------------------------
!------------------------------------------------------------------
!----------------------------------------------------------------------
!------------------------------------------------------------------
Tauphc(1:NEQ)=0.0d+00
Tauphd(1:NEQ)=0.0d+00
!Tauphcmat=0.0d+00
!------------------------------------------------------------------
!$OMP PARALLEL SHARED(p4int,Yvertex) , &
!$OMP& PRIVATE(n1,n2,n3,n4,p1,p2,p3,i,p4,p3red) , &
!$OMP& SHARED(iphmat,ippmat,p4phmat,p4ppmat,Yphcj1,Yphcj2,Yphcj3,Yphcj4,Yppj1,Yppj2) ,&
!$OMP& SHARED(Yphcj1_im,Yphcj2_im,Yphcj3_im,Yphcj4_im,Yppj1_im,Yppj2_im) 

num_threads = omp_get_num_threads()

!print*,'Total number of threads are=',num_threads

Do p3red=1,npatch/num_threads

!Do p3=1,npatch

p3=p3red+(OMP_GET_THREAD_NUM())*(npatch/num_threads)

!print*, 'I am number',p3

Do p2=1,npatch
Do p1=1,npatch

Do n4=1,nband
Do n3=1,nband
Do n2=1,nband
Do n1=1,nband

i=imat(n1,n2,n3,n4,p1,p2,p3)

!i=(n1-1)*((nband**3)*(npatch**3))+(n2-1)*((nband**2)*(npatch**3))&
!     +(n3-1)*((nband)*(npatch**3))+(n4-1)*(npatch**3)&
!     +(p1-1)*(npatch**2)+(p2-1)*(npatch)+p3
!--------------------------------------------------------------------
! ph-c term:
!--------------------------------------------------------------------

Yphcj1(n4,n2,p2,n3,n1,p3,p1)=Yvertex(i)


Yphcj2(n4,n2,p2,p3,n3,n1,p1)=Yvertex(i)


Yphcj3(n4,n2,p2,n1,n3,p1,p3)=Yvertex(i)


Yphcj4(n4,n2,p2,p1,n1,n3,p3)=Yvertex(i)

!--------------------------------------------------------------------
if(vertex_flag=='comp')then

Yphcj1_im(n4,n2,p2,n3,n1,p3,p1)=Yvertex(i+(NEQ/2))


Yphcj2_im(n4,n2,p2,p3,n3,n1,p1)=Yvertex(i+(NEQ/2))


Yphcj3_im(n4,n2,p2,n1,n3,p1,p3)=Yvertex(i+(NEQ/2))


Yphcj4_im(n4,n2,p2,p1,n1,n3,p3)=Yvertex(i+(NEQ/2))

elseif(vertex_flag=='real')then

Yphcj1_im(n4,n2,p2,n3,n1,p3,p1)=0.0d+00


Yphcj2_im(n4,n2,p2,p3,n3,n1,p1)=0.0d+00


Yphcj3_im(n4,n2,p2,n1,n3,p1,p3)=0.0d+00


Yphcj4_im(n4,n2,p2,p1,n1,n3,p3)=0.0d+00

endif

!--------------------------------------------------------------------

Yppj1(n4,n3,p3,n2,n1,p2,p1)=Yvertex(i)

Yppj2(n4,n3,p3,p2,n2,n1,p1)=Yvertex(i)

!--------------------------------------------------------------------
if(vertex_flag=='comp')then

Yppj1_im(n4,n3,p3,n2,n1,p2,p1)=Yvertex(i+(NEQ/2))
Yppj2_im(n4,n3,p3,p2,n2,n1,p1)=Yvertex(i+(NEQ/2))

elseif(vertex_flag=='real')then

Yppj1_im(n4,n3,p3,n2,n1,p2,p1)=0.0d+00
Yppj2_im(n4,n3,p3,p2,n2,n1,p1)=0.0d+00

endif
!--------------------------------------------------------------------

!iphmat(n2,n4,p2,n3,n1,p3,p1)=i

iphmat(n4,n2,p2,n1,n3,p1,p3)=i

ippmat(n4,n3,p3,n2,n1,p2,p1)=i

p4=p4int(n1,n2,n3,n4,p1,p2,p3)

!p4phmat(n2,n4,p2,n3,n1,p3,p1)=p4

!p4phmat(n4,n2,p2,n1,n3,p1,p3)=p4

!p4ppmat(n4,n3,p3,n2,n1,p2,p1)=p4

end do
end do
end do
end do


end do
end do

end do !p3red

!end do !p3

!----------------------------------------------------------------------
!$OMP END PARALLEL
!----------------------------------------------------------------------
!------------------------------------------------------------------
! Matsumph
!------------------------------------------------------------------
CALL SYSTEM_CLOCK(count_i_ph, count_rate_ph)
!------------------------------------------------------------------
!$OMP PARALLEL SHARED(stepsize,x,eF,sgn,tauphc,tauphd,nref,nrefgap) , &
!$OMP& PRIVATE(n1,p1,n2,p2,l1,l2,p3max,energy,kpx,kpy,anglekp,radkp,step,sint,k,anglek,i) , &
!$OMP& SHARED(thetakmat,nrootsmat,gvelFSmat,thetap,KFSpmat,Yvertex,p4int,KFSmat) , &
!$OMP& PRIVATE(Matsumph,angle1,KFS1,angle2,KFS2,j1,j2,j3,j4,p3,n3,n4,p4,p1p,p2p,p3p,p1red) , &
!$OMP& PRIVATE(Matsumph1,Matsumph2,Matsumph3,Matsumph4,dispersioncall,dispersioncalll1,n1r,n3r,l1r,l2r) &
!$OMP& SHARED(Yphcj1,Yphcj2,Yphcj3,Yphcj4,p4phmat,iphmat,Yphcj1_im,Yphcj2_im,Yphcj3_im,Yphcj4_im)
!----------------------------------------------------------------------
!----------------------------------------------------------------------

Do p3=1,nprot!npatch

!Do p1=1,npatch

!print*, 'num_threads=',num_threads

Do p1red=1,npatch/num_threads

p1=p1red+(OMP_GET_THREAD_NUM())*(npatch/num_threads)

!print*, 'I am number',p1,'out of',num_threads*(npatch/num_threads)

Do n3=1,nband

Do n1=1,nband

angle1=thetap(n1,p1)
KFS1=KFSpmat(n1,p1)

angle2=thetap(n3,p3)
KFS2=KFSpmat(n3,p3)

!------------------------------------------------------------------

Do sint=1,2

energy=eF+(sgn(sint)*x)

Do k=1,npatch

Do l1=1,nband,nrefgap

anglek=thetakmat(sint,l1,k)

if (nrootsmat(sint,l1,k)>0)then

dispersioncalll1=dispersion(l1,KFSmat(sint,l1,k),thetakmat(sint,l1,k),eF)

Do l2=1,nband,nrefgap

!------------------------------------------------------------------ 
!------------------------------------------------------------------
! 2.ph channel : p1-p3 - Matsumph1 and Matsumph2
!------------------------------------------------------------------

kpx=KFS1*cos(angle1)-KFS2*cos(angle2)+KFSmat(sint,l1,k)*cos(anglek)
kpy=KFS1*sin(angle1)-KFS2*sin(angle2)+KFSmat(sint,l1,k)*sin(anglek)

CALL FOLDING(kpx,kpy,kpx,kpy)

anglekp=atan2(kpy,kpx)
anglekp=(-SIGN(1.0d+00,anglekp)+1.0d+00)*(pi)+anglekp

radkp=sqrt(((kpx)**2)+((kpy)**2))


dispersioncall=dispersion(l2,radkp,anglekp,eF)
!------------------------------------------------------------------
! Regulator:
!------------------------------------------------------------------
if(abs(dispersioncall)>scalecommon+stepsize)then
step=1.0d+00
elseif(abs(dispersioncall)<scalecommon-stepsize)then
step=0.0d+00
else
step=0.5d+00
endif
!------------------------------------------------------------------ 
!------------------------------------------------------------------
!------------------------------------------------------------------
!new:

if(energy-eF < 0.0d+00 .and. dispersioncall<0.0d+00)then
Matsumph=0.0d+00
elseif(energy-eF .GE. 0.0d+00 .and. dispersioncall .GE. 0.0d+00)then
Matsumph=0.0d+00
elseif(energy-eF < 0.0d+00)then
Matsumph=1.0d+00/(dispersioncalll1-dispersioncall)
elseif(energy-eF .GE. 0.0d+00)then
Matsumph=-1.0d+00/(dispersioncalll1-dispersioncall)
endif

if(abs(dispersioncalll1-dispersioncall)<10.0d+00**(-14.0d+00))then

!print*,'careful! devided by zero!'

Matsumph=(beta/4.0d+00)*(((tanh(beta*(energy-eF)/2.0d+00))**2)-1.0d+00)

!print*,n1,n2,l1,l2,p1,p2,k,Matsumph,((tanh(beta*(energy-eF)/2.0d+00))**2),(beta/4.0d+00),energy-eF

endif

!------------------------------------------------------------------

Matsumph=step*Matsumph*((KFSmat(sint,l1,k))/(gvelFSmat(sint,l1,k)))

!------------------------------------------------------------------

if(sint==1)then
Matsumph1(l2,l1,k)=Matsumph

!if(p3==1 .and. p1red==1 .and. n3==6 .and. n1==3 .and. sint==1 .and. k==6 .and. l1==5 .and. l2==3)then
!print*,'Matsumph not zero!',p3,p1red,p1,n3,n1,sint,k,l1,l2,matsumph,step,KFSmat(sint,l1,k),gvelFSmat(sint,l1,k)&
!      ,dispersioncalll1,dispersioncall,l2,radkp,anglekp,eF,KFS1,angle1,KFS2,angle2
!endif


elseif(sint==2)then
Matsumph2(l2,l1,k)=Matsumph
endif

!------------------------------------------------------------------
! 3 and 4: p3-p1 - Matsumph3 and Matsumph4
!------------------------------------------------------------------

!kpx=-KFSpmat(n1,p1)*cos(thetap(n1,p1))+KFSpmat(n3,p3)*cos(thetap(n3,p3))+KFSmat(sint,l1,k)*cos(thetakmat(sint,l1,k))
!kpy=-KFSpmat(n1,p1)*sin(thetap(n1,p1))+KFSpmat(n3,p3)*sin(thetap(n3,p3))+KFSmat(sint,l1,k)*sin(thetakmat(sint,l1,k))

kpx=-KFS1*cos(angle1)+KFS2*cos(angle2)+KFSmat(sint,l1,k)*cos(anglek)
kpy=-KFS1*sin(angle1)+KFS2*sin(angle2)+KFSmat(sint,l1,k)*sin(anglek)

anglekp=atan2(kpy,kpx)
anglekp=(-SIGN(1.0d+00,anglekp)+1.0d+00)*(pi)+anglekp

radkp=sqrt(((kpx)**2)+((kpy)**2))

dispersioncall=dispersion(l2,radkp,anglekp,eF)
!------------------------------------------------------------------
! Regulator:
!------------------------------------------------------------------
if(abs(dispersioncall)>scalecommon+stepsize)then
step=1.0d+00
elseif(abs(dispersioncall)<scalecommon-stepsize)then
step=0.0d+00
else
step=0.5d+00
endif
!------------------------------------------------------------------ 
!------------------------------------------------------------------
!------------------------------------------------------------------
!new:

if(energy-eF < 0.0d+00 .and. dispersioncall<0.0d+00)then
Matsumph=0.0d+00
elseif(energy-eF .GE. 0.0d+00 .and. dispersioncall .GE. 0.0d+00)then
Matsumph=0.0d+00
elseif(energy-eF < 0.0d+00)then
Matsumph=1.0d+00/(dispersioncalll1-dispersioncall)
elseif(energy-eF .GE. 0.0d+00)then
Matsumph=-1.0d+00/(dispersioncalll1-dispersioncall)
endif

if(abs(dispersioncalll1-dispersioncall)<10.0d+00**(-14.0d+00))then

!print*,'careful! devided by zero!'

Matsumph=(beta/4.0d+00)*(((tanh(beta*(energy-eF)/2.0d+00))**2)-1.0d+00)

!print*,n1,n2,l1,l2,p1,p2,k,Matsumph,((tanh(beta*(energy-eF)/2.0d+00))**2),(beta/4.0d+00),energy-eF

endif

!------------------------------------------------------------------

Matsumph=step*Matsumph*((KFSmat(sint,l1,k))/(gvelFSmat(sint,l1,k)))

!------------------------------------------------------------------

if(sint==1)then
Matsumph3(l2,l1,k)=Matsumph
elseif(sint==2)then
Matsumph4(l2,l1,k)=Matsumph
endif
!------------------------------------------------------------------
!------------------------------------------------------------------ 
!------------------------------------------------------------------ 
!------------------------------------------------------------------
!------------------------------------------------------------------ 
end do !l2
!------------------------------------------------------------------

!------------------------------------------------------------------
end if !k on FS
end do !k
end do !l1

end do !sint
!------------------------------------------------------------------
Do p2=1,npatch

Do n2=1,nband
Do n4=1,nband,2

p4=p4int(n1,n2,n3,n4,p1,p2,p3)

!i=imat(n1,n2,n3,n4,p1,p2,p3)

!p4=p4phmat(n4,n2,p2,n1,n3,p1,p3)

i=iphmat(n4,n2,p2,n1,n3,p1,p3)

if(p4==0)then

Tauphc(i)=0.0d+00

elseif(p4>0)then

!------------------------------------------------------------------
Do sint=1,2

Do k=1,npatch

Do l1=1,nband

if (nrootsmat(sint,l1,k)>0)then

Do l2=1,nband 
!------------------------------------------------------------------

l1r=nref(l1)
l2r=nref(l2)


!j1=imat(n3,l1,n1,l2,p3,k,p1)  !term4
!j2=imat(n2,l1,n4,l2,p2,k,p4)

!j3=imat(n1,l1,n3,l2,p1,k,p3)  !term2
!j4=imat(n4,l1,n2,l2,p4,k,p2)


if(sint==1)then

!Tauphc(i)=Tauphc(i)-1.0d+00*(1.0d+00/((1.0d+00)*(npatch)))&
!            *(((Yphcj1(l2,l1,k,n1,n3,p1,p3)*Yphcj2(l2,l1,k,p4,n4,n2,p2))*(Matsumph3(l2r,l1r,k)))&
!             +((Yphcj3(l2,l1,k,n1,n3,p1,p3)*Yphcj4(l2,l1,k,p4,n4,n2,p2))*(Matsumph1(l2r,l1r,k)))) 


Tauphc(i)=Tauphc(i)-1.0d+00*(1.0d+00/((1.0d+00)*(npatch)))&
            *((((Yphcj1(l2,l1,k,n1,n3,p1,p3)+ic*Yphcj1_im(l2,l1,k,n1,n3,p1,p3))&
               *(Yphcj2(l2,l1,k,p4,n4,n2,p2)+ic*Yphcj2_im(l2,l1,k,p4,n4,n2,p2)))*(Matsumph3(l2r,l1r,k)))&
             +(((Yphcj3(l2,l1,k,n1,n3,p1,p3)+ic*Yphcj3_im(l2,l1,k,n1,n3,p1,p3))&
               *(Yphcj4(l2,l1,k,p4,n4,n2,p2)+ic*Yphcj4_im(l2,l1,k,p4,n4,n2,p2)))*(Matsumph1(l2r,l1r,k))))

elseif(sint==2)then 

!Tauphc(i)=Tauphc(i)-1.0d+00*(1.0d+00/((1.0d+00)*(npatch)))&
!            *(((Yphcj1(l2,l1,k,n1,n3,p1,p3)*Yphcj2(l2,l1,k,p4,n4,n2,p2))*(Matsumph4(l2r,l1r,k)))&
!             +((Yphcj3(l2,l1,k,n1,n3,p1,p3)*Yphcj4(l2,l1,k,p4,n4,n2,p2))*(Matsumph2(l2r,l1r,k)))) 

Tauphc(i)=Tauphc(i)-1.0d+00*(1.0d+00/((1.0d+00)*(npatch)))&
            *((((Yphcj1(l2,l1,k,n1,n3,p1,p3)+ic*Yphcj1_im(l2,l1,k,n1,n3,p1,p3))&
               *(Yphcj2(l2,l1,k,p4,n4,n2,p2)+ic*Yphcj2_im(l2,l1,k,p4,n4,n2,p2)))*(Matsumph4(l2r,l1r,k)))&
             +(((Yphcj3(l2,l1,k,n1,n3,p1,p3)+ic*Yphcj3_im(l2,l1,k,n1,n3,p1,p3))&
               *(Yphcj4(l2,l1,k,p4,n4,n2,p2)+ic*Yphcj4_im(l2,l1,k,p4,n4,n2,p2)))*(Matsumph2(l2r,l1r,k))))


endif

!------------------------------------------------------------------
end do !l2
end if!
end do !l1
end do !k
end do !sint

!i=iphmat(n2,n4,p2,n3,n1,p3,p1)

!Tauphc(i)=tauphcmat(n2,n4,p2,n3,n1,p3,p1)

j1=imat(n2,n1,n3,n4,p2,p1,p3)
Tauphd(j1)=-Tauphc(i)
!------------------------------------------------------------------
end if !p4small
!------------------------------------------------------------------
!j1=imat(nbar_mat(n1),nbar_mat(n2),nbar_mat(n3),nbar_mat(n4),p1,p2,p3)
!Tauphc(j1)=Tauphc(i)

!j2=imat(nbar_mat(n2),nbar_mat(n1),nbar_mat(n3),nbar_mat(n4),p2,p1,p3)
!Tauphd(j2)=-Tauphc(j1)
!------------------------------------------------------------------
end do !n2
end do !n4
end do !p2
!------------------------------------------------------------------
end do !n3
end do !n1

end do !p1red

end do !p3
!------------------------------------------------------------------
!----------------------------------------------------------------------
!$OMP END PARALLEL
!----------------------------------------------------------------------
CALL SYSTEM_CLOCK(count_f_ph, count_rate_ph)
!----------------------------------------------------------------------
!print*,'time spent on calculating p-h terms=',real((count_f_ph-count_i_ph)/count_rate_ph)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
Taupp(1:NEQ)=cmplx(0.0d+00,0.0d+00)
YDvertex(1:NEQ)=0.0d+00
nrep(1:Ntot)=0.0d+00
!Tauppmat=0.0d+00
!------------------------------------------------------------------
CALL SYSTEM_CLOCK(count_i_pp, count_rate_pp)
!------------------------------------------------------------------
! Matsum-pp
!------------------------------------------------------------------
!$OMP PARALLEL SHARED(stepsize,x,eF,sgn,taupp,tauphc,tauphd,YDvertex,nrep,nref,nrefgap) , &
!$OMP& PRIVATE(n1,p1,n2,p2,l1,l2,p3max,energy,kpx,kpy,anglekp,radkp,step,sint,k,anglek,i,iasym,irot,dispersioncalll1) , &
!$OMP& SHARED(thetakmat,nrootsmat,gvelFSmat,thetap,KFSpmat,Yvertex,p4int,KFSmat)&
!$OMP& PRIVATE(Matsumpp,angle1,KFS1,angle2,KFS2,j1,j2,p3,n3,n4,p4,p1p,p2p,p3p,Matsumpp1,Matsumpp2,dispersioncall,l1r,l2r) &
!$OMP& SHARED(Yppj1,Yppj2,p4ppmat,ippmat,Yppj1_im,Yppj2_im) , &
!$OMP& PRIVATE(p1red)
!----------------------------------------------------------------------
!------------------------------------------------------------------
! Matsum-pp
!------------------------------------------------------------------
!------------------------------------------------------------------

Do p1red=1,npatch/num_threads

p1=p1red+(OMP_GET_THREAD_NUM())*(npatch/num_threads)

!print*, 'I am number',p1,'out of',num_threads*(npatch/num_threads)

Do p2=1,npatch

Do n1=1,nband,1

angle1=thetap(n1,p1)
KFS1=KFSpmat(n1,p1)

Do n2=1,n1

angle2=thetap(n2,p2)
KFS2=KFSpmat(n2,p2)


Do sint=1,2

energy=eF+(sgn(sint)*x)

Do k=1,npatch

Do l1=1,nband,nrefgap

anglek=thetakmat(sint,l1,k)

if (nrootsmat(sint,l1,k)>0)then

dispersioncalll1=dispersion(l1,KFSmat(sint,l1,k),thetakmat(sint,l1,k),eF)

!------------------------------------------------------------------ 

kpx=KFS1*cos(angle1)+KFS2*cos(angle2)-KFSmat(sint,l1,k)*cos(anglek)
kpy=KFS1*sin(angle1)+KFS2*sin(angle2)-KFSmat(sint,l1,k)*sin(anglek)

anglekp=atan2(kpy,kpx)
anglekp=(-SIGN(1.0d+00,anglekp)+1.0d+00)*(pi)+anglekp

radkp=sqrt(((kpx)**2)+((kpy)**2))

!------------------------------------------------------------------

Do l2=1,nband,nrefgap

dispersioncall=dispersion(l2,radkp,anglekp,eF)
!------------------------------------------------------------------ 
! Regulator:
!------------------------------------------------------------------
if(abs(dispersioncall)>x+stepsize)then
step=1.0d+00
elseif(abs(dispersioncall)<x-stepsize)then
step=0.0d+00
else
step=0.5d+00
endif
!------------------------------------------------------------------ 
!------------------------------------------------------------------
!------------------------------------------------------------------
!new:

if(energy-eF < 0.0d+00 .and. dispersioncall>0.0d+00)then
Matsumpp=0.0d+00
elseif(energy-eF .GE. 0.0d+00 .and. dispersioncall .LE. 0.0d+00)then
Matsumpp=0.0d+00
elseif(energy-eF < 0.0d+00)then
Matsumpp=1.0d+00/(dispersioncalll1+dispersioncall)
elseif(energy-eF .GE. 0.0d+00)then
Matsumpp=-1.0d+00/(dispersioncalll1+dispersioncall)
endif


if(abs(dispersioncalll1+dispersioncall)<10.0d+00**(-14.0d+00))then
!print*,'careful! devided by zero!'
Matsumpp=(beta/4.0d+00)*(((tanh(beta*(energy-eF)/2.0d+00))**2)-1.0d+00)
endif
!------------------------------------------------------------------

Matsumpp=step*Matsumpp*(((KFSmat(sint,l1,k))/(gvelFSmat(sint,l1,k))))

!------------------------------------------------------------------

if(sint==1)then
Matsumpp1(l2,l1,k)=Matsumpp
elseif(sint==2)then
Matsumpp2(l2,l1,k)=Matsumpp
endif

!------------------------------------------------------------------
!------------------------------------------------------------------
!------------------------------------------------------------------ 
end do !l2
!------------------------------------------------------------------
!------------------------------------------------------------------
!------------------------------------------------------------------
end if !k on FS
end do !k
end do !l1

end do !sint
!------------------------------------------------------------------
Do p3=1,nprot

Do n3=1,nband

Do n4=1,nband,2

p4=p4int(n1,n2,n3,n4,p1,p2,p3)

!i=imat(n1,n2,n3,n4,p1,p2,p3)

i=ippmat(n4,n3,p3,n2,n1,p2,p1)

!p4=p4ppmat(n4,n3,p3,n2,n1,p2,p1)

if(p4==0)then

Taupp(i)=cmplx(0.0d+00,0.0d+00)

elseif(p4>0)then

!------------------------------------------------------------------
Do sint=1,2

Do k=1,npatch

Do l1=1,nband

if (nrootsmat(sint,l1,k)>0)then

Do l2=1,nband 

l1r=nref(l1)
l2r=nref(l2)

if(sint==1)then

!Taupp(i)=Taupp(i)-(1.0d+00/((1.0d+00)*(npatch)))*Yppj1(l2,l1,k,n2,n1,p2,p1)*Yppj2(l2,l1,k,p4,n4,n3,p3)&
!                 *Matsumpp1(l2r,l1r,k)


Taupp(i)=Taupp(i)-(1.0d+00/((1.0d+00)*(npatch)))*(Yppj1(l2,l1,k,n2,n1,p2,p1)+ic*Yppj1_im(l2,l1,k,n2,n1,p2,p1))&
                                                *(Yppj2(l2,l1,k,p4,n4,n3,p3)+ic*Yppj2_im(l2,l1,k,p4,n4,n3,p3))&
                                                *(Matsumpp1(l2r,l1r,k))

elseif(sint==2)then

!Taupp(i)=Taupp(i)-(1.0d+00/((1.0d+00)*(npatch)))*Yppj1(l2,l1,k,n2,n1,p2,p1)*Yppj2(l2,l1,k,p4,n4,n3,p3)&
!                 *Matsumpp2(l2r,l1r,k)

Taupp(i)=Taupp(i)-(1.0d+00/((1.0d+00)*(npatch)))*(Yppj1(l2,l1,k,n2,n1,p2,p1)+ic*Yppj1_im(l2,l1,k,n2,n1,p2,p1))&
                                                *(Yppj2(l2,l1,k,p4,n4,n3,p3)+ic*Yppj2_im(l2,l1,k,p4,n4,n3,p3))&
                                                *(Matsumpp2(l2r,l1r,k))


endif


end do  !l2

end if !nroots

end do  !k

end do  !l1

end do !sint
!------------------------------------------------------------------
endif
!------------------------------------------------------------------
!i=ippmat(n3,n4,p3,n2,n1,p1,p2)
!Taupp(i)=tauppmat(n3,n4,p3,n2,n1,p1,p2)

!j1=imat(nbar_mat(n1),nbar_mat(n2),nbar_mat(n3),nbar_mat(n4),p1,p2,p3)
!Taupp(j1)=Taupp(i)
!------------------------------------------------------------------
end do !p3

end do !n4
end do !n3
!------------------------------------------------------------------
end do !n2
end do !n1

end do !p2
end do !p1red
!----------------------------------------------------------------------
!$OMP END PARALLEL
!----------------------------------------------------------------------
!----------------------------------------------------------------------
CALL SYSTEM_CLOCK(count_f_pp, count_rate_pp)
!----------------------------------------------------------------------
! Summing all terms:
!------------------------------------------------------------------
!$OMP PARALLEL SHARED(stepsize,x,eF,sgn,taupp,tauphc,tauphd,YDvertex,nrep,nref,nrefgap) , &
!$OMP& PRIVATE(n1,p1,n2,p2,l1,l2,p3max,energy,kpx,kpy,anglekp,radkp,step,sint,k,anglek,i,iasym,irot,dispersioncalll1) , &
!$OMP& SHARED(thetakmat,nrootsmat,gvelFSmat,thetap,KFSpmat,Yvertex,p4int,KFSmat)&
!$OMP& PRIVATE(Matsumpp,angle1,KFS1,angle2,KFS2,j1,j2,p3,n3,n4,p4,p1p,p2p,p3p,Matsumpp1,Matsumpp2,dispersioncall,l1r,l2r) &
!$OMP& SHARED(Yppj1,Yppj2,p4ppmat,ippmat,nbar_mat) , &
!$OMP& PRIVATE(p1red)
!----------------------------------------------------------------------
!------------------------------------------------------------------

Do p1red=1,npatch/num_threads

p1=p1red+(OMP_GET_THREAD_NUM())*(npatch/num_threads)

Do p3=1,nprot

Do p2=1,npatch

Do n4=1,nband,2


Do n3=1,nband


Do n1=1,nband,1


Do n2=1,n1!nband,1

p4=p4int(n1,n2,n3,n4,p1,p2,p3)

!i=imat(n1,n2,n3,n4,p1,p2,p3)

i=ippmat(n4,n3,p3,n2,n1,p2,p1)

!p4=p4ppmat(n4,n3,p3,n2,n1,p2,p1)

if(p4==0)then

YDvertex(i)=0.0d+00

elseif(p4>0)then

!-----------------------------------------------------------------
!------------------------------------------------------------------
!------------------------------------------------------------------
!i=ippmat(n3,n4,p3,n2,n1,p1,p2)
!Taupp(i)=tauppmat(n3,n4,p3,n2,n1,p1,p2)

!j1=imat(nbar_mat(n1),nbar_mat(n2),nbar_mat(n3),nbar_mat(n4),p1,p2,p3)
!Taupp(j1)=Taupp(i)


if(nrep(i)==0)then

YDvertex(i)=real(Taupp(i)+Tauphc(i)+Tauphd(i))/(2.0d+00*pi)

if(vertex_flag=='comp')then

YDvertex(i+(NEQ/2))=aimag(Taupp(i)+Tauphc(i)+Tauphd(i))/(2.0d+00*pi)

endif

nrep(i)=1

!------------------------------------------------------------------
! recnstruction using SU(2) sym in spin-sector:
!------------------------------------------------------------------

j2=imat(nbar_mat(n1),nbar_mat(n2),nbar_mat(n3),nbar_mat(n4),p1,p2,p3)
YDvertex(j2)=YDvertex(i)

if(vertex_flag=='comp')then
YDvertex(j2+(NEQ/2))=YDvertex(i+(NEQ/2))
endif

nrep(j2)=1

!------------------------------------------------------------------
!------------------------------------------------------------------
! vertex recostruction employing anti-symmetric properties:
!------------------------------------------------------------------
if(Idantisym==1)then

iasym=imat(n2,n1,n3,n4,p2,p1,p3)

if(nrep(iasym)==0 .and. n1>n2)then

YDvertex(iasym)=-YDvertex(i)                      ! real part reconstruction 

if(vertex_flag=='comp')then
YDvertex(iasym+(NEQ/2))=-YDvertex(i+(NEQ/2))      ! imaginary part reconstruction 
endif

nrep(iasym)=1


!------------------------------------------------------------------
! recnstruction using SU(2) sym in spin-sector:
!------------------------------------------------------------------

j2=imat(nbar_mat(n2),nbar_mat(n1),nbar_mat(n3),nbar_mat(n4),p2,p1,p3)

YDvertex(j2)=YDvertex(iasym)

if(vertex_flag=='comp')then
YDvertex(j2+(NEQ/2))=YDvertex(iasym+(NEQ/2))
endif

nrep(j2)=1

!------------------------------------------------------------------

endif

endif

!----------------------------------------------------------------------
!----------------------------------------------------------------------
! vertex recostruction due to the rot-sym of Fermi surface!
!----------------------------------------------------------------------

if(idFSrotsym==1)then

Do j1=1,nprotpatch-1

p1p=p1+j1*nprot-npatch*int((p1+(j1*nprot)-10.0d+00**(-10.0d+00))/(1.0d+00*npatch))
p2p=p2+j1*nprot-npatch*int((p2+(j1*nprot)-10.0d+00**(-10.0d+00))/(1.0d+00*npatch))
p3p=p3+j1*nprot-npatch*int((p3+(j1*nprot)-10.0d+00**(-10.0d+00))/(1.0d+00*npatch))

irot=imat(n1,n2,n3,n4,p1p,p2p,p3p)

if(nrep(irot)==0)then

YDvertex(irot)=YDvertex(i)                   ! real part reconstruction

if(vertex_flag=='comp')then
YDvertex(irot+(NEQ/2))=YDvertex(i+(NEQ/2))   ! imaginary part reconstruction
endif

nrep(irot)=1


!------------------------------------------------------------------
! recnstruction using SU(2) sym in spin-sector:
!------------------------------------------------------------------

j2=imat(nbar_mat(n1),nbar_mat(n2),nbar_mat(n3),nbar_mat(n4),p1p,p2p,p3p)
YDvertex(j2)=YDvertex(irot)

if(vertex_flag=='comp')then
YDvertex(j2+(NEQ/2))=YDvertex(irot+(NEQ/2))
endif

nrep(j2)=1

!------------------------------------------------------------------


endif

if(Idantisym==1)then

iasym=imat(n2,n1,n3,n4,p2p,p1p,p3p)

if(nrep(iasym)==0)then

YDvertex(iasym)=-YDvertex(i)                   ! real part reconstruction

if(vertex_flag=='comp')then
YDvertex(iasym+(NEQ/2))=-YDvertex(i+(NEQ/2))   ! imaginary part reconstruction
endif

nrep(iasym)=1


!------------------------------------------------------------------
! recnstruction using SU(2) sym in spin-sector:
!------------------------------------------------------------------

j2=imat(nbar_mat(n2),nbar_mat(n1),nbar_mat(n3),nbar_mat(n4),p2p,p1p,p3p)
YDvertex(j2)=YDvertex(iasym)

if(vertex_flag=='comp')then
YDvertex(j2+(NEQ/2))=YDvertex(iasym+(NEQ/2))
endif

nrep(j2)=1

!------------------------------------------------------------------


endif

endif


!end if ! p1,p2,p3 inside BZ!

end do !j1

endif
!----------------------------------------------------------------------
endif !if (nrep(i)==0)
endif !p4 small
!------------------------------------------------------------------
if(vertex_flag=='comp')then
maxvertex_n=sqrt((Yvertex(i))**2+(Yvertex(i+(NEQ/2)))**2)
else
maxvertex_n=Yvertex(i)
endif
!------------------------------------------------------------------
if(i==1)then
maxvertex=maxvertex_n
imax=1

elseif(abs(maxvertex_n)>abs(maxvertex))then
maxvertex=maxvertex_n
imax=i
endif

if(i==1)then
minvertex=maxvertex_n
imin=1
elseif(abs(maxvertex_n)<abs(minvertex))then
minvertex=maxvertex_n
imin=i
endif

!-----------------------------------------------------------------
end do !p3

end do !n4
end do !n3
!------------------------------------------------------------------
end do !n2
end do !n1

end do !p2
end do !p1red
!----------------------------------------------------------------------
!$OMP END PARALLEL
!----------------------------------------------------------------------
!print*,'time spent for pp term and summing up everything is=',real((count_f_pp-count_i_pp)/count_rate_pp)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
xcount=0
!----------------------------------------------------------------------
!----------------------------------------------------------------------
end if !nrootstotal>0
!----------------------------------------------------------------------
!----------------------------------------------------------------------
print*,'icall=',icall,'scale/t=',scalecommon/t,'maxV=',maxvertex,'minV=',minvertex
if(icall>100 .and. abs(maxvertex)>max_vertex_required)then
print*,'vertex is diverging!, i.e., I am setting divflag to 1'
divflag=1
YDvertex=0.0d+00
endif
!----------------------------------------------------------------------
!----------------------------------------------------------------------
DEALLOCATE(nrep,Matsumpp1,Matsumpp2,Matsumph1,Matsumph2)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
end subroutine Fflow
!----------------------------------------------------------------------
!----------------------------------------------------------------------
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
DOUBLE PRECISION function dispersionana(blab,mom,theta,energy)

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

dispersionana=0.0d+00

if(blab .LE. (nband/2.0d+00))then        

dispersionana=1.0*t*sqrt(3.0d+00+2.0d+00*(cos(Lconst*kx*(sqrt(3.0d+00))))&
                +4.0d+00*(Cos(Lconst*kx*(sqrt(3.0d+00))/2.0d+00))*(cos(Lconst*1.5d+00*ky)))-energy

else 

dispersionana=-1.0d+00*t*sqrt(3.0d+00+2.0d+00*(cos(Lconst*kx*(sqrt(3.0d+00))))&
                +4.0d+00*(Cos(Lconst*kx*(sqrt(3.0d+00))/2.0d+00))*(cos(Lconst*1.5d+00*ky)))-energy


endif
!----------------------------------------------------------------------
!----------------------------------------------------------------------
end function dispersionana
!----------------------------------------------------------------------
!----------------------------------------------------------------------
DOUBLE PRECISION function dispersion(blab,mom,theta,energy)

use commons

implicit none

common eF,kxsamp,kysamp,dissamp

DOUBLE PRECISION,INTENT(IN)::mom,theta,energy
INTEGER, INTENT(IN)::blab
DOUBLE PRECISION::eF,scalecommon,kx,ky,eF_min,eF_max
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
INTEGER::kxint,kyint,xint,yint,signnx,signny,id,countx,county,idp1,idp2,idp3,idp4,idp0,&
         flag4,flag0,i,kyintmin,xintp2fold,yintp2fold,xintp3fold,yintp3fold
DOUBLE PRECISION,DIMENSION(1:Ndis)::kxsamp,kysamp
DOUBLE PRECISION,DIMENSION(1:nband,1:Ndis,1:Ndis)::dissamp

INTEGER,DIMENSION(1:Ndis)::kyintmaxmat,kyintminmat

DOUBLE PRECISION::px0,px1,px2,px3,px4,px2fold,px3fold,&
                  py0,py1,py2,py3,py4,py2fold,py3fold,&
                  pz0,pz1,pz2,pz3,pz4,&
                  nx,ny,nz,Ax,Ay,Az,Bx,By,Bz,nl,&
                  anglep1,anglep2,anglep3,KFSp1,KFSp2,KFSp3,dispersionana,dispersion_test&
                  ,kpx,kpy

!----------------------------------------------------------------------
!Do kxint=1,Ndis
!Do kyint=1,Ndis
!write(2,*)kxsamp(kxint),kysamp(kyint),dissamp(1,kxint,kyint)
!end do
!end do
!----------------------------------------------------------------------
! for Hexagonal lattice:
!----------------------------------------------------------------------
kx=mom*cos(theta)
ky=mom*sin(theta)

CALL FOLDING(kx,ky,kx,ky)

!----------------------------------------------------------------------
if(disp_flag=='ana')then

dispersion=dispersionana(blab,mom,theta,energy)

!----------------------------------------------------------------------
elseif(disp_flag=='try')then

dispersion=dispersion_test(blab,mom,theta,energy)

!----------------------------------------------------------------------
!----------------------------------------------------------------------
elseif(disp_flag=='num')then
!----------------------------------------------------------------------
!xint=1
!countx=0
!------------------------------------------------------------------
!------------------------------------------------------------------
!if(kyintmaxmat(xint)-kyintmin>1)then
!------------------------------------------------------------------

CALL Grid_search(kx,ky,xint,yint) 

!print*,'inside function:',xint,yint,kx,ky,kxsamp(xint),kysamp(yint),dissamp(1,xint,yint)

!------------------------------------------------------------------
!------------------------------------------------------------------

px2=kxsamp(xint+1)
py2=kysamp(yint)

pz2=dissamp(blab,xint+1,yint)!-energy

!------------------------------------------------------------------

px3=kxsamp(xint)
py3=kysamp(yint+1)

pz3=dissamp(blab,xint,yint+1)!-energy

!------------------------------------------------------------------

if(ky .LE. ((((kysamp(yint+1)-kysamp(yint))/(kxsamp(xint)-kxsamp(xint+1)))&
             *(kx-kxsamp(xint)))+(kysamp(yint+1))))then

!print*,'triangle1'

px1=kxsamp(xint)
py1=kysamp(yint)
pz1=dissamp(blab,xint,yint)!-energy

flag0=1
flag4=0

!print*,'here-0'

else

px1=kxsamp(xint+1)
py1=kysamp(yint+1)
pz1=dissamp(blab,xint+1,yint+1)!-energy

!print*,'here-4'

flag0=0
flag4=1

!print*,'triangle2'

endif

!------------------------------------------------------------------
!-------------------------------------

Ax=px2-px1
Ay=py2-py1
Az=pz2-pz1

Bx=px3-px1
By=py3-py1
Bz=pz3-pz1


nx=Ay*Bz-Az*By
ny=-Ax*Bz+Az*Bx
nz=Ax*By-Ay*Bx

!nl=(sqrt(nx**2+ny**2+nz**2))

!nx=nx/nl
!ny=ny/nl
!nz=nz/nl

!------------------------------------------------------------------
! Linear interpolation to estimate the dispersion:
!------------------------------------------------------------------

if(abs(nz)>10.0d+00**(-8.0d+00))then

dispersion=pz1-(nx/nz)*(kx-px1)-(ny/nz)*(ky-py1)-energy

else


print*,'ERROR!!!! from dispersion function',kx/BZconst,ky/BZconst,Ax,By,Ay,Bx,xint,yint!,kyintmaxmat(xint)

!print*,'checking=',px1/BZconst,px2/BZconst,px3/BZconst,kxsamp(xint),kxsamp(xint+1)

endif
!----------------------------------------------------------------------
!----------------------------------------------------------------------
end if !disp_flag
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!print*,'I am in dispersion with flag=',disp_flag
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
!--------------------------------------
! Numerical:
!--------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Numerical calculation of the group velocity:
!----------------------------------------------------------------------
!----------------------------------------------------------------------

ek0=dispersion(blab,mom,theta,eF)

!--------------------------------------
!--------------------------------------

deltakx=1.0d+00*pi/(10.0d+00**(5.0d+00))
deltaky=1.0d+00*pi/(10.0d+00**(5.0d+00))


kx=(deltakx/2.0d+00)+mom*cos(theta)
ky=mom*sin(theta)

kradkxp=sqrt(((kx)**2)+((ky)**2))

anglekxp=atan2(ky,kx)
anglekxp=(-SIGN(1.0d+00,anglekxp)+1.0d+00)*(pi)+anglekxp

ekxp=dispersion(blab,kradkxp,anglekxp,eF)

!---------------

kx=(-deltakx/2.0d+00)+mom*cos(theta)
ky=mom*sin(theta)

kradkxn=sqrt(((kx)**2)+((ky)**2))

anglekxn=atan2(ky,kx)
anglekxn=(-SIGN(1.0d+00,anglekxn)+1.0d+00)*(pi)+anglekxn

ekxn=dispersion(blab,kradkxn,anglekxn,eF)


groupvelkx=(ekxp-ekxn)/(deltakx)

!--------------------------------------

kx=mom*cos(theta)
ky=(deltaky/2.0d+00)+mom*sin(theta)

kradkyp=sqrt(((kx)**2)+((ky)**2))

anglekyp=atan2(ky,kx)
anglekyp=(-SIGN(1.0d+00,anglekyp)+1.0d+00)*(pi)+anglekyp

ekyp=dispersion(blab,kradkyp,anglekyp,eF)

!---------------
kx=mom*cos(theta)
ky=(-deltaky/2.0d+00)+mom*sin(theta)

kradkyn=sqrt(((kx)**2)+((ky)**2))

anglekyn=atan2(ky,kx)
anglekyn=(-SIGN(1.0d+00,anglekyn)+1.0d+00)*(pi)+anglekyn

ekyn=dispersion(blab,kradkyn,anglekyn,eF)


groupvelky=(ekyp-ekyn)/(deltaky)

!---------------------------------------------------------------------- 

grouprad=(groupvelkx)*(cos(theta))+(groupvelky)*(sin(theta))
groupangle=(groupvelky)*(cos(theta))-(groupvelkx)*(sin(theta))

!---------------------------------------------------------------------- 
!--------------------------------------
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
DOUBLE PRECISION:: fermi,dispersion,groupvel,dispersionana
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
kzerop=t*((10.0d+00)**(-5.0d+00))

zeroer = error_roots_disp!t*(10.0d+00**(-4.0d+00))
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
if(((dispersionana(blabb,zeromin,angle,energy))*dispersionana(blabb,zeromax,angle,energy))<0.0d+00)then

Do

call zero_rc ( zeromin, zeromax, zeroer, zeroarg, zerostatus, zerovalue )

if ( zerostatus == 0 ) then
exit 
end if

zerovalue=dispersionana(blabb,zeroarg,angle,energy)


end do

kFS=zeroarg
nroots=1

elseif(abs(dispersionana(blabb,zeromin,angle,energy))< zeroer)then

kFS=zeromin
nroots=1

elseif(abs(dispersionana(blabb,zeromax,angle,energy))< zeroer)then

kFS=zeromax
nroots=1

elseif(((dispersionana(blabb,zeromin,angle,energy))*(dispersionana(blabb,zeromax,angle,energy)))>0.0d+00)then

nroots=0


endif

!----------------------------------------------------------------------

if(nroots>0 .and. abs(dispersionana(blabb,kFS,angle,energy))>10.0d+00*zeroer)then
print*,'Note: roots of dispersion has not been found properly!',dispersionana(blabb,kFS,angle,energy),angle/pi,eF
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
subroutine patch_lab_Ge(flag,theta,mom,blab,theta_min_mat,theta_max_mat,nplab,anglek)

use commons

implicit none


CHARACTER*3, INTENT (IN)::flag
DOUBLE PRECISION, INTENT(IN)::theta,mom
INTEGER, INTENT (IN)::blab

DOUBLE PRECISION, DIMENSION(1:nband), INTENT(IN)::theta_min_mat,theta_max_mat

INTEGER, INTENT(OUT)::nplab
DOUBLE PRECISION, INTENT(OUT)::anglek


DOUBLE PRECISION::kx,ky,kxprime,kyprime,base,KFSmax,kxmax,kymax,thetar,kred!,patcherr1,patcherr2
INTEGER::nplab0,nplabm,nplabp,BZid,np1,np2
!----------------------------------------------------------------------
!patcherr1=10.0d+00**(-15.0d+00)
!patcherr2=10.0d+00**(-14.0d+00)
!----------------------------------------------------------------------


if(flag=='mul')then

np1=int((theta-patcherr2)/(pi/3.0d+00))

thetar=theta-(pi/3.0d+00)*np1-patcherr2

!print*,'patch-check:','thetar=',thetar,'thetamin',theta_min_mat(blab),'thetamax',theta_max_mat(blab)

np2=1

if(npatch==12)then

CALL kred_generator('test',num_theta_ired,blab,thetar,kred)

if(mom>kred)then

np2=np2+1

endif


else


if(thetar<theta_min_mat(blab))then

np2=1


elseif(thetar .GE. theta_max_mat(blab))then

np2=npatch/6

elseif((npatch/6)-2>0)then

np2=2+2*int((thetar-(theta_min_mat(blab)))/((theta_max_mat(blab)-theta_min_mat(blab))/((npatch/12)-1)))

CALL kred_generator('test',num_theta_ired,blab,thetar,kred)


if(mom>kred)then

np2=np2+1

endif

endif

endif

nplab=np1*(npatch/6)+np2

!print*,'checkin-patch:',np1,np2,nplab


!print*,'from patch_lab_Ge:',theta/(pi/3.0d+00),thetar/(pi/3.0d+00),np1,np2,nplab

elseif(flag=='sim')then

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
elseif (flag=='com')then 

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
!----------------------------------------------------------------------
INTEGER function imat(n1,n2,n3,n4,p1,p2,p3)

use commons

implicit none

integer::n1,n2,n3,n4,p1,p2,p3

imat=(n1-1)*((nband**3)*(npatch**3))+(n2-1)*((nband**2)*(npatch**3))&
     +(n3-1)*((nband)*(npatch**3))+(n4-1)*(npatch**3)&
     +(p1-1)*(npatch**2)+(p2-1)*(npatch)+p3

end function imat
!----------------------------------------------------------------------
!--------------------------------------------------------------
Subroutine num_unique(Ntot,v_test,error,ucount)

implicit none

integer,intent(in)::Ntot
double precision,dimension(1:Ntot),intent(in) :: v_test
double precision,intent(in)::error
integer,intent(out)::ucount

integer::dim_red
double precision,dimension(:),allocatable::vec,vec_red

ucount=1
dim_red=size(v_test)

Do

if(ucount==1)then
allocate(vec(1:size(v_test)))
vec=v_test
else
deallocate(vec)
allocate(vec(1:size(vec_red)))
vec=vec_red
deallocate(vec_red)
endif

dim_red=size(pack(vec,abs(abs(vec)-abs(vec(1)))>error))


allocate(vec_red(1:dim_red))

vec_red=pack(vec,abs(abs(vec)-abs(vec(1)))>error)

if(dim_red==0)exit

ucount=ucount+1

end do

end subroutine num_unique
!----------------------------------------------------------------------
SUBROUTINE Grid_search(kx,ky,xint,yint)
!----------------------------------------------------------------------
use commons

implicit none

common eF,kxsamp,kysamp,dissamp

!----------------------------------------------------------------------
DOUBLE PRECISION,INTENT(IN)::kx,ky
INTEGER,INTENT(OUT)::xint,yint
!----------------------------------------------------------------------
DOUBLE PRECISION,DIMENSION(1:Ndis)::kxsamp,kysamp
DOUBLE PRECISION,DIMENSION(1:nband,1:Ndis,1:Ndis)::dissamp
INTEGER,DIMENSION(1:Ndis)::kyintmaxmat,kyintminmat
!----------------------------------------------------------------------
INTEGER::kyintmin,countx,county
DOUBLE PRECISION::eF,eF_min,eF_max
!----------------------------------------------------------------------
!----------------------------------------------------------------------

xint=1
countx=0

Do 

if(xint==Ndis-1)then

exit

elseif(kx .GE. kxsamp(xint) .and. kx< kxsamp(xint+1)) then

if(kx<0)then
xint=xint+1
endif

exit

elseif(xint<Ndis)then

xint=xint+1

endif

countx=countx+1

if(countx>Ndis)then
print*,'I ma stuck in loop-x!'
endif

end do

!------------------------------------------------------------------

yint=1
county=0

Do 

if(yint==Ndis-1)then

exit

elseif(ky .GE. kysamp(yint) .and. ky< kysamp(yint+1)) then

if(ky<0)then
yint=yint+1
endif

exit

elseif(yint<Ndis)then

yint=yint+1

endif

county=county+1

if(county>Ndis)then
print*,'I ma stuck in loop-y!'
endif

end do

!----------------------------------------------------------------------
!----------------------------------------------------------------------
END SUBROUTINE 
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Subroutine to find the roots of the dipersion:
!----------------------------------------------------------------------
subroutine roots_dispersion_theta(blab,angle,energy,KFSmin,KFSmax,nroots,KFS)

use commons

implicit none

common eF

!------------------------------------------------------------------ 
! Parameter declaration:
!------------------------------------------------------------------ 
DOUBLE PRECISION,INTENT(IN)::angle,energy
INTEGER,INTENT(IN)::blab
DOUBLE PRECISION,INTENT(IN)::KFSmin,KFSmax
DOUBLE PRECISION,INTENT(OUT)::KFS
INTEGER,INTENT(OUT)::nroots

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
DOUBLE PRECISION::gvsign,kxr,kyr,kxbound,kybound,zeroer
DOUBLE PRECISION,DIMENSION(1:20)::KFSext,KFSmat

INTEGER::flag
DOUBLE PRECISION::dispersion_theta

!------------------------------------------------------------------
!------------------------------------------------------------------
kzerop=t*((10.0d+00)**(-5.0d+00))

zeroer = error_roots_disp!t*(10.0d+00**(-3.0d+00))
!----------------------------------------------------------------------            
!Initialization:
!----------------------------------------------------------------------

blabb=blab! important checking!!!!

!print*,'blab=',blab

!----------------------------------------------------------------------
! Initialization:
!----------------------------------------------------------------------

zeromin = KFSmin
zeromax = KFSmax

zerostatus=0
nroots=-1
KFS=0.0d+00

icountersp=0

!----------------------------------------------------------------------
if(((dispersion_theta(blabb,zeromin,angle,energy))*dispersion_theta(blabb,zeromax,angle,energy))<0.0d+00)then

Do

call zero_rc ( zeromin, zeromax, zeroer, zeroarg, zerostatus, zerovalue )

if ( zerostatus == 0 ) then
exit 
end if

zerovalue=dispersion_theta(blabb,zeroarg,angle,energy)


end do

kFS=zeroarg
nroots=1

elseif(abs(dispersion_theta(blabb,zeromin,angle,energy))< zeroer)then

kFS=zeromin
nroots=1

elseif(abs(dispersion_theta(blabb,zeromax,angle,energy))< zeroer)then

kFS=zeromax
nroots=1

elseif(((dispersion_theta(blabb,zeromin,angle,energy))*(dispersion_theta(blabb,zeromax,angle,energy)))>0.0d+00)then

nroots=0


endif

!----------------------------------------------------------------------

if(nroots>0 .and. abs(dispersion_theta(blabb,kFS,angle,energy))>10.0d+00*zeroer)then
print*,'Note: roots of dispersion has not been found properly!',dispersion_theta(blabb,kFS,angle,energy),angle/pi,eF
endif

if(nroots .LE. 0)then
KFS=0.0d+00
end if

!print*, 'Inside:', dispersion(blabb,KFS,angle,energy),blabb,KFS,angle

!-----------------------------------------------------------------------
end subroutine roots_dispersion_theta
!-----------------------------------------------------------------------
!----------------------------------------------------------------------
! Subroutine to find the roots of the dipersion:
!----------------------------------------------------------------------
subroutine roots_dispersion_theta2(blab,angle,energy,KFSmin,KFSmax,nroots,KFS)

use commons

implicit none

common eF

!------------------------------------------------------------------ 
! Parameter declaration:
!------------------------------------------------------------------ 
DOUBLE PRECISION,INTENT(IN)::angle,energy
INTEGER,INTENT(IN)::blab
DOUBLE PRECISION,INTENT(IN)::KFSmin,KFSmax
DOUBLE PRECISION,INTENT(OUT)::KFS
INTEGER,INTENT(OUT)::nroots

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
DOUBLE PRECISION::gvsign,kxr,kyr,kxbound,kybound,zeroer
DOUBLE PRECISION,DIMENSION(1:20)::KFSext,KFSmat

INTEGER::flag
DOUBLE PRECISION::dispersion_theta

!------------------------------------------------------------------
!------------------------------------------------------------------
kzerop=t*((10.0d+00)**(-5.0d+00))

zeroer = error_roots_disp!t*(10.0d+00**(-3.0d+00))
!----------------------------------------------------------------------            
!Initialization:
!----------------------------------------------------------------------

blabb=blab! important checking!!!!

!print*,'blab=',blab

!----------------------------------------------------------------------
! Initialization:
!----------------------------------------------------------------------

zeromin = KFSmin
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
end subroutine roots_dispersion_theta2
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
Subroutine Extermum_theta(blab,angle,dimwork,nsign,KFSext,flag,e_ext,e_bound)

use commons

implicit none


INTEGER,INTENT(IN)::blab,dimwork
DOUBLE PRECISION,INTENT(IN)::angle
INTEGER,INTENT(OUT)::nsign
DOUBLE PRECISION,DIMENSION(1:10)::KFSext,dispext
INTEGER,INTENT(OUT)::flag
DOUBLE PRECISION,INTENT(OUT)::e_ext,e_bound

DOUBLE PRECISION,DIMENSION(1:dimwork)::dispmat_ex,ddispmat_ex

DOUBLE PRECISION::dispersion_theta,KFS,KFSmax,dispersion,angle_min,angle_max
INTEGER::i,ncrit,j,openstatus,k,unitm,i_min,i_max

DOUBLE PRECISION,DIMENSION(1:Nnumdisp)::kxdftmat,kydftmat,dispdftmat
DOUBLE PRECISION,DIMENSION(1:Nnumdisp)::kpxmat_interp,kpymat_interp,dispmat_interp
DOUBLE PRECISION::p_interp,dispcheck,kpxmin,kpx,kpy
!------------------------------------------------------------------

i_min=0
i_max=0

ncrit=0
nsign=1

CALL BZ_Edge(angle,KFSmax)

Do i=1,50

KFS=KFSmax*(i)/50.0d+00

!dispmat_ex(i)=dispersion(blab,KFS,angle,0.0d+00)

!p_interp=4.0d+00

!kpxmat_interp(1)=KFS*(cos(angle))
!kpymat_interp(1)=KFS*(sin(angle))

!CALL shepard_interp_2d(Nnumdisp,kxdftmat,kydftmat,dispdftmat,p_interp,1,kpxmat_interp,kpymat_interp,dispmat_interp)

dispmat_ex(i)=dispersion_theta(blab,KFS,angle,0.0d+00)!dispmat_interp(1)


if(i>1)then

ddispmat_ex(i)=(dispmat_ex(i)-dispmat_ex(i-1))/(KFSmax/50.0d+00)

!write(3,*)kpxmat_interp(1)/BZconst,kpymat_interp(1)/BZconst,dispmat_ex(i),angle/(2.0d+00),KFS/BZconst,ddispmat_ex(i),nsign

endif

end do

!------------------------------------------------------------------
Do i=2,49

if(i>1)then

KFS=KFSmax*i/50.0d+00

!write(2,*)KFS/BZconst,dispmat_ex(i),ddispmat_ex(i),ddispmat_ex(i)*ddispmat_ex(i+1)

if(ddispmat_ex(i)*ddispmat_ex(i+1)<0.0d+00)then

Do j=1,2

if(ddispmat_ex(i)*ddispmat_ex(i+1+j)<0.0d+00)then

ncrit=ncrit+1

endif  !changing sign with j add


if(ncrit>1)then

KFSext(nsign)=KFS
dispext(nsign)=dispmat_ex(i)

nsign=nsign+1

endif  !ncrit

enddo

else

ncrit=0

endif  !changing sign

endif !i>1


if(nsign>5)then
flag=0
else
flag=1
end if


if(nsign==1)then

CALL BZ_edge(angle,KFS)

KFSext(1)=KFS
e_ext=dispmat_ex(50)
e_bound=dispmat_ex(50)

elseif(nsign==2)then

KFSext(2)=KFS
e_ext=dispext(1)
e_bound=dispmat_ex(50)

elseif(nsign>2)then

e_ext=dispext(1)
e_bound=dispext(nsign-1)

endif

end do !i


if(abs(e_ext-e_bound)>10.0d+00**(-10.0d+00) .and. i_min==0)then

angle_min=angle
i_min=1

elseif(abs(e_ext-e_bound)<10.0d+00**(-10.0d+00) .and. i_min==1 .and. i_max==0)then

angle_max=angle
i_max=1

end if

kpx=KFSext(1)*(cos(angle))
kpy=KFSext(1)*(sin(angle))

!------------------------------------------------------------------
end Subroutine Extermum_theta
!----------------------------------------------------------------------
!-----------------------------------------------------------------------
Subroutine Extermum_theta_bunch(blab,dimwork,n_a,flag,angle_mat,kred_mat,e_ext_mat,e_bound_mat)

use commons

implicit none

!-----------------------------------------------------------------------
INTEGER,INTENT(IN)::blab
INTEGER,INTENT(IN)::dimwork
INTEGER,INTENT(IN)::n_a

INTEGER,INTENT(OUT)::flag
DOUBLE PRECISION,DIMENSION(1:n_a),INTENT(OUT)::angle_mat
DOUBLE PRECISION,DIMENSION(1:n_a),INTENT(OUT)::kred_mat
DOUBLE PRECISION,DIMENSION(1:n_a),INTENT(OUT)::e_ext_mat
DOUBLE PRECISION,DIMENSION(1:n_a),INTENT(OUT)::e_bound_mat
!-----------------------------------------------------------------------
DOUBLE PRECISION,DIMENSION(1:10)::KFSext
DOUBLE PRECISION::e_ext
DOUBLE PRECISION::e_bound

DOUBLE PRECISION::angle
DOUBLE PRECISION,DIMENSION(1:10)::dispext

DOUBLE PRECISION,DIMENSION(1:dimwork)::dispmat_ex,ddispmat_ex

DOUBLE PRECISION::dispersion_theta,KFS,KFSmax,dispersion,angle_min,angle_max
INTEGER::i,ncrit,j,openstatus,k,unitm,i_min,i_max,nsign

DOUBLE PRECISION,DIMENSION(1:Nnumdisp)::kxdftmat,kydftmat,dispdftmat
DOUBLE PRECISION,DIMENSION(1:Nnumdisp)::kpxmat_interp,kpymat_interp,dispmat_interp
DOUBLE PRECISION::p_interp,dispcheck,kpxmin,kpx,kpy,dispmin,dispmax
!------------------------------------------------------------------
!------------------------------------------------------------------

CALL READ_DATA(blab,Nnumdisp,kxdftmat,kydftmat,dispdftmat,dispmin,dispmax)

kxdftmat(1:Nnumdisp)=BZconst*kxdftmat(1:Nnumdisp)/(abs(kpxmin))!(0.03855004306107383d+00) 
kydftmat(1:Nnumdisp)=BZconst*kydftmat(1:Nnumdisp)/(abs(kpxmin))!(0.03855004306107383d+00)


!------------------------------------------------------------------

i_min=0
i_max=0

Do k=1,n_a

angle=(k-1.0d+00)*(pi/3.0d+00)/n_a


ncrit=0
nsign=1

CALL BZ_Edge(angle,KFSmax)

Do i=1,dimwork

KFS=KFSmax*(i)/dimwork


!p_interp=4.0d+00

!kpxmat_interp(1)=KFS*(cos(angle))
!kpymat_interp(1)=KFS*(sin(angle))

!CALL shepard_interp_2d(Nnumdisp,kxdftmat,kydftmat,dispdftmat,p_interp,1,kpxmat_interp,kpymat_interp,dispmat_interp)


!CALL READ_DATA(blab,Nnumdisp,kxdftmat,kydftmat,dispdftmat,dispmin,dispmax)

!kxdftmat(1:Nnumdisp)=BZconst*kxdftmat(1:Nnumdisp)/(abs(kpxmin))!(0.03855004306107383d+00) 
!kydftmat(1:Nnumdisp)=BZconst*kydftmat(1:Nnumdisp)/(abs(kpxmin))!(0.03855004306107383d+00)


p_interp=16.0d+00

kpxmat_interp(1)=KFS*(cos(angle))
kpymat_interp(1)=KFS*(sin(angle))

CALL shepard_interp_2d(Nnumdisp,kxdftmat,kydftmat,dispdftmat,p_interp,1,kpxmat_interp,kpymat_interp,dispmat_interp)


dispmat_ex(i)=dispmat_interp(1)


!dispmat_ex(i)=dispersion_theta(blab,KFS,angle,0.0d+00)

if(i>1)then

ddispmat_ex(i)=(dispmat_ex(i)-dispmat_ex(i-1))/(KFSmax/dimwork)

!write(3,*)kpxmat_interp(1)/BZconst,kpymat_interp(1)/BZconst,dispmat_ex(i),angle/(2.0d+00),KFS/BZconst,ddispmat_ex(i),nsign

endif

end do

!------------------------------------------------------------------
Do i=2,dimwork-1

if(i>1)then

KFS=KFSmax*i/dimwork

!write(2,*)KFS/BZconst,dispmat_ex(i),ddispmat_ex(i),ddispmat_ex(i)*ddispmat_ex(i+1)

if(ddispmat_ex(i)*ddispmat_ex(i+1)<0.0d+00)then

Do j=1,2

if(ddispmat_ex(i)*ddispmat_ex(i+1+j)<0.0d+00)then

ncrit=ncrit+1

endif  !changing sign with j add


if(ncrit>1)then

KFSext(nsign)=KFS
dispext(nsign)=dispmat_ex(i)

nsign=nsign+1

endif  !ncrit

enddo

else

ncrit=0

endif  !changing sign

endif !i>1


if(nsign>5)then
flag=0
else
flag=1
end if


if(nsign==1)then

CALL BZ_edge(angle,KFS)

KFSext(1)=KFS
e_ext=dispmat_ex(dimwork)
e_bound=dispmat_ex(dimwork)

elseif(nsign==2)then

KFSext(2)=KFS
e_ext=dispext(1)
e_bound=dispmat_ex(dimwork)

elseif(nsign>2)then

e_ext=dispext(1)
e_bound=dispext(nsign-1)

endif

end do !i


if(abs(e_ext-e_bound)>10.0d+00**(-10.0d+00) .and. i_min==0)then

angle_min=angle
i_min=1

elseif(abs(e_ext-e_bound)<10.0d+00**(-10.0d+00) .and. i_min==1 .and. i_max==0)then

angle_max=angle
i_max=1

end if

kpx=KFSext(1)*(cos(angle))
kpy=KFSext(1)*(sin(angle))

angle_mat(k)=angle/(2.0d+00*pi)
kred_mat(k)=KFSext(1)/BZconst
e_ext_mat(k)=e_ext
e_bound_mat(k)=e_bound


!write(2,*)angle/(2.0d+00*pi),KFSext(1)/BZconst,e_ext,e_bound,kpx/BZconst,kpy/BZconst!,ex_min_f,ex_max_f


end do 

!----------------------------------------------------------------------
end Subroutine Extermum_theta_bunch
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
DOUBLE PRECISION function dispersion_theta(blab,mom,theta,energy)

use commons

implicit none

common eF!,eF_min,eF_max

DOUBLE PRECISION,INTENT(IN)::mom,theta,energy
INTEGER, INTENT(IN)::blab
DOUBLE PRECISION::eF,kx,ky,eF_min,eF_max

INTEGER*4,parameter::N_interp=1
REAL*8,DIMENSION(1:N_interp)::kpxmat_interp,kpymat_interp,dispmat_interp
REAL*8::p_interp

INTEGER::openstatus,i,unitm

DOUBLE PRECISION,DIMENSION(1:Nnumdisp)::kxdftmat,kydftmat,dispdftmat
DOUBLE PRECISION::dispmin,dispmax,kpxmin,kpymin,kpx,kpy,dispcheck

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!------------------------------------------------------------------
!open(unit=41,file="dispersion_band_1.txt",status="old",action="read",position="rewind",IOSTAT=openstatus)
!open(unit=42,file="dispersion_band_2.txt",status="old",action="read",position="rewind",IOSTAT=openstatus)
!open(unit=43,file="dispersion_band_3.txt",status="old",action="read",position="rewind",IOSTAT=openstatus)
!open(unit=44,file="dispersion_band_4.txt",status="old",action="read",position="rewind",IOSTAT=openstatus)
!-----------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------

!unitm=40+blab

!print*,'unitm is',unitm

!Do i=1,Nnumdisp

!read(unitm,*)kpx,kpy,dispcheck

!kxdftmat(i)=kpx
!kydftmat(i)=kpy

!dispdftmat(i)=dispcheck*100.0d+00

!--------------------------------
! Kpxmin:
!--------------------------------
!if(i==1)then
!kpxmin=kxdftmat(i)
!endif

!if(kxdftmat(i)<kpxmin)then
!kpxmin=kxdftmat(i)
!endif
!--------------------------------

!end do

!print*,'kpxmin=',kpxmin


CALL READ_DATA(blab,Nnumdisp,kxdftmat,kydftmat,dispdftmat,dispmin,dispmax)


!kxdftmat(1:Nnumdisp)=BZconst*kxdftmat(1:Nnumdisp)/(abs(kpxmin))!(0.03855004306107383d+00) 
!kydftmat(1:Nnumdisp)=BZconst*kydftmat(1:Nnumdisp)/(abs(kpxmin))!(0.03855004306107383d+00)


p_interp=16.0d+00

kpxmat_interp(1)=mom*(cos(theta))
kpymat_interp(1)=mom*(sin(theta))

CALL shepard_interp_2d(Nnumdisp,kxdftmat,kydftmat,dispdftmat,p_interp,1,kpxmat_interp,kpymat_interp,dispmat_interp)

dispersion_theta=dispmat_interp(1)-energy

!----------------------------------------------------------------------
end function dispersion_theta
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! This subroutine generates patches for a given chemical potential
!-------------------------------
! Inputs:
!-------------------------------
! 1. eF (real) :chemical potential
!-------------------------------
! Outputs:
!-------------------------------
!1. thetapn (two dimensional real array of size (nband,npatch)) : angles of the patches 
!2. KFSpmatn (two dimensional real array of size (nband,npatch)): k values on the Fermi surfaces
!3. nrootspmatn (two dimensional integer array of size (nband,npatch)) : 1 if the Fermi surface has been found and 0 otherwise
!----------------------------------------------------------------------
!----------------------------------------------------------------------
Subroutine patch_generator(eF,thetapn,KFSpmatn,nrootspmatn)

use commons

implicit none

DOUBLE PRECISION,INTENT(IN)::eF                                     !chemical potential            
DOUBLE PRECISION,DIMENSION(1:nband,1:npatch),INTENT(OUT)::thetapn   !angles of the patches 
DOUBLE PRECISION,DIMENSION(1:nband,1:npatch),INTENT(OUT)::KFSpmatn  !k values on the Fermi surfaces
INTEGER,DIMENSION(1:nband,1:npatch),INTENT(OUT)::nrootspmatn        !1 if the Fermi surface has been found and 0 otherwise
!----------------------------------------
! Dummy variables:
!----------------------------------------
INTEGER::n_angles
DOUBLE PRECISION,DIMENSION(1:nband,1:npatch)::angle_mat
INTEGER,DIMENSION(1:nband)::nFS_max_mat

INTEGER::l1,i_min,i_max,icounter,unitm,flag_ex,i,k,j,nsign,nFS &
         ,nroots,ns,nFS_max,BZid,nroots_real
DOUBLE PRECISION::angle,KFS,ex_min,ex_max,e_min,e_max,energy1,energy2&
                 ,angle_min,angle_max,kpx,kpy,KFSmax

DOUBLE PRECISION,DIMENSION(1:nband,1:npatch)::thetap,KFSpmat,kmin_mat,kmax_mat
INTEGER,DIMENSION(1:nband,1:npatch)::nrootspmat

DOUBLE PRECISION,DIMENSION(1:10)::KFSext_ex
!----------------------------------------------------------------------

if(disp_flag=='num')then

!----------------------------------------------------------------------

Do l1=1,nband/2

CALL Bands_scan('test',num_theta_ired,l1,angle_min,angle_max,e_min,e_max,nFS_max)

nFS_max_mat((2*l1-1):(2*l1))=nFS_max

!print*,'band',l1,'nFS_max=',nFS_max,'angle_min=',angle_min,'angle_max=',angle_max,'e_min=',e_min,'e_max=',e_max

!----------------------------------------------------------------------

if(nFS_max==1 .or. eF<e_min .or. eF>e_max)then

nFS=1
n_angles=npatch

else

nFS=2
n_angles=(npatch+12)/2

endif

!----------------------------------------------------------------------

!print*,'n_angles=',n_angles,'nFS=',nFS

k=0

!--------------------------------

Do j=1,n_angles/6

if(nFS==1)then
angle=(2.0d+00*pi)*(j-0.5d+00)/n_angles

CALL BZ_edge(angle,KFSmax)

thetap(l1,j)=angle

kmin_mat(l1,j)=0.0d+00

!print*,'l1=',l1,j,kmin_mat(l1,j)

kmax_mat(l1,j)=KFSmax

else

if(n_angles/6==1)then

angle=angle_min+(angle_max-angle_min)/(2.0d+00)

else

if(j==1)then
angle=angle_min/2.0d+00
ns=1
elseif(j==n_angles/6)then
angle=angle_max+(((pi/3.0d+00)-angle_max)/2.0d+00)
ns=1
else
angle=angle_min+(angle_max-angle_min)*(j-1.5d+00)/((n_angles/6)-2)
ns=2
endif

endif

CALL BZ_edge(angle,KFSmax)

Do i=1,ns

k=k+1
thetap(l1,k)=angle

if(ns==1)then
kmin_mat(l1,k)=0.0d+00
kmax_mat(l1,k)=KFSmax
else

if(i==1)then
CALL Extermum_theta(l1,angle,50,nsign,KFSext_ex,flag_ex,ex_max,ex_min)

kmin_mat(l1,k)=0.0d+00
kmax_mat(l1,k)=KFSext_ex(1)

else

kmin_mat(l1,k)=KFSext_ex(1)+10.0d+00**(-5.0d+00)
kmax_mat(l1,k)=KFSext_ex(2)

endif

endif

end do  !i

endif

angle_mat((2*l1-1):2*l1,j)=angle

end do !j

!print*,'finding the roots:'

nroots_real=0

Do k=1,npatch/6

angle=thetap(l1,k)

CALL roots_dispersion_theta2(2*l1,angle,eF,kmin_mat(l1,k),kmax_mat(l1,k),nroots,KFS)

!print*,'k=',k,'nroots=',nroots,kmin_mat(l1,k),kmax_mat(l1,k)

nroots_real=nroots_real+nroots


if(nroots==0)then
KFS=kmax_mat(l1,k)
nroots=1
endif


KFSpmat(l1,k)=KFS
nrootspmat(l1,k)=nroots


Do j=1,5

thetap(l1,k+j*(npatch/6))=angle+j*(pi/3.0d+00)
KFSpmat(l1,k+j*(npatch/6))=KFS
nrootspmat(l1,k+j*(npatch/6))=nroots

end do !j

end do !k

if(nroots_real>0)then
print*, 'bands', 2*l1-1,'and', 2*l1, 'are relevant!'
endif


Do k=1,npatch

angle=thetap(l1,k)
KFS=KFSpmat(l1,k)
nroots=nrootspmat(l1,k)

thetapn((2*l1-1):2*l1,k)=angle
KFSpmatn((2*l1-1):2*l1,k)=KFS
nrootspmatn((2*l1-1):2*l1,k)=nroots

!print*,'reconstruction check',l1,2*l1-1,2*l1

kpx=(KFS)*cos(angle)
kpy=(KFS)*sin(angle)

!print*,kpx/BZconst,kpy/BZconst,nroots

!write(2,*)kpx/BZconst,kpy/BZconst

end do
!----------------------------------------------------------------------
end do
!----------------------------------------------------------------------
elseif(disp_flag=='ana')then

Do l1=1,nband

nroots_real=0

Do k=1,npatch

angle=(2.0d+00*pi/npatch)*(k-0.5d+00)

CALL roots_dispersion(l1,angle,eF,nroots,KFS)

nroots_real=nroots_real+nroots

if(nroots<1)then
CALL BZ_edge(angle,KFS)
nroots=1
endif

CALL Hexagonal(kpx,kpy,BZid)

if(BZid==0 .and. nroots==1)then

CALL Folding(kpx,kpy,kpx,kpy)

KFS=sqrt((kpx**2)+(kpy**2))

angle=atan2(kpy,kpx)
angle=(-SIGN(1.0d+00,angle)+1.0d+00)*(pi)+angle

endif

thetapn(l1,k)=angle
KFSpmatn(l1,k)=KFS
nrootspmatn(l1,k)=nroots

if(l1==1 .or. l1==5)then
print*,'band=',l1,'k=',k,'nroots=',nroots,KFS/BZconst
endif

end do

if(nroots_real>0)then
print*,'band',l1,'is relevant'
endif

end do

end if
!----------------------------------------------------------------------
!----------------------------------------------------------------------
end Subroutine patch_generator
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! This subroutines gives the following information about band indices
!  that are independant of the choice of chemical potential
!--------------------------------
! Inputs:
!--------------------------------
! 1. flag (character) : 'test' or 'real'
! 2. N (integer) : number of angles in [0,pi/3], to find extermum of dipersion
! 3. l1 (integer) : band index
!--------------------------------
! Outputs:
!--------------------------------
! 1. angle_min (real) : minimum angle in [0,pi/3] to have at max two Fermi surfaces
! 2. angle_max (real) : maximum angle in [0,pi/3] to have at max two Fermi surfaces
! 3. e_min (real) : minimum chemical potential to have at max two Fermi surfaces
! 4. e_max (real) : maximum chemical potential to have at max two Fermi surfaces
! 5. nFS_max (integer) : maximum number of Fermi Surfaces
!----------------------------------------------------------------------
Subroutine Bands_scan(flag,N,l1,angle_min,angle_max,e_min,e_max,nFS_max)

use commons

implicit none

CHARACTER*4,INTENT(IN)::flag              !'test' or 'real'                           
INTEGER,INTENT(IN)::N                     !number of angles in [0,pi/3], to find extermum of dipersion
INTEGER,INTENT(IN)::l1                    !band index 
DOUBLE PRECISION,INTENT(OUT)::angle_min   !minimum angle in [0,pi/3] to have at max two Fermi surfaces
DOUBLE PRECISION,INTENT(OUT)::angle_max   !maximum angle in [0,pi/3] to have at max two Fermi surfaces
DOUBLE PRECISION,INTENT(OUT)::e_min       !minimum chemical potential to have at max two Fermi surfaces
DOUBLE PRECISION,INTENT(OUT)::e_max       !maximum chemical potential to have at max two Fermi surfaces
INTEGER,INTENT(OUT)::nFS_max

!-------------------------------------
! Dummy variables:
!-------------------------------------
INTEGER::i,unitm,i_min,i_max,icounter,openstatus
DOUBLE PRECISION::angle,KFS,energy1,energy2,ex_min,ex_max

INTEGER::dimwork_e,flag_e
DOUBLE PRECISION,DIMENSION(1:N)::angle_mat_e,kred_mat_e,e_ext_mat_e,e_bound_mat_e

!----------------------------------------------------------------------
if(flag=='test')then
!-----------------------------------------------------------------------------------------------------
open(unit=71,file="band1_kred.txt",status="old",action="read",position="rewind",IOSTAT=openstatus)
open(unit=72,file="band2_kred.txt",status="old",action="read",position="rewind",IOSTAT=openstatus)
open(unit=73,file="band3_kred.txt",status="old",action="read",position="rewind",IOSTAT=openstatus)
open(unit=74,file="band4_kred.txt",status="old",action="read",position="rewind",IOSTAT=openstatus)
!-----------------------------------------------------------------------------------------------------
elseif(flag=='real')then
!-----------------------------------------------------------------------------------------------------
open(unit=71,file="band1_kred.txt",status="new",action="write",position="rewind",IOSTAT=openstatus)
open(unit=72,file="band2_kred.txt",status="new",action="write",position="rewind",IOSTAT=openstatus)
open(unit=73,file="band3_kred.txt",status="new",action="write",position="rewind",IOSTAT=openstatus)
open(unit=74,file="band4_kred.txt",status="new",action="write",position="rewind",IOSTAT=openstatus)
!-----------------------------------------------------------------------------------------------------
endif
!----------------------------------------------------------------------

print*,'I am in Bands_scan: band_index=',l1

i_min=0
i_max=0
icounter=0

e_min=0.0d+00
e_max=0.0d+00

!--------------------------------
unitm=70+l1

if(unitm==71)then
open(unit=unitm,file='band1_kred.txt')
elseif(unitm==72)then
open(unit=unitm,file='band2_kred.txt')
elseif(unitm==73)then
open(unit=unitm,file='band3_kred.txt')
elseif(unitm==74)then
open(unit=unitm,file='band4_kred.txt')
endif


if(flag=='real')then
CALL Extermum_theta_bunch(l1,num_kvec_thf,N,flag_e,angle_mat_e,kred_mat_e,e_ext_mat_e,e_bound_mat_e)
endif


Do i=1,N

!-----------------------------------------------------------------------
if(flag=='test')then

read(unitm,*)angle,KFS,ex_max,ex_min

elseif(flag=='real')then

angle=angle_mat_e(i)
KFS=kred_mat_e(i)
ex_max=e_ext_mat_e(i)
ex_min=e_bound_mat_e(i)

write(unitm,*)angle,KFS,ex_max,ex_min

endif

!-----------------------------------------------------------------------

energy1=min(ex_min,ex_max)
energy2=max(ex_min,ex_max)

ex_min=energy1
ex_max=energy2


if(abs(ex_max-ex_min)>10.0d+00**(-10.0d+00))then

icounter=icounter+1


!write(3,*)angle,ex_max,ex_min,abs(ex_max-ex_min)


if(icounter==1)then
e_min=ex_min
e_max=ex_max
endif

if(ex_min<e_min)then
e_min=ex_min
endif

if(ex_max>e_max)then
e_max=ex_max
endif

endif

if(abs(ex_max-ex_min)>10.0d+00**(-10.0d+00) .and. i_min==0)then
angle_min=angle*2.0d+00*pi
i_min=1

elseif(abs(ex_max-ex_min)<10.0d+00**(-10.0d+00) .and. i_min==1 .and. i_max==0 .and. angle>1.0d+00/12.0d+00)then
angle_max=angle*2.0d+00*pi
i_max=1
endif

end do !i

CLOSE(unit=unitm)

if(icounter>0)then
nFS_max=2
else
nFS_max=1
endif

!print*,'check:',l1,nFS_max,e_min,e_max,icounter

!----------------------------------------------------------------------
end Subroutine Bands_scan
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! This subroutine finds for a given angle the k value for which there 
! is a unique Fermi Surface.
! Inputs :
! 1. flag (character): 'flag=test' to read from file and 'flag=real' to find from dispersion
! 2. N (integer) : number of angles in case flag=test
! 3. l1 (integer) : band-index can take 1,2,3,4
! 4. angle (real) : angle in [0,pi/3]
! Outputs:
! 1.kred (real) : the k value for which we have a unique Fermi surface
!----------------------------------------------------------------------
!----------------------------------------------------------------------
Subroutine kred_generator(flag,N,l1,angle,kred)

use commons

implicit none

CHARACTER*4,INTENT(IN)::flag           !'flag=test' to read from file and 'flag=real' to find from dispersion
INTEGER,INTENT(IN)::N                  !number of angles in case flag=test
INTEGER,INTENT(IN)::l1                 !band-index can take 1,2,3,4
DOUBLE PRECISION,INTENT(IN)::angle     !angle in [0,pi/3]
DOUBLE PRECISION,INTENT(OUT)::kred     !the k value for which we have a unique Fermi surface

!----------------------------------------------------------------------
! Dummy variables
!----------------------------------------------------------------------

DOUBLE PRECISION::theta,KFS
INTEGER::i,istar,unitm
DOUBLE PRECISION,DIMENSION(1:N)::theta_mat,kred_mat
!----------------------------------------------------------------------

!print*,'I am in kred_generator!','l1=',l1,'angle/(2*pi)',angle/(2.0d+00*pi)

if(angle>pi/3.0d+00)then
print*,'Give me angles in the range [0,pi/3]'
endif


unitm=70+l1

if(unitm==71)then
open(unit=unitm,file='band1_kred.txt')
elseif(unitm==72)then
open(unit=unitm,file='band2_kred.txt')
elseif(unitm==73)then
open(unit=unitm,file='band3_kred.txt')
elseif(unitm==74)then
open(unit=unitm,file='band4_kred.txt')
endif

Do i=1,N

if(flag=='test')then
read(unitm,*)theta,KFS
endif

theta_mat(i)=theta*2.0d+00*pi
Kred_mat(i)=KFS*BZconst

end do

CLOSE(unit=unitm)

Do i=1,N

if(angle .LE. theta_mat(1))then
istar=1
elseif(angle .GE. theta_mat(N))then
istar=N
elseif(i>1 .and. angle > theta_mat(i-1) .and. angle .LE. theta_mat(i))then
istar=i
endif

end do




if(istar==1)then

kred=kred_mat(istar)+(angle-theta_mat(istar))&
                    *((kred_mat(istar+1)-kred_mat(istar))/(theta_mat(istar+1)-theta_mat(istar)))

else

kred=kred_mat(istar)+(angle-theta_mat(istar))&
                    *((kred_mat(istar)-kred_mat(istar-1))/(theta_mat(istar)-theta_mat(istar-1)))


endif


!print*,'from kred_generator:','istar=',istar,kred_mat(istar)/BZconst,kred/BZconst


end Subroutine kred_generator
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! This subroutine provides information about the bands
!----------------------------------
! Inputs:
!----------------------------------
!1. N (integer): number of angles to find extermums of dispersion 
!----------------------------------
! Outputs:
!----------------------------------
!1.angle_min_mat (one dimensional real array of size nband) : minimum angle in [0,pi/3] to have at max two Fermi surfaces
!2.angle_max_mat (one dimensional real of size nband) : maximum angle in [0,pi/3] to have at max two Fermi surfaces
!3.e_min_mat (one dimensional real of size nband) : minimum chemical potential to have at max two Fermi surfaces
!4.e_max_mat (one dimensional real of size nband) : maximum chemical potential to have at max two Fermi surfaces
!5.nFS_max_mat (one dimensional integer array of size nband) : maximum number of Fermi Surfaces
!----------------------------------------------------------------------
!----------------------------------------------------------------------
Subroutine Bands_Info(flag,N,angle_min_mat,angle_max_mat,e_min_mat,e_max_mat,nFS_max_mat)

use commons

implicit none

CHARACTER*4::flag                                               !it can be real or test to call band_scan subroutine
INTEGER,INTENT(IN)::N                                           !number of angles to find extermums of dispersion 
DOUBLE PRECISION,DIMENSION(1:nband),INTENT(OUT)::angle_min_mat  !minimum angle in [0,pi/3] to have at max two Fermi surfaces
DOUBLE PRECISION,DIMENSION(1:nband),INTENT(OUT)::angle_max_mat  !maximum angle in [0,pi/3] to have at max two Fermi surfaces
DOUBLE PRECISION,DIMENSION(1:nband),INTENT(OUT)::e_min_mat      !minimum chemical potential to have at max two Fermi surfaces
DOUBLE PRECISION,DIMENSION(1:nband),INTENT(OUT)::e_max_mat      !maximum chemical potential to have at max two Fermi surfaces
INTEGER,DIMENSION(1:nband),INTENT(OUT)::nFS_max_mat             !maximum number of Fermi Surfaces

!----------------------------------
! Dummy variables:
!----------------------------------
INTEGER::l1,nFS_max
DOUBLE PRECISION::angle_min,angle_max,e_min,e_max
!----------------------------------

print*,'I am in Bands_Info'


Do l1=1,nband/2

CALL Bands_scan(flag,N,l1,angle_min,angle_max,e_min,e_max,nFS_max)

angle_min_mat((2*l1-1):(2*l1))=angle_min
angle_max_mat((2*l1-1):(2*l1))=angle_max
e_min_mat((2*l1-1):(2*l1))=e_min
e_max_mat((2*l1-1):(2*l1))=e_max
nFS_max_mat((2*l1-1):(2*l1))=nFS_max

!print*,l1,(2*l1-1),(2*l1),nFS_max

end do
!----------------------------------------------------------------------

end Subroutine Bands_Info
!----------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------
! This subroutine helps to introduce a patching scheme for a given chemical potential
!-----------------------
!Inputs:
!-----------------------
!1. N (integer) : Number of angles for finding the chemical potential that has a unique Fermi surface
!2. Nk (integer) : Number of k points for a given angle to find extermums of dispersion
!-----------------------
!Outputs:
!-----------------------
!1. n_angles_mat (2 dimensional integer array of size (nband,2)) : number of angles for various bands for 
!   a given npatch for the case of one or two Fermi surfaces
!2. angles_mat (3 dimensional real array of size (nband,2,n_angles_mat)) : the angles for variou patches
!   for the case of one or two Fermi surfaces
!3. nFS_mat (3 dimensional integer array of size (nband,2,n_angles_mat)) : the number of Fermi surfaces for a given angle
!4. kmin_mat (4 dimensional real array of size (nband,2,n_angles_mat,2)) : minimum value of vec(k) such that dispersion is monotonic
!5. kmax_mat (4 dimensional real array of size (nband,2,n_angles_mat,2)) : maximum value of vec(k) such that dispersion is monotonic
!6. kred_mat (2 dimensional real array of size (nband,patch)) :magnitute of vec(k) for which there is a unique Fermi surface 
!----------------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------------
Subroutine patch_assist(N,Nk,n_angles_mat,angles_mat,nFS_mat,kmin_mat,kmax_mat,kred_mat)

use commons

implicit none

INTEGER,INTENT(IN)::N                                                             
INTEGER,INTENT(IN)::Nk
INTEGER,DIMENSION(1:nband,1:2),INTENT(OUT)::n_angles_mat                                !number of angles for various bands for a given npatch for t
DOUBLE PRECISION,DIMENSION(1:nband,1:2,1:npatch),INTENT(OUT)::angles_mat                !the angles for variou patches
INTEGER,DIMENSION(1:nband,1:2,1:npatch),INTENT(OUT)::nFS_mat                            !the number of Fermi surfaces for a given angle
DOUBLE PRECISION,DIMENSION(1:nband,1:2,1:npatch,1:2),INTENT(OUT)::kmin_mat              !minimum value of vec(k) such that dispersion is monotonic
DOUBLE PRECISION,DIMENSION(1:nband,1:2,1:npatch,1:2),INTENT(OUT)::kmax_mat              !maximum value of vec(k) such that dispersion is monotonic
DOUBLE PRECISION,DIMENSION(1:nband,1:npatch),INTENT(OUT)::kred_mat                      !magnitute of vec(k) for which there is a unique Fermi surface            
!--------------------------------------------
! Dummy variables:
!--------------------------------------------
INTEGER::k,j,n_angles,l1,nsign,flag_ex,l2
DOUBLE PRECISION::angle,KFSmax,ex_max,ex_min,kred

DOUBLE PRECISION,DIMENSION(1:nband)::angle_min_mat,angle_max_mat,e_min_mat,e_max_mat
INTEGER,DIMENSION(1:nband)::nFS_max_mat

DOUBLE PRECISION,DIMENSION(1:10)::KFSext_ex
!----------------------------------------------------------------------
!----------------------------------------------------------------------

print*,'I am in patch_assist'


CALL Bands_Info('test',N,angle_min_mat,angle_max_mat,e_min_mat,e_max_mat,nFS_max_mat)

!----------------------------------------------------------------------
Do l1=1,nband
!----------------------------------------------------------------------
print*,'l1=',l1,'nFS=',nFS_max_mat(l1)
!---------------------------------

n_angles=npatch

n_angles_mat(l1,1)=n_angles

Do j=1,n_angles/6

angle=(2.0d+00*pi)*(j-0.5d+00)/n_angles

angles_mat(l1,1,j)=angle

CALL BZ_edge(angle,KFSmax)

kmin_mat(l1,1,j,1)=0.0d+00

kmax_mat(l1,1,j,1)=KFSmax

kred_mat(l1,j)=KFSmax

nFS_mat(l1,1,j)=1

end do  !j
!---------------------------------
!---------------------------------

if(nFS_max_mat(l1)==2)then

n_angles=(npatch+12)/2

n_angles_mat(l1,2)=n_angles

!---------------------------------
Do j=1,n_angles/6

if(n_angles/6==1)then
angle=angle_min_mat(l1)+(angle_max_mat(l1)-angle_min_mat(l1))/(2.0d+00)
nFS_mat(l1,2,j)=2
else

if(j==1)then
angle=angle_min_mat(l1)/2.0d+00
nFS_mat(l1,2,j)=1
elseif(j==n_angles/6)then
angle=angle_max_mat(l1)+(((pi/3.0d+00)-angle_max_mat(l1))/2.0d+00)
nFS_mat(l1,2,j)=1
else
angle=angle_min_mat(l1)+(angle_max_mat(l1)-angle_min_mat(l1))*(j-1.5d+00)/((n_angles/6)-2)
nFS_mat(l1,2,j)=2
endif

endif

angles_mat(l1,2,j)=angle

l2=1+int((l1-0.1d+00)/2.0d+00)

!print*,'l1 is',l1,'l2 is',l2,'angle=',angle,angle_min_mat(l1),angle_max_mat(l1),nFS_max_mat(l1)

CALL Extermum_theta(l2,angle,Nk,nsign,KFSext_ex,flag_ex,ex_max,ex_min)

kmin_mat(l1,2,j,1)=0.0d+00
kmax_mat(l1,2,j,1)=KFSext_ex(1)

kmin_mat(l1,2,j,2)=KFSext_ex(1)+10.0d+00**(-5.0d+00)
kmax_mat(l1,2,j,2)=KFSext_ex(2)


CALL kred_generator('test',N,l2,angle,kred)

kred_mat(l1,j)=kred

end do  !j
!---------------------------------
endif  !nFS_max_mat(l1)==2
!---------------------------------
!----------------------------------------------------------------------
end do
!----------------------------------------------------------------------
end Subroutine patch_assist
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! This subroutine calculates the patch index for 
!      vec(p4)=vex(p1)+vec(p2)-vec(p3): 
!---------------------------------------------------------
!Inputs:
!---------------------------------------------------------
!1. N (integer) : number of patches
!2. thetap_mat (two dimensional array of size (nband,npatch)) : angles of patches
!3. KFSp_mat (two dimensional array of size (nband,npatch)) : corresponding magnitute of vectors on the Fermi surface
!4. energy (real) :chemical potential
!5. bands_info_mat (2 dimensional array of size (5,nband))
   !bands_info_mat(1,1:nband)=1.0d+00*nFS_max_mat(1:nband)   !maximum number of Fermi Surfaces
   !bands_info_mat(2,1:nband)=angle_min_mat(1:nband)         !minimum angle in [0,pi/3] for which we might get two Fermi Surfaces
   !bands_info_mat(3,1:nband)=angle_max_mat(1:nband)         !maximum angle in [0,pi/3] for which we might get two Fermi Surfaces
   !bands_info_mat(4,1:nband)=e_min_mat(1:nband)             !minimum dispersion for which we have two Fermi surfaces
   !bands_info_mat(5,1:nband)=e_max_mat(1:nband)             !maximum dispersion for which we have two Fermi surfaces 
!---------------------------------------------------------
!Outputs:
!---------------------------------------------------------
!1. p4vec (7 dimensional real array with size
!         (nband,nband,nband,npatch,npatch,npatch,2))        !kx and ky value of vec(p4)
!2. p4int  (7 dimensional integer array of size 
!          (nband,nband,nband,nband,npatch,npatch,npatch))   !integer associated to vec(p4)
!---------------------------------------------------------
!----------------------------------------------------------------------
subroutine p4_generator(N,thetap_mat,KFSp_mat,energy,bands_info_mat,p4vec,p4int)
!----------------------------------------------------------------------

use commons

IMPLICIT NONE

INTEGER,INTENT(IN)::N                                                                   !number of patches
DOUBLE PRECISION,DIMENSION(1:nband,1:N),INTENT(IN)::thetap_mat                          !chosen angles associated to patches
DOUBLE PRECISION,DIMENSION(1:nband,1:N),INTENT(IN)::KFSp_mat                            !corresponding magnitute of vectors on the Fermi surface

DOUBLE PRECISION,INTENT(IN)::energy                                                     !chemical potential
DOUBLE PRECISION,DIMENSION(1:5,1:nband),INTENT(IN)::bands_info_mat                      !Information about bands (see above)

DOUBLE PRECISION,DIMENSION(1:nband,1:nband,1:nband,1:N,1:N,1:N,1:2),INTENT(OUT)::p4vec  !the kx and ky componant of the vec(p4)
INTEGER,DIMENSION(1:nband,1:nband,1:nband,1:nband,1:N,1:N,1:N),INTENT(OUT)::p4int       !patch index of vec(p4)

!----------------------------------------------------------------------
! Dummy variables for the evaluation:
!----------------------------------------------------------------------


INTEGER::n1,n2,n3,n4,p1,p2,p3,BZid,k,nFS
DOUBLE PRECISION::angle1,angle2,angle3,KFS1,KFS2,KFS3,kpx,kpy &
                 ,kpxr,kpyr,anglek,anglered,KFSred,KFS

INTEGER,DIMENSION(1:nband)::nFS_mat
DOUBLE PRECISION,DIMENSION(1:nband)::KFSred_mat,angle_min_mat,angle_max_mat

CHARACTER*3::patch_flag
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Vec-add on the Fermi Surface!
!----------------------------------------------------------------------
Do n1=1,nband

Do p1=1,npatch

angle1=thetap_mat(n1,p1)
KFS1=KFSp_mat(n1,p1)

Do n2=1,nband

Do p2=1,npatch

angle2=thetap_mat(n2,p2)
KFS2=KFSp_mat(n2,p2)

Do n3=1,nband

Do p3=1,npatch

angle3=thetap_mat(n3,p3)
KFS3=KFSp_mat(n3,p3)

kpx=KFS1*(cos(angle1))+KFS2*(cos(angle2))-KFS3*(cos(angle3))
kpy=KFS1*(sin(angle1))+KFS2*(sin(angle2))-KFS3*(sin(angle3))


CALL Hexagonal(kpx,kpy,BZid)


if(BZid==0)then

CALL Folding(kpx,kpy,kpxr,kpyr)

else

kpxr=kpx
kpyr=kpy

endif

CALL Hexagonal(kpxr,kpyr,BZid)

if(BZid==0)then
print*,'careful p4 not on the FS!',kpxr/BZconst,kpyr/BZconst,p1,p2,p3,n1,n2,n3
end if

anglek=atan2(kpyr,kpxr)
anglek=(-SIGN(1.0d+00,anglek)+1.0d+00)*(pi)+anglek


p4vec(n1,n2,n3,p1,p2,p3,1)=kpxr
p4vec(n1,n2,n3,p1,p2,p3,2)=kpyr


if(sqrt((kpxr)**2+(kpyr)**2).LE. 10.0d+00**(-4.0d+00))then

p4int(n1,n2,n3,1:nband,p1,p2,p3)=0

else

!------------------------------------------------------------------
! Finding the corresponding patch index:
!------------------------------------------------------------------
Do n4=1,nband

angle_min_mat(1:nband)=bands_info_mat(2,1:nband)
angle_max_mat(1:nband)=bands_info_mat(3,1:nband)

if(bands_info_mat(1,n4)<2.0d+00 .or. energy<bands_info_mat(4,n4) .or. energy>bands_info_mat(5,n4))then

patch_flag='sim'

else

patch_flag='mul'

endif


CALL patch_lab_Ge(patch_flag,anglek,sqrt(((kpxr)**2)+((kpyr)**2)),n4,angle_min_mat,angle_max_mat,k,anglered)
p4int(n1,n2,n3,n4,p1,p2,p3)=k

if(k<1 .or. k>npatch)then
print*,'patch-lab-Ge-not-working-properly-here is p4'
endif

end do !n4
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
end subroutine p4_generator
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
! This subroutine initialize the vertex function
!---------------------------------
! Inputs: 
!---------------------------------
!1. flag (character) : at the moment it can be only 'local'
!2. thetap (two dimensional real array (nband,npatch)) : angles for various bands and patches 
!3. KFSpmat (two dimensional real array (nband,npatch)) : corresponding k value on the Fermi surface
!---------------------------------
! Outputs: 
!---------------------------------------------------------------------------------
!2. Yvertexcomplex (one dimensional complex array of size (nband**4)*(npatch**3))
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
Subroutine Initial_vertex(flag,thetap,KFSpmat,Yvertexcomplex)

use commons

IMPLICIT NONE

CHARACTER*5,INTENT(IN)::flag                          ! input: flag='local' corresponds to local Coulomb interaction
DOUBLE PRECISION,DIMENSION(1:nband,1:npatch)::thetap  ! input: real the angles corresponding to patches
DOUBLE PRECISION,DIMENSION(1:nband,1:npatch)::KFSpmat ! input: real the magnitute of vector on the Fermi surface 
COMPLEX,DIMENSION(1:Ntot),INTENT(OUT)::Yvertexcomplex ! output: Intial vertex 

!-----------------------------------------------------------------------
! Rest of dummy variables:
!-----------------------------------------------------------------------


INTEGER::p1,p2,p3,n1orb,n2orb,n3orb,n4orb,n1spin,n2spin,n3spin,n4spin &
         ,n1type,n2type,n3type,n4type,n1,n2,n3,n4,i,iprime &
         ,n1prime,n2prime,n3prime,n4prime,k,n1p,n2p,n3p,n4p,imat_fun,imat

DOUBLE PRECISION::Yverdummy,YverdummyV,Yverdummyu,kpx,kpy,kdelta,eF,eF_min,eF_max,maxvertex

COMPLEX*16,DIMENSION(1:2,1:2)::Umatp4

DOUBLE PRECISION,DIMENSION(1:Ntot)::YvertexABinitial
COMPLEX*16,DIMENSION(1:2,1:2,1:npatch)::Umattype

INTEGER,DIMENSION(1:ntype,1:norbital,1:nspin)::nbandmat
INTEGER,DIMENSION(1:nband)::nspinmat,norbmat,ntypemat

DOUBLE PRECISION,DIMENSION(1:nband,1:nband,1:nband,1:npatch,1:npatch,1:npatch,1:2)::p4vec
INTEGER,DIMENSION(1:nband,1:nband,1:nband,1:nband,1:npatch,1:npatch,1:npatch)::p4int

DOUBLE PRECISION,DIMENSION(1:nband,1:npatch)::thetap_mat,KFSp_mat

DOUBLE PRECISION,DIMENSION(1:Ndis)::kxsamp,kysamp
DOUBLE PRECISION,DIMENSION(1:nband,1:Ndis,1:Ndis)::dissamp

INTEGER,DIMENSION(1:Ndis)::kyintmaxmat,kyintminmat

!------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
Do k=1,npatch

Umattype(1,1,k)=(1.0d+00/sqrt(2.0d+00))*cmplx(-1.0d+00,0.0d+00)
Umattype(1,2,k)=(1.0d+00/sqrt(2.0d+00))*cmplx(1.0d+00,0.0d+00)
Umattype(2,1,k)=Umattype(1,2,k)
Umattype(2,2,k)=Umattype(1,2,k)

end do !k
!----------------------------------------------------------------------
!----------------------------------------------------------------------
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

end do !n1spin
end do !n1orb
end do !n1type

!-----------------------------------------------------------------------

Yvertexcomplex(1:Ntot)=cmplx(0.0d+00,0.0d+00)

!-----------------------------------------------------------------------
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


i=imat(n1,n2,n3,n4,p1,p2,p3)


if(p1==1 .and. p2==1 .and. p3==1)then


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

iprime=imat(n1prime,n2prime,n3prime,n4prime,p1,p2,p3)



YvertexABinitial(iprime)=((kdelta(n1p,n2p))*(kdelta(n3p,n4p))*(kdelta(n1p,n3p)))*&
            ((u)*(kdelta(n1orb,n2orb)*kdelta(n3orb,n4orb)*kdelta(n1orb,n3orb))&
              *(((kdelta(n1spin,n3spin))*(kdelta(n2spin,n4spin)))-((kdelta(n1spin,n4spin))*(kdelta(n2spin,n3spin))))&
            +(u)*(kdelta(n1spin,n2spin)*kdelta(n3spin,n4spin)*kdelta(n1spin,n3spin))&                          !(uprime-Jcoup)
              *(((kdelta(n1orb,n3orb))*(kdelta(n2orb,n4orb)))-((kdelta(n1orb,n4orb))*(kdelta(n2orb,n3orb))))&
            +(u)*(1.0d+00-kdelta(n2orb,n1orb))*(1.0d+00-kdelta(n1spin,n2spin))&
                     *((kdelta(n1spin,n3spin))*(kdelta(n2spin,n4spin))*(kdelta(n1orb,n3orb))*(kdelta(n2orb,n4orb))&
                      -(kdelta(n2spin,n3spin))*(kdelta(n1spin,n4spin))*(kdelta(n2orb,n3orb))*(kdelta(n1orb,n4orb))))


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
else

iprime=imat(n1,n2,n3,n4,1,1,1)

Yvertexcomplex(i)=Yvertexcomplex(iprime)

endif !(p1=p2=p3=1)
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
end subroutine initial_vertex
!------------------------------------------------------------------
!------------------------------------------------------------------
! This subroutine converts the 7-dimensional array 
! (n1,n2,n3,n4,p1,p2,p3) into 1 dimensional
!------------------------------
! Inputs:
!------------------------------
!1. nb (integer) :number of bands (8)
!2. np (integer) :number of patches
!------------------------------
! Output:
!------------------------------
!1. imat (seven dimensional integer array of size (nb**4)*(np**3))
!2.n1mat (one dimensional integer array of size (nb**4)*(np**3)) : gives the first band index corresponding to 7 dimensional array
!3.n2mat (one dimensional integer array of size (nb**4)*(np**3)) : gives the second band index corresponding to 7 dimensional array
!4.n3mat (one dimensional integer array of size (nb**4)*(np**3)) : gives the third band index corresponding to 7 dimensional array
!5.n4mat (one dimensional integer array of size (nb**4)*(np**3)) : gives the fourth band index corresponding to 7 dimensional array
!------------------------------------------------------------------
subroutine imat_generator(nb,np,imat,n1mat,n2mat,n3mat,n4mat)

use commons

IMPLICIT NONE

INTEGER,INTENT(IN)::nb                                                   ! number of bands
INTEGER,INTENT(IN)::np                                                   ! number of patches
INTEGER,DIMENSION(1:nb,1:nb,1:nb,1:nb,1:np,1:np,1:np),INTENT(OUT)::imat  ! 7 dimensional output array
INTEGER,DIMENSION(1:Ntot),INTENT(OUT)::n1mat                             ! 1 dimensional array of first band index
INTEGER,DIMENSION(1:Ntot),INTENT(OUT)::n2mat                             ! 1 dimensional array of second band index
INTEGER,DIMENSION(1:Ntot),INTENT(OUT)::n3mat                             ! 1 dimensional array of third band index
INTEGER,DIMENSION(1:Ntot),INTENT(OUT)::n4mat                             ! 1 dimensional array of fourth band index

INTEGER::i,n1,n2,n3,n4,p1,p2,p3

i=0

Do n1=1,nb
Do n2=1,nb
Do n3=1,nb
Do n4=1,nb

Do p1=1,np
Do p2=1,np
Do p3=1,np


i=i+1


imat(n1,n2,n3,n4,p1,p2,p3)=i

n1mat(i)=n1
n2mat(i)=n2
n3mat(i)=n3
n4mat(i)=n4

end do
end do
end do

end do
end do
end do
end do

!------------------------------------------------------------------
end subroutine imat_generator
!------------------------------------------------------------------
!------------------------------------------------------------------
! This subroutine is just plotting the Brillioun Zone
! for Hexagonal lattice
!------------------------
! Inputs:
!------------------------
! 1.N (integer) : number of angles
!------------------------------------------------------------------
!------------------------------------------------------------------
Subroutine Plot_BZ(N)

use commons

IMPLICIT NONE

INTEGER,INTENT(IN)::N     !number of angles in interval [0,2*pi]

!------------------------------------------
! Dummy variables:
!------------------------------------------

DOUBLE PRECISION::angle,KFSmax,kpx,kpy  
INTEGER::k,openstatus

!-----------------------------------------------------------------------------------------------------
open(unit=5,file="Sep121BZ.txt",status="new",action="write",position="rewind",IOSTAT=openstatus)
!-----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
Do k=1,N

angle=2.0d+00*pi*((k-1.0d+00)/(N*1.0d+00))

CALL BZ_edge(angle,KFSmax)

kpx=KFSmax*cos(angle)
kpy=KFSmax*sin(angle)

write(5,*)kpx/BZconst,kpy/BZconst

end do

end subroutine Plot_BZ
!------------------------------------------------------------------
!------------------------------------------------------------------
! This subroutines define the regular grid:
!----------------------------------------------
! Inputs:
!----------------------------------------------
!1. N (integer) : number of kx/ky points
!----------------------------------------------
! Outputs:
!----------------------------------------------
!1. kxsamp (one dimensonal real array of size N) : kx values of regular grid
!2. kysamp (one dimensonal real array of size N) : ky values of regular grid 
!3. kxsamp_interp (one dimensonal real array of size N*N) : kx values of regular grid inside of Brilloune Zone
!4. kysamp_interp (one dimensonal real array of size N*N) : ky values of regular grid outside of Brilloune Zone
!------------------------------------------------------------------
Subroutine Sample_points(N,kxsamp,kysamp,kpxmat_interp,kpymat_interp)

use commons

IMPLICIT NONE

INTEGER,INTENT(IN)::N                                                      !number of kx/ky points
DOUBLE PRECISION,DIMENSION(1:N),INTENT(OUT)::kxsamp                        !kx values of regular grid  
DOUBLE PRECISION,DIMENSION(1:N),INTENT(OUT)::kysamp                        !ky values of regular grid                       
DOUBLE PRECISION,DIMENSION(1:(N*N)),INTENT(OUT)::kpxmat_interp             !kx values of regular grid inside of Brilloune Zone          
DOUBLE PRECISION,DIMENSION(1:(N*N)),INTENT(OUT)::kpymat_interp             !ky values of regular grid inside of Brilloune Zone
!-----------------------------------------
! Dummy variables:
!-----------------------------------------
DOUBLE PRECISION::kpx,kpy
INTEGER::j,kxint,kyint
!------------------------------------------------------------------
!------------------------------------------------------------------
! Sample points:
!------------------------------------------------------------------

j=1

Do kyint=1,N
Do kxint=1,N

kxsamp(kxint)=-BZconst+(2.0d+00*BZconst)*((kxint-1.0d+00)/(N-1.0d+00))
kysamp(kyint)=-((sqrt(3.0d+00))/2.0d+00)*BZconst+(2.0d+00*BZconst*(sqrt(3.0d+00))/2.0d+00)*((kyint-1.0d+00)/(N-1.0d+00))

kpx=kxsamp(kxint)
kpy=kysamp(kyint)

!----------------------------------------------

CALL Folding(kpx,kpy,kpx,kpy)

kpxmat_interp(j)=kpx
kpymat_interp(j)=kpy

!dissamp(l1,kxint,kyint)=dispersionana(l1,kxsamp(kxint),kysamp(kyint),0.0d+00)

j=j+1

end do
end do

end subroutine Sample_points
!------------------------------------------------------------------
!This subroutine writes the vertex tensor in a txt-file:
!----------------
!Inputs:
!----------------
!1. NEQ (integer) : dimension of the vertex
!2. Yvertexre (one dimensional array) : Vertex
!3. n1max (integer) : the first band index to which maximum vertex belongs to
!4. n1max (integer) : the second band index to which maximum vertex belongs to
!5. n1max (integer) : the third band index to which maximum vertex belongs to
!6. n1max (integer) : the forth band index to which maximum vertex belongs to
!------------------------------------------------------------------
Subroutine write_vertex(NEQ,Yvertex,n1max,n2max,n3max,n4max)

use commons

IMPLICIT NONE

!------------------------------------------------------------------
INTEGER,INTENT(IN)::NEQ                                  !dimension of vertex
DOUBLE PRECISION,DIMENSION(1:NEQ),INTENT(IN)::Yvertex    !one dimentional vertex
INTEGER,INTENT(IN)::n1max                                !the first band index to which maximum vertex belongs to
INTEGER,INTENT(IN)::n2max                                !the second band index to which maximum vertex belongs to
INTEGER,INTENT(IN)::n3max                                !the third band index to which maximum vertex belongs to
INTEGER,INTENT(IN)::n4max                                !the forth band index to which maximum vertex belongs to
!------------------------------------------------------------------
INTEGER::p1,p2,p3,n1,n2,n3,n4,i,imat,openstatus
!------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------
open(unit=1,file="vertex_test1_num_ana.txt",status="new",action="write",position="rewind",IOSTAT=openstatus)
!-----------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------
open(unit=11,file="vertex_maxvertex_sector_num_ana.txt",status="new",action="write",position="rewind",IOSTAT=openstatus)
!-----------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------
open(unit=12,file="vertex_intraband_spin_singlet_num_ana.txt",status="new",action="write",position="rewind",IOSTAT=openstatus)
!-----------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------
open(unit=13,file="vertex_intraband_spin_triplet_num_ana.txt",status="new",action="write",position="rewind",IOSTAT=openstatus)
!-----------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------
open(unit=14,file="vertex_interband_spin_singlet_num_ana.txt",status="new",action="write",position="rewind",IOSTAT=openstatus)
!-----------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------
open(unit=15,file="vertex_interband_spin_triplet_num_ana.txt",status="new",action="write",position="rewind",IOSTAT=openstatus)
!-----------------------------------------------------------------------------------------------------

Do n4=1,nband
Do n3=1,nband
Do n2=1,nband
Do n1=1,nband

Do p3=1,npatch
Do p2=1,npatch
Do p1=1,npatch

i=imat(n1,n2,n3,n4,p1,p2,p3)

write(1,*)Yvertex(i),p1,p2,p3,n1,n2,n3,n4


if(n1==n1max .and. n2==n2max .and. n3==n3max  .and. n4==n4max )then
write(11,*)Yvertex(i),p1,p2,p3
!print*,Yvertex(i),p1,p2,p3
endif

if(n1==n3 .and. n2==n4)then

if(2*int(n1/2) .ne. n1 .and. n1==n2)then
write(13,*)Yvertex(i),p1,p2,p3,n1,n2,n3,n4
!print*,'13:',n1,n2,n3,n4
endif

if(2*int(n1/2) .ne. n1 .and. n1==n2-1)then
write(12,*)Yvertex(i),p1,p2,p3,n1,n2,n3,n4
!print*,'12:',n1,n2,n3,n4
endif

if(2*int(n1/2) .ne. n1 .and. 2*int(n2/2) .ne. n2 .and. n2>n1)then
write(15,*)Yvertex(i),p1,p2,p3,n1,n2,n3,n4
!print*,'15:',n1,n2,n3,n4
endif

if(2*int(n1/2) .ne. n1 .and. 2*int(n2/2)==n2 .and. n2-n1>1)then
write(14,*)Yvertex(i),p1,p2,p3,n1,n2,n3,n4
!print*,'14:',n1,n2,n3,n4
endif

endif


end do !n4
end do !n3
end do !n2
end do !n1

end do !p3
end do !p2
end do !p1
!-------------------------------------------------------------------------------------------
end subroutine write_vertex
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
! This subroutine interpolates the DFT data to produce data on a 
! Regular grid
!---------------
! Inputs:
!---------------
! 1. flag (character) : 'real' does the interpolation, 'test' reads from file
! 2. N (integer) : number of DFT data points
! 3. Nsamp (integer) : square root of number of regular grid points
!---------------
! Outputs:
!---------------
! 1.dispmat_ex (3 dimensional double precision array) : the interpolated dispersion 
!   the first dimension specifies the band index 1,2,3,4,
!   the second and third dimension specifies kx and ky, respectively
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
subroutine Transform_regular_grid(flag,N,Nsamp,dissamp)
!-------------------------------------------------------------------------------------------

use commons

IMPLICIT NONE

CHARACTER*4,INTENT(IN)::flag                                                ! 'real' does interpolation, 'test' reads the interpolated data
INTEGER,INTENT(IN)::N                                                       ! number of dft k points
INTEGER,INTENT(IN)::Nsamp                                                   ! number of grid dimension

DOUBLE PRECISION,DIMENSION(1:nband,1:Nsamp,1:Nsamp),INTENT(OUT)::dissamp    ! Interpolated dispersion

!-------------------------------------------------------------------------
! Rest of Dummy variables:
!-------------------------------------------------------------------------
DOUBLE PRECISION,DIMENSION(1:N)::kpxmat,kpymat
DOUBLE PRECISION,DIMENSION(1:N)::dispmat

INTEGER::N_interp,openstatus,j,kxint,kyint,units,l1
DOUBLE PRECISION,DIMENSION(1:(Nsamp*Nsamp))::kpxmat_interp,kpymat_interp,dispmat_interp &
                                             ,kpxmat_ex,kpymat_ex,dispmat_ex
DOUBLE PRECISION,DIMENSION(1:Nsamp)::kxsamp,kysamp

DOUBLE PRECISION::p_interp,kpx,kpy,discheck,dispmin,dispmax,mom,angle,dispersionana

!------------------------------------------------------------------

if(flag=='real')then

open(unit=61,file="band1_interp.txt",status="new",action="write",position="rewind",IOSTAT=openstatus)
open(unit=62,file="band2_interp.txt",status="new",action="write",position="rewind",IOSTAT=openstatus)
open(unit=63,file="band3_interp.txt",status="new",action="write",position="rewind",IOSTAT=openstatus)
open(unit=64,file="band4_interp.txt",status="new",action="write",position="rewind",IOSTAT=openstatus)

elseif(flag=='test')then

open(unit=61,file="band1_interp.txt",status="old",action="read",position="rewind",IOSTAT=openstatus)
open(unit=62,file="band2_interp.txt",status="old",action="read",position="rewind",IOSTAT=openstatus)
open(unit=63,file="band3_interp.txt",status="old",action="read",position="rewind",IOSTAT=openstatus)
open(unit=64,file="band4_interp.txt",status="old",action="read",position="rewind",IOSTAT=openstatus)

endif

!-----------------------------------------------------------------------------------------------------
!open(unit=2,file="Sep12test2.txt",status="new",action="write",position="rewind",IOSTAT=openstatus)
!-----------------------------------------------------------------------------------------------------

CALL Sample_points(Nsamp,kxsamp,kysamp,kpxmat_interp,kpymat_interp)

N_interp=Nsamp*Nsamp

Do l1=1,nband/2
!------------------------------------------------------------------
! transforming to regular grid:
!------------------------------------------------------------------
if(flag=='real')then

CALL READ_DATA(l1,N,kpxmat,kpymat,dispmat,dispmin,dispmax)

p_interp=16.0d+00

CALL shepard_interp_2d(N,kpxmat,kpymat,dispmat,p_interp,N_interp,kpxmat_interp,kpymat_interp,dispmat_interp)

!-------------------------------------------------------------
Do j=1,N_interp

mom=sqrt(((kpxmat_interp(j))**2)+((kpymat_interp(j))**2))

angle=atan2(kpymat_interp(j),kpxmat_interp(j))
angle=(-SIGN(1.0d+00,angle)+1.0d+00)*(pi)+angle

!dispmat_interp(j)=dispersionana(l1*2,mom,angle,0.0d+00)


end do !j
!-------------------------------------------------------------

dispmat_ex(1:N_interp)=dispmat_interp(1:N_interp)

units=60+l1

Do j=1,N_interp

write(units,*)dispmat_interp(j)

end do !j

!------------------------------------------------------------------
elseif(flag=='test')then

units=60+l1

Do j=1,N_interp

read(units,*)discheck  

kpxmat_ex(j)=kpxmat_interp(j)
kpymat_ex(j)=kpymat_interp(j)

dispmat_ex(j)=discheck

end do

endif
!------------------------------------------------------------------
!------------------------------------------------------------------
j=1
Do kyint=1,Nsamp
Do kxint=1,Nsamp

kpxmat_ex(j)=kxsamp(kxint)
kpymat_ex(j)=kysamp(kyint)

kpx=kpxmat_ex(j)
kpy=kpymat_ex(j)

discheck=dispmat_ex(j)

dissamp((2*l1-1):(2*l1),kxint,kyint)=discheck

j=j+1

end do
end do
!------------------------------------------------------------------
end do
!------------------------------------------------------------------
end Subroutine Transform_regular_grid
!------------------------------------------------------------------
!---------------------------------------------------------------------
! This subroutine reads the DFT-data
!-----------------------
!  Inputs:
!-----------------------
!  1.nbl (integer) : band-index in [1,4]
!  2.Ndata (integer) :number of DFT-data points
!-----------------------
!  Outputs:
!-----------------------
!  1.kpxmat (one dimensional array of size Ndata) : kx values
!  2.kpymat (one dimensional array of size Ndata) : ky values
!  3.dispmat (one dimensional array of size Ndata) : dispersion values
!  4.dispmin (real) : minimum value of dispersion 
!  5.dispmax (real) : maximum value of dispersion 
!---------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------

Subroutine READ_DATA(nbl,Ndata,kpxmat,kpymat,dispmat,dispmin,dispmax)

use commons

IMPLICIT NONE

!--------------------------------------------------------------------------------------
INTEGER,INTENT(IN)::nbl                                    !band index in [1,4]
INTEGER,INTENT(IN)::Ndata                                  ! number of K points

DOUBLE PRECISION,DIMENSION(1:Ndata),INTENT(OUT)::kpxmat    ! kx values
DOUBLE PRECISION,DIMENSION(1:Ndata),INTENT(OUT)::kpymat    ! ky values
DOUBLE PRECISION,DIMENSION(1:Ndata),INTENT(OUT)::dispmat   ! dispersion
DOUBLE PRECISION,INTENT(OUT)::dispmin                      ! min-value of dispersion
DOUBLE PRECISION,INTENT(OUT)::dispmax                      ! max-value of dispersion
!--------------------------------------------------------------------------------------
! Dummy variables:
!--------------------------------------------------------------------------------------

INTEGER::openstatus,im
DOUBLE PRECISION::kpxm,kpym,discheckm,kpxmin,kpxmax,kpymin,kpymax,mom,angle,dispersionana

INTEGER::unitm               ! Index to read files

unitm=40+nbl

!-----------------------------------------------------------------------------------------------------
! Files of num-dispersion:
!-----------------------------------------------------------------------------------------------------
open(unit=40,file="ALL_kpoints.dat",status="old",action="read",position="rewind",IOSTAT=openstatus)
!-----------------------------------------------------------------------------------------------------
open(unit=41,file="ALL_band_1.dat",status="old",action="read",position="rewind",IOSTAT=openstatus)
open(unit=42,file="ALL_band_2.dat",status="old",action="read",position="rewind",IOSTAT=openstatus)
open(unit=43,file="ALL_band_3.dat",status="old",action="read",position="rewind",IOSTAT=openstatus)
open(unit=44,file="ALL_band_4.dat",status="old",action="read",position="rewind",IOSTAT=openstatus)
!-----------------------------------------------------------------------------------------------------

Do im=1,Ndata

read(40,*)kpxm,kpym
read(unitm,*)discheckm

discheckm=discheckm*100.0d+00

kpxmat(im)=kpxm
kpymat(im)=kpym

dispmat(im)=discheckm

if(im==1)then
kpxmin=kpxm
kpxmax=kpxm

kpymin=kpym
kpymax=kpym

dispmin=discheckm
dispmax=discheckm

endif

!--------------------------------
! Kpxmin:
!--------------------------------
if(kpxmat(im)<kpxmin)then
kpxmin=kpxmat(im)
endif
!--------------------------------
! kpxmax:
!--------------------------------
if(kpxmat(im)>kpxmax)then
kpxmax=kpxmat(im)
endif
!--------------------------------
! kpymin:
!--------------------------------
if(kpymat(im)<kpymin)then
kpymin=kpymat(im)
endif
!--------------------------------
! kpymax:
!--------------------------------
if(kpymat(im)>kpymax)then
kpymax=kpymat(im)
endif
!--------------------------------
! dispmin:
!--------------------------------
if(dispmat(im)<dispmin)then
dispmin=dispmat(im)
endif
!--------------------------------
! dispmax:
!--------------------------------
if(dispmat(im)>dispmax)then
dispmax=dispmat(im)
endif
!--------------------------------
end do
!------------------------------------------------------------------
!------------------------------------------------------------------
! Making convensions consistent:
!------------------------------------------------------------------
kpxmat(1:Ndata)=kpxmat(1:Ndata)*(BZconst)/abs(kpxmin)
kpymat(1:Ndata)=kpymat(1:Ndata)*((sqrt(3.0d+00))*BZconst/2.0d+00)/abs(kpymin)


Do im=1,Ndata


mom=sqrt(((kpxmat(im))**2)+((kpymat(im))**2))

angle=atan2(kpymat(im),kpxmat(im))
angle=(-SIGN(1.0d+00,angle)+1.0d+00)*(pi)+angle

dispmat(im)=dispersionana(nbl*2,mom,angle,0.0d+00)

end do


end Subroutine READ_DATA
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------
! This subroutine checks the interpolated dispersion with DFT data
!---------------------------------------------------------------------
!  Inputs:
!---------------------------------------------------------------------
!  1.Ndata (integer) :number of DFT-data points
!---------------------------------------------------------------------
!  Outputs:
!---------------------------------------------------------------------
!  4 files 'comparison_band1-4.txt'
!---------------------------------------------------------------------
!--------------------------------------------------------------------------------------
Subroutine Comp_interp_DFT(N)

use commons

IMPLICIT NONE

INTEGER,INTENT(IN)::N

INTEGER::i,l1,unitm,openstatus
DOUBLE PRECISION::kx,ky,dispersion_theta,dispersion,dispmin,dispmax&
                 ,disp_dft,disp_interp,angle,mom
DOUBLE PRECISION,DIMENSION(1:N)::kpxmat,kpymat,dispmat

!--------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------
open(unit=91,file="comparison_band1.txt",status="new",action="write",position="rewind",IOSTAT=openstatus)
open(unit=92,file="comparison_band2.txt",status="new",action="write",position="rewind",IOSTAT=openstatus)
open(unit=93,file="comparison_band3.txt",status="new",action="write",position="rewind",IOSTAT=openstatus)
open(unit=94,file="comparison_band4.txt",status="new",action="write",position="rewind",IOSTAT=openstatus)
!-----------------------------------------------------------------------------------------------------

Do l1=1,4

CALL READ_DATA(l1,N,kpxmat,kpymat,dispmat,dispmin,dispmax)

unitm=90+l1

Do i=1,N

kx=kpxmat(i)
ky=kpymat(i)

disp_dft=dispmat(i)

mom=sqrt((kx**2)+(ky**2))

angle=atan2(ky,kx)
angle=(-SIGN(1.0d+00,angle)+1.0d+00)*(pi)+angle

disp_interp=dispersion(2*l1,mom,angle,0.0d+00)

write(unitm,*)kx/BZconst,ky/BZconst,disp_dft,disp_interp,disp_dft-disp_interp

end do
end do

end Subroutine Comp_interp_DFT
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
! This subroutine puts the graphene dispersion on a square grid 
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!  Inputs:
!--------------------------------------------------------------------------------------
!  1. N : number of k-points
!--------------------------------------------------------------------------------------
!  Outputs:
!--------------------------------------------------------------------------------------
!  one file 'ALL_kpoints.dat' that contains the k points with format kx and ky
!  and four other files 'ALL_band_1-4.dat' including the eigenvalues (of H_0)
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
Subroutine GENERATE_DATA(N)

use commons

IMPLICIT NONE

INTEGER,INTENT(IN)::N

!--------------------------------------------------------------
! Dummy variables:
!--------------------------------------------------------------
INTEGER::kxint,kyint,l1,unitm,openstatus,j
DOUBLE PRECISION::kx,ky,dispersion,dispersionana,mom,angle
DOUBLE PRECISION,DIMENSION(1:N)::kxsamp,kysamp
DOUBLE PRECISION,DIMENSION(1:N*N)::kpxmat_interp,kpymat_interp
!-----------------------------------------------------------------------------------------------------
open(unit=40,file="ALL_kpoints.dat",status="new",action="readwrite",position="rewind",IOSTAT=openstatus)
!-----------------------------------------------------------------------------------------------------
open(unit=41,file="ALL_band_1.dat",status="new",action="readwrite",position="rewind",IOSTAT=openstatus)
open(unit=42,file="ALL_band_2.dat",status="new",action="readwrite",position="rewind",IOSTAT=openstatus)
open(unit=43,file="ALL_band_3.dat",status="new",action="readwrite",position="rewind",IOSTAT=openstatus)
open(unit=44,file="ALL_band_4.dat",status="new",action="readwrite",position="rewind",IOSTAT=openstatus)
!-----------------------------------------------------------------------------------------------------

CALL Sample_points(N,kxsamp,kysamp,kpxmat_interp,kpymat_interp)

open(unit=40,file='ALL_kpoints.dat')

open(unit=41,file='ALL_band_1.dat')
open(unit=42,file='ALL_band_2.dat')
open(unit=43,file='ALL_band_3.dat')
open(unit=44,file='ALL_band_4.dat')


Do j=1,N*N

kx=kpxmat_interp(j)
ky=kpymat_interp(j)

!----------------------------------------------
!----------------------------------------------

write(40,*)kx,ky

mom=sqrt((kx**2)+(ky**2))

angle=atan2(ky,kx)
angle=(-SIGN(1.0d+00,angle)+1.0d+00)*(pi)+angle

Do l1=1,nband/2

dispersion=0.01d+00*dispersionana(l1*2,mom,angle,0.0d+00)

unitm=40+l1

write(unitm,*)dispersion

end do !l1

end do !j


CLOSE(40)
CLOSE(41)
CLOSE(42)
CLOSE(43)
CLOSE(44)
!--------------------------------------------------------------
end Subroutine GENERATE_DATA
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
! This subroutine puts the graphene dispersion on a given grid
!           which will be imported as .dat file 
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!  Inputs:
!--------------------------------------------------------------------------------------
!  1. N : number of k-points
!--------------------------------------------------------------------------------------
!  Outputs:
!--------------------------------------------------------------------------------------
!  one file 'ALL_kpoints.dat' that contains the k points with format kx and ky
!  and four other files 'ALL_band_1-4.dat' including the eigenvalues (of H_0)
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
Subroutine GENERATE_DATA_GIVEN_GRID(N)

use commons

IMPLICIT NONE

INTEGER,INTENT(IN)::N

!--------------------------------------------------------------
! Dummy variables:
!--------------------------------------------------------------
INTEGER::kxint,kyint,l1,unitm,openstatus,j
DOUBLE PRECISION::kx,ky,dispersion,dispersionana,mom,angle,kpxmin,kpymin
DOUBLE PRECISION,DIMENSION(1:N)::kpxmat_re,kpymat_re
!-----------------------------------------------------------------------------------------------------
open(unit=40,file="ALL_kpoints.dat",status="old",action="read",position="rewind",IOSTAT=openstatus)
!-----------------------------------------------------------------------------------------------------
open(unit=41,file="ALL_band_1.dat",status="new",action="readwrite",position="rewind",IOSTAT=openstatus)
open(unit=42,file="ALL_band_2.dat",status="new",action="readwrite",position="rewind",IOSTAT=openstatus)
open(unit=43,file="ALL_band_3.dat",status="new",action="readwrite",position="rewind",IOSTAT=openstatus)
open(unit=44,file="ALL_band_4.dat",status="new",action="readwrite",position="rewind",IOSTAT=openstatus)
!-----------------------------------------------------------------------------------------------------

open(unit=40,file='ALL_kpoints.dat')

open(unit=41,file='ALL_band_1.dat')
open(unit=42,file='ALL_band_2.dat')
open(unit=43,file='ALL_band_3.dat')
open(unit=44,file='ALL_band_4.dat')

!----------------------------------------------
Do j=1,N
!----------------------------------------------

read(40,*)kx,ky

kpxmat_re(j)=kx
kpymat_re(j)=ky

if(j==1)then
kpxmin=kpxmat_re(1)
kpymin=kpymat_re(1)
endif

if(kpxmat_re(j) < kpxmin)then
kpxmin=kpxmat_re(j)
endif

if(kpymat_re(j) < kpymin)then
kpymin=kpymat_re(j)
endif

end do
!----------------------------------------------
!print*,'kpxmin=',kpxmin
!print*,'kpymin=',kpymin
!print*,'kpymin/kpxmin=',kpymin/kpxmin,sqrt(3.0d+00)/2.0d+00
!----------------------------------------------

kpxmat_re(1:N*N)=kpxmat_re(1:N*N)*(BZconst)/abs(kpxmin)
kpymat_re(1:N*N)=kpymat_re(1:N*N)*((sqrt(3.0d+00))*BZconst/2.0d+00)/abs(kpymin)

Do j=1,N

kx=kpxmat_re(j)
ky=kpymat_re(j)

mom=sqrt((kx**2)+(ky**2))

angle=atan2(ky,kx)
angle=(-SIGN(1.0d+00,angle)+1.0d+00)*(pi)+angle

Do l1=1,nband/2

dispersion=0.01d+00*dispersionana(l1*2,mom,angle,0.0d+00)

unitm=40+l1

write(unitm,*)dispersion

end do !l1

end do !j


CLOSE(40)
CLOSE(41)
CLOSE(42)
CLOSE(43)
CLOSE(44)
!--------------------------------------------------------------
end Subroutine GENERATE_DATA_GIVEN_GRID
!--------------------------------------------------------------------------------------
!-----------------------------------------------------------------------
! This subroutine provides the information,i.e.,
! the number of times the dipsersion changes sign,
! and the extermums of dispersion
! in order to find the Fermi Surfaces!
!--------------------------------------------------------
! Inputs:
!--------------------------------------------------------
! 1. l1 : band index
! 2. dimwork : number of k values for a given angle
! 3. N_a : number of angles
!--------------------------------------------------------
! Outputs:
!--------------------------------------------------------
! 1. numsign : number of times dispersion changes sign for a given angle
! 2. KFSext : extermums of the dispersion 
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
Subroutine Extermum_Ge(blab,dimwork,N_a,nsign_mat,KFSext)

use commons

implicit none

!---------------------------------------------------------------
INTEGER,INTENT(IN)::blab 
INTEGER,INTENT(IN)::dimwork
INTEGER,INTENT(IN)::N_a


INTEGER,DIMENSION(1:N_a),INTENT(OUT)::nsign_mat
DOUBLE PRECISION,DIMENSION(1:10,1:N_a),INTENT(OUT)::KFSext

!---------------------------------------------------------------
! Dummy variables:
!---------------------------------------------------------------

DOUBLE PRECISION::angle
INTEGER::nsign

DOUBLE PRECISION::e_ext,e_bound

DOUBLE PRECISION,DIMENSION(1:dimwork)::dispmat_ex,ddispmat_ex

DOUBLE PRECISION::dispersion_theta,KFS,KFSmax,dispersion,angle_min,angle_max
INTEGER::i,ncrit,j,openstatus,k,unitm,i_min,i_max

DOUBLE PRECISION,DIMENSION(1:Nnumdisp)::kxdftmat,kydftmat,dispdftmat
DOUBLE PRECISION,DIMENSION(1:Nnumdisp)::kpxmat_interp,kpymat_interp,dispmat_interp
DOUBLE PRECISION::p_interp,dispcheck,kpxmin,kpx,kpy
!------------------------------------------------------------------

Do k=1,N_a

angle=(k-1.0d+00)*(pi/3.0d+00)/N_a

ncrit=0
nsign=1

CALL BZ_Edge(angle,KFSmax)

Do i=1,50

KFS=KFSmax*(i)/50.0d+00


dispmat_ex(i)=dispersion_theta(blab,KFS,angle,0.0d+00)!dispmat_interp(1)


if(i>1)then

ddispmat_ex(i)=(dispmat_ex(i)-dispmat_ex(i-1))/(KFSmax/50.0d+00)


endif

end do

!------------------------------------------------------------------
Do i=2,49

if(i>1)then

KFS=KFSmax*i/50.0d+00


if(ddispmat_ex(i)*ddispmat_ex(i+1)<0.0d+00)then

Do j=1,2

if(ddispmat_ex(i)*ddispmat_ex(i+1+j)<0.0d+00)then

ncrit=ncrit+1

endif  !changing sign with j add


if(ncrit>1)then

KFSext(nsign,k)=KFS

nsign=nsign+1

endif  !ncrit

enddo

else

ncrit=0

endif  !changing sign

endif !i>1


if(nsign>5)then
print*,'careful! more than five Fermi Surfaces!'
end if


if(nsign==1)then

CALL BZ_edge(angle,KFS)

KFSext(1,k)=KFS

else

KFSext(2,k)=KFS

endif

end do ! i


end do ! k

!----------------------------------------------------------------------
end Subroutine Extermum_Ge
!----------------------------------------------------------------------
