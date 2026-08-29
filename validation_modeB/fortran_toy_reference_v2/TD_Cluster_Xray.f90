program main
  implicit none
  include 'mpif.h'
  integer, parameter :: root=0
  integer, parameter :: ncellsx=1
  integer, parameter :: ncellsy=1
  integer, parameter :: ncells=ncellsx*ncellsy
  !  integer, parameter :: ncells=8 ! uncomment to use Betts 8A cluster. Also comment out above line and set flag_B ne 1
  integer, parameter :: twistx=1    ! number of twists in x direction 
  integer, parameter :: twisty=1    ! number of twists in y direction
  integer, parameter :: norbstm=5
  integer, parameter :: norbsl=6   ! number of ligand orbitals
  !  integer, parameter :: norbsl=0    ! if this is set to zero there will be warning messages when compiling but it can be ignored
  integer, parameter :: norbspc=norbstm+norbsl
  integer, parameter :: norbs=ncells*norbspc
  integer, parameter :: norbspc_max=11
  integer, parameter :: norbs_max=ncells*norbspc_max
  integer, parameter :: nholes=2
  integer, parameter :: nd=2              ! gives the target configuration d^(nd) in hole language
  integer, parameter :: qnumbers=ncells*(nholes+1)     ! if dlso=0 - no spin-orbit coupling in valence shell
  !  integer, parameter :: qnumbers=ncells                ! if dlso/=0 - spin-orbit coupling in valence shell
  integer, parameter :: max_array=2147483647
!  integer, parameter :: nfreq=1600
!  integer, parameter :: ifreq=180
  integer, parameter :: nfreq=16
  integer, parameter :: ifreq=16
  integer, parameter :: flag_w=0    ! flag_w eq 1 turns on writing out all eigenvalues for Nholes
  integer, parameter :: flag_B=1    ! flag_B ne 1 turns on Betts cluster
  integer, parameter :: flag_R=2    ! flag_R ne 1 turns on XAS/RIXS calculation
  integer, parameter :: flag_X=1    ! flag_X ne 1 turns on XPS calculation
  integer, parameter :: flag_skip_XAS = 1 ! flag_skip_XAS ne 1 will skip XAS when doing a XAS/RIXS calculation
  integer, parameter :: stepmax=2   ! 1 diagonalizes Nup + Ndown, Norbs Hamiltonian
                                    ! 2 diagonalizes 1 plus Nup-1 + Ndown, Norbs hole removal Hamiltonian (XAS+RIXS)
                                    ! 3 diagonalizes N+1 hole addition Hamiltonian (XPS) - THIS IS CURRENTLY NOT WORKING FOR THIS CODE
  integer, parameter :: edge=1      ! 1 is for Tm L-edge, 2 ligand K-edge
  integer, parameter :: ncore=6     ! 6 for Tm L-edge, 4 for ligand K-edge
  integer, parameter :: usearpack=3 ! set usearpack to > 2 to switch ARPACK off
                                    ! 0 (just valence)
                                    ! 1 (both valence and core)
                                    ! 2 (just core) states
   integer, parameter :: timesteps=80
!  integer, parameter :: timesteps=80
!  integer, parameter :: timesteps=6
  integer, parameter :: nb=1
  integer, parameter :: nt=1
  real, parameter :: small=0.000001
  ! CuO parameters
!  double precision, parameter :: e2=2.2d0  ! d_z^2 site energy   
!  double precision, parameter :: e3=1.0d0  ! d_xy   
!  double precision, parameter :: e4=1.5d0  ! d_xz  
!  double precision, parameter :: e5=1.5d0  ! d_yz    
  ! NiO parameters
  double precision, parameter :: e2=0.0d0  ! d_z^2 site energy   
  double precision, parameter :: e3=1.0d0  ! d_xy   
  double precision, parameter :: e4=1.0d0  ! d_xz  
  double precision, parameter :: e5=1.0d0  ! d_yz    
  double precision, parameter :: f0=6.0d0                 ! valence Coulomb matrix elements
  double precision, parameter :: f2=0.12857d0             ! F_2
  double precision, parameter :: f4=0.025d0               ! F_4
  double precision, parameter :: raca=f0-dfloat(49)*f4
  double precision, parameter :: racb=f2-dfloat(5)*f4
  double precision, parameter :: racc=dfloat(35)*f4
  double precision, parameter :: aveu=raca-dfloat(14)*racb/dfloat(9)+dfloat(7)*racc/dfloat(9)
  
  double precision, parameter :: cf0=f0         ! core hole Coulomb matrix elements
  double precision, parameter :: cf2=6.18d0     ! F^2
  double precision, parameter :: cf3=4.635d0    ! G^1
  double precision, parameter :: cf4=2.63d0     ! G^3
  
  double precision, parameter :: upp0=1.0d0      ! ligand F0
  double precision, parameter :: upp2=0.1d0      ! ligand F2
  double precision, parameter :: uppf=upp0       ! ligand F0 1s-2p interaction
  double precision, parameter :: uppg=0.25d0     ! ligand G1 1s-2p interaction
  double precision, parameter :: lso=11.5d0      ! core spin-orbit coupling
  double precision, parameter :: dlso=0.0d0      ! valence spin-orbit coupling -> requires larger matrix diagonalization
  ! just for CuO
!  double precision, parameter :: delxx=2.8d0 ! Chunjing - 2.522
!  double precision, parameter :: delxy=2.8d0 ! Chunjing - 0.948
!  double precision, parameter :: delxz=2.8d0 ! Chunjing - 1.093
  ! for all others
  double precision, parameter :: delxx=nd*aveu 
  double precision, parameter :: delxy=nd*aveu  
  double precision, parameter :: delxz=nd*aveu 
  double precision, parameter :: delyx=delxy
  double precision, parameter :: delyy=delxx
  double precision, parameter :: delyz=delxz
  double precision, parameter :: e6=delxx      ! ligand x p_x
  double precision, parameter :: e7=delxy      ! ligand x p_y
  double precision, parameter :: e8=delxz      ! ligand x p_z
  double precision, parameter :: e9=delyx      ! ligand y p_x
  double precision, parameter :: e10=delyy     ! ligand y p_y
  double precision, parameter :: e11=delyz     ! ligand y p_z
  double precision, parameter :: vpds=2.0d0/dsqrt(3.0d0)
  double precision, parameter :: vpdp=0.225d0
!  double precision, parameter :: vpdp=0.0d0
  double precision, parameter :: vpps=0.0d0
  double precision, parameter :: vppp=0.25d0
!  double precision, parameter :: vppp=0.0d0
!  double precision, parameter :: tpd=0.0d0     ! Chunjing - 1.2 
  double precision, parameter :: tpd=dsqrt(3.0d0)*vpds/2.0d0                ! hybridization between d_x^2-y^2 and ligand p_xx and p_yy
  double precision, parameter :: tpdz=tpd/4.0d0   ! Chunjing - 0.3          ! hybridization between d_z^2 and ligand p_xx and p_yy
!  double precision, parameter :: tpdz=vpds/2.0d0
!  double precision, parameter :: tpdxy=0.225d0  ! Chunjing - 0.8            ! hybridization between d_xy and ligand p_xy and p_yx
  double precision, parameter :: tpdxy=vpdp
!  double precision, parameter :: tpdxz=0.72d0  ! Chunjing - 0.8
  double precision, parameter :: tpdxz=vpdp                                 ! hybridization between d_xz and ligand p_xz
!  double precision, parameter :: tpdyz=0.72d0  ! Chunjing - 0.8
  double precision, parameter :: tpdyz=vpdp                                 ! hybridization between d_yz and ligand p_yz
  double precision, parameter :: i4=-1.0d0      ! should be set to -1 if crystal field energies are used as input and not calculated from point charge model
  double precision, parameter :: theta_sc=149.5d0/dfloat(2*180) ! in units of pi 
  double precision, parameter :: ecore=100.0d0
  double precision, parameter :: pump_freq=0.5d0
  double precision, parameter :: pump_ampx=10.0d0
  double precision, parameter :: pump_ampy=0.0d0
  double precision, parameter :: pump_sigma=1.5d0
  
  integer momx(ncells),momy(ncells),cstates(ncells),cdim_max,box(2*norbs),hNEV,countj(2*norbs+ncore*ncells),l1,l2,l3,l4
  integer, allocatable, dimension(:) :: mlabel,slabel,uflag,mflag,sflag,nl,cl,nf,cf,fsign,cmlabel,cbinlabel,hmlabel,hslabel,gsign,lflag,mlflag, allowed_tm
  integer, allocatable, dimension(:) :: IndexI,IndexJ,tIndexI,tIndexJ
  integer, allocatable, dimension(:,:) :: vect,basis,lap,clap,dh,ph,flap,degen,sz1,sz2,sb4,binlabel,cbasis,hbasis,hbinlabel
  integer, allocatable, dimension(:,:,:) :: cpairs,cothers,dpairs,dothers
  integer, allocatable, dimension(:,:,:,:) :: hpairs,hothers
  integer aaa(2*norbs),bbb(2*norbs),alpha,gsec(qnumbers),statemax,iwork,ccc(2*norbs+ncore*ncells),ddd(2*norbs+ncore*ncells),tmom,tspin,hstatemax
  integer state,i,j,k,l,ip,ndimu,ndimd,n,m,p,count,ndim_max,nu,ws,nc,nevnmax,size_flag,dflag,aa(nholes-1),bb(nholes-1),gam,del,tlabel
  integer info,ii,jj,kk,ll,prod,factor,count2,count4,count6,nev,nevn,ndim,nmax,iup,fup,ncv,gstatedegen,twist1,twist2,beta,bella,reid
  integer momgs,wsp,gdegen(qnumbers),inx,iny,inz,fnx,fny,fnz,istate(qnumbers),spr_size,mtot(qnumbers),sztot(qnumbers),qqq,qstates,ppp
  integer mdiff(ncells,ncells),madd(ncells,ncells),symm,npairs,pp,factor1,factor2,dstates(ncells,nholes+1),binstates(qnumbers),qlist(qnumbers)
  integer hdim_max,hstates(ncells),hnpairs,hmtot(ncells),hqlist(ncells),hbinstates(ncells),hqstates
  integer haa(nholes),hbb(nholes),extk,xpslist1(qnumbers*ncells),xpslist2(qnumbers*ncells),xps1(2*norbs+ncore*ncells,1000),xps2(2*norbs+ncore*ncells,1000),xps3(2*norbs+ncore*ncells,1000)
  integer u1,u2,u3,u4,outout,count3
  integer i_time,ti,tj,qq,tm,i_time1,i_time2,t1,t2,t3,t4,tm_i
  integer(8) cmax,mpimax
  integer siz, myrank, numproc, mpistat
  real test,bob,time1,time2,time3,time4,Ylm,nchecks(16),time5,time6,smax
  real, allocatable, dimension(:,:) :: sz,sp,sm,spin,scorb
  real, allocatable, dimension(:,:,:) :: sorb
  double precision pi,omega,omegai,gamma,sum2,gammax,sump,a(qnumbers),orb(5),xpsmatrix,hint(0:2,-2:2,0:2,-2:2,0:2,-2:2,0:2,-2:2),omegai_max
  double precision kx,ky,egs,emax,cint(norbspc+ncore/2,norbspc+ncore/2,norbspc+ncore/2,norbspc+ncore/2),binomial,qx,qy,px,py,sum_max,XAS_max,term,pol_in(3),pol_out(3)
  double precision gaunt,xtal(2*norbs+ncore*ncells,2*norbs+ncore*ncells),irixs(ncells,nfreq),coree,hmin,hmax
  double precision sum,lor,gint(norbspc,norbspc,norbspc,norbspc),tmatx(norbspc_max,norbspc_max),tmaty(norbspc_max,norbspc_max),dso(2*norbspc,2*norbspc),phi,theta_in,theta_out,ephi,etheta
  double precision t(2*norbs_max+ncore*ncells,2*norbs_max+ncore*ncells),rixs,omegai_min,XAS_maxx,tmato(norbspc_max,norbspc_max)
  double precision core(6,6),xcef(norbspc_max,norbspc_max),pgsymm(5,5),spectra(qnumbers,3,3,qnumbers,ifreq,nfreq),xps(nfreq),occ(norbspc,norbspc)
  double precision p_time,vec_px(-nb*timesteps:nt*timesteps),vec_py(-nb*timesteps:nt*timesteps),pp_time,p1_time,p2_time,p3_time,p4_time
  double precision o1,o2,o3,o4,o5,o6,o7,o8,o9,o10,o11,epsilon,gauss,sigma,ttest,s1,s2,tmeas,s3,s4,omegas,sum4,tch,l12

  double precision, allocatable, dimension(:) :: w,work,sparseH,tsparseH,dn,DD,edn
  double precision, allocatable, dimension(:,:) :: vn,h,ZZ,toto1,toto2,toto3,deig,hdeig,evn,edeig
  double precision, allocatable, dimension(:,:,:) :: v,vh,ev
  double complex, allocatable, dimension(:) :: ioverlap,foverlap,d_mat,u11
  double complex, allocatable, dimension(:,:) :: t1_evol,t2_evol,tkern,biga,bigb,bigc,bigd,bige
  double complex, allocatable, dimension(:,:,:) :: fvert,ivert,evolution,tvn,tevolution
  double complex, allocatable, dimension(:,:,:,:) :: tv,tfvert,tivert
  double complex, allocatable, dimension(:,:,:,:,:) :: tpivert
  double complex, allocatable, dimension(:,:,:,:,:,:,:) :: ckern
  double complex csum,ksum,b_occ(norbspc,norbspc),csum4
  complex :: eye,gammakm(4,-4:4),ccsum
  character jobz,uplo
  character*20 :: filename
  !
  !  This code evaluates XAS/RIXS in one shot. It has three parts.
  !  It first calculates the ground and excited states of the valence Hamiltonian, labeled by quantum numbers (1)total momentum and (2)total S_z spin
  !  Then it calculates the core hole intermediate states for the core Hamiltonian
  !  Then it assembles the XAS/RIXS spectra using Kramers-Heisenberg formula
  !  Twisted average boundary conditions can be employed
  !
  !  to compile on Unix:  gfortran-mp-8 -O2 -o code -L/usr/lib -llapack name.f90
  !
  !  to compile on a Mac: mpif90 -O2 -o code -ffree-line-length-none -fcheck=all -L/usr/local/lib -larpack -framework accelerate TD_Cluster_Xray.f90
  !
  ! error trapping
  !
  if(ncells*nd.ne.nholes)then
     write(6,*)'Configuration set up error: Ncells x nd =/ nholes'
     goto 999
!     write(6,*)'But we are proceeding anyway...'
  end if
 
  if(raca.lt.0.0d0)then
     write(6,*)'Error: Racah A is negative'
     goto 999
  end if

  if(racb.lt.0.0d0)then
     write(6,*)'Error: Racah B is negative'
     goto 999
  end if
  
  if(dlso.ne.0.0d0)then
     if(qnumbers.ne.ncells)then
        write(6,*)'This is not set up correctly for valence spin-orbit coupling. Sz is not a good quantum number.'
        goto 999
     end if
  end if
  
  if(edge.eq.1)then
     if(ncore.ne.6)then
        write(6,*)'Error: edge and ncore not set properly'
        goto 999
     end if
  end if
  if(edge.eq.2)then
     if(ncore.ne.4)then
        write(6,*)'Error: edge and ncore not set properly'
        goto 999
     end if
  end if

  if(edge.gt.2.or.edge.lt.1)then
     write(6,*)'Error: edge not set properly'
     goto 999
  end if

  if(flag_B.ne.1)then
     write(6,*)'Error: Betts cluster not coded for this app'
     goto 999
  end if

  if(ncells.eq.1)then
     if(nholes.ne.nd)then
        write(6,*)'Not set up for proper n_hole configuration: n_hole is not equal to nd'
        goto 999
     end if
  end if

  if(flag_R.ne.1.and.stepmax.ne.2)then
     write(6,*)'Error: Not set properly for XAS/RIXS calculation'
     goto 999
  end if
  
  if(flag_X.ne.1.and.stepmax.ne.3)then
     write(6,*)'Error: Not set properly for XPS calculation'
     goto 999
  end if

  !
  !         CODE SET UP STARTS HERE
  !
  call CPU_time(time1)
  call mpi_init(mpistat)
  call mpi_comm_rank(mpi_comm_world,myrank,mpistat)
  call mpi_comm_size(mpi_comm_world,numproc,mpistat)
  print*, 'Process', myrank, ' of ', numproc, ' started.'

  if(flag_R.eq.1.and.numproc.gt.1)then
     if(myrank.eq.0)then
        write(6,*)'Parallel code only is invoked when calculating XAS & RIXS.',numproc
        write(6,*)'Run again with numproc=1.'
     end if
     goto 999
  end if

  pi=2.0d0*dasin(1.0d0)

  ephi=0.0d0
  etheta=0.0d0

  call gethint(f0,f2,f4,cf0,cf2,cf3,cf4,upp0,upp2,uppf,uppg,hint)

  if(nholes.gt.1)then
     allocate(fsign((nholes)*(nholes-1)/2))
     call getfsign(nholes,fsign)
  else
     allocate(fsign(1))
     fsign(1)=0
  end if

  eye=dcmplx(0.0d0,1.0d0)
  jobz='V'
  uplo='U'

  if(myrank.eq.0)write(6,*)
  open(unit=16,file='output.dat')

  core=dfloat(0)

  if(ncore.eq.6)then
     do i=1,6
        if(i.le.3)then
           core(i,i)=-0.5d0*lso*dfloat(i-2)
        else
           core(i,i)=0.5d0*lso*dfloat(i-5)
        end if
     end do
     
     core(1,5)=-0.5d0*lso*dsqrt(dfloat(2))
     core(5,1)=core(1,5)
     core(2,6)=core(1,5)
     core(6,2)=core(1,5)
  end if
  
!  inx=1
!  iny=0
!  inz=0
  !  write(6,*)'enter incident light polarization (x,y,z)'                                                    
  !  read(5,*)inx
  !  read(5,*)iny
  !  read(5,*)inz
!  fnx=1
!  fny=0
!  fnz=0
  fup=1
  !  write(6,*)'enter 1 for unpolarized scattered light'
  !  read(5,*)fup
  !  if(fup.ne.1)then
  !     write(6,*)'enter scattered light polarization for RIXS (x,y,z)'                          
  !     read(5,*)fnx
  !     read(5,*)fny
  !     read(5,*)fnz
  !  end if
  !
  !  symm=2
  !  call getpgsymm(symm,pgsymm)
  !
  !  do i=1,5
  !     do j=1,5
  !        if(abs(pgsymm(i,j)).gt.small)write(6,*)i,j,pgsymm(i,j)
  !     end do
  !  end do

  ! 3d TM at center
  xcef=dfloat(0)
  xcef(1,1)=e3/dfloat(2)
  xcef(5,5)=xcef(1,1)
  xcef(1,5)=-xcef(1,1)
  xcef(5,1)=xcef(1,5)
  xcef(3,3)=e2
  xcef(2,2)=(e4+e5)/dfloat(2)
  xcef(4,4)=xcef(2,2)
  xcef(2,4)=(-e4+e5)/dfloat(2)
  xcef(4,2)=xcef(2,4)

  ! 2p ligand at x
  xcef(6,6)=0.5d0*(e6+e7)       !ml=-1
  xcef(7,7)=e8                  !ml=0, pz orbital
  xcef(8,8)=0.5d0*(e6+e7)       !ml=1
  xcef(6,8)=0.5d0*(e6-e7)
  xcef(8,6)=xcef(6,8)
  ! 2p ligand at y
  xcef(9,9)=0.5d0*(e9+e10)      !ml=-1
  xcef(10,10)=e11               !ml=0
  xcef(11,11)=0.5d0*(e9+e10)    !ml=1
  xcef(9,11)=0.5d0*(e9-e10)
  xcef(11,9)=xcef(9,11)
  
  tmatx=dfloat(0)
  tmaty=dfloat(0)
  tmato=dfloat(0)

  tmatx(1,6)=0.5d0*(tpd-tpdxy)           ! (2,-2) <-> (1,-1)
  tmatx(6,1)=tmatx(1,6)
  tmatx(5,6)=0.5d0*(tpd+tpdxy)           ! (2,2)  <-> (1,-1)
  tmatx(6,5)=tmatx(5,6)
  tmatx(1,8)=0.5d0*(tpd+tpdxy)           ! (2,-2) <-> (1,1)
  tmatx(8,1)=tmatx(1,8)
  tmatx(5,8)=0.5d0*(tpd-tpdxy)           ! (2,2)  <-> (1,1)
  tmatx(8,5)=tmatx(5,8)
  tmatx(3,6)=tpdz/dsqrt(dfloat(2))       ! (2,0)  <-> (1,-1)
  tmatx(6,3)=tmatx(3,6)
  tmatx(3,8)=tpdz/dsqrt(dfloat(2))       ! (2,0)  <-> (1,1)
  tmatx(8,3)=tmatx(3,8)
  tmatx(2,7)=-tpdxz/dsqrt(dfloat(2))     ! (2,-1) <-> (1,0)
  tmatx(7,2)=tmatx(2,7)
  tmatx(4,7)=tpdxz/dsqrt(dfloat(2))      ! (2,1)  <-> (1,0)
  tmatx(7,4)=tmatx(4,7)
  
  tmaty(1,9)=-0.5d0*(tpd-tpdxy)          ! (2,-2) <-> (1,-1)
  tmaty(9,1)=tmaty(1,9)
  tmaty(5,9)=-0.5d0*(tpd+tpdxy)          ! (2,2)  <-> (1,-1)
  tmaty(9,5)=tmaty(5,9)
  tmaty(1,11)=0.5d0*(tpd+tpdxy)          ! (2,-2) <-> (1,1)
  tmaty(11,1)=tmaty(1,11)
  tmaty(5,11)=0.5d0*(tpd-tpdxy)          ! (2,2)  <-> (1,1)
  tmaty(11,5)=tmaty(5,11)
  tmaty(3,9)=tpdz/dsqrt(dfloat(2))       ! (2,0)  <-> (1,-1)
  tmaty(9,3)=tmaty(3,9)
  tmaty(3,11)=-tpdz/dsqrt(dfloat(2))     ! (2,0)  <-> (1,1)
  tmaty(11,3)=tmaty(3,11)
  tmaty(2,10)=tpdyz/dsqrt(dfloat(2))     ! (2,-1) <-> (1,0)
  tmaty(10,2)=tmaty(2,10)
  tmaty(4,10)=tpdyz/dsqrt(dfloat(2))     ! (2,1)  <-> (1,0)
  tmaty(10,4)=tmaty(4,10)
  
!  tmato(6,11)=0.5d0*(vpps-vppp)          ! (1,-1) <-> (1,1)
  tmato(6,11)=-vppp          ! (1,-1) <-> (1,1)
  tmato(11,6)=tmato(6,11)
!  tmato(8,9)=-0.5d0*(vpps-vppp)          ! (1,1)  <-< (1,-1)
  tmato(8,9)=vppp          ! (1,1)  <-< (1,-1)
  tmato(9,8)=tmato(8,9)
                                         ! to keep the matrix real [xx<->yx,xy<->yy,xz<->yz] neglected and [xx<->yy, xy<->yx] kept.
  dso=dfloat(0)

  do i=1,5
     dso(i,i)=-0.5d0*dlso*dfloat(i-3)
     dso(i+norbspc,i+norbspc)=0.5d0*dlso*dfloat(i-3)
  end do
  
  dso(1,norbspc+2)=-dlso
  dso(norbspc+2,1)=dso(1,norbspc+2)
  dso(2,norbspc+3)=-dlso*dsqrt(3.0d0)/dsqrt(2.0d0)
  dso(norbspc+3,2)=dso(2,norbspc+3)
  dso(3,norbspc+4)=-dlso*dsqrt(3.0d0)/dsqrt(2.0d0)
  dso(norbspc+4,3)=dso(3,norbspc+4)
  dso(4,norbspc+5)=-dlso
  dso(norbspc+5,4)=dso(4,norbspc+5)
!
! Define the make-up of the unit cell
!
  emax=-100.0d0
  
  if(flag_B.eq.1)then
     !for a square or rectangular cluster
     call setup_s(ncellsx,ncellsy,ncells,momx,momy,mdiff,madd)
  else
     !8 site Betts cluster
     call setup_B(ncells,momx,momy,mdiff,madd)
  end if
  !
  ! Set up is done. Code starts here
  ! First do valence states, then get core hole states
  !
  npairs=nholes*(nholes-1)/2

  if(stepmax.eq.1)then
     allocate(uflag(2*norbs))
     allocate(mflag(2*norbs))
     allocate(sflag(2*norbs))
     allocate(lflag(2*norbs))
     allocate(mlflag(2*norbs))
  end if     
  if(stepmax.eq.2.or.stepmax.eq.3)then
     allocate(uflag(2*norbs+ncore*ncells))
     allocate(mflag(2*norbs+ncore*ncells))
     allocate(sflag(2*norbs+ncore*ncells))
     allocate(lflag(2*norbs+ncore*ncells))
     allocate(mlflag(2*norbs+ncore*ncells))
  end if
  uflag=0
  mflag=0
  sflag=0
  do i=1,ncells
     do j=0,norbspc-1
        nu=1+norbspc*(i-1)+j
        uflag(nu)=1+j
        uflag(nu+norbs)=1+j
        mflag(nu)=i
        mflag(nu+norbs)=i
        sflag(nu)=1
        sflag(nu+norbs)=-1
        if(j.lt.5)then
           lflag(nu)=2
           lflag(nu+norbs)=2
           mlflag(nu)=j-2
           mlflag(nu+norbs)=j-2
        else
           lflag(nu)=1
           lflag(nu+norbs)=1
           if(j.le.7)then
              mlflag(nu)=j-6
              mlflag(nu+norbs)=j-6
           else
              mlflag(nu)=j-9
              mlflag(nu+norbs)=j-9
           end if
        end if
     end do
  end do
  if(stepmax.eq.2.or.stepmax.eq.3)then
     if(edge.eq.1)then
        do j=1,3
           do i=0,ncells-1
              uflag(2*norbs+j+i*3)=j+norbspc
              mflag(2*norbs+j+i*3)=i+1
              sflag(2*norbs+j+i*3)=1
              lflag(2*norbs+j+i*3)=1
              mlflag(2*norbs+j+i*3)=j-2
              uflag(2*norbs+j+i*3+3*ncells)=j+norbspc
              mflag(2*norbs+j+i*3+3*ncells)=i+1
              sflag(2*norbs+j+i*3+3*ncells)=-1
              lflag(2*norbs+j+i*3+3*ncells)=1
              mlflag(2*norbs+j+i*3+3*ncells)=j-2
           end do
        end do
     end if
     if(edge.eq.2)then
        do j=1,2
           do i=0,ncells-1
              uflag(2*norbs+j+i*2)=j+norbspc
              mflag(2*norbs+j+i*2)=i+1
              sflag(2*norbs+j+i*2)=1
              lflag(2*norbs+j+i*2)=0
              mlflag(2*norbs+j+i*2)=0
              uflag(2*norbs+j+i*2+2*ncells)=j+norbspc
              mflag(2*norbs+j+i*2+2*ncells)=i+1
              sflag(2*norbs+j+i*2+2*ncells)=-1
              lflag(2*norbs+j+i*2+2*ncells)=0
              mlflag(2*norbs+j+i*2+2*ncells)=0
           end do
        end do
     end if
  end if

  !this explains the mapping. suppose we have 5 orbitals per unit cell, and we have 2 unit cells.
  !                               the numbering is (2*norbspc+6)*ncells
  !
  !       N=1(k=0)          N=2(k=pi)
  !       --------------    --------------
  !       1,2,3,4,5         6,7,8,9,10          -> the spin up locations
  !       11,12,13,14,15    16,17,18,19,20      -> the spin down locations
  !       21,22,23          24,25,26            -> the core spin up locations
  !       27,28,29          30,31,32            -> the core spin down locations
  !
  !
     
  siz=1
  ndim_max=nint(binomial(2*norbs,nholes))

  !  write(6,*)binomial(2*norbs,nholes),idnint(binomial(2*norbs,nholes))
  
  !  An NxN matrix requires (sizeof(TYPE))*N*N) bytes.
  !  If the available memory is X bytes and double is used [size of(double) is 8 bytes].
  !  Then the largest dimension of N is X = 8 * N*N ==> N = sqrt(X/8).
  !  With 16GB of ram, X=16*(1024^3), so N_max = 46K
  
  bob=ndim_max/1024
  test=8.0*bob*bob/1024
  if(myrank.eq.0)then
     write(6,*)'Hilbert space without core hole -',2*norbs,'choose',nholes,'-:',ndim_max
     if(nint(test).ne.0)write(6,*)'This code requires at least',nint(test),'GB RAM for full ED. Ndim_max=46K for 16 GB'
     write(6,*)
  end if
  !
  !         CODE SET UP ENDS HERE
  !
  !       First step is to set up the basis states for N_holes/N_orbitals and organize them in terms of quantum numbers
  !
  allocate(basis(ndim_max,nholes))
  allocate(mlabel(ndim_max))
  allocate(slabel(ndim_max))
  basis=0
  state=0
  box=0

  call enumerate(nholes,2*norbs,2*norbs,ndim_max,basis,box,state)
  
  dstates=0
  
  do i=1,ndim_max
     count=0
     count2=0
     count4=0
     do j=1,nholes
        ip=mflag(basis(i,j))
        k=sflag(basis(i,j))
        count=count+momx(ip)
        count2=count2+momy(ip)
        count4=count4+k
        if(flag_B.eq.1)then
           if(count.ge.ncellsx)count=count-ncellsx
           if(count2.ge.ncellsy)count2=count2-ncellsy
           if(count.lt.0)count=count+ncellsx
           if(count2.lt.0)count2=count2+ncellsy
        else
           if(count.gt.2)count=count-4
           if(count.lt.-1)count=count+4
           if(count2.gt.2)count2=count2-4
           if(count2.lt.-1)count2=count2+4
        end if
     end do
     do ii=1,ncells
        if(count.eq.momx(ii).and.count2.eq.momy(ii))jj=ii
     end do
     mlabel(i)=jj
     if(dlso.eq.0.0d0)then
        do k=1,nholes+1
           if(count4.eq.-nholes+(k-1)*2)then
              ! the state has been now identified of having total momentum jj and total Sz spin k
              slabel(i)=k
              dstates(jj,k)=dstates(jj,k)+1
           end if
        end do
     else
        dstates(jj,1)=dstates(jj,1)+1
     end if
  end do
  statemax=1
  do i=1,ncells
     if(dlso.eq.0.0d0)then
        do j=1,nholes+1
           if(dstates(i,j).gt.statemax)statemax=dstates(i,j)
           if(myrank.eq.0)write(6,*)i,j,dstates(i,j)
        end do
     else
        if(dstates(i,1).gt.statemax)statemax=dstates(i,1)
     end if
  end do

  if(myrank.eq.0)write(6,*)'States partitioned. Statemax:',statemax

  allocate(binlabel(qnumbers,statemax))
  
  count=0
  if(dlso.eq.0.0d0)then
     do k=1,ncells
        do l=1,nholes+1
           count=count+1
           count2=0
           do n=1,ndim_max
              if(mlabel(n).eq.k.and.slabel(n).eq.l)then
                 count2=count2+1
                 binlabel(count,count2)=n
              end if
           end do
           binstates(count)=count2
           mtot(count)=k
           sztot(count)=-nholes+(l-1)*2
        end do
     end do
  else
     do k=1,ncells
        count=k
        count2=0
        do n=1,ndim_max
           if(mlabel(n).eq.k)then
              count2=count2+1
              binlabel(count,count2)=n
           end if
        end do
        binstates(count)=count2
        mtot(count)=k
     end do
  end if

  qlist=0
  count=0
  do i=1,qnumbers
     if(binstates(i).ne.0)then
        count=count+1
        qlist(count)=i
        if(myrank.eq.0)write(6,*)'qnumbers, number of states:',i,binstates(i)
     end if
  end do

  qstates=count
  if(myrank.eq.0)write(6,*)'total quantum numbers, largest Hilbert subspace:',qstates,statemax
  
  !
  !       DONE SETTING up the basis states for N_holes/N_orbitals
  !
  !
  ! Populate the N_hole Hamiltonian matrix in this basis, partitioned in total momentum sectors
  !
  ! do this for each momentum sector individually since they do not mix
  !
  ! setting up for a sum over twisted boundary conditions
  !
  call CPU_time(time2)

  if(myrank.eq.0)then
     write(16,*)'time initializing',time2-time1
     write(6,*)'time initializing',time2-time1
  end if
  
  vec_px=dfloat(0)
  vec_py=dfloat(0)
  if(myrank.eq.0)then
     open(unit=20,file='Pump_Pulse.dat')
     open(unit=21,file='Occ_Pulse.dat')
  end if
  do i=-nb*timesteps,nt*timesteps
     p_time=pi*dfloat(2*i)/dfloat(timesteps)
     vec_px(i)=pump_ampx*dsin(2.0d0*pi*pump_freq*p_time+0.0d0*pi)*dexp(-0.5d0*p_time*p_time/pump_sigma)
     vec_py(i)=pump_ampy*dsin(2.0d0*pi*pump_freq*p_time+0.0d0*pi)*dexp(-0.5d0*p_time*p_time/pump_sigma)
     if(myrank.eq.0)write(20,*)p_time,vec_px(i),vec_py(i)
  end do
  if(myrank.eq.0)close(20)

  NEV=statemax
  if(usearpack.eq.1)NEV=400
  do twist1=0,twistx-1
     do twist2=0,twisty-1
        !
        ! IF TWISTS ARE NOT EMPLOYED CODE BEGINS HERE
        !
        allocate(degen(qnumbers,NEV))
        allocate(deig(qnumbers,NEV))
        allocate(v(qnumbers,statemax,NEV))
        allocate(edeig(qnumbers,NEV))
        allocate(ev(qnumbers,statemax,NEV))
        allocate(tv(qnumbers,statemax,NEV,-nb*timesteps:nt*timesteps))
        
        do qqq=1,qstates
           ws=qlist(qqq)
           ndim=binstates(ws)
           if(dlso.eq.0.0d0)then
              if(myrank.eq.0)write(6,*)'Sector:         ','K:',momx(mtot(ws)),momy(mtot(ws)),'2*Sz:',sztot(ws),'Hilbert space:',binstates(ws)
           else
              if(myrank.eq.0)write(6,*)'Sector:         ','K:',momx(mtot(ws)),momy(mtot(ws)),'Hilbert space:',binstates(ws)
           end if
!           write(6,*)'twist',twist1,twist2,twistx-1,twisty-1
           if(myrank.eq.0)then
              write(6,*)'------------------------------------------------------------------------------------------------------'
              write(6,*)'writing valence Hamiltonian...'
           end if
           if(usearpack.gt.1)then
              NEV=ndim
              allocate(h(ndim,ndim))
              h=dfloat(0)
              allocate(d_mat(ndim))
              allocate(t1_evol(ndim,ndim))
              allocate(t2_evol(ndim,ndim))
              allocate(evolution(-nb*timesteps:nt*timesteps,ndim,ndim))
              evolution=dcmplx(0.0d0,0.0d0)
           else
              NEV=min(400,ndim/2)
              NCV=2*NEV
              allocate(tIndexI(min(max_array,ndim*ndim)))
              allocate(tIndexJ(min(max_array,ndim*ndim)))
              allocate(tsparseH(min(max_array,ndim*ndim)))
              count=0
           end if

           allocate(dpairs(binstates(ws),npairs,2))
           allocate(dothers(binstates(ws),npairs,nholes-2))
           dpairs=0
           dothers=0

           do n=1,binstates(ws)
              count4=0
              do i=1,nholes-1
                 do j=i+1,nholes
                    count4=count4+1
                    dpairs(n,count4,1)=basis(binlabel(ws,n),i)
                    dpairs(n,count4,2)=basis(binlabel(ws,n),j)
                    if(nholes.gt.2)then
                       count6=0
                       do pp=1,nholes
                          if(pp.ne.i.and.pp.ne.j)then
                             count6=count6+1
                             dothers(n,count4,count6)=basis(binlabel(ws,n),pp)
                          end if
                       end do
                    end if
                 end do
              end do
           end do
              
           do i_time=-nb*timesteps,nt*timesteps
              h=dfloat(0)
              t=dfloat(0)
              xtal=dfloat(0)
              do i=1,2*norbs
                 do j=1,2*norbs
                    if(mflag(i).eq.mflag(j))then
                       k=mflag(i)
                       ii=uflag(i)
                       jj=uflag(j)
                       if(sflag(i)*sflag(j).eq.1)then
                          xtal(i,j)=xcef(ii,jj)
                          if(flag_B.eq.1)then
                             kx=pi*(dfloat(2*momx(k))/dfloat(ncellsx)+dfloat(2*twist1)/dfloat(twistx))+pi
                             t(i,j)=t(i,j)-dfloat(2)*tmatx(ii,jj)*dsin((kx+vec_px(i_time))/dfloat(2))
                             ky=pi*(dfloat(2*momy(k))/dfloat(ncellsy)+dfloat(2*twist2)/dfloat(twisty))+pi
                             t(i,j)=t(i,j)-dfloat(2)*tmaty(ii,jj)*dsin((ky+vec_py(i_time))/dfloat(2))
                             t(i,j)=t(i,j)-dfloat(4)*tmato(ii,jj)*dsin((kx+vec_px(i_time))/2.0d0)*dsin((ky+vec_py(i_time))/2.0d0)
                          else
                             !                       kx=pi*dfloat(4*momx(k))/dfloat(ncells)
                             !                       ky=pi*dfloat(4*momy(k))/dfloat(ncells)                                                                                              
                             !                       t(i,j)=t(i,j)-dfloat(2)*tmaty(ii,jj)*(dcos(kx)+dcos(ky))                                                    
                          end if
                       end if
                       if(sflag(i).eq.1.and.sflag(j).eq.1)xtal(i,j)=xtal(i,j)+dso(ii,jj)
                       if(sflag(i).eq.-1.and.sflag(j).eq.-1)xtal(i,j)=xtal(i,j)+dso(ii+norbspc,jj+norbspc)
                       if(sflag(i).eq.1.and.sflag(j).eq.-1)xtal(i,j)=xtal(i,j)+dso(ii,jj+norbspc)
                       if(sflag(j).eq.1.and.sflag(i).eq.-1)xtal(i,j)=xtal(i,j)+dso(ii+norbspc,jj)
                    end if
                 end do
              end do

              do n=1,binstates(ws)
                 do m=n,binstates(ws)
                    sum=dfloat(0)
                    do i=1,npairs
                       factor1=fsign(i)
                       do j=1,npairs
                          factor2=fsign(j)
                          factor=factor1*factor2
                          prod=1
                          do pp=1,nholes-2
                             if(dothers(n,i,pp).ne.dothers(m,j,pp))prod=0
                          end do
                          if(prod.eq.1)then
                             u1=uflag(dpairs(n,i,1))
                             u2=uflag(dpairs(n,i,2))
                             u3=uflag(dpairs(m,j,1))
                             u4=uflag(dpairs(m,j,2))
                             call checkint(u1,u2,u3,u4,norbspc,norbstm,ncore,outout)
                             if(outout.eq.1)then
                                l1=lflag(dpairs(n,i,1)) !c
                                ii=mlflag(dpairs(n,i,1)) !c
                                l2=lflag(dpairs(n,i,2)) !d
                                jj=mlflag(dpairs(n,i,2)) !d
                                l3=lflag(dpairs(m,j,1)) !a
                                kk=mlflag(dpairs(m,j,1)) !a
                                l4=lflag(dpairs(m,j,2)) !b 
                                ll=mlflag(dpairs(m,j,2)) !b 
                                !                             ! spin a=spin c and spin b = spin d
                                if(sflag(dpairs(n,i,2))*sflag(dpairs(m,j,2)).eq.1)then
                                   if(sflag(dpairs(n,i,1))*sflag(dpairs(m,j,1)).eq.1)then
                                      !                                   ! pa-pc=pd-pb
                                      k=mdiff(mflag(dpairs(n,i,1)),mflag(dpairs(m,j,1)))
                                      l=mdiff(mflag(dpairs(m,j,2)),mflag(dpairs(n,i,2)))
                                      if(k.eq.l)then
                                         if(usearpack.gt.1)then
                                            h(n,m)=h(n,m)+factor*0.5d0*(hint(l1,ii,l2,jj,l4,ll,l3,kk)+hint(l2,jj,l1,ii,l3,kk,l4,ll))/dfloat(ncells)
                                         else
                                            sum=sum+factor*0.5d0*(hint(l1,ii,l2,jj,l4,ll,l3,kk)+hint(l2,jj,l1,ii,l3,kk,l4,ll))/dfloat(ncells)
                                         end if
                                      end if
                                   end if
                                end if
                                !                          ! spin b=spin c and spin a = spin d
                                if(sflag(dpairs(n,i,2))*sflag(dpairs(m,j,1)).eq.1)then
                                   if(sflag(dpairs(n,i,1))*sflag(dpairs(m,j,2)).eq.1)then
                                      !                                   ! pa-pc=pd-pb
                                      k=mdiff(mflag(dpairs(n,i,2)),mflag(dpairs(m,j,1)))
                                      l=mdiff(mflag(dpairs(m,j,2)),mflag(dpairs(n,i,1)))
                                      if(k.eq.l)then
                                         if(usearpack.gt.1)then
                                            h(n,m)=h(n,m)-factor*0.5d0*(hint(l1,ii,l2,jj,l3,kk,l4,ll)+hint(l2,jj,l1,ii,l4,ll,l3,kk))/dfloat(ncells)
                                         else
                                            sum=sum-factor*0.5d0*(hint(l1,ii,l2,jj,l3,kk,l4,ll)+hint(l2,jj,l1,ii,l4,ll,l3,kk))/dfloat(ncells)
                                         end if
                                      end if
                                   end if
                                end if
                             end if
                          end if
                       end do
                    end do

                    do i=1,nholes
                       ii=basis(binlabel(ws,n),i)
                       factor1=(-1)**(i+1)
                       if(nholes.gt.1)then
                          count4=0
                          do ip=1,nholes
                             if(ip.ne.i)then
                                count4=count4+1
                                aa(count4)=basis(binlabel(ws,n),ip)
                             end if
                          end do
                       end if
                       
                       do j=1,nholes
                          jj=basis(binlabel(ws,m),j)
                          factor2=(-1)**(j+1)
                          if(nholes.gt.1)then
                             count4=0
                             do ip=1,nholes
                                if(ip.ne.j)then
                                   count4=count4+1
                                   bb(count4)=basis(binlabel(ws,m),ip)
                                end if
                             end do
                          end if
                          count4=0
                          if(nholes.gt.1)then
                             do k=1,nholes-1
                                do l=1,nholes-1
                                   if(aa(k).eq.bb(l))count4=count4+1
                                end do
                             end do
                          end if
                          if(count4.eq.nholes-1)then
                             factor=factor1*factor2
                             if(mflag(ii).eq.mflag(jj))then
                                if(usearpack.gt.1)then
                                   h(n,m)=h(n,m)+factor*(t(ii,jj)-i4*xtal(ii,jj))
                                else
                                   sum=sum+factor*(t(ii,jj)-i4*xtal(ii,jj))
                                end if
                             end if
                          end if
                       end do
                    end do
                    
                    if(usearpack.le.1)then
                       if(abs(sum).gt.small)then
                          count=count+1
                          tIndexI(count)=n
                          tIndexJ(count)=m
                          tsparseH(count)=sum
                          if(n.ne.m)then
                             count=count+1
                             tIndexI(count)=m
                             tIndexJ(count)=n
                             tsparseH(count)=sum
                          end if
                       end if
                    else
                       h(m,n)=h(n,m)
                    end if
                    
                 end do
              end do

              if(i_time.eq.-nb*timesteps)then
                 if(myrank.eq.0)then
                    if(usearpack.le.1)write(6,*)'count,max_array,ndim*ndim',count,max_array,ndim*ndim
                 end if
              end if
              if(i_time.eq.nt*timesteps)then
                 deallocate(dpairs)
                 deallocate(dothers)
              end if
                 
!              write(6,*)'valence Hamiltonian diagonalizing...'
              if(usearpack.gt.1)then
                 if(i_time.eq.-nb*timesteps)then
                    allocate(w(ndim))
                    iwork=3*ndim-1
                    allocate(work(iwork))
                    degen(ws,1)=1
                    gdegen(ws)=1
                 end if
                 call dsyev(jobz,uplo,ndim,h,ndim,w,work,iwork,info)
                 do i=1,ndim
                    deig(ws,i)=w(i)
                    if(i_time.eq.-nb*timesteps)then
                       if(i.ne.1)then
                          if(abs(sngl(deig(ws,i))-sngl(deig(ws,1))).lt.small)gdegen(ws)=gdegen(ws)+1
                          if(abs(deig(ws,i)-deig(ws,i-1)).lt.small)then
                             degen(ws,i)=degen(ws,i-1)+1
                          else
                             degen(ws,i)=1
                          end if
                       end if
                    end if
                    do j=1,ndim
                       v(ws,j,i)=h(j,i)
                    end do
                 end do
                 if(i_time.eq.nt*timesteps)then
                    deallocate(h)
                    deallocate(w)
                    deallocate(work)
                 end if
              end if
              if(usearpack.le.1)then
                 spr_size=count
                 if(myrank.eq.0)write(6,*)'spr_size:',count
                 allocate(IndexI(spr_size))
                 allocate(IndexJ(spr_size))
                 allocate(sparseH(spr_size))
                 do i=1,spr_size
                    IndexI(i)=tIndexI(i)
                    IndexJ(i)=tIndexJ(i)
                    sparseH(i)=tsparseH(i)
                 end do
                 deallocate(tIndexI)
                 deallocate(tIndexJ)
                 deallocate(tsparseH)
                 allocate(DD(NCV))
                 allocate(ZZ(ndim,NCV))
                 call arpackdiag(spr_size,ndim,NEV,NCV,IndexI,IndexJ,sparseH,DD,ZZ)

                 if(i_time.eq.-nb*timesteps)then
                    degen(ws,1)=1
                    gdegen(ws)=1
                 end if

                 do i=1,NEV
                    deig(ws,i)=DD(i)
                    if(i_time.eq.-nb*timesteps)then
                       if(i.ne.1)then
                          if(abs(sngl(deig(ws,i))-sngl(deig(ws,1))).lt.small)gdegen(ws)=gdegen(ws)+1
                          if(abs(deig(ws,i)-deig(ws,i-1)).lt.small)then
                             degen(ws,i)=degen(ws,i-1)+1
                          else
                             degen(ws,i)=1
                          end if
                       end if
                    end if
                    do j=1,ndim
                       v(ws,j,i)=ZZ(j,i)
                    end do
                 end do

                 if(i_time.eq.nt*timesteps)then
                    deallocate(IndexI)
                    deallocate(IndexJ)
                    deallocate(sparseH)
                    deallocate(DD)
                    deallocate(ZZ)
                 end if
              end if
              
              if(i_time.eq.-nb*timesteps)then
                 if(deig(ws,NEV).gt.emax)emax=deig(ws,NEV)
                 ! for clarity's sake, v(ws,j,i) is the jth eigenvector element of the ith eigenvalue for momentum sector ws                                                           
                 !                                                                                                                                                                       
                 if(myrank.eq.0)then
                    write(6,*)'--------------------------------------------------------------------------------------------------'
                    if(dlso.eq.0.0d0)then
                       write(6,*)'Quantum number, momentum [units pi/Nx,pi/Ny], spin, eigenvalue (lower set of states), degeneracy'
                       write(16,*)'Quantum number, momentum [units pi/Nx,pi/Ny], spin, eigenvalue (lower set of states), degeneracy'
                    else
                       write(6,*)'Quantum number, momentum [units pi/Nx,pi/Ny], eigenvalue (lower set of states), degeneracy'
                       write(16,*)'Quantum number, momentum [units pi/Nx,pi/Ny], eigenvalue (lower set of states), degeneracy'
                    end if
                    write(6,*)'--------------------------------------------------------------------------------------------------'
                    write(16,*)
                    do n=1,min(NEV,40)
                       if(dlso.eq.0.0d0)then
                          write(6,*)ws,momx(mtot(ws)),momy(mtot(ws)),sztot(ws),n,sngl(deig(ws,n)),degen(ws,n)
                          write(16,*)ws,momx(mtot(ws)),momy(mtot(ws)),sztot(ws),n,sngl(deig(ws,n)),degen(ws,n)
                       else
                          write(6,*)ws,momx(mtot(ws)),momy(mtot(ws)),n,sngl(deig(ws,n)),degen(ws,n)
                          write(16,*)ws,momx(mtot(ws)),momy(mtot(ws)),n,sngl(deig(ws,n)),degen(ws,n)
                       end if
                    end do
                    write(6,*)
                    write(16,*)
                 end if
                 do i=1,NEV
                    edeig(ws,i)=deig(ws,i)
                    do j=1,ndim
                       ev(ws,j,i)=v(ws,j,i)
                    end do
                 end do
              end if

              !eigenstate evaluation is now done for a given time step                                                                                                                   
              ! evaluate exp(-i*H(t)dt)=U^T(t) * exp(-i deig(t) dt) * U(t)                                                                                                               
              do i=1,ndim
                 d_mat(i)=cdexp(-pi*eye*dfloat(2)*deig(ws,i)/dfloat(timesteps))
              end do
              do i=1,ndim
                 do j=1,ndim
                    t1_evol(i,j)=d_mat(i)*v(ws,j,i)
                 end do
              end do
              do i=1,ndim
                 do j=1,ndim
                    csum=dcmplx(0.0d0,0.0d0)
                    do k=1,ndim
                       csum=csum+v(ws,i,k)*t1_evol(k,j)
                    end do
                    t2_evol(i,j)=csum
                 end do
              end do
              do i=1,ndim
                 do j=1,ndim
                    evolution(i_time,i,j)=t2_evol(i,j)
                 end do
              end do
              if(i_time.ne.-nb*timesteps)then
                 do n=1,ndim
                    do m=1,ndim
                       csum=dcmplx(0.0d0,0.0d0)
                       do j=1,ndim
                          csum=csum+evolution(i_time,n,j)*tv(ws,j,m,i_time-1)
                       end do
                       tv(ws,n,m,i_time)=csum
                    end do
                 end do
              else
                 do n=1,ndim
                    do m=1,ndim
                       tv(ws,n,m,i_time)=ev(ws,n,m)
                    end do
                 end do
              end if
              deig=dfloat(0)
              v=dfloat(0)
           end do !i_time loop      
           deallocate(d_mat)
           deallocate(t1_evol)
           deallocate(t2_evol)
           deallocate(evolution)
        end do
        
        !This completes the sum over momentum sectors for the states without the core hole
        !
        ! measure some things in the ground state at t=-infinity
        !
        egs=1000000.0d0
        if(myrank.eq.0)then
           write(6,*)'ground state energies of each sector:'
           write(16,*)'ground state energies of each sector:'
        end if
        i=0
        do qqq=1,qstates
           ws=qlist(qqq)
           i=i+1
           if(myrank.eq.0)then
              write(6,*)ws,sngl(edeig(ws,1))
              write(16,*)ws,sngl(edeig(ws,1))
           end if
           a(i)=sngl(edeig(ws,1))
        end do
        if(myrank.eq.0)then
           write(6,*)
           write(16,*)
        end if
        
        call quicksort(a, 1, i)
        count=0
        gsec=0
        do qqq=1,qstates
           ws=qlist(qqq)
           if(sngl(edeig(ws,1)).eq.a(1))then
              if(myrank.eq.0)then
                 write(6,*)ws,'belongs to the ground state sector and has degeneracy',gdegen(ws),'ground state energy:',sngl(a(1)),'momentum:',momx(mtot(ws)),momy(mtot(ws)),'2*S_z:',sztot(ws)
                 write(16,*)ws,'belongs to the ground state sector and has degeneracy',gdegen(ws),'ground state energy:',sngl(a(1)),'momentum:',momx(mtot(ws)),momy(mtot(ws)),'2*S_z:',sztot(ws)
              end if
              if(count.eq.0)then
                 momgs=mtot(ws)
              else
                 if(mtot(ws).ne.momgs)then
                    write(6,*)'Ground state is degenerate for more than 1 total momentum. Code cannot proceed.'
                    goto 999
                 end if
              end if
              count=count+gdegen(ws)
              gsec(ws)=1
              !determine the ground state hole occupations
              if(myrank.eq.0)write(6,*)'Ground state occupations'

              do wsp=1,gdegen(ws)
                 occ=dfloat(0)
                 do n=1,binstates(ws)
                    do i=1,nholes
                       ii=basis(binlabel(ws,n),i)
                       factor1=(-1)**(i+1)
                       if(nholes.gt.1)then
                          count4=0
                          do ip=1,nholes
                             if(ip.ne.i)then
                                count4=count4+1
                                aa(count4)=basis(binlabel(ws,n),ip)
                             end if
                          end do
                       end if
                       
                       do m=1,binstates(ws)
                          do j=1,nholes
                             jj=basis(binlabel(ws,m),j)
                             factor2=(-1)**(j+1)
                             if(nholes.gt.1)then
                                count4=0
                                do ip=1,nholes
                                   if(ip.ne.j)then
                                      count4=count4+1
                                      bb(count4)=basis(binlabel(ws,m),ip)
                                   end if
                                end do
                             end if
                             count4=0
                             if(nholes.gt.1)then
                                do k=1,nholes-1
                                   do l=1,nholes-1
                                      if(aa(k).eq.bb(l))count4=count4+1
                                   end do
                                end do
                             end if
                             if(count4.eq.nholes-1)then
                                factor=factor1*factor2
                                if(mflag(ii).eq.mflag(jj))then
                                   occ(uflag(ii),uflag(jj))=occ(uflag(ii),uflag(jj))+factor*ev(ws,n,wsp)*ev(ws,m,wsp)
                                end if
                             end if
                          end do
                       end do
                    end do
                 end do
                 if(myrank.eq.0)then
                    write(6,*)'x^2-y^2 occupation:',(occ(1,1)+occ(5,5)+occ(1,5)+occ(5,1))/dfloat(2)
                    write(6,*)'z^2 occupation:',occ(3,3)
                    write(6,*)'xy occupation:',(occ(1,1)+occ(5,5)-occ(1,5)-occ(5,1))/dfloat(2)
                    write(6,*)'xz occupation:',(occ(2,2)+occ(4,4)+occ(2,4)+occ(4,2))/dfloat(2)
                    write(6,*)'yz occupation:',(occ(2,2)+occ(4,4)-occ(2,4)-occ(4,2))/dfloat(2)
                    write(16,*)'x^2-y^2 occupation:',(occ(1,1)+occ(5,5)+occ(1,5)+occ(5,1))/dfloat(2)
                    write(16,*)'z^2 occupation:',occ(3,3)
                    write(16,*)'xy occupation:',(occ(1,1)+occ(5,5)-occ(1,5)-occ(5,1))/dfloat(2)
                    write(16,*)'xz occupation:',(occ(2,2)+occ(4,4)+occ(2,4)+occ(4,2))/dfloat(2)
                    write(16,*)'yz occupation:',(occ(2,2)+occ(4,4)-occ(2,4)-occ(4,2))/dfloat(2)
                    if(norbsl.ne.0)then
                       write(6,*)'ligand occupation:'
                       write(6,*)'xx occupation:',(occ(6,6)+occ(8,8)+occ(6,8)+occ(8,6))/dfloat(2)
                       write(6,*)'xy occupation:',(occ(6,6)+occ(8,8)-occ(6,8)-occ(8,6))/dfloat(2)
                       write(6,*)'xz occupation:',occ(7,7)
                       write(6,*)'yx occupation:',(occ(9,9)+occ(11,11)+occ(9,11)+occ(11,9))/dfloat(2)
                       write(6,*)'yy occupation:',(occ(9,9)+occ(11,11)-occ(9,11)-occ(11,9))/dfloat(2)
                       write(6,*)'yz occupation:',occ(10,10)
                       write(16,*)'ligand occupation:'
                       write(16,*)'xx occupation:',(occ(6,6)+occ(8,8)+occ(6,8)+occ(8,6))/dfloat(2)
                       write(16,*)'xy occupation:',(occ(6,6)+occ(8,8)-occ(6,8)-occ(8,6))/dfloat(2)
                       write(16,*)'xz occupation:',occ(7,7)
                       write(16,*)'yx occupation:',(occ(9,9)+occ(11,11)+occ(9,11)+occ(11,9))/dfloat(2)
                       write(16,*)'yy occupation:',(occ(9,9)+occ(11,11)-occ(9,11)-occ(11,9))/dfloat(2)
                       write(16,*)'yz occupation:',occ(10,10)
                    end if
                    write(6,*)
                    write(16,*)
                 end if
              end do
           end if
        end do
        
        gstatedegen=count
        if(myrank.eq.0)then
           write(6,*)'ground state is',gstatedegen,'-fold degenerate'
           write(16,*)'ground state is',gstatedegen,'-fold degenerate'
        end if
        gamma=0.01d0

        if(myrank.eq.0)then
           if(twist1.eq.0.and.twist2.eq.0)open(unit=55,file='DOS.dat')
        end if
        if(usearpack.gt.1)then
           do i=1,nfreq
              omega=-1.0d0+a(1)+(emax+2.0d0-a(1))*dfloat(i)/dfloat(nfreq)
              sum=0.0d0
              do qqq=1,qstates
                 ws=qlist(qqq)
                 !              do j=1,min(NEV,binstates(ws))
                 do j=1,binstates(ws)
                    sum=sum+(1.0d0/pi)*gamma/(gamma*gamma+(omega-edeig(ws,j))**2)/dfloat(ncells*twistx*twisty)
                 end do
              end do
              if(myrank.eq.0)write(55,*)omega,sum
           end do
        else
           do i=1,nfreq
              omega=-1.0d0+a(1)+(emax+2.0d0-a(1))*dfloat(i)/dfloat(nfreq)
              sum=0.0d0
              do qqq=1,qstates
                 ws=qlist(qqq)
                 do j=1,NEV
                    sum=sum+(1.0d0/pi)*gamma/(gamma*gamma+(omega-edeig(ws,j))**2)/dfloat(ncells*twistx*twisty)
                 end do
              end do
              if(myrank.eq.0)write(55,*)omega,sum
           end do
        end if
        if(myrank.eq.0)then
           close(55)
           write(6,*)
        end if
        
        if(ncells.eq.1)then
           if(abs(i4).lt.small)then
              if(nholes.gt.1.and.nholes.le.5)then
                 write(6,*)'Checks for',nholes,'holes without crystal fields'
                 call checks(nholes,raca,racb,racc,nchecks)     
                 if(nholes.eq.2)then
                    write(6,*)(nchecks(i),i=1,5)
                 end if
                 if(nholes.eq.3)then
                    write(6,*)(nchecks(i),i=1,7)
                 end if
                 if(nholes.gt.3)then
                    do i=1,16
                       write(6,*)i,nchecks(i)
                    end do
                 end if
              end if
           end if
        end if
        
        !write out momentum dependent eigenstates
        if(flag_w.eq.1)then
           open(unit=8,file='eigenvalues.dat')
           do i=0,ncellsx-1
              do qqq=1,qstates
                 ws=qlist(qqq)
                 if(momx(mtot(ws)).eq.i)then
                    kx=pi*dfloat(2*momx(mtot(ws)))/dfloat(ncellsx)
                    do j=0,ncellsy-1
                       if(momy(mtot(ws)).eq.j)then
                          ky=pi*dfloat(2*momy(mtot(ws)))/dfloat(ncellsy)
                          do n=1,binstates(ws)
                             if(myrank.eq.0)write(8,*)kx,ky,(deig(ws,k),k=1,binstates(ws))
                          end do
                       end if
                    end do
                 end if
              end do
              if(myrank.eq.0)write(8,*)
           end do
           if(myrank.eq.0)close(8)
        end if

        call CPU_time(time3)
        
        if(myrank.eq.0)then
           write(6,*)'-----------------------------------------------------'
           write(6,*)'-------End of non-core hole calculation--------------'
           write(6,*)'-----------------------------------------------------'
           write(6,*)
           write(6,*)'time writing and diagonalizing valence Hamiltonian matrix',time3-time2
           write(16,*)'time writing and diagonalizing valence Hamiltonian matrix',time3-time2
        end if
        !calculating time dependent observables                                                                                                                                          
        !                                                                                                                                      
        do qqq=1,qstates
           ws=qlist(qqq)
           if(gsec(ws).eq.1)then    !ground state                                                                                                                                        
              allocate(u11(binstates(ws)))
              do wsp=1,gdegen(ws)
                 do ti=-nb*timesteps,nt*timesteps
                    p_time=pi*dfloat(2*ti)/dfloat(timesteps)
                    tj=-nb*timesteps
                    do i=1,binstates(ws)
                       u11(i)=tv(ws,i,wsp,ti)
                    end do
                    b_occ=dcmplx(0.0d0,0.0d0)
                    do n=1,binstates(ws)
                       do i=1,nholes
                          ii=basis(binlabel(ws,n),i)
                          factor1=(-1)**(i+1)
                          if(nholes.gt.1)then
                             count4=0
                             do ip=1,nholes
                                if(ip.ne.i)then
                                   count4=count4+1
                                   aa(count4)=basis(binlabel(ws,n),ip)
                                end if
                             end do
                          end if
                          do m=1,binstates(ws)
                             do j=1,nholes
                                jj=basis(binlabel(ws,m),j)
                                factor2=(-1)**(j+1)
                                if(nholes.gt.1)then
                                   count4=0
                                   do ip=1,nholes
                                      if(ip.ne.j)then
                                         count4=count4+1
                                         bb(count4)=basis(binlabel(ws,m),ip)
                                      end if
                                   end do
                                end if
                                count4=0
                                if(nholes.gt.1)then
                                   do k=1,nholes-1
                                      do l=1,nholes-1
                                         if(aa(k).eq.bb(l))count4=count4+1
                                      end do
                                   end do
                                end if
                                if(count4.eq.nholes-1)then
                                   factor=factor1*factor2
                                   if(mflag(ii).eq.mflag(jj))then
                                      b_occ(uflag(ii),uflag(jj))=b_occ(uflag(ii),uflag(jj))+factor*dconjg(u11(m))*u11(n)
                                   end if
                                end if
                             end do
                          end do
                       end do
                    end do
                    o1=cdabs((b_occ(1,1)+b_occ(5,5)+b_occ(1,5)+b_occ(5,1))/dfloat(2))
                    o2=cdabs(b_occ(3,3))
                    o3=cdabs((b_occ(1,1)+b_occ(5,5)-b_occ(1,5)-b_occ(5,1))/dfloat(2))
                    o4=cdabs((b_occ(2,2)+b_occ(4,4)+b_occ(2,4)+b_occ(4,2))/dfloat(2))
                    o5=cdabs((b_occ(2,2)+b_occ(4,4)-b_occ(2,4)-b_occ(4,2))/dfloat(2))
                    if(norbsl.ne.0)then
                       o6=cdabs((b_occ(6,6)+b_occ(8,8)+b_occ(6,8)+b_occ(8,6))/dfloat(2))
                       o7=cdabs((b_occ(6,6)+b_occ(8,8)-b_occ(6,8)-b_occ(8,6))/dfloat(2))
                       o8=cdabs(b_occ(7,7))
                       o9=cdabs((b_occ(9,9)+b_occ(11,11)+b_occ(9,11)+b_occ(11,9))/dfloat(2))
                       o10=cdabs((b_occ(9,9)+b_occ(11,11)-b_occ(9,11)-b_occ(11,9))/dfloat(2))
                       o11=cdabs(b_occ(10,10))
                    end if
                    if(norbsl.ne.0)then
                       if(myrank.eq.0)write(21,*)p_time,o1,o2,o3,o4,o5,o6,o7,o8,o9,o10,o11
                    else
                       if(myrank.eq.0)write(21,*)p_time,o1,o2,o3,o4,o5
                    end if
                 end do
              end do
              deallocate(u11)
           end if
        end do
        close(21)

!        do qqq=1,qstates
!           ws=qlist(qqq)
!           do i_time=-nb*timesteps,nt*timesteps
!              do n=1,binstates(ws)
!                 csum=dcmplx(0.0d0,0.0d0)
!                 do m=1,binstates(ws)
!                    csum=csum+ev(ws,m,n)*tv(ws,m,n,i_time)
!                 end do
!                 write(6,*)csum,cdexp(-eye*pi*edeig(ws,n)*dfloat(2*(i_time+nb*timesteps))/dfloat(timesteps))
!              end do
!           end do
!        end do
!        goto 666
        !
        ! NOTE: |psi_m(t)> = tv(ws,n,m,t) is the eigenvector component n for eigenstate m at time t 
        !
        ! the check above shows that <psi_m(t)|psi_m(t=t_min)>=exp(-i(t-t_min)E_m)
        !
        if(stepmax.eq.2)then
           !
           !     NOW DOING CORE HOLE CALCULATION - SET UP
           !     spin orbit is included so only quantum numbers are set by total momentum
           !     since the core hole can have any momentum at the same energy it is sufficient to calculate...
           !
           if(twist1.eq.0.and.twist2.eq.0)then
              ndimu=nint(binomial(2*norbs,nholes-1))
              cdim_max=ncore*ncells*ndimu
              if(myrank.eq.0)then
                 write(6,*)'edge (1=TM L, 2=Ligand K):',edge
                 write(6,*)'Hilbert space for core hole states:',cdim_max  
              end if
                 
              allocate(cbasis(cdim_max,nholes))
              allocate(cmlabel(cdim_max))
              allocate(vect(ndimu,nholes-1))
              cbasis=0
              state=0
              box=0
              
              call enumerate(nholes-1,2*norbs,2*norbs,ndimu,vect,box,state)
              
              cstates=0
              state=1
              
              do i=1,ndimu
                 do l=1,ncore*ncells
                    do k=1,nholes-1
                       cbasis(state,k)=vect(i,k)
                    end do
                    cbasis(state,nholes)=2*norbs+l
                    state=state+1
                 end do
              end do
              
              cstates=0
              
              if(myrank.eq.0)write(6,*)'Core states partitioned'
              do i=1,cdim_max
                 count=0
                 count2=0
                 do j=1,nholes
                    ip=mflag(cbasis(i,j))
                    count=count+momx(ip)
                    count2=count2+momy(ip)
                    if(uflag(cbasis(i,j)).le.2*norbs)count4=count4+sflag(cbasis(i,j))
                    if(flag_B.eq.1)then
                       if(count.ge.ncellsx)count=count-ncellsx
                       if(count2.ge.ncellsy)count2=count2-ncellsy
                       if(count.lt.0)count=count+ncellsx
                       if(count2.lt.0)count2=count2+ncellsy
                    else
                       if(count.gt.2)count=count-4
                       if(count.lt.-1)count=count+4
                       if(count2.gt.2)count2=count2-4
                       if(count2.lt.-1)count2=count2+4
                    end if
                 end do
                 do ii=1,ncells
                    if(count.eq.momx(ii).and.count2.eq.momy(ii))jj=ii
                 end do
                 ! the state has been now identified of having total momentum jj
                 cstates(jj)=cstates(jj)+1
                 cmlabel(i)=jj
              end do
              
              allocate(cbinlabel(cstates(momgs)))
              allocate(cpairs(cstates(momgs),npairs,2))
              allocate(cothers(cstates(momgs),npairs,nholes-2))
              cpairs=0
              cothers=0
              
              count=0
              do n=1,cdim_max
                 if(cmlabel(n).eq.momgs)then
                    count=count+1
                    count4=0
                    cbinlabel(count)=n
                    do i=1,nholes-1
                       do j=i+1,nholes
                          count4=count4+1
                          cpairs(count,count4,1)=cbasis(n,i)
                          cpairs(count,count4,2)=cbasis(n,j)
                          if(nholes.gt.2)then
                             count6=0
                             do pp=1,nholes
                                if(pp.ne.i.and.pp.ne.j)then
                                   count6=count6+1
                                   cothers(count,count4,count6)=cbasis(n,pp)
                                end if
                             end do
                          end if
                       end do
                    end do
                 end if
              end do
              
              !diagonalize the (N_holes-1) + (1) core-hole matrix
              if(myrank.eq.0)write(6,*)'constructing core hole Hamiltonian...'
              ndim=cstates(momgs)              
           end if

           if(usearpack.ge.1.and.usearpack.le.2)then
              NEVN=min(800,cstates(momgs)/2)
              NCV=2*NEVN
              allocate(dn(NEVN))
              allocate(vn(cstates(momgs),NEVN))
              allocate(tIndexI(10000000))
              allocate(tIndexJ(10000000))
              allocate(tsparseH(10000000))
              count=0
              allocate(d_mat(ndim))
              allocate(t1_evol(ndim,ndim))
              allocate(t2_evol(ndim,ndim))
              allocate(evolution(-nb*timesteps:nt*timesteps,ndim,ndim))
           else
              NEVN=cstates(momgs)
              allocate(dn(cstates(momgs)))
              allocate(vn(cstates(momgs),cstates(momgs)))
              allocate(h(cstates(momgs),cstates(momgs)))
              h=dfloat(0)
              allocate(d_mat(ndim))
              allocate(t1_evol(ndim,ndim))
              allocate(t2_evol(ndim,ndim))
              allocate(evolution(-nb*timesteps:nt*timesteps,ndim,ndim))
           end if
           
           allocate(edn(NEVN))
           allocate(evn(cstates(momgs),NEVN))
           allocate(tvn(cstates(momgs),NEVN,-nb*timesteps:nt*timesteps))
           
           do i_time=-nb*timesteps,nt*timesteps
              h=dfloat(0)
              t=dfloat(0)
              xtal=dfloat(0)
              do i=1,2*norbs
                 do j=1,2*norbs
                    if(mflag(i).eq.mflag(j))then
                       k=mflag(i)
                       ii=uflag(i)
                       jj=uflag(j)
                       if(sflag(i)*sflag(j).eq.1)then
                          xtal(i,j)=xcef(ii,jj)
                          if(flag_B.eq.1)then
                             kx=pi*(dfloat(2*momx(k))/dfloat(ncellsx)+dfloat(2*twist1)/dfloat(twistx))+pi
                             t(i,j)=t(i,j)+dfloat(2)*tmatx(ii,jj)*dsin((kx+vec_px(i_time))/dfloat(2))
                             ky=pi*(dfloat(2*momy(k))/dfloat(ncellsy)+dfloat(2*twist2)/dfloat(twisty))+pi
                             t(i,j)=t(i,j)+dfloat(2)*tmaty(ii,jj)*dsin((ky+vec_py(i_time))/dfloat(2))
                             t(i,j)=t(i,j)+dfloat(4)*tmato(ii,jj)*dsin((kx+vec_px(i_time))/2.0d0)*dsin((ky+vec_py(i_time))/2.0d0)
                          else
                             !                       kx=pi*dfloat(4*momx(k))/dfloat(ncells)                                               
                             !                       ky=pi*dfloat(4*momy(k))/dfloat(ncells)                                               
                             !                       t(i,j)=t(i,j)-dfloat(2)*tmaty(ii,jj)*(dcos(kx)+dcos(ky))                             
                          end if
                       end if
                       if(sflag(i).eq.1.and.sflag(j).eq.1)xtal(i,j)=xtal(i,j)+dso(ii,jj)
                       if(sflag(i).eq.-1.and.sflag(j).eq.-1)xtal(i,j)=xtal(i,j)+dso(ii+norbspc,jj+norbspc)
                       if(sflag(i).eq.1.and.sflag(j).eq.-1)xtal(i,j)=xtal(i,j)+dso(ii,jj+norbspc)
                       if(sflag(j).eq.1.and.sflag(i).eq.-1)xtal(i,j)=xtal(i,j)+dso(ii+norbspc,jj)
                    end if
                 end do
              end do
              
              do n=1,cstates(momgs)
                 do m=n,cstates(momgs)
                    sum=dfloat(0)
                    do i=1,nholes
                       ii=cbasis(cbinlabel(n),i)
                       factor1=(-1)**(i+1)
                       if(nholes.gt.1)then
                          count4=0
                          do ip=1,nholes
                             if(ip.ne.i)then
                                count4=count4+1
                                aa(count4)=cbasis(cbinlabel(n),ip)
                             end if
                          end do
                       end if
                       
                       do j=1,nholes
                          jj=cbasis(cbinlabel(m),j)
                          factor2=(-1)**(j+1)
                          if(nholes.gt.1)then
                             count4=0
                             do ip=1,nholes
                                if(ip.ne.j)then
                                   count4=count4+1
                                   bb(count4)=cbasis(cbinlabel(m),ip)
                                end if
                             end do
                          end if
                          count4=0
                          if(nholes.gt.1)then
                             do k=1,nholes-1
                                do l=1,nholes-1
                                   if(aa(k).eq.bb(l))count4=count4+1
                                end do
                             end do
                          end if
                          if(count4.eq.nholes-1)then
                             factor=factor1*factor2
                             if(mflag(ii).eq.mflag(jj))then
                                if(uflag(ii).le.norbspc.and.uflag(jj).le.norbspc)then
                                   if(usearpack.ge.1.and.usearpack.le.2)then
                                      sum=sum+factor*(t(ii,jj)-i4*xtal(ii,jj))
                                   else
                                      h(n,m)=h(n,m)+factor*(t(ii,jj)-i4*xtal(ii,jj))
                                   end if
                                end if
                             end if
                             if(uflag(ii).gt.norbspc.and.uflag(jj).gt.norbspc)then
                                k=uflag(ii)-norbspc
                                l=uflag(jj)-norbspc
                                if(sflag(ii).eq.-1)k=k+3
                                if(sflag(jj).eq.-1)l=l+3
                                if(usearpack.ge.1.and.usearpack.le.2)then
                                   sum=sum+core(k,l)*factor
                                else
                                   h(n,m)=h(n,m)+core(k,l)*factor
                                end if
                             end if
                          end if
                       end do
                    end do
                 
                    do i=1,npairs
                       factor1=fsign(i)
                       do j=1,npairs
                          factor2=fsign(j)
                          factor=factor1*factor2
                          prod=1
                          do pp=1,nholes-2
                             if(cothers(n,i,pp).ne.cothers(m,j,pp))prod=0
                          end do
                          if(prod.eq.1)then
                             u1=uflag(cpairs(n,i,1)) !c
                             l1=lflag(cpairs(n,i,1)) !c
                             ii=mlflag(cpairs(n,i,1)) !c
                             u2=uflag(cpairs(n,i,2)) !d
                             l2=lflag(cpairs(n,i,2)) !d
                             jj=mlflag(cpairs(n,i,2)) !d
                             u3=uflag(cpairs(m,j,1)) !a
                             l3=lflag(cpairs(m,j,1)) !a
                             kk=mlflag(cpairs(m,j,1)) !a
                             u4=uflag(cpairs(m,j,2)) !b 
                             l4=lflag(cpairs(m,j,2)) !b 
                             ll=mlflag(cpairs(m,j,2)) !b 
                             call checkint(u1,u2,u3,u4,norbspc,norbstm,ncore,outout)
                             if(outout.eq.1)then
                                !                             ! spin a=spin c and spin b = spin d
                                if(sflag(cpairs(n,i,2))*sflag(cpairs(m,j,2)).eq.1)then
                                   if(sflag(cpairs(n,i,1))*sflag(cpairs(m,j,1)).eq.1)then
                                      !                                   ! pa-pc=pd-pb
                                      k=mdiff(mflag(cpairs(n,i,1)),mflag(cpairs(m,j,1)))
                                      l=mdiff(mflag(cpairs(m,j,2)),mflag(cpairs(n,i,2)))
                                      if(k.eq.l)then
                                         if(usearpack.ge.1.and.usearpack.le.2)then
                                            sum=sum+factor*0.5d0*(hint(l1,ii,l2,jj,l4,ll,l3,kk)+hint(l2,jj,l1,ii,l3,kk,l4,ll))/dfloat(ncells)
                                         else
                                            h(n,m)=h(n,m)+factor*0.5d0*(hint(l1,ii,l2,jj,l4,ll,l3,kk)+hint(l2,jj,l1,ii,l3,kk,l4,ll))/dfloat(ncells)
                                         end if
                                      end if
                                   end if
                                end if
                                !                          ! spin b=spin c and spin a = spin d
                                if(sflag(cpairs(n,i,2))*sflag(cpairs(m,j,1)).eq.1)then
                                   if(sflag(cpairs(n,i,1))*sflag(cpairs(m,j,2)).eq.1)then
                                      !                                   ! pa-pc=pd-pb
                                      k=mdiff(mflag(cpairs(n,i,2)),mflag(cpairs(m,j,1)))
                                      l=mdiff(mflag(cpairs(m,j,2)),mflag(cpairs(n,i,1)))
                                      if(k.eq.l)then
                                         if(usearpack.ge.1.and.usearpack.le.2)then
                                            sum=sum-factor*0.5d0*(hint(l1,ii,l2,jj,l3,kk,l4,ll)+hint(l2,jj,l1,ii,l4,ll,l3,kk))/dfloat(ncells)
                                         else
                                            h(n,m)=h(n,m)-factor*0.5d0*(hint(l1,ii,l2,jj,l3,kk,l4,ll)+hint(l2,jj,l1,ii,l4,ll,l3,kk))/dfloat(ncells)
                                         end if
                                      end if
                                   end if
                                end if
                             end if
                          end if
                       end do
                    end do
                    if(usearpack.ge.1.and.usearpack.le.2)then
                       if(abs(sum).gt.small)then
                          count=count+1
                          tIndexI(count)=n
                          tIndexJ(count)=m
                          tsparseH(count)=sum
                          if(n.ne.m)then
                             count=count+1
                             tIndexI(count)=m
                             tIndexJ(count)=n
                             tsparseH(count)=sum
                          end if
                       end if
                    else
                       h(m,n)=h(n,m)
                    end if
                 end do
              end do

              if(myrank.eq.0)then
                 if(i_time.eq.-nb*timesteps)write(6,*)'Core Hamiltonian diagonalizing...twist1,twist2',twist1,twist2
              end if
              
              if(usearpack.ge.1.and.usearpack.le.2)then
                 spr_size=count
                 if(myrank.eq.0)write(6,*)'spr_size:',count
                 allocate(IndexI(spr_size))
                 allocate(IndexJ(spr_size))
                 allocate(sparseH(spr_size))
                 do i=1,spr_size
                    IndexI(i)=tIndexI(i)
                    IndexJ(i)=tIndexJ(i)
                    sparseH(i)=tsparseH(i)
                 end do
                 deallocate(tIndexI)
                 deallocate(tIndexJ)
                 deallocate(tsparseH)
                 allocate(DD(NCV))
                 allocate(ZZ(cstates(momgs),NCV))
                 call arpackdiag(spr_size,cstates(momgs),NEVN,NCV,IndexI,IndexJ,sparseH,DD,ZZ)
              
                 do i=1,NEVN
                    dn(i)=DD(i)
                    do j=1,cstates(momgs)
                       vn(j,i)=ZZ(j,i)
                    end do
                 end do

                 if(i_time.eq.nt*timesteps)then
                    deallocate(IndexI)
                    deallocate(IndexJ)
                    deallocate(sparseH)
                    deallocate(DD)
                    deallocate(ZZ)
                 end if
              else
                 if(i_time.eq.-nb*timesteps)then
                    allocate(w(cstates(momgs)))
                    allocate(work(3*cstates(momgs)-1))
                 end if
                 call dsyev(jobz,uplo,cstates(momgs),h,cstates(momgs),w,work,3*cstates(momgs)-1,info)
                 do i=1,cstates(momgs)
                    dn(i)=w(i)
                    do j=1,cstates(momgs)
                       vn(j,i)=h(j,i)
                    end do
                 end do
                 if(i_time.eq.-nb*timesteps)then
                    do i=1,NEVN
                       edn(i)=w(i)
                       do j=1,cstates(momgs)
                          evn(j,i)=h(j,i)
                       end do
                    end do
                 end if
                 if(i_time.eq.nt*timesteps)then
                    deallocate(h)
                    deallocate(w)
                    deallocate(work)
                 end if
              end if
           
           !
           ! for clarity's sake, v(ws,j,i) is the jth eigenvector element of the ith eigenvalue for momentum sector ws
           !
           if(i_time.eq.-nb*timesteps)then
              if(twist1.eq.0.and.twist2.eq.0)then
                 if(myrank.eq.0)open(unit=57,file='cDOS.dat')
              end if
              do i=1,nfreq
                 omega=-1.0d0+dn(1)+dfloat(i)*(2.0d0-dn(1)+dn(NEVN))/dfloat(nfreq)
                 sum=0.0d0
                 do j=1,min(NEVN,cstates(momgs))
                    sum=sum+(1.0d0/pi)*gamma/(gamma*gamma+(omega-dn(j))**2)/dfloat(ncells*twistx*twisty)
                 end do
                 if(myrank.eq.0)write(57,*)omega,sum
              end do
              if(myrank.eq.0)close(57)
           end if
           call CPU_time(time4)  
           if(myrank.eq.0)then
              if(i_time.eq.nt*timesteps)then
                 write(6,*)'-----------------------------------------------------'
                 write(6,*)'-------End of core hole calculation--------------'
                 write(6,*)'-----------------------------------------------------'
                 write(6,*)'time writing and diagonalizing core Hamiltonian matrix',time4-time3
                 write(16,*)'time writing and diagonalizing core Hamiltonian matrix',time4-time3
              end if
           end if

           !eigenstate evaluation is now done for a given time step                                                                    

           ! evaluate exp(-i*H(t)dt)=U^T(t) * exp(-i deig(t) dt) * U(t)                                                                
           do i=1,ndim
              d_mat(i)=cdexp(-dfloat(2)*pi*eye*dn(i)/dfloat(timesteps))
           end do
           do i=1,ndim
              do j=1,ndim
                 t1_evol(i,j)=d_mat(i)*vn(j,i)
              end do
           end do
           do i=1,ndim
              do j=1,ndim
                 csum=dcmplx(0.0d0,0.0d0)
                 do k=1,ndim
                    csum=csum+vn(i,k)*t1_evol(k,j)
                 end do
                 t2_evol(i,j)=csum
              end do
           end do
           do i=1,ndim
              do j=1,ndim
                 evolution(i_time,i,j)=t2_evol(i,j)
              end do
           end do
           if(i_time.ne.-nb*timesteps)then
              do n=1,ndim
                 do m=1,ndim
                    csum=dcmplx(0.0d0,0.0d0)
                    do j=1,ndim
                       csum=csum+evolution(i_time,n,j)*tvn(j,m,i_time-1)
                    end do
                    tvn(n,m,i_time)=csum
                 end do
              end do
           else
              do n=1,ndim
                 do m=1,ndim
                    tvn(n,m,i_time)=evn(n,m)
                 end do
              end do
           end if
           dn=dfloat(0)
           vn=dfloat(0)

        end do   !time loop i_time                                                                                                     
        
        deallocate(d_mat)
        deallocate(t1_evol)
        deallocate(t2_evol)
!        deallocate(evolution)
        call CPU_time(time5)
     end if
     
     if(stepmax.eq.3)then
           !
           !     NOW DOING HOLE ADDITION CALCULATION
           !
        hnpairs=nholes*(nholes+1)/2

        if(twist1.eq.0.and.twist2.eq.0)then
           ndimu=nint(binomial(2*norbs,nholes+1))
           ndimd=nint(binomial(2*norbs,nholes))
           hdim_max=ncore*ncells*ndimd+ndimu
           write(6,*)'edge (1=TM L, 2=Ligand K):',edge
           write(6,*)'Hilbert space for hole addition states:',hdim_max  
           
           allocate(hbasis(hdim_max,nholes+1))
           allocate(hmlabel(hdim_max))
           allocate(hslabel(hdim_max))
           allocate(vect(ndimu,nholes+1))
           hbasis=0
           state=0
           box=0
           
           call enumerate(nholes+1,2*norbs,2*norbs,ndimu,vect,box,state)
           
           do i=1,ndimu
              do k=1,nholes+1
                 hbasis(i,k)=vect(i,k)
              end do
           end do
           
           deallocate(vect)
           allocate(vect(ndimd,nholes))
           state=0
           box=0
           
           call enumerate(nholes,2*norbs,2*norbs,ndimd,vect,box,state)
           state=ndimu

           do i=1,ndimd
              do l=1,ncore*ncells
                 state=state+1
                 do k=1,nholes
                    hbasis(state,k)=vect(i,k)
                 end do
                 hbasis(state,nholes+1)=2*norbs+l
              end do
           end do
           
           do i=1,hdim_max
              if(i.eq.ndimu)write(6,*)'N+1 basis states done'
           end do
           
           hstates=0
           
           do i=1,hdim_max
              count=0
              count2=0
              do j=1,nholes+1
                 ip=mflag(hbasis(i,j))
                 count=count+momx(ip)
                 count2=count2+momy(ip)
                 if(flag_B.eq.1)then
                    if(count.ge.ncellsx)count=count-ncellsx
                    if(count2.ge.ncellsy)count2=count2-ncellsy
                    if(count.lt.0)count=count+ncellsx
                    if(count2.lt.0)count2=count2+ncellsy
                 else
                    if(count.gt.2)count=count-4
                    if(count.lt.-1)count=count+4
                    if(count2.gt.2)count2=count2-4
                    if(count2.lt.-1)count2=count2+4
                 end if
              end do
              do ii=1,ncells
                 if(count.eq.momx(ii).and.count2.eq.momy(ii))jj=ii
              end do
              hstates(jj)=hstates(jj)+1
              hmlabel(i)=jj
           end do
           
           hstatemax=1
           do i=1,ncells
              if(hstates(i).gt.hstatemax)hstatemax=hstates(i)
           end do
           
           write(6,*)'Hole addition states partitioned'
           
           allocate(hbinlabel(ncells,hstatemax))
           allocate(hpairs(ncells,hstatemax,hnpairs,2))
           allocate(hothers(ncells,hstatemax,hnpairs,nholes-1))
           hpairs=0
           hothers=0
           
           do k=1,ncells
              count2=0
              do n=1,hdim_max
                 if(hmlabel(n).eq.k)then
                    count2=count2+1
                    hbinlabel(k,count2)=n
                    count4=0
                    do i=1,nholes
                       do j=i+1,nholes+1
                          count4=count4+1
                          hpairs(k,count2,count4,1)=hbasis(n,i)
                          hpairs(k,count2,count4,2)=hbasis(n,j)
                          if(nholes+1.gt.2)then
                             count6=0
                             do pp=1,nholes+1
                                if(pp.ne.i.and.pp.ne.j)then
                                   count6=count6+1
                                   hothers(k,count2,count4,count6)=hbasis(n,pp)
                                end if
                             end do
                          end if
                       end do
                    end do
                 end if
              end do
              hbinstates(k)=count2
              hmtot(k)=k
           end do
           
           hqlist=0
           count=0
           do i=1,ncells
              if(hbinstates(i).ne.0)then
                 count=count+1
                 hqlist(count)=i
                 write(6,*)'momentum, number of states:',i,hbinstates(i)
              end if
           end do
           
           hqstates=count
           write(6,*)hqstates,hstatemax
           
           if(nholes.ge.1)then
              allocate(gsign((nholes)*(nholes+1)/2))
              call getfsign(nholes+1,gsign)
           else
              allocate(gsign(1))
              gsign(1)=0
           end if
           
        end if

        hNEV=hstatemax
        allocate(hdeig(ncells,hNEV))
        allocate(vh(ncells,hstatemax,hNEV))
        do qqq=1,hqstates
           ws=hqlist(qqq)
           ndim=hbinstates(ws)
           write(6,*)'Sector:         ','K:',momx(hmtot(ws)),momy(hmtot(ws)),'Hilbert space:',hbinstates(ws)
           write(6,*)'twist',twist1,twist2,twistx-1,twisty-1
           write(6,*)'------------------------------------------------------------------------------------------------------'
           write(6,*)'writing hole addition Hamiltonian...'
           if(usearpack.gt.1)then
              hNEV=ndim
              allocate(h(ndim,ndim))
              h=dfloat(0)
           else
              hNEV=min(100,ndim/2)
              NCV=2*hNEV
              allocate(tIndexI(min(max_array,ndim*ndim)))
              allocate(tIndexJ(min(max_array,ndim*ndim)))
              allocate(tsparseH(min(max_array,ndim*ndim)))
              count=0
           end if
           do n=1,hbinstates(ws)
              do m=n,hbinstates(ws)
                 sum=dfloat(0)
                 
                 do i=1,nholes+1
                    ii=hbasis(hbinlabel(ws,n),i)
                    factor1=(-1)**(i+1)
                    count4=0
                    do ip=1,nholes+1
                       if(ip.ne.i)then
                          count4=count4+1
                          haa(count4)=hbasis(hbinlabel(ws,n),ip)
                       end if
                    end do
                    
                    do j=1,nholes+1
                       jj=hbasis(hbinlabel(ws,m),j)
                       factor2=(-1)**(j+1)
                       count4=0
                       do ip=1,nholes+1
                          if(ip.ne.j)then
                             count4=count4+1
                             hbb(count4)=hbasis(hbinlabel(ws,m),ip)
                          end if
                       end do
                       
                       count4=0
                       if(nholes.gt.0)then
                          do k=1,nholes
                             do l=1,nholes
                                if(haa(k).eq.hbb(l))count4=count4+1
                             end do
                          end do
                       end if
                       if(count4.eq.nholes)then
                          factor=factor1*factor2
                          if(mflag(ii).eq.mflag(jj))then
                             if(uflag(ii).le.norbspc.and.uflag(jj).le.norbspc)then
                                if(usearpack.gt.1)then
                                   h(n,m)=h(n,m)+factor*(t(ii,jj)-i4*xtal(ii,jj))
                                else
                                   sum=sum+factor*(t(ii,jj)-i4*xtal(ii,jj))
                                end if
                             end if
                          end if
                          if(uflag(ii).gt.norbspc.and.uflag(jj).gt.norbspc)then
                             k=uflag(ii)-norbspc
                             l=uflag(jj)-norbspc
                             if(k.eq.l.and.sflag(ii).eq.sflag(jj))then
                                coree=ecore
                             else
                                coree=dfloat(0)
                             end if
                             if(sflag(ii).eq.-1)k=k+3
                             if(sflag(jj).eq.-1)l=l+3
                             if(usearpack.ge.1.and.usearpack.le.2)then
                                sum=sum+(core(k,l)+coree)*factor
                             else
                                h(n,m)=h(n,m)+(core(k,l)+coree)*factor
                             end if
                          end if
                       end if
                    end do
                 end do
                 
                 do i=1,hnpairs
                    factor1=gsign(i)
                    do j=1,hnpairs
                       factor2=gsign(j)
                       factor=factor1*factor2
                       prod=1
                       do pp=1,nholes-1
                          if(hothers(ws,n,i,pp).ne.hothers(ws,m,j,pp))prod=0
                       end do
                       if(prod.eq.1)then
                          l1=lflag(hpairs(ws,n,i,1)) !c
                          ii=mlflag(hpairs(ws,n,i,1)) !c
                          l2=lflag(hpairs(ws,n,i,2)) !d
                          jj=mlflag(hpairs(ws,n,i,2)) !d
                          l3=lflag(hpairs(ws,m,j,1)) !a
                          kk=mlflag(hpairs(ws,m,j,1)) !a
                          l4=lflag(hpairs(ws,m,j,2)) !b 
                          ll=mlflag(hpairs(ws,m,j,2)) !b 
                          !                             ! spin a=spin c and spin b = spin d
                          if(sflag(hpairs(ws,n,i,2))*sflag(hpairs(ws,m,j,2)).eq.1)then
                             if(sflag(hpairs(ws,n,i,1))*sflag(hpairs(ws,m,j,1)).eq.1)then
                                !                                   ! pa-pc=pd-pb
                                k=mdiff(mflag(hpairs(ws,n,i,1)),mflag(hpairs(ws,m,j,1)))
                                l=mdiff(mflag(hpairs(ws,m,j,2)),mflag(hpairs(ws,n,i,2)))
                                if(k.eq.l)then
                                   if(usearpack.ge.1.and.usearpack.le.2)then
                                      sum=sum+factor*0.5d0*(hint(l1,ii,l2,jj,l4,ll,l3,kk)+hint(l2,jj,l1,ii,l3,kk,l4,ll))/dfloat(ncells)
                                   else
                                      h(n,m)=h(n,m)+factor*0.5d0*(hint(l1,ii,l2,jj,l4,ll,l3,kk)+hint(l2,jj,l1,ii,l3,kk,l4,ll))/dfloat(ncells)
                                   end if
                                end if
                             end if
                          end if
                          !                          ! spin b=spin c and spin a = spin d
                          if(sflag(hpairs(ws,n,i,2))*sflag(hpairs(ws,m,j,1)).eq.1)then
                             if(sflag(hpairs(ws,n,i,1))*sflag(hpairs(ws,m,j,2)).eq.1)then
                                !                                   ! pa-pc=pd-pb
                                k=mdiff(mflag(hpairs(ws,n,i,2)),mflag(hpairs(ws,m,j,1)))
                                l=mdiff(mflag(hpairs(ws,m,j,2)),mflag(hpairs(ws,n,i,1)))
                                if(k.eq.l)then
                                   if(usearpack.ge.1.and.usearpack.le.2)then
                                      sum=sum-factor*0.5d0*(hint(l1,ii,l2,jj,l3,kk,l4,ll)+hint(l2,jj,l1,ii,l4,ll,l3,kk))/dfloat(ncells)
                                   else
                                      h(n,m)=h(n,m)-factor*0.5d0*(hint(l1,ii,l2,jj,l3,kk,l4,ll)+hint(l2,jj,l1,ii,l4,ll,l3,kk))/dfloat(ncells)
                                   end if
                                end if
                             end if
                          end if
                       end if
                    end do
                 end do

                 if(usearpack.ge.1.and.usearpack.le.2)then
                    if(abs(sum).gt.small)then
                       count=count+1
                       tIndexI(count)=n
                       tIndexJ(count)=m
                       tsparseH(count)=sum
                       if(n.ne.m)then
                          count=count+1
                          tIndexI(count)=m
                          tIndexJ(count)=n
                          tsparseH(count)=sum
                       end if
                    end if
                 else
                    h(m,n)=h(n,m)
                 end if
                 
              end do
           end do
           
           if(usearpack.le.1)write(6,*)count
           
           write(6,*)'hole addition Hamiltonian diagonalizing...'
           write(6,*)
           if(usearpack.gt.1)then
              allocate(w(ndim))
              iwork=3*ndim-1
              allocate(work(iwork))
              call dsyev(jobz,uplo,ndim,h,ndim,w,work,iwork,info)
              do i=1,ndim
                 hdeig(ws,i)=w(i)
                 do j=1,ndim
                    vh(ws,j,i)=h(j,i)
                 end do
              end do
              deallocate(h)
              deallocate(w)
              deallocate(work)
           end if
           if(usearpack.le.1)then
              spr_size=count
              write(6,*)'spr_size:',count
              allocate(IndexI(spr_size))
              allocate(IndexJ(spr_size))
              allocate(sparseH(spr_size))
              do i=1,spr_size
                 IndexI(i)=tIndexI(i)
                 IndexJ(i)=tIndexJ(i)
                 sparseH(i)=tsparseH(i)
              end do
              deallocate(tIndexI)
              deallocate(tIndexJ)
              deallocate(tsparseH)
              allocate(DD(NCV))
              allocate(ZZ(ndim,NCV))
              call arpackdiag(spr_size,ndim,NEV,NCV,IndexI,IndexJ,sparseH,DD,ZZ)
              
              do i=1,hNEV
                 hdeig(ws,i)=DD(i)
                 do j=1,ndim
                    vh(ws,j,i)=ZZ(j,i)
                 end do
              end do
              
              deallocate(IndexI)
              deallocate(IndexJ)
              deallocate(sparseH)
              deallocate(DD)
              deallocate(ZZ)
           end if
           
           ! for clarity's sake, v(ws,j,i) is the jth eigenvector element of the ith eigenvalue for momentum sector ws
           !
           write(6,*)'--------------------------------------------------------------------------------------------------'
           write(6,*)'Quantum number, momentum [units pi/Nx,pi/Ny], spin, eigenvalue (lower set of states)'
           write(6,*)'--------------------------------------------------------------------------------------------------'
           write(16,*)'Quantum number, momentum units [pi/Nx,pi/Ny], spin, eigenvalue (lower set of states)'
           write(6,*)
           write(16,*)
           do n=1,min(40,hNEV)
              write(6,*)ws,momx(hmtot(ws)),momy(hmtot(ws)),n,sngl(hdeig(ws,n))
              write(16,*)ws,momx(hmtot(ws)),momy(hmtot(ws)),n,sngl(hdeig(ws,n))
           end do
           write(6,*)
           write(16,*)
        end do
        
        hmin=100.0d0
        hmax=-100.0d0
        do qqq=1,hqstates
           ws=hqlist(qqq)
           do j=1,hbinstates(ws)
              if(hdeig(ws,j).lt.hmin)hmin=hdeig(ws,j)
              if(hdeig(ws,j).gt.hmax)hmax=hdeig(ws,j)
           end do
        end do
        
        gamma=0.01d0
        sum2=dfloat(0)
        if(twist1.eq.0.and.twist2.eq.0)open(unit=58,file='hDOS.dat')
        do i=1,100*nfreq
           omega=-1.0d0+hmin+dfloat(i)*(2.0d0-hmin+hmax)/dfloat(100*nfreq)
           sum=0.0d0
           do qqq=1,hqstates
              ws=hqlist(qqq)
              if(usearpack.gt.1)then
                 do j=1,hbinstates(ws)
                    sum=sum+(1.0d0/pi)*gamma/(gamma*gamma+(omega-hdeig(ws,j))**2)/dfloat(ncells*twistx*twisty)
                 end do
              else 
                 do j=1,hNEV
                    sum=sum+(1.0d0/pi)*gamma/(gamma*gamma+(omega-hdeig(ws,j))**2)/dfloat(ncells*twistx*twisty)
                 end do
              end if
           end do
           write(58,*)omega,sum
           sum2=sum2+sum*(2.0d0-hmin+hmax)/dfloat(100*nfreq)
        end do
        close(58)
        write(6,*)'hole addition DOS calculated, weight:',sum2
     end if
        
     if(flag_R.eq.1) goto 888        
     !
     ! This part calculates XAS & RIXS
     !
     if(myrank.eq.0)then
        write(6,*)
        write(6,*)'Finding XAS & RIXS overlap states between valence-hole and core-hole sectors...(twist1,twist2)',twist1,twist2
        write(6,*)
     end if
     
     if(twist1.eq.0.and.twist2.eq.0)then
        cmax=min(cstates(momgs)*statemax,max_array)
        if(cmax.lt.0)cmax=max_array
        if(cmax.lt.0)then
           write(6,*)'cmax cannot convert'
           goto 999
        end if
        allocate(ivert(qnumbers,gstatedegen,nevn))
        allocate(fvert(qnumbers,nevn,nev))
        allocate(lap(qnumbers,cmax))
        allocate(clap(qnumbers,cmax))
        allocate(flap(qnumbers,cmax))
        allocate(dh(qnumbers,cmax))
        allocate(ph(qnumbers,cmax))
        allocate(tivert(qnumbers,gstatedegen,nevn,-nb*timesteps:nt*timesteps))
        allocate(tpivert(qnumbers,gstatedegen,nevn,-nb*timesteps:nt*timesteps,-nb*timesteps:nt*timesteps))
        allocate(tfvert(qnumbers,nev,nevn,-nb*timesteps:nt*timesteps))
        allocate(tkern(-nb*timesteps:nt*timesteps,-nb*timesteps:nt*timesteps))
        spectra=dfloat(0)
        
        ivert=dcmplx(0.0d0,0.0d0)
        fvert=dcmplx(0.0d0,0.0d0)
        tivert=dcmplx(0.0d0,0.0d0)
        tfvert=dcmplx(0.0d0,0.0d0)
        tpivert=dcmplx(0.0d0,0.0d0)
        
        ! identifies the core and valence holes involved in the transition between states core states n and valence states ws,m
        do qqq=1,qstates
           ws=qlist(qqq)
           count=0
           do m=1,binstates(ws)
              do n=1,cstates(momgs)
                 prod=0
                 do i=1,nholes
                    count2=1
                    do j=1,nholes-1
                       if(basis(binlabel(ws,m),i).eq.cbasis(cbinlabel(n),j))then
                          prod=prod+1
                          count2=0
                       end if
                    end do
                    if(count2.ne.0)k=basis(binlabel(ws,m),i)
                 end do
                 if(prod.eq.nholes-1)then
                    factor1=1
                    do i=1,nholes
                       if(basis(binlabel(ws,m),i).lt.k)then
                          factor1=-factor1
                       else
                          goto 5555
                       end if
                    end do
5555                continue
                    jj=cbasis(cbinlabel(n),nholes)    ! core hole location
                    count=count+1
                    flap(ws,count)=factor1*(-1)**(nholes-1)
                    lap(ws,count)=m           !valence basis state
                    clap(ws,count)=n          !core basis state
                    ph(ws,count)=cbasis(cbinlabel(n),nholes)    ! core hole location
                    dh(ws,count)=k                              ! valence hole eliminated
                 end if
              end do
           end do
           istate(ws)=count
           if(myrank.eq.0)write(6,*)'total # overlap states:',count,'in sector:',ws,cstates(momgs),binstates(ws)
        end do
     end if
177  continue

     if(myrank.eq.0)then
        write(6,*)
        write(6,*)'determining polarization dependence...'
        write(6,*)
     end if
     if(twist1.eq.0.and.twist2.eq.0)then
        omegai_min=-3.0d0+edn(1)-a(1)
        omegai_max=5.0d0-edn(1)+edn(NEVN)
        allocate(ckern(3,3,qnumbers,NEV,ifreq,qnumbers,gstatedegen))
        ckern=dcmplx(dfloat(0),dfloat(0))
     end if
     
     XAS_maxx=dfloat(0)
     do alpha=1,1
!     do alpha=1,1
        inx=0
        iny=0
        inz=0
        if(alpha.eq.1)inx=1
        if(alpha.eq.2)iny=1
        if(alpha.eq.3)inz=1
        do beta=2,3
!        do beta=1,1
           fnx=0
           fny=0
           fnz=0
           if(beta.eq.1)fnx=1
           if(beta.eq.2)fny=1
           if(beta.eq.3)fnz=1
           if(myrank.eq.0)write(6,*)'incident, scattered photon dipole matrix elements: nx,ny,nz=',inx,iny,inz,fnx,fny,fnz
           ! constructs the matrix elements for incident (ivert) and scattered (fvert) photon polarizations
           do qqq=1,qstates
              ws=qlist(qqq)
              allocate(ioverlap(istate(ws)))
              allocate(foverlap(istate(ws)))
              allocate(cl(istate(ws)))
              allocate(nl(istate(ws)))
              allocate(cf(istate(ws)))
              allocate(nf(istate(ws)))
              ioverlap=dcmplx(0.0d0,0.0d0)
              foverlap=dcmplx(0.0d0,0.0d0)
              
              count2=0
              count4=0
              do i=1,istate(ws)
                 if(sflag(ph(ws,i))*sflag(dh(ws,i)).eq.1)then
                    k=lflag(dh(ws,i))
                    kk=mlflag(dh(ws,i))
                    j=lflag(ph(ws,i))
                    jj=mlflag(ph(ws,i))
                    u1=uflag(dh(ws,i))
                    u2=uflag(ph(ws,i))
                    count=0
                    count3=0
                    if(edge.eq.1)then
                       if(u1.le.5.and.u2.ge.1+norbspc)count=1
                    end if
                    if(edge.eq.2)then
                       if(u2.eq.1+norbspc)then
                          if(u1.ge.6.and.u1.le.8)count=1
                       end if
                       if(u2.eq.2+norbspc)then
                          if(u1.ge.9.and.u1.le.11)then
                             count=1
                             count3=1
                          end if
                       end if
                    end if
                    if(count.eq.1)then
                       if(abs(gaunt(k,kk,j,jj,1)).gt.small)then
                          if(gsec(ws).eq.1)then
                             csum=dcmplx(dfloat(0),dfloat(0))
                             if(kk-jj.eq.1)csum=dcmplx(-dfloat(inx),dfloat(iny))/dsqrt(2.0d0)
                             if(kk-jj.eq.0)csum=dcmplx(dfloat(inz),0.0d0)
                             if(kk-jj.eq.-1)csum=dcmplx(dfloat(inx),dfloat(iny))/dsqrt(2.0d0)
                             if(cdabs(csum).gt.small)then
                                count2=count2+1
                                ioverlap(count2)=csum*gaunt(k,kk,j,jj,1)*flap(ws,i)
                                !                    lap(ws,count)=m           !valence basis state
                                !                    clap(ws,count)=n          !core basis state
                                nl(count2)=lap(ws,i)
                                cl(count2)=clap(ws,i)
                             end if
                          end if
                          csum=dcmplx(dfloat(0),dfloat(0))
                          if(kk-jj.eq.1)csum=dcmplx(-dfloat(fnx),dfloat(fny))/dsqrt(2.0d0)
                          if(kk-jj.eq.0)csum=dcmplx(dfloat(fnz),0.0d0)
                          if(kk-jj.eq.-1)csum=dcmplx(dfloat(fnx),dfloat(fny))/dsqrt(2.0d0)
                          if(cdabs(csum).gt.small)then
                             count4=count4+1
                             foverlap(count4)=csum*gaunt(k,kk,j,jj,1)*flap(ws,i)
                             nf(count4)=lap(ws,i)
                             cf(count4)=clap(ws,i)
                          end if
                       end if
                    end if
                 end if
              end do
              
              do n=1,min(nevn,cstates(momgs))
                 if(gsec(ws).eq.1)then
                    do j=1,gdegen(ws)
                       csum=dcmplx(0.0d0,0.0d0)
                       do i=1,count2
                          csum=csum+ioverlap(i)*ev(ws,nl(i),j)*evn(cl(i),n)
                       end do
                       ivert(ws,j,n)=csum
                    end do
                 end if
                 do m=1,min(nev,binstates(ws))
                    csum=dcmplx(0.0d0,0.0d0)
                    do i=1,count4
                       csum=csum+foverlap(i)*ev(ws,nf(i),m)*evn(cf(i),n)
                    end do
                    fvert(ws,n,m)=csum
                 end do
              end do
              do ti=-nb*timesteps,nt*timesteps
                 do n=1,min(nevn,cstates(momgs))
                    !
                    ! NOTE:
                    ! tivert(ws,j,n,ti) is the matrix element ioverlap between the time-evolved ground state (ws,j) at time ti and the time evolved core hole state n at time ti
                    ! = <Phi_n(t)|dipole|Psi_gs(t)>
                    !
                    if(gsec(ws).eq.1)then
                       do j=1,gdegen(ws)
                          csum=dcmplx(0.0d0,0.0d0)
                          do i=1,count2
                             csum=csum+ioverlap(i)*tv(ws,nl(i),j,ti)*evn(cl(i),n)
                          end do
                          tivert(ws,j,n,ti)=csum
                       end do
                    end if
                    !
                    ! NOTE:
                    ! tfvert(ws,m,n,ti) is the matrix element foverlap between the time-evolved state (ws,m) at time ti and the time evolved core hole state n at time ti
                    ! = <Phi_n(t)|dipole|Psi_m(t)>
                    !
                    do m=1,min(nev,binstates(ws))
                       csum=dcmplx(0.0d0,0.0d0)
                       do i=1,count4
                          csum=csum+foverlap(i)*tv(ws,nf(i),m,ti)*evn(cf(i),n)
                       end do
                       tfvert(ws,m,n,ti)=csum
                    end do
                 end do
              end do

              deallocate(ioverlap)
              deallocate(foverlap)
              deallocate(cl)
              deallocate(nl)
              deallocate(cf)
              deallocate(nf)

              !
              !confirmed tivert(ti)~exp(-i (ti-t0) E_gs)
              !confirmed tfvert(n,ti)~exp(-i (ti-t0) E_n)
              !
              !confirmed tpivert(tj,ti)~exp(-i(tj-ti)*E_CH(m)-i (ti-t0) E_gs)
              !
              ! NOTE:
              ! tpivert(ws,j,m,tj,ti) propagates the core hole state n at time ti (reached from ground state ws, j) to state m at time tj
              ! summed over all core hole intermediate states n
              ! = sum_n U_{m,n}(tj,ti) <Phi_n(ti)|dipole|Psi_j(ti)>
              !
              !

              if(gsec(ws).eq.1)then
                 do j=1,gdegen(ws)
                    if(myrank.eq.0)then
                       write(6,*)'Calculating full evolution propagator for core hole states needed for RIXS for ground state ws,j...',ws,j
                    end if
                    allocate(tevolution(-nb*timesteps:nt*timesteps,ndim,ndim))
              
                    allocate(biga(ndim,ndim))
                    allocate(bigb(ndim,ndim))
                    allocate(bigc(ndim,ndim))
                    allocate(bigd(ndim,ndim))
                    allocate(bige(ndim,ndim))
                    do i_time1=-nb*timesteps,nt*timesteps
                       if(myrank.eq.0)write(6,*)i_time1,nt*timesteps
                       do i_time2=i_time1+1,nt*timesteps
                          if(i_time2.eq.i_time1+1)then
                             tevolution=dcmplx(0.0d0,0.0d0)
                             do n=1,ndim
                                do i_time=-nb*timesteps,nt*timesteps
                                   tevolution(i_time,n,n)=dcmplx(1.0d0,0.0d0)
                                end do
                             end do
                          end if
                          do n=1,ndim
                             do m=1,ndim
                                biga(n,m)=evolution(i_time2,n,m)
                                bigb(n,m)=tevolution(i_time2-1,n,m)
                             end do
                          end do
                          bigc=matmul(biga,bigb)
                          do n=1,ndim
                             do m=1,ndim
                                tevolution(i_time2,n,m)=bigc(n,m)
                             end do
                          end do
                          
                          bigd=matmul(bigc,evn)
                          bige=matmul(transpose(evn),bigd)

                          do m=1,min(nevn,cstates(momgs))
                             csum=dcmplx(0.0d0,0.0d0)
                             do n=1,min(nevn,cstates(momgs))
                                csum=csum+tivert(ws,j,n,i_time1)*bige(n,m)
                             end do
                             tpivert(ws,j,m,i_time2,i_time1)=csum
                          end do
                       end do
                    end do

                    deallocate(tevolution)
                    deallocate(biga)
                    deallocate(bigb)
                    deallocate(bigc)
                    deallocate(bigd)
                    deallocate(bige)
                 end do
              
                 if(myrank.eq.0)write(6,*)'Done'
                 
              end if
           end do

!           deallocate(evolution)
           gammax=0.1d0
           
           if(beta.eq.1)then ! WE ARE CHANGING TO THIS TO ALSO GET Y XAS
              if(twist1.eq.0.and.twist2.eq.0)then
                 sum_max=dfloat(0)
                 if(myrank.eq.0)then
                    write(6,*)'Calculating equilibrium XAS...twist1,twist2',twist1,twist2
                    if(alpha.eq.1.and.flag_skip_XAS.eq.1)then
                       bella=6666
                       reid=78
                       open(unit=bella,file='XAS_X.dat')
                       open(unit=reid,file='trXAS_X.dat')
                    end if
                    if(alpha.eq.2.and.flag_skip_XAS.eq.1)then
                       bella=6667
                       reid=79
                       open(unit=bella,file='XAS_Y.dat')
                       open(unit=reid,file='trXAS_Y.dat')
                    end if
                    if(alpha.eq.3.and.flag_skip_XAS.eq.1)then
                       bella=6668
                       reid=80
                       open(unit=bella,file='XAS_Z.dat')
                       open(unit=reid,file='trXAS_Z.dat')
                    end if
                 end if
              end if
              do m=1,ifreq
                 omegai=omegai_min+omegai_max*dfloat(m)/dfloat(ifreq)
                 sum=dfloat(0)
                 sum2=dfloat(0)
                 do i=1,min(nevn,cstates(momgs))
                    lor=gammax/(gammax*gammax+(omegai-(edn(i)-a(1)))**2)
                    do qqq=1,qstates
                       ws=qlist(qqq)
                       if(gsec(ws).eq.1)then
                          do j=1,gdegen(ws)
                             sum2=sum2+lor*(ivert(ws,j,i)*dconjg(ivert(ws,j,i)))
                          end do
                          sum=sum+lor
                       end if
                    end do
                 end do
                 if(sum2.gt.sum_max)then
                    sum_max=sum2
                    XAS_max=omegai
                 end if
                 if(myrank.eq.0)write(bella,*)omegai,sum2,sum,XAS_max
              end do
              if(myrank.eq.0)write(6,*)alpha,'XAS peak at',XAS_max,'value',sum_max,'for twist',twist1,twist2
              if(XAS_max.gt.XAS_maxx)XAS_maxx=XAS_max
              ! this does the time-domain part
              ! make the time-dependent kernal tkern(ti,tj)
              if(myrank.eq.0)then
                 write(6,*)'Calculating time-domain XAS...twist1,twist2',twist1,twist2
              end if
              tkern=dcmplx(0.0d0,0.0d0)
              do ti=-nb*timesteps,nt*timesteps
                 do tj=-nb*timesteps,nt*timesteps
                    do qqq=1,qstates
                       ws=qlist(qqq)
                       if(gsec(ws).eq.1)then
                          do j=1,min(nevn,cstates(momgs))
                             tkern(ti,tj)=tkern(ti,tj)+dconjg(tivert(ws,k,j,tj))*tpivert(ws,k,j,tj,ti)
                             ! Restric caclulation to be in the same sector of the groundstate. If you see it it must
                             ! a dd without spin flip.
                          end do
                       end if
                    end do
                 end do
              end do
              !
              ! CALCULATES TIME-RESOLVED XAS
              !
              !              gauss=dfloat(120)*pi/dfloat(120*(nt+nb))
!              gauss=dfloat(640)*pi/dfloat(timesteps*(nt+nb))
              gauss=dfloat(160)*pi/dfloat(320*8)
!              do tm=-nb*timesteps/2,nt*timesteps,6
               if(flag_skip_XAS.eq.1)then
                  do tm=-2*timesteps,4*timesteps,6
      !                 if(myrank.eq.0)write(6,*)tm
                     tmeas=pi*dfloat(2*tm)/dfloat(timesteps)
                     if(myrank.eq.0)write(6,*)tmeas
                     do m=1,ifreq
                        omegai=omegai_min+omegai_max*dfloat(m)/dfloat(ifreq)
                        csum=dcmplx(0.0d0,0.0d0)
                        do ti=-nb*timesteps+myrank,nt*timesteps,numproc
                           pp_time=pi*dfloat(2*ti)/dfloat(timesteps)
                           s1=gauss*dexp(-0.5d0*gauss*gauss*(pp_time-tmeas)**2)/dsqrt(2.0d0*pi)
                           !                       s1=1.0d0
                           do tj=ti+1,nt*timesteps
                              p_time=pi*dfloat(2*tj)/dfloat(timesteps)
                              s2=gauss*dexp(-0.5d0*gauss*gauss*(p_time-tmeas)**2)/dsqrt(2.0d0*pi)
                              !                          s2=1.0d0
                              csum=csum+4.0d0*pi*pi*s1*s2*tkern(ti,tj)*cdexp(eye*omegai*(p_time-pp_time))/dfloat(timesteps*timesteps)
                           end do
                        end do
                        call MPI_REDUCE(csum,csum4,siz,MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,mpistat)
                        csum=csum4
                        if(myrank.eq.0)write(reid,*)tmeas,omegai,dreal(csum),dimag(csum)
                     end do
                     if(myrank.eq.0)write(reid,*)
                  end do
               end if
345           continue
              if(myrank.eq.0.and.flag_skip_XAS.eq.1)then
                 close(reid)
                 close(bella)
              end if
           end if
           
!           goto 853
           
           !
           ! CALCULATES TIME-RESOLVED RIXS
           !
           if(beta.eq.1)then
              if(myrank.eq.0)write(6,*)'Calculating time-evolved RIXS - twist1,twist2:',twist1,twist2
           end if
           
           count2=0
           !           gauss=dfloat(240)*pi/dfloat(timesteps*(nt+nb))
           !           gauss=dfloat(60)*pi/dfloat(timesteps*(nt+nb))
           !           gauss=dfloat(240)*pi/dfloat(320*8)
           ! v2 diagnostic test: narrowed probe (sigma=1.0, gauss=1/sigma) with delays kept well
           ! inside the existing nb=nt=1 window, so the same (t1,t2) grid/tpivert build is reused
           ! at no extra compute cost, but every tested delay sits >=3 probe-sigma from either edge.
           gauss=1.0d0
           tch=dfloat(2)*pi/1.5d0
           do qqq=1,qstates
              ws=qlist(qqq)
              if(gsec(ws).eq.1)then
                 if(count2.eq.0)then
                    count2=count2+1
                    if(myrank.eq.0)write(6,*)'DEBUG which ws processed: ws=',ws,'2Sz=',sztot(ws),'binstates=',binstates(ws)
                    do ppp=1,qstates
                       wsp=qlist(ppp)
                       tmom=mdiff(mtot(ws),mtot(wsp))
                       tspin=0
                       if(dlso.eq.0.0d0)tspin=sztot(ws)-sztot(wsp)
                       if(tspin.eq.0)then
                          count=0
                          tmeas=pi*dfloat(2*tm)/dfloat(timesteps)
                          !                          do tm=-nb*timesteps/2,nt*timesteps,60
                          ! FIX: this allocate sits inside the ppp (wsp) loop, and with SOC on
                          ! (dlso != 0) the "if(dlso.eq.0.0d0)tspin=..." line above never fires, so
                          ! tspin stays 0 and "tspin.eq.0" is trivially true for EVERY ppp -- the
                          ! original code reallocates allowed_tm on every matching ppp without ever
                          ! deallocating, crashing under -fcheck=all on the second match (confirmed:
                          ! the original toy run did the exact same thing and only ever completed the
                          ! FIRST matching wsp channel before crashing -- see
                          ! /private/tmp/fortran_ref_build/run.log and fortran_ref_build3/run.log).
                          if (allocated(allowed_tm)) deallocate(allowed_tm)
                          allocate(allowed_tm(5))
                          allowed_tm = [-40,-20,0,20,40]
                          do tm=-2*timesteps,4*timesteps,2
                              do tm_i = 1, size(allowed_tm)
                                 if (tm == allowed_tm(tm_i)) then
                                    tmeas=pi*dfloat(2*tm)/dfloat(timesteps)
                                    if(dlso.eq.0.0d0)then
                                       tlabel=100*count+10*tmom+abs(tspin)
                                    else
                                       tlabel=100*count+10*tmom
                                    end if
                                    ! v2 diagnostic fix: tlabel alone collides across distinct wsp
                                    ! (ppp) channels that happen to share the same tmom/tspin -- this
                                    ! silently overwrote earlier channels' files (status='unknown').
                                    ! Add a ppp-dependent offset so every wsp gets its own file set.
                                    tlabel=tlabel+1000*ppp
                                    if(myrank.eq.0)then
                                       write(6,*)'Polarization:',alpha,beta
                                       write(6,*)'Quantum numbers (time, Delta S, Delta Q:',count,tspin,tmom
                                       if(alpha.eq.1.and.beta.eq.1)then
                                          write(filename,'("trSPEC_XX_",I4.4,".dat")')tlabel
                                          open(unit=tlabel+321,file=filename,status='unknown')
                                       end if
                                       if(alpha.eq.1.and.beta.eq.2)then
                                          write(filename,'("trSPEC_XY_",I4.4,".dat")')tlabel
                                          open(unit=tlabel+321,file=filename,status='unknown')
                                       end if
                                       if(alpha.eq.1.and.beta.eq.3)then
                                          write(filename,'("trSPEC_XZ_",I4.4,".dat")')tlabel
                                          open(unit=tlabel+321,file=filename,status='unknown')
                                       end if
                                       if(alpha.eq.2.and.beta.eq.1)then
                                          write(filename,'("trSPEC_YX_",I4.4,".dat")')tlabel
                                          open(unit=tlabel+321,file=filename,status='unknown')
                                       end if
                                       if(alpha.eq.2.and.beta.eq.2)then
                                          write(filename,'("trSPEC_YY_",I4.4,".dat")')tlabel
                                          open(unit=tlabel+321,file=filename,status='unknown')
                                       end if
                                       if(alpha.eq.2.and.beta.eq.3)then
                                          write(filename,'("trSPEC_YZ_",I4.4,".dat")')tlabel
                                          open(unit=tlabel+321,file=filename,status='unknown')
                                       end if
                                       if(alpha.eq.3.and.beta.eq.1)then
                                          write(filename,'("trSPEC_ZX_",I4.4,".dat")')tlabel
                                          open(unit=tlabel+321,file=filename,status='unknown')
                                       end if
                                       if(alpha.eq.3.and.beta.eq.2)then
                                          write(filename,'("trSPEC_ZY_",I4.4,".dat")')tlabel
                                          open(unit=tlabel+321,file=filename,status='unknown')
                                       end if
                                       if(alpha.eq.3.and.beta.eq.3)then
                                          write(filename,'("trSPEC_ZZ_",I4.4,".dat")')tlabel
                                          open(unit=tlabel+321,file=filename,status='unknown')
                                       end if
                                    end if
                                    count=count+1
                                    ! get the time-dependent RIXS kernal
                                    tkern=dcmplx(0.0d0,0.0d0)
                                    do t1=-nb*timesteps,nt*timesteps
                                       do t2=t1+1,nt*timesteps
                                          csum=dcmplx(0.0d0,0.0d0)
                                          do k=1,gdegen(ws)
                                             do l=1,min(nev,binstates(wsp))
                                                do i=1,min(nevn,cstates(momgs))
                                                   csum=csum+dconjg(tfvert(wsp,l,i,t2))*tpivert(ws,k,i,t2,t1)
                                                end do
                                             end do
                                          end do
                                          tkern(t1,t2)=csum
                                       end do
                                    end do
                                    do m=1,ifreq
                                       omegai=omegai_min+omegai_max*dfloat(m)/dfloat(ifreq)
                                       if(myrank.eq.0)write(6,*)tm,m,ifreq
                                       do n=1,nfreq
                                          omega=dfloat(8*n)/dfloat(nfreq)-2.0d0
                                          omegas=omegai-omega
                                          csum=dcmplx(0.0d0,0.0d0)
                                          do t1=-nb*timesteps+myrank,nt*timesteps,numproc
                                             p1_time=pi*dfloat(2*t1)/dfloat(timesteps)
                                             s1=gauss*dexp(-0.5d0*gauss*gauss*(p1_time-tmeas)**2)/dsqrt(2.0d0*pi)
                                             !                                s1=dfloat(1)
                                             do t2=t1+1,nt*timesteps
                                                p2_time=pi*dfloat(2*t2)/dfloat(timesteps)
                                                s2=gauss*dexp(-0.5d0*gauss*gauss*(p2_time-tmeas)**2)/dsqrt(2.0d0*pi)
                                                !                                   s2=dfloat(1)
                                                l12=dexp(-dabs(p1_time-p2_time)/tch)
                                                csum=csum+l12*s1*s2*tkern(t1,t2)*cdexp(-eye*(omegai*p1_time-omegas*p2_time))*(dfloat(2)*pi/dfloat(timesteps))**2
                                             end do
                                          end do
                                          call MPI_REDUCE(csum,csum4,siz,MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,mpistat)
                                          !                                call MPI_REDUCE(sum2,sum4,siz,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,mpistat)
                                          if(myrank.eq.0)write(tlabel+321,*)omegai,omega,dreal(csum4*dconjg(csum4))
                                       end do
                                       if(myrank.eq.0)write(tlabel+321,*)
                                    end do
                                    if(myrank.eq.0)close(tlabel+321)
                                 end if
                              end do
                          end do
                       end if
                    end do
                 end if
              end if
           end do
           
           do qqq=1,qstates
              ws=qlist(qqq)
              if(gsec(ws).eq.1)then
                 do m=1,gdegen(ws)
                    do ppp=1,qstates
                       wsp=qlist(ppp)
                       do n=1,ifreq
                          omegai=omegai_min+omegai_max*dfloat(n)/dfloat(ifreq)
                          do i=1,NEV
                             csum=dcmplx(0.0d0,0.0d0)
                             do j=1,min(nevn,cstates(momgs))
                                csum=csum+ivert(ws,m,j)*dconjg(fvert(wsp,j,i))/(omegai-(dn(j)-a(1))+eye*gammax)
                             end do
                             ckern(alpha,beta,wsp,i,n,ws,m)=csum
                          end do
                       end do
                    end do
                 end do
              end if
           end do
853        continue
        end do
     end do
     
     !     goto 888
     
     if(myrank.eq.0)then
        write(6,*)'---------------------'
        write(6,*)'Assembling equilibrium spectra...'
        write(6,*)'---------------------'
     end if
     gamma=0.01d0
     do alpha=1,3
        do beta=1,3
           count2=0
           do qqq=1,qstates
              ws=qlist(qqq)
              if(gsec(ws).eq.1)then
                 if(count2.eq.0)then
                    count2=count2+1
                    do ppp=1,qstates
                       wsp=qlist(ppp)
                       do m=1,nfreq
                          omega=dfloat(8*m)/dfloat(nfreq)-2.0d0
                          do n=1,ifreq
                             omegai=omegai_min+omegai_max*dfloat(n)/dfloat(ifreq)
                             rixs=dfloat(0)
                             do j=1,gdegen(ws)
                                do i=1,min(nev,binstates(wsp))
                                   rixs=rixs+dreal(ckern(alpha,beta,wsp,i,n,ws,j)*dconjg(ckern(alpha,beta,wsp,i,n,ws,j)))*gamma/((omega-(deig(wsp,i)-a(1)))**2+gamma**2)
                                end do
                             end do
                             spectra(ws,alpha,beta,wsp,n,m)=spectra(ws,alpha,beta,wsp,n,m)+rixs/dfloat(twistx*twisty)
                          end do
                       end do
                    end do
                 end if
              end if
           end do
        end do
     end do
     
     deallocate(dn)
     deallocate(vn)
     
888  continue
  
     if(flag_X.eq.1)goto 887
     !
     ! This part is from the old code for equilbrium stuff. DO NOTE USE!!!!!!!
     !
     open(unit=166,file='XPS.dat')
     do extk=1,ncells
        !determine which states can be connected for a given momentum of photoemitted electrons
        count=0
        xpslist1=0
        xpslist2=0
        do n=1,qstates
           ws=qlist(n)
           p=mtot(ws)
           do m=1,ncells
              if(extk.eq.mdiff(m,p))then
                 count=count+1
                 xpslist1(count)=n
                 xpslist2(count)=m
              end if
           end do
        end do
        ! determine overlap states for the quantum number state pairs
        xps=dfloat(0)
        do i=1,count
           ws=qlist(xpslist1(i))
           wsp=xpslist2(i)
           do j=1,2*norbs+ncore*ncells
              count2=0
              if(mflag(j).eq.extk)then
                 do n=1,binstates(ws)
                    prod=1
                    do k=1,nholes
                       if(basis(binlabel(ws,n),k).eq.j)prod=0
                    end do
                    if(prod.eq.1)then
                       !                          write(6,*)i,j,(basis(binlabel(ws,n),k),k=1,nholes)
                       do m=1,hbinstates(wsp)
                          count=0
                          prod=0
                          do k=1,nholes+1
                             if((hbasis(hbinlabel(wsp,m),k)).eq.j)then
                                prod=1
                                factor=(-1)**k
                             end if
                          end do
                          if(prod.eq.1)then
                             count=0
                             do k=1,nholes+1
                                prod=0
                                do p=1,nholes
                                   if(basis(binlabel(ws,n),p).eq.hbasis(hbinlabel(wsp,m),k))prod=1
                                end do
                                count=count+prod
                             end do
                             if(count.eq.nholes)then
                                count2=count2+1
                                xps1(j,count2)=n
                                xps2(j,count2)=m
                                xps3(j,count2)=factor
                             end if
                          end if
                       end do
                    end if
                 end do
                 countj(j)=count2
              end if
           end do
           do n=1,binstates(ws)
              do m=1,hbinstates(wsp)
                 xpsmatrix=0
                 do j=1,2*norbs+ncore*ncells
                    sum=dfloat(0)
                    do qqq=1,countj(j)
                       ii=xps1(j,qqq)
                       jj=xps2(j,qqq)
                       sum=sum+v(ws,ii,n)*vh(wsp,jj,m)*dfloat(xps3(j,qqq))
                    end do
                 end do
                 do j=1,nfreq
                    omega=dfloat(j)*(hmax-a(1))/dfloat(nfreq)
                    xps(j)=xps(j)+xpsmatrix*(gamma/pi)/(gamma*gamma+(omega-(hdeig(wsp,m)-deig(ws,n)))**2)
                 end do
              end do
           end do
        end do
        do j=1,nfreq
           omega=dfloat(j)*(hmax-a(1))/dfloat(nfreq)
           if(myrank.eq.0)write(166,*)extk,omega,xps(j)
        end do
        if(myrank.eq.0)write(166,*)
     end do
     if(myrank.eq.0)close(166)
     
887  continue
     
     deallocate(degen)
     deallocate(v)
     
     if(twistx*twisty.ne.1)then
        write(6,*)'incrementing twist',twist1,twist2,'(',twistx,twisty,')'
        write(6,*)
     end if
  end do
end do

if(flag_R.ne.1)then
   !     goto 555
   if(myrank.eq.0)write(6,*)'Writing RIXS to files...'
   do alpha=1,3
        do beta=1,3
           !   do alpha=1,1
           !        do beta=1,1
           count2=0
           do qqq=1,qstates
              ws=qlist(qqq)
              if(gsec(ws).eq.1)then
                 if(count2.eq.0)then
                    count2=count2+1
                    do ppp=1,qstates
                       wsp=qlist(ppp)
                       tmom=mdiff(mtot(ws),mtot(wsp))
                       if(dlso.eq.0.0d0)then
                          tspin=sztot(ws)-sztot(wsp)
                          tlabel=10*tmom+abs(tspin)
                       else
                          tlabel=10*tmom
                       end if
                       if(myrank.eq.0)then
                          if(alpha.eq.1.and.beta.eq.1)then
                             write(filename,'("SPEC_XX_",I2.2,".dat")')tlabel
                             open(unit=tlabel+654,file=filename,status='unknown')
                          end if
                          if(alpha.eq.1.and.beta.eq.2)then
                             write(filename,'("SPEC_XY_",I2.2,".dat")')tlabel
                             open(unit=tlabel+654,file=filename,status='unknown')
                          end if
                          if(alpha.eq.1.and.beta.eq.3)then
                             write(filename,'("SPEC_XZ_",I2.2,".dat")')tlabel
                             open(unit=tlabel+654,file=filename,status='unknown')
                          end if
                          if(alpha.eq.2.and.beta.eq.1)then
                             write(filename,'("SPEC_YX_",I2.2,".dat")')tlabel
                             open(unit=tlabel+654,file=filename,status='unknown')
                          end if
                          if(alpha.eq.2.and.beta.eq.2)then
                             write(filename,'("SPEC_YY_",I2.2,".dat")')tlabel
                             open(unit=tlabel+654,file=filename,status='unknown')
                          end if
                          if(alpha.eq.2.and.beta.eq.3)then
                             write(filename,'("SPEC_YZ_",I2.2,".dat")')tlabel
                             open(unit=tlabel+654,file=filename,status='unknown')
                          end if
                          if(alpha.eq.3.and.beta.eq.1)then
                             write(filename,'("SPEC_ZX_",I2.2,".dat")')tlabel
                             open(unit=tlabel+654,file=filename,status='unknown')
                          end if
                          if(alpha.eq.3.and.beta.eq.2)then
                             write(filename,'("SPEC_ZY_",I2.2,".dat")')tlabel
                             open(unit=tlabel+654,file=filename,status='unknown')
                          end if
                          if(alpha.eq.3.and.beta.eq.3)then
                             write(filename,'("SPEC_ZZ_",I2.2,".dat")')tlabel
                             open(unit=tlabel+654,file=filename,status='unknown')
                          end if
                       end if
                       do m=1,nfreq
                          omega=dfloat(8*m)/dfloat(nfreq)-2.0d0
                          do n=1,ifreq
                             omegai=omegai_min+omegai_max*dfloat(n)/dfloat(ifreq)
                             if(myrank.eq.0)write(654+tlabel,*)omega,omegai,spectra(ws,alpha,beta,wsp,n,m)
                          end do
                          if(myrank.eq.0)write(654+tlabel,*)
                       end do
                       if(myrank.eq.0)close(654+tlabel)
                    end do
                 end if
              end if
           end do
        end do
     end do

555  continue

     goto 666
     !
     ! here is the part where I construct sigma,pi cross sections for arbitrary q
     !
     if(myrank.eq.0)write(6,*)'Calculating the Milano spectra...'
     do l=1,4
        irixs=dfloat(0)
        ! l=1   sigma, sigma:
        ! l=2   sigma, pi:
        ! l=3   pi, sigma:
        ! l=4   pi, pi:
        do qqq=1,qstates
           ws=qlist(qqq)
           if(gsec(ws).eq.1)then
              do ppp=1,qstates
                 wsp=qlist(ppp)
                 tmom=mdiff(mtot(ws),mtot(wsp))
                 if(dlso.eq.0.0d0)then
                    tspin=sztot(ws)-sztot(wsp)
                    tlabel=10*tmom+abs(tspin)
                 else
                    tlabel=10*tmom
                 end if
                 if(tmom.eq.1.or.tmom.eq.2)then
                    pol_in=dfloat(0)
                    pol_out=dfloat(0)
                    if(l.eq.1)then
                       pol_in(2)=dfloat(1)
                       pol_out(2)=dfloat(1)
                    end if
                    if(l.eq.2)then
                       pol_in(2)=dfloat(1)
                       if(tmom.eq.1)then
                          pol_out(1)=0.9648d0
                          pol_out(3)=0.263d0
                       else
                          pol_out(1)=0.1946d0
                          pol_out(3)=0.9809d0
                       end if
                    end if
                    if(l.eq.3)then
                       pol_out(2)=dfloat(1)
                       if(tmom.eq.1)then
                          pol_in(1)=0.9648d0
                          pol_in(3)=-0.263d0
                       else
                          pol_in(1)=0.6655d0
                          pol_in(3)=0.7464d0
                       end if
                    end if
                    if(l.eq.4)then
                       if(tmom.eq.1)then
                          pol_in(1)=0.9648d0
                          pol_in(3)=-0.263d0
                          pol_out(1)=0.9648d0
                          pol_out(3)=0.263d0
                       else
                          pol_in(1)=0.6655d0
                          pol_in(3)=0.7464d0
                          pol_out(1)=0.1946d0
                          pol_out(3)=0.9809d0
                       end if
                    end if
                    do m=1,nfreq
                       omega=dfloat(8*m)/dfloat(nfreq)
                       gamma=0.05d0+0.05d0*(omega/2.5d0)**4.0d0
                       !                                   if(omega.gt.2.4d0)then
                       !                                      gamma=0.02d0+(omega-2.4d0)/dfloat(3)
                       !                                   else
                       !                                      gamma=0.02d0
                       !                                   end if
!                       gamma=0.1d0
                       do n=1,ifreq
                          omegai=omegai_min+omegai_max*dfloat(n)/dfloat(ifreq)
                          if(omegai.lt.-0.5d0)then
                             rixs=dfloat(0)
                             do j=1,gdegen(ws)
                                do i=1,min(nev,binstates(wsp))
                                   ksum=dcmplx(0.0d0,0.0d0)
                                   do alpha=1,3
                                      do beta=1,3
!                                         ksum=ksum+pol_in(beta)*pol_out(alpha)*ckern(alpha,beta,wsp,i,n,ws,j)
                                         ksum=ksum+pol_in(alpha)*pol_out(beta)*ckern(alpha,beta,wsp,i,n,ws,j)
                                      end do
                                   end do
                                   rixs=rixs+dreal(ksum*dconjg(ksum))*gamma/((omega-(deig(wsp,i)-a(1)))**2+gamma**2)
                                end do
                             end do
                             !                    if(dabs(omegai-XAS_maxx).lt.0.1d0.and.omegai.le.XAS_maxx)then
                             !                    if(dabs(omegai-XAS_maxx).lt.0.25d0)then
                             irixs(tmom,m)=irixs(tmom,m)+rixs
                             !                    write(7+k+5*l,*)omega,rixs,omegai-XAS_maxx,XAS_maxx,k
                          end if
                       end do
                    end do
                 end if
              end do
           end if
        end do
        do k=1,ncells
           if(myrank.eq.0)then
              if(l.eq.1)then
                 write(filename,'("Sigma_Sigma_",I2.2,".dat")')k
                 open(unit=7+k+5*l,file=filename,status='unknown')
              end if
              if(l.eq.2)then
                 write(filename,'("Sigma_Pi_",I2.2,".dat")')k
                 open(unit=7+k+5*l,file=filename,status='unknown')
              end if
              if(l.eq.3)then
                 write(filename,'("Pi_Sigma_",I2.2,".dat")')k
                 open(unit=7+k+5*l,file=filename,status='unknown')
              end if
              if(l.eq.4)then
                 write(filename,'("Pi_Pi_",I2.2,".dat")')k
                 open(unit=7+k+5*l,file=filename,status='unknown')
              end if
           end if
           do m=1,nfreq
              omega=dfloat(8*m)/dfloat(nfreq)
              if(myrank.eq.0)write(7+k+5*l,*)omega,1.5d0*irixs(k,m)/dfloat(10)
           end do
           if(myrank.eq.0)close(7+k+5*l)
        end do
     end do
  end if
     
666 continue
  
  call CPU_time(time6)
  if(myrank.eq.0)then
     write(6,*)'---------------------------------------------------'
     write(6,*)'---------------------------------------------------'
     write(6,*)'Succesful End'
     write(6,*)'Number of orbitals:',norbs
     write(6,*)'Number of holes:',nholes
     write(16,*)'---------------------------------------------------'
     write(16,*)'---------------------------------------------------'
     write(16,*)'Succesful End'
     write(16,*)'Number of orbitals:',norbs
     write(16,*)'Number of holes:',nholes
     if(flag_B.eq.1)then
        write(6,*)'Ncellsx=',ncellsx,'Ncellsy=',ncellsy
        write(16,*)'Ncellsx=',ncellsx,'Ncellsy=',ncellsy
     else
        write(6,*)'Betts cluster 8A'
        write(16,*)'Betts cluster 8A'
     end if
     write(6,*)'Number of twisted boundary conditions averaged over:',twistx*twisty
     write(6,*)'Total Hilbert space:',2*norbs,'choose',nholes,'=',ndim_max
     write(16,*)'Number of twisted boundary conditions averaged over:',twistx*twisty
     write(16,*)'Total Hilbert space:',2*norbs,'choose',nholes,'=',ndim_max
     write(6,*)'Valence spin-orbit coupling:',sngl(dlso)
     write(16,*)'Valence spin-orbit coupling:',sngl(dlso)
     if(dlso.eq.0.0d0)then
        write(6,*)'Sector:        total momentum [units pi/Nx,pi/Ny]   Total Sz spin       Number of States in Sector'
        write(16,*)'Sector:        total momentum [units pi/Nx,pi/Ny]   Total Sz spin       Number of States in Sector'
     else
        write(6,*)'Sector:        total momentum [units pi/Nx,pi/Ny]   Number of States in Sector'
        write(16,*)'Sector:        total momentum [units pi/Nx,pi/Ny]   Number of States in Sector'
     end if
     count=0
     do i=1,ncells
        if(dlso.eq.0.0d0)then
           do j=1,nholes+1
              count=count+1
              write(6,*)count,':  K = (',momx(i),momy(i),'),','     2Sz = ',-nholes+2*(j-1),':',dstates(i,j)
              write(16,*)count,':  K = (',momx(i),momy(i),'),','     2Sz = ',-nholes+2*(j-1),':',dstates(i,j)
           end do
        else
           write(6,*)i,':  K = (',momx(i),momy(i),'):',dstates(i,1)
           write(16,*)i,':  K = (',momx(i),momy(i),'):',dstates(i,1)
        end if
     end do
     write(6,*)'Projected Hilbert space, states used:'
     write(16,*)'Projected Hilbert space, states used:'
     if(dlso.eq.0.0d0)then
        do i=1,ncells*(nholes+1)
           if(usearpack.ge.2)write(6,*)i,binstates(i),binstates(i)
           if(usearpack.lt.2)write(6,*)i,binstates(i),NEV
           if(usearpack.ge.2)write(16,*)i,binstates(i),binstates(i)
           if(usearpack.lt.2)write(16,*)i,binstates(i),NEV
        end do
     else
        do i=1,ncells
           if(usearpack.ge.2)write(6,*)i,binstates(i),binstates(i)
           if(usearpack.lt.2)write(6,*)i,binstates(i),NEV
           if(usearpack.ge.2)write(16,*)i,binstates(i),binstates(i)
           if(usearpack.lt.2)write(16,*)i,binstates(i),NEV
        end do
     end if
     if(stepmax.eq.2)then
        write(6,*)'Core hole total Hilbert space: (',2*norbs,'choose',nholes-1,')times (6 times',ncells,') =',cdim_max
        write(6,*)'Core hole partitioned Hilbert space:',cstates(momgs)
        if(usearpack.eq.1.or.usearpack.eq.2)write(6,*)'CH eigenstates used',min(nevn,cstates(momgs))
        if(usearpack.eq.0.or.usearpack.gt.2)write(6,*)'CH eigenstates used',cstates(momgs)
        write(16,*)'Core hole total Hilbert space: (',2*norbs,'choose',nholes-1,')times (6 times',ncells,') =',cdim_max
        write(16,*)'Core hole partitioned Hilbert space:',cstates(momgs)
        if(usearpack.eq.1.or.usearpack.eq.2)write(16,*)'CH eigenstates used',min(nevn,cstates(momgs))
        if(usearpack.eq.0.or.usearpack.gt.2)write(16,*)'CH eigenstates used',cstates(momgs)
     end if
     if(stepmax.eq.3)then
        write(6,*)'Hole addition total Hilbert space:',hdim_max
        write(6,*)'Sector:        total momentum [units pi/Nx,pi/Ny]   Number of States in Sector'
        write(16,*)'Hole addition total Hilbert space:',hdim_max
        write(16,*)'Sector:        total momentum [units pi/Nx,pi/Ny]   Number of States in Sector'
        do i=1,ncells
           write(6,*)i,':  K = (',momx(i),momy(i),'),',':',hstates(i)
           write(16,*)i,':  K = (',momx(i),momy(i),'),',':',hstates(i)
        end do
     end if
     write(6,*)'---------------------------------------------------'
     write(6,*)'---------------------------------------------------'
     write(6,*)'U:',sngl(f0+4.0*f2+36.0*f4)
     write(6,*)'Uprime,J_H:'
     write(6,*)sngl(f0+4*f2-34*f4),sngl(35*f4),'for dx^2-y^2/dxy orbitals'
     write(6,*)sngl(f0-2*f2-4*f4),sngl(3*f2+20*f4),'for dx^2-y^2/dxz/yz orbitals'
     write(6,*)sngl(f0-4*f2+6*f4),sngl(4*f2+15*f4),'for dx^2-y^2/dz^2 orbitals'
     write(6,*)'Average U:',aveu
     write(6,*)
     write(6,*)'Condon-Shortly parameters used for 3d-3d [eV]:'
     write(6,*)'F0=',sngl(f0),'F_2=',sngl(f2),'F_4=',sngl(f4)
     write(6,*)
     write(6,*)'Condon-Shortly parameters used for 3d-2p [eV]:'
     write(6,*)'F0=',sngl(cf0),'F^2=',sngl(cf2),'G^1=',sngl(cf3),'G^3=',sngl(cf4)
     write(6,*)
     write(6,*)'spin exchange (4t^2/U) is',sngl(4*tpd*tpd/(f0+4*f2+36*f4))
     write(6,*)'---------------------------------------------------'
     write(6,*)'Total time for run:',time6-time1,'seconds.',(time6-time1)/60,'minutes'
     write(6,*)'time initializing:',time2-time1,'seconds.'
     write(6,*)'time for valence diagonalization:',time3-time2,'seconds.'
     if(stepmax.eq.2)then
        write(6,*)'time for core hole diagonalization:',time4-time3,'seconds.'
        write(6,*)'time for time evolution:',time5-time4,'seconds.'
        if(flag_R.ne.1)write(6,*)'time for XAS/RIXS:',time6-time5,'seconds.'
     end if
     
     write(16,*)'---------------------------------------------------'
     write(16,*)'---------------------------------------------------'
     write(16,*)'U:',sngl(f0+4.0*f2+36.0*f4)
     write(16,*)'Uprime,J_H:'
     write(16,*)sngl(f0+4*f2-34*f4),sngl(35*f4),'for dx^2-y^2/dxy orbitals'
     write(16,*)sngl(f0-2*f2-4*f4),sngl(3*f2+20*f4),'for dx^2-y^2/dxz/yz orbitals'
     write(16,*)sngl(f0-4*f2+6*f4),sngl(4*f2+15*f4),'for dx^2-y^2/dz^2 orbitals'
     write(16,*)'Average U:',aveu
     write(16,*)
     write(16,*)'Condon-Shortly parameters used for 3d-3d [eV]:'
     write(16,*)'F0=',sngl(f0),'F_2=',sngl(f2),'F_4=',sngl(f4)
     write(16,*)
     write(16,*)'Condon-Shortly parameters used for 3d-2p [eV]:'
     write(16,*)'F0=',sngl(cf0),'F^2=',sngl(cf2),'G^1=',sngl(cf3),'G^3=',sngl(cf4)
     write(16,*)
     write(16,*)'spin exchange (4t^2/U) is',sngl(4*tpd*tpd/(f0+4*f2+36*f4))
     write(16,*)'---------------------------------------------------'
     write(16,*)'Total time for run:',time6-time1,'seconds.',(time6-time1)/60,'minutes'
     write(16,*)'time initializing:',time2-time1,'seconds.'
     write(16,*)'time for valence diagonalization:',time3-time2,'seconds.'
     if(stepmax.eq.2)then
        write(16,*)'time for core hole diagonalization:',time4-time3,'seconds.'
        write(16,*)'time for time evolution:',time5-time4,'seconds.'
        if(flag_R.ne.1)write(16,*)'time for XAS/RIXS:',time6-time5,'seconds.'
     end if
     close(16)
  end if
  
999 continue
  call mpi_barrier(MPI_COMM_WORLD,mpistat)
  call mpi_finalize(mpistat)
        
  stop
end program main

recursive subroutine enumerate(balls,boxes,norbs,ndim_max,vect,box,state)
  implicit none
  integer, intent(in) :: balls,boxes,norbs,ndim_max
  integer i,count
  integer, intent(inout) :: vect(ndim_max,norbs),box(norbs),state
  
  if(boxes.le.0)then
     write(6,*)'error'
     return
  end if
  if(balls.gt.boxes)then
     return
  end if
  if(balls.lt.0)return
  if(boxes.eq.1)then
     box(1)=balls
     state=state+1
     count=0
     do i=1,norbs
        if(box(i).eq.1)then
           count=count+1
           vect(state,count)=i
        end if
     end do
     return
  else
     do i=0,1
        box(boxes)=i
        call enumerate(balls-i,boxes-1,norbs,ndim_max,vect,box,state)
     end do
  end if
  return
end subroutine enumerate

double precision function binomial(n,k)
  implicit none
  double precision factoid
  integer, intent(IN) :: n,k
  integer i
  factoid=1.0
  do i=1,k
     factoid=factoid*dfloat(n-k+i)/dfloat(i)
  end do
  binomial=factoid
end function binomial

integer function fact(n)
  implicit none
  integer :: p,i
  integer, intent(IN) :: n
  p=1
  do i=1,n
     p=p*i
  end do
  fact=p
end function fact

double precision function gaunt(l1,m1,l2,m2,kk)
  implicit none
  integer :: l1,m1,l2,m2,kk

  gaunt=dfloat(0)

  if(abs(m1).gt.l1)goto 911
  if(abs(m2).gt.l2)goto 911

  if(l1.eq.2.and.l2.eq.2)then
     if(kk.eq.0)then
        if(m1*m2.eq.4)gaunt=dfloat(1)
        if(m1*m2.eq.1)gaunt=dfloat(1)
        if(m1.eq.0.and.m2.eq.0)gaunt=dfloat(1)
     end if
     if(kk.eq.2)then
        if(m1*m2.eq.4)gaunt=-dfloat(2)/dfloat(7)
        if(m1*m2.eq.2)then
           if(abs(m1).gt.abs(m2))gaunt=dsqrt(dfloat(6))/dfloat(7)
           if(abs(m2).gt.abs(m1))gaunt=-dsqrt(dfloat(6))/dfloat(7)
        end if
        if(abs(m1).eq.2.and.m2.eq.0)gaunt=-dfloat(2)/dfloat(7)
        if(abs(m2).eq.2.and.m1.eq.0)gaunt=-dfloat(2)/dfloat(7)
        if(m1*m2.eq.1)gaunt=dfloat(1)/dfloat(7)
        if(abs(m1).eq.1.and.m2.eq.0)gaunt=dfloat(1)/dfloat(7)
        if(abs(m2).eq.1.and.m1.eq.0)gaunt=-dfloat(1)/dfloat(7)
        if(m1.eq.0.and.m2.eq.0)gaunt=dfloat(2)/dfloat(7)
        if(m1*m2.eq.-1)gaunt=-dsqrt(dfloat(6))/dfloat(7)
     end if
     if(kk.eq.4)then
        if(m1*m2.eq.4)gaunt=dfloat(1)/dfloat(21)
        if(m1*m2.eq.2)then
           if(abs(m1).gt.abs(m2))gaunt=-dsqrt(dfloat(5))/dfloat(21)
           if(abs(m2).gt.abs(m1))gaunt=dsqrt(dfloat(5))/dfloat(21)
        end if
        if(abs(m1).eq.2.and.m2.eq.0)gaunt=dsqrt(dfloat(15))/dfloat(21)
        if(abs(m2).eq.2.and.m1.eq.0)gaunt=dsqrt(dfloat(15))/dfloat(21)
        if(m1*m2.eq.1)gaunt=-dfloat(4)/dfloat(21)
        if(abs(m1).eq.1.and.m2.eq.0)gaunt=dsqrt(dfloat(30))/dfloat(21)
        if(abs(m2).eq.1.and.m1.eq.0)gaunt=-dsqrt(dfloat(30))/dfloat(21)
        if(m1.eq.0.and.m2.eq.0)gaunt=dfloat(6)/dfloat(21)
        if(m1*m2.eq.-4)gaunt=dsqrt(dfloat(70))/dfloat(21)
        if(m1*m2.eq.-2)then
           if(abs(m1).gt.abs(m2))gaunt=-dsqrt(dfloat(35))/dfloat(21)
           if(abs(m2).gt.abs(m1))gaunt=dsqrt(dfloat(35))/dfloat(21)
        end if
        if(m1*m2.eq.-1)gaunt=-dsqrt(dfloat(40))/dfloat(21)
     end if
  end if
  if(l1.eq.1.and.l2.eq.1)then
     if(kk.eq.0)then
        if(m1*m2.eq.1)gaunt=dfloat(1)
        if(m1.eq.0.and.m2.eq.0)gaunt=dfloat(1)
     end if
     if(kk.eq.2)then
        if(m1*m2.eq.1)gaunt=-dfloat(1)/dfloat(5)
        if(abs(m1).eq.1.and.m2.eq.0)gaunt=dsqrt(dfloat(3))/dfloat(5)
        if(abs(m2).eq.1.and.m1.eq.0)gaunt=-dsqrt(dfloat(3))/dfloat(5)
        if(m1.eq.0.and.m2.eq.0)gaunt=dfloat(2)/dfloat(5)
        if(m1*m2.eq.-1)gaunt=-dsqrt(dfloat(6))/dfloat(5)
     end if
  end if
  if(l1.eq.2.and.l2.eq.1)then
     if(abs(m2).le.l2)then
        if(kk.eq.1)then
           if(m1*m2.eq.2)gaunt=dsqrt(dfloat(6))/dsqrt(dfloat(15))
           if(m1*m2.eq.1)gaunt=dsqrt(dfloat(3))/dsqrt(dfloat(15))
           if(m1.eq.0.and.abs(m2).eq.1)gaunt=dfloat(1)/dsqrt(dfloat(15))
           if(m2.eq.0.and.abs(m1).eq.1)gaunt=dsqrt(dfloat(3))/dsqrt(dfloat(15))
           if(m2.eq.0.and.m1.eq.0)gaunt=dfloat(2)/dsqrt(dfloat(15))
        end if
        if(kk.eq.3)then
           if(m1*m2.eq.2)gaunt=-dsqrt(dfloat(3))/dsqrt(dfloat(245))
           if(m1*m2.eq.1)gaunt=-dfloat(3)/dsqrt(dfloat(245))
           if(m1.eq.0.and.abs(m2).eq.1)gaunt=-dsqrt(dfloat(18))/dsqrt(dfloat(245))
           if(abs(m1).eq.2.and.m2.eq.0)gaunt=dsqrt(dfloat(15))/dsqrt(dfloat(245))
           if(abs(m1).eq.1.and.m2.eq.0)gaunt=dsqrt(dfloat(24))/dsqrt(dfloat(245))
           if(m1.eq.0.and.m2.eq.0)gaunt=dsqrt(dfloat(27))/dsqrt(dfloat(245))
           if(m1*m2.eq.-2)gaunt=-dsqrt(dfloat(45))/dsqrt(dfloat(245))
           if(m1*m2.eq.-1)gaunt=-dsqrt(dfloat(30))/dsqrt(dfloat(245))
        end if
     end if
  end if
  if(l1.eq.1.and.l2.eq.2)then
     if(abs(m1).le.l1)then
        if(kk.eq.1)then
           if(m1*m2.eq.2)gaunt=-dsqrt(dfloat(6))/dsqrt(dfloat(15))
           if(m1*m2.eq.1)gaunt=dsqrt(dfloat(3))/dsqrt(dfloat(15))
           if(m2.eq.0.and.abs(m1).eq.1)gaunt=-dfloat(1)/dsqrt(dfloat(15))
           if(m1.eq.0.and.abs(m2).eq.1)gaunt=-dsqrt(dfloat(3))/dsqrt(dfloat(15))
           if(m2.eq.0.and.m1.eq.0)gaunt=dfloat(2)/dsqrt(dfloat(15))
        end if
        if(kk.eq.3)then
           if(m1*m2.eq.2)gaunt=dsqrt(dfloat(3))/dsqrt(dfloat(245))
           if(m1*m2.eq.1)gaunt=-dfloat(3)/dsqrt(dfloat(245))
           if(m2.eq.0.and.abs(m1).eq.1)gaunt=dsqrt(dfloat(18))/dsqrt(dfloat(245))
           if(abs(m2).eq.2.and.m1.eq.0)gaunt=dsqrt(dfloat(15))/dsqrt(dfloat(245))
           if(abs(m2).eq.1.and.m1.eq.0)gaunt=-dsqrt(dfloat(24))/dsqrt(dfloat(245))
           if(m1.eq.0.and.m2.eq.0)gaunt=dsqrt(dfloat(27))/dsqrt(dfloat(245))
           if(m1*m2.eq.-2)gaunt=dsqrt(dfloat(45))/dsqrt(dfloat(245))
           if(m1*m2.eq.-1)gaunt=-dsqrt(dfloat(30))/dsqrt(dfloat(245))
        end if
     end if
  end if
  if(l1.eq.0.and.l2.eq.0)then
     if(kk.eq.0)gaunt=dfloat(1)
  end if
  if(l1.eq.1.and.l2.eq.0)then
     if(kk.eq.1)then
        if(m1.eq.0)gaunt=-dfloat(1)/dsqrt(3.0d0)
        if(m1.eq.1.or.m1.eq.-1)gaunt=dfloat(1)/dsqrt(3.0d0)
     end if
  end if
  if(l1.eq.0.and.l2.eq.1)then
     if(kk.eq.1)then
        if(m2.eq.0)gaunt=dfloat(1)/dsqrt(3.0d0)
        if(m2.eq.1.or.m2.eq.-1)gaunt=-dfloat(1)/dsqrt(3.0d0)
     end if
  end if
911 continue
  return
end function gaunt

function Ylm(L,M,the,phi)
! spherical harmonic function for 0 <= L <= 4, -L <= M <= L                                                                  
  implicit none
  complex:: Ylm            ! spherical harmonic function                                                                     
  complex:: eye
  integer, intent(in):: L,M  ! indices                                                                                       
  real, intent(in):: the ! polar angle, theta (angle w/ positive z-axis)                                                     
  real, intent(in):: phi ! azimuthal angle, phi (angle of projection on x-y plane w/ pos. x-axis)                            
  real :: pi

  pi=2.0*asin(1.0)

  eye=cmplx(0.0,1.0)

  select case(L)
  case(0)
     select case(M)
     case(0)
        Ylm=0.5/sqrt(pi)
     case default
!        call stop_error("Ylm error: invalid M")                                                                             
     end select
  case(1)
     select case(M)
     case(0)
        Ylm=sqrt(3/pi)*cos(the)/2
     case(1,-1)
        Ylm=-exp(eye*phi)*sqrt(3/(2*pi))*sin(the)/2
     case default
!        call stop_error("Ylm error: invalid M")                                                                             
     end select
  case(2)
     select case(M)
     case(0)
        Ylm=sqrt(5/pi)*(-1 + 3*cos(the)**2)/4
     case(1,-1)
        Ylm=-exp(eye*phi)*sqrt(15/(2*pi))*cos(the)*sin(the)/2
     case(2,-2)
        Ylm=exp(2*eye*phi)*sqrt(15/(2*pi))*sin(the)**2/4
     case default
!        call stop_error("Ylm error: invalid M")                                                                             
     end select
  case(3)
     select case(M)
     case(0)
        Ylm=sqrt(7/pi)*(-3*cos(the) + 5*cos(the)**3)/4
     case(1,-1)
        Ylm=-exp(eye*phi)*sqrt(21/pi)*(-1 + 5*cos(the)**2)*sin(the)/8
     case(2,-2)
        Ylm=exp(2*eye*phi)*sqrt(105/(2*pi))*cos(the)*sin(the)**2/4
     case(3,-3)
        Ylm=-exp(3*eye*phi)*sqrt(35/pi)*sin(the)**3/8
     case default
!        call stop_error("Ylm error: invalid M")                                                                             
     end select
  case(4)
     select case(M)
     case(0)
        Ylm=3*sqrt(1/pi)*(35*cos(the)**4 - 30*cos(the)**2+3)/16
     case(1,-1)
        Ylm=3*exp(eye*phi)*sqrt(5/pi)*(-3*cos(the) + 7*cos(the)**3)*sin(the)/8
     case(2,-2)
        Ylm=3*exp(2*eye*phi)*sqrt(5/(2*pi))*(7*cos(the)**2-1)*sin(the)**2/8
     case(3,-3)
        Ylm=3*exp(3*eye*phi)*sqrt(35/pi)*cos(the)*sin(the)**3/8
     case(4,-4)
        Ylm=3*exp(4*eye*phi)*sqrt(35/(2*pi))*sin(the)**4/16
     case default
!        call stop_error("Ylm error: invalid M")                                                                             
     end select
  case default
!     call stop_error("Ylm error: invalid L")                                                                                
  end select
  if (M < 0) Ylm=(-1)**(-M)*conjg(Ylm)
end function Ylm

subroutine getfsign(N,fsign)
  implicit none
  integer :: N,i,j,k,fsign(N*(N-1)/2)
  k=1
  do i=N-1,1,-1
     do j=1,i
        fsign(k)=(-1)**(1+j)
        k=k+1
     end do
  end do
  return
end subroutine getfsign

subroutine getgammakm(Q,x,y,z,nc,gammakm)
  implicit none
  complex :: Ylm,sum,gammakm(4,-4:4)
  real :: Qave,Rave,phi,theta,pi,rrr,small
  integer :: nc,i,k,m
  real :: Q(nc),x(nc),y(nc),z(nc)

  gammakm=0.0

  Qave=dfloat(0)
  Rave=dfloat(0)
  do i=1,nc
     rrr=sqrt(x(i)**2+y(i)**2+z(i)**2)
     Qave=Qave+Q(i)/float(nc)
     Rave=Rave+rrr/float(nc)
  end do

  small=0.0000001

  pi=2.0*asin(1.0)
  do k=1,4
     do m=-k,k
        sum=cmplx(0)
        do i=1,nc
           rrr=sqrt(x(i)**2+y(i)**2+z(i)**2)
           theta=acos(z(i)/rrr)
           if(x(i).eq.0.0)then
              phi=pi/2.0
           else
              phi=atan(y(i)/x(i))
           end if
           sum=sum+sqrt(4*pi/(2*k+1))*(Q(i)/Qave)*Ylm(k,m,theta,phi)*(Rave/rrr)**(1+k)
        end do
        if(cabs(sum).lt.small)sum=cmplx(0.0,0.0)
        gammakm(k,m)=sum
     end do
  end do
  return
end subroutine getgammakm

subroutine getpgsymm(symm,pgsymm)
  implicit none
  double precision, parameter :: small=0.0001d0
  integer :: symm,nc,i,j,k
  real, allocatable, dimension(:) :: Q,x,y,z
  double precision :: gaunt,pgsymm(5,5)
  complex :: gammakm(4,-4:4),ccsum

  write(6,*)'point group symmetry'                                                               
  write(6,*)'1 for octahedral, 2 for square planar, 3 for square pyramidal, 4 to enter generic ligand valence and generic position around transition metal'      
  write(6,*)'point group symmetry =',symm
  write(6,*)
  
  if(symm.eq.4)then
     write(6,*)'Enter number of nearest ligands to get xtal field from point charge model'
     read(5,*)nc

     allocate(Q(nc))
     allocate(x(nc))
     allocate(y(nc))
     allocate(z(nc))

     do i=1,nc
        write(6,*)'enter valence charge'
        read(5,*)Q(i)
        write(6,*)'position from transition metal ion (x,y,z)'
        read(5,*)x(i),y(i),z(i)
     end do
  end if

  if(symm.eq.1)then
     nc=6
     allocate(Q(nc))
     allocate(x(nc))
     allocate(y(nc))
     allocate(z(nc))
     do i=1,nc
        Q(i)=-2.0
     end do
     x=0.0
     y=0.0
     z=0.0
     x(1)=1
     y(2)=1
     z(3)=1
     x(4)=-1
     y(5)=-1
     z(6)=-1
  end if

  if(symm.eq.2)then
     nc=4
     allocate(Q(nc))
     allocate(x(nc))
     allocate(y(nc))
     allocate(z(nc))
     do i=1,nc
        Q(i)=-2.0
     end do
     x=0.0
     y=0.0
     z=0.0
     x(1)=1
     y(2)=1
     x(3)=-1
     y(4)=-1
  end if

  if(symm.eq.3)then
     nc=5
     allocate(Q(nc))
     allocate(x(nc))
     allocate(y(nc))
     allocate(z(nc))
     do i=1,nc
        Q(i)=-2.0
     end do
     x=0.0
     y=0.0
     z=0.0
     x(1)=1
     y(2)=1
     x(3)=-1
     y(4)=-1
     z(5)=1
  end if

  call getgammakm(Q,x,y,z,nc,gammakm)

  pgsymm=0.0d0

  do i=-2,2
     do j=-2,2
        ccsum=cmplx(0,0)
        do k=2,4,2
           if(k.eq.3)then
              ccsum=ccsum+gammakm(k+1,i-j)*gaunt(2,i,2,j,k)
           else
              ccsum=ccsum+gammakm(k,i-j)*gaunt(2,i,2,j,k)
           end if
        end do
        if(abs(real(ccsum)).gt.small)pgsymm(i+3,j+3)=real(ccsum)
     end do
  end do
  return
end subroutine getpgsymm

subroutine gethint(f0,f2,f4,cf0,cf2,cf3,cf4,fp,fp2,gpf,gpg,hint)
  implicit none
  double precision :: sum,gaunt,f0,f2,f4,cf0,cf2,cf3,cf4,fp,fp2,gpf,gpg,hint(0:2,-2:2,0:2,-2:2,0:2,-2:2,0:2,-2:2),rk(0:4,0:2,0:2,0:2,0:2)
  integer :: l1,l2,l3,l4,i,j,k,l,r
  hint=dfloat(0)
  rk=dfloat(0)
  rk(0,1,1,1,1)=fp     ! U_pp F0
  rk(1,1,0,0,1)=gpf     ! U_(1s-2p) L core F0
  rk(1,0,1,1,0)=gpf     ! U_(1s-2p) L core F0
  rk(1,1,0,1,0)=gpg     ! U_(1s-2p) L core G1
  rk(1,0,1,0,1)=gpg     ! U_(1s-2p) L core G1
  rk(2,1,1,1,1)=fp2    ! U_pp F2
  rk(0,2,2,2,2)=f0     ! U_dd
  rk(2,2,2,2,2)=f2*dfloat(49)   ! U_dd exchange
  rk(4,2,2,2,2)=f4*dfloat(441)  ! U_dd exchange
  rk(0,1,2,2,1)=cf0    ! U_pd Tm core F_0
  rk(0,2,1,1,2)=cf0    ! U_pd Tm core F_0
  rk(2,1,2,2,1)=cf2    ! U_pd Tm core F_2
  rk(2,2,1,1,2)=cf2    ! U_pd Tm core F_2
  rk(1,1,2,1,2)=cf3    ! U_pd Tm core G_1
  rk(1,2,1,2,1)=cf3    ! U_pd Tm core G_1
  rk(3,1,2,1,2)=cf4    ! U_pd Tm core G_3
  rk(3,2,1,2,1)=cf4    ! U_pd Tm core G_3
  do l1=0,2
     do l2=0,2
        do l3=0,2
           do l4=0,2
              do i=-l1,l1
                 do j=-l2,l2
                    do k=-l3,l3
                       do l=-l4,l4
                          sum=dfloat(0)
                          if(i+j.eq.k+l)then
                             do r=0,4
                                sum=sum+gaunt(l1,i,l4,l,r)*gaunt(l3,k,l2,j,r)*rk(r,l1,l2,l3,l4)
                             end do
                          end if
                          hint(l1,i,l2,j,l3,k,l4,l)=sum
                       end do
                    end do
                 end do
              end do
           end do
        end do
     end do
  end do
  return
end subroutine gethint

subroutine getgint(f0,f2,f4,norbspc,gint)
  implicit none
  integer :: norbspc,i,j,k,l,r,pp
  double precision :: sum,gaunt,arr(0:4),f0,f2,f4
  double precision :: gint(norbspc,norbspc,norbspc,norbspc)

  arr=dfloat(0)
  arr(0)=f0
  arr(2)=f2*dfloat(7*7)
  arr(4)=f4*dfloat(21*21)
  gint=dfloat(0)
  do i=-2,2
     do j=-2,2
        do k=-2,2
           do l=-2,2
              sum=dfloat(0)
              if(i+j.eq.k+l)then
                 do pp=0,4,2
                    sum=sum+gaunt(2,i,2,l,pp)*gaunt(2,k,2,j,pp)*arr(pp)
                 end do
              end if
              gint(i+3,j+3,k+3,l+3)=sum
           end do
        end do
     end do
  end do
  return
end subroutine getgint

subroutine checkint(u1,u2,u3,u4,norbspc,norbstm,ncore,outout)
  implicit none
  integer :: u1,u2,u3,u4,norbspc,outout,norbstm,ncore,pout
  !make sure that we do not include interaction terms that are unphysical such as core TM-ligand interaction                            
  !make sure that only interactions on the same atom are considered                                                                     
  outout=1
  if(u1.le.norbstm)then   ! u1 is TM                                                                                                    
     if(u2.gt.norbstm.and.u2.le.norbspc)outout=0  ! u2 is a ligand                                                                      
     if(u3.gt.norbstm.and.u3.le.norbspc)outout=0  ! u3 is a ligand                                                                      
     if(u4.gt.norbstm.and.u4.le.norbspc)outout=0  ! u4 is a ligand                                                                      
  end if
  if(u2.le.norbstm)then   ! u2 is TM                                                                                                    
     if(u1.gt.norbstm.and.u1.le.norbspc)outout=0  ! u1 is a ligand                                                                      
     if(u3.gt.norbstm.and.u3.le.norbspc)outout=0  ! u3 is a ligand                                                                      
     if(u4.gt.norbstm.and.u4.le.norbspc)outout=0  ! u4 is a ligand                                                                      
  end if
  if(u3.le.norbstm)then   ! u3 is TM                                                                                                    
     if(u2.gt.norbstm.and.u2.le.norbspc)outout=0  ! u2 is a ligand                                                                      
     if(u1.gt.norbstm.and.u1.le.norbspc)outout=0  ! u1 is a ligand                                                                      
     if(u4.gt.norbstm.and.u4.le.norbspc)outout=0  ! u4 is a ligand                                                                      
  end if
  if(u4.le.norbstm)then   ! u4 is TM                                                                                                    
     if(u2.gt.norbstm.and.u2.le.norbspc)outout=0  ! u2 is a ligand                                                                      
     if(u3.gt.norbstm.and.u3.le.norbspc)outout=0  ! u3 is a ligand                                                                      
     if(u1.gt.norbstm.and.u1.le.norbspc)outout=0  ! u1 is a ligand                                                                      
  end if
  if(u1.gt.norbstm.and.u1.le.norbspc-3)then   ! u1 is ox ligand                                                                         
     pout=1
     if(u2.gt.norbstm.and.u2.le.norbspc-3)pout=pout+1  ! u2 is ox ligand                                                                
     if(u3.gt.norbstm.and.u3.le.norbspc-3)pout=pout+1  ! u3 is ox ligand                                                                
     if(u4.gt.norbstm.and.u4.le.norbspc-3)pout=pout+1  ! u4 is ox ligand                                                                
     if(pout.ne.4)outout=0
  end if
  if(u1.gt.norbspc-3.and.u1.le.norbspc)then   ! u1 is oy ligand                                                                         
     pout=1
     if(u2.gt.norbspc-3.and.u2.le.norbspc)pout=pout+1  ! u2 is oy ligand                                                                
     if(u3.gt.norbspc-3.and.u3.le.norbspc)pout=pout+1  ! u3 is oy ligand                                                                
     if(u4.gt.norbspc-3.and.u4.le.norbspc)pout=pout+1  ! u4 is oy ligand                                                                
     if(pout.ne.4)outout=0
  end if
  if(ncore.eq.6)then
     if(u1.gt.norbspc)then  ! u1 is a TM core                                                                                              
        if(u2.gt.norbstm.and.u2.lt.norbspc)outout=0 ! u2 is a ligand                                                                    
        if(u3.gt.norbstm.and.u3.lt.norbspc)outout=0 ! u3 is a ligand                                                                    
        if(u4.gt.norbstm.and.u4.lt.norbspc)outout=0 ! u4 is a ligand                                                                    
     end if
     if(u2.gt.norbspc)then  ! u2 is a TM core                                                                                              
        if(u1.gt.norbstm.and.u1.lt.norbspc)outout=0 ! u1 is a ligand                                                                    
        if(u3.gt.norbstm.and.u3.lt.norbspc)outout=0 ! u3 is a ligand                                                                    
        if(u4.gt.norbstm.and.u4.lt.norbspc)outout=0 ! u4 is a ligand                                                                    
     end if
     if(u3.gt.norbspc)then  ! u3 is a TM core                                                                                              
        if(u1.gt.norbstm.and.u1.lt.norbspc)outout=0 ! u1 is a ligand                                                                    
        if(u2.gt.norbstm.and.u2.lt.norbspc)outout=0 ! u2 is a ligand                                                                    
        if(u4.gt.norbstm.and.u4.lt.norbspc)outout=0 ! u4 is a ligand                                                                    
     end if
     if(u4.gt.norbspc)then  ! u4 is a TM core                                                                                              
        if(u1.gt.norbstm.and.u1.lt.norbspc)outout=0 ! u1 is a ligand                                                                    
        if(u2.gt.norbstm.and.u2.lt.norbspc)outout=0 ! u2 is a ligand                                                                    
        if(u3.gt.norbstm.and.u3.lt.norbspc)outout=0 ! u3 is a ligand                                                                    
     end if
  end if
  if(ncore.eq.4)then
     if(u1.gt.norbspc)then  ! u1 is a ligand core                                                                                              
        if(u2.le.norbstm)outout=0 ! u2 is a TM                                                                                          
        if(u3.le.norbstm)outout=0 ! u3 is a TM                                                                                          
        if(u4.le.norbstm)outout=0 ! u4 is a TM                                                                                          
        if(u1.gt.norbspc.and.u1.le.norbspc+2)then !u1 is ox ligand core
           pout=1
           if((u2.gt.norbstm.and.u2.le.norbspc-3).or.(u2.gt.norbspc.and.u2.le.norbspc+2))pout=pout+1 ! u2 is ox
           if((u3.gt.norbstm.and.u3.le.norbspc-3).or.(u3.gt.norbspc.and.u3.le.norbspc+2))pout=pout+1 ! u3 is ox
           if((u4.gt.norbstm.and.u4.le.norbspc-3).or.(u4.gt.norbspc.and.u4.le.norbspc+2))pout=pout+1 ! u4 is ox
           if(pout.ne.4)outout=0
        end if
        if(u1.gt.norbspc+2)then !u1 is oy ligand core
           pout=1
           if((u2.gt.norbspc-3.and.u2.le.norbspc).or.(u2.gt.norbspc+2))pout=pout+1 ! u2 is oy 
           if((u3.gt.norbspc-3.and.u3.le.norbspc).or.(u3.gt.norbspc+2))pout=pout+1 ! u3 is oy 
           if((u4.gt.norbspc-3.and.u4.le.norbspc).or.(u4.gt.norbspc+2))pout=pout+1 ! u4 is oy 
           if(pout.ne.4)outout=0
        end if
     end if
     if(u2.gt.norbspc)then  ! u2 is a ligand core                                                                                              
        if(u1.le.norbstm)outout=0 ! u1 is a TM                                                                                          
        if(u3.le.norbstm)outout=0 ! u3 is a TM                                                                                          
        if(u4.le.norbstm)outout=0 ! u4 is a TM                                                                                          
        if(u2.gt.norbspc.and.u2.le.norbspc+2)then !u2 is ox ligand core
           pout=1
           if((u1.gt.norbstm.and.u1.le.norbspc-3).or.(u1.gt.norbspc.and.u1.le.norbspc+2))pout=pout+1 ! u1 is ox
           if((u3.gt.norbstm.and.u3.le.norbspc-3).or.(u3.gt.norbspc.and.u3.le.norbspc+2))pout=pout+1 ! u3 is ox
           if((u4.gt.norbstm.and.u4.le.norbspc-3).or.(u4.gt.norbspc.and.u4.le.norbspc+2))pout=pout+1 ! u4 is ox
           if(pout.ne.4)outout=0
        end if
        if(u2.gt.norbspc+2)then !u2 is oy ligand core
           pout=1
           if((u1.gt.norbspc-3.and.u1.le.norbspc).or.(u2.gt.norbspc+2))pout=pout+1 ! u1 is oy 
           if((u3.gt.norbspc-3.and.u3.le.norbspc).or.(u3.gt.norbspc+2))pout=pout+1 ! u3 is oy 
           if((u4.gt.norbspc-3.and.u4.le.norbspc).or.(u4.gt.norbspc+2))pout=pout+1 ! u4 is oy 
           if(pout.ne.4)outout=0
        end if
     end if
     if(u3.gt.norbspc)then  ! u3 is a ligand core                                                                                              
        if(u1.le.norbstm)outout=0 ! u1 is a TM                                                                                          
        if(u2.le.norbstm)outout=0 ! u2 is a TM                                                                                          
        if(u4.le.norbstm)outout=0 ! u4 is a TM                                                                                          
        if(u3.gt.norbspc.and.u3.le.norbspc+2)then !u3 is ox ligand core
           pout=1
           if((u2.gt.norbstm.and.u2.le.norbspc-3).or.(u2.gt.norbspc.and.u2.le.norbspc+2))pout=pout+1 ! u2 is ox
           if((u1.gt.norbstm.and.u1.le.norbspc-3).or.(u1.gt.norbspc.and.u1.le.norbspc+2))pout=pout+1 ! u1 is ox
           if((u4.gt.norbstm.and.u4.le.norbspc-3).or.(u4.gt.norbspc.and.u4.le.norbspc+2))pout=pout+1 ! u4 is ox
           if(pout.ne.4)outout=0
        end if
        if(u3.gt.norbspc+2)then !u3 is oy ligand core
           pout=1
           if((u2.gt.norbspc-3.and.u2.le.norbspc).or.(u2.gt.norbspc+2))pout=pout+1 ! u2 is oy 
           if((u1.gt.norbspc-3.and.u1.le.norbspc).or.(u1.gt.norbspc+2))pout=pout+1 ! u1 is oy 
           if((u4.gt.norbspc-3.and.u4.le.norbspc).or.(u4.gt.norbspc+2))pout=pout+1 ! u4 is oy 
           if(pout.ne.4)outout=0
        end if
     end if
     if(u4.gt.norbspc)then  ! u4 is a ligand core                                                                                              
        if(u1.le.norbstm)outout=0 ! u1 is a TM                                                                                          
        if(u2.le.norbstm)outout=0 ! u2 is a TM                                                                                          
        if(u3.le.norbstm)outout=0 ! u3 is a TM                                                                                          
        if(u4.gt.norbspc.and.u4.le.norbspc+2)then !u4 is ox ligand core
           pout=1
           if((u1.gt.norbstm.and.u1.le.norbspc-3).or.(u1.gt.norbspc.and.u1.le.norbspc+2))pout=pout+1 ! u1 is ox
           if((u2.gt.norbstm.and.u2.le.norbspc-3).or.(u2.gt.norbspc.and.u2.le.norbspc+2))pout=pout+1 ! u2 is ox
           if((u3.gt.norbstm.and.u3.le.norbspc-3).or.(u3.gt.norbspc.and.u3.le.norbspc+2))pout=pout+1 ! u3 is ox
           if(pout.ne.4)outout=0
        end if
        if(u4.gt.norbspc+2)then !u4 is oy ligand core
           pout=1
           if((u1.gt.norbspc-3.and.u1.le.norbspc).or.(u1.gt.norbspc+2))pout=pout+1 ! u1 is oy 
           if((u2.gt.norbspc-3.and.u2.le.norbspc).or.(u2.gt.norbspc+2))pout=pout+1 ! u2 is oy 
           if((u3.gt.norbspc-3.and.u3.le.norbspc).or.(u3.gt.norbspc+2))pout=pout+1 ! u3 is oy 
           if(pout.ne.4)outout=0
        end if
     end if
  end if
  return
end subroutine checkint
  
subroutine getcint(f0,f2,f4,cf0,cf2,cf3,cf4,norbspc,cint)
  implicit none
  double precision :: sum,gaunt,cf0,cf2,cf3,cf4,arrc1(0:2),arrc2(3),arr(0:4),f0,f2,f4
  integer :: norbspc,l1,l2,l3,l4,i,j,k,l,r
  double precision :: cint(norbspc+3,norbspc+3,norbspc+3,norbspc+3)

  arr=dfloat(0)
  arr(0)=f0
  arr(2)=f2*dfloat(7*7)
  arr(4)=f4*dfloat(21*21)

  arrc1=dfloat(0)
  arrc1(0)=cf0
  arrc1(2)=cf2
  arrc2=dfloat(0)
  arrc2(1)=cf3
  arrc2(3)=cf4

  cint=dfloat(0)
  
  do l1=1,2
     do l2=1,2
        do l3=1,2
           do l4=1,2
              do i=-l1,l1
                 do j=-l2,l2
                    do k=-l3,l3
                       do l=-l4,l4
                          sum=dfloat(0)
                          if(i+j.eq.k+l)then
                             if(l1+l2+l3+l4.ne.8)then
                                do r=0,2,2
                                   if(l2.eq.l3.and.l1.eq.l4)sum=sum+gaunt(l1,i,l4,l,r)*gaunt(l3,k,l2,j,r)*arrc1(r)
                                   if(l2.eq.l4.and.l1.eq.l3)sum=sum+gaunt(l1,i,l4,l,r+1)*gaunt(l3,k,l2,j,r+1)*arrc2(r+1)
                                end do
                             else
                                do r=0,4,2
                                   sum=sum+gaunt(l1,i,l4,l,r)*gaunt(l3,k,l2,j,r)*arr(r)
                                end do
                             end if
                          end if
                          if(l1+l2+l3+l4.eq.8)then
                             cint(i+3,j+3,k+3,l+3)=sum
                          end if
                          if(l1+l2+l3+l4.eq.4)then
                             cint(i+norbspc+2,j+norbspc+2,k+norbspc+2,l+norbspc+2)=sum
                          end if
                          if(l1+l3.eq.2.and.l2+l4.ne.2)then
                             cint(i+norbspc+2,j+3,k+norbspc+2,l+3)=sum
                          end if
                          if(l1+l4.eq.2.and.l2+l3.ne.2)then
                             cint(i+norbspc+2,j+3,k+3,l+norbspc+2)=sum
                          end if
                          if(l2+l3.eq.2.and.l1+l4.ne.2)then
                             cint(i+3,j+norbspc+2,k+norbspc+2,l+3)=sum
                          end if
                          if(l2+l4.eq.2.and.l1+l3.ne.2)then
                             cint(i+3,j+norbspc+2,k+3,l+norbspc+2)=sum
                          end if
                       end do
                    end do
                 end do
              end do
           end do
        end do
     end do
  end do
  return
end subroutine getcint

subroutine setup_s(ncellsx,ncellsy,ncells,momx,momy,mdiff,madd)
  implicit none
  integer :: ncellsx,ncellsy,ncells,count,i,j,ii
  integer :: momx(ncells),momy(ncells)
  integer :: mdiff(ncells,ncells),madd(ncells,ncells)
  integer :: kkx,kky,ppx,ppy,diffx,diffy,addx,addy
  
  count=1
  do i=0,ncellsx-1
     do j=0,ncellsy-1
        momx(count)=i
        momy(count)=j
        count=count+1
     end do
  end do

  if(ncellsx.eq.1)momx(1)=0
  if(ncellsy.eq.1)momy(1)=0

  do i=1,ncells
     kkx=momx(i)
     kky=momy(i)
     do j=1,ncells
        ppx=momx(j)
        ppy=momy(j)
        diffx=kkx-ppx
        if(diffx.lt.0)diffx=diffx+ncellsx
        diffy=kky-ppy
        if(diffy.lt.0)diffy=diffy+ncellsy
        addx=kkx+ppx
        if(addx.ge.ncellsx)addx=addx-ncellsx
        addy=kky+ppy
        if(addy.ge.ncellsy)addy=addy-ncellsy
        do ii=1,ncells
           if(momx(ii).eq.diffx.and.momy(ii).eq.diffy)then
              mdiff(i,j)=ii
           end if
           if(momx(ii).eq.addx.and.momy(ii).eq.addy)then
              madd(i,j)=ii
           end if
        end do
     end do
  end do
  return
end subroutine setup_s
  
subroutine setup_B(ncells,momx,momy,mdiff,madd)
  implicit none
  integer :: ncells,i,j,ii
  integer :: momx(ncells),momy(ncells)
  integer :: mdiff(ncells,ncells),madd(ncells,ncells)
  integer :: kkx,kky,ppx,ppy,diffx,diffy,addx,addy
  
  momx(1)=0
  momy(1)=0
  momx(2)=1
  momy(2)=1
  momx(3)=-1
  momy(3)=1
  momx(4)=1
  momy(4)=-1
  momx(5)=-1
  momy(5)=-1
  momx(6)=2
  momy(6)=0
  momx(7)=0
  momy(7)=2
  momx(8)=2
  momy(8)=2

  do i=1,ncells
     kkx=momx(i)
     kky=momy(i)
     do j=1,ncells
        ppx=momx(j)
        ppy=momy(j)
        diffx=kkx-ppx
        diffy=kky-ppy
        if(diffx.lt.-1)diffx=diffx+4
        if(diffy.lt.-1)diffy=diffy+4
        if(diffx.gt.2)diffx=diffx-4
        if(diffy.gt.2)diffy=diffy-4
        addx=kkx+ppx
        if(addx.gt.2)addx=addx-4
        addy=kky+ppy
        if(addy.ge.2)addy=addy-4
        if(addx.lt.-1)addx=addx+4
        if(addy.lt.-1)addy=addy+4
        do ii=1,ncells
           if(momx(ii).eq.diffx.and.momy(ii).eq.diffy)then
              mdiff(i,j)=ii
           end if
           if(momx(ii).eq.addx.and.momy(ii).eq.addy)then
              madd(i,j)=ii
           end if
        end do
     end do
  end do

  return
end subroutine setup_B

subroutine checks(n,raca,racb,racc,nchecks)
  implicit none
  integer n
  double precision :: raca,racb,racc
  real :: nchecks(16)

  nchecks=0.0
  if(n.eq.2)then
     nchecks(1)=raca+14*racb+7*racc
     nchecks(2)=raca-3*racb+2*racc
     nchecks(3)=raca+4*racb+2*racc
     nchecks(4)=raca+7*racb
     nchecks(5)=raca-8*racb
  end if
  if(n.eq.3)then
     nchecks(1)=3*raca+5*racb+5*racc+dsqrt((193*racb**2+8*racb*racc+4*racc**2))
     nchecks(2)=3*raca+5*racb+5*racc-dsqrt((193*racb**2+8*racb*racc+4*racc**2))
     nchecks(3)=3*raca+9*racb+3*racc
     nchecks(4)=3*raca-11*racb+3*racc
     nchecks(5)=3*raca-6*racb+3*racc
     nchecks(6)=3*raca
     nchecks(7)=3*raca-15*racb
  end if
  if(n.eq.4)then
     nchecks(1)=6*raca+10*racb+10*racc+2*dsqrt(193*racb*racb+8*racb*racc+4*racc*racc)
     nchecks(2)=6*raca+10*racb+10*racc-2*dsqrt(193*racb*racb+8*racb*racc+4*racc*racc)
     nchecks(3)=6*raca+9*racb+7.5*racc+1.5*dsqrt(144*racb*racb+8*racb*racc+racc*racc)
     nchecks(4)=6*raca+9*racb+7.5*racc-1.5*dsqrt(144*racb*racb+8*racb*racc+racc*racc)
     nchecks(5)=6*raca+6*racc
     nchecks(6)=6*raca-5*racb+7.5*racc+0.5*dsqrt(708*racb*racb-12*racb*racc+9*racc*racc)
     nchecks(7)=6*raca-5*racb+7.5*racc-0.5*dsqrt(708*racb*racb-12*racb*racc+9*racc*racc)
     nchecks(8)=6*raca-15*racb+6*racc
     nchecks(9)=6*raca-5*racb+5.5*racc+0.5*dsqrt(912*racb*racb-24*racb*racc+9*racc*racc)
     nchecks(10)=6*raca-5*racb+5.5*racc-0.5*dsqrt(912*racb*racb-24*racb*racc+9*racc*racc)
     nchecks(11)=6*raca-5*racb+4*racc
     nchecks(12)=6*raca-5*racb+5.5*racc+1.5*dsqrt(68*racb*racb+4*racb*racc+racc*racc)
     nchecks(13)=6*raca-5*racb+5.5*racc-1.5*dsqrt(68*racb*racb+4*racb*racc+racc*racc)
     nchecks(14)=6*raca-12*racb+4*racc
     nchecks(15)=6*raca-17*racb+4*racc
     nchecks(16)=6*raca-21*racb
  end if
  if(n.eq.5)then
     nchecks(1)=10*raca-3*racb+8*racc
     nchecks(2)=10*raca+20*racb+10*racc
     nchecks(3)=10*raca-3*racb+11*racc+3*dsqrt(57*racb*racb+2*racb*racc+racc*racc)
     nchecks(4)=10*raca-3*racb+11*racc-3*dsqrt(57*racb*racb+2*racb*racc+racc*racc)
     nchecks(5)=10*raca-4*racb+10*racc
     nchecks(6)=10*raca-25*racb+10*racc
     nchecks(7)=10*raca-9*racb+8*racc
     nchecks(8)=10*raca+3*racb+10*racc
     nchecks(9)=10*raca-13*racb+8*racc
     nchecks(10)=10*raca-22*racb+10*racc
     nchecks(11)=10*raca-24*racb+8*racc
     nchecks(12)=10*raca-28*racb+7*racc
     nchecks(13)=10*raca-18*racb+5*racc
     nchecks(14)=10*raca-13*racb+7*racc
     nchecks(15)=10*raca-25*racb+5*racc
     nchecks(16)=10*raca-35*racb
  end if
  return
end subroutine checks

recursive subroutine quicksort(a, first, last)
  implicit none
  double precision ::  a(*), x, t
  integer first, last
  integer i, j
  
  x = a( (first+last) / 2 )
  i = first
  j = last
  do
     do while (a(i) < x)
        i=i+1
     end do
     do while (x < a(j))
        j=j-1
     end do
     if (i >= j) exit
     t = a(i);  a(i) = a(j);  a(j) = t
     i=i+1
     j=j-1
  end do
  if (first < i-1) call quicksort(a, first, i-1)
  if (j+1 < last)  call quicksort(a, j+1, last)
end subroutine quicksort

subroutine arpackdiag(spr_size,Hsize,nev,ncv,IndexI,IndexJ,sparseH,d,z)
  implicit none

  integer :: nev,ncv,maxncv,spr_size,ii,jj
  integer     ::   Hsize,ldv,ido,lworkl,info,ierr,j,maxitr,mode,ishfts,i
  integer     ::   iparam(11), ipntr(11), IndexI(spr_size), IndexJ(spr_size), ldz
  double precision :: tol, sigma,workl(ncv*(ncv+8)),small 
  double precision :: v(Hsize,ncv),workd(3*Hsize),d(ncv),resid(Hsize),sparseH(spr_size),z(Hsize,ncv)
  
  character        bmat*1, which*2
  logical     ::   rvec,selection(Hsize)

  maxncv=ncv
  bmat = 'I'
  which = 'SA'!Require the smallest algebraic eigen-values.                                                                                                             
  ldv=Hsize
  ldz=Hsize
  small=0.000001d0
  
  !     %--------------------------------------------------%                                                                                                              
  !     | The work array WORKL is used in DSAUPD as        |   
  !     | workspace.  Its dimension LWORKL is set as       |
  !     | illustrated below.  The parameter TOL determines |
  !     | the stopping criterion.  If TOL<=0, machine      |
  !     | precision is used.  The variable IDO is used for |
  !     | reverse communication and is initially set to 0. |
  !     | Setting INFO=0 indicates that a random vector is |
  !     | generated in DSAUPD to start the Arnoldi         |
  !     | iteration.                                       |
  !     %--------------------------------------------------%
  
  tol = 1.0d-14!The tolerance factor for the convergence of Lanczos processes. If set to be zero, the error can be < 1.0d-16.
  info = 0
  ido = 0
  
  !     %---------------------------------------------------%                                                                                                             
  !     | This program uses exact shifts with respect to    |                                                                                                             
  !     | the current Hessenberg matrix (IPARAM(1) = 1).    |                                                                                                             
  !     | IPARAM(3) specifies the maximum number of Arnoldi |                                                                                                             
  !     | iterations allowed.  Mode 1 of DSAUPD is used     |                                                                                                             
  !     | (IPARAM(7) = 1).  All these options may be        |                                                                                                             
  !     | changed by the user. For details, see the         |                 
  !     | documentation in DSAUPD.                          |
  !     %---------------------------------------------------%
  
  ishfts = 1
  maxitr = 300 !Max numver of iterations of the Lanczos processes.
  mode   = 1
  
  iparam = 0
  ipntr = 0
  iparam(1) = ishfts
  iparam(3) = maxitr
  iparam(7) = mode
  lworkl = ncv*(ncv+8)
  
  !*****************************************************************************
  
  !     %-------------------------------------------%
  !     | M A I N   L O O P (Reverse communication) |
  !     %-------------------------------------------%
  
10 continue
  
  !        %---------------------------------------------%                                                                                                                
  !        | Repeatedly call the routine DSAUPD and take |                                                                                                                
  !        | actions indicated by parameter IDO until    |                                                                                                                
  !        | either convergence is indicated or maxitr   |                                                                                                                
  !        | has been exceeded.                          |                                                                                                                
  !        %---------------------------------------------%                                                                                                                
  
  call dsaupd(ido,bmat,Hsize,which,nev,tol,resid,ncv,v,ldv,iparam,ipntr,workd,workl,lworkl,info)
  
  if (ido .eq. -1 .or. ido .eq. 1) then
     
     !           %--------------------------------------%
     !           | Perform matrix vector multiplication |
     !           |              y <--- OP*x             |
     !           | The user should supply his/her own   |
     !           | matrix vector multiplication routine |
     !           | here that takes workd(ipntr(1)) as   |
     !           | the input, and return the result to  |
     !           | workd(ipntr(2)).                     |
     !           %--------------------------------------%
     
     !*********************************************************************************
     
     workd(ipntr(2):ipntr(2)+Hsize-1)=0
     do ii=1, spr_size
        workd(ipntr(2)+IndexI(ii)-1)=workd(ipntr(2)+IndexI(ii)-1)+workd(ipntr(1)+IndexJ(ii)-1)*sparseH(ii)
     end do
     
     !*********************************************************************************
     
     !           %-----------------------------------------%                                                                                                            
     !           | L O O P   B A C K to call DSAUPD again. |                                                                                                            
     !           %-----------------------------------------%                                                                                                            
     
     go to 10
     
  end if

  !     %----------------------------------------%                                                                                                                        
  !     | Either we have convergence or there is |                                                                                                                        
  !     | an error.                              |                                                                                                                        
  !     %----------------------------------------%                                                                                                                        
  
  
  if ( info .lt. 0 ) then
     
     !        %--------------------------%                                                                                                                                   
     !        | Error message. Check the |                                                                                                                                   
     !        | documentation in DSAUPD. |                                                                                                                                   
     !        %--------------------------%                                                                                                                                   
     
     print *, ' '
     print *, ' Error with _saupd, info = ', info
     print *, ' Check documentation in _saupd '
     print *, ' '
     stop
     
  else
     
     !        %-------------------------------------------%                                                                                                                  
     !        | No fatal errors occurred.                 |                                                                                                                  
     !        | Post-Process using DSEUPD.                |                                                                                                                  
     !        |                                           |                                                                                                                  
     !        | Computed eigenvalues may be extracted.    |                                                                                                                  
     !        |                                           |                                                                                                                  
     !        | Eigenvectors may also be computed now if  |                                                                                                                  
     !        | desired.  (indicated by rvec = .true.)    |                                                                                                                  
     !        %-------------------------------------------%                                                                                                                  
     
     rvec = .true.
     
     call dseupd(rvec,'All',selection,d,z,ldz,sigma,bmat,Hsize,which,nev,tol,resid,ncv,v,ldv,iparam,ipntr,workd,workl,lworkl,ierr)
     
     !        %----------------------------------------------%
     !        | Eigenvalues are returned in the first column |
     !        | of the two dimensional array D and the       |
     !        | corresponding eigenvectors are returned in   |
     !        | the first NEV columns of the two dimensional |
     !        | array V if requested.  Otherwise, an         |
     !        | orthogonal basis for the invariant subspace  |
     !        | corresponding to the eigenvalues in D is     |
     !        | returned in V.                               |
     !        %----------------------------------------------%
     
     do ii=1, nev
        !      write(6,*)ii, d(ii)                                                                                                                                              
        do jj=1, Hsize
           z(jj,ii)=v(jj,ii)
           !         if(abs(v(jj,ii)).gt.small)write(6,*)jj,v(jj,ii)                                                                                                               
        end do
        !      write(6,*)';;;;;;;;'                                                                                                                                             
     end do
     !   close(17)                                                                                                                                                           
     
     if ( ierr .ne. 0) then
        !            %------------------------------------%                                                                                                                     
        !            | Error condition:                   |                                                                                                                     
        !            | Check the documentation of DSEUPD. |                                                                                                                     
        !            %------------------------------------%                                                                                                                     
        print *, ' '
        print *, ' Error with _seupd, info = ', ierr
        print *, ' Check the documentation of _seupd. '
        print *, ' '
        stop
     else
20      continue
     !            %-------------------------------%                                                                                                                          
     !            | Display computed residuals    |                                                                                                                          
     !            %-------------------------------%
     end if
  
  !        %------------------------------------------%                                                                                                                   
  !        | Print additional convergence information |                                                                                                                   
  !        %------------------------------------------%                                                                                                                   
  
     if ( info .eq. 1) then
        print *, ' '
        print *, ' Maximum number of iterations reached.'
        print *, ' '
     else if ( info .eq. 3) then
        print *, ' '
        print *, ' No shifts could be applied during implicit',&
             ' Arnoldi update, try increasing NCV.'
        print *, ' '
     end if
     
     print *, ' '
     print *, ' _SDRV1 '
     print *, ' ====== '
     print *, ' '
     print *, ' Size of the matrix is ', Hsize
     print *, ' The number of Ritz values requested is ', nev
     print *, ' The number of Arnoldi vectors generated (NCV) is ', ncv
     print *, ' What portion of the spectrum: ', which
     print *, ' The number of Implicit Arnoldi update iterations taken is ', iparam(3)
     print *, ' The number of OP*x is ', iparam(9)
     print *, ' The convergence criterion is ', tol
     print *, ' '
     
  end if

  !     %---------------------------%                                                                                                                                     
  !     | Done with program arpack. |                                                                                                                                     
  !     %---------------------------%                                                                                                                                     

  return
end subroutine arpackdiag

subroutine EX_ARPACK(N,NEV,NCV,MATRICE,D,Z)
  IMPLICIT NONE
  INTEGER :: N
  INTEGER :: NEV
  INTEGER :: NCV
  INTEGER :: IDO, LDV, LDZ
  INTEGER :: IPARAM(11), LWORKL, INFO
  INTEGER :: IPNTR(11), PURSUE
  LOGICAL :: SELECTION(N)
  REAL(KIND=8) :: TOL, RESID(N), V(N,N), MATRICE(N,N)
  REAL(KIND=8) :: WORKD(3*N), WORKL(NCV*(NCV+8))
  REAL(KIND=8) :: Z(N,NEV), D(NEV), SIGMA
!---------------------------------
 
  LDV=N
  LDZ=N
  IDO=0
  PURSUE=1
  INFO=0
  IPARAM=0
  IPNTR=0
  TOL=0D0
!---------------------------------
!ENTER THE MATRIX VALUES
!---------------------------------
!  MATRICE=0D0
!  MATRICE(1,1)=-dfloat(1)/dfloat(6)+7
!  MATRICE(5,5)=-dfloat(1)/dfloat(6)
!  MATRICE(1,5)=-dfloat(5)/dfloat(6)
!  MATRICE(5,1)=MATRICE(1,5)
!  MATRICE(2,2)=dfloat(4)/dfloat(6)
!  MATRICE(4,4)=MATRICE(2,2)+3
!  MATRICE(3,3)=-dfloat(1)
!----------------------------------
  
!----------------------------------
! CALCULATE THE EIGENVALUES AND EIGENVECTORS 
!----------------------------------
  IPARAM(1)=1  
  IPARAM(3)=25
  IPARAM(7)=1
  LWORKL=NCV*(NCV+8)
  DO WHILE (PURSUE==1) 
     CALL DSAUPD(IDO,'I',N,'SA',NEV,TOL,RESID,NCV,V, &
          & LDV, IPARAM, IPNTR, WORKD, WORKL, LWORKL, INFO)
     IF (IDO==1 .OR. IDO==-1) THEN
        CALL MAT_VEC(N,MATRICE,WORKD(IPNTR(1)),WORKD(IPNTR(2)))
     ELSE IF (IDO==99) THEN
        PURSUE=0
     END IF
  END DO
  !------------------------------------
  
 
  !------------------------------------
  !POST PROCESSING WITH DSEUPD
  !------------------------------------
  CALL DSEUPD(.TRUE.,'A',SELECTION,D,Z,LDZ,SIGMA, &
       & 'I', N, 'SM', NEV, TOL, RESID, NCV, V, &
       & LDV, IPARAM, IPNTR, WORKD, WORKL, LWORKL, INFO)
  !------------------------------------
  !****************************************************************
  
  
  !****************************************************************
  !**************        DISPLAY THE RESULTS            ************
  !****************************************************************
!  PRINT*, '#1 : '
!  PRINT*, '    EIGENVALUE  -> ', D(1)
!  PRINT*, '    EIGENVECTOR -> ', Z(:,1)
!  PRINT*, '#2 : '
!  PRINT*, '    EIGENVALUE -> ', D(2)
!  PRINT*, '    EIGENVECTOR -> ', Z(:,2)
  !****************************************************************
  
  
  !****************************************************************
  !******************      END OF SUBROUTINE      *****************
  !****************************************************************
  return
end subroutine EX_ARPACK

SUBROUTINE MAT_VEC(N,MATRICE,X,Y)
  !--------------------------------------------
  
  !--------------------------------------------
  ! VARIABLE DECLARATION 
  !--------------------------------------------
  IMPLICIT NONE
  INTEGER :: I, J, N
  REAL(KIND=8) :: X(N), Y(N), MATRICE(N,N)
  !--------------------------------------------
 
  !--------------------------------------------
  ! MATRIX-VECTOR MULTIPLICATION 
  !--------------------------------------------
  Y=0D0
  DO I=1,N
     DO J=1,N
        Y(I)=Y(I)+MATRICE(I,J)*X(J)
     END DO
  END DO
  !--------------------------------------------
  
  !--------------------------------------------
  ! END OF THE SUBROUTINE 
  !--------------------------------------------
  return
END SUBROUTINE MAT_VEC

