MODULE tomog_precision
  USE number_types
  INTEGER,PARAMETER::prec = dp;
END MODULE tomog_precision

!====================================================
! This sets all the data structure and methods for 
! the tomography model 
!====================================================
MODULE tomog_model

  USE tomog_precision

!!!!!!!!!!!!!!!!!!
  IMPLICIT NONE
!!!!!!!!!!!!!!!!!!

!=============================================================
! Data Structure for holding the tomogram : TmgModel
!
!
! For Elliptical Tomography:
! Np - # of parameters in each cell
! grid size is Nx X Nz
! x0,z0 - upper left corner coordinates
! dx,dz - cell dimensions
! val(1:Np,1:Nz,1:Nx) - holds model parameters values
!
  TYPE TmgModel
     INTEGER :: Np,Nx,Nz
     REAL(PREC) :: x0,z0,dx,dz
     REAL(PREC),DIMENSION(:,:,:),POINTER :: val=>null()
  END TYPE TmgModel
!==============================================================

CONTAINS

  FUNCTION Set_TmgModel(Np,Nx,Nz,x0,z0,dx,dz,Model) RESULT(ierr)
    IMPLICIT NONE
    INTEGER,INTENT(in)           :: Np,Nx,Nz
    REAL(PREC),INTENT(in)        :: x0,z0,dx,dz
    TYPE(TmgModel),INTENT(INOUT) :: Model
    INTEGER::ierr

    Model%Np = Np;
    Model%Nx = Nx;
    Model%Nz = Nz;
    Model%x0 = x0;Model%z0 = z0;
    Model%dx = dx;Model%dz = dz;


    IF(ASSOCIATED(Model%val)) THEN
       DEALLOCATE(Model%val,STAT= ierr)
       IF(ierr) RETURN
       NULLIFY(Model%val)
    END IF

    ALLOCATE(Model%val(3,Nx,Nz),STAT= ierr)
    IF(ierr /= 0) STOP 'Allocation Failed for Slowness Model'

  END FUNCTION Set_TmgModel

  FUNCTION deallocate_TmgModel(Model) RESULT(ierr)
    IMPLICIT NONE
    TYPE(TmgModel),INTENT(INOUT) :: Model
    INTEGER::ierr

    IF(ASSOCIATED(Model%val)) THEN
       DEALLOCATE(Model%val,STAT= ierr)
       IF(ierr) RETURN
       NULLIFY(Model%val)
    END IF

  END FUNCTION deallocate_TmgModel

END MODULE tomog_model


MODULE tomog_geometry
!====================================================
! This sets all the data structure and methods for
! for tomography acquisition  geometry
!====================================================

  USE tomog_precision

!!!!!!!!!!!!!!!!!!!
  IMPLICIT NONE
!!!!!!!!!!!!!!!!!!!

!=============================================================
! Data Structure for holding acquisition geometry: TmgGeom
!
!
! For Ellitical Tomography:
! Nrays         - # of rays 
! xrec(1:Nrays) - receivers coordinates
! zrec(1:Nrays) - receivers coordinates
! xsou(1:Nrays) - sources coordinates
! zsou(1:Nrays) - sources coordinates
!
!==============================================================
  TYPE TmgGeom
     INTEGER :: Nrays

     REAL(PREC),DIMENSION(:),POINTER::xrec=>null()
     REAL(PREC),DIMENSION(:),POINTER::zrec=>null()

     REAL(PREC),DIMENSION(:),POINTER::xsou=>null()
     REAL(PREC),DIMENSION(:),POINTER::zsou=>null()

  END TYPE TmgGeom

  INTERFACE Allocate_TmgGeom
     MODULE PROCEDURE Set_TmgGeom
  END INTERFACE

CONTAINS

  FUNCTION Set_TmgGeom(Ndata,Geometry) RESULT(ierr)
    IMPLICIT NONE
    INTEGER,INTENT(in)           :: Ndata
    TYPE(TmgGEOM),INTENT(INOUT)  :: Geometry
    INTEGER::ierr

    Geometry%Nrays = Ndata;

    IF(ASSOCIATED(Geometry%xsou)) THEN
       DEALLOCATE(Geometry%xsou,STAT= ierr)
       IF(ierr) RETURN
       NULLIFY(Geometry%xsou)
    END IF


    IF(ASSOCIATED(Geometry%zsou)) THEN
       DEALLOCATE(Geometry%zsou,STAT= ierr)
       IF(ierr) RETURN
       NULLIFY(Geometry%zsou)
    END IF

    IF(ASSOCIATED(Geometry%xrec)) THEN
       DEALLOCATE(Geometry%xrec,STAT= ierr)
       IF(ierr) RETURN
       NULLIFY(Geometry%xrec)
    END IF

    IF(ASSOCIATED(Geometry%zrec)) THEN
       DEALLOCATE(Geometry%zrec,STAT= ierr)
       IF(ierr) RETURN
       NULLIFY(Geometry%zrec)
    END IF

    ALLOCATE(Geometry%xsou(Ndata),STAT= ierr)
    IF(ierr) RETURN
    ALLOCATE(Geometry%zsou(Ndata),STAT= ierr)
    IF(ierr) RETURN
    ALLOCATE(Geometry%xrec(Ndata),STAT= ierr)
    IF(ierr) RETURN
    ALLOCATE(Geometry%zrec(Ndata),STAT= ierr)
    IF(ierr) RETURN

  END FUNCTION Set_TmgGeom


  FUNCTION deallocate_TmgGeom(Geometry) RESULT(ierr)
    IMPLICIT NONE
    TYPE(TmgGEOM),INTENT(INOUT) :: Geometry
    INTEGER::ierr

    IF(ASSOCIATED(Geometry%xsou)) THEN
       DEALLOCATE(Geometry%xsou,STAT= ierr)
       IF(ierr == 0) RETURN
       NULLIFY(Geometry%xsou)
    END IF

    IF(ASSOCIATED(Geometry%zsou)) THEN
       DEALLOCATE(Geometry%zsou,STAT= ierr)
       IF(ierr == 0) RETURN
       NULLIFY(Geometry%zsou)
    END IF

    IF(ASSOCIATED(Geometry%xrec)) THEN
       DEALLOCATE(Geometry%xrec,STAT= ierr)
       IF(ierr == 0) RETURN
       NULLIFY(Geometry%xrec)
    END IF

    IF(ASSOCIATED(Geometry%zrec)) THEN
       DEALLOCATE(Geometry%zrec,STAT= ierr)
       IF(ierr == 0) RETURN
       NULLIFY(Geometry%zrec)
    END IF

  END FUNCTION deallocate_TmgGeom


END MODULE tomog_geometry


MODULE tomog_Modeling
!==========================================================!$
! This should assemble all the tomography
! together, hopefully in the future and with your help 
! it will get better and better
!==========================================================!$

  USE tomog_model
  USE tomog_geometry

  USE sparse
  USE LINSYS_CG
  USE LINSYS_LSQR
  USE LINSYS_SIRT
  USE LINSYS_ART

!  USE raysmooth_2d


!!!!!!!!!!!!!!!!!!!
  IMPLICIT NONE
!!!!!!!!!!!!!!!!!!!

  ! Tomogram :

  TYPE(TmgModel),PRIVATE,POINTER    :: TSlow => Null()         ! Model
  TYPE(TmgGeom), PRIVATE,POINTER    :: Tgeom => Null()         ! Geometry

  REAL(PREC),PRIVATE, POINTER :: Tobs(:) => Null()     ! Traveltime data

  ! Linear system solver:

  REAL(PREC),PRIVATE, POINTER :: Trhs(:) => Null()  
  REAL(PREC),PRIVATE, POINTER :: Tsol(:) => Null()  
  TYPE(sparseMatrix) , PRIVATE :: Tmatrix 

  REAL(PREC),PRIVATE :: w_data,compliance_min,compliance_max

  ! Useful variables to checking if geometry fits the model grid
  REAL(PREC) :: xsmin,xsmax,zsmin,zsmax
  REAL(PREC) :: xrmin,xrmax,zrmin,zrmax

  ! Useful variables to measure data fitting 
  REAL(PREC) :: tmg_data_mean,tmg_data_std,tmg_rsd_mean,tmg_rsd_std

  !Raytracing Parameters
  INTEGER,PARAMETER :: ODE_DIM      = 5;    ! order of the ODE system for raytracing
  INTEGER,PARAMETER :: FAN_SIZE     = 360;  ! number of rays in the rays fan for each source 
  INTEGER,PARAMETER :: MAXRAY_STEPS = 1000; ! maximum number of steps along a ray for Frechet derivative computation 
  REAL(dp) :: rpath(ODE_DIM,MAXRAY_STEPS)   ! raypath coordinates


CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SETUP routines
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE Tmg_SetModel(Np,Nx,Nz,x0,z0,dx,dz)
    INTEGER    ::Np,Nx,Nz,ierr
    REAL(PREC)::dx,dz,x0,z0
    REAL(PREC) :: xs(2),xr(2),dist,vel_avg
    INTEGER :: i

    ALLOCATE(Tslow);
    ierr = Set_TmgModel(Np,Nx,Nz,x0,z0,dx,dz,Tslow);

    vel_avg = 0.d0;

    DO i=1,Tgeom%nrays

       xs = (/Tgeom%xsou(i),Tgeom%zsou(i)/);
       xr = (/Tgeom%xrec(i),Tgeom%zrec(i)/);

       dist = SQRT(DOT_PRODUCT(xr-xs,xr-xs));

       vel_avg = vel_avg + dist/Tobs(i)
    END DO

    vel_avg = vel_avg/float(Tgeom%nrays);

    Tslow%val(:,:,:) = 0.d0;
    DO i=1,2
       Tslow%val(i,:,:) = vel_avg**2;
    END DO


  END SUBROUTINE Tmg_SetModel


  SUBROUTINE TmgSetup(fmodel,fdata)
    CHARACTER*(*),INTENT(in) :: fmodel,fdata
    INTEGER::ierr
    INTEGER               ::Nrow,Ncol,Ndata
    INTEGER,PARAMETER::fid_m = 77;
    INTEGER,PARAMETER::fid_d = 67;
    INTEGER :: i

    OPEN(unit=fid_d,file=fdata,status='old')
    READ(fid_d,*,err=1000) ndata

    ierr = Set_TmgGeom(Ndata,TGeom);

    ! testar status
    ALLOCATE(Tobs(ndata),STAT=ierr)
    IF(ierr) THEN
       STOP 'No memmory available for Tobs'
    END IF

    DO i = 1,ndata
       READ(fid_d,*,err=1000) TGeom%xsou(i),TGeom%zsou(i), &
            TGeom%xrec(i),TGeom%zrec(i), &
            Tobs(i)
    END DO
1000 CLOSE(fid_d)

    xsmin = MINVAL(TGeom%xsou,dim=1);
    xsmax = MAXVAL(TGeom%xsou,dim=1);

    zsmin = MINVAL(TGeom%zsou,dim=1);
    zsmax = MAXVAL(TGeom%zsou,dim=1);

    xrmin = MINVAL(TGeom%xrec,dim=1);
    xrmax = MAXVAL(TGeom%xrec,dim=1);

    zrmin = MINVAL(TGeom%zrec,dim=1);
    zrmax = MAXVAL(TGeom%zrec,dim=1);

  END SUBROUTINE TmgSetup

  
  SUBROUTINE Tmg_setlinsys()

    INTEGER :: Ncol,Nrow,NzeroMax,ierr

    Ncol     = Tslow%np*Tslow%nx*Tslow%nz;
    Nrow     = Tgeom%nrays + ncol + ncol + ncol + ncol;
    NzeroMax = Tgeom%nrays*(12*Tslow%np*MAX(Tslow%nx,Tslow%nz)) + ncol + 8*ncol;

    ierr = SparseMatrix_Allocate(Nrow,Ncol,NzeroMax, Tmatrix)
    IF(ierr) THEN
       STOP 'No memmory available for TMatrix'
    END IF

    IF(.NOT.ASSOCIATED(Tsol)) THEN
       ALLOCATE(Tsol(Ncol),STAT = ierr);
       IF(ierr) THEN
          STOP 'No memmory available for Tsol'
       END IF
    ELSE
       DEALLOCATE(Tsol,STAT = ierr);
       Tsol => null();
       ALLOCATE(Tsol(Ncol),STAT = ierr);
       IF(ierr) THEN
          STOP 'No memory available for Tsol'
       END IF
    END IF

    IF(.NOT.ASSOCIATED(Trhs)) THEN
       ALLOCATE(Trhs(Nrow),STAT = ierr);
       IF(ierr) THEN
          STOP 'No memory available for Trhs'
       END IF
    ELSE
       DEALLOCATE(Trhs,STAT = ierr);
       Trhs => null();
       ALLOCATE(Trhs(Nrow),STAT = ierr);
       IF(ierr) THEN
          STOP 'No memory available for Trhs'
       END IF
    END IF

  END SUBROUTINE Tmg_setlinsys
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! I/O Routines
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Tmg_read_Data(fdata)
    CHARACTER*(*),INTENT(in) :: fdata
    INTEGER::i,ierr,ndata
    INTEGER,PARAMETER::fid_d = 67;
    REAL::xs,zs,xr,zr,tvt

    OPEN(unit=fid_d,file=fdata,status='old')
    ndata=0
    DO
       READ(fid_d,*,err=1000,END=3000) xs,zs,xr,zr,tvt
       ndata = ndata + 1
    END DO
3000 CLOSE(fid_d)

    ALLOCATE(TGeom,STAT=ierr)
    IF(ierr) THEN
       STOP 'No memmory available for TGeom'
    END IF

    ierr = Set_TmgGeom(Ndata,TGeom);
    IF(ierr) THEN
       STOP 'Allocation failed for TmgGeom'
    END IF

    ! testar status
    ALLOCATE(Tobs(ndata),STAT=ierr)
    IF(ierr) THEN
       STOP 'No memmory available for Tobs'
    END IF

    OPEN(unit=fid_d,file=fdata,status='old')
    DO i = 1,ndata
       READ(fid_d,*,err=1000) TGeom%xsou(i),TGeom%zsou(i), &
            TGeom%xrec(i),TGeom%zrec(i), &
            Tobs(i)
    END DO
1000 CLOSE(fid_d)

    xsmin = MINVAL(TGeom%xsou,dim=1);
    xsmax = MAXVAL(TGeom%xsou,dim=1);

    zsmin = MINVAL(TGeom%zsou,dim=1);
    zsmax = MAXVAL(TGeom%zsou,dim=1);

    xrmin = MINVAL(TGeom%xrec,dim=1);
    xrmax = MAXVAL(TGeom%xrec,dim=1);

    zrmin = MINVAL(TGeom%zrec,dim=1);
    zrmax = MAXVAL(TGeom%zrec,dim=1);

! Data statistics:

! mean
    tmg_data_mean =  SUM(Tobs(1:ndata))/float(ndata);
! standard deviation
    tmg_data_std  =  SQRT(DOT_PRODUCT((Tobs(1:ndata)-tmg_data_mean), &
         (Tobs(1:ndata)-tmg_data_mean)  &
         )/float(ndata-1)               &
         );

  END SUBROUTINE Tmg_read_Data

  SUBROUTINE Tmg_read_model(fmodel)
    CHARACTER*(*),INTENT(in) :: fmodel
    INTEGER::ierr
    INTEGER,PARAMETER::fid_m = 77;
    INTEGER :: np,nx,nz,ip,ix,iz
    REAL(PREC) :: x0,z0,dx,dz
 
    REAL :: value
    CHARACTER(LEN=200) :: fbin

     OPEN(unit=fid_m,file=fmodel,status='old')
     READ(fid_m,*,err=1000) np,nx,nz,x0,z0,dx,dz
     READ(fid_m,'(a)') fbin
1000 CLOSE(fid_m)

    PRINT *,np,nx,nz,x0,z0,dx,dz

    ALLOCATE(Tslow);
    ierr = Set_TmgModel(Np,Nx,Nz,x0,z0,dx,dz,Tslow)
 
    OPEN(unit=fid_m,file=fbin,form='unformatted',recordtype='stream')
    DO ip=1,Np
       DO iz=1,Nz
          DO ix=1,Nx
             READ(fid_m) value 
             Tslow%val(ip,ix,iz) = value;
          END DO
       END DO
    END DO


    CALL tomog_output('init_tomog')

  END SUBROUTINE Tmg_read_model


  SUBROUTINE tomog_output(fname)
    CHARACTER*(*) :: fname
    CHARACTER(len=200) :: fout
    REAL(PREC)::media,desv
    INTEGER :: Np,Nx,Nz
    INTEGER :: ip,ix,iz

    Np=Tslow%Np;Nx= Tslow%Nx;Nz=Tslow%Nz;

    fout=TRIM(ADJUSTL(fname))//'.hdr'

    OPEN(36,File = fout,status= 'unknown');
    WRITE(36,'(i3,x,i3,x,i3,x,4(f7.2))') Tslow%Np,Tslow%Nx,Tslow%Nz,Tslow%x0,Tslow%z0,Tslow%dx,Tslow%dz
    WRITE(36,'(a)') TRIM(ADJUSTL(fname))//'.bin'
    CLOSE(36)

    fout=TRIM(ADJUSTL(fname))//'.bin'
    OPEN(39,file=fout,form='unformatted',recordtype='stream',status='unknown')

    DO ip=1,Np
       DO iz=1,Nz
          WRITE(39) (REAL(Tslow%val(ip,ix,iz)), ix = 1,Nx)
       END DO
    END DO

    CLOSE(39)

  END SUBROUTINE tomog_output

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  RAYTRACING AND NON-LINEAR TOMOGRAPHY ITERATIONS  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Tmg_Display_Limits()

    WRITE(*,'("Source   : "/,"xmin = ",f10.5,"  xmax = ",f10.5,/,   &
         "zmin = ",f10.5,"  zmax = ",f10.5,//,  &
         "Receiver : "/,"xmin = ",f10.5,"  xmax = ",f10.5,/,   &
         "zmin = ",f10.5,"  zmax = ",f10.5,/)') &
         xsmin,xsmax,zsmin,zsmax,xrmin,xrmax,zrmin,zrmax

  END SUBROUTINE Tmg_Display_Limits


  SUBROUTINE tomog_solve(method,lambda_1,lambda_2,lambda_3)

    REAL(PREC),PARAMETER     :: TOL=1.0e-3
    INTEGER,INTENT(in)       :: method
    REAL(PREC),INTENT(INOUT) :: lambda_1,lambda_2,lambda_3

    REAL(PREC),ALLOCATABLE :: cnorm(:)
    REAL(PREC)             :: res, damping_homog, damping_iso,damping_horiz, cnorm_max
    REAL(PREC)             :: xs(2),xr(2),dist,vel_avg,tmg_rsd_std_0,alpha,dalpha,tsol_norm
    INTEGER                :: niter,nitmax,nrays,ncol,nzero,np,nx,nz,n_iter
    INTEGER                :: ierr,istop,iray,irow,index,ip,ix,iz,icol,irtrace
    INTEGER,PARAMETER      :: n_iter_max = 50
    INTEGER,PARAMETER      :: imodel = 3
    CHARACTER(LEN=32) :: fout
    CHARACTER(LEN=3)  :: fidx

    CALL Tmg_setLinsys(); ! allocate memory for linear system

    ALLOCATE(cnorm(Tmatrix%Ncol),STAT = ierr)
    IF(ierr) THEN
       STOP 'No memmory available for cnorm'
    END IF

    ! Always make sure allocation status has been checked

    ! Start from a homogeneous isotropic model:


    ! Possible Initial Homogeneous Models
    !
!!$    SELECT CASE(imodel)
!!$
!!$    CASE(1) ! Minimum slowness :
!!$
!!$       vel_avg=1.0E32
!!$       DO iray=1,Tgeom%nrays
!!$
!!$          xs = (/Tgeom%xsou(iray),Tgeom%zsou(iray)/);
!!$          xr = (/Tgeom%xrec(iray),Tgeom%zrec(iray)/);
!!$
!!$          dist = SQRT(DOT_PRODUCT(xr-xs,xr-xs));
!!$
!!$          vel_avg = MIN(vel_avg,dist/Tobs(iray))
!!$       END DO
!!$
!!$    CASE(2)
!!$       !
!!$       !  use  average slowness to weight the data:
!!$       !
!!$       vel_avg=0.d0;
!!$       DO iray=1,Tgeom%nrays
!!$
!!$          xs = (/Tgeom%xsou(iray),Tgeom%zsou(iray)/);
!!$          xr = (/Tgeom%xrec(iray),Tgeom%zrec(iray)/);
!!$
!!$          dist = SQRT(DOT_PRODUCT(xr-xs,xr-xs));
!!$
!!$          vel_avg = vel_avg + (dist/Tobs(iray))**2
!!$       END DO
!!$
!!$       vel_avg = SQRT(vel_avg/float(Tgeom%nrays));
!!$
!!$    CASE(3)
!!$       !
!!$       !  use  average slowness to weight the data:
!!$       !
!!$       vel_avg=0.d0;
!!$       DO iray=1,Tgeom%nrays
!!$
!!$          xs = (/Tgeom%xsou(iray),Tgeom%zsou(iray)/);
!!$          xr = (/Tgeom%xrec(iray),Tgeom%zrec(iray)/);
!!$
!!$          dist = SQRT(DOT_PRODUCT(xr-xs,xr-xs));
!!$
!!$          vel_avg = vel_avg + dist/Tobs(iray)
!!$       END DO
!!$
!!$       vel_avg = vel_avg/float(Tgeom%nrays);
!!$
!!$    CASE DEFAULT
!!$
!!$       !
!!$       !  use  average slowness to weight the data:
!!$       !
!!$       vel_avg=0.d0;
!!$       DO iray=1,Tgeom%nrays
!!$
!!$          xs = (/Tgeom%xsou(iray),Tgeom%zsou(iray)/);
!!$          xr = (/Tgeom%xrec(iray),Tgeom%zrec(iray)/);
!!$
!!$          dist = SQRT(DOT_PRODUCT(xr-xs,xr-xs));
!!$
!!$          vel_avg = vel_avg + dist/Tobs(iray)
!!$       END DO
!!$
!!$       vel_avg = vel_avg/float(Tgeom%nrays);
!!$
!!$
!!$    END SELECT
!!$
!!$    w_data = vel_avg;
!!$
!!$    Tslow%val(:,:,:) = 0.d0;
!!$    DO ip=1,2
!!$       Tslow%val(ip,:,:) = vel_avg**2;
!!$    END DO

    CALL set_aijmodel(tslow%np,tslow%nx,tslow%nz,tslow%x0,tslow%z0,tslow%dx,tslow%dz,Tslow%val)

    CALL ray_fan_init(FAN_SIZE,1) ! Allocate memory for the ray fan

    ! Lets give a try to continuation method for regularization (might be turned of in the future)
    !    Cont_loop:DO 

    !       IF(lambda_1 < 0.01.AND.method /= 1) EXIT Cont_loop

    tmg_rsd_std_0 = tmg_data_std;
    n_iter    = 0

    ! Step length control: the ideia is not to allow
    ! larger steps in the slowness model at the beginning of
    ! the non-linear iterations 
    dalpha    = 0.10_PREC
    alpha     = 0.0_PREC

    NL_tomog_loop: DO    ! Main Loop for nonlinear tomography

       n_iter    = n_iter + 1

       alpha = MIN(alpha + dalpha, 1.0_prec)

       CALL raytracing();

       CALL tomog_residuals()

       WRITE(fidx,'(I3)') 100+n_iter
       fout = "partial_"//fidx(2:3)
       CALL tomog_output(fout)

       ! Stop Criteria:

       IF(N_iter > N_iter_max) THEN
          istop = 3
          EXIT NL_tomog_loop
       END IF

       ! if residuals do not change much between iterations 
       ! or if number of iterations reach maximum allowed value 
       IF( ABS(tmg_rsd_std-tmg_rsd_std_0) < TOL*tmg_rsd_std_0 ) THEN
          istop = 1;
          EXIT NL_tomog_loop
       END IF

       ! if data is fitted under the specified accuracy  
       IF( tmg_rsd_std < TOL * tmg_data_std ) THEN
          istop = 2;
          EXIT NL_tomog_loop
       END IF

       tmg_rsd_std_0 =  tmg_rsd_std

       !=============================================================
       !
       ! Solution od the linear system :
       !
       !=============================================================


       ! If regularization is required compute column norm  
       ! in order to normalize damping values

       cnorm = 0.0_PREC;
       nrays = Tmatrix%nrow;
       ncol  = Tmatrix%Ncol;
       nzero = TMatrix%nzero;

       np = Tslow%np;
       nx = Tslow%nx;
       nz = Tslow%nz;

       DO index = 1,Nzero
          cnorm(TMatrix%icol(index)) = cnorm(TMatrix%icol(index)) + &
               (TMatrix%elem(index))**2
       END DO

       cnorm = SQRT(cnorm);
       !======================================================================
       !
       ! Regularization :
       !
       ! 1-damping * Identity Matrix  * x = 0
       !   tries to determine the solution closest to initial model that fits the data
       !
       ! 2-make a damping in favor of isotropy this is necessary where ray coverage is
       !   not enough to recover anisotropy
       !
       !
       !        |     |     |   |
       !        |  A  | x   | t |
       !        |     |   = |   |
       !        |  I  |     | 0 |
       !        |     |     |   |
       !        | Dx  |     | 0 |
       !        |     |     |   |
       !
       !======================================================================
       cnorm_max = MAXVAL(cnorm,dim=1);
       irow    = Tmatrix%Nrow  
       PRINT *,'  Nrows = ',tmatrix%nrow


!!$       IF(lambda_1 > 0) THEN ! this makes the model close to the initial model
!!$          damping_homog = lambda_1 * cnorm_max;
!!$          DO index = 1,Tmatrix%Ncol
!!$             iray  = iray  + 1;
!!$             nzero = nzero + 1;
!!$             TRHS(iray)          = 0.0_PREC;
!!$             Tmatrix%irow(iray)  = nzero;
!!$
!!$             Tmatrix%icol(nzero) = index;
!!$             Tmatrix%elem(nzero) = damping_homog; 
!!$          END DO
!!$       END IF

       IF(lambda_1 > 0) THEN ! this makes the model close to the initial model
          damping_homog = lambda_1 * cnorm_max;

          DO iz=1,Tslow%Nz
             DO ix=1,Tslow%Nx
                DO ip=1,TSlow%Np
                   irow  = irow  + 1;
                   TRHS(irow)          = -damping_homog * Tslow%val(ip,ix,iz);
                   Tmatrix%irow(irow)  = nzero+1;

                   nzero = nzero + 1;
                   Tmatrix%icol(nzero) = ip + ((ix-1)+(iz-1)*Nx)*Np;
                   Tmatrix%elem(nzero) = damping_homog; 
                END DO
             END DO
          END DO

       END IF


       IF(lambda_2 > 0) THEN !this  makes the model isotropic 
          damping_iso = lambda_2 * cnorm_max;

          DO iz = 1,Nz
             DO ix = 1,Nx         

                irow                  = irow  + 1;
                TRHS(irow)            = damping_iso*(Tslow%val(2,ix,iz)-Tslow%val(1,ix,iz));
                nzero                 = nzero + 1;
                Tmatrix%irow(irow)    = nzero;

                Tmatrix%icol(nzero) = ((ix-1) + (iz-1)*Nx)*Np + 1;
                Tmatrix%elem(nzero) = damping_iso; 

                nzero                 = nzero + 1;
                Tmatrix%icol(nzero) = ((ix-1) + (iz-1)*Nx)*Np + 2;
                Tmatrix%elem(nzero) =-damping_iso; 

             END DO
          END DO
       END IF

       IF(lambda_3 > 0) THEN
          damping_horiz = lambda_3 * cnorm_max;

          DO iz = 1,Nz
             DO ix = 1,Nx-1		
                DO ip=1,Np
                   irow                  = irow  + 1;
                   nzero                 = nzero + 1;
                   TRHS(irow)            = damping_horiz*(Tslow%val(ip,ix+1,iz)-Tslow%val(ip,ix,iz)); !0.0_PREC;
                   Tmatrix%irow(irow)    = nzero;

                   Tmatrix%icol(nzero) = ((ix-1) + (iz-1)*Nx)*Np + ip;
                   Tmatrix%elem(nzero) = damping_horiz; 

                   nzero = nzero + 1;
                   Tmatrix%icol(nzero) = ((ix-1) + (iz-1)*Nx)*Np + ip + Np;
                   Tmatrix%elem(nzero) =-damping_horiz; 

                END DO
             END DO
          END DO


          DO iz = 1,Nz-1
             DO ix = 1,Nx		
                DO ip=1,Np
                   irow                  = irow  + 1;
                   nzero                 = nzero + 1;
                   TRHS(irow)            = damping_horiz*(Tslow%val(ip,ix,iz+1)-Tslow%val(ip,ix,iz)); !0.0_PREC;
                   Tmatrix%irow(irow)    = nzero;

                   Tmatrix%icol(nzero) = ((ix-1) + (iz-1)*Nx)*Np + ip;
                   Tmatrix%elem(nzero) = damping_horiz; 

                   nzero = nzero + 1;
                   Tmatrix%icol(nzero) = ((ix-1) + (iz-1)*Nx)*Np + ip + Nx*Np;
                   Tmatrix%elem(nzero) =-damping_horiz; 

                END DO
             END DO
          END DO


       END IF

       Tmatrix%nrow         = irow;
       Tmatrix%nzero        = nzero;
       Tmatrix%irow(irow+1) = nzero+1;

       linsys_solver: SELECT CASE(method)

       CASE (1)

          Nitmax = 15*Ncol;  ! SIRT ain't famous for its efficiency
          CALL sirt_solver(TMatrix,Tsol,TRHS,res,nitmax,niter)

       CASE (2)

          Nitmax = 10*Ncol;  ! This choice is very  conservative
          CALL cg_solver(TMatrix,Tsol,TRHS,res,nitmax,niter)

       CASE (3)

          Nitmax = 10*Ncol;  ! This choice is very  conservative
          CALL lsqr_solver(TMatrix,Tsol,TRHS,res,nitmax,niter)

       CASE default

          Nitmax = 10*Ncol;  ! This choice is very conservative
          CALL cg_solver(TMatrix,Tsol,TRHS,res,nitmax,niter)

       END SELECT linsys_solver

       IF(niter >= nitmax) WRITE(*,*) 'Maximum number of Iterations reached in lin_sys_SOLVER res = ',res


       PRINT *, 'L1 = ',SUM(ABS(Tsol)),'    res = ',tmg_rsd_std_0 
       Tsol = alpha * Tsol
       CALL add_perturbation(Tsol)

       CALL get_parameters(Tsol)

       ! Add the perturbation to the current model


       DO iz=1,Tslow%Nz
          DO ix=1,Tslow%Nx
             DO ip=1,Tslow%Np
                icol = ip + ((ix-1) + (iz-1)*Nx)*Np
                Tslow%val(ip,ix,iz) = Tsol(icol)
             END DO
          END DO
       END DO

    END DO NL_tomog_loop

    CALL get_model(Tsol)

    ! Add the perturbation to the current model

    DO ip=1,Tslow%Np
       DO iz=1,Tslow%Nz
          DO ix=1,Tslow%Nx
             icol = ip + ((ix-1) + (iz-1)*Nx)*Np
             Tslow%val(ip,ix,iz) = Tsol(icol)
          END DO
       END DO
    END DO



    ! END DO Cont_loop

    WRITE(*,*) 'Program StraightRay_Tomog finished with NO ERRORS '
    WRITE(*,*) 'Root Mean Square of Residuals  :: ',tmg_rsd_std;
    WRITE(*,*) 'Number of Nonlinear iterations :: ',n_iter
    WRITE(*,*) 'STOPING criteria               :: ',istop

  END SUBROUTINE tomog_solve

  SUBROUTINE tomog_residuals()
    INTEGER :: ierr

    tmg_rsd_mean =  SUM(TRHS(1:Tgeom%nrays))/float(Tgeom%nrays);

    tmg_rsd_std  =  SQRT(DOT_PRODUCT((TRHS(1:Tgeom%nrays)-tmg_rsd_mean), &
         (TRHS(1:Tgeom%nrays)-tmg_rsd_mean)  &
         )/float(Tgeom%nrays-1)               &
         );

  END SUBROUTINE tomog_residuals


  SUBROUTINE raytracing()
    USE sparse
    USE raysmooth_2d

    REAL(PREC),DIMENSION(:),ALLOCATABLE::Trow
    REAL(PREC) ::xs0,zs0,xr,zr,xs,zs,length,x,z,s,ds,dsx,dsz,tvt
    REAL(PREC) ::phi,phi_min,phi_max
    INTEGER    ::np,nx,nz
    INTEGER    ::ierr,isou,irec,ip,ix,iz,ixp,izp,icol,iray,index,nzero
    INTEGER    ::nsteps,irow

    ! Local Variables:

    REAL(PREC):: eps,tol,dtvt,xmin,zmin,xmax,zmax,dx,dz
    REAL(PREC):: x0(2),x1(2),n(2),xv(4),zv(4);
    INTEGER    :: is,ir,ixs,izs,ixr,izr

    ALLOCATE(Trow(Tmatrix%Ncol),STAT= ierr);
    IF(ierr /= 0) THEN
       STOP 'raytracing::No memmory available for Trow'
    END IF

    np = TSlow%Np;
    nx = TSlow%Nx;
    nz = TSlow%Nz;
    dx = TSlow%dx;
    dz = TSlow%dz;

    xmin = TSlow%x0;
    zmin = TSlow%z0;
    xmax = xmin + FLOAT(nx-1)*dx;
    zmax = zmin + FLOAT(nz-1)*dz;

    !Tolerances : 

    eps  = 0.001_prec*MIN(dx,dz);
    tol  = 10.0_prec*eps;
    ds   = 0.02_prec*min(dx,dz);

    nzero = 0;

    IF(np == 3 ) THEN ! for tilted elliptical medium

       xs0=-100000000.d0;
       zs0=-100000000.d0;

       irow = 0

       loop_on_rays_1 : DO iray=1,Tgeom%Nrays

          xs = Tgeom%xsou(iray);zs = Tgeom%zsou(iray);
          xr = Tgeom%xrec(iray);zr = Tgeom%zrec(iray);

          ! check if source and receiver are inside the grid

          tvt    = 0.0_PREC;
          rpath  = 0.0_PREC;
          Trow   = 0.0_PREC;

          IF(ABS(xs-xs0) > 0.1d0*dx .OR. ABS(zs-zs0) > 0.1d0*dz )  THEN

             xs0 = xs;
             zs0 = zs;

             IF(ABS(xs0 - xmin) < 0.05 * dx) THEN
                phi_min =-90.d0;
                phi_max = 90.d0;
                xs0 = xmin+0.01 * dx ! move source stricly inside the model 
             ELSEIF(ABS(xs0 - xmax) < 0.05 * dx) THEN
                phi_min = 90.d0;
                phi_max =270.d0;
                xs0 = xmax-0.01 * dx
             ELSEIF(ABS(zs0 - zmin) < 0.05 * dz) THEN
                phi_min = 0.d0;
                phi_max =180.d0;
                zs0 = zmin+0.01 * dz
             ELSEIF(ABS(zs0 - zmax) < 0.05 * dz) THEN
                phi_min = 180.d0;
                phi_max = 360.d0;
                zs0 = zmax-0.01 * dz
             END IF

             CALL ray_fan(xs0,zs0,phi_min,phi_max);

          END IF

          CALL interp_fan_data(xr,zr,phi,tvt,ierr)


          IF(ierr == 0 ) THEN
             irow = irow + 1
             CALL raypath(xs0,zs0,phi,ds,nsteps,rpath)
             CALL frechet_deriv(nsteps,ds,rpath,tmatrix%ncol,Trow)

             ! fill in  Tmatrix:

             Trhs(irow)       = Tobs(iray) - tvt  ;
             Tmatrix%irow(irow) = nzero+1;
             DO index = 1,Tmatrix%Ncol
                IF(ABS(Trow(index)) > 0.0) THEN
                   nzero = nzero + 1;
                   Tmatrix%icol(nzero) = index;
                   Tmatrix%elem(nzero) = Trow(index);
                END IF
             END DO

          END IF

       END DO loop_on_rays_1

       Tmatrix%nrow         = irow
       Tmatrix%nzero        = nzero
       Tmatrix%irow(irow+1) = nzero+1;

    ELSE

       xs0=-100000000.d0;
       zs0=-100000000.d0;

       irow = 0

       loop_on_rays_2 : DO iray=1,Tgeom%Nrays

          xs = Tgeom%xsou(iray);zs = Tgeom%zsou(iray);
          xr = Tgeom%xrec(iray);zr = Tgeom%zrec(iray);

          ! check if source and receiver are inside the grid

          tvt    = 0.0_PREC;
          rpath  = 0.0_PREC;
          Trow   = 0.0_PREC;

          IF(ABS(xs-xs0) > 0.1d0*dx .OR. ABS(zs-zs0) > 0.1d0*dz )  THEN

             xs0 = xs;
             zs0 = zs;

             IF(ABS(xs0 - xmin) < 0.05 * dx) THEN
                phi_min =-90.d0;
                phi_max = 90.d0;
                xs0 = xmin + 0.01 * dx ! move source stricly inside the model 
             ELSEIF(ABS(xs0 - xmax) < 0.05 * dx) THEN
                phi_min = 90.d0;
                phi_max =270.d0;
                xs0 = xmax-0.01 * dx
             ELSEIF(ABS(zs0 - zmin) < 0.05 * dz) THEN
                phi_min = 0.d0;
                phi_max =180.d0;
                zs0 = zmin+0.01 * dz
             ELSEIF(ABS(zs0 - zmax) < 0.05 * dz) THEN
                phi_min = 180.d0;
                phi_max = 360.d0;
                zs0 = zmax-0.01 * dz
             END IF

             CALL ray_fan(xs0,zs0,phi_min,phi_max);

          END IF

          CALL interp_fan_data(xr,zr,phi,tvt,ierr)

          IF(ierr == 0 ) THEN
             irow = irow + 1
             CALL raypath(xs0,zs0,phi,ds,nsteps,rpath)
             CALL frechet_deriv(nsteps,ds,rpath,tmatrix%ncol,Trow)

             ! fill in  Tmatrix:

             Trhs(irow)       = Tobs(iray) - tvt  ;
             Tmatrix%irow(irow) = nzero+1;
             DO index = 1,Tmatrix%Ncol
                IF(ABS(Trow(index)) > 0.0) THEN
                   nzero = nzero + 1;
                   Tmatrix%icol(nzero) = index;
                   Tmatrix%elem(nzero) = Trow(index);
                END IF
             END DO

          END IF

       END DO loop_on_rays_2

       Tmatrix%nrow         = irow
       Tmatrix%nzero        = nzero
       Tmatrix%irow(irow+1) = nzero+1;


    END IF

    DEALLOCATE(Trow,STAT= ierr);
    IF(ierr) THEN
       STOP 'raytracing::Memory deallocation failed  for Trow'
    END IF


  END SUBROUTINE raytracing


  ! Release Memory:
  SUBROUTINE tomog_cleanup()
    INTEGER :: ierr
    DEALLOCATE(TRHS,STAT = ierr);
    IF(ierr) THEN
       STOP 'deallocation failed for TRHS'
    END IF

    DEALLOCATE(Tsol,STAT= ierr);
    IF(ierr) THEN
       STOP 'deallocation failed for Tsol'
    END IF

    ierr = SM_Deallocate(TMatrix);
    IF(ierr) THEN
       STOP 'deallocation failed for SM_Deallocate'
    END IF

    NULLIFY(Tslow)
    NULLIFY(Tgeom)

  END SUBROUTINE tomog_cleanup



END MODULE tomog_Modeling
