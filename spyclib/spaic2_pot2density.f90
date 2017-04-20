! Written by Ulf Saalmann

MODULE spaic2_pot2density
IMPLICIT NONE

  include "spaic2_defs.f90"

  ! ALLOCATABLES
  double precision, allocatable :: parapot(:)
  double precision, allocatable :: pararho(:)
  double precision, allocatable :: bsplinepot(:, :)
  double precision, allocatable :: bsplinerho(:,:)
  double precision, allocatable :: rbump(:)
  double precision, allocatable :: dbump(:)
  double precision, allocatable :: rr(:)
  double precision, allocatable :: woodsaxon_potential(:)
  
  ! INTERFACE

  ! SUBROUTINES
  public :: linsys
  public :: initialize
  public :: diagnostics
  public :: set_woodsaxon_params
  public :: deinitialize
  public :: free_all_memory
  public :: compute_density


CONTAINS

  SUBROUTINE set_woodsaxon_params(params, n)
    integer, intent(in) :: n
    double precision, dimension(n), intent(in) :: params
    if(allocated(parapot)) then
      deallocate(parapot)
    endif
    allocate(parapot(1:n))
    parapot(:) = 2 * params(:) - 1
    if(num_para .NE. n) then
      num_para = n
    endif
    call initialize()
  END SUBROUTINE set_woodsaxon_params


  SUBROUTINE deinitialize()
  IMPLICIT NONE
    if(allocated(pararho)) then
      deallocate(pararho)
    endif
    if(allocated(bsplinepot)) then
      deallocate(bsplinepot)
    endif
    if(allocated(bsplinerho)) then
      deallocate(bsplinerho)
    endif
    if(allocated(rbump)) then
      deallocate(rbump)
    endif
    if(allocated(dbump)) then
      deallocate(dbump)
    endif
    if(allocated(rr)) then
      deallocate(rr)
    endif
    if(allocated(woodsaxon_potential)) then
      deallocate(woodsaxon_potential)
    endif
  END SUBROUTINE deinitialize


  SUBROUTINE free_all_memory()
  IMPLICIT NONE
    if(allocated(parapot)) then
      deallocate(parapot)
    endif
    call deinitialize()
  END SUBROUTINE free_all_memory


  SUBROUTINE initialize()
  IMPLICIT NONE
    integer :: i,ii
    double precision :: r

    call deinitialize()
    allocate(pararho(1:num_para))
    allocate(bsplinepot(0:nr, num_para), bsplinerho(0:nr, num_para))
    allocate(rbump(num_para), dbump(num_para))
    allocate(rr(0:nr))
    allocate(woodsaxon_potential(0:nr))

    ! compute all the b-spline coefficients
    do ii=0,nr
      r = ii*dr + 1d-8
      rr(ii) = r
      do i=1,num_para
        bsplinerho(ii,i) = bspline3norm(i,r,rpot/(num_para+2))
        bsplinepot(ii,i) = bspline3intnorm(i,r,rpot/(num_para+2))
      end do
    end do

    ! bumpy Wood-Saxon mit Coulomb tail
    do i=1,num_para
      rbump(i) = (i+3)*((rclus+drclus)/num_para)
      dbump(i) = 2*((rclus+2*drclus)/num_para)
    end do

  END SUBROUTINE initialize

  SUBROUTINE diagnostics()
    write (*, *), ''
    write (*, *), '######################'
    write (*, *), 'DIAGNOSTIC INFORMATION'
    write (*, *), '######################'
    write (*, *), 'num_para', num_para
    write (*, *), 'nr', nr
    write (*, *), 'nr', nr
    if(allocated(parapot)) then
      write (*, *), 'parapot allocated'
    else
      write (*, *), 'parapot not allocated!'
    endif
    if(allocated(pararho)) then
      write (*, *), 'pararho allocated'
    else
      write (*, *), 'parapho not allocated!'
    endif
    if(allocated(bsplinepot)) then
      write (*, *), 'bsplinepot allocated'
    else
      write (*, *), 'bsplinepot not allocated!'
    endif
    if(allocated(bsplinerho)) then
      write (*, *), 'bsplinerho allocated'
    else
      write (*, *), 'bsplinerho not allocated!'
    endif
    if(allocated(rbump)) then
      write (*, *), 'rbump allocated'
    else
      write (*, *), 'rbump not allocated!'
    endif
    write (*, *), '######################'
    write (*, *), 'DIAGNOSTIC INFORMATION'
    write (*, *), '######################'
    write (*, *), ''
  END SUBROUTINE diagnostics


  SUBROUTINE compute_density()
    integer :: ii, i1, i2
  
    double precision :: rhowosabump(0:nr)
    double precision, allocatable :: amat(:,:), bvec(:), caux(:)
    double precision :: r

    do ii=0,nr
      r = ii*dr
      rhowosabump(ii) = rhowosa(r) &
           + vclus*fluc*sum(parapot(:)*(4*r*(r-rbump(:))**2+dbump(:)**2*(-6*r+4*rbump(:))) &
           / (dbump(:)**4*exp((r-rbump(:))**2/dbump(:)**2)*r))
    end do
  
    allocate(amat(num_para+1,num_para+1), bvec(num_para+1), caux(num_para+1))
    bvec(:) = 1
    amat(:,:) = 1
    amat(num_para+1,num_para+1) = 0
    do i1=1,num_para
      bvec(i1) = sum(rr(1:)*rhowosabump(1:)*bsplinerho(1:,i1))
      do i2 = i1, num_para
        amat(i1,i2) = sum(rr(1:)*bsplinerho(1:,i1)*bsplinerho(1:,i2))
        amat(i2,i1) = amat(i1,i2)
      end do
    end do
    call linsys(amat, bvec, caux)
    pararho(1:num_para) = caux(1:num_para)

    ! Compute Wood-Saxon-based potential
    woodsaxon_potential(:) = 0
    do ii = 0, nr
      r = ii*dr
      woodsaxon_potential(ii) = potwosa(r) + vclus*fluc*sum(parapot(:)*exp(-(r-rbump(:))**2/dbump(:)**2))
    end do

    deallocate(amat)
    deallocate(bvec)
    deallocate(caux)

  END SUBROUTINE compute_density


  FUNCTION potwosa(r)
  IMPLICIT NONE
    double precision, intent(in) :: r
    double precision :: potwosa
    potwosa = -vclus/(1+exp((r-rclus)/drclus)) - (1/r)/(1+exp((rclus-r)/drclus))
  END FUNCTION potwosa

  
  FUNCTION rhowosa(r)
  IMPLICIT NONE
    double precision, intent(in) :: r
    double precision :: rhowosa
    rhowosa = (1/cosh((r-rclus)/(2*drclus))**2*(2*drclus*vclus &
        +(1-r*vclus)*tanh((r-rclus)/(2*drclus))))/(4*drclus**2*r)
  END FUNCTION rhowosa
  
  
  FUNCTION bspline3(j,r,dr)
    implicit none
    integer, intent(in) :: j
    double precision, intent(in) :: r,dr
    double precision :: bspline3
    double precision :: raux
    raux = (r - (j-1)*dr)/dr
    bspline3 = 0
    if (0<raux .and. raux<3) then
      select case (nint(raux+5d-1))
      case(1)
        bspline3 = raux**2/2
      case(2)
        bspline3 = (2*raux*(3-raux)-3)/2
      case(3)
        bspline3 = (raux-3)**2/2
      end select
    end if
  END FUNCTION bspline3


  FUNCTION bspline3int(j,r,dr)
    implicit none
    integer, intent(in) :: j
    double precision, intent(in) :: r,dr
    double precision :: bspline3int
    double precision :: raux,pota,potb
    raux = (r - (j-1)*dr)/dr
    bspline3int = 0
    select case(nint(raux+5d-1))
    case(:0)
      pota = 0
      potb = 12 + 24*j
    case(1)
      pota = raux**3*(10*(-1+j)**2+15*(-1+j)*raux+6*raux**2)
      potb = 12+(4-3*raux)*raux**3-4*j*(-6+raux**3)
    case(2)
      pota = 1+j*(-5+j*(10+j*(-10+j*(5+2*j)))) &
           -10*(-1+2*j*(1+j))*(-1+j+raux)**3+15*(1+2*j)*(-1+j+raux)**4-12*(-1+j+raux)**5
      potb = 21+4*j*(3+(-3+raux)*raux*(-3+2*raux))+2*(-3+raux)*raux*(6+raux*(-7+3*raux))
    case(3)
      pota = -2-j*(20+j*(20+j*(40+j*(10+j)))) &
           +10*(2+j)**2*(-1+j+raux)**3-15*(2+j)*(-1+j+raux)**4+6*(-1+j+raux)**5
      potb = -((-3+raux)**3*(-1+4*j+3*raux))
    case(4:)
      pota = 30*(1+2*j*(1+j))
      potb = 0
    end select
    bspline3int = - ((dr**2/(60*(-1+j+raux))) * pota + (dr**2/24) * potb)
  END FUNCTION bspline3int

  
  FUNCTION bspline3norm(j,r,dr)
    implicit none
    integer, intent(in) :: j
    double precision, intent(in) :: r,dr
    double precision :: bspline3norm
    bspline3norm = bspline3(j,r,dr) / ((dr**3/2) * (1+2*j*(1+j)))
  END FUNCTION bspline3norm


  FUNCTION bspline3intnorm(j,r,dr)
    implicit none
    integer, intent(in) :: j
    double precision, intent(in) :: r,dr
    double precision :: bspline3intnorm
    bspline3intnorm = bspline3int(j,r,dr) / ((dr**3/2) * (1+2*j*(1+j)))
  END FUNCTION bspline3intnorm


  
  SUBROUTINE linsys(a,b,x)
    implicit none
    double precision, dimension(:,:), intent(in) :: a
    double precision, dimension(:), intent(in) :: b
    double precision, dimension(:), intent(out) :: x

    integer     ::  info,n,nrhs,ipiv(size(a,1))
    double precision, dimension(size(a,1),size(a,2))    ::  aa
    double precision, dimension(size(b),1)    ::  bb

    n = size(b,1)
    nrhs = 1

    aa = a
    bb(:,1) = b
    call dgesv(n,nrhs,aa,n,ipiv,bb,n,info)

    x = bb(:,1)
  END SUBROUTINE linsys

END MODULE spaic2_pot2density
