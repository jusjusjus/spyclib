! Written by Ulf Saalmann

MODULE spaic
  IMPLICIT NONE

  ! Global parameters
  double precision, parameter :: pi = 3.141592653589793238462
  integer :: npara, jr, nr
  double precision :: rclus = 14.478
  double precision :: drclus = 0.9
  double precision :: vclus = 0.226
  double precision :: fluc = 0.3
  integer :: lin = 0 ! s
  integer :: kin = 1 ! 1
  integer :: dl  = 1 ! l -> l+1
  integer :: nw = 100 ! num frequencies
  double precision :: wmin = 0.2
  double precision :: wmax = 2.2
  double precision, dimension(:), allocatable :: para

  ! hard-wired parameter
  integer :: jpot=1
  double precision :: dr = 0.05, rmax = 200

  ! Internal parameters
  double precision, allocatable :: knots(:), pots(:,:), pot0(:), pot1(:), pot(:), bpot(:), cpot(:)
  double precision, allocatable :: rr(:), bpsi(:), cpsi(:), frequencies(:), spectrum(:)

  ! Function and Subroutines
  public :: set_params
  public :: bspline
  public :: cnumerov
  public :: bnumerov
  public :: integ
  public :: wkbratio
  public :: compute_crosssection
  public :: initialize
  public :: deinitialize
  public :: free_all_memory


CONTAINS

  ! SUBROUTINE deallocate_dv(mem)
  !   double precision, dimension(:), allocatable, intent(inout) :: mem
  !   if(allocated(mem)) then
  !     deallocate(mem)
  !   endif
  ! END SUBROUTINE deallocate_dv

  SUBROUTINE set_params(params, n)
    integer, intent(in) :: n
    double precision, dimension(n), intent(in) :: params
    double precision, allocatable :: para(:)
    if(allocated(para)) then
      deallocate(para)
    endif
    allocate(para(1:n))
    para(:) = params(:)
    if(npara .NE. n) then
      npara = n
      call initialize()
    endif
  END SUBROUTINE set_params


  SUBROUTINE compute_crosssection()
    call crossection(nr)
  END SUBROUTINE compute_crosssection


  subroutine crossection(nr)
    implicit none
    integer, intent(in) :: nr
    integer :: i
    double precision :: w, eb

    if (lin+dl<0) then
      write (*, *) "negative angular momentum: lin+dl"
      stop
    end if

    ! potential
    pot(:) = pot0(:)
    do i=1,npara
      pot(:) = pot(:) +  (2*para(i)-1)*pots(:,i)
    end do
    bpot(:) = pot(:) + lin*(lin+1)*pot1(:)
    cpot(:) = pot(:) + (lin+dl)*(lin+dl+1)*pot1(:)

    ! bound state
    call bnumerov(lin,kin,bpot(:),eb,bpsi, nr)

    ! continuum states
    do i=1,nr
      if (abs(pot(i)+(lin+dl)*(lin+dl+1)*pot1(i))>1d-3)  jr=i
    end do
    if (nr-jr<10) then
      write (*, *) "too few asymptotic points: ", nr-jr
      stop
    end if

    if (eb+wmin<0) then
      write (*, *) "freqeuncy too small: ", nr-jr, eb+wmin
      stop
    end if

    do i=0,nw
      w = wmin + i*((wmax-wmin)/nw)
      call cnumerov(lin+dl, eb+w,cpot(0:), cpsi, nr)
      frequencies(i) = w
      spectrum(i) = (sum(bpsi(:)*rr(:)*cpsi(:)))**2
    end do

  end subroutine crossection


  SUBROUTINE deinitialize()
  IMPLICIT NONE
    if(allocated(knots)) then
      deallocate(knots)
    endif
    if(allocated(pots)) then
      deallocate(pots)
    endif
    if(allocated(rr)) then
      deallocate(rr)
    endif
    if(allocated(bpot)) then
      deallocate(bpot)
    endif
    if(allocated(cpot)) then
      deallocate(cpot)
    endif
    if(allocated(bpsi)) then
      deallocate(bpsi)
    endif
    if(allocated(cpsi)) then
      deallocate(cpsi)
    endif
    if(allocated(pot)) then
      deallocate(pot)
    endif
    if(allocated(pot0)) then
      deallocate(pot0)
    endif
    if(allocated(pot1)) then
      deallocate(pot1)
    endif
    if(allocated(frequencies)) then
      deallocate(frequencies)
    endif
    if(allocated(spectrum)) then
      deallocate(spectrum)
    endif
  END SUBROUTINE deinitialize


  SUBROUTINE free_all_memory()
  IMPLICIT NONE
    call deinitialize()
    if(allocated(para)) then
      deallocate(para)
    endif
  END SUBROUTINE free_all_memory

  subroutine initialize()
    implicit none
    integer :: i,ii
    double precision :: r

    nr = nint(rmax/dr)

    call deinitialize()
    allocate(knots(1:npara+1))   
    allocate(pots(0:nr,1:npara))
    allocate(rr(0:nr),bpot(0:nr),cpot(0:nr),bpsi(0:nr),cpsi(0:nr))
    allocate(pot(0:nr),pot0(0:nr),pot1(0:nr))
    allocate(frequencies(0:nw))
    allocate(spectrum(0:nw))

    ! define the knots for b-splines
    select case(jpot)
    case(1)
      do i=1,npara+1
        knots(i) = (i-1)*((rclus+2*drclus)/npara)
      end do
    end select

    ! calculate all potentials
    do ii=0,nr
      r = ii*dr
      rr(ii) = r
      pot0(ii) = -vclus/(1+exp((r-rclus)/drclus))
      if (ii>0) &
        pot1(ii) = 1/(2*r**2)
      do i=1,npara
        pots(ii,i) = fluc*vclus*bspline(i,r,npara+2,knots(:))
      end do
    end do
    pot1(0) = pot1(1)

  end subroutine initialize

  ! =========
  ! B-splines
  ! =========

  function bspline(j, x, n, xx)
    implicit none
    integer, intent(in) :: j,n
    double precision, intent(in) :: x,xx(1:n-1)
    double precision :: bspline

    bspline = 0
    if (j==1) then
      if (xx(1)<=x .and. x<xx(2)) &
        bspline = ((x-xx(2))/(xx(1)-xx(2)))**2
    end if
    if (j==2) then
      if (xx(1)<=x .and. x<xx(2)) &
        bspline = -(((x - xx(1))*(x*(2*xx(1) - xx(2) - xx(3)) &
        + 2*xx(2)*xx(3) - xx(1)*(xx(2) + xx(3)))) &
        / ((xx(1) - xx(2))**2*(xx(1) - xx(3))))
      if (xx(2)<=x .and. x<=xx(3)) &
        bspline = (x - xx(3))**2 &
        / ((-xx(1) + xx(3))*(-xx(2) + xx(3)))
    end if
    if (2<j .and. j<n-1) then
      if (xx(j-2)<=x .and. x<xx(j-1)) &
        bspline = (x - xx(j-2))**2/((xx(j-2) - xx(j-1))*(xx(j-2) - xx(j)))
      if (xx(j-1)<=x .and. x<xx(j)) &
        bspline = (((x - xx(j-2))*(x - xx(j)))/(xx(j-2) - xx(j)) &
        + ((x - xx(j-1))*(x - xx(j+1)))/(xx(j-1) - xx(j+1)))/(-xx(j-1) + xx(j))
      if (xx(j)<=x .and. x<xx(j+1)) &
        bspline = (x - xx(j+1))**2/((-xx(j-1) + xx(j+1))*(-xx(j) + xx(j+1)))
    end if
    if (j==n-1) then
      if (xx(n-3)<=x .and. x<xx(n-2)) &
        bspline = (x - xx(n-3))**2/((xx(n-3) - xx(n-1))*(xx(n-3) - xx(n-2)))
      if (xx(n-2)<=x .and. x<xx(n-1)) &
        bspline = -(((x - xx(n-1))*(xx(n-3)*(xx(n-1) - 2*xx(n-2)) &
        + xx(n-1)*xx(n-2) + x*(xx(n-3) - 2*xx(n-1) + xx(n-2)))) &
        / ((xx(n-3) - xx(n-1))*(xx(n-1) - xx(n-2))**2))
    end if
    if (j==n) then
      if (xx(n-2)<=x .and. x<=xx(n-1)) &
        bspline = (x - xx(n-2))**2/(xx(n-1) - xx(n-2))**2
    end if

  end function bspline



  ! ==============
  ! numerov method
  ! ==============

  subroutine cnumerov(l,e,pot,psi, nr)
    implicit none
    integer, intent(in) :: l, nr
    double precision, intent(in) :: e,pot(0:nr)
    double precision, intent(out) :: psi(0:nr)

    integer :: j
    double precision :: norm

    call integ(l,e,pot,j,psi, nr)

    norm = maxval(abs(psi(jr:nr)))
    psi(:) = psi(:)/norm
    if (psi(1)<0)  psi(:) = -psi(:)

  end subroutine cnumerov



  subroutine bnumerov(l,k,pot,e,psi, nr)
    implicit none
    integer, intent(in) :: l,k, nr
    double precision, intent(in) :: pot(0:nr)
    double precision, intent(out) :: e,psi(0:nr)

    integer :: j1,j2,j,i,nzeros
    double precision :: e1,e2

    nzeros = k-1
    if (nzeros<0) then
      write (*, *) "bnumerov: nzeros<0"
      stop 
    end if

    ! lower limit for bracketing
    e1 = -0.1
    j1 = nzeros+1
    do while(j1>nzeros)
      e1 = 2*e1
      call integ(l,e1,pot(0:nr),j1,psi, nr)
    end do

    ! upper limit for bracketing
    e2 = 1
    j2 = nzeros
    do while(j2.le.nzeros)
      e2 = 2*e2
      call integ(l,e2,pot,j2,psi, nr)
    end do

    ! naive bracketing
    do i=1,33
      ! could be improved
      e = (e1+e2)/2
      call integ(l,e,pot,j,psi, nr)
      if (j>nzeros) then
        e2 = e
      else
        e1 = e
      end if
    end do

    if (e>0 .and. sqrt(2*e)*dr>0.5) then
      write (*, *) "bnumerov: too large grid spacing? momentum:", sqrt(2*e)
    end if

  end subroutine bnumerov


  subroutine integ(l,e,pot,nzeros,psi,nr)
    implicit none
    integer, intent(in) :: l, nr
    double precision, intent(in) :: e
    double precision, intent(in) :: pot(0:nr)
    integer, intent(out) :: nzeros
    double precision, intent(out) :: psi(0:nr)

    integer :: i,j0,flag,j
    double precision :: r,rtp,drq
    double precision :: raux(0:2)

    flag = 1
    nzeros = 0
    drq = dr**2
    psi(:) = 0

    ! outer classical turning point
    j = nr
    do while (e<pot(j) .and. j>0)
      j = j-1
    end do
    rtp = j*dr

    if (pot(1)<e) then
      ! no inner classical forbidden region
      j0 = 1
      r = dr
      psi(1) = r**(l+1)*(1-r/(l+1))
      raux(2) = drq * (pot(1)-e)
      r = 2*dr
      psi(2) = r**(l+1)*(1-r/(l+1))
      raux(1) = drq * (pot(2)-e)
    else
      ! with inner classical forbidden region
      ! j0*dr is the radius where propagation should be possible
      j0 = max(1,nint(sqrt(l*(l+1)/(2*e+1/drq))/dr))
      r = j0*dr
      psi(j0) = 1d-10
      raux(2) = drq * (pot(j0)-e)
      r = (j0+1)*dr 
      ! use WKB
      psi(j0+1) = psi(j0) * wkbratio(pot(j0),pot(j0+1),dr,e)
      raux(1) = drq * (pot(j0+1)-e)

      do i=j0-1,1,-1
        psi(i) = psi(i+1)/wkbratio(pot(i),pot(i+1),dr,e)
      end do
    end if


    do i=j0+2,nr
      r = i*dr
      raux(0) = drq * (pot(i)-e)
      psi(i) = ( (12+10*raux(1))*psi(i-1) - (6-raux(2))*psi(i-2) ) / (6-raux(0))
      raux(2) = raux(1)
      raux(1) = raux(0)
      if (flag>0 .and. psi(i)<0) then
        flag = -flag 
        nzeros = nzeros+1
      end if
      if (flag<0 .and. psi(i)>0) then
        flag = -flag
        nzeros = nzeros+1
      end if
      if (e<0 .and. r>rtp+1 .and. abs(psi(i))>abs(psi(i-1))) then
        psi(i) = psi(i-1)/2
        exit
      end if
    end do


    ! normalize wave-functions
    psi(:) = psi(:)/sqrt(sum(psi(:)**2))
  end subroutine integ


  function wkbratio(pot1,pot2,dr,e)
    implicit none
    double precision, intent(in) :: pot1,pot2,dr,e
    double precision :: wkbratio
    double precision :: raux1,raux2,raux3
    raux1 = sqrt(2*(pot1-e))
    raux2 = sqrt(2*((pot1+pot2)/2-e))
    raux3 = sqrt(2*(pot2-e))
    wkbratio =  exp((dr/6)*(raux1+4*raux2+raux3))
  end function wkbratio

END MODULE spaic
