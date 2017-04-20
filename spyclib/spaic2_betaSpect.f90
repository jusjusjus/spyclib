! Written by Ulf Saalmann
MODULE spaic2_betaSpect
IMPLICIT NONE

  ! GLOBAL PARAMETERS
  include "spaic2_defs.f90"

  integer :: num_qns

  double precision :: freqmin=0.25d0, freqmax=1.25d0
  integer, parameter :: nfreq=100

  ! ALLOCATABLES
  double precision, allocatable :: knots(:),pots(:,:),pot0(:),pot1(:),pot(:),bpot(:),cpot(:)
  double precision, allocatable :: rr(:),bpsi(:),cpsi(:), potcha(:),rhos(:,:)
  double precision, allocatable :: freq(:)
  double precision, allocatable :: spec(:, :), beta(:, :)
  double precision, allocatable :: para(:)
  integer, allocatable :: quantum_numbers(:, :)

  ! INTERFACE

  ! SUBROUTINES
  public :: initialize
  public :: crossection
  public :: deinitialize
  public :: compute_beta
  public :: free_all_memory
  public :: set_potential_and_qns


CONTAINS

  SUBROUTINE free_all_memory()
  IMPLICIT NONE
    call deinitialize()
    if(allocated(para)) then
      deallocate(para)
    endif
  END SUBROUTINE free_all_memory


  SUBROUTINE set_potential_and_qns(params, params_size, qns, qns_size)
    integer, intent(in) :: params_size, qns_size
    double precision, dimension(params_size), intent(in) :: params
    integer, dimension(qns_size, 2), intent(in) :: qns

    ! write (*, *), qns(:, :)
    ! para = params
    if(allocated(para)) then
      deallocate(para)
    endif
    allocate(para(1:params_size))
    para(:) = params(:)
    num_para = params_size

    ! quantum_numbers = qns
    if(allocated(quantum_numbers)) then
      deallocate(quantum_numbers)
    endif
    allocate(quantum_numbers(1:qns_size, 1:2))
    quantum_numbers(:, :) = qns(:, :)
    num_qns = qns_size

    call initialize()
  END SUBROUTINE set_potential_and_qns


  SUBROUTINE compute_beta()
    integer :: ii, l, n
    double precision :: econt(0:nfreq), dipm(0:nfreq), dipp(0:nfreq), &
      phim(0:nfreq), phip(0:nfreq)

    ! calculate differential cross section, beta parameter, etc.
    do ii=1,num_qns
      l = quantum_numbers(ii, 1)
      n = quantum_numbers(ii, 2)
      if (l==0) then
        dipm(:) = 0
        call crossection(l, n, l+1, freq, nfreq, econt,dipp,phip)
      else
        call crossection(l, n, l-1, freq, nfreq, econt,dipm,phim)
        call crossection(l, n, l+1, freq, nfreq, econt,dipp,phip)
      end if
      
      beta(:, ii) = (l*(l-1) * dipm(:)**2 + (l+1)*(l+2) * dipp(:)**2 &
           - 6*l*(l+1) * dipm(:)*dipp(:) * cos(phip(:)-phim(:))) &
           / ((2*l+1) * (l*dipm(:)**2 + (l+1)*dipp(:)**2) + 1d-99)
      spec(:, ii) = dipm(:)**2 + dipp(:)**2

    end do

  END SUBROUTINE compute_beta

  include "spaic2_misc.f90"
  include "spaic2_cwf.f90"

  SUBROUTINE crossection(lin, nin, lfin, freq, nfreq, econt, dip, phi)
  IMPLICIT NONE
    integer, intent(in) :: lin,nin,lfin, nfreq
    double precision, intent(in) :: freq(0:nfreq) 
    double precision, intent(out) :: econt(0:nfreq),dip(0:nfreq),phi(0:nfreq) 
    integer :: i, j
    double precision :: eb, kc, angular,s,c
 
    j = max(lin,lfin)
    angular = j/sqrt(4*j*j-1d0)

    bpot(:) = pot(:) + lin*(lin+1)*pot1(:)
    cpot(:) = pot(:) + lfin*(lfin+1)*pot1(:)

    ! calculate bound state
    call bdiag(nin, bpot(:), nr,eb,bpsi)


    if (eb+freq(0)<ecmin) &
         call warning("some frequencies are too small to ionize the sysmtem")

    ! calculate continum states and dipole matrix element
    do i=0, nfreq
      econt(i) = eb+freq(i)
      kc = sqrt(2*econt(i))
      if (econt(i)<ecmin) then
        dip(i) = 0
        phi(i) = 0
        cpsi(:) = 0
        s = 1
        c = 0
      else
        call cnumerov(lfin,econt(i),cpot(0:), nr,cpsi,s,c)
        !! write(0,'(3i1,i4,f8.3,2f8.4)') lin,nin,lfin,i,27.21*econt(i),s,c
        phi(i) = atan2(c,s)
        dip(i) = angular * sum(bpsi(:)*rr(:)*cpsi(:))/sqrt(s**2+c**2)
      end if
      
    end do

  END SUBROUTINE crossection


  subroutine initialize()
    implicit none
    integer :: i,ii
    double precision :: r,drbspline

    call deinitialize()

    allocate(pots(0:nr, 1:num_para))
    allocate(rr(0:nr), bpot(0:nr), cpot(0:nr), bpsi(0:nr), cpsi(0:nr))
    allocate(pot(0:nr), pot0(0:nr), pot1(0:nr))
    allocate(spec(0:nfreq, 1:num_qns), beta(0:nfreq, 1:num_qns))
    allocate(potcha(1:num_para), rhos(0:nr, 1:num_para))
    allocate(freq(0:nfreq))

    ! Compute frequency axis
    do i = 0, nfreq
      freq(i) = freqmin + i * ((freqmax-freqmin)/nfreq)
    end do
    
    ! calculate all potentials
    do ii=0,nr
      r = ii*dr + 1d-6
      rr(ii) = r
      if (ii>0) pot1(ii) = 1/(2*r**2)
    end do
    pot1(0) = pot1(1)

    ! make the knot distance that the last b-spline vanishes at rpot
    drbspline = rpot/(num_para+2)
    do ii = 0, nr
      r = rr(ii) + 1d-10
      pot0(ii) =0
      do i=1,num_para
        pots(ii,i) = bspline3intnorm(i,r,drbspline)
      end do
    end do

    pot0(:) = 0
    pot(:)  = 0
    do i=1,num_para
      ! shift the parameter's c.o.m. to guarantee sum(para(:)) = 1
      pot(:) = pot(:) + (para(i)-(sum(para(:))-1)/num_para)*pots(:,i)
      ! scale all parameters to guarantee sum(para(:)) = 1
      ! pot(:) = pot(:) + (para(i)/sum(para(:)))*pots(:,i)
    end do

  END SUBROUTINE initialize


  SUBROUTINE deinitialize()
  IMPLICIT NONE
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
    if(allocated(potcha)) then
      deallocate(potcha)
    endif
    if(allocated(rhos)) then
      deallocate(rhos)
    endif
    if(allocated(spec)) then
      deallocate(spec)
    endif
    if(allocated(beta)) then
      deallocate(beta)
    endif
    if(allocated(freq)) then
      deallocate(freq)
    endif
  END SUBROUTINE deinitialize


END MODULE spaic2_betaSpect
