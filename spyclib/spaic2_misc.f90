! Written by Ulf Saalmann

SUBROUTINE warning(s)
  ! issue a warning, but do not stop calcuations
  implicit none
  character*(*), intent(in) :: s
  write(*,'(a)') s
END SUBROUTINE warning

! =========
! B-splines
! =========

! 3rd-order splines with an equi-distant knot sequence: x_j=(j-1)*dr
function bspline3(j,r,dr)
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
end function bspline3

! integration of the radial Poisson equation with a charge distribution
! given by 3rd-order splines with an equi-distant knot sequence: x_j=(j-1)*dr
! i.e.  (1/r) int_0^r ds s^2 b(s) + int_r^infty ds s b(s) 
function bspline3int(j,r,dr)
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
end function bspline3int

! same as bspline3 but normalized: int dr r^2 b(r) = 1
function bspline3norm(j,r,dr)
  implicit none
  integer, intent(in) :: j
  double precision, intent(in) :: r,dr
  double precision :: bspline3norm
  bspline3norm = bspline3(j,r,dr) / ((dr**3/2) * (1+2*j*(1+j)))
end function bspline3norm

! same as bspline3int but normalized: int dr r^2 b(r) = 1
function bspline3intnorm(j,r,dr)
  implicit none
  integer, intent(in) :: j
  double precision, intent(in) :: r,dr
  double precision :: bspline3intnorm
  bspline3intnorm = bspline3int(j,r,dr) / ((dr**3/2) * (1+2*j*(1+j)))
end function bspline3intnorm


! ==============
! numerov method
! ==============

subroutine cnumerov(l,e,pot, nr,psi,s,c)
  implicit none
  integer, intent(in) :: l, nr
  double precision, intent(in) :: e,pot(0:nr)
  double precision, intent(out) :: psi(0:nr),s,c

  integer :: i, j
  double precision :: r
  double precision :: ps,pc,ss,sc,cc
  double precision, allocatable, save :: psisin(:),psicos(:)

  logical, save :: first=.true.

  if (first) then
    allocate(psicos(nasymp:nr),psisin(nasymp:nr))
    first = .false.
  end if

  call integ(l,e,pot, nr,j,psi,.false.)

  ! fit to asymtoptic solutions (for normalization and phase shift)
  if (sqrt(2*e)*((nr-nasymp)*dr)<pi) &
       call warning("fitting region is rather small")

  do i=nasymp,nr
    r = i*dr
    psisin(i) = asympsin(l,e,r)
    psicos(i) = asympcos(l,e,r)
  end do

  ps = sum(psisin(:)*psi(nasymp:nr))
  pc = sum(psicos(:)*psi(nasymp:nr))
  ss = sum(psisin(:)*psisin(:))
  sc = sum(psisin(:)*psicos(:))
  cc = sum(psicos(:)*psicos(:))
  s = (cc*ps-sc*pc)/(cc*ss-sc**2)
  c = (ss*pc-sc*ps)/(cc*ss-sc**2)
  
end subroutine cnumerov


! sin-like asymptotic function
function asympsin(l,e,r)
  implicit none
  integer, intent(in) :: l
  double precision, intent(in) :: e,r
  double precision :: asympsin
  if (nint(r/dr)<nasymp) then
    asympsin = 0
    return
  end if
  asympsin = coulombf(l,e,r)
end function asympsin

! cos-like asymptotic function
function asympcos(l,e,r)
  implicit none
  integer, intent(in) :: l
  double precision, intent(in) :: e,r
  double precision :: asympcos
  if (nint(r/dr)<nasymp) then
    asympcos = 0
    return
  end if
  asympcos = coulombg(l,e,r)
end function asympcos


! ======================
! Coulomb wave functions
! ======================

function coulombf(l,e,r)
  ! regular
  implicit none
  integer, intent(in) :: l
  double precision, intent(in) :: e,r
  double precision :: coulombf

  integer, save :: lx=-1
  double precision, save :: ex=0
  double precision, save :: cwf(nasymp:nr)

  integer :: i, j
  double precision :: drq,raux(0:2)

  if (l.ne.lx .or. abs(e-ex)>1d-6) then
    cwf(nasymp) = cf1(l,e)
    cwf(nasymp+1) = cf2(l,e)

    drq = dr**2
    raux(2) = drq * (-1/rr(nasymp)+l*(l+1)*pot1(nasymp) - e)
    raux(1) = drq * (-1/rr(nasymp+1)+l*(l+1)*pot1(nasymp+1) - e)
    do i=nasymp+2,nr
      raux(0) = drq * (-1/rr(i)+l*(l+1)*pot1(i) - e)
      cwf(i) = ( (12+10*raux(1))*cwf(i-1) - (6-raux(2))*cwf(i-2) ) / (6-raux(0))
      raux(2) = raux(1)
      raux(1) = raux(0)
    end do
  end if

  j = nint(r/dr)
  if (j<nasymp .or. j>nr) then
    write(0,'("coulombf: r is outside allowed regin")')
    stop
  end if
  coulombf = cwf(j)

end function coulombf


FUNCTION coulombg(l,e,r)
  ! irregular
IMPLICIT NONE
  integer, intent(in) :: l
  double precision, intent(in) :: e,r
  double precision :: coulombg

  integer, save :: lx=-1
  double precision, save :: ex=0
  double precision, save :: cwf(nasymp:nr)

  integer :: i, j
  double precision :: drq,raux(0:2)

  if (l.ne.lx .or. abs(e-ex)>1d-6) then
    cwf(nasymp) = cg1(l,e)
    cwf(nasymp+1) = cg2(l,e)

    drq = dr**2
    raux(2) = drq * (-1/rr(nasymp)+l*(l+1)*pot1(nasymp) - e)
    raux(1) = drq * (-1/rr(nasymp+1)+l*(l+1)*pot1(nasymp+1) - e)
    do i = nasymp+2, nr
      raux(0) = drq * (-1/rr(i)+l*(l+1)*pot1(i) - e)
      cwf(i) = ( (12+10*raux(1))*cwf(i-1) - (6-raux(2))*cwf(i-2) ) / (6-raux(0))
      raux(2) = raux(1)
      raux(1) = raux(0)
    end do
  end if

  j = nint(r/dr)
  if (j<nasymp .or. j>nr) then
    write(0,'("coulombg: r is outside allowed regin")')
    stop
  end if
  coulombg = cwf(j)

END FUNCTION coulombg


! ==========================
! spherical Bessel functions
! ==========================

function sbj(k,x)
  ! spherical bessel j
  implicit none
  integer, intent(in) :: k
  double precision, intent(in) :: x
  double precision :: sbj
  select case(k)
  case(0)
    sbj = sin(x)/x
  case(1) 
    sbj = (-(x*cos(x)) + sin(x))/x**2
  case(2) 
    sbj = -((3*x*cos(x) + (-3 + x**2)*sin(x))/x**3)
  case(3) 
    sbj = (x*(-15 + x**2)*cos(x) + 3*(5 - 2*x**2)*sin(x))/x**4
  case(4)
    sbj = (5*x*(-21 + 2*x**2)*cos(x) + (105 - 45*x**2 + x**4)*sin(x))/x**5
  case(5)
    sbj = (-(x*(945 - 105*x**2 + x**4)*cos(x)) + 15*(63 - 28*x**2 + x**4)*sin(x))/x**6
  case default
    write(0,'("sbj: k>5")')
  end select
end function sbj


FUNCTION sby(k,x)
  ! spherical bessel j
  implicit none
  integer, intent(in) :: k
  double precision, intent(in) :: x
  double precision :: sby
  select case(k)
  case(0)
    sby = -(cos(x)/x)
  case(1) 
    sby = -((cos(x) + x*sin(x))/x**2)
  case(2) 
    sby = ((-3 + x**2)*cos(x) - 3*x*sin(x))/x**3
  case(3) 
    sby = (3*(-5 + 2*x**2)*cos(x) + x*(-15 + x**2)*sin(x))/x**4
  case(4) 
    sby = -(((105 - 45*x**2 + x**4)*cos(x) + 5*x*(21 - 2*x**2)*sin(x))/x**5)
  case(5)
    sby = -((15*(63 - 28*x**2 + x**4)*cos(x) + x*(945 - 105*x**2 + x**4)*sin(x))/x**6)
  case default
    sby = 0.0d0
    write(0,'("sby: k>5")')
  end select
END FUNCTION sby



! integration for Numerov method
SUBROUTINE integ(l,e,pot, nr,nzeros,psi,qprint)
IMPLICIT NONE
  integer, intent(in) :: l, nr
  double precision, intent(in) :: e
  double precision, intent(in) :: pot(0:nr)
  integer, intent(out) :: nzeros
  double precision, intent(out) :: psi(0:nr)
  logical, intent(in) :: qprint

  integer :: i,j0,flag,j
  double precision :: r, rtp, drq
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
  if (qprint) write(*,'("# e rtp",f8.4,f8.2)') e ,rtp

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
    ! write(*,'("# j0",i3)') j0
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

  ! write(*,'(f6.2,2f10.6)') (i*dr,pot(i),psi(i),i=0,nr); stop

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


! diagonalization for bound state
SUBROUTINE bdiag(k, pot, nr,e,psi)
IMPLICIT NONE
  integer, intent(in) :: k, nr
  double precision, intent(in) :: pot(0:nr)
  double precision, intent(out) :: e,psi(0:nr)

  double precision :: raux
  double precision :: hmat(0:3,1:nasymp),vaux(1:nasymp)

  raux = 1/(360*dr**2)
  hmat(0,1:nasymp) =  490*raux + pot(1:nasymp)    
  hmat(1,1:nasymp) = -270*raux
  hmat(2,1:nasymp) =   27*raux
  hmat(3,1:nasymp) =   -2*raux

  call bandeiglevs(nasymp,3,hmat,vaux)
  ! write(*,'("# lev",99f8.4)') vaux(1:min(10,nasymp))
  e = vaux(k)
  call bandeigvec(nasymp,3,hmat,e,vaux)
  if (vaux(2)>0) then
    psi(1:nasymp) = vaux(1:nasymp)
  else
    psi(1:nasymp) = -vaux(1:nasymp)
  end if
  psi(nasymp+1:) = 0

end subroutine bdiag




! ==============
! linear algebra
! ==============

subroutine bandeiglevs(n,ns,a,lev)
  implicit none
  integer,intent(in)::n,ns
  double precision,intent(in)::a(0:ns,n)
  double precision,intent(out)::lev(n)

  integer :: INFo
  double precision::M(ns+1,n),W(n),Z(1,1),WORK(3*n)

  M=a

  call DSBEV('N','L',n,ns,M,ns+1,W,Z,1,WORK,INFo)
  if(INFo<0)then
    write(0,*)"bandeiglevs|DSBEV: an argument had an illegal value"
    stop
  end if
  if(INFo>0)then
    write(0,*)"bandeiglevs|DSBEV: algorithm failed to converge, i=",INFo
    stop
  end if
  lev=W
end subroutine bandeiglevs



subroutine bandeigvec(n,ns,a,lev,vec)
  implicit none
  integer,intent(in)::n,ns
  double precision,intent(in)::a(0:ns,n)
  double precision,intent(in)::lev
  double precision,intent(out)::vec(n)
  logical::check=.FALSE.

  integer::i,j,ii, INFo,IPIV(n)
  double precision::AB(2*ns+ns+1,n),B(n,1),C(n),raux

  B=1d0/n

  do ii=1,3
    AB=0
    do j=1,n
      do i=max(1,j-ns),min(n,j+ns)
        AB(ns+ns+1+i-j,j)=a(abs(i-j),j)
        if(i==j)AB(ns+ns+1+i-j,j)=AB(ns+ns+1+i-j,j)-lev
      end do
    end do
    if(check)C(1:n)=B(1:n,1)
    call DGBSV(n,ns,ns,1,AB,2*ns+ns+1,IPIV,B,n,INFo)
    raux=1d0/sqrt(dot_product(B(1:n,1),B(1:n,1)))
    B(1:n,1)=raux*B(1:n,1)

    if(check)then
      raux=1/dot_product(B(1:n,1),C(1:n))
      C=0
      do j=1,n
        do i=max(1,j-ns),min(n,j+ns)
          C(j)=C(j)+a(abs(i-j),j)*B(i,1)
        end do
      end do
      print*,"XXX",lev,dot_product(B(1:n,1),C(1:n)),1+raux
    end if

  end do

  vec(1:n)=B(1:n,1)
end subroutine bandeigvec

